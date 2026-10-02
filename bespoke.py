import glob
import os
import select
import subprocess
import pysam
import tempfile
import time
from contextlib import redirect_stdout, redirect_stderr
import concurrent.futures
from reporting import safe_report


def merge_bed_files(bed_paths, outdir):
    merged_bed_path = tempfile.mktemp(prefix=".temp_", suffix=".bed", dir=outdir)
    with open(merged_bed_path, 'w') as merged_bed_file:
        for bed_path in bed_paths:
            with open(bed_path, 'r') as bed_file:
                for line in bed_file:
                    merged_bed_file.write(line)
    return merged_bed_path


def intersect_bam_with_bed(bam_path, merged_bed_path, outdir, reporter):
    bam_basename = os.path.splitext(os.path.basename(bam_path))[0]
    output_bam_path = os.path.join(outdir, f"roi.{bam_basename}.bam")

    reporter.info("Splitting BAM for parallel processing...")
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        contigs = bam_file.references

    def fetch_and_intersect(contig):
        temp_bam_path = tempfile.mktemp(prefix=f".temp_{contig}_", suffix=".bam", dir=outdir)
        bed_intervals = []
        with open(merged_bed_path, 'r') as bed_file:
            for line in bed_file:
                bed_contig, start, end = line.strip().split()[:3]
                if bed_contig == contig:
                    bed_intervals.append((int(start), int(end)))

        with open(os.devnull, 'w') as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                with pysam.AlignmentFile(bam_path, "rb") as bam_file, \
                    pysam.AlignmentFile(temp_bam_path, "wb", template=bam_file) as temp_bam_file:
                    written_reads = set()
                    for read in bam_file.fetch(contig):
                        if read.query_name not in written_reads:
                            for start, end in bed_intervals:
                                if start <= read.reference_start <= end:
                                    temp_bam_file.write(read)
                                    written_reads.add(read.query_name)
                                    break
        return temp_bam_path

    temp_bam_paths = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_contig = {executor.submit(fetch_and_intersect, contig): contig for contig in contigs}
        for future in concurrent.futures.as_completed(future_to_contig):
            contig = future_to_contig[future]
            try:
                temp_bam_path = future.result()
                temp_bam_paths.append(temp_bam_path)
                reporter.progress(f"Processing contig: {contig}")
            except Exception as exc:
                reporter.error(f"Contig {contig} generated an exception: {exc}")

    reporter.info("Merging temporary BAM files...")
    with open(os.devnull, 'w') as devnull:
        with redirect_stdout(devnull), redirect_stderr(devnull):
            with pysam.AlignmentFile(output_bam_path, "wb", template=pysam.AlignmentFile(bam_path, "rb")) as out_bam_file:
                for temp_bam_path in temp_bam_paths:
                    reporter.progress(f"Merging {temp_bam_path}...")
                    with pysam.AlignmentFile(temp_bam_path, "rb") as temp_bam_file:
                        for read in temp_bam_file:
                            out_bam_file.write(read)
                    os.remove(temp_bam_path)

    reporter.info("Sorting the output BAM file...")
    sorted_bam_path = output_bam_path.replace(".bam", ".sorted.bam")
    pysam.sort("-@", str(os.cpu_count()), "-o", sorted_bam_path, output_bam_path)
    os.remove(output_bam_path)

    reporter.info("Indexing the sorted BAM file...")
    pysam.index("-@", str(os.cpu_count()), sorted_bam_path)

    reporter.info(f"Output BAM: {sorted_bam_path}", color='g')


def bespoke(args, reporter):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if args.logfile:
        logfile_path = os.path.join(args.outdir, args.logfile)
    else:
        current_time = time.strftime("%Y%m%d-%H%M%S")
        logfile_path = os.path.join(args.outdir, f"{args.command}-{current_time}.log")
    reporter.set_logfile(logfile_path)

    reporter.info("Intersecting BAM with BED regions")
    merged_bed_path = merge_bed_files(args.bed, args.outdir)
    intersect_bam_with_bed(args.bam, merged_bed_path, args.outdir, reporter)
    os.remove(merged_bed_path)


def _stream_process(proc, reporter, status_message):
    pending = ["|", "/", "-", "\\"]
    while True:
        ready_fds, _, _ = select.select([proc.stdout, proc.stderr], [], [], 0.1)

        if proc.stdout in ready_fds:
            output = proc.stdout.readline()
            if output:
                reporter.stdout(output.strip())

        if proc.stderr in ready_fds:
            error = proc.stderr.readline()
            if error:
                reporter.stderr(error.strip())

        if proc.poll() is not None:
            break

        reporter.status(f"{status_message} {pending[int(time.time() * 10) % 4]}")

    for output in proc.stdout:
        reporter.stdout(output.strip())
    for error in proc.stderr:
        reporter.stderr(error.strip())


def bespoke2(args, reporter):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    merged_bed_path = None
    try:
        if args.logfile:
            logfile_path = os.path.join(args.outdir, args.logfile)
        else:
            current_time = time.strftime("%Y%m%d-%H%M%S")
            logfile_path = os.path.join(args.outdir, f"{args.command}-{current_time}.log")
        reporter.set_logfile(logfile_path)

        reporter.info("Merging BED files")
        merged_bed_path = merge_bed_files(args.bed, args.outdir)
        cmd = ["bedtools", "sort", "-i", merged_bed_path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        _stream_process(proc, reporter, "Sorting...")

        if proc.returncode != 0:
            safe_report(reporter, "error", "Intersecting BAM with BED regions failed. Exiting.")
            raise Exception("Intersecting BAM with BED regions. Exiting.")

        bam_basename = os.path.splitext(os.path.basename(args.bam))[0]
        output_bam_path = os.path.join(args.outdir, f"roi.{bam_basename}.bam")
        temp_bam_path = tempfile.mktemp(prefix=".temp_", suffix=".bam", dir=args.outdir)

        reporter.info("Intersecting BAM with BED regions")
        cmd = [
            "samtools", "view", "-bh",
            "-@", str(os.cpu_count()),
            "-L", merged_bed_path,
            # -L alone keeps only reads whose own start falls in the BED, so a
            # mate outside the target intervals is dropped. BAMSurgeon realigns
            # complete fragments, so losing mates dilutes the requested VAF and
            # pushes selected reads onto its orphan path. --fetch-pairs
            # recovers the mates; it also enables samtools' multi-region
            # iterator, which requires an indexed input BAM.
            "--fetch-pairs",
        ]
        cmd.extend(["-o", temp_bam_path, args.bam])

        reporter.command(cmd)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        _stream_process(proc, reporter, "Extracting region of interest...")

        if proc.returncode != 0:
            safe_report(reporter, "error", "Intersecting BAM with BED regions failed. Exiting.")
            raise Exception("Intersecting BAM with BED regions. Exiting.")

        reporter.info("Sorting the output BAM file...")
        pysam.sort("-@", str(os.cpu_count()), "-o", output_bam_path, temp_bam_path)
        os.remove(temp_bam_path)

        reporter.info("Indexing the sorted BAM file...")
        pysam.index("-@", str(os.cpu_count()), output_bam_path)

        reporter.info(f"Completed. Results saved to {output_bam_path}", color='g')

    except Exception as e:
        safe_report(reporter, "error", f"Error: {str(e)}")
        raise

    finally:
        safe_report(reporter, "info", "Cleaning up...")
        pattern = os.path.join(args.outdir, '*.temp*')
        files_to_delete = glob.glob(pattern)
        for file_path in files_to_delete:
            if os.path.isfile(file_path):
                try:
                    os.remove(file_path)
                except OSError as exc:
                    safe_report(reporter, "warning", f"Cleanup could not remove {file_path}: {exc}")
        if merged_bed_path and os.path.exists(merged_bed_path):
            try:
                os.remove(merged_bed_path)
            except OSError as exc:
                safe_report(reporter, "warning", f"Cleanup could not remove {merged_bed_path}: {exc}")
