import glob
import os
from pathlib import Path
import select
import subprocess
import tempfile
import time
from reporting import safe_report


def build_deepvariant_command(args):
    """Build a nested-Docker command with paths valid on the host and in the container.

    The Oracode container is launched with data directories mounted at identical
    absolute paths. Docker therefore receives one common host directory and the
    DeepVariant container sees the same files below /workdir.
    """

    paths = [
        Path(args.reference).expanduser().resolve(),
        Path(args.bam).expanduser().resolve(),
        Path(args.bed).expanduser().resolve(),
        Path(args.outdir).expanduser().resolve(),
    ]
    workspace = Path(os.path.commonpath([str(path) for path in paths]))
    if not workspace.is_dir():
        workspace = workspace.parent

    def mounted(path):
        return Path("/workdir") / path.relative_to(workspace)

    reference, bam, bed, outdir = paths
    command = [
        "docker", "run", "--rm",
        "--gpus", "all",
        "-v", f"{workspace}:/workdir",
        "google/deepvariant",
        "/opt/deepvariant/bin/run_deepvariant",
        "--model_type=WES",
        f"--ref={mounted(reference)}",
        f"--reads={mounted(bam)}",
        f"--regions={mounted(bed)}",
        f"--output_vcf={mounted(outdir) / 'mutated.vcf.gz'}",
        f"--num_shards={str(os.cpu_count())}",
    ]
    min_candidate_reads = int(getattr(args, "min_mutated_reads", 3))
    if min_candidate_reads < 1:
        raise ValueError("DeepVariant candidate alt-read count must be at least 1")
    if hasattr(args, "min_mutated_reads"):
        command.append(
            "--make_examples_extra_args="
            f"vsc_min_count_snps={min_candidate_reads},"
            f"vsc_min_count_indels={min_candidate_reads}"
        )
    return command


def educe_variant(args, reporter):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    #temp_dir = os.path.join(args.outdir, "tmp")
    #if not os.path.exists(temp_dir):
    #    os.makedirs(temp_dir)
    
    try:
        if args.logfile:
            logfile_path = os.path.join(args.outdir, args.logfile)
        else:
            current_time = time.strftime("%Y%m%d-%H%M%S")
            logfile_path = os.path.join(args.outdir, f"{args.command}-{current_time}.log")
        reporter.set_logfile(logfile_path)

        vc_cmd = []

        if args.tool == "deepvariant":
            reporter.info("Starting variant calling using DeepVariant...", color='c')
            vc_cmd = build_deepvariant_command(args)
        elif args.tool == "tvc":
            reporter.info("Starting variant calling using Torrent Variant Caller...", color='c')
            # Step 1: Prepare the command
            vc_cmd = [
                "tvc",
                "--num-threads", str(os.cpu_count()),
                "--input-bam", f"{args.bam}",
                "--target-file", f"{args.bed}",
                "--output-dir", f"{args.outdir}",
                "-r", f"{args.reference}"
            ]
            if hasattr(args, "mindepth"):
                min_depth = int(args.mindepth)
                if min_depth < 2:
                    raise ValueError("TVC minimum coverage must be at least 2")
                reference_index = vc_cmd.index("-r")
                vc_cmd[reference_index:reference_index] = [
                    "--gen-min-coverage", str(min_depth),
                    "--snp-min-coverage", str(min_depth),
                    "--mnp-min-coverage", str(min_depth),
                    "--indel-min-coverage", str(min_depth),
                    "--hotspot-min-coverage", str(min_depth),
                ]
        else:
            safe_report(reporter, "error", "Invalid variant calling tool. Only DeepVariant (--tool deepvariant) and Torrent Variant Caller (--tool tvc) are supported. Exiting.")
            raise Exception("Invalid variant calling tool. Only DeepVariant (--tool deepvariant) and Torrent Variant Caller (--tool tvc) are supported. Exiting.")

        reporter.command(vc_cmd)

        # Step 2: Run the command
        pending = ["|", "/", "-", "\\"]
        proc = subprocess.Popen(vc_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        while True:
            ready_fds, _, _ = select.select([proc.stdout, proc.stderr], [], [], 0.1)  # Faster spinner update
            
            if proc.stdout in ready_fds:
                output = proc.stdout.readline()
                if output:
                    reporter.stdout(output.strip())
            
            if proc.stderr in ready_fds:
                error = proc.stderr.readline()
                if error:
                    reporter.stderr(error.strip())
            
            # Check if process has terminated
            if proc.poll() is not None:
                break
            
            reporter.status(f"Calling variants... {pending[int(time.time() * 10) % 4]}")

        # Ensure all remaining output is read
        for output in proc.stdout:
            reporter.stdout(output.strip())
        for error in proc.stderr:
            reporter.stderr(error.strip())

        if proc.returncode != 0:
            safe_report(reporter, "error", "Variant caller execution failed. Exiting.")
            raise Exception("Variant caller execution failed. Exiting.")
        
        if args.tool == "tvc":            
            reporter.info("Merging VCF files...", color='c')
            merged_vcf_path = os.path.join(args.outdir, "mutated.vcf.gz")    
            temp_merged_vcf_path = tempfile.mktemp(prefix=".temp_reflow_", suffix=".vcf", dir=args.outdir)
            temp_compressed_snv_vcf_path = tempfile.mktemp(prefix=".temp_reflow_", suffix=".vcf.gz", dir=args.outdir)
            temp_compressed_indel_vcf_path = tempfile.mktemp(prefix=".temp_reflow_", suffix=".vcf.gz", dir=args.outdir)

            reporter.info("Compressing the sorted VCF file...", color='c')
            with open(temp_compressed_snv_vcf_path, 'wb') as output_handle:
                subprocess.run(["bgzip", "-c", os.path.join(args.outdir, "small_variants.vcf")], stdout=output_handle, check=True)
            reporter.info("Indexing the compressed VCF file...", color='c')
            subprocess.run(["bcftools", "index", temp_compressed_snv_vcf_path], check=True)

            reporter.info("Compressing the sorted VCF file...", color='c')
            with open(temp_compressed_indel_vcf_path, 'wb') as output_handle:
                subprocess.run(["bgzip", "-c", os.path.join(args.outdir, "indel_assembly.vcf")], stdout=output_handle, check=True)
            reporter.info("Indexing the compressed VCF file...", color='c')
            subprocess.run(["bcftools", "index", temp_compressed_indel_vcf_path], check=True)

            subprocess.run(["bcftools", "concat", "-a", "-O", "v", "-o", temp_merged_vcf_path,
                            temp_compressed_snv_vcf_path, 
                            temp_compressed_indel_vcf_path], check=True)     
            
            reporter.info("Sorting the merged VCF file...", color='c')
            subprocess.run(["bcftools", "sort", "-O", "z", "-o", merged_vcf_path, temp_merged_vcf_path], check=True)

            reporter.info("Indexing the compressed VCF file...", color='c')
            subprocess.run(["bcftools", "index", merged_vcf_path], check=True)
            
            #os.remove(temp_compressed_indel_vcf_path)
            #os.remove(temp_compressed_snv_vcf_path)
            #os.remove(temp_merged_vcf_path)
            #[os.remove(f) for f in glob.glob('.temp*')]            

        reporter.info(f"Variant calling completed. Results saved to {args.outdir}/mutated.vcf.gz", color='g')

    finally:
        # Clean up
        safe_report(reporter, "info", "Cleaning up...", color='d')
        pattern = os.path.join(args.outdir, '*.temp*')
        files_to_delete = glob.glob(pattern)
        for file_path in files_to_delete:
            if os.path.isfile(file_path):
                try:
                    os.remove(file_path)
                except OSError as exc:
                    safe_report(reporter, "warning", f"Cleanup could not remove {file_path}: {exc}")
        #if os.path.exists(temp_dir):
        #    os.rmdir(temp_dir)
