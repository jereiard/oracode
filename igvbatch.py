import os
import time
import pysam
from reporting import safe_report


def extract_common_positions(reporter, vcf_path, mut_path):
    common_positions = []
    vcf = pysam.VariantFile(vcf_path)
    with open(mut_path, 'r') as mut_file:
        for line in mut_file:
            contig, start, end = line.strip().split()[:3]
            reporter.info(f"Checking {contig}:{start}-{end}...")
            start, end = int(start), int(end)
            for record in vcf.fetch(contig, start - 2, end + 2):
                if start - 2 <= record.pos <= end + 2:
                    common_positions.append((contig, start, end, record))
    return common_positions


def create_igv_batch_script(common_positions, mut_path, output_dir):
    mut_filename = os.path.splitext(os.path.basename(mut_path))[0]
    script_path = os.path.join(output_dir, f"igvbatch.{mut_filename}.script")
    with open(script_path, 'w') as script_file:
        script_file.write("preference SAM.SHOW_CENTER_LINE true\n")
        script_file.writelines("preference SAM.SHOW_SOFT_CLIPPED true\n")
        script_file.writelines("preference FLANKING_REGION\t100\n")
        script_file.write(f"snapshotDirectory screenshots_{mut_filename}\n")
        screenshots_dir = os.path.join(output_dir, f"screenshots_{mut_filename}")
        os.mkdir(screenshots_dir) if not os.path.exists(screenshots_dir) else None
        script_file.write("viewaspairs\n")
        script_file.write("collapse\n")
        script_file.write("colorBy READ_STRAND\n")
        for contig, start, end, record in common_positions:
            ref = record.ref
            alt = record.alts[0]
            genotype = record.samples[0]['GT']
            genotype_str = "hetero" if genotype == (0, 1) else "homo"
            filename = f"{contig}-{start}-{end}-{ref}-{alt}-{genotype_str}.png"
            script_file.write(f"goto {contig}:{int(start)-80}-{int(end)+80}\n")
            script_file.write(f"sort BASE {contig}:{(int(start)+int(end))//2}\n")
            script_file.write(f"snapshot {filename}\n")
    return script_path


def save_stats(mut_path, common_positions, output_dir):
    mut_filename = os.path.splitext(os.path.basename(mut_path))[0]
    stats_path = os.path.join(output_dir, f"stats.{mut_filename}.tsv")
    total_mutations = sum(1 for _ in open(mut_path))
    common_mutations = len(common_positions)
    mutagenesis_rate = (common_mutations / total_mutations) * 100 if total_mutations else 0.0
    with open(stats_path, 'w') as stats_file:
        stats_file.write(f"Mutagenesis rate: {mutagenesis_rate:.2f}%\n")
        with open(mut_path, 'r') as mut_file:
            for line in mut_file:
                contig, start, end = line.strip().split()[:3]
                is_common = any(contig == pos[0] and int(start) == pos[1] and int(end) == pos[2] for pos in common_positions)
                stats_file.write(f"{line.strip()}\t{'common' if is_common else 'not common'}\n")
    return stats_path


def igvbatch(args, reporter):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if args.logfile:
        logfile_path = os.path.join(args.outdir, args.logfile)
    else:
        current_time = time.strftime("%Y%m%d-%H%M%S")
        logfile_path = os.path.join(args.outdir, f"{args.command}-{current_time}.log")
    reporter.set_logfile(logfile_path)

    try:
        reporter.info("Extracting common genomic positions from VCF and mut files...")
        common_positions = extract_common_positions(reporter, args.vcf, args.mut)
        reporter.info(f"Found {len(common_positions)} common positions.")

        reporter.info("Creating IGV batch script...")
        script_path = create_igv_batch_script(common_positions, args.mut, args.outdir)
        reporter.info(f"IGV batch script saved to {script_path}.")

        reporter.info("Saving stats file...")
        stats_path = save_stats(args.mut, common_positions, args.outdir)
        reporter.info(f"Stats file saved to {stats_path}.")

        reporter.info("IGV batch processing completed.", color='g')
    except Exception as e:
        safe_report(reporter, "error", f"Error processing IGV batch: {e}")
        raise
