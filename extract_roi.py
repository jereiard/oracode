import os
import time

indent = 4
half_indent = 2


def read_gene_symbols(gene_file):
    with open(gene_file, 'r') as f:
        return [line.strip() for line in f]


def process_refseq_db(gene_symbols, outdir, reporter, genelist_basename, refseq_db):
    for refseq_file in [refseq_db]:
        total_genes = len(gene_symbols)
        processed_genes = 0
        progress = 0.0
        with open(refseq_file, 'r') as f:
            reporter.info(f"Reading RefSeq DB: {refseq_file}", half_indent)
            lines = f.readlines()
            reporter.info(f"Lines Read: {len(lines)}", indent)
        reporter.info(f"Extracting ROI regions: {gene_symbols}", indent)

        output_lines = []
        for gene in gene_symbols:
            len_output_lines_before = len(output_lines)
            for line in lines:
                columns = line.strip().split('\t')
                if columns[4] == gene:
                    output_lines.append(f"{columns[1]}\t{columns[2]}\t{columns[3]}\t{columns[4]}\n")
                    processed_genes += 1
                    progress = (processed_genes / total_genes) * 100 if total_genes else 100.0

            if len(output_lines) == len_output_lines_before:
                reporter.error(f"{gene} was not found in the RefSeq database. Skipping...", indent)
            reporter.progress(f"Progress: {processed_genes}/{total_genes} ({progress:.2f}%)")

        reporter.info(f"Completed: {processed_genes} of {total_genes} genes included.", indent, color='g')

        if output_lines:
            output_file = os.path.join(outdir, f"roi.{genelist_basename}.bed")
            with open(output_file, 'w') as f:
                f.writelines(output_lines)


def extract_roi(args, reporter):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if args.logfile:
        logfile_path = os.path.join(args.outdir, args.logfile)
    else:
        current_time = time.strftime("%Y%m%d-%H%M%S")
        logfile_path = os.path.join(args.outdir, f"{args.command}-{current_time}.log")
    reporter.set_logfile(logfile_path)

    reporter.info("Extracting regions of interest based on gene symbols")
    genelist_basename = os.path.splitext(os.path.basename(args.genelist))[0]

    reporter.info(f"Reading gene list: {args.genelist}", half_indent)
    gene_symbols = read_gene_symbols(args.genelist)
    process_refseq_db(gene_symbols, args.outdir, reporter, genelist_basename, args.refseq_db)
