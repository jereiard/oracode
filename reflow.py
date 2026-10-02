import argparse
import sys

from extract_roi import extract_roi
from bespoke import bespoke2
from oracode import oracode
from educe_variant import educe_variant
from igvbatch import igvbatch
from remede import remede
from proficiency import PROFILE_NAMES, proficiency
from reporting import build_reporter, safe_report


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description="Process gene symbols and refseq database.")
    parser.add_argument(
        "--ui",
        choices=["textual", "plain", "curses"],
        default="textual",
        help="Runtime UI reporter (default: textual; falls back to plain if unavailable)",
    )
    parser.add_argument("--log-level", default="INFO", help="Runtime log level (default: INFO)")
    parser.add_argument(
        "--ui-hold",
        choices=["never", "on-error", "always"],
        default="never",
        help="Keep Textual UI open after completion: never, on-error, or always (default: never)",
    )
    subparsers = parser.add_subparsers(dest="command")

    extract_parser = subparsers.add_parser("extractROI", help="Extract regions of interest based on gene symbols")
    extract_parser.add_argument("--genelist", required=True, help="Path to the gene symbol file")
    extract_parser.add_argument("--refseq-db", required=True, help="Path to a specific RefSeq database file")
    extract_parser.add_argument("--outdir", required=True, help="Output directory")
    extract_parser.add_argument("--logfile", help="Filename to save reporter output as a logfile (default: subcommand-date-time.log in the output directory)")

    bespoke_parser = subparsers.add_parser("bespoke", help="Intersect BAM with BED regions")
    bespoke_parser.add_argument("--bam", required=True, help="Path to the BAM file")
    bespoke_parser.add_argument("--bed", required=True, nargs='+', help="Paths to the BED files")
    bespoke_parser.add_argument("--outdir", required=True, help="Output directory")
    bespoke_parser.add_argument("--logfile", help="Filename to save reporter output as a logfile (default: subcommand-date-time.log in the output directory)")

    oracode_parser = subparsers.add_parser("oracode", help="Run bamsurgeon to add SNVs and indels to BAM")
    oracode_parser.add_argument("--bam", required=True, help="Path to the BAM file")
    oracode_parser.add_argument("--snv", help="Path to the SNV file")
    oracode_parser.add_argument("--reference", required=True, help="Path to the reference genome file")
    oracode_parser.add_argument("--indel", help="Path to the indel file")
    oracode_parser.add_argument("--outdir", required=True, help="Output directory")
    oracode_parser.add_argument("--aligner", required=True, help="Aligner to use (e.g., bwa, bowtie2)")
    oracode_parser.add_argument("--logfile", help="Filename to save reporter output as a logfile (default: subcommand-date-time.log in the output directory)")

    educe_variant_parser = subparsers.add_parser("educeVariant", help="Call SNV/INDEL mutation using DeepVariant")
    educe_variant_parser.add_argument("--bam", required=True, help="Path to the BAM file")
    educe_variant_parser.add_argument("--bed", required=True, help="Path to the BED file")
    educe_variant_parser.add_argument("--reference", required=True, help="Path to the reference genome file")
    educe_variant_parser.add_argument("--tool", required=True, help="Name of the variant calling tool (e.g., deepvariant, tvc)")
    educe_variant_parser.add_argument("--outdir", required=True, help="Output directory")
    educe_variant_parser.add_argument("--logfile", help="Filename to save reporter output as a logfile (default: subcommand-date-time.log in the output directory)")

    igvbatch_parser = subparsers.add_parser("igvbatch", help="Create IGV batch script for common genomic positions")
    igvbatch_parser.add_argument("--vcf", required=True, help="Path to the VCF file")
    igvbatch_parser.add_argument("--mut", required=True, help="Path to the mut file")
    igvbatch_parser.add_argument("--outdir", required=True, help="Output directory")
    igvbatch_parser.add_argument("--logfile", help="Filename to save reporter output as a logfile (default: subcommand-date-time.log in the output directory)")

    remede_parser = subparsers.add_parser("remede", help="Fetch and modify reads in BAM file")
    remede_parser.add_argument("--bam", required=True, help="Path to the BAM file")
    remede_parser.add_argument("--muts", required=True, nargs='+', help="Paths to the mutation files")
    remede_parser.add_argument("--outdir", required=True, help="Output directory")
    remede_parser.add_argument("--floworder", default="TCAG", help="Flow order (default: TCAG)")
    remede_parser.add_argument("--barcode", default="CAGATCCATCGAT", help="Barcode (default: CAGATCCATCGAT)")
    remede_parser.add_argument("--library-key", default="TACGTACGTCTGAGCATCGATCGATGTACAGC", help="Library key (default: TACGTACGTCTGAGCATCGATCGATGTACAGC)")
    remede_parser.add_argument("--except-reads-by-name", nargs='*', help="Read names to exclude (optional)")
    remede_parser.add_argument("--logfile", help="Filename to save reporter output as a logfile (default: subcommand-date-time.log in the output directory)")

    proficiency_parser = subparsers.add_parser(
        "proficiency",
        help="Run the four-profile proficiency workflow from one variants CSV",
    )
    proficiency_parser.add_argument(
        "variants",
        nargs="?",
        default="materials/variants.csv",
        help="Variants CSV (default: materials/variants.csv)",
    )
    proficiency_parser.add_argument(
        "--data-root",
        help="Directory containing materials/, sources/, and references/ (auto-detected by default)",
    )
    proficiency_parser.add_argument(
        "--outdir",
        help="Output directory (default: results/<variants CSV stem>)",
    )
    proficiency_parser.add_argument(
        "--profiles",
        nargs="+",
        choices=["all", *PROFILE_NAMES],
        default=["all"],
        help="Profiles to run (default: all four profiles)",
    )
    proficiency_parser.add_argument(
        "--vaf",
        type=float,
        default=0.5,
        help="Default variant allele fraction when CSV has no vaf column (default: 0.5)",
    )
    proficiency_parser.add_argument(
        "--mindepth",
        type=int,
        default=2,
        help=(
            "Minimum BAMSurgeon site depth and TVC total coverage threshold "
            "(minimum/default: 2)"
        ),
    )
    proficiency_parser.add_argument(
        "--min-mutated-reads",
        type=int,
        default=1,
        help=(
            "Minimum BAMSurgeon mutated reads and DeepVariant candidate alt-read count "
            "(default: 1)"
        ),
    )
    proficiency_parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Validate inputs and create target/mutation files without modifying BAMs",
    )
    proficiency_parser.add_argument("--resume", action="store_true", help="Reuse completed BAM/VCF outputs")
    proficiency_parser.add_argument("--fail-fast", action="store_true", help="Stop after the first failed profile")
    proficiency_parser.add_argument("--skip-igv", action="store_true", help="Generate IGV batch files but do not launch IGV")
    proficiency_parser.add_argument("--igv-command", help="IGV or igv.sh executable path")
    proficiency_parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate BAMs after the standardized mutated.bam is created",
    )
    proficiency_parser.add_argument(
        "--logfile",
        help="Pipeline log filename inside the output directory (default: pipeline.log)",
    )
    return parser.parse_args(argv)


def show_banner():
    print("================================================")
    print("             ORACODE in GENECEPTION             ")
    print("                version 0.1.0                   ")
    print("         Copyright 2024, by Joowon Jang         ")
    print("================================================")
    print("""
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
""")


def show_help():
    print("Usage: python reflow.py [--ui textual|plain|curses] [subcommand] [options]\r")
    print("Subcommands:\r")
    print("  extractROI --genelist <path> --refseq-db <path> --outdir <path>  Extract regions of interest based on gene symbols\r")
    print("  bespoke    --bam <path> --bed <path...> --outdir <path>          Intersect BAM with BED regions\r")
    print("  oracode    --bam <path> --reference <path> --aligner <name>      Add SNVs/INDELs using bamsurgeon\r")
    print("  educeVariant --bam <path> --bed <path> --reference <path>        Call variants\r")
    print("  igvbatch   --vcf <path> --mut <path> --outdir <path>             Create IGV batch scripts\r")
    print("  remede     --bam <path> --muts <path...> --outdir <path>         Repair IonTorrent read metadata\r")
    print("  proficiency [materials/variants.csv] [--data-root <path>]        Run hg19/hg38 x Illumina/IonTorrent workflow\r")


def dispatch(args, reporter):
    if args.command == "extractROI":
        return extract_roi(args, reporter)
    if args.command == "bespoke":
        return bespoke2(args, reporter)
    if args.command == "oracode":
        return oracode(args, reporter)
    if args.command == "educeVariant":
        return educe_variant(args, reporter)
    if args.command == "igvbatch":
        return igvbatch(args, reporter)
    if args.command == "remede":
        return remede(args, reporter)
    if args.command == "proficiency":
        return proficiency(args, reporter)
    show_help()
    return None


def _report_exception(reporter, exc):
    safe_report(reporter, "error", str(exc), exc_info=True)


def _close_reporter(reporter):
    try:
        reporter.close()
    except Exception as close_exc:
        if sys.exc_info()[0] is not None:
            print(f"WARN: Reporter cleanup failed: {close_exc}", file=sys.stderr)
        else:
            raise


def main(argv=None):
    show_banner()
    args = parse_arguments(argv)
    if args.command == "help" or args.command is None:
        show_help()
        return

    if args.ui == "textual":
        try:
            from textual_ui import TextualDashboardUnavailable, run_textual_dashboard

            return run_textual_dashboard(args, dispatch, hold=args.ui_hold, log_level=args.log_level)
        except TextualDashboardUnavailable as exc:
            reporter = build_reporter("plain", log_level=args.log_level)
            reporter.warning(f"Textual UI unavailable; falling back to plain reporter: {exc}")
        else:  # pragma: no cover - return above exits on success
            return
    else:
        reporter = build_reporter(args.ui, log_level=args.log_level)

    try:
        return dispatch(args, reporter)
    except Exception as exc:
        _report_exception(reporter, exc)
        raise
    finally:
        _close_reporter(reporter)


if __name__ == "__main__":
    main()
