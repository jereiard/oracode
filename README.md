# Oracode

> 한국어: [README.ko.md](README.ko.md)

Oracode generates in silico proficiency-testing materials for germline NGS.
It extracts target regions from a public alignment BAM, introduces defined
variants with BAMSurgeon, rewrites the Ion Torrent flow signal so the edited
reads are callable, runs variant calling (TVC for Ion Torrent, DeepVariant for
Illumina), and produces per-variant IGV snapshots and reports.

The defining step is Ion Torrent flow-signal regeneration: a BAM whose query
sequence is edited without updating the ZM/FO/KS tags is not interpretable by
Torrent Variant Caller, so introduced variants are not called. Oracode
recalculates the ZM array from the flow order and a deterministic per-read seed.

## Contents

| Path | Purpose |
|---|---|
| `reflow.py` | CLI entry point (`proficiency`, `extractROI`, `bespoke`, `oracode`, `remede`, `educeVariant`, `igvbatch`) |
| `proficiency.py` | End-to-end four-profile workflow (hg19/hg38 x Illumina/IonTorrent) |
| `oracode.py` | Sequence-level variant introduction (BAMSurgeon driver) |
| `remede.py` | Ion Torrent flow-signal (ZM) regeneration |
| `bespoke.py` | Region-of-interest extraction from the source BAM |
| `educe_variant.py` | Variant calling (TVC / DeepVariant) and VCF merging |
| `igvbatch.py`, `igv-headless.sh` | Headless IGV snapshot generation |
| `extract_roi.py` | Target-region definition from a gene list |
| `reporting.py`, `textual_ui.py`, `curses_utils.py`, `logfile_manager.py` | Console reporting and UI |
| `Dockerfile` | Image build (base + final image) |
| `conda-explicit.txt` | Pinned conda package list for the base |
| `patches/` | BAMSurgeon patches |
| `oracode`, `console`, `dist.sh` | Host launcher and image build scripts |

## Install

See [INSTALL.md](INSTALL.md). Build the image from this repository. The only
external input is the Ion Torrent toolchain (`tmap`/`tvc`); supply it via
`TS_IMAGE`, built from the [Torrent Suite sources](https://github.com/iontorrent/TS).

## Quick start

```bash
./oracode --ui plain proficiency materials/variants.csv --prepare-only
./oracode proficiency materials/variants.csv
```

## Inputs

Two inputs are required: a variants CSV and a data root. The data root is
auto-detected as the directory containing `materials/`, `sources/`, and
`references/`; use `--data-root` to set it explicitly.

### Variants CSV

A header row with these columns (a bundled example is `materials/variants.csv`):

| Column | Example | Meaning |
|---|---|---|
| `contig` | `chr17` | Chromosome |
| `pos_hg19` | `7123443` | 1-based position on hg19 |
| `pos_hg38` | `7220124` | 1-based position on hg38 |
| `ref` | `C` | Reference allele |
| `alt` | `A` | Alternate allele |
| `gene` | `ACADVL` | Gene symbol |
| `type` | `SNV` | `SNV`, `INS`, or `DEL` |

Two columns are optional:

| Column | Example | Meaning |
|---|---|---|
| `vaf` | `0.25` | Variant allele fraction to introduce, `0 < vaf <= 1`. A blank cell or an absent column falls back to `--vaf` (default `0.5`). |
| `id` | `CVSET1-042` | Variant identifier recorded in `variants.results.csv`. When the column is absent or a cell is blank, it is derived from the locus as `gene:contig:pos_hg19:ref>alt`, and a row that lists the same locus again gets `#2`, `#3`, ... appended (e.g. `ACADVL:chr17:7123443:C>A#2`). When specifying it manually, the identifier must be unique. |

`vaf` is written into the per-variant BAMSurgeon input files, so each variant
is spiked in at that fraction rather than at one global rate; out-of-range or
non-numeric values abort the run.

### Data root layout

```
materials/
  variants.csv                              # default input (any path works)
  gencode.v49lift37.annotation.sorted.gtf.gz   # hg19 targets (proficiency)
  gencode.v49.annotation.sorted.gtf.gz         # hg38 targets (proficiency)
  <gene list>                               # extractROI only (optional)
  <refseq>.db                               # extractROI only (optional)
references/
  hg19.chrs.fa  (+ .fai, .bwt/.pac/.sa/.amb/.ann, .tmap.*)   # bwa + tmap indexes
  hg38.chrs.fa  (+ .fai, .bwt/.pac/.sa/.amb/.ann, .tmap.*)
sources/
  illumina/hg19/HG001/HG001.bam  (+ .bai)
  illumina/hg19/HG001/HG001.bed
  illumina/hg38/HG001/HG001.bam  (+ .bai)
  illumina/hg38/HG001/HG001.bed
  iontorrent/hg19/HG001/HG001.bam (+ .bai)
  iontorrent/hg19/HG001/HG001.bed
  iontorrent/hg38/HG001/HG001.bam (+ .bai)
  iontorrent/hg38/HG001/HG001.bed
```

Only the GENCODE GTF (materials/) and the source BAM/BED (sources/) are read
by `proficiency`. The gene list and RefSeq database are **not used by
`proficiency`**; they belong to the separate `extractROI` subcommand.

Notes:

- The reference FASTA must carry both the `bwa` index (Illumina `--aligner mem`)
  and the `tmap` index (Ion Torrent `--aligner tmap`).
- The `sources/<platform>/<genome>/HG001/HG001.bed` is the panel/ROI definition
  for that platform.
- Ion Torrent source BAMs must be realigned Torrent Suite BAMs carrying the
  `ZM`/`FO`/`KS` tags; Illumina source BAMs are standard paired-end BAMs.
- Only the profiles you pass with `--profiles` need their files present, e.g.
  `--profiles iontorrent-hg38` requires only that one source directory.

## License

Apache License 2.0. See [LICENSE](LICENSE).
