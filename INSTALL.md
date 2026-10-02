# Install

> 한국어: [INSTALL.ko.md](INSTALL.ko.md)

Oracode runs as a Docker image built from this repository.

## Requirements

- Linux host with Docker (the launcher mounts the current directory and the
  Docker socket).
- ~20 GB free disk for the image and working files.
- No GPU is required; the launcher passes `--gpus all` for the optional
  DeepVariant path and falls back gracefully when none is present.
- **Ion Torrent toolchain image (required).** Ion Torrent variant calling uses
  `tmap`, `tvc`, `tvcutils`, and `tvcassembly`, which are not redistributable.
  You must supply your own image built from the Torrent Suite sources
  (https://github.com/iontorrent/TS) following that repository's build
  instructions (`buildTools/BUILD.txt`; the Analysis module is built with
  `MODULES=Analysis ./buildTools/build.sh`). Point the build at it with
  `--build-arg TS_IMAGE=<your-tsuite-image>` (default
  `jereiard/tsuite:5.18.1`). Everything else is built from public sources.

## Build the image

Build the image from this repository. The Dockerfile is self-contained: the
base (Ubuntu 20.04 + Miniforge + samtools/bcftools + BAMSurgeon + Picard) and
the final image are both built here from public sources. Only the Ion Torrent
toolchain (`TS_IMAGE`) is an external input. Network access is required to
fetch the base packages, BAMSurgeon, Picard, and IGV; `dist.sh` forwards any
proxy environment variables.

```bash
# Build (forwards HTTP(S)_PROXY / NO_PROXY when set).
./dist.sh

# Run against the built tag.
ORACODE_IMAGE=oracode ./oracode --help
```

Build arguments:

| Argument | Default | Purpose |
|---|---|---|
| `TS_IMAGE` | `jereiard/tsuite:5.18.1` | Supplies `tmap`, `tvc`, `tvcutils`, `tvcassembly`, and TVC parameter sets. Required for Ion Torrent variant calling. |
| `BAMSURGEON_REPO` | `https://www.github.com/jereiard/bamsurgeon` | BAMSurgeon source. |
| `BAMSURGEON_COMMIT` | `50b14964…` | Pinned BAMSurgeon commit. |
| `PICARD_VERSION` | `1.131` | Picard tools release. |
| `IGV_VERSION` / `IGV_SHA256` | `2.19.8` / pinned | IGV download and checksum. |

Verify the build:

```bash
docker run --rm oracode --help
docker run --rm oracode remede --help
```

## Running

`./oracode` wraps `docker run` and forwards all arguments:

```bash
# Validate inputs and build target/BAMSurgeon files only.
./oracode --ui plain proficiency materials/variants.csv --prepare-only

# Run all four profiles (hg19/hg38 x Illumina/IonTorrent).
./oracode proficiency materials/variants.csv

# Point at a data root explicitly when it is not auto-detected.
./oracode proficiency materials/variants.csv --data-root /path/to/data-root
```

Open a shell inside the image:

```bash
./console
```

## Data layout

The data root must contain `materials/`, `sources/`, and `references/`.
The exact directory tree (variants CSV columns, reference FASTA indexes,
per-platform source BAM/BED paths) is specified under **Inputs** in
[README.md](README.md) / [README.ko.md](README.ko.md).

## Offline / air-gapped hosts

IGV performs a Google OAuth lookup at startup that stalls headless batch runs
on air-gapped hosts. The image ships an offline OAuth config and
`igv-headless`, and `run_igv` points IGV at it, so no network call is made.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `tvc: command not found` | The TS tools were not copied in. Check the `TS_IMAGE` build arg. |
| IGV stalls at startup | OAuth lookup; confirm `igv-headless` and the offline OAuth config are in the image. |
| Proxy build failure | Export `HTTP_PROXY`/`HTTPS_PROXY`/`NO_PROXY` and use `./dist.sh`. |
| Output BAM unreadable by TVC | Flow signal was not regenerated; ensure `remede` runs before variant calling. |
