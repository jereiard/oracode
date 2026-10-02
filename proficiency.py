"""End-to-end proficiency-test workflow driven by a single variants CSV.

The legacy commands remain available, but this module gives them a stable
orchestration layer and a deterministic output contract for the four default
HG001 profiles (hg19/hg38 x Illumina/IonTorrent).
"""

from __future__ import annotations

import csv
import gzip
import html
import json
import os
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import pysam

from reporting import safe_report


GENOMES = ("hg19", "hg38")
PLATFORMS = ("illumina", "iontorrent")
PROFILE_NAMES = tuple(f"{platform}-{genome}" for genome in GENOMES for platform in PLATFORMS)
REQUIRED_COLUMNS = {"contig", "pos_hg19", "pos_hg38", "ref", "alt", "gene", "type"}


class PipelineConfigurationError(ValueError):
    """Raised when inputs or default materials do not satisfy the workflow contract."""


class PipelineExecutionError(RuntimeError):
    """Raised after reports are written when one or more profiles failed."""


@dataclass(frozen=True)
class Variant:
    variant_id: str
    contig: str
    pos_hg19: int
    pos_hg38: int
    ref: str
    alt: str
    gene: str
    kind: str
    vaf: float

    def position(self, genome: str) -> int:
        if genome == "hg19":
            return self.pos_hg19
        if genome == "hg38":
            return self.pos_hg38
        raise PipelineConfigurationError(f"Unsupported genome: {genome}")

    def normalized_key(self, genome: str) -> Tuple[str, int, str, str]:
        pos, ref, alt = normalize_alleles(self.position(genome), self.ref, self.alt)
        return self.contig, pos, ref, alt


@dataclass(frozen=True)
class Profile:
    name: str
    genome: str
    platform: str
    bam: Path
    source_bed: Path
    reference: Path
    gtf: Path
    aligner: str
    caller: str


@dataclass(frozen=True)
class PreparedInputs:
    genome: str
    directory: Path
    snv_path: Path
    indel_path: Path
    target_bed: Path


@dataclass
class ProfileResult:
    profile: str
    genome: str
    platform: str
    status: str = "pending"
    bam: Optional[str] = None
    validation_vcf: Optional[str] = None
    igv_batch: Optional[str] = None
    screenshots_dir: Optional[str] = None
    screenshot_count: int = 0
    roi_header_normalized: bool = False
    igv_status: str = "not-run"
    detected_variant_ids: List[str] = field(default_factory=list)
    error: Optional[str] = None


def normalize_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Return a simple representation-independent variant key.

    Trimming shared sequence makes VCF-style anchored indels comparable with
    callers that emit an equivalent minimal representation.
    """

    ref = ref.upper()
    alt = alt.upper()
    while ref and alt and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return pos, ref or "-", alt or "-"


def load_variants(path: Path, default_vaf: float = 0.5) -> List[Variant]:
    if not path.is_file():
        raise PipelineConfigurationError(f"Variants CSV not found: {path}")
    if not 0 < default_vaf <= 1:
        raise PipelineConfigurationError("Default VAF must be greater than 0 and at most 1")

    variants: List[Variant] = []
    seen_ids: Set[str] = set()
    generated_id_counts: Dict[str, int] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise PipelineConfigurationError(
                f"Variants CSV is missing required columns: {', '.join(sorted(missing))}"
            )
        for row_number, row in enumerate(reader, start=2):
            try:
                contig = row["contig"].strip()
                pos_hg19 = int(row["pos_hg19"])
                pos_hg38 = int(row["pos_hg38"])
                ref = row["ref"].strip().upper()
                alt = row["alt"].strip().upper()
                gene = row["gene"].strip()
                kind = row["type"].strip().upper()
                vaf = float((row.get("vaf") or default_vaf))
            except (TypeError, ValueError) as exc:
                raise PipelineConfigurationError(f"Invalid variants CSV row {row_number}: {exc}") from exc

            if not contig or not gene or pos_hg19 < 1 or pos_hg38 < 1:
                raise PipelineConfigurationError(f"Invalid coordinate/gene at CSV row {row_number}")
            if not ref or not alt or any(base not in "ACGTN" for base in ref + alt):
                raise PipelineConfigurationError(f"Invalid REF/ALT at CSV row {row_number}: {ref}>{alt}")
            if not 0 < vaf <= 1:
                raise PipelineConfigurationError(f"Invalid VAF at CSV row {row_number}: {vaf}")
            _validate_variant_shape(kind, ref, alt, row_number)

            supplied_id = (row.get("id") or "").strip()
            if supplied_id:
                variant_id = supplied_id
            else:
                base_id = f"{gene}:{contig}:{pos_hg19}:{ref}>{alt}"
                occurrence = generated_id_counts.get(base_id, 0) + 1
                generated_id_counts[base_id] = occurrence
                variant_id = base_id if occurrence == 1 else f"{base_id}#{occurrence}"
            if variant_id in seen_ids:
                raise PipelineConfigurationError(
                    f"Duplicate variant ID at CSV row {row_number}: {variant_id}"
                )
            seen_ids.add(variant_id)
            variants.append(
                Variant(variant_id, contig, pos_hg19, pos_hg38, ref, alt, gene, kind, vaf)
            )

    if not variants:
        raise PipelineConfigurationError("Variants CSV does not contain any variants")
    return variants


def _validate_variant_shape(kind: str, ref: str, alt: str, row_number: int) -> None:
    if kind == "SNV" and (len(ref) != 1 or len(alt) != 1):
        raise PipelineConfigurationError(f"SNV must have one-base REF/ALT at CSV row {row_number}")
    if kind == "INS" and not (len(alt) > len(ref) and alt.startswith(ref)):
        raise PipelineConfigurationError(
            f"INS must use anchored VCF alleles (ALT starts with REF) at CSV row {row_number}"
        )
    if kind == "DEL" and not (len(ref) > len(alt) and ref.startswith(alt)):
        raise PipelineConfigurationError(
            f"DEL must use anchored VCF alleles (REF starts with ALT) at CSV row {row_number}"
        )
    if kind not in {"SNV", "INS", "DEL"}:
        raise PipelineConfigurationError(f"Unsupported variant type at CSV row {row_number}: {kind}")


def find_data_root(explicit: Optional[str] = None) -> Path:
    candidates: List[Path] = []
    if explicit:
        candidates.append(Path(explicit))
    if os.environ.get("ORACODE_DATA_ROOT"):
        candidates.append(Path(os.environ["ORACODE_DATA_ROOT"]))
    candidates.extend([Path.cwd(), Path(__file__).resolve().parent])

    seen: Set[Path] = set()
    for candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if candidate in seen:
            continue
        seen.add(candidate)
        if all((candidate / name).is_dir() for name in ("materials", "sources", "references")):
            return candidate
    checked = ", ".join(str(path) for path in seen)
    raise PipelineConfigurationError(
        "Could not locate materials/, sources/, and references/. "
        f"Set --data-root or ORACODE_DATA_ROOT. Checked: {checked}"
    )


def resolve_profiles(data_root: Path, selected: Sequence[str]) -> List[Profile]:
    names = list(PROFILE_NAMES if not selected or "all" in selected else selected)
    invalid = sorted(set(names) - set(PROFILE_NAMES))
    if invalid:
        raise PipelineConfigurationError(f"Unsupported profiles: {', '.join(invalid)}")

    profiles: List[Profile] = []
    for name in PROFILE_NAMES:
        if name not in names:
            continue
        platform, genome = name.split("-", 1)
        sample_dir = data_root / "sources" / platform / genome / "HG001"
        reference = data_root / "references" / f"{genome}.chrs.fa"
        gtf_name = (
            "gencode.v49lift37.annotation.sorted.gtf.gz"
            if genome == "hg19"
            else "gencode.v49.annotation.sorted.gtf.gz"
        )
        profile = Profile(
            name=name,
            genome=genome,
            platform=platform,
            bam=sample_dir / "HG001.bam",
            source_bed=sample_dir / "HG001.bed",
            reference=reference,
            gtf=data_root / "materials" / gtf_name,
            aligner="mem" if platform == "illumina" else "tmap",
            caller="deepvariant" if platform == "illumina" else "tvc",
        )
        _validate_profile_files(profile)
        profiles.append(profile)
    return profiles


def _validate_profile_files(profile: Profile) -> None:
    required = [
        profile.bam,
        Path(f"{profile.bam}.bai"),
        profile.source_bed,
        profile.reference,
        Path(f"{profile.reference}.fai"),
        profile.gtf,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise PipelineConfigurationError(
            f"Profile {profile.name} is missing required files: {', '.join(missing)}"
        )


def validate_reference_alleles(variants: Sequence[Variant], genome: str, reference: Path) -> None:
    mismatches: List[str] = []
    fasta = pysam.FastaFile(str(reference))
    try:
        references = set(fasta.references)
        for variant in variants:
            pos = variant.position(genome)
            if variant.contig not in references:
                mismatches.append(f"{variant.variant_id}: contig {variant.contig} missing")
                continue
            observed = fasta.fetch(variant.contig, pos - 1, pos - 1 + len(variant.ref)).upper()
            if observed != variant.ref:
                mismatches.append(
                    f"{variant.variant_id}: {variant.contig}:{pos} expected {variant.ref}, reference has {observed}"
                )
    finally:
        fasta.close()
    if mismatches:
        preview = "; ".join(mismatches[:10])
        suffix = f" (+{len(mismatches) - 10} more)" if len(mismatches) > 10 else ""
        raise PipelineConfigurationError(f"Reference allele validation failed for {genome}: {preview}{suffix}")


def write_bamsurgeon_inputs(
    variants: Sequence[Variant], genome: str, output_dir: Path
) -> Tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    snv_path = output_dir / "bamsurgeon.snv.tsv"
    indel_path = output_dir / "bamsurgeon.indel.tsv"
    snv_lines: List[str] = []
    indel_lines: List[str] = []
    seen_snv_lines: Set[str] = set()
    seen_indel_lines: Set[str] = set()

    for variant in sorted(variants, key=lambda item: (item.contig, item.position(genome), item.variant_id)):
        pos = variant.position(genome)
        if variant.kind == "SNV":
            line = f"{variant.contig}\t{pos}\t{pos}\t{variant.vaf:g}\t{variant.alt}\n"
            if line not in seen_snv_lines:
                seen_snv_lines.add(line)
                snv_lines.append(line)
        elif variant.kind == "INS":
            inserted = variant.alt[len(variant.ref) :]
            line = f"{variant.contig}\t{pos}\t{pos + 1}\t{variant.vaf:g}\tINS\t{inserted}\n"
            if line not in seen_indel_lines:
                seen_indel_lines.add(line)
                indel_lines.append(line)
        else:
            deletion_length = len(variant.ref) - len(variant.alt)
            line = (
                f"{variant.contig}\t{pos}\t{pos + deletion_length}\t"
                f"{variant.vaf:g}\tDEL\t*\n"
            )
            if line not in seen_indel_lines:
                seen_indel_lines.add(line)
                indel_lines.append(line)

    snv_path.write_text("".join(snv_lines), encoding="utf-8")
    indel_path.write_text("".join(indel_lines), encoding="utf-8")
    return snv_path, indel_path


def _parse_gtf_attributes(attributes: str) -> Dict[str, str]:
    parsed: Dict[str, str] = {}
    for item in attributes.rstrip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def _merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def load_gene_intervals(gtf_path: Path, genes: Set[str]) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    found: Set[str] = set()
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    with opener(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9 or columns[2] != "gene":
                continue
            attrs = _parse_gtf_attributes(columns[8])
            gene = attrs.get("gene_name")
            if gene not in genes:
                continue
            # GTF is 1-based inclusive; BED is 0-based half-open.
            start, end = int(columns[3]) - 1, int(columns[4])
            intervals.setdefault(columns[0], []).append((start, end))
            found.add(gene)

    missing = sorted(genes - found)
    if missing:
        raise PipelineConfigurationError(
            f"Genes missing from {gtf_path.name}: {', '.join(missing)}"
        )
    return {contig: _merge_intervals(values) for contig, values in intervals.items()}


def intersect_source_bed(
    source_bed: Path, gene_intervals: Mapping[str, Sequence[Tuple[int, int]]]
) -> Dict[str, List[Tuple[int, int]]]:
    intersections: Dict[str, List[Tuple[int, int]]] = {}
    cursors: Dict[str, int] = {contig: 0 for contig in gene_intervals}
    with source_bed.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            columns = line.split()
            if len(columns) < 3 or columns[0] not in gene_intervals:
                continue
            contig, source_start, source_end = columns[0], int(columns[1]), int(columns[2])
            targets = gene_intervals[contig]
            cursor = cursors[contig]
            while cursor < len(targets) and targets[cursor][1] <= source_start:
                cursor += 1
            cursors[contig] = cursor
            index = cursor
            while index < len(targets) and targets[index][0] < source_end:
                start = max(source_start, targets[index][0])
                end = min(source_end, targets[index][1])
                if start < end:
                    intersections.setdefault(contig, []).append((start, end))
                index += 1
    return {contig: _merge_intervals(values) for contig, values in intersections.items()}


def write_target_bed(
    variants: Sequence[Variant], genome: str, gtf: Path, source_bed: Path, output_path: Path
) -> Path:
    genes = {variant.gene for variant in variants}
    # The proficiency CSV can intentionally contain intronic or panel-edge
    # positions that are outside the vendor capture BED.  Reducing to the
    # complete GENCODE gene intervals preserves any reads at those positions;
    # intersecting with the capture BED would silently make those mutations
    # impossible.  The source BED is still a required profile material and is
    # retained for provenance, but does not narrow the GENCODE target here.
    target_intervals = load_gene_intervals(gtf, genes)
    if not target_intervals:
        raise PipelineConfigurationError(
            f"No GENCODE intervals found for target genes in {genome}: {gtf}"
        )

    uncovered: List[str] = []
    for variant in variants:
        pos0 = variant.position(genome) - 1
        if not any(start <= pos0 < end for start, end in target_intervals.get(variant.contig, [])):
            uncovered.append(f"{variant.variant_id}@{variant.contig}:{variant.position(genome)}")
    if uncovered:
        preview = ", ".join(uncovered[:10])
        suffix = f" (+{len(uncovered) - 10} more)" if len(uncovered) > 10 else ""
        raise PipelineConfigurationError(
            f"Variants outside their selected GENCODE gene intervals for {genome}: {preview}{suffix}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for contig in sorted(target_intervals, key=_contig_sort_key):
            for start, end in target_intervals[contig]:
                handle.write(f"{contig}\t{start}\t{end}\n")
    return output_path


def _contig_sort_key(contig: str) -> Tuple[int, str]:
    token = contig.removeprefix("chr")
    if token.isdigit():
        return int(token), ""
    return 10_000, token


def prepare_inputs(
    variants: Sequence[Variant], profiles: Sequence[Profile], output_root: Path
) -> Dict[str, PreparedInputs]:
    prepared: Dict[str, PreparedInputs] = {}
    by_genome = {profile.genome: profile for profile in profiles}
    for genome in GENOMES:
        if genome not in by_genome:
            continue
        representative = by_genome[genome]
        directory = output_root / "inputs" / genome
        validate_reference_alleles(variants, genome, representative.reference)
        snv_path, indel_path = write_bamsurgeon_inputs(variants, genome, directory)

        # Illumina and IonTorrent source BEDs are expected to describe the same
        # capture design for a genome. Validate each selected profile and keep a
        # profile-specific BED if they ever diverge.
        profile_beds: Dict[str, Path] = {}
        for profile in [item for item in profiles if item.genome == genome]:
            target_path = directory / f"targets.{profile.platform}.bed"
            profile_beds[profile.platform] = write_target_bed(
                variants, genome, profile.gtf, profile.source_bed, target_path
            )
        canonical = profile_beds.get("illumina") or profile_beds["iontorrent"]
        prepared[genome] = PreparedInputs(genome, directory, snv_path, indel_path, canonical)
    return prepared


def detect_variants(vcf_path: Path, variants: Sequence[Variant], genome: str) -> Set[str]:
    expected: Dict[Tuple[str, int, str, str], str] = {
        variant.normalized_key(genome): variant.variant_id for variant in variants
    }
    detected: Set[str] = set()
    variant_file = pysam.VariantFile(str(vcf_path))
    try:
        for record in variant_file:
            if set(record.filter.keys()) != {"PASS"}:
                continue
            if not any(
                allele is not None and allele > 0
                for sample in record.samples.values()
                for allele in (sample.get("GT") or ())
            ):
                continue
            for alt in record.alts or ():
                pos, ref, normalized_alt = normalize_alleles(record.pos, record.ref, alt)
                variant_id = expected.get((record.contig, pos, ref, normalized_alt))
                if variant_id:
                    detected.add(variant_id)
    finally:
        variant_file.close()
    return detected


def create_igv_batch(
    profile: Profile,
    variants: Sequence[Variant],
    bam_path: Path,
    output_dir: Path,
) -> Tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    screenshots_dir = output_dir / "screenshots"
    screenshots_dir.mkdir(parents=True, exist_ok=True)
    batch_path = output_dir / "igv.batch"
    with batch_path.open("w", encoding="utf-8") as handle:
        # "genome" must be the FIRST command. IGV's batch runner loads the
        # user's default genome (hg38) over the network before executing the
        # first command unless that command is itself "genome"
        # (BatchRunner.runWithDefaultGenome). On an air-gapped host that lookup
        # times out and stalls the run. Emitting no "new" avoids the same eager
        # reload; loading our local FASTA directly keeps the batch offline.
        handle.write(f"genome {profile.reference}\n")
        handle.write("setSleepInterval 100\n")
        handle.write("preference SAM.SHOW_CENTER_LINE true\n")
        handle.write("preference SAM.SHOW_SOFT_CLIPPED true\n")
        handle.write("preference FLANKING_REGION 100\n")
        handle.write(f"snapshotDirectory {screenshots_dir}\n")
        handle.write(f"load {bam_path}\n")
        handle.write("collapse\n")
        handle.write("colorBy READ_STRAND\n")
        if profile.platform == "illumina":
            handle.write("viewaspairs\n")
        for variant in variants:
            pos = variant.position(profile.genome)
            start, end = max(1, pos - 80), pos + max(len(variant.ref), len(variant.alt)) + 80
            filename = igv_screenshot_filename(variant)
            handle.write(f"goto {variant.contig}:{start}-{end}\n")
            handle.write(f"sort BASE {variant.contig}:{pos}\n")
            handle.write(f"snapshot {filename}\n")
        handle.write("exit\n")
    return batch_path, screenshots_dir


def igv_screenshot_filename(variant: Variant) -> str:
    return f"{variant.variant_id.replace(':', '_').replace('>', '-')}.png"


def resolve_igv_offline_directory(executable: Optional[str]) -> Optional[Path]:
    """Return an IGV user directory that disables the Google OAuth lookup.

    IGV opens an HTTPS connection to PROVISIONING_URL_DEFAULT while building its
    menu bar, before it is running as a batch. On an air-gapped host that
    connection fails (proxy 407 or timeout) and IGV raises a modal dialog that
    stalls the headless run. A user directory containing a valid (though
    unusable) oauth-config.json makes ``fetchOauthConfigs`` register a Google
    provider, so ``getGoogleProvider`` never reaches the network.

    The directory is searched next to the IGV executable (``<igv>/oauth-config.json``)
    and at ``IGV_OFFLINE_DIR``; the first that already has an ``oauth-config.json``
    wins, otherwise ``None`` keeps the prior behaviour.
    """
    candidates: List[Path] = []
    if executable:
        candidates.append(Path(executable).resolve().parent / "offline")
    override = os.environ.get("IGV_OFFLINE_DIR")
    if override:
        candidates.append(Path(override))
    for candidate in candidates:
        if (candidate / "oauth-config.json").is_file():
            return candidate
    return None


def run_igv(batch_path: Path, command: Optional[str], reporter) -> str:
    snapshot_directory: Optional[Path] = None
    genome: Optional[str] = None
    expected: List[str] = []
    for line in batch_path.read_text(encoding="utf-8").splitlines():
        command_name, _, argument = line.partition(" ")
        if command_name == "genome":
            genome = argument
        elif command_name == "snapshotDirectory":
            snapshot_directory = Path(argument)
        elif command_name == "snapshot":
            expected.append(argument)

    executable = command or os.environ.get("IGV_EXECUTABLE")
    if executable:
        resolved = shutil.which(executable) or executable
    else:
        resolved = shutil.which("igv.sh") or shutil.which("igv")
    if not resolved:
        reporter.warning(
            f"IGV executable not found; batch script generated at {batch_path}. "
            "Set --igv-command or IGV_EXECUTABLE to create screenshots."
        )
        return "batch-only"
    cmd = [resolved]
    if genome:
        cmd.extend(["--genome", genome])
    offline_directory = resolve_igv_offline_directory(resolved)
    if offline_directory is not None:
        cmd.extend(["--igvDirectory", str(offline_directory)])
    cmd.extend(["-b", str(batch_path)])
    if not os.environ.get("DISPLAY"):
        headless = shutil.which("igv-headless")
        if headless:
            cmd = [headless, *cmd]
        elif shutil.which("xvfb-run"):
            cmd = ["xvfb-run", "-a", *cmd]
    reporter.command(cmd)
    subprocess.run(cmd, check=True)
    if snapshot_directory is None:
        raise PipelineExecutionError(f"IGV batch has no snapshotDirectory: {batch_path}")
    missing = [
        name
        for name in expected
        if not (snapshot_directory / name).is_file()
        or (snapshot_directory / name).stat().st_size == 0
    ]
    if missing:
        sample = ", ".join(missing[:3])
        raise PipelineExecutionError(
            f"IGV did not create {len(missing)}/{len(expected)} screenshots: {sample}"
        )
    return "completed"


def _link_final_bam(source: Path, destination: Path) -> Path:
    if source.resolve() == destination.resolve():
        return destination
    for path in (destination, Path(f"{destination}.bai")):
        if path.exists() or path.is_symlink():
            path.unlink()
    os.link(source, destination)
    source_index = Path(f"{source}.bai")
    if source_index.is_file():
        os.link(source_index, Path(f"{destination}.bai"))
    else:
        pysam.index(str(destination))
    return destination


def bam_has_reads(path: Path) -> bool:
    if not path.is_file():
        return False
    with pysam.AlignmentFile(str(path), "rb") as alignment:
        return next(alignment.fetch(until_eof=True), None) is not None


def normalize_roi_header_for_reference(
    bam_path: Path, reference_path: Path, reporter
) -> bool:
    """Normalize a derived ROI BAM when its SQ dictionary differs from FASTA.

    BamSurgeon requires the original and realigned donor BAMs to use identical
    numeric reference IDs. Some hg38 source BAMs include alternate contigs and
    use a different SQ order from the supplied primary-contig FASTA. The source
    BAM is never modified; only the already reduced ROI copy is rewritten.
    """

    with pysam.FastaFile(str(reference_path)) as fasta:
        reference_names = tuple(fasta.references)
        reference_lengths = tuple(fasta.lengths)
    reference_id = {name: index for index, name in enumerate(reference_names)}

    with pysam.AlignmentFile(str(bam_path), "rb") as source:
        observed = tuple(zip(source.references, source.lengths))
        expected = tuple(zip(reference_names, reference_lengths))
        if observed == expected:
            return False

        reporter.warning(
            f"Normalizing derived ROI BAM SQ order to match reference: {bam_path}"
        )
        header = source.header.to_dict()
        header["SQ"] = [
            {"SN": name, "LN": length}
            for name, length in zip(reference_names, reference_lengths)
        ]
        header.setdefault("HD", {})["SO"] = "unsorted"
        unsorted_path = bam_path.with_suffix(".canonical.unsorted.bam")
        sorted_path = bam_path.with_suffix(".canonical.sorted.bam")
        with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as destination:
            for read in source.fetch(until_eof=True):
                if read.is_unmapped or read.reference_id < 0:
                    continue
                source_reference = source.get_reference_name(read.reference_id)
                if source_reference not in reference_id:
                    continue
                mate_reference = (
                    source.get_reference_name(read.next_reference_id)
                    if read.next_reference_id >= 0
                    else None
                )
                read.reference_id = reference_id[source_reference]
                if mate_reference in reference_id:
                    read.next_reference_id = reference_id[mate_reference]
                elif read.next_reference_id >= 0:
                    read.next_reference_id = -1
                    read.next_reference_start = -1
                    read.mate_is_unmapped = True
                    read.template_length = 0
                destination.write(read)

    try:
        pysam.sort(
            "-@",
            str(max(1, min(os.cpu_count() or 1, 16))),
            "-o",
            str(sorted_path),
            str(unsorted_path),
        )
        sorted_path.replace(bam_path)
        index_path = Path(f"{bam_path}.bai")
        if index_path.exists():
            index_path.unlink()
        pysam.index(str(bam_path))
    finally:
        for temporary in (unsorted_path, sorted_path):
            if temporary.exists():
                temporary.unlink()
    if not bam_has_reads(bam_path):
        raise PipelineExecutionError(
            f"ROI BAM contains no primary-reference reads after header normalization: {bam_path}"
        )
    return True


def _remove_intermediate_bams(directory: Path, final_bam: Path) -> None:
    for path in directory.glob("*.bam*"):
        if path == final_bam or path == Path(f"{final_bam}.bai"):
            continue
        if path.is_file() or path.is_symlink():
            path.unlink()


class ProficiencyPipeline:
    def __init__(self, args, reporter):
        self.args = args
        self.reporter = reporter
        self.data_root = find_data_root(args.data_root)
        variants_path = Path(args.variants).expanduser()
        if not variants_path.is_absolute() and not variants_path.is_file():
            variants_path = self.data_root / variants_path
        self.variants_path = variants_path.resolve()
        self.output_root = (
            Path(args.outdir).expanduser().resolve()
            if args.outdir
            else (Path.cwd() / "results" / self.variants_path.stem).resolve()
        )

    def run(self) -> Mapping[str, object]:
        self.output_root.mkdir(parents=True, exist_ok=True)
        if self.args.mindepth < 2:
            raise PipelineConfigurationError("--mindepth must be at least 2")
        if self.args.min_mutated_reads < 1:
            raise PipelineConfigurationError("--min-mutated-reads must be at least 1")
        logfile = self.output_root / (self.args.logfile or "pipeline.log")
        self.reporter.set_logfile(str(logfile))
        self.reporter.info(f"Proficiency workflow input: {self.variants_path}")
        self.reporter.info(f"Data root: {self.data_root}")
        self.reporter.info(f"Output root: {self.output_root}")

        variants = load_variants(self.variants_path, self.args.vaf)
        profiles = resolve_profiles(self.data_root, self.args.profiles)
        prepared = prepare_inputs(variants, profiles, self.output_root)
        self.reporter.info(
            f"Prepared {len(variants)} variants for {len(profiles)} profiles", color="g"
        )

        results: List[ProfileResult] = []
        if self.args.prepare_only:
            results = [
                ProfileResult(profile.name, profile.genome, profile.platform, status="prepared")
                for profile in profiles
            ]
        else:
            for index, profile in enumerate(profiles, start=1):
                self.reporter.progress(f"Profile {index}/{len(profiles)}: {profile.name}")
                try:
                    result = self._run_profile(profile, prepared[profile.genome], variants)
                except Exception as exc:
                    safe_report(self.reporter, "error", f"Profile {profile.name} failed: {exc}")
                    result = ProfileResult(
                        profile.name,
                        profile.genome,
                        profile.platform,
                        status="failed",
                        error=str(exc),
                    )
                    results.append(result)
                    if self.args.fail_fast:
                        break
                    continue
                results.append(result)

        artifacts = write_reports(self.output_root, self.variants_path, variants, profiles, results)
        manifest = {
            "input": str(self.variants_path),
            "data_root": str(self.data_root),
            "output_root": str(self.output_root),
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "prepare_only": bool(self.args.prepare_only),
            "profiles": [asdict(result) for result in results],
            "artifacts": {key: str(value) for key, value in artifacts.items()},
        }
        manifest_path = self.output_root / "manifest.json"
        manifest["artifacts"]["manifest"] = str(manifest_path)
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

        failures = [result for result in results if result.status == "failed"]
        if failures:
            raise PipelineExecutionError(
                f"{len(failures)} profile(s) failed; partial report: {artifacts['html_report']}"
            )
        self.reporter.info(f"Workflow completed: {artifacts['html_report']}", color="g")
        return manifest

    def _run_profile(
        self, profile: Profile, prepared: PreparedInputs, variants: Sequence[Variant]
    ) -> ProfileResult:
        from bespoke import bespoke2
        from educe_variant import educe_variant
        from oracode import oracode
        from remede import remede

        directory = self.output_root / profile.genome / profile.platform
        directory.mkdir(parents=True, exist_ok=True)
        final_bam = directory / "mutated.bam"
        validation_dir = directory / "validation"
        validation_vcf = validation_dir / "mutated.vcf.gz"
        pipeline_log = str(self.output_root / (self.args.logfile or "pipeline.log"))
        roi_header_normalized = False

        if not (self.args.resume and bam_has_reads(final_bam)):
            target_bed = prepared.directory / f"targets.{profile.platform}.bed"
            bespoke2(
                SimpleNamespace(
                    command="bespoke",
                    bam=str(profile.bam),
                    bed=[str(target_bed)],
                    outdir=str(directory),
                    logfile=pipeline_log,
                ),
                self.reporter,
            )
            roi_bam = directory / f"roi.{profile.bam.stem}.bam"
            if not roi_bam.is_file():
                raise PipelineExecutionError(f"ROI BAM was not created: {roi_bam}")
            roi_header_normalized = normalize_roi_header_for_reference(
                roi_bam, profile.reference, self.reporter
            )

            oracode(
                SimpleNamespace(
                    command="oracode",
                    bam=str(roi_bam),
                    snv=str(prepared.snv_path),
                    indel=str(prepared.indel_path),
                    reference=str(profile.reference),
                    outdir=str(directory),
                    aligner=profile.aligner,
                    allow_partial=True,
                    require_paired=False,
                    force=True,
                    mindepth=self.args.mindepth,
                    min_mutated_reads=self.args.min_mutated_reads,
                    tag_reads=True,
                    logfile=pipeline_log,
                ),
                self.reporter,
            )
            bamsurgeon_bam = directory / f"bs.{profile.aligner}.{roi_bam.stem}.bam"
            mutation_bam = bamsurgeon_bam
            if profile.platform == "iontorrent":
                # HG001 IonTorrent BAM headers already contain the corrected chrM
                # length, so no reheader step belongs in this workflow.
                remede(
                    SimpleNamespace(
                        command="remede",
                        bam=str(bamsurgeon_bam),
                        muts=[str(prepared.snv_path), str(prepared.indel_path)],
                        outdir=str(directory),
                        logfile=pipeline_log,
                        floworder="TCAG",
                        barcode="CAGATCCATCGAT",
                        library_key="TACGTACGTCTGAGCATCGATCGATGTACAGC",
                        except_reads_by_name=None,
                    ),
                    self.reporter,
                )
                mutation_bam = directory / f"reflow.{bamsurgeon_bam.stem}.bam"
            if not mutation_bam.is_file():
                raise PipelineExecutionError(f"Mutated BAM was not created: {mutation_bam}")
            _link_final_bam(mutation_bam, final_bam)

        if not (self.args.resume and validation_vcf.is_file()):
            validation_dir.mkdir(parents=True, exist_ok=True)
            target_bed = prepared.directory / f"targets.{profile.platform}.bed"
            educe_variant(
                SimpleNamespace(
                    command="educeVariant",
                    bam=str(final_bam),
                    bed=str(target_bed),
                    reference=str(profile.reference),
                    tool=profile.caller,
                    mindepth=self.args.mindepth,
                    min_mutated_reads=self.args.min_mutated_reads,
                    outdir=str(validation_dir),
                    logfile=pipeline_log,
                ),
                self.reporter,
            )
        if not validation_vcf.is_file():
            raise PipelineExecutionError(f"Validation VCF was not created: {validation_vcf}")

        detected = detect_variants(validation_vcf, variants, profile.genome)
        batch_path, screenshots_dir = create_igv_batch(
            profile, variants, final_bam, directory / "igv"
        )
        igv_status = "skipped" if self.args.skip_igv else run_igv(
            batch_path, self.args.igv_command, self.reporter
        )
        screenshot_count = sum(
            1 for path in screenshots_dir.glob("*.png") if path.stat().st_size > 0
        )

        if not self.args.keep_intermediates:
            _remove_intermediate_bams(directory, final_bam)

        return ProfileResult(
            profile.name,
            profile.genome,
            profile.platform,
            status="completed",
            bam=str(final_bam),
            validation_vcf=str(validation_vcf),
            igv_batch=str(batch_path),
            screenshots_dir=str(screenshots_dir),
            screenshot_count=screenshot_count,
            roi_header_normalized=roi_header_normalized,
            igv_status=igv_status,
            detected_variant_ids=sorted(detected),
        )


def write_reports(
    output_root: Path,
    variants_path: Path,
    variants: Sequence[Variant],
    profiles: Sequence[Profile],
    results: Sequence[ProfileResult],
) -> Dict[str, Path]:
    result_by_profile = {result.profile: result for result in results}
    csv_path = output_root / "variants.results.csv"
    profile_columns = [profile.name for profile in profiles]
    completed_profiles = {
        result.profile for result in results if result.status == "completed"
    }
    all_outputs_completed = set(profile_columns) == completed_profiles

    rows: List[Dict[str, object]] = []
    for variant in variants:
        row: Dict[str, object] = {
            "variant_id": variant.variant_id,
            "contig": variant.contig,
            "pos_hg19": variant.pos_hg19,
            "pos_hg38": variant.pos_hg38,
            "ref": variant.ref,
            "alt": variant.alt,
            "gene": variant.gene,
            "type": variant.kind,
            "vaf": variant.vaf,
        }
        detected_count = 0
        for profile_name in profile_columns:
            result = result_by_profile.get(profile_name)
            if result is None or result.status != "completed":
                state = "not_run"
            elif variant.variant_id in result.detected_variant_ids:
                state = "detected"
                detected_count += 1
            else:
                state = "missed"
            row[profile_name] = state
        row["detected_profiles"] = detected_count
        row["shared_all_outputs"] = (
            all_outputs_completed and detected_count == len(profile_columns)
        )
        rows.append(row)

    fieldnames = [
        "variant_id",
        "contig",
        "pos_hg19",
        "pos_hg38",
        "ref",
        "alt",
        "gene",
        "type",
        "vaf",
        *profile_columns,
        "detected_profiles",
        "shared_all_outputs",
    ]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    html_path = output_root / "report.html"
    _write_html_report(html_path, variants_path, rows, profiles, results)
    pdf_path = output_root / "report.pdf"
    _write_minimal_pdf(pdf_path, rows, results)
    return {
        "results_csv": csv_path,
        "html_report": html_path,
        "pdf_report": pdf_path,
    }


def _write_html_report(
    path: Path,
    variants_path: Path,
    rows: Sequence[Mapping[str, object]],
    profiles: Sequence[Profile],
    results: Sequence[ProfileResult],
) -> None:
    profile_names = [profile.name for profile in profiles]
    result_by_profile = {result.profile: result for result in results}
    summary_cards = []
    for profile in profiles:
        result = result_by_profile.get(profile.name)
        status = result.status if result else "not-run"
        count = len(result.detected_variant_ids) if result else 0
        summary_cards.append(
            f"<div class='card'><h3>{html.escape(profile.name)}</h3>"
            f"<p>Status: <strong>{html.escape(status)}</strong></p>"
            f"<p>Detected: {count}/{len(rows)}</p>"
            f"<p>IGV PNG: {result.screenshot_count if result else 0}/{len(rows)}</p></div>"
        )

    header = "".join(f"<th>{html.escape(name)}</th>" for name in profile_names)
    table_rows = []
    for row in rows:
        state_cells = []
        for name in profile_names:
            state = html.escape(str(row[name]))
            result = result_by_profile.get(name)
            value = state
            if result and result.screenshots_dir:
                variant = next(item for item in rows if item["variant_id"] == row["variant_id"])
                screenshot = Path(result.screenshots_dir) / (
                    str(variant["variant_id"]).replace(":", "_").replace(">", "-") + ".png"
                )
                if screenshot.is_file() and screenshot.stat().st_size > 0:
                    href = html.escape(os.path.relpath(screenshot, path.parent))
                    value = f"<a href='{href}'>{state}</a>"
            state_cells.append(f"<td class='{state}'>{value}</td>")
        states = "".join(state_cells)
        table_rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['variant_id']))}</td>"
            f"<td>{html.escape(str(row['gene']))}</td>"
            f"<td>{html.escape(str(row['contig']))}</td>"
            f"<td>{row['pos_hg19']}</td><td>{row['pos_hg38']}</td>"
            f"<td>{html.escape(str(row['ref']))}&gt;{html.escape(str(row['alt']))}</td>"
            f"<td>{html.escape(str(row['type']))}</td>{states}"
            f"<td>{'yes' if row['shared_all_outputs'] else 'no'}</td>"
            "</tr>"
        )

    path.write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>Oracode proficiency report</title><style>"
        "body{font-family:Arial,sans-serif;margin:2rem;color:#172033}"
        ".cards{display:flex;gap:1rem;flex-wrap:wrap}.card{border:1px solid #ccd3df;"
        "border-radius:8px;padding:0.75rem 1rem;min-width:180px}"
        "table{border-collapse:collapse;width:100%;font-size:12px;margin-top:1rem}"
        "th,td{border:1px solid #d8dee9;padding:5px;text-align:left}th{background:#eef2f7}"
        ".detected{background:#d9f6df}.missed{background:#ffe0e0}.not_run{background:#eee}"
        "</style></head><body>"
        "<h1>Oracode proficiency-test report</h1>"
        f"<p>Input: <code>{html.escape(str(variants_path))}</code></p>"
        f"<div class='cards'>{''.join(summary_cards)}</div>"
        "<h2>Variant introduction validation</h2><table><thead><tr>"
        "<th>ID</th><th>Gene</th><th>Contig</th><th>hg19</th><th>hg38</th>"
        f"<th>Allele</th><th>Type</th>{header}<th>Shared by all outputs</th>"
        f"</tr></thead><tbody>{''.join(table_rows)}</tbody></table>"
        "</body></html>",
        encoding="utf-8",
    )


def _pdf_escape(value: str) -> str:
    ascii_value = value.encode("ascii", "replace").decode("ascii")
    return ascii_value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _write_minimal_pdf(
    path: Path, rows: Sequence[Mapping[str, object]], results: Sequence[ProfileResult]
) -> None:
    lines = ["Oracode proficiency-test report", f"Variants: {len(rows)}"]
    for result in results:
        lines.append(
            f"{result.profile}: {result.status}, detected {len(result.detected_variant_ids)}/{len(rows)}, "
            f"IGV PNG {result.screenshot_count}/{len(rows)}"
        )
    shared = sum(1 for row in rows if row["shared_all_outputs"])
    lines.append(f"Shared by all outputs: {shared}/{len(rows)}")
    commands = ["BT", "/F1 11 Tf", "50 790 Td"]
    for index, line in enumerate(lines):
        if index:
            commands.append("0 -16 Td")
        commands.append(f"({_pdf_escape(line)}) Tj")
    commands.append("ET")
    stream = "\n".join(commands).encode("ascii")

    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] "
        b"/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        b"<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n" + stream + b"\nendstream",
    ]
    output = bytearray(b"%PDF-1.4\n")
    offsets = [0]
    for index, obj in enumerate(objects, start=1):
        offsets.append(len(output))
        output.extend(f"{index} 0 obj\n".encode("ascii"))
        output.extend(obj)
        output.extend(b"\nendobj\n")
    xref = len(output)
    output.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    output.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        output.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
    output.extend(
        f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
        f"startxref\n{xref}\n%%EOF\n".encode("ascii")
    )
    path.write_bytes(output)


def proficiency(args, reporter):
    return ProficiencyPipeline(args, reporter).run()
