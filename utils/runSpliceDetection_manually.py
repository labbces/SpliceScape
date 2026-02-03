#!/usr/bin/env python3
"""
SRA -> (download) -> BBDuk -> STAR -> (MAJIQ + SGSeq) with resume + cleanup of intermediates.

Key behavior (as requested):
- After producing cleaned fastq.gz, delete original raw fastq.gz
- After producing BAM and validating it, delete cleaned fastq.gz
- After finishing BOTH MAJIQ and SGSeq, delete BAM

External tools expected on PATH (depending on chosen options):
- curl or wget (URL mode)
- prefetch + fasterq-dump (SRA toolkit, SRA mode)
- bbduk.sh (BBTools)
- STAR
- majiq
- Rscript (+ SGSeq installed in R)
- samtools (optional but recommended for BAM validation)

Example usage:
  # SRA toolkit download (paired inferred by SRA):
  ./pipeline.py --sra SRR12345678 --download-mode sra \
    --outdir results \
    --adapters adapters.fa \
    --star-index /path/to/STARindex \
    --gtf /path/to/annotation.gtf \
    --majiq-config /path/to/majiq_config.ini \
    --threads 16

#TODO: Update URL, only the base URL must be provided.
  # Direct URL download (paired):
  ./pipeline.py --sra SRR123 --download-mode url \
    --url-r1 https://.../sample_1.fastq.gz \
    --url-r2 https://.../sample_2.fastq.gz \
    --outdir results \
    --adapters adapters.fa \
    --star-index /path/to/STARindex \
    --gtf /path/to/annotation.gtf \
    --majiq-config /path/to/majiq_config.ini \
    --threads 16
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence
from urllib.parse import urljoin
import urllib.request
from urllib.error import URLError, HTTPError
import json
import time

LOG_LEVEL: int = 1

# ---------------------------
# Helpers
# ---------------------------

def set_log_level(level: int) -> None:
    """Ajusta o nível global de log. Use 1=INFO, 2=WARN, 3=DEBUG."""
    global LOG_LEVEL
    try:
        LOG_LEVEL = int(level)
    except Exception:
        # ignore invalid values and keep previous level
        pass


def log(level: int, msg: str, *args) -> None:
    """Imprime mensagens dependendo do nível atual de log.

    Parâmetros:
      level: 1 (INFO), 2 (WARN), 3 (DEBUG). Mensagens com level > LOG_LEVEL são ignoradas.
      msg: string, pode conter formatação estilo printf se args forem fornecidos.
    """
    try:
        if level > LOG_LEVEL:
            return
    except Exception:
        # se LOG_LEVEL for inválido por algum motivo, recupera comportamento seguro
        return

    labels = {1: "INFO", 2: "WARN", 3: "DEBUG"}
    label = labels.get(level, "INFO")
    try:
        message = msg % args if args else msg
    except Exception:
        # fallback simples se formatação falhar
        message = "{}".format(msg) + (" " + " ".join(map(str, args)) if args else "")

    print(f"[{label}] {message}", flush=True)

def star_index_dir(outdir: Path, species: str) -> Path:
    """
    Shared STAR index directory, species-specific.
    Example: <outdir>/STAR_index_hs
    """
    safe = species.strip().lower()
    if not safe.isalnum():
        raise RuntimeError(f"Invalid species abbreviation: {species}")
    return outdir / f"STAR_index_{safe}"

def star_index_marker(idx_dir: Path) -> Path:
    return idx_dir / ".star_index.done"

def star_index_meta(idx_dir: Path) -> Path:
    return idx_dir / ".star_index.meta.json"

def file_fingerprint(path: Path) -> dict:
    """
    Fast fingerprint (no hashing): size + mtime + resolved path.
    Good enough to detect most mismatches without expensive hashing.
    """
    st = path.stat()
    return {
        "path": str(path.resolve()),
        "size": int(st.st_size),
        "mtime": float(st.st_mtime),
    }


def ensure_dir(d: Path) -> None:
    d.mkdir(parents=True, exist_ok=True)

def file_ok(p: Path, min_bytes: int = 1) -> bool:
    try:
        return p.is_file() and p.stat().st_size >= min_bytes
    except FileNotFoundError:
        return False

def safe_unlink(p: Optional[Path]) -> None:
    if not p:
        return
    try:
        p.unlink()
    except FileNotFoundError:
        pass

def which_or_die(exe: str) -> str:
    path = shutil.which(exe)
    if not path:
        raise RuntimeError(f"Required executable not found on PATH: {exe}")
    return path

def run_cmd(
    cmd: Sequence[str],
    *,
    cwd: Optional[Path] = None,
    stdout_path: Optional[Path] = None,
    stderr_path: Optional[Path] = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    if stdout_path:
        ensure_dir(stdout_path.parent)
    if stderr_path:
        ensure_dir(stderr_path.parent)
    log(3, "run_cmd cmd: %s", cmd)
    stdout_handle = open(stdout_path, "wb") if stdout_path else None
    stderr_handle = open(stderr_path, "wb") if stderr_path else None
    try:
        return subprocess.run(
            list(cmd),
            cwd=str(cwd) if cwd else None,
            stdout=stdout_handle if stdout_handle else None,
            stderr=stderr_handle if stderr_handle else None,
            check=check,
        )
    finally:
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle:
            stderr_handle.close()

def touch(path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text("ok\n")

def step_done(marker: Path) -> bool:
    return file_ok(marker, 2)

def gzip_ok(path: Path, nbytes: int = 1024 * 1024) -> bool:
    if not file_ok(path, 10):
        return False
    try:
        with gzip.open(path, "rb") as fh:
            fh.read(nbytes)
        return True
    except Exception:
        return False

def bam_ok(bam: Path) -> bool:
    """Prefer samtools quickcheck if available; otherwise do a size check."""
    if not file_ok(bam, 1000):
        return False
    samtools = shutil.which("samtools")
    if samtools:
        cp = subprocess.run([samtools, "quickcheck", "-v", str(bam)], capture_output=True)
        return cp.returncode == 0
    return True


# ---------------------------
# Paths
# ---------------------------

@dataclass
class SamplePaths:
    sample_id: str
    work_dir: Path

    # raw fastq.gz
    raw_r1: Path
    raw_r2: Optional[Path]

    # cleaned fastq.gz
    clean_r1: Path
    clean_r2: Optional[Path]
    bbduk_log: Path

    # STAR
    star_dir: Path
    star_prefix: str
    star_log: Path
    bam: Path

    # MAJIQ
    majiq_dir: Path
    majiq_log: Path

    # SGSeq
    sgseq_dir: Path
    sgseq_log: Path

    # markers
    mk_download: Path
    mk_bbduk: Path
    mk_star: Path
    mk_majiq: Path
    mk_sgseq: Path

def build_paths(sample_id: str, out_dir: Path, paired: bool, species: str) -> SamplePaths:
    work = out_dir / sample_id
    ensure_dir(work)

    raw_dir = work / "01_raw"
    clean_dir = work / "02_bbduk"
    star_dir = work / "03_star"
    majiq_dir = work / "04_majiq"
    sgseq_dir = work / "05_sgseq"
    markers = work / "markers"

    for d in (raw_dir, clean_dir, star_dir, majiq_dir, sgseq_dir, markers):
        ensure_dir(d)

    raw_r1 = raw_dir / f"{sample_id}_1.fastq.gz"
    raw_r2 = (raw_dir / f"{sample_id}_2.fastq.gz") if paired else None

    clean_r1 = clean_dir / f"{sample_id}.trimmed.R1.fastq.gz"
    clean_r2 = (clean_dir / f"{sample_id}.trimmed.R2.fastq.gz") if paired else None
    bbduk_log = clean_dir / f"{sample_id}.bbduk.log.txt"

    # STAR writes files under prefix; final BAM is {prefix}.Aligned.sortedByCoord.out.bam
    star_prefix = f"{species.lower()}_{sample_id}."
    star_log = star_dir / f"{species.lower()}_{sample_id}.STAR.log.txt"
    bam = star_dir / f"{star_prefix}Aligned.sortedByCoord.out.bam"

    majiq_log = majiq_dir / f"{sample_id}.majiq.log.txt"
    sgseq_log = sgseq_dir / f"{sample_id}.sgseq.log.txt"

    return SamplePaths(
        sample_id=sample_id,
        work_dir=work,
        raw_r1=raw_r1,
        raw_r2=raw_r2,
        clean_r1=clean_r1,
        clean_r2=clean_r2,
        bbduk_log=bbduk_log,
        star_dir=star_dir,
        star_prefix=star_prefix,
        star_log=star_log,
        bam=bam,
        majiq_dir=majiq_dir,
        majiq_log=majiq_log,
        sgseq_dir=sgseq_dir,
        sgseq_log=sgseq_log,
        mk_download=markers / "01_download.done",
        mk_bbduk=markers / "02_bbduk.done",
        mk_star=markers / "03_star.done",
        mk_majiq=markers / "04_majiq.done",
        mk_sgseq=markers / "05_sgseq.done",
    )


# ---------------------------
# Step 1: Download
# ---------------------------

def download_url_to_file(
    url: str,
    out: Path,
    *,
    user: Optional[str] = None,
    password: Optional[str] = None,
    force: bool = False,
) -> None:
    """Download URL -> out, supporting basic auth. Uses curl/wget if available; falls back to urllib."""
    ensure_dir(out.parent)
    if not force and file_ok(out, 10):
        return

    curl = shutil.which("curl")
    wget = shutil.which("wget")

    if curl:
        cmd = ["curl", "-L", "--retry", "5", "--retry-delay", "2", "-o", str(out)]
        if user is not None and password is not None:
            cmd.extend(["-u", f"{user}:{password}"])
        cmd.append(url)
        run_cmd(cmd)
        return

    if wget:
        cmd = ["wget", "-O", str(out)]
        if user is not None and password is not None:
            cmd.extend(["--user", user, "--password", password])
        cmd.append(url)
        run_cmd(cmd)
        return

    # urllib fallback
    req = urllib.request.Request(url)
    if user is not None and password is not None:
        import base64
        token = f"{user}:{password}".encode("utf-8")
        auth = b"Basic " + base64.b64encode(token)
        req.add_header("Authorization", auth.decode("ascii"))

    try:
        with urllib.request.urlopen(req) as r, open(out, "wb") as w:
            shutil.copyfileobj(r, w)
    except HTTPError as e:
        raise RuntimeError(f"HTTP error {e.code} downloading {url}") from e
    except URLError as e:
        raise RuntimeError(f"URL error downloading {url}: {e.reason}") from e

def build_fastq_urls_from_base(
    url_base: str,
    sra_id: str,
    paired: bool,
    *,
    paired_suffix_r1: str = "_1.fastq.gz",
    paired_suffix_r2: str = "_2.fastq.gz",
    single_suffix: str = ".fastq.gz",
) -> tuple[str, Optional[str]]:
    """
    Given:
      url_base='https://my.server.com/data/'
      sra_id='SRR010101'

    Returns:
      R1='https://my.server.com/data/SRR010101_1.fastq.gz'
      R2='https://my.server.com/data/SRR010101_2.fastq.gz'  (if paired)

    Uses urljoin to avoid double slashes issues.
    """
    if not url_base.endswith("/"):
        url_base = url_base + "/"

    if paired:
        url_r1 = urljoin(url_base, f"{sra_id}{paired_suffix_r1}")
        url_r2 = urljoin(url_base, f"{sra_id}{paired_suffix_r2}")
        return url_r1, url_r2
    else:
        url_r1 = urljoin(url_base, f"{sra_id}{single_suffix}")
        return url_r1, None

def download_from_url_base(
    *,
    url_base: str,
    sra_id: str,
    paired: bool,
    paths: "SamplePaths",
    force: bool,
    user: Optional[str] = None,
    password: Optional[str] = None,
) -> None:
    """Download R1/R2 from a base URL using the SRA ID to form filenames."""
    if not force and step_done(paths.mk_download):
        return

    url_r1, url_r2 = build_fastq_urls_from_base(url_base, sra_id, paired)

    download_url_to_file(url_r1, paths.raw_r1, user=user, password=password, force=force)

    if paired:
        if paths.raw_r2 is None:
            raise RuntimeError("Pipeline paths are single-end but paired=True was requested.")
        download_url_to_file(url_r2, paths.raw_r2, user=user, password=password, force=force)

    # sanity check gzip readability
    if not gzip_ok(paths.raw_r1):
        raise RuntimeError(f"Downloaded R1 is not a readable gzip: {paths.raw_r1}")
    if paths.raw_r2 and not gzip_ok(paths.raw_r2):
        raise RuntimeError(f"Downloaded R2 is not a readable gzip: {paths.raw_r2}")

    touch(paths.mk_download)

def download_from_sra(sra_id: str, paths: SamplePaths, threads: int, force: bool) -> None:
    if not force and step_done(paths.mk_download):
        return

    which_or_die("prefetch")
    which_or_die("fasterq-dump")

    sra_cache = paths.work_dir / "sra_cache"
    tmp_fq = paths.work_dir / "tmp_fasterq"
    ensure_dir(sra_cache)
    ensure_dir(tmp_fq)

    # 1) prefetch
    run_cmd(["prefetch", sra_id, "--output-directory", str(sra_cache)])

    # 2) fasterq-dump
    run_cmd([
        "fasterq-dump", sra_id,
        "--outdir", str(tmp_fq),
        "--threads", str(threads),
        "--skip-technical",
        "--split-files",
    ])

    fq1 = tmp_fq / f"{sra_id}_1.fastq"
    fq2 = tmp_fq / f"{sra_id}_2.fastq"
    fq_single = tmp_fq / f"{sra_id}.fastq"

    pigz = shutil.which("pigz")

    def compress_move(src_fastq: Path, dst_gz: Path) -> None:
        ensure_dir(dst_gz.parent)
        if not force and file_ok(dst_gz, 10):
            return
        if pigz:
            # pigz -p threads -c src > dst
            with open(dst_gz, "wb") as out:
                cp = subprocess.run([pigz, "-p", str(threads), "-c", str(src_fastq)], stdout=out)
                if cp.returncode != 0:
                    raise RuntimeError(f"pigz failed on {src_fastq}")
        else:
            with open(src_fastq, "rb") as f_in, gzip.open(dst_gz, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    if fq1.exists():
        # SRA produced split files (_1, _2). Treat as paired in source.
        compress_move(fq1, paths.raw_r1)

        if paths.raw_r2 is not None:
            # User requested paired: require fq2
            if not fq2.exists():
                raise RuntimeError("Requested paired-end (--paired) but SRA did not produce _2.fastq.")
            compress_move(fq2, paths.raw_r2)
        else:
            # User requested single-end: ignore fq2 if present
            pass

    elif fq_single.exists():
        # SRA produced a single fastq
        if paths.raw_r2 is not None:
            raise RuntimeError("Requested paired-end (--paired) but SRA produced single-end FASTQ.")
        compress_move(fq_single, paths.raw_r1)

    else:
        raise RuntimeError(f"Could not find fasterq-dump outputs in {tmp_fq}")

    # sanity checks
    if not gzip_ok(paths.raw_r1):
        raise RuntimeError(f"SRA-derived file not readable gzip: {paths.raw_r1}")
    if paths.raw_r2 and not gzip_ok(paths.raw_r2):
        raise RuntimeError(f"SRA-derived file not readable gzip: {paths.raw_r2}")

    # cleanup tmp fasterq output
    shutil.rmtree(tmp_fq, ignore_errors=True)

    touch(paths.mk_download)


# ---------------------------
# Step 2: BBDuk (and delete raw)
# ---------------------------

def run_bbduk(
    paths: SamplePaths,
    adapters_fa: Path,
    threads: int,
    force: bool,
    delete_raw: bool = True,
    extra_args: Optional[list[str]] = None,
) -> None:
    which_or_die("bbduk.sh")

    if not force and step_done(paths.mk_bbduk):
        return

    if not file_ok(adapters_fa, 10):
        raise RuntimeError(f"Adapters file not found or empty: {adapters_fa}")

    cmd = [
        "bbduk.sh",
        f"in1={paths.raw_r1}",
        f"out1={paths.clean_r1}",
        f"ref={adapters_fa}",
        "ktrim=r",
        "k=23",
        "trimpolya=10",
        "minlen=25",
        "mink=11",
        "hdist=1",
        "qtrim=rl",
        "trimq=10",
        "tbo=t",
        f"threads={threads}",
        "showspeed=t",
    ]
    if paths.raw_r2 and paths.clean_r2:
        cmd.extend([f"in2={paths.raw_r2}", f"out2={paths.clean_r2}"])
    if extra_args:
        cmd.extend(extra_args)

    # capture both stdout and stderr into the same log
    run_cmd(cmd, stdout_path=paths.bbduk_log, stderr_path=paths.bbduk_log, check=True)

    # validate outputs before cleanup
    if not gzip_ok(paths.clean_r1):
        raise RuntimeError(f"BBDuk produced invalid gzip: {paths.clean_r1}")
    if paths.clean_r2 and not gzip_ok(paths.clean_r2):
        raise RuntimeError(f"BBDuk produced invalid gzip: {paths.clean_r2}")

    touch(paths.mk_bbduk)

    # delete raw fastqs
    if delete_raw:
        safe_unlink(paths.raw_r1)
        safe_unlink(paths.raw_r2)


# ---------------------------
# Step 3: STAR (and delete cleaned fastqs after BAM OK)
# ---------------------------

def build_or_reuse_star_index(
    *,
    outdir: Path,
    species: str,
    genome_fasta: Path,
    annotation: Path,
    threads: int,
    sjdb_overhang: int,
    force: bool,
    extra_args: Optional[list[str]] = None,
) -> Path:
    """
    Shared STAR index living under outdir/subdir.

    Reuse if:
      - .star_index.done exists
      - genomeParameters.txt exists
      - meta matches (genome, annotation, sjdbOverhang, STAR version, extra_args)

    Otherwise rebuild in-place.
    """
    which_or_die("STAR")

    if not file_ok(genome_fasta, 10):
        raise RuntimeError(f"Genome FASTA not found/empty: {genome_fasta}")
    if not file_ok(annotation, 10):
        raise RuntimeError(f"Annotation not found/empty: {annotation}")

    idx_dir = star_index_dir(outdir, species)
    ensure_dir(idx_dir)

    marker = star_index_marker(idx_dir)
    meta_path = star_index_meta(idx_dir)
    genome_params = idx_dir / "genomeParameters.txt"
    idx_log = idx_dir / "STAR.genomeGenerate.log.txt"

    # STAR version (best-effort)
    star_ver = None
    try:
        cp = subprocess.run(["STAR", "--version"], capture_output=True, text=True)
        if cp.returncode == 0:
            star_ver = (cp.stdout.strip() or cp.stderr.strip() or None)
    except Exception:
        pass

    desired = {
        "species": species.lower(),
        "genome_fasta": file_fingerprint(genome_fasta),
        "annotation": file_fingerprint(annotation),
        "sjdb_overhang": int(sjdb_overhang),
        "extra_args": extra_args or [],
        "star_version": star_ver,
    }

    # Reuse if ok
    if not force and file_ok(marker, 2) and file_ok(genome_params, 10) and file_ok(meta_path, 10):
        try:
            existing = json.loads(meta_path.read_text())
        except Exception:
            existing = {}

        mismatches = []
        for k in ("species", "genome_fasta", "annotation", "sjdb_overhang", "extra_args", "star_version"):
            if existing.get(k) != desired.get(k):
                mismatches.append(k)

        if not mismatches:
            log(1, "Reusing STAR index: %s", idx_dir)
            return idx_dir

        raise RuntimeError(
            "STAR index exists but does not match current inputs/params. "
            f"Mismatched: {', '.join(mismatches)}. "
            "Use a different --outdir / --star-index-subdir or rerun with --force."
        )

    # If marker/meta missing or incomplete index, rebuild
    if not force and file_ok(marker, 2) and not file_ok(genome_params, 10):
        log(2, "STAR index marker exists but genomeParameters.txt is missing/empty. Rebuilding: %s", idx_dir)

    log(1, "Building STAR index in: %s", idx_dir)

    cmd = [
        "STAR",
        "--runMode", "genomeGenerate",
        "--runThreadN", str(threads),
        "--genomeDir", str(idx_dir),
        "--genomeFastaFiles", str(genome_fasta),
        "--sjdbGTFfile", str(annotation),
        "--sjdbOverhang", str(sjdb_overhang),
    ]
    if extra_args:
        cmd.extend(extra_args)

    run_cmd(cmd, stdout_path=idx_log, stderr_path=idx_log, check=True)

    if not file_ok(genome_params, 10):
        raise RuntimeError(f"STAR index generation finished but missing/empty: {genome_params}")

    # Write meta + marker
    meta_path.write_text(json.dumps(desired, indent=2) + "\n")
    touch(marker)

    log(1, "STAR index ready: %s", idx_dir)
    return idx_dir


def run_star(
    paths: SamplePaths,
    star_index_dir: Path,
    threads: int,
    force: bool,
    delete_clean: bool = True,
    extra_args: Optional[list[str]] = None,
) -> None:
    which_or_die("STAR")

    if not force and step_done(paths.mk_star):
        return

    if not star_index_dir.is_dir():
        raise RuntimeError(f"STAR index dir not found: {star_index_dir}")

    reads = [str(paths.clean_r1)]
    if paths.clean_r2:
        reads.append(str(paths.clean_r2))

    cmd = [
        "STAR",
        "--genomeDir", str(star_index_dir),
        "--runThreadN", str(threads),
        "--readFilesIn", *reads,
        "--readFilesCommand", "zcat",
        "--outFileNamePrefix", str(paths.star_dir / paths.star_prefix),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outSAMstrandField", "intronMotif",
        "--limitBAMsortRAM", "0",
        "--twopassMode", "Basic",
    ]
    if extra_args:
        cmd.extend(extra_args)

    run_cmd(cmd, stdout_path=paths.star_log, stderr_path=paths.star_log, check=True)

    # validate BAM
    if not bam_ok(paths.bam):
        raise RuntimeError(f"STAR produced BAM that failed validation: {paths.bam}")

    touch(paths.mk_star)

    # delete cleaned fastq.gz
    if delete_clean:
        safe_unlink(paths.clean_r1)
        safe_unlink(paths.clean_r2)


# ---------------------------
# Step 4: MAJIQ
# ---------------------------

def run_majiq(
    paths: SamplePaths,
    gtf: Path,
    majiq_config: Path,
    threads: int,
    force: bool,
    extra_args_build: Optional[list[str]] = None,
    extra_args_psi: Optional[list[str]] = None,
) -> None:
    which_or_die("majiq")

    if not force and step_done(paths.mk_majiq):
        return

    if not file_ok(paths.bam, 1000):
        raise RuntimeError(f"BAM not found for MAJIQ step: {paths.bam}")
    if not file_ok(gtf, 10):
        raise RuntimeError(f"GTF not found: {gtf}")
    if not file_ok(majiq_config, 10):
        raise RuntimeError(f"MAJIQ config not found: {majiq_config}")

    build_dir = paths.majiq_dir / "build"
    psi_dir = paths.majiq_dir / "psi"
    ensure_dir(build_dir)
    ensure_dir(psi_dir)

    # MAJIQ build
    cmd_build = [
        "majiq", "build",
        str(gtf),
        "-c", str(majiq_config),
        "-j", str(threads),
        "-o", str(build_dir),
    ]
    if extra_args_build:
        cmd_build.extend(extra_args_build)

    run_cmd(cmd_build, stdout_path=paths.majiq_log, stderr_path=paths.majiq_log, check=True)

    # MAJIQ psi (expects *.majiq in build_dir)
    cmd_psi = [
        "majiq", "psi",
        str(build_dir),
        "-j", str(threads),
        "-o", str(psi_dir),
    ]
    if extra_args_psi:
        cmd_psi.extend(extra_args_psi)

    run_cmd(cmd_psi, stdout_path=paths.majiq_log, stderr_path=paths.majiq_log, check=True)

    # minimal expected outputs (varies by MAJIQ version; keep marker as truth)
    touch(paths.mk_majiq)


# ---------------------------
# Step 5: SGSeq (via generated R script)
# ---------------------------

SGSEQ_R_SCRIPT = r'''
suppressPackageStartupMessages({
  library(SGSeq)
  library(GenomicFeatures)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly=TRUE)
bam <- args[1]
gtf <- args[2]
outdir <- args[3]

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# Create TxDb from GTF
txdb <- makeTxDbFromGFF(gtf, format="gtf")

# Define sample info for SGSeq
si <- data.frame(sample_name="sample1", file_bam=bam, stringsAsFactors=FALSE)

# Prepare annotations and count features
# Note: You may want to tune which features/filters depending on your needs.
annot <- getFeatures(txdb)
feat <- analyzeFeatures(si, features=annot)

saveRDS(feat, file=file.path(outdir, "sgseq_features.rds"))

# Export a simple table of splice junctions (if present)
# Different SGSeq objects can be coerced; here we try a common extractor.
try({
  sj <- as.data.frame(getSpliceJunctions(feat))
  write.table(sj, file=file.path(outdir, "splice_junctions.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}, silent=TRUE)

cat("SGSeq finished.\n")
'''

def run_sgseq(
    paths: SamplePaths,
    gtf: Path,
    force: bool,
) -> None:
    which_or_die("Rscript")

    if not force and step_done(paths.mk_sgseq):
        return

    if not file_ok(paths.bam, 1000):
        raise RuntimeError(f"BAM not found for SGSeq step: {paths.bam}")
    if not file_ok(gtf, 10):
        raise RuntimeError(f"GTF not found: {gtf}")

    r_script_path = paths.sgseq_dir / "run_sgseq.R"
    ensure_dir(paths.sgseq_dir)
    r_script_path.write_text(SGSEQ_R_SCRIPT)

    cmd = ["Rscript", str(r_script_path), str(paths.bam), str(gtf), str(paths.sgseq_dir)]
    run_cmd(cmd, stdout_path=paths.sgseq_log, stderr_path=paths.sgseq_log, check=True)

    touch(paths.mk_sgseq)


# ---------------------------
# Final cleanup: delete BAM only if BOTH MAJIQ and SGSeq are done
# ---------------------------

def cleanup_bam_if_done(paths: SamplePaths, delete_bam: bool = True) -> None:
    if not delete_bam:
        return
    if step_done(paths.mk_majiq) and step_done(paths.mk_sgseq):
        safe_unlink(paths.bam)
        # common index names
        safe_unlink(paths.bam.with_suffix(paths.bam.suffix + ".bai"))  # .bam.bai
        safe_unlink(paths.bam.with_suffix(".bai"))                     # .bai


# ---------------------------
# Main
# ---------------------------

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="SRA/URL -> BBDuk -> STAR -> MAJIQ + SGSeq (resumable, with cleanup).")

    ap.add_argument("--sra", required=True, help="SRA run identifier (e.g., SRR12345678). Used as sample ID.")
    ap.add_argument("--download-mode", choices=["sra", "url"], required=True, help="Download from SRA Toolkit or from direct URLs.")
    ap.add_argument("--url-base", default=None,
                    help="Base URL (directory/prefix) where FASTQs live, e.g. 'https://my.server.com/data/'. "
                    "The script will append <SRA>_1.fastq.gz (and _2 if paired).")
    ap.add_argument("--paired", action="store_true",
                    help="If set, download paired-end: <SRA>_1.fastq.gz and <SRA>_2.fastq.gz. "
                    "Otherwise download single-end: <SRA>.fastq.gz (or customize in code).")
    ap.add_argument("--url-user", default=None, help="Username for basic auth (optional).")
    ap.add_argument("--url-pass", default=None, help="Password for basic auth (optional).")

    ap.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    ap.add_argument("--threads", type=int, default=8, help="Threads for tools that support it.")

    ap.add_argument("--adapters", required=True, type=Path, help="FASTA of adapters/contaminants for BBDuk ref=...")
    ap.add_argument("--genome-fasta", required=True, type=Path, help="Genome FASTA used to build STAR index.")

    ap.add_argument("--annotation", required=True, type=Path, help="Annotation file (GTF preferred; GFF3 also supported by STAR in many cases).")    # ap.add_argument("--majiq-config", required=True, type=Path, help="MAJIQ config ini (your existing working config).")

    ap.add_argument("--sjdb-overhang", type=int, default=149, help="STAR sjdbOverhang (usually readLength-1). Default 149 for 150bp reads.")

    ap.add_argument("--species", required=True, help="Species abbreviation used to name the shared STAR index (e.g. hs, mm, ath, zm).")

    ap.add_argument("--force", action="store_true", help="Force rerun steps even if markers exist.")
    ap.add_argument("--keep-raw", action="store_true", help="Do not delete raw fastq.gz after BBDuk.")
    ap.add_argument("--keep-clean", action="store_true", help="Do not delete cleaned fastq.gz after STAR.")
    ap.add_argument("--keep-bam", action="store_true", help="Do not delete BAM after MAJIQ+SGSeq are done.")

    return ap.parse_args()

def main() -> None:
    args = parse_args()

    sample_id = args.sra

    paired = args.paired
    paths = build_paths(sample_id, args.outdir, paired=paired, species=args.species)

    # 1) Download
    if args.download_mode == "url":
        if not args.url_base:
            raise RuntimeError("In URL mode you must provide --url-base, e.g. https://my.server.com/data/")

        download_from_url_base(
            url_base=args.url_base,
            sra_id=sample_id,
            paired=paired,
            paths=paths,
            force=args.force,
            user=args.url_user,
            password=args.url_pass,
        )
    else:
        # For SRA: we don't know single/paired until fasterq-dump runs.
        # To avoid mismatch, we try paired layout first; if it turns out single-end, we'll just produce raw_r1.
        download_from_sra(sample_id, paths, threads=args.threads, force=args.force)

    # 2) BBDuk (delete raw afterwards unless keep-raw)
    run_bbduk(
        paths,
        adapters_fa=args.adapters,
        threads=args.threads,
        force=args.force,
        delete_raw=(not args.keep_raw),
    )

    # 3) STAR (delete cleaned afterwards unless keep-clean)
    idx_dir = build_or_reuse_star_index(
        outdir=args.outdir,
        species=args.species,
        genome_fasta=args.genome_fasta,
        annotation=args.annotation,
        threads=args.threads,
        sjdb_overhang=args.sjdb_overhang,
        force=args.force,
    )

    run_star(
        paths,
        star_index_dir=idx_dir,
        threads=args.threads,
        force=args.force,
        delete_clean=(not args.keep_clean),
    )

    # # 4) MAJIQ
    # run_majiq(
    #     paths,
    #     gtf=args.gtf,
    #     majiq_config=args.majiq_config,
    #     threads=args.threads,
    #     force=args.force,
    # )

    # # 5) SGSeq
    # run_sgseq(
    #     paths,
    #     gtf=args.gtf,
    #     force=args.force,
    # )

    # Final cleanup (delete BAM only if BOTH done, unless keep-bam)
    cleanup_bam_if_done(paths, delete_bam=(not args.keep_bam))

if __name__ == "__main__":
    main()
