# SpliceScape

**Version 2.1**

SpliceScape is a Nextflow pipeline for large-scale identification, quantification, and characterization of alternative splicing events from paired-end RNA-seq data.

The workflow integrates **BBDuk** (read cleaning), **STAR** (splice-aware alignment), **SGSeq** and **MAJIQ v3 / VOILA** (splicing analysis), and **SQLite** (integrated results and sample-level quality-control statistics). SpliceScape supports configurable concurrency and sample-level recovery for selected processing failures.

![SpliceScape workflow](images/GeneralViewPipeline.png)

## 1. Installation

### 1.1 Clone the repository

```bash
git clone https://github.com/labbces/SpliceScape.git
cd SpliceScape
```

### 1.2 Create the software environment

The repository provides [`envs/environment_majiqv3.yml`](envs/environment_majiqv3.yml), including Python, R/Bioconductor (SGSeq), BBDuk/BBMap, STAR, samtools, and Nextflow 25.04.8.

**Before creating the environment:** the YAML currently includes installation-specific `name` and `prefix` entries (`/Storage/progs/SpliceScape_env_p312`). Remove those two lines from your local copy of the YAML when installing outside the original development system.

```bash
# Make a local, portable copy of the supplied specification.
sed '/^name:/d; /^prefix:/d' envs/environment_majiqv3.yml > environment_majiqv3.local.yml

conda env create \
  --prefix /path/to/SpliceScape_env \
  --file environment_majiqv3.local.yml

conda activate /path/to/SpliceScape_env
```

Replace `/path/to/SpliceScape_env` with an absolute, writable installation path. Mamba may be used instead of Conda.

Check core tools:

```bash
python --version
nextflow -version
STAR --version
samtools --version
Rscript --version
```

### 1.3 Install MAJIQ v3 separately

**MAJIQ v3 is not included in the Conda YAML** and must be installed according to the MAJIQ developers' instructions. Obtain the appropriate MAJIQ license and follow the official installation procedure:

- [MAJIQ official website](https://majiq.biociphers.org/)
- [MAJIQ download and licensing](https://majiq.biociphers.org/app_download/)
- [MAJIQ installation documentation](https://biociphers.bitbucket.io/majiq-docs/getting-started-guide/installing.html)

Install MAJIQ/VOILA into the environment used to run SpliceScape and confirm that the MAJIQ command is accessible:

```bash
majiq --version
```

Keep the license file in a secure location. Its path is provided in the run-specific JSON configuration.

## 2. Input files

The pipeline requires:

1. **Genome FASTA** — uncompressed `.fa` or `.fasta`.
2. **Genome annotation** — uncompressed GFF3/GTF with gene and exon annotations suitable for the splicing tools.
3. **Sample manifest** — text file with one accession per line.

Example `samples.txt`:

```text
SRR17687198
SRR23406859
SRR12548574
```

The currently active workflow uses **`WGET_DOWNLOADER`**, which retrieves paired FASTQ files from the server configured in the JSON. Provide a URL from which files can be accessed using the sample accession, with corresponding `_1.fastq.gz` and `_2.fastq.gz` read files. **The active workflow does not download directly from NCBI based solely on SRA accessions.**

For background on reference-file preparation, see the [SpliceScape Wiki](https://github.com/labbces/SpliceScape/wiki).

## 3. JSON configuration

SpliceScape v2.1 uses two complementary configurations:

- A **JSON file** supplies sample inputs, reference paths, biological parameters, output locations, cleanup options, credentials, and concurrency.
- [`reads_processing/splicescape-v2.config`](reads_processing/splicescape-v2.config) defines Nextflow executor profiles and CPU/memory allocations. Adjust its SGE queue and resources for your installation.

Create a JSON file, for example `splicescape.json`:

```json
{
  "species": "Athaliana_447",
  "species_name": "Arabidopsis thaliana",
  "reads_file": "/data/project/samples.txt",
  "genomeFASTA": "/data/genomes/Athaliana_447_TAIR10.fa",
  "genomeGFF": "/data/genomes/Athaliana_447_Araport11.gene_exons.gff3",

  "outdir": "/data/project/results",
  "backbonedir": "/data/project/reference",
  "db": "/data/project/results/db/splicescape.sqlite",

  "max_active_samples": 3,
  "threads": 12,
  "sgseq_cores": 8,
  "majiq_cores": 12,

  "bbduk": "/path/to/SpliceScape_env/bin/bbduk.sh",
  "minlength": 60,
  "trimq": 20,
  "k": 27,
  "rref": "/data/references/adapters.fa",
  "maxmem": "20g",

  "r_libs": "/path/to/R/library",

  "url": "https://example.org/rnaseq/",
  "user": "download_username",
  "password": "download_password",

  "rm_reads_scratch_raw": true,
  "rm_reads_scratch_clean": true,
  "rm_reads_inspect_bbduk_logs": false,
  "rm_mapping_inspect_bamfiles": false,
  "rm_majiq_scratch": false,

  "environment": {
    "conda_env": "/path/to/SpliceScape_env",
    "majiq_license_file": "/path/to/majiq_license.lic"
  }
}
```

**Replace the illustrative paths, URL, credentials, and filenames.** The `bbduk` field must point to the actual BBDuk executable installed on your system; its name may vary. Prefer absolute paths and ensure output and database directories are writable. Protect configuration files containing credentials, for example with `chmod 600 splicescape.json`. Do not commit credentials or license files.

### 3.1 Input and output parameters

| Parameter | Description |
|---|---|
| `species` | Reference identifier used in output names, indexing, and database records (e.g. `Athaliana_447`). |
| `species_name` | Descriptive organism name stored with the database records. |
| `reads_file` | File with one paired-end sample accession per line. |
| `genomeFASTA` | Path to uncompressed reference FASTA. |
| `genomeGFF` | Path to uncompressed reference annotation. |
| `outdir` | Main directory for outputs and intermediate files. |
| `backbonedir` | Directory for shared reference resources such as the STAR index. |
| `db` | SQLite database path. |

### 3.2 Concurrency and processing resources

| Parameter | Description |
|---|---|
| `max_active_samples` | Maximum number of admitted samples in their complete processing lifecycle; default **2**. |
| `threads` | Threads passed to STAR genome indexing and alignment. |
| `sgseq_cores` | Core-count argument for SGSeq. |
| `majiq_cores` | Core-count argument for MAJIQ/VOILA stages. |

For example, `"max_active_samples": 3` allows at most three concurrently admitted samples. When a sample finishes or reaches a **supported** terminal failure state, its slot is returned and another sample may start. Actual task scheduling is also limited by the execution profile and available cluster resources.

### 3.3 BBDuk parameters

| Parameter | Description |
|---|---|
| `bbduk` | Full path to the BBDuk executable or launch script. |
| `minlength` | Minimum length retained after trimming. |
| `trimq` | Quality threshold used for read trimming. |
| `k` | K-mer size for matching the contaminant/reference FASTA. |
| `rref` | Adapter or contaminant reference FASTA. |
| `maxmem` | Java memory allocation for BBDuk, e.g. `20g`. |

### 3.4 Environment and download parameters

| Parameter | Description |
|---|---|
| `r_libs` | Directory containing the required R libraries. |
| `environment.conda_env` | Conda environment prefix; its `bin` directory is prepended to `PATH` by the launcher. |
| `environment.majiq_license_file` | Path to the MAJIQ license; the launcher exports `MAJIQ_LICENSE_FILE`. |
| `url` | Base URL for paired FASTQ download via `WGET_DOWNLOADER`. |
| `user` | Download server username. |
| `password` | Download server password. |

### 3.5 Cleanup parameters

| Parameter | Description |
|---|---|
| `rm_reads_scratch_raw` | Delete raw downloaded FASTQs after successful trimming. |
| `rm_reads_scratch_clean` | Delete cleaned FASTQs after successful mapping. |
| `rm_reads_inspect_bbduk_logs` | Control removal of BBDuk inspection logs in the applicable pipeline path. |
| `rm_mapping_inspect_bamfiles` | Delete intermediate mapped BAMs during successful sample finalization. |
| `rm_majiq_scratch` | Control removal of MAJIQ scratch files. |

All cleanup settings are JSON booleans (`true` or `false`). Retaining intermediate files is useful for troubleshooting. BBDuk and STAR summary logs must remain available until QC import has completed.

## 4. Running SpliceScape

The supported entry point is [`reads_processing/run_splicescape.py`](reads_processing/run_splicescape.py). Run these commands from the repository root.

### 4.1 SGE cluster

```bash
python3 reads_processing/run_splicescape.py \
  --params-file /absolute/path/to/splicescape.json \
  -profile sge
```

Before running, edit the SGE queue and per-process resource requests in [`reads_processing/splicescape-v2.config`](reads_processing/splicescape-v2.config) if needed. The supplied SGE profile currently uses `introns.q`.

### 4.2 Local execution

```bash
python3 reads_processing/run_splicescape.py \
  --params-file /absolute/path/to/splicescape.json \
  -profile local
```

Ensure the local machine has enough RAM and CPU capacity for STAR, SGSeq, and MAJIQ tasks.

### 4.3 Resume a previous execution

```bash
python3 reads_processing/run_splicescape.py \
  --params-file /absolute/path/to/splicescape.json \
  -profile sge \
  --resume
```

**Resume caveat:** Some intermediate files are kept outside the Nextflow work directories and may have been removed by cleanup options. Verify that required external artifacts are still present before resuming. For an independent validation run, use a fresh output directory without `--resume`.

## 5. Results and QC statistics

SpliceScape writes BBDuk preprocessing outputs, STAR mapping outputs, SGSeq and MAJIQ results, and an integrated SQLite database.

For samples whose **SGSeq and MAJIQ analyses and database import complete successfully**, the SQLite `sample_qc_metrics` table contains numeric BBDuk and STAR summary statistics:

```sql
SELECT sra_id, tool, metric, numeric_value, text_value
FROM sample_qc_metrics
WHERE sra_id = 'SRR17687198'
ORDER BY tool, metric;
```

QC metrics are stored using an idempotent UPSERT keyed by sample accession, genome version, tool, and metric. Re-importing the same metrics updates existing rows without duplicating them.

## 6. Failure handling and known limitations

Version 2.1 has sample-level failure detection for `WGET_DOWNLOADER`, `ALTERNATIVE_RUN_BBDUK`, `MAPPING_STAR`, `SGSEQ`, and `MAJIQ_SJ`. A supported sample failure is recorded without necessarily aborting the remaining batch.

Recovery for later MAJIQ stages, `POPULATE_DB`, `FINALIZE_SAMPLE`, scheduler-level failures, and some resume scenarios is **not yet implemented**. See [`reads_processing/TODO.md`](reads_processing/TODO.md).

## 7. Documentation

- [SpliceScape Wiki](https://github.com/labbces/SpliceScape/wiki)
- [MAJIQ documentation](https://biociphers.bitbucket.io/majiq-docs/)
- [MAJIQ download and licensing](https://majiq.biociphers.org/app_download/)
