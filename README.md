# VCF Filter – Comprehensive Somatic Variant Filtering & Splitter

A robust, production‑ready Python command‑line tool for filtering and splitting VCF (Variant Call Format) files from somatic variant callers (such as Mutect2, Strelka2, or VarScan2). It supports fine‑grained quality control, variant‑type classification, and multiple output modes—all with **no external dependencies**.

> **Repository:** [https://github.com/ahmadsam66/VCF_filter](https://github.com/ahmadsam66/VCF_filter)

---

## Table of Contents

- [Installation](#installation)
- [Domain Background](#domain-background)
- [Tool Functionality](#tool-functionality)
- [Command‑Line Options](#commandline-options)
- [Usage](#usage)
  - [Command‑Line Interface (CLI)](#commandline-interface-cli)
  - [Interactive Terminal Mode](#interactive-terminal-mode)
- [Example Workflows](#example-workflows)
- [License](#license)

---

## Installation

### Clone the repository (recommended)

```bash
git clone https://github.com/ahmadsam66/VCF_filter.git
cd VCF_filter
```

### Make the script executable

```bash
chmod +x vcf_filter.py
```

### Verify Python version

The script requires **Python 3.6+**. Check your version:

```bash
python3 --version
```

If you have an older version, install Python 3 from your system’s package manager (e.g., `apt install python3` on Ubuntu, `brew install python` on macOS).

### Optional: Add to PATH

To run the script from anywhere, you can add it to your `PATH` or create a symbolic link:

```bash
# Example: link to ~/bin (create if it doesn't exist)
mkdir -p ~/bin
ln -s $(pwd)/vcf_filter.py ~/bin/vcf_filter
export PATH="$HOME/bin:$PATH"   # add this to your .bashrc/.zshrc
```

Now you can simply run `vcf_filter` from any directory.

### Direct download (alternative)

If you prefer not to clone the entire repository:

```bash
curl -O https://raw.githubusercontent.com/ahmadsam66/VCF_filter/main/vcf_filter.py
chmod +x vcf_filter.py
```

---

## Domain Background

**VCF (Variant Call Format)** is the standard file format for storing genetic variation data. A VCF file consists of a header (lines starting with `#`) and data lines with the following required columns:

| Column | Name     | Description                                                                 |
|--------|----------|-----------------------------------------------------------------------------|
| 1      | CHROM    | Chromosome or contig name.                                                  |
| 2      | POS      | 1‑based position of the variant.                                            |
| 3      | ID       | Variant identifier (often `.`).                                             |
| 4      | REF      | Reference allele at that position.                                          |
| 5      | ALT      | Alternate allele(s) (comma‑separated for multi‑allelic sites).              |
| 6      | QUAL     | Phred‑scaled quality score for the variant call.                            |
| 7      | FILTER   | Filter status (`PASS` if accepted, otherwise a flag or `;`‑separated list). |
| 8      | INFO     | Additional information (key=value pairs).                                   |
| 9      | FORMAT   | Format specification for sample columns (e.g. `GT:AD:DP:AF`).               |
| 10+    | SAMPLES  | One or more sample columns, each with values corresponding to FORMAT.       |

The **INFO** and **FORMAT** fields contain many quality metrics used by somatic callers. This tool allows you to filter on those metrics (depth, mapping quality, allele frequency, strand bias, etc.) to obtain a high‑confidence set of somatic variants.

---

## Tool Functionality

- **Variant Classification** – automatically classifies each record as:
  - `SNV` (single‑nucleotide variant)
  - `Indel` (insertion/deletion)
  - `MNV` (multi‑nucleotide variant, e.g., two adjacent substitutions)
- **Quality and Depth Filters**
  - `QUAL` score (column 6)
  - Total read depth from `INFO/DP`
  - Sample‑level depth from `FORMAT/DP`
- **Mapping and Base Quality**
  - `INFO/MQ` – mapping quality
  - `INFO/MMQ` – median mapping quality
  - `INFO/MBQ` – median base quality
  - `INFO/MFRL` – median fragment length
  - `INFO/MPOS` – median distance from read end
- **Mutect2‑Specific Filters**
  - `TLOD` – tumor log‑odds
  - `NLOD` – normal log‑odds
  - `NALOD` – normal allele log‑odds
  - `ECNT` – event count (number of alt alleles in the haplotype)
  - Exclude short tandem repeats (STR)
- **Sample‑Level Allele Metrics**
  - `FORMAT/AF` – allele frequency in tumor
  - `FORMAT/AD` – allelic depth (total and alt read counts)
  - `FORMAT/F1R2` and `F2R1` – read counts supporting the variant on forward/reverse strands
- **Phasing Support** – optionally require that a variant is phased (presence of `PGT` or `|` in GT)
- **Output Flexibility**
  - Write all passing variants to a single VCF
  - Split passing variants into separate files for SNVs, Indels, and MNVs
  - Filter only a specific variant type (e.g., only SNVs)
- **Interactive Mode** – for users who prefer a guided setup

---

## Command‑Line Options

All options are available both as command‑line arguments and (in a subset) via the interactive mode.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-i`, `--input` | path | **required** | Input VCF file path. |
| `-o`, `--output` | path | auto‑generated | Exact output VCF filename (only when not splitting). |
| `-t`, `--tag` | string | `filtered` | Suffix tag for auto‑generated output filename (e.g., `my_vcf_filtered.vcf`). |
| `--var-type` | `{all, snv, indel, mnv}` | `all` | Keep only the specified variant type. |
| `--split-types` | flag | `False` | Split output into separate files: `_snvs.vcf`, `_indels.vcf`, `_mnvs.vcf`. |
| `--chrom` | string | `None` | Keep variants only on this chromosome (e.g., `chr1`). |
| `--pass-only` | flag | `False` | Keep only records with `FILTER == PASS`. |
| `--min-qual` | float | `None` | Minimum QUAL score (column 6). |
| `--min-mq` | float | `None` | Minimum mapping quality (INFO/MQ). |
| `--min-mmq` | float | `None` | Minimum median mapping quality (INFO/MMQ). |
| `--min-mbq` | float | `None` | Minimum median base quality (INFO/MBQ). |
| `--min-info-dp` | float | `None` | Minimum total depth from INFO/DP. |
| `--min-sample-dp` | int | `None` | Minimum tumor sample depth (FORMAT/DP). |
| `--min-tlod` | float | `None` | Minimum tumor log‑odds (INFO/TLOD, Mutect2). |
| `--min-nlod` | float | `None` | Minimum normal log‑odds (INFO/NLOD, Mutect2). |
| `--min-nalod` | float | `None` | Minimum normal allele log‑odds (INFO/NALOD, Mutect2). |
| `--max-ecnt` | int | `None` | Maximum event count (INFO/ECNT, Mutect2). |
| `--min-mfrl` | float | `None` | Minimum median fragment length (INFO/MFRL). |
| `--min-mpos` | float | `None` | Minimum median distance from read end (INFO/MPOS). |
| `--exclude-str` | flag | `False` | Exclude variants flagged as STR (short tandem repeats). |
| `--min-af` | float | `None` | Minimum tumor allele frequency (FORMAT/AF). |
| `--min-alt-ad` | int | `None` | Minimum alternate allele read count (FORMAT/AD). |
| `--min-f1r2` | int | `None` | Minimum F1R2 read count (FORMAT/F1R2). |
| `--min-f2r1` | int | `None` | Minimum F2R1 read count (FORMAT/F2R1). |
| `--require-phased` | flag | `False` | Keep only variants that are phased (PGT present or `|` in GT). |

---

## Usage

### Command‑Line Interface (CLI)

Run the script from the terminal with the input file and any desired filters.

**Basic syntax:**

```bash
./vcf_filter.py -i <input.vcf> [options]
```

**Example:** filter by quality and PASS flag:

```bash
./vcf_filter.py -i somatic.vcf --pass-only --min-qual 30
```

If you have a **gzipped** VCF, you can decompress on‑the‑fly and pipe it:

```bash
zcat somatic.vcf.gz | ./vcf_filter.py -i /dev/stdin --pass-only
```

### Interactive Terminal Mode

Run the script **without any arguments** to start an interactive session. It will guide you through the most common filtering parameters step by step.

```bash
./vcf_filter.py
```

You will be prompted for:
- Input VCF path
- Variant type selection (all, SNV only, Indel only, or split into separate files)
- Output filename or tag
- Each individual filter (press ENTER to skip any)

This mode is ideal for beginners or one‑off analyses.

---

## Example Workflows

Below are 10 practical examples, ranging from simple to complex.

### 1. Basic PASS filter

Keep only variants that passed all filters and are of high quality (QUAL >= 30).

```bash
./vcf_filter.py -i raw_somatic.vcf --pass-only --min-qual 30
```

### 2. Filter by chromosome

Extract all PASS variants on chromosome 2.

```bash
./vcf_filter.py -i raw_somatic.vcf --chrom chr2 --pass-only
```

### 3. Depth filters only

Require total INFO depth >= 20 and tumor sample depth >= 10.

```bash
./vcf_filter.py -i raw_somatic.vcf --min-info-dp 20 --min-sample-dp 10
```

### 4. Mutect2‑specific LOD filters

Apply recommended Mutect2 thresholds: TLOD >= 6, NLOD >= 3, NALOD >= 2, and ECNT <= 3.

```bash
./vcf_filter.py -i mutect2.vcf --min-tlod 6.0 --min-nlod 3.0 --min-nalod 2.0 --max-ecnt 3
```

### 5. Strand bias and read position filters

Use high‑stringency filters to reduce strand bias artifacts.

```bash
./vcf_filter.py -i somatic.vcf --min-mbq 30 --min-mfrl 50 --min-mpos 10
```

### 6. Allele frequency and alternate read count

Keep only variants with tumor AF >= 5% and at least 5 alternate reads.

```bash
./vcf_filter.py -i somatic.vcf --min-af 0.05 --min-alt-ad 5
```

### 7. Only SNVs

Extract only single‑nucleotide variants (ignore indels and MNVs).

```bash
./vcf_filter.py -i somatic.vcf --var-type snv
```

### 8. Split into separate SNV, Indel, and MNV files

This creates three output files: `somatic_snvs.vcf`, `somatic_indels.vcf`, and `somatic_mnvs.vcf`.

```bash
./vcf_filter.py -i somatic.vcf --split-types
```

### 9. Full production filter (combined)

Apply a comprehensive set of filters often used in somatic variant calling pipelines.

```bash
./vcf_filter.py -i tumor_normal.vcf \
  --pass-only \
  --min-qual 30 \
  --min-info-dp 20 \
  --min-sample-dp 10 \
  --min-tlod 6.0 \
  --min-nlod 3.0 \
  --min-nalod 2.0 \
  --max-ecnt 3 \
  --min-mbq 30 \
  --min-mfrl 50 \
  --min-mpos 10 \
  --min-af 0.05 \
  --min-alt-ad 5 \
  --exclude-str \
  --var-type all
```

### 10. Interactive mode

Launch the interactive terminal and follow the prompts to set up a custom filtering session.

```bash
./vcf_filter.py
```

---

## License

This project is distributed under the **MIT License**. See below for the full license text:

```
MIT License

Copyright (c) 2025 [Your Name or Organization]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
