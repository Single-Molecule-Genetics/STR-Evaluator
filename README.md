# STR-Evaluator

**STR-Evaluator** is a Python pipeline for detecting and evaluating **short tandem repeat (STR) variants** from paired-end (PE) sequencing data, leveraging duplex consensus sequence (DCS) tags for higher confidence. It assigns evidence **tiers** to STR calls, computes site-level frequencies, applies optional Bayesian filtering, and produces both **machine-readable outputs** (pickles, Excel summaries) and **visual quality control plots** (PDFs of stutter/frequency vs tier distributions).

---

## Purpose

- Detect STR events from PE reads, restrict to tags seen in DCS for robustness.
- Resolve ambiguous calls (“ties”) per tag × position with user-selectable policies.
- Quantify per-tag evidence across mates/strands.
- Assign **tiers** reflecting confidence and evidence patterns.
- Summarize coverage, allele frequencies, and heterology at STR sites.
- Visualize stutter vs tier distributions per STR locus.

---

## Inputs

- **STR reference pickle** (`chromosome_STRdict`):
  - Maps genomic STR spans to annotation strings (chr-start-end-motif).
- **DCS BAM**:
  - Duplex consensus sequences, used to define the set of trusted read tags.
- **PE BAM**:
  - Paired-end reads to scan for STR motifs.
- **Regions**:
  - Comma-separated chromosome list to process (e.g. `chr1,chr2,chr3`).
- **Optional cross/genotype Excel**:
  - Provides expected STRs for engineered strains.
  - Requires `--strain1` and `--strain2` arguments.
- **Parameters**:
  - Minimum sequence length (default: 60 bp).
  - Tie resolution policy (`favor`, `random`, `remove`).

---

## Workflow

### 1. Setup
- Creates a `results/` directory with subfolders for plots.
- Optionally loads cross/genotype Excel to collect “exception” STRs.

### 2. How STRs are found
- **Reference STR catalog**:  
  Uses a precomputed pickle (`chromosome_STRdict`) listing known STR loci across the genome.
- **Read scanning**:  
  Each PE read is scanned for repeat motifs with regex thresholds:  
  - mononucleotide ≥ 5 repeats,  
  - dinucleotide ≥ 6 repeats,  
  - trinucleotide ≥ 4 repeats.  
- **Mapping to reference**:  
  Candidate repeat runs are mapped to genomic positions via the BAM alignment. Only repeats overlapping loci in `chromosome_STRdict` are kept.  
- **Normalization**:  
  STR calls are standardized into an **id** string:  
  `chr-pos-REF-ALT`  
  where REF = reference motif, ALT = observed motif and length.  
- **Anchoring with DCS tags**:  
  Only PE read tags also found in the DCS BAM are retained, ensuring higher-confidence support.

### 3. Tie resolution
- Groups by `(start_pos, tag)` to identify cases where one tag supports multiple STR calls.
- Applies `clean_tied_data` policy:
  - **favor** expected alleles (from cross table),
  - **random** select one,
  - **remove** drop tied events.
- Produces a tie-free dataframe.

### 4. Tag-level aggregation
- For each tag × STR id:
  - Counts evidence in each mate/strand code (`ab.1`, `ba.2`, `ab.2`, `ba.1`).
  - Computes within-strand percentages.
- Collects per-tag totals across all orientations.

### 5. Tier assignment
- Each tag × STR id is classified into a **tier** that reflects:
  - Strand/mate balance,
  - Family size,
  - Frequency proportions.
- Tiers are ordered consistently and color-coded for plotting.
- A `tier_score` provides a monotonic ranking.

### 6. Site-level metrics
- Computes position-level metrics:
  - **Coverage** (number of tags),
  - **VAF** (allele fraction per id),
  - **Heterology** flag (mid-range VAF),
  - **Repeat size/type** (short/long).
- Optional **Bayesian filtering**:
  - Combines tier score, coverage, VAF, and stutter priors to estimate `P(True STR)`.

### 7. Summaries
- **Raw tie-free pickle**: per-tag STR calls with alt, VAF, heterology.
- **Summary pickle**: detailed table with tiers, counts, frequencies, coverage.
- **Compressed Excel**: per-locus summaries with AFs and tier distributions.

### 8. Visualization
- Per STR locus:
  - Tabulates family size distance (FS) vs tier contributions.
  - Produces stacked bar plots showing tier breakdown.
- Outputs individual PDFs per locus and merged region-level PDFs.

---

## Outputs

- **Pickles (`results/*.pkl`)**
  - Tie-free raw calls.
  - Main summary dataframe.
  - Cached intermediate DataFrames.
- **Excel summaries**
  - Compressed per-locus frequency tables with conditional formatting by tier.
- **PDF plots**
  - Stutter vs tier visualizations per STR locus.
  - Optionally merged PDFs per region.

---

## Key Concepts

- **id format**: `chr-start-ref-altMotif` (encodes locus and observed STR).
- **tag subtypes**: `ab.1`, `ba.2`, `ab.2`, `ba.1` (mates × strands).
- **tag_nr**: `1` = duplex-consistent, `2` = other patterns.
- **VAF**: allele fraction from distinct tag calls, not read depth.
- **Heterology**: positions with VAF ≈ 0.4–0.7 flagged as heterozygous-like.
- **Tiers**: evidence grading system encoding strength and symmetry.

---
## Example command

```bash
python3 STR-Evaluator.py \
  --str_finder_pkl STRs_index.pkl \
  --dcsbam DCS.bam \
  --pebam PE.bam \
  --regions chr1,chr2,chr3 \
  --minseqlength 60
```

## Performance

- **Parallelized** using `ProcessPoolExecutor`:
  - Region-level scanning,
  - Tie resolution,
  - Tag aggregation.
- **Memory-conscious**:
  - Uses staged concatenation (`fast_concat`).
  - Frees intermediates explicitly (`del_variable`).

---
## License
This project is licensed under the **MIT License**. You are free to use, modify, and distribute it for research and academic purposes.  

## Contributions
We welcome contributions! To get started:  
1. Fork the repository.  
2. Create a feature branch.  
3. Submit a pull request with a clear description of your changes.  

## About Us  
This tool was developed by **Shehab Moukbel ALi Aldawla** as part of the **Single Molecule Genetics (SMG) team** at JKU.  
For more information, visit our [team page](https://www.jku.at/institut-fuer-biophysik/ueber-uns/team/single-molecule-genetics-irene-tiemann-boege/) or join our community discussions. 

## Citation
If you use **STR-Evaluator** in your research, please cite it as:

> Shehab Moukbel Ali Aldawla. *STR-Evaluator: Detection, tier-based classification, and visualization of short tandem repeat variants from paired-end sequencing*.  
> GitHub repository: https://github.com/Single-Molecule-Genetics/STR-Evaluator  
> Version: 1.0.0 (2025)