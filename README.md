# TE CpG Density

`annotate_cpg_density.py` calculates three sequence-based CpG metrics for
genomic intervals. It was developed for transposable-element (TE) loci, but it
can be used with any set of genomic regions represented as a three-column BED
file.

For every input interval, the script retrieves the corresponding reference
sequence through `pybedtools`/`bedtools nuc`, counts cytosines, guanines, and
case-insensitive `CG` dinucleotides, and writes three BED-like output files.

## What the script calculates

For an interval of length *L*, let:

- *N<sub>CpG</sub>* be the observed number of `CG` dinucleotides;
- *N<sub>C</sub>* be the number of cytosines; and
- *N<sub>G</sub>* be the number of guanines.

The reported metrics are:

| Output suffix | Metric | Formula |
| --- | --- | --- |
| `_CpG-per-base.bed` | CpG per interval base | *N<sub>CpG</sub> / L* |
| `_CpG-per-GC.bed` | CpG per C or G base | *N<sub>CpG</sub> / (N<sub>C</sub> + N<sub>G</sub>)* |
| `_CpG-OE.bed` | CpG observed/expected ratio | *N<sub>CpG</sub> x L / (N<sub>C</sub> x N<sub>G</sub>)* |

The first metric is reported as `0.000000` for a zero-length interval. The
second is `NA` when an interval contains neither C nor G. Observed/expected is
`NA` when either the C count or G count is zero. Numeric values are written to
six decimal places.

## Requirements

- Python 3
- [BEDTools](https://bedtools.readthedocs.io/) available on `PATH`
- [pybedtools](https://daler.github.io/pybedtools/)
- Optional: [SAMtools](https://www.htslib.org/) to create a FASTA index

One reproducible installation route is Conda/Bioconda:

```bash
conda create --name te-cpg-density \
  --channel conda-forge --channel bioconda \
  python=3.10 pybedtools=0.12.0 bedtools=2.31.0 samtools
conda activate te-cpg-density
```

Alternatively, install the Python dependency with:

```bash
python -m pip install -r requirements.txt
```

BEDTools is a command-line dependency and must still be installed separately
and discoverable on `PATH`.

## Input

The command accepts one reference genome followed by one or more BED files:

```text
python annotate_cpg_density.py <genome.fa> <TE1.bed> [<TE2.bed> ...]
```

### Reference FASTA

Use the same genome assembly and chromosome naming convention used to generate
the BED coordinates. An accompanying FASTA index is recommended:

```bash
samtools faidx genome.fa
```

If `genome.fa.fai` is absent, the script prints a warning and continues.

### BED intervals

Input must be a headerless, tab-delimited, three-column BED file:

```text
chr1    100000    100250
chr1    205000    205430
chrX    75000     75400
```

BED coordinates are zero-based and end-exclusive. Intervals are processed
independently, including intervals that overlap.

> **Important:** the current parser uses the fixed `bedtools nuc` column
> positions produced from a three-column BED file. BED4, BED6, or wider input
> shifts those positions and is not supported by this version.

## Usage

From the repository directory:

```bash
python annotate_cpg_density.py genome.fa repeatmasker_TE_loci.bed
```

Multiple BED files can be processed in one command:

```bash
python annotate_cpg_density.py genome.fa \
  LINE_loci.bed SINE_loci.bed LTR_loci.bed
```

For each BED file, outputs are written to the current working directory using
the input basename. For example, `LINE_loci.bed` produces:

```text
LINE_loci_CpG-per-base.bed
LINE_loci_CpG-per-GC.bed
LINE_loci_CpG-OE.bed
```

Each output has four tab-delimited columns and no header:

```text
chrom    start    end    value
```

Example:

```text
chr1    100000    100250    0.024000
chr1    205000    205430    0.011628
```

## Interpretation

- **CpG per base** measures absolute CpG density across the whole interval.
- **CpG per C+G** normalizes the CpG count by GC-bearing bases and can help
  separate CpG abundance from interval length.
- **CpG observed/expected** compares observed CpGs with the count expected from
  the interval's marginal C and G composition. Values below 1 indicate CpG
  depletion and values above 1 indicate enrichment under this simple model.

CpG composition is reference-sequence based; it does not measure DNA
methylation. For TE analyses, interpret these values in the context of TE
family, locus age, genomic context, assembly version, and interval length.

## Behavioral details and limitations

- Output files are created in the directory from which the command is run, not
  beside the input BED file.
- Inputs from different directories that share a basename generate the same
  output names and can overwrite one another.
- Chromosome names in the BED and FASTA must match (for example, `chr1` versus
  `1`).
- Ambiguous reference bases contribute to interval length but not to C or G
  counts.
- Invalid or out-of-range intervals are handled by BEDTools; inspect BEDTools
  messages when an interval cannot be evaluated.
- Malformed `bedtools nuc` rows and rows without the expected columns are
  skipped, matching the original script behavior.

## Validation

The repository includes a synthetic FASTA/BED integration test with intervals
covering CpG-rich sequence, ambiguous bases, and an interval without C or G:

```bash
python -m unittest discover -s tests -v
```

The test requires both `pybedtools` and the `bedtools` executable.

## Acknowledgment

The script retains a simple positional command-line interface and the
comment-led procedural organization used by the Python utilities in
[TEProf2Paper](https://github.com/twlab/TEProf2Paper), while remaining a Python
3 program. Its CpG calculations and input/output behavior are unchanged from
the original `annotate_cpg_density.py` implementation.

## Citation

If this script contributes to a publication, cite the repository and the core
tools it calls:

- Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing
  genomic features. *Bioinformatics*. 2010;26(6):841-842.
  <https://doi.org/10.1093/bioinformatics/btq033>
- Dale RK, Pedersen BS, Quinlan AR. Pybedtools: a flexible Python library for
  manipulating genomic datasets and annotations. *Bioinformatics*.
  2011;27(24):3423-3424.
  <https://doi.org/10.1093/bioinformatics/btr539>
