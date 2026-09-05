#!/usr/bin/env python3
# purpose: Annotate genomic intervals with three CpG-density measurements.
# usage: python annotate_cpg_density.py <genome.fa> <TE1.bed> [<TE2.bed> ...]
# argument structure: reference genome FASTA followed by one or more BED files

import os
import sys

from pybedtools import BedTool


if len(sys.argv) < 3:
    sys.stderr.write(
        "Usage: python annotate_cpg_density.py "
        "<genome.fa> <TE1.bed> [<TE2.bed> ...]\n"
    )
    sys.exit(1)

genome_fasta = sys.argv[1]
bed_files = sys.argv[2:]

# Verify that the genome FASTA index (.fai) exists, warn if not
if not os.path.isfile(genome_fasta):
    sys.stderr.write(
        "Error: Genome FASTA file '{}' not found.\n".format(genome_fasta)
    )
    sys.exit(1)

if not os.path.isfile(genome_fasta + ".fai"):
    sys.stderr.write(
        "Warning: No index (.fai) found for {}. It's recommended to index the "
        "FASTA (using samtools faidx) for better performance.\n".format(
            genome_fasta
        )
    )

# Process each BED file independently
for bed_path in bed_files:
    # Prepare output file names
    base_name = os.path.splitext(os.path.basename(bed_path))[0]
    out_file1 = "{}_CpG-per-base.bed".format(base_name)
    out_file2 = "{}_CpG-per-GC.bed".format(base_name)
    out_file3 = "{}_CpG-OE.bed".format(base_name)

    # Load BED intervals
    bed = BedTool(bed_path)

    # Use bedtools nuc via pybedtools to get nucleotide content and CpG count
    # -pattern "CG" counts CpG occurrences; -C makes the search case-insensitive
    nuc_result = bed.nucleotide_content(
        fi=genome_fasta,
        pattern="CG",
        C=True,
    )

    # Open output files for writing
    with open(out_file1, "w") as fout1, open(
        out_file2, "w"
    ) as fout2, open(out_file3, "w") as fout3:
        # Each line in nuc_result corresponds to one interval. bedtools nuc
        # outputs a header as the first line, so skip it if present.
        with open(nuc_result.fn) as nuc_file:
            for i, line in enumerate(nuc_file):
                if i == 0:
                    # Detect a bedtools nuc header beginning with "chrom"
                    columns = line.strip().split()
                    if columns[0].lower() == "chrom" or columns[0].lower().startswith(
                        "chromosome"
                    ):
                        continue

                # Process a data line
                fields = line.strip().split()
                if len(fields) < 12:
                    # Skip lines without the expected nucleotide-content columns
                    continue

                # Original BED fields (assuming at least chr, start, end)
                chrom = fields[0]
                start = fields[1]
                end = fields[2]

                # bedtools nuc fields for a three-column BED input:
                # fields[5] = number of A bases
                # fields[6] = number of C bases
                # fields[7] = number of G bases
                # fields[8] = number of T bases
                # fields[9] = number of N bases
                # fields[11] = interval length
                # fields[12] = pattern occurrences when -pattern is used
                # These indices shift if the input BED has extra columns.
                try:
                    c_count = int(fields[6])
                    g_count = int(fields[7])
                    length = int(fields[11])
                except ValueError:
                    # Skip lines with an unexpected bedtools nuc format
                    continue

                # With -pattern and without -seq, the pattern count is expected
                # in the last column. Fall back to the preceding column if needed.
                cpg_count = 0
                if fields[-1].isdigit():
                    cpg_count = int(fields[-1])
                else:
                    try:
                        cpg_count = int(fields[-2])
                    except Exception:
                        cpg_count = 0

                # Calculate CpG per base
                if length > 0:
                    cpg_per_base = cpg_count / length
                else:
                    cpg_per_base = 0.0

                # Calculate CpG per C+G base
                if (c_count + g_count) > 0:
                    cpg_per_gc = cpg_count / (c_count + g_count)
                else:
                    # No C or G in the sequence: the metric is not applicable
                    cpg_per_gc = None

                # Calculate the observed/expected CpG ratio
                if c_count > 0 and g_count > 0:
                    cpg_oe = (cpg_count * float(length)) / (c_count * g_count)
                else:
                    # If either count is zero, the expected CpG value is zero
                    cpg_oe = None

                # Format values to six decimal places; undefined values are NA
                val1 = "{:.6f}".format(cpg_per_base)
                val2 = (
                    "{:.6f}".format(cpg_per_gc)
                    if cpg_per_gc is not None
                    else "NA"
                )
                val3 = "{:.6f}".format(cpg_oe) if cpg_oe is not None else "NA"

                # Write chromosome, start, end, and metric value
                fout1.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, val1))
                fout2.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, val2))
                fout3.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, val3))

    print(
        "Processed {}: outputs -> {}, {}, {}".format(
            bed_path,
            out_file1,
            out_file2,
            out_file3,
        )
    )
