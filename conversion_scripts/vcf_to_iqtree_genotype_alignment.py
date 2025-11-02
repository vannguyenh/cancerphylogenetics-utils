#!/usr/bin/env python3
"""
vcf_to_iqtree_genotype_alignment.py

Convert a VCF into a SNP alignment using IQ-TREE’s genotype symbols
(based on the provided conversion table). Supports phased (|) and
unphased (/) GT encodings and writes FASTA and/or PHYLIP (relaxed).

Assumptions/Rules
- Keeps ONLY biallelic SNPs (REF and single ALT, both length 1).
- Uses the IQ-TREE genotype mapping from the user's table:
  * Phased: A|C->M, C|A->!, T|G->&, etc.
  * Unphased: A/C->M, C/T->Y, G/T->K, etc.
  * Homozygotes: A/A or A|A -> 'A' (likewise for C,G,T).
- Missing/unsupported -> '?'.
- PHYLIP output uses relaxed names (no strict 10-char truncation).

Usage examples
--------------
# Both outputs (default), autodetect gz by extension:
python vcf_to_iqtree_genotype_alignment.py -i input.vcf.gz -o outprefix

# FASTA only
python vcf_to_iqtree_genotype_alignment.py -i input.vcf -o out --fasta-only

# PHYLIP only
python vcf_to_iqtree_genotype_alignment.py -i input.vcf -o out --phylip-only

"""

import sys
import os
import argparse
import gzip
from collections import defaultdict

PHASED_MAP = {
    "A|A":"A","C|C":"C","G|G":"G","T|T":"T",
    "A|C":"M","A|G":"R","A|T":"W","C|G":"S","C|T":"Y","G|T":"K",
    "C|A":"!","G|A":'"',"T|A":"#","G|C":"$","T|C":"%","T|G":"&",
}

UNPHASED_MAP = {
    "A/A":"A","C/C":"C","G/G":"G","T/T":"T",
    "A/C":"M","C/A":"M",
    "A/G":"R","G/A":"R",
    "A/T":"W","T/A":"W",
    "C/G":"S","G/C":"S",
    "C/T":"Y","T/C":"Y",
    "G/T":"K","T/G":"K",
}

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def convert_gt_to_symbol(gt, ref, alt, missing_char="?"):
    """Map a genotype string to the IQ-TREE symbol according to the table."""
    if gt in (".", "./.", ".|.", None):
        return missing_char

    # Decide phased vs unphased
    if "|" in gt:
        alle_idx = gt.split("|")
        phased = True
    elif "/" in gt:
        alle_idx = gt.split("/")
        phased = False
    else:
        # Haploid or single allele e.g. "0"
        alle_idx = [gt]
        phased = False

    # Translate allele indices (0->REF, 1->ALT); others unsupported
    alle_bases = []
    for a in alle_idx:
        if a == ".":
            alle_bases.append(missing_char)
        else:
            try:
                ai = int(a)
            except ValueError:
                return missing_char
            if ai == 0:
                alle_bases.append(ref.upper())
            elif ai == 1:
                alle_bases.append(alt.upper())
            else:
                return missing_char  # multiallelic not supported by table

    if any(b == missing_char for b in alle_bases):
        return missing_char

    if len(alle_bases) == 1:
        b = alle_bases[0]
        return b if b in ("A","C","G","T") else missing_char

    a1, a2 = alle_bases[0], alle_bases[1]
    key = f"{a1}|{a2}" if phased else f"{a1}/{a2}"
    return (PHASED_MAP if phased else UNPHASED_MAP).get(key, missing_char)

def write_fasta(out_path, names, seqs, wrap=80):
    with open(out_path, "w") as fh:
        for n in names:
            s = seqs[n]
            fh.write(f">{n}\n")
            if wrap and wrap > 0:
                for i in range(0, len(s), wrap):
                    fh.write(s[i:i+wrap] + "\n")
            else:
                fh.write(s + "\n")

def write_phylip_relaxed(out_path, names, seqs):
    with open(out_path, "w") as fh:
        fh.write(f"{len(names)} {len(seqs[names[0]]) if names else 0}\n")
        for n in names:
            fh.write(f"{n} {seqs[n]}\n")

def main():
    ap = argparse.ArgumentParser(description="VCF -> IQ-TREE genotype alignment (FASTA/PHYLIP).")
    ap.add_argument("-i", "--vcf", required=True, help="Input VCF (optionally .gz)")
    ap.add_argument("-o", "--outprefix", required=True, help="Output prefix (writes .fasta/.phy)")
    ap.add_argument("--missing-char", default="?", help="Character for missing/unsupported [default: ?]")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--fasta-only", action="store_true", help="Write only FASTA")
    g.add_argument("--phylip-only", action="store_true", help="Write only PHYLIP (relaxed)")
    args = ap.parse_args()

    fasta_out = args.outprefix + ".fasta"
    phy_out   = args.outprefix + ".phy"

    with open_maybe_gzip(args.vcf) as f:
        sample_names = []
        seq_lists = defaultdict(list)
        total_variants = 0
        kept_sites = 0

        for line in f:
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                sample_names = parts[9:]
                for s in sample_names:
                    seq_lists[s] = []
                continue

            total_variants += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            fmt = parts[8]
            gts = parts[9:]

            # Keep only biallelic SNPs
            if "," in alt:
                continue
            if len(ref) != 1 or len(alt) != 1:
                continue

            fmt_keys = fmt.split(":")
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                gt_idx = None

            kept_sites += 1
            for si, s in enumerate(sample_names):
                if gt_idx is None or si >= len(gts):
                    sym = args.missing_char
                else:
                    fields = gts[si].split(":")
                    gt_field = fields[gt_idx] if gt_idx < len(fields) else "."
                    sym = convert_gt_to_symbol(gt_field, ref, alt, missing_char=args.missing_char)
                seq_lists[s].append(sym)

    # Build sequences
    seqs = {s: "".join(seq_lists[s]) for s in sample_names}

    # Write outputs
    if not args.phylip_only:
        write_fasta(fasta_out, sample_names, seqs)
    if not args.fasta_only:
        write_phylip_relaxed(phy_out, sample_names, seqs)

    # Summary
    aln_len = len(seqs[sample_names[0]]) if sample_names else 0
    sys.stderr.write(
        f"[vcf_to_iqtree] Samples: {len(sample_names)} | Sites kept: {kept_sites} | "
        f"Alignment length: {aln_len}\n"
        f"[vcf_to_iqtree] Wrote: "
        f"{fasta_out if not args.phylip_only else '(FASTA skipped)'} "
        f"{phy_out if not args.fasta_only else '(PHYLIP skipped)'}\n"
    )

if __name__ == "__main__":
    main()
