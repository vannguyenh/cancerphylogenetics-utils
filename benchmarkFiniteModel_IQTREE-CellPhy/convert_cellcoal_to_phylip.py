#!/usr/bin/env python3
"""
Convert CellCoal snv_gen output to GT16 PHYLIP format for IQ-TREE or CellPhy.

IQ-TREE and CellPhy use DIFFERENT encodings for GT16 phased heterozygotes:

  Genotype    IQ-TREE   CellPhy
  --------    -------   -------
  A|C           M         1
  A|G           R         2
  A|T           W         3
  C|G           S         4
  C|T           Y         5
  G|T           K         6
  C|A           !         !
  G|A           "         "
  T|A           #         #
  G|C           $         $
  T|C           %         %
  T|G           &         &

For GT10 (unphased), both tools use the same IUPAC encoding (snv_hap files).

USAGE:
  # Convert replicate 1 SNV-only for IQ-TREE GT16 (default):
  python3 convert_cellcoal_to_phylip.py --input_dir results_S1_ADO0.00_ERR0.00 --replicate 1 --tool iqtree

  # Convert replicate 1 full genotypes for IQ-TREE GT16:
  python3 convert_cellcoal_to_phylip.py --input_dir results_S1_ADO0.00_ERR0.00 --replicate 1 --tool iqtree --sites full

  # Convert all 20 replicates, both tools, both site types:
  for i in $(seq 1 20); do
      for tool in iqtree cellphy; do
          for sites in snv full; do
              python3 convert_cellcoal_to_phylip.py --input_dir results_S1_ADO0.00_ERR0.00 --replicate $i --tool $tool --sites $sites
          done
      done
  done
"""

import argparse
import os
import sys

# ============================================================
# GT16 encoding: phased diploid genotypes -> single character
# In CellCoal snv_gen files, each genotype is a two-letter pair
# where first letter = maternal allele, second = paternal allele.
# ============================================================

# IQ-TREE GT16 encoding: IUPAC letters for standard-order heterozygotes
GT16_MAP_IQTREE = {
    'AA': 'A', 'CC': 'C', 'GG': 'G', 'TT': 'T',
    'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K',
    'CA': '!', 'GA': '"', 'TA': '#', 'GC': '$', 'TC': '%', 'TG': '&',
}

# CellPhy GT16 encoding: numeric 1-6 for standard-order heterozygotes
GT16_MAP_CELLPHY = {
    'AA': 'A', 'CC': 'C', 'GG': 'G', 'TT': 'T',
    'AC': '1', 'AG': '2', 'AT': '3', 'CG': '4', 'CT': '5', 'GT': '6',
    'CA': '!', 'GA': '"', 'TA': '#', 'GC': '$', 'TC': '%', 'TG': '&',
}


def encode_gt16(genotype_pair, gt16_map):
    """Encode a two-letter CellCoal genotype to a GT16 single char."""
    if '?' in genotype_pair or 'N' in genotype_pair or '-' in genotype_pair:
        return '?'
    geno = genotype_pair.upper()
    if geno in gt16_map:
        return gt16_map[geno]
    else:
        print(f"WARNING: Unknown genotype '{genotype_pair}', encoding as '?'",
              file=sys.stderr)
        return '?'


def parse_genotype_file(filepath, has_site_indices=True):
    """
    Parse CellCoal genotype file (snv_gen or full_gen).
    Format: first line = num_taxa num_sites
            second line = site indices (snv_gen only)
            remaining lines = taxon_name followed by space-separated 2-char genotypes

    Args:
        filepath: path to the genotype file
        has_site_indices: True for snv_gen (has site index line), False for full_gen

    Returns: (num_taxa, num_sites, taxa_dict)
             where taxa_dict = {taxon_name: [genotype_pairs]}
    """
    taxa = {}
    with open(filepath, 'r') as f:
        header = f.readline().strip().split()
        num_taxa = int(header[0])
        num_sites = int(header[1])

        # snv_gen has a second line with site indices; full_gen does not
        if has_site_indices:
            _site_indices = f.readline().strip()

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            taxon_name = parts[0]
            genotypes = parts[1:]
            taxa[taxon_name] = genotypes

    return num_taxa, num_sites, taxa


def write_phylip(filepath, taxa_sequences, num_sites):
    """Write PHYLIP format file."""
    num_taxa = len(taxa_sequences)
    with open(filepath, 'w') as f:
        f.write(f"{num_taxa} {num_sites}\n")
        for taxon_name, sequence in taxa_sequences:
            f.write(f"{taxon_name:<20s}{sequence}\n")


def convert_gt16(input_dir, replicate, output_dir, tool, sites):
    """
    Convert CellCoal genotype file to GT16 PHYLIP format.

    Args:
        input_dir: CellCoal results directory
        replicate: replicate number (1-based)
        output_dir: where to write output
        tool: 'iqtree' or 'cellphy' (determines encoding)
        sites: 'snv' (variable sites only) or 'full' (all sites)
    """
    gt16_map = GT16_MAP_IQTREE if tool == 'iqtree' else GT16_MAP_CELLPHY

    if sites == 'snv':
        gen_file = os.path.join(input_dir, 'snv_genotypes_dir',
                                f'snv_gen.{replicate:04d}')
        has_site_indices = True
    else:
        gen_file = os.path.join(input_dir, 'full_genotypes_dir',
                                f'full_gen.{replicate:04d}')
        has_site_indices = False

    if not os.path.exists(gen_file):
        print(f"ERROR: File not found: {gen_file}", file=sys.stderr)
        return None

    num_taxa, num_sites, taxa = parse_genotype_file(gen_file, has_site_indices)

    taxa_sequences = []
    for taxon_name in sorted(taxa.keys()):
        genotypes = taxa[taxon_name]
        encoded = ''.join(encode_gt16(g, gt16_map) for g in genotypes)
        taxa_sequences.append((taxon_name, encoded))

    out_file = os.path.join(output_dir, f'rep{replicate:04d}_GT16_{tool}_{sites}.phy')
    write_phylip(out_file, taxa_sequences, num_sites)
    print(f"GT16 ({tool}, {sites}): {out_file} ({num_taxa} taxa, {num_sites} sites)")
    return out_file


def main():
    parser = argparse.ArgumentParser(
        description='Convert CellCoal snv_gen to GT16 PHYLIP format.\n'
                    'Supports both IQ-TREE and CellPhy GT16 encodings.\n'
                    'NOTE: For GT10, use snv_hap files directly '
                    '(IUPAC-encoded with -x option, works for both tools).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input_dir', required=True,
                        help='CellCoal results directory')
    parser.add_argument('--replicate', type=int, required=True,
                        help='Replicate number (1-based)')
    parser.add_argument('--tool', choices=['iqtree', 'cellphy'], default='iqtree',
                        help='Target tool encoding (default: iqtree). '
                             'iqtree: M,R,W,S,Y,K for het. '
                             'cellphy: 1,2,3,4,5,6 for het.')
    parser.add_argument('--sites', choices=['snv', 'full'], default='snv',
                        help='Site type (default: snv). '
                             'snv: variable sites only (from snv_genotypes_dir). '
                             'full: all sites (from full_genotypes_dir).')
    parser.add_argument('--output_dir', default=None,
                        help='Output directory (default: input_dir/phylip_gt16_<tool>_<sites>)')

    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.path.join(args.input_dir, f'phylip_gt16_{args.tool}_{args.sites}')
    os.makedirs(args.output_dir, exist_ok=True)

    convert_gt16(args.input_dir, args.replicate, args.output_dir, args.tool, args.sites)


if __name__ == '__main__':
    main()