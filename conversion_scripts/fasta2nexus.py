#!/usr/bin/env python3
"""
Convert an aligned FASTA to NEXUS (DNA).
- Preserves names; wraps in single quotes if needed.
- Converts U->T; leaves gaps '-' and missing '?' alone.
Usage:
  python fasta2nexus.py -i input.fasta -o output.nex
"""
import argparse, re
from pathlib import Path

def parse_fasta(fp: Path):
    if not fp.exists():
        raise FileNotFoundError(fp)
    recs = []
    name = None
    seq_chunks = []
    with fp.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seq = ''.join(seq_chunks)
                    recs.append((name, seq))
                name = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            seq = ''.join(seq_chunks)
            recs.append((name, seq))
    if not recs:
        raise ValueError('No sequences found')
    lengths = {len(s) for _, s in recs}
    if len(lengths) != 1:
        raise ValueError(f'Sequences not equal length: {sorted(lengths)}')
    return recs

def quote_taxon(t):
    if re.search(r"[^A-Za-z0-9_.-]", t):
        t2 = t.replace("'", "''")
        return "'" + t2 + "'"
    return t

def write_nexus(records, outpath: Path):
    ntax = len(records)
    nchar = len(records[0][1])
    lines = []
    lines.append('#NEXUS')
    lines.append('Begin data;')
    lines.append(f'    Dimensions ntax={ntax} nchar={nchar};')
    lines.append('    Format datatype=DNA missing=? gap=-;')
    lines.append('    Matrix')
    for name, seq in records:
        lines.append(f'    {quote_taxon(name):<30} {seq}')
    lines.append('    ;')
    lines.append('End;')
    outpath.write_text('\n'.join(lines))

def main():
    ap = argparse.ArgumentParser(description='Convert aligned FASTA to NEXUS (DNA).')
    ap.add_argument('-i','--input', required=True, help='Input FASTA')
    ap.add_argument('-o','--output', required=True, help='Output NEXUS')
    args = ap.parse_args()
    records = parse_fasta(Path(args.input))
    write_nexus(records, Path(args.output))
    print(f'Wrote {args.output} (ntax={len(records)}, nchar={len(records[0][1])})')

if __name__ == '__main__':
    main()
