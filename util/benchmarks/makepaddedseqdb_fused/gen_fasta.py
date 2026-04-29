#!/usr/bin/env python3
"""Generate a synthetic nucleotide FASTA with long sequences.

Sequences are long enough (default 50 KB) to trigger splitting with
the default --max-seq-len 10000.  Used by the benchmark/correctness suite.
"""

import argparse
import random
import sys


def parse_size(s):
    """Parse a human-readable size string (e.g. '100M', '5G') to bytes."""
    s = s.strip().upper()
    multipliers = {'B': 1, 'K': 1024, 'M': 1024**2, 'G': 1024**3, 'T': 1024**4}
    if s[-1] in multipliers:
        return int(float(s[:-1]) * multipliers[s[-1]])
    return int(s)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', default='100M',
                        help='Target total FASTA size (e.g. 100M, 1G)')
    parser.add_argument('--seq-len', type=int, default=50000,
                        help='Length of each sequence in bases (default: 50000)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility')
    parser.add_argument('-o', '--output', default='/dev/stdout',
                        help='Output file (default: stdout)')
    args = parser.parse_args()

    target_bytes = parse_size(args.size)
    seq_len = args.seq_len
    rng = random.Random(args.seed)
    bases = 'ACGT'
    line_width = 80

    written = 0
    seq_id = 0

    out = open(args.output, 'w') if args.output != '/dev/stdout' else sys.stdout
    try:
        while written < target_bytes:
            header = f'>synth_{seq_id} len={seq_len}\n'
            out.write(header)
            written += len(header)

            remaining = seq_len
            while remaining > 0:
                chunk = min(remaining, line_width)
                line = ''.join(rng.choice(bases) for _ in range(chunk))
                out.write(line)
                out.write('\n')
                written += chunk + 1
                remaining -= chunk

            seq_id += 1
    finally:
        if out is not sys.stdout:
            out.close()

    print(f'Generated {seq_id} sequences, ~{written} bytes', file=sys.stderr)


if __name__ == '__main__':
    main()
