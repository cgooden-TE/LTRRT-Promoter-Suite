#!/usr/bin/env python3
"""Score IsoClassifier per-read assignments against the origins encoded in simulated read names."""
import sys
from collections import Counter, defaultdict

import pandas as pd


def load(prefix):
    parts = []
    for suf in ('_3p_softclip_per_read_spliced.tsv', '_3p_softclip_per_read_nonspliced.tsv'):
        parts.append(pd.read_csv(prefix + suf, sep='\t', dtype=str))
    df = pd.concat(parts, ignore_index=True)
    truth = df['Read'].str.split('|', expand=True)
    df['scen'], df['true_feat'], df['true_orient'], df['kind'] = truth[0], truth[1], truth[2], truth[4 - 1]
    df['ambiguous'] = df['Nested_In'].fillna('').ne('') | (
        df.get('Assign_Method', pd.Series('', index=df.index)).fillna('') == 'em')
    df['correct'] = (df['Feature'] == df['true_feat']) & (df['Orientation'] == df['true_orient'])
    return df


def report(prefix, label):
    df = load(prefix)
    print(f'== {label} ==')
    amb = df[df['ambiguous']]
    print(f'reads with a 5\' end inside >1 feature: {len(amb)}; '
          f'correctly attributed: {amb["correct"].sum()} ({100 * amb["correct"].mean():.1f}%)')
    for scen, g in amb.groupby('scen'):
        print(f'  scenario {scen}: {g["correct"].sum()}/{len(g)} correct')
    # per-unit totals: assigned vs true, restricted to reads that were assignable at all
    assigned = Counter(zip(df['Feature'], df['Orientation']))
    true = Counter(zip(df['true_feat'], df['true_orient']))
    keys = sorted(set(assigned) | set(true))
    em = {}
    try:
        u = pd.read_csv(prefix + '_em_units.tsv', sep='\t')
        for f, o, n, uq in zip(u['Feature'], u['Orientation'], u['EM_Reads'], u['Unique_Reads']):
            em[(f, o)] = n
    except FileNotFoundError:
        pass
    print(f'  {"unit":34s} {"true":>6s} {"assigned":>9s}' + (f' {"EM_expected":>12s}' if em else ''))
    for k in keys:
        extra = ''
        if em:
            extra = f' {em[k]:12.1f}' if k in em else f' {"(unique)":>12s}'
        print(f'  {k[0] + ":" + k[1]:34s} {true.get(k, 0):6d} {assigned.get(k, 0):9d}' + extra)
    return df


if __name__ == '__main__':
    for arg in sys.argv[1:]:
        prefix, _, label = arg.partition('=')
        report(prefix, label or prefix)
