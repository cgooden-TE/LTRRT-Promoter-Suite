# Inferential replicates for `EM_Reads`

Status: **design only, nothing implemented.** Written 2026-09-22 against
`isoclassifier_nested_em` at the Stage C commit.

## Why

`_em_units.tsv` reports `EM_Reads` as a point estimate, and DECLTR hands it to DESeq2 as though it
were an observed count. Two rows carrying 1000 reads are treated identically whether that 1000 was
counted or apportioned.

The reliability is not uniform, and it is worst exactly where the biology is interesting. Stable
rows are well-anchored genes; contested rows are nested TEs, TEs in genes, and overlapping gene
pairs. On chr1, 44% of units in contested clusters have `Unique_Reads` at or above 95% of their
`EM_Reads` (ordinary counting noise), and the other 56% are apportioned to some degree.

The sharpest argument is this repository's own history. The same reads, under four versions of the
model, gave `Zm00001eb049250` an `EM_Reads` of 108.6 and then 1,426.9, and `TE_homo_137237` 1,487
and then 169. Identical data. A point estimate with no spread invites conclusions that will not
survive the next methodological revision.

## What a replicate is

Resample the reads of a cluster with replacement, re-run the EM on the resampled set, and record
each unit's `EM_Reads`. Repeat B times. The spread across replicates is that unit's inferential
uncertainty. This is what Salmon's `--numBootstraps` does; `--numGibbsSamples` is the posterior
sampling alternative.

The resampling has to be over **reads**, not over per-read posteriors. Sampling each read's
assignment from its stored posterior is much cheaper but holds the fitted parameters fixed, so it
captures only which-read-went-where and understates the real spread.

Full resampling captures the thing that matters most here: **peak-calling instability propagates**.
A TSS peak that barely clears the Poisson test or `tss_trunc_ratio` is called in some replicates and
not others, and the unit's count swings with it. Since a peak call is the most consequential single
event in this model, that is precisely the uncertainty worth surfacing.

## What it does not do

Bootstrapping measures **precision, not accuracy**. A confidently wrong assignment is wrong in every
replicate and comes back with a tight interval.

This is not hypothetical. The antisense calls the junction audit proves wrong have a median
posterior of 1.000; they would bootstrap beautifully tight. Anyone reading a narrow interval as
validation would be misled in the place this pipeline is weakest.

So the two measurements are complementary and neither substitutes for the other:

- `junction_audit.py` says whether an assignment is **right**, on the ~50% of contested reads
  carrying junctions;
- bootstrap replicates say how **stable** a count is, on all of them.

Any report that carries one should say which.

## Cost, and where it actually goes

Bootstrap cost scales with **reads**, not clusters. This kills the obvious shortcut: on chr2, 61% of
the contested read mass sits in 10 of 943 clusters, so restricting the bootstrap to the heavy
clusters costs roughly as much as bootstrapping everything, and skipping the long tail (69% of
chr1's clusters) saves about 10% of the compute.

Concentration, for reference:

| clusters covering | chr1 (1,134 clusters) | chr2 (943 clusters) |
|---|---|---|
| 50% of read mass | 46 | 4 |
| 75% | 153 | 41 |
| 90% | 350 | 171 |
| 99% | 749 | 529 |

The savings therefore have to come from elsewhere, in this order:

1. **Re-run only `resolve()`, in process.** The BAM pass, GFF parsing and unit construction are
   per-run costs that a bootstrap has no reason to repeat. This is the single biggest lever and it
   is the whole reason for implementing this inside IsoClassifier rather than as a wrapper script.
2. **B = 20-30, not 100.** Enough for a coefficient of variation; linear in cost.
3. **m-out-of-n for the heavy clusters.** For a 218,000-read cluster, resample m << n reads per
   replicate and scale the variance by the subsample fraction. Standard technique, known caveats
   (the scaling assumes the estimator is asymptotically normal in n, which is worth checking against
   one fully bootstrapped cluster).
4. **Skip clusters with no contest.** Units where `Unique_Reads` ~ `EM_Reads` need no replicates.
5. **Parallelise by cluster**, reusing the existing `--threads` machinery unchanged. Clusters are
   connected components, so replicates are embarrassingly parallel and need no new concurrency code.

Order of magnitude on chr1: the EM stage is ~250 s of a 313 s run, so 20 in-process replicates is
roughly 80 minutes per chromosome before any of 3-5 are applied.

## First-pass implementation

Deliberately minimal: no m-out-of-n, no variance model, no sampling design. Get honest numbers for
contested clusters first, then decide what to approximate.

### 1. Flag

```
--em-bootstrap N        number of inferential replicates (default 0, off)
--em-bootstrap-seed S   base seed, so a run is reproducible (default 0)
```

### 2. Where it hooks in

`em_assign_deferred()` (IsoClassifier.py) already builds `records` and calls `nested_em.resolve()`
once. The bootstrap loop goes around that call, reusing the same `records` list and the same
`builder` closure, so nothing is re-read or re-parsed:

```python
point = nested_em.resolve(records, builder, p, pool=pool)
reps = []
for b in range(args.em_bootstrap):
    reps.append(nested_em.resolve(records, builder, p, pool=pool,
                                  resample_seed=args.em_bootstrap_seed + b,
                                  log=lambda *_: None))
```

The point estimate stays exactly what it is today. Replicates are additional, never a replacement,
so existing output is unchanged when the flag is off.

### 3. What `resolve()` needs

One new argument. When `resample_seed` is not None, within each cluster draw `len(ridx)` read
indices with replacement from `ridx` before the passes begin, and run the existing machinery on that
multiset. Everything downstream -- peak calling, the Kaplan-Meier fits, `em_cluster` -- then
operates on resampled reads with no further change, which is what makes peak instability show up in
the spread.

Two details that are easy to get wrong:

- **Resample within cluster, not globally.** Clusters are independent; a global resample would
  change cluster membership and make replicates non-comparable.
- **Anchors resample too.** They are what pin a unit's shape, so holding them fixed would
  understate the spread. They are in `records` already, so this is automatic as long as the
  resample is over all of a cluster's record indices rather than only the deferred ones.

Duplicated reads are handled naturally: a read drawn twice contributes twice to the M-step sums. The
leave-one-out subtractions stay per-record and remain correct.

### 4. Output

A fourth table, `<output>_em_units_boot.tsv`, one row per unit:

| column | meaning |
|---|---|
| `Feature`, `Orientation` | join key to `_em_units.tsv` |
| `EM_Reads` | the point estimate, repeated for convenience |
| `Boot_Mean`, `Boot_SD` | across replicates |
| `Boot_CV` | `Boot_SD / max(Boot_Mean, 1)` -- the column to sort by |
| `Boot_P05`, `Boot_P95` | 5th and 95th percentiles |
| `Boot_Zero_Frac` | share of replicates giving this unit under 1 read; the "does it exist at all" number |
| `N_Replicates` | B actually completed for this unit |

`Boot_Zero_Frac` is the one to watch for phantoms: a unit whose reads depend on a marginal peak will
lose them entirely in some replicates.

Raw per-replicate counts go to `<output>_em_boot_raw.tsv` (`Feature`, `Orientation`, `Replicate`,
`EM_Reads`) so downstream tools that want the full distribution -- swish/fishpond expect exactly
this shape -- can have it without re-running anything.

### 5. Validation before trusting it

- **B-stability.** Run B = 10, 20, 40 on one chromosome and check that `Boot_SD` has converged. If
  it has not, the whole design needs more replicates, not fewer.
- **Anchored units are a null control.** A unit with `Unique_Reads` ~ `EM_Reads` should come back
  with `Boot_CV` near the Poisson expectation, ~`1/sqrt(n)`. If those units show a large spread,
  the resampling is wrong, not the biology.
- **Known-contested units are a positive control.** `Zm00001eb049250` / `TE_homo_137237` on chr1
  and the cluster-721 gene pair on chr2 should show visibly larger spread than their neighbours.
- **Do not compare `Boot_SD` against the junction audit** and conclude anything about accuracy.
  They measure different things; see the section above.

### 6. What to add only after the first pass runs

- m-out-of-n for heavy clusters, once the cost is actually measured rather than estimated.
- A variance model (`Boot_CV` regressed on `Unique_Reads/EM_Reads`, posterior spread,
  `Cluster_Units`, `TSS_Peaks_In_Zone`, `Promoter_Anchor`) to impute uncertainty for units that were
  not bootstrapped, in the spirit of a DESeq2 dispersion trend. This is what lets a downstream tool
  have a value for every feature.
- Gibbs sampling as an alternative to resampling, if the bootstrap turns out to be poorly behaved
  for units with very few reads.

## Open questions

- **Does the annotation belong in the resample?** Peak instability is captured; annotation error is
  not. `Promoter_Anchor = annotated` units (32-37% of gene sense units) carry a systematic
  uncertainty that no read-level bootstrap will show. Perturbing the promoter window would capture
  it, and that is a much bigger design.
- **What does DESeq2 do with this?** Passing replicates through DECLTR is not free; swish is built
  for it, DESeq2 is not. Deciding the downstream path should probably come before implementing the
  upstream one.
