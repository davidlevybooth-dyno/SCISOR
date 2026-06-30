# Evaluation: scoring stack and rational-design comparison (WS0, WS1)

Goal: turn shrunk-sequence FASTAs into a ranked, auditable scorecard, and compare
SCISOR's deletions against hand/rational designs on the same axes. SCISOR optimizes
naturalness only; this is the layer that tells us whether a construct is actually
foldable, functional, and deliverable.

## 1. Scoring axes

| Axis | Metric | Tool | Question it answers |
|---|---|---|---|
| Structure (global) | mean pLDDT, pTM | ESMFold / Boltz | Does it fold confidently? |
| Structure (vs original) | templated RMSD / TM-score per retained domain | structure align (e.g. US-align/TM-align) | Are kept domains preserved? |
| Motif retention | presence + position of catalytic/ligand motifs | sequence/feature check | Did we keep the functional residues? |
| Interface (where relevant) | co-fold iPAE | Boltz / AF-multimer-style | Is the binding interface intact (e.g. SYNGAP1 + PSD-95, TSC2 + Hamartin)? |
| AAV-fit | CDS bp, fits ≤3825 bp / ≤1274 aa | arithmetic | Will it package? |
| Naturalness | SCISOR log-likelihood / per-position deletion prob | SCISOR | How natural is the shrunken sequence? |

The "dual-depressor" caveat applies when scoring linker-containing constructs:
disordered linkers (pLDDT ~25–34) and chimeric junctions both lower confidence
without indicating a broken design. Score domains and junctions separately; do not
penalize a construct for low linker pLDDT alone.

## 2. Comparators

Score these side-by-side per target:
- SCISOR baseline samples (current 32-sample sets at 1,274 aa).
- SCISOR constrained samples (WS2, once locks exist).
- Rational/manual designs (MT9 for TSC2; SHANK3 and SYNGAP1 minigenes).
- External baselines where available (Raygun, ProGen2 length-conditioned).

### Rational-design ingestion (input contract)
Rational designs are supplied later as FASTA. Each record:
```
>TSC2_MT9-H1|kind rational|target TSC2|notes (GGGGS)x3 junctions
M... (the designed sequence) ...
```
The harness aligns each rational design back to the wild-type target to recover its
implied deletion map (for overlaying against SCISOR's), then scores it identically.

## 3. Harness: `scisor_score.py`

A new script (sibling of [scisor_shrink.py](../scisor_shrink.py)) that:
1. Reads one or more FASTAs (SCISOR output headers carry the deletion map; rational
   designs are aligned to WT to derive theirs).
2. Runs structure prediction (pluggable backend: ESMFold local, or Boltz; A100 for
   throughput).
3. Computes the axes in section 1 into a tidy table (CSV/Parquet): one row per
   (target, source, sample), columns = metrics.
4. Emits per-residue **deletion-frequency maps** across the sample set (how often each
   WT position is deleted) for the tolerance overlay (M1.2).
5. Ranks constructs per target by a configurable weighted score; flags AAV-fit pass/fail.

Suggested layout:
```
scisor_score.py            # CLI
SCISOR/scoring/
  fold.py                  # ESMFold/Boltz backends -> structure + pLDDT/pTM
  align.py                 # WT alignment -> deletion map; templated RMSD/TM
  motifs.py                # frozen-set / motif retention checks (shares schema with WS2)
  aav.py                   # length -> bp, budget check
  report.py                # tables, deletion-frequency maps, ranking
```

CLI sketch:
```sh
python scisor_score.py \
    --inputs results/aav/TSC2_aav1274_n32.fasta rational/TSC2_MT9.fasta \
    --wt targets/TSC2.fasta \
    --fold-backend esmfold \
    --motifs configs/frozen/TSC2.json \
    --out results/scores/TSC2.csv
```

## 4. Deletion-tolerance map (M1.2)

Two complementary views per target:
- **Sample frequency map:** across the 32 baseline samples, fraction of samples that
  delete each WT residue. High-frequency = SCISOR considers it expendable.
- **M-sweep:** shrink by d% for d in {1,3,5,10,20,30,40}; track which positions start
  getting deleted as pressure increases (Fig. 3b style).

Overlay both on the domain map. Headline question (M1.3): does SCISOR's "keep" signal
light up the frozen sets (SYNGAP1 C2/RasGAP/CC/QTRV; SHANK3 SH3/PDZ/SAM/ANK; TSC2
HEAT/GAP) while its "delete" signal lands on the scaffold the rational designs already
cut (SYNGAP1 ~730–1188; SHANK3 proline-rich; TSC2 middle)? Agreement validates the
rational priors; disagreement is a concrete flag.

## 5. Milestones / outcomes
- **Harness v1** runs end-to-end on the three current baselines → ranked table.
- **Tolerance maps** for TSC2 / SHANK3 / SYNGAP1.
- **Baseline-vs-rational comparison report** with the agreement metric and the
  micro-dystrophin calibration check.

## 6. Open questions
- Folding backend of record (ESMFold for speed vs Boltz for interfaces) — likely both,
  ESMFold for screening and Boltz for interface/co-fold on shortlisted constructs.
- Where templated RMSD uses the experimental structure vs an AF2/ESMFold model of WT.
- Exact weighting for the composite rank (set per target with the team).
