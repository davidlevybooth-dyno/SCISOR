# Constrained diffusion: domain locking (WS2)

Goal: let the user declare regions that **must not be deleted** (catalytic residues,
ligand motifs, structured domains) and guarantee SCISOR never removes them, while
still hitting the target length. This is the first and cheapest form of "constrained
diffusion" — it changes only the sampling step, not the model.

## 1. Where it hooks in

SCISOR already zeroes deletion probability at special tokens inside `p_sample`. From
[SCISOR/SCISOR/shortening_scud.py](../SCISOR/shortening_scud.py) (lines ~127–132):

```python
model_probs = (
    model_probs
    * (x != self.pad_token_id)
    * (x != self.tokenizer.cls_token_id)
    * (x != self.tokenizer.eos_token_id)
)
```

A keep-mask is the same operation, applied to user-specified positions. Proposed
addition (hard lock):

```python
if getattr(self, "keep_mask", None) is not None:
    # keep_mask: (B, L) bool/0-1, 1 = protected (never delete)
    model_probs = model_probs * (1 - self.keep_mask.to(model_probs.dtype))
```

Because deletions happen one residue per step and `p_sample` re-reads `model_probs`
each step, the mask must be tracked in the **shrinking coordinate frame**: after each
deletion, drop the deleted column from the mask so it stays aligned with `x`. This is
the main implementation subtlety (`shrink_sequence` already maintains
`preserved_indices`; the mask follows the same bookkeeping).

### Feasibility guard
If `(#deletable positions) < num_deletions`, the target length is infeasible under the
locks. The runner must detect this up front and report it (e.g. "TSC2 → 1274 needs 533
deletions but only 410 positions are unlocked"), rather than silently under-deleting.

## 2. Constraint variants

- **Hard lock (default):** probability set to 0 at protected positions. Guarantees
  zero deletions there.
- **Soft penalty:** multiply protected-position probability by a small factor (or add a
  log-prob penalty) instead of zero — discourages but allows, useful for "prefer to
  keep" regions.
- **Oversample-and-filter (no code change):** draw extra samples and reject any that
  touched a frozen position. Cheap to prototype, wasteful at high lock fractions;
  useful as a correctness check against the masked implementation.

## 3. Region specification schema

Per-target JSON, resolved to a 0/1 mask over WT coordinates. Sources: UniProt features
(domains, active sites, binding sites) plus manual entries for motifs SCISOR/UniProt
won't flag.

```json
{
  "target": "SYNGAP1",
  "uniprot": "Q96PV0",
  "length": 1343,
  "keep": [
    {"name": "C2", "start": 150, "end": 261, "reason": "domain"},
    {"name": "RasGAP", "start": 362, "end": 553, "reason": "catalytic domain"},
    {"name": "catalytic_R", "pos": [485], "reason": "GAP arginine finger"},
    {"name": "coiled_coil", "start": 1189, "end": 1262, "reason": "trimerization"},
    {"name": "PDZ_ligand_QTRV", "start": 1340, "end": 1343, "reason": "PSD-95 binding"}
  ]
}
```
(Indices illustrative — to be finalized per target with the team; 0- vs 1-based and
inclusive/exclusive conventions documented in the schema header.)

Frozen sets to encode (from the rational playbook):
- **SYNGAP1:** C2 / RasGAP / coiled-coil / C-terminal -QTRV.
- **SHANK3:** SH3 / PDZ / SAM / ANK.
- **TSC2:** HEAT / GAP (and the dimerization interfaces retained in MT9).

This schema is shared with the evaluation harness (motif retention) and with WS3
(block deletion operates on the *complement* of the keep set).

## 4. Runner integration

Add `--keep-mask configs/frozen/<TARGET>.json` to [scisor_shrink.py](../scisor_shrink.py):
- Build the 0/1 mask in WT coordinates, expand per sample, attach as `model.keep_mask`.
- Run the feasibility guard before sampling.
- Record the active constraint set in the output FASTA header for auditability.

```sh
python scisor_shrink.py --input ../targets/SYNGAP1.fasta \
    --output ../results/aav/SYNGAP1_aav1274_locked_n32.fasta \
    --target-length 1274 --num-samples 32 --batch-size 16 \
    --keep-mask configs/frozen/SYNGAP1.json
```

## 5. Milestones / outcomes
- **M2.1:** keep-mask in `p_sample` with coordinate-tracked bookkeeping; unit test that
  no frozen position is ever deleted (cross-checked against oversample-and-filter).
- **M2.2:** frozen-set schema + `--keep-mask` flag; constraint sets for the 3 targets.
- **M2.3:** constrained 32-sample baselines; report the naturalness/score cost of locks
  vs unconstrained, and confirm target length still met (or feasibility flagged).

## 6. Limitations (motivates WS3)
Locking prevents deleting the wrong residues, but scattered single-residue deletions
*between* locked domains can still disorganize structure and leave domains tethered by
nothing. Preserving tertiary organization — contiguous block deletion plus engineered
linkers — is [structural-integrity-linkers.md](structural-integrity-linkers.md).
