# Performance: making SCISOR fast (WS-P)

A deliberately **separate improvement path** from quality. SCISOR is slow because it
deletes **one residue per model forward pass, serially**. This workstream cuts that cost
without changing what the model represents. Every change here is measured on the speed
**and** quality scoreboards in [benchmarking.md](benchmarking.md) against the frozen
baseline — a speedup that quietly degrades shrinks is not a win.

This track depends only on WS0 (the harness) and runs in parallel with the quality work
(WS2 locks). Landing it early makes every later folding/iteration loop — including WS3 —
cheaper.

## 1. The bottleneck (verified in code)

[../SCISOR/shortening_scud.py](../SCISOR/shortening_scud.py), `shrink_sequence`:

```python
while S.sum() > 0:
    dels_this_interval = (S > 0).int()   # exactly ONE deletion per forward pass
    x, deleted_indices = self.p_sample(x, None, None, S, dels_this_interval, ...)
    ...
    S = S - dels_this_interval
```

`dels_this_interval = (S > 0).int()` deletes exactly one residue per step. A
~533-deletion target (TSC2/SHANK3 → 1274 aa) therefore runs ~533 sequential forward
passes per batch — the ~60 min/target cost. SYNGAP1 (69 deletions) is minutes.

Crucially, the machinery for multiple deletions already exists:
- `p_sample` deletes *k* at once: `torch.multinomial(tempered_probs, int(num_dels[i]))`.
- `sample_sequence(n_T=...)` already groups deletions into discretization intervals.

The only blocker in `shrink_sequence` is its preserved-index bookkeeping, hardcoded to one
deletion/step:

```python
preserved_indices = [
    l[:i] + l[i + 1 :]
    for l, i in zip(preserved_indices, [d[0] if d else 10000 for d in deleted_indices])
]
```

## 2. Primary lever — deletions-per-step (the headline)

This is the order-of-magnitude win and a small, contained change.

1. **Fix bookkeeping** in `shrink_sequence` to drop **all** deleted indices per step (not
   just `d[0]`), so k>1 deletions/step are tracked correctly. Keep the reconstruction
   assert in [../scisor_shrink.py](../scisor_shrink.py) green.
2. **Expose a knob** — `--dels-per-step` (equivalently a denoising-step count / `n_T`) on
   [../scisor_shrink.py](../scisor_shrink.py), reusing the interval-grouping already in
   `sample_sequence`.
3. **Sweep and find the knee** — dels-per-step ∈ {1 (baseline), 2, 4, 8, 16, 32}. Speed
   scoreboard for throughput; quality scoreboard (naturalness + motif retention on all;
   folding on the shortlist) for fidelity. Pick the **largest speedup whose quality stays
   within tolerance of the exact 1/step baseline**.

Why it can cost quality: deleting k positions in one step draws them from a single
distribution instead of re-evaluating the model after each deletion. The approximation
error grows with k — quantifying it on the scoreboard is the entire point.

## 3. Windowing wildcard

`window_size=2048` ([../configs/uniref90.yaml](../configs/uniref90.yaml)). The three
primary targets (<2048 aa) run full attention each step. Dystrophin (3685 aa) exceeds the
window and triggers **random-window sampling each step** (`predict_with_windows` /
`get_window_indices`). Report dystrophin as a separate calibration line; do not pool its
speed/quality numbers with the three targets.

## 4. Secondary levers (only if the knee leaves headroom)

These reduce per-forward cost, not the count:

- **FlashAttention on A100/H100.** The path exists but is forced off for the T4 in
  [../SCISOR/faesm.py](../SCISOR/faesm.py) (`use_fa=False`); `esm.py` already gates on
  `flash_attn_installed`. Install `flash_attn` and enable `use_fa=True` where the GPU
  supports it.
- **`torch.compile`** on the denoiser.
- **Batch-size scaling** on A100/H100 memory vs the T4's 16 GB.
- **Dependency refresh** (carefully): the `transformers==4.46.3` pin and the fresh ESM2
  tokenizer exist because newer `transformers` removed symbols `faesm` imports. Any bump
  must keep the ProteinGym guard green.

## 5. Milestones / outcomes

- **M-P1:** `shrink_sequence` multi-deletion bookkeeping fixed + `--dels-per-step` knob;
  reconstruction assert holds.
- **M-P2:** dels-per-step sweep with the chosen knee documented on both scoreboards;
  target: an order-of-magnitude reduction in forward passes on the ~530-deletion targets
  at equal quality, with the dystrophin caveat recorded.
- **M-P3 (optional):** secondary-lever micro-benchmarks (FA, compile, batch scaling).

## 6. Risks

- Deletions-per-step is an approximation; the sweep is the mitigation, not a guess.
- Dystrophin windowing behaves differently — keep it separate.
- Dependency bumps can silently change model behavior; the ProteinGym guard is the gate.
