# Benchmarking: frozen baseline + two scoreboards (WS0 operational)

This is the operational companion to [evaluation.md](evaluation.md). `evaluation.md`
defines *what* to score (the quality axes); this doc defines *how we benchmark
rigorously*: a frozen baseline, two independent scoreboards, a fold-free model-quality
guard, and the reproducibility rules every change must follow.

Rationale: SCISOR has two distinct kinds of improvement — **quality** (better real-world
shrinks of SYNGAP1/SHANK3/TSC2) and **speed** ([performance.md](performance.md)). These
must be measured on separate scoreboards so a speedup is never mistaken for a quality
change and vice versa. Nothing in the improvement program starts until this harness and a
committed baseline exist.

## 1. Compute placement

Hybrid, decided for this program:

- **SCISOR sampling** runs on a **local A100** (the denoiser is a 35M ESM2; sampling is
  light). A single **T4** run is kept only as the slow-baseline timing reference that
  matches the historical ~60 min/target number.
- **Folding-based quality eval** (the heavy tier) runs on the existing **`~/phi-api` k8s
  H100 runners** — `esmfold2` and `boltz2` already exist there
  (`infrastructure/cloudbuild/runners/{esmfold2,boltz2}.yaml`, weights jobs under
  `infrastructure/k8s/`). We call them as a service; we do **not** stand up folding from
  scratch. (`phi-api` also already has an `esmc` runner — a head start for
  [esmc-migration.md](esmc-migration.md).)

## 2. Speed scoreboard — `scisor_bench.py`

A new sibling of [scisor_shrink.py](../scisor_shrink.py) that wraps a shrink run and
records, per (target, config, GPU):

| Metric | Why |
|---|---|
| Wall-clock per target | Headline cost |
| Model forward-pass count | The serial bottleneck; the thing perf work reduces |
| Sequences/sec | Throughput at a given batch size |
| Peak GPU memory | Batch-size / GPU-fit ceiling |
| GPU-hours | Cost accounting across T4 / A100 / H100 |

Runs with fixed seeds and prints a tidy row per config so a sweep is one table.

## 3. Quality scoreboard — `scisor_score.py` + `SCISOR/scoring/`

Implements the layout in [evaluation.md](evaluation.md) §3, in two tiers:

- **Cheap tier (local A100, minutes, no folding):**
  - SCISOR naturalness / NLL of the shrunk sequence and per-position deletion probability.
  - Motif retention against the frozen-set schema (shared with
    [constrained-diffusion.md](constrained-diffusion.md)).
  - AAV-fit (length → bp; ≤ 3825 bp / ≤ 1274 aa).
  - Per-residue **deletion-frequency maps** (M1.2 tolerance overlay).
- **Heavy tier (phi-api H100, shortlists only):**
  - A thin folding client posts sequences to the `esmfold2` / `boltz2` runners and parses
    mean pLDDT / pTM and templated TM-score / RMSD per retained domain.
  - Apply the dual-depressor rule from [evaluation.md](evaluation.md) §1 (do not penalize
    low linker pLDDT per se).

The cheap tier gates everything; the heavy tier is reserved for ranked shortlists to keep
H100 spend bounded.

## 4. Model-quality regression guard — `scisor_proteingym.py`

A **fold-free correctness gate** so any speed change or refactor can be validated in
minutes, independent of the folding loop:

- Score deletion variants with SCISOR's predicted deletion log-prob.
- Compute **Spearman vs the ProteinGym indel/deletion DMS set** (the paper's SOTA claim).
- The baseline Spearman is committed; any change that moves it beyond tolerance is a
  regression, full stop.

Sub-task: acquire and pin the ProteinGym indel/deletion data. Until it lands,
naturalness / NLL is the interim fold-free guard.

## 5. Frozen baseline protocol

The reference everything is measured against:

1. Regenerate the **32-sample, T=1, one-deletion-per-step** baselines for
   **SYNGAP1, SHANK3, TSC2, and dystrophin** with **pinned seeds**.
2. Commit a **metrics manifest** (CSV/JSON) + the exact commands + seeds + GPU-hours under
   a tracked `benchmarks/` directory. Note `results/` and `*.fasta` are gitignored, so we
   commit **metrics, not raw FASTAs**.
3. Record the baseline **ProteinGym Spearman** alongside.

Dystrophin is a separate calibration line, never pooled with the three primary targets:
at 3685 aa it exceeds `window_size=2048` ([configs/uniref90.yaml](../configs/uniref90.yaml))
and triggers random-window sampling each step (see [performance.md](performance.md)).

## 6. Acceptance rules (apply to every change)

- Report deltas on **both** scoreboards vs the frozen baseline — never one in isolation.
- ProteinGym parity must **not regress** beyond tolerance.
- Pin dependencies, fix seeds, record the exact command and GPU-hours for reproducibility.
- A speed change must show its quality delta; a quality change must show its speed cost.

## 7. Milestones / outcomes

- **Harness v1:** one command → speed scoreboard; one command → ranked quality table.
- **Frozen baseline:** committed metrics + ProteinGym Spearman for the 3 targets +
  dystrophin.
- **Guard:** `scisor_proteingym.py` reproduces the paper-level deletion-effect parity.
