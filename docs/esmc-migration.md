# ESM2 → ESMC migration (WS4)

Goal: replace SCISOR's ESM2 backbone with **ESMC** (ESM Cambrian), which outperforms
ESM2 in our hands. This requires retraining the denoiser and is the most
infrastructure-heavy workstream; treat it as parallelizable R&D with a clear go/no-go.

> Head start: `~/phi-api` already runs an **`esmc` runner** on k8s H100s
> (`infrastructure/cloudbuild/runners/esmc.yaml`, `populate-esmc-weights-job.yaml`), so
> ESMC weights/inference are already provisioned — reuse it for the embedding-layer sweep
> (M4.1) instead of bootstrapping ESMC from scratch. The ProteinGym guard in
> [benchmarking.md](benchmarking.md) is the M4.2 parity gate.

## 1. Where ESM2 lives today

SCISOR's denoiser is an ESM2 model with two additions: a per-layer FiLM conditioning on
the diffusion variable, and a linear head projecting to a per-position deletion logit.

- [SCISOR/SCISOR/faesm.py](../SCISOR/faesm.py) — `FAESM_Base` wraps
  `FAEsmForMaskedLM.from_pretrained("facebook/esm2_t12_35M_UR50D", ...)`, takes the
  `last_hidden_state`, and applies `self.proj = nn.Linear(embed_dim, 1)` to get the
  deletion prediction. Conditioning (`t` or `S`) is passed in as `conditioning`.
- [SCISOR/SCISOR/esm.py](../SCISOR/esm.py) — vendored FAESM with `FAEsmEncoder` adding a
  FiLM layer after each transformer block (`film_layers`, `conditioning_embedder`),
  zero-initialized so training starts from the pretrained ESM2.
- The diffusion math (`shortening_scud.py`, `shortening_diffusion.py`) is
  backbone-agnostic: it consumes a per-position score from `x0_model`, so swapping the
  backbone is contained to the two files above plus retraining.

## 2. What changes with ESMC

ESMC is not a drop-in for the HuggingFace ESM2 classes:
- **API/SDK.** ESMC ships via the ESM SDK (local weights) and/or the Forge API, not the
  `transformers` `EsmModel` path FAESM subclasses. The FiLM-after-each-layer hook must
  be re-implemented against ESMC's transformer stack (or via forward hooks if internals
  aren't subclassable).
- **Multi-layer embeddings.** ESMC exposes richer per-layer representations that encode
  different aspects of protein information. Which layer(s) to feed the denoiser is an
  open choice (last layer vs a mid layer vs a learned combination). This needs an
  **embedding-layer sweep** — a small study selecting the layer(s) that best predict
  deletion effects.
- **Tokenizer/vocab.** Confirm ESMC's vocabulary and special tokens; the deletion logic
  in `p_sample` masks cls/eos/pad by id, so those ids must be remapped.
- **No FAESM flash wrapper.** The T4/SDPA workaround (`use_fa=False`) is ESM2/FAESM
  specific; ESMC has its own attention implementation. On Ampere (A100) flash is
  available; on T4 confirm a non-flash path exists for inference.
- **Licensing.** Confirm ESMC license terms for our use before committing.

## 3. Work plan

1. **Adapter (M4.1).** New backbone module mirroring `FAESM_Base`'s interface
   (`forward(x, t, input_mask, S) -> per-position logits`) but built on ESMC. Re-wire
   FiLM conditioning (zero-init, as today) and the linear deletion head. Decide the
   conditioned embedding layer(s) (start: last layer; then sweep).
2. **Embedding-layer sweep (M4.1).** On a held-out set (ProteinGym deletion effects),
   compare which ESMC layer(s) give the best deletion-effect prediction; fix the choice.
3. **Retrain SCISOR-C (M4.2).** Train the denoiser on UniRef (UniRef90 to match U90).
   Reuse the existing training entry (`train.py` + configs) with a new model config.
   *Compute:* A100 — the paper trained S/M ~1 week on 2× A100, L ~4 days on 4× H100;
   budget accordingly. Watch `val_l01` down, as upstream notes.
4. **Re-benchmark (M4.2).** ProteinGym deletion-effect parity vs ESM2 SCISOR first
   (correctness gate), then quality.
5. **Targets (M4.3).** Re-run baselines + WS2 constrained runs with SCISOR-C on
   TSC2 / SHANK3 / SYNGAP1; compare on the WS0 scorecard.

## 4. Milestones / outcomes
- **M4.1:** ESMC adapter with FiLM + head; forward pass and a short training run
  converge; embedding-layer choice fixed.
- **M4.2:** SCISOR-C trained; **ProteinGym parity** with ESM2 SCISOR (go/no-go gate).
- **M4.3:** SCISOR-C scored on our targets; documented improvement or a clear no-go.

## 5. Risks / decision points
- **Cost before signal.** Don't commit A100-weeks until the adapter + sweep show ESMC
  embeddings predict deletions at least as well as ESM2 (decision gate at M4.1).
- **Conditioning placement.** FiLM-per-layer worked for ESM2; ESMC's depth/width may
  need a different conditioning scheme.
- **Reproducibility.** Pin ESMC SDK version and weights; the inference patches
  (transformers pin, tokenizer) are ESM2-specific and won't carry over unchanged.
