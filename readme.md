# SCISOR: A Diffusion Model to Shrink Proteins While Maintaining Their Function

[Ethan Baron](https://baronet2.github.io/)\*, [Alan N. Amin](https://alannawzadamin.github.io)\*, [Ruben Weitzman](https://rubenweitzman.github.io/), [Debora Susan Marks](https://www.deboramarkslab.com/deboramarks), [Andrew Gordon Wilson](https://cims.nyu.edu/~andrewgw/). * equal contribution

Paper: https://openreview.net/pdf?id=YqQoNJWY22

![SCISOR learns to reverse a noising process of random insertions](images/concept.png)

## Abstract
Many proteins useful in modern medicine or bioengineering are challenging to make in the lab, fuse with other proteins in cells, or deliver to tissues in the body because their sequences are too long. Shortening these sequences typically involves costly, time-consuming experimental campaigns. Ideally, we could instead use modern models of massive databases of sequences from nature to learn how to propose shrunken proteins that resemble sequences found in nature. Unfortunately, these models struggle to efficiently search the combinatorial space of all deletions, and are not trained with inductive biases to learn how to delete.

To address this gap, we propose SCISOR, a novel discrete diffusion model that deletes letters from sequences to generate protein samples that resemble those found in nature. To do so, SCISOR trains a de-noiser to reverse a forward noising process that adds random insertions to natural sequences.

As a generative model, SCISOR fits evolutionary sequence data competitively with previous large models. In evaluation, SCISOR achieves state-of-the-art predictions of the functional effects of deletions on ProteinGym. Finally, we use the SCISOR de-noiser to shrink long protein sequences, and show that its suggested deletions result in significantly more realistic proteins and more often preserve functional motifs than previous models of evolutionary sequences.

## Citation

If you use SCISOR, please cite our paper:

```
@inproceedings{SCISOR,
  title        = {A Diffusion Model to Shrink Proteins While Maintaining their Function},
  author       = {Baron, Ethan and Amin, Alan\,N. and Weitzman, Ruben and Marks, Debora\,S. and Wilson, Andrew\,G.},
  booktitle    = {ICML Workshop on Exploration in AI Today},
  year         = {2025},
  month        = {Jul},
  url          = {https://openreview.net/pdf?id=YqQoNJWY22},
}
```

## Usage

We provide the trained SCISOR models on HuggingFace at https://huggingface.co/SCISOR/SCISOR/tree/main.

### Installation

Install dependencies by running `pip install .` with a recent version of Python.

### Generating Unconditional Samples with SCISOR

Running `python sample_sequences.py` will produce a file `sampled_sequences.fasta`.

### Shrinking Proteins with SCISOR

Running `python shrink_sequences.py` takes in a `fasta` file like `sampled_sequences.fasta`, and produces a file `shrunk_sequences.fasta`.

### Shrinking real (long) proteins for inference — `scisor_shrink.py`

`shrink_sequences.py` is geared at the paper's benchmark and has a couple of sharp
edges for applied use. `scisor_shrink.py` is a drop-in inference runner that fixes them:

- **No silent length filter.** Upstream's `read_fasta_to_df` ends with
  `.head(100).query("Length <= 1000")`, which silently drops every sequence longer
  than 1000 aa. `scisor_shrink.py` keeps all sequences.
- **Local checkpoints.** Loads from `weights/SCISOR_U90_S.ckpt` (configurable via
  `--ckpt`) instead of re-downloading from HuggingFace each run.
- **Absolute size budgets.** `--target-length N` deletes down to exactly N residues
  (e.g. AAV packaging limits), in addition to `--shrink-pct`.
- **Stochastic sampling.** `--num-samples K` draws K independent shrinks per input
  (SCISOR samples deletions via multinomial); `--temperature 0` gives a single
  deterministic top-k shrink.

```sh
python scisor_shrink.py --input targets.fasta --output shrunk.fasta \
    --target-length 1274 --num-samples 32 --batch-size 16
```

**Running on non-Ampere GPUs / without FlashAttention.** FlashAttention-2 does not
support pre-Ampere GPUs (e.g. Tesla T4, sm_75), and `faesm`'s flash rotary path
hard-requires `flash_attn`. The ESM backbone is therefore constructed with
`use_fa=False` (see `SCISOR/faesm.py`), which keeps HuggingFace's pure-PyTorch
RotaryEmbedding and runs attention via Torch SDPA — the same path the model already
falls back to at forward time when `flash_attn` is unavailable. For inference you do
**not** need the training-only dependencies (`evodiff`, `faiss-cpu`, `hydra-core`);
`pip install -e . --no-deps` plus `torch torchvision pytorch_lightning transformers
faesm pandas tqdm wandb huggingface_hub einops omegaconf` is sufficient. Note that
recent `transformers` (5.x) removed symbols `faesm` imports from the ESM module;
pin `transformers==4.46.3`.

## Improvement program (planned work)

This fork is being evolved from a naturalness-only deletion model into a constrained,
structure-aware protein-minimization platform, with rigorous benchmarking. If you are
picking this up on a fresh machine, start with the docs below — they are the source of
truth for the plan.

- **[docs/ROADMAP.md](docs/ROADMAP.md)** — top-level plan, workstreams (WS0–WS5 + WS-P),
  phased timeline, compute plan.
- **[docs/benchmarking.md](docs/benchmarking.md)** — how we measure rigorously: frozen
  baseline, two scoreboards (speed vs quality), the fold-free ProteinGym guard.
- **[docs/performance.md](docs/performance.md)** — the separate **speed** path
  (deletions-per-step; SCISOR currently deletes one residue per forward pass, serially).
- **[docs/evaluation.md](docs/evaluation.md)** — quality scoring axes and rational-design
  comparison.
- **[docs/constrained-diffusion.md](docs/constrained-diffusion.md)** — keep-mask domain
  locking (the first **quality** improvement).
- **[docs/structural-integrity-linkers.md](docs/structural-integrity-linkers.md)** — block
  deletion + linker repair.
- **[docs/esmc-migration.md](docs/esmc-migration.md)** — ESM2 → ESMC backbone swap.

**Start here (first improvements).** Phase 0 is a gate: build the harness and freeze a
baseline. Then two independent experiments run in parallel:

- **Path A (speed)** — deletions-per-step in [docs/performance.md](docs/performance.md).
- **Path B (quality)** — keep-masks in [docs/constrained-diffusion.md](docs/constrained-diffusion.md).

**Compute.** Hybrid: SCISOR sampling on a local A100; folding-based quality eval via the
existing `~/phi-api` k8s **H100** runners (`esmfold2`, `boltz2`); one T4 run kept only as
the slow-baseline timing reference.

**Planned scripts / flags (not yet implemented — see the docs):** `scisor_bench.py`
(speed scoreboard), `scisor_score.py` + `SCISOR/scoring/` (quality scoreboard),
`scisor_proteingym.py` (model-quality guard), and new
[scisor_shrink.py](scisor_shrink.py) flags `--dels-per-step` (speed) and `--keep-mask`
(locks).

### Training SCISOR with Uniref50

To train SCISOR on Uniref50, start by downloading the data and pre-processing it into batches:

```sh
wget -O uniref50.tar.gz "https://zenodo.org/records/6564798/files/uniref50.tar.gz?download=1"
tar -xvzf uniref50.tar.gz -C data
python data/preprocess_uniref50.py
```

Currently, `preprocess_uniref50.py` is set up to process 50 batches of 256 sequences, for demonstration purposes. If you wish to train on more data, increase `num_training_batches`.

Then run:
```sh
python train.py  --config-name=uniref50
```

If using `wandb`, the key is to watch `val_l01` go down! We suggest training with an A100 GPU.

The config parameters can be updated to specify the dataset path, model architecture, training hyperparameters, sampling settings, etc.

### Training SCISOR with Uniref90

To download the Uniref90 dataset, use:
```sh
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz
```

We provide the code we used to pre-process Uniref90. Due to this dataset's size, we follow the following steps:
1. Shuffle the dataset and save the shuffled dataset in separate shards (see `data/shuffle_uniref90.py`)
2. Preprocess each shard individually using `data/preprocess_uniref90_shard.py`

After preprocessing, the model can be trained with:
```sh
python train.py  --config-name=uniref90
```