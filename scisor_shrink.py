#!/usr/bin/env python
"""SCISOR shrink runner (patched for inference on real, long proteins).

Faithful adaptation of upstream ``shrink_sequences.py`` with the changes needed
to actually run SCISOR on long therapeutic targets and to draw multiple
stochastic samples:

1. ``read_fasta_to_df`` no longer applies the upstream ``.head(100).query("Length <= 1000")``
   filter, which silently dropped every sequence longer than 1000 aa.
2. The checkpoint is loaded from a local path (default: ``weights/SCISOR_U90_S.ckpt``
   next to this script) instead of re-downloading from HuggingFace on every run.
3. Shrinking can be driven by a percentage (``--shrink-pct``) or an explicit
   ``--target-length`` (deletions = max(0, L - target_length)). The latter is the
   right mode for hard size budgets (e.g. AAV packaging).
4. ``--num-samples`` draws N independent stochastic shrinks per input sequence
   (SCISOR samples deletions via multinomial at ``temperature>0``); use
   ``--temperature 0`` for a single deterministic (top-k) shrink.

Output FASTA keeps upstream-style deletion-map headers (deletion positions in
original coordinates) and adds a ``|sample k`` field.
"""
import argparse
import os

import pandas as pd
import torch
import torch.nn.functional as F
from tqdm import tqdm
from transformers import AutoTokenizer

from SCISOR.shortening_scud import ShorteningSCUD

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ESM_TOKENIZER = "facebook/esm2_t6_8M_UR50D"


def read_fasta_to_df(fasta_file):
    with open(fasta_file, "r") as file:
        content = file.read()
    sequences = []
    for entry in content.strip().split(">"):
        if not entry:
            continue
        lines = entry.strip().split("\n")
        header = lines[0]
        sequence = "".join(lines[1:]).replace(" ", "").replace("\r", "")
        sequences.append(
            {"Header": header, "Sequence": sequence, "Length": len(sequence)}
        )
    # NOTE: upstream applied .head(100).query("Length <= 1000") here; removed so
    # long real targets (SYNGAP1/SHANK3/TSC2/dystrophin) are not silently dropped.
    return pd.DataFrame(sequences).drop_duplicates().reset_index(drop=True)


def untokenize(seq, tokenizer):
    return (
        tokenizer.decode(seq)
        .replace(" ", "")
        .replace("<cls>", "")
        .replace("<eos>", "")
        .replace("<pad>", "")
    )


def save_sequences_to_fasta(fasta_file, seqs, headers):
    if os.path.dirname(fasta_file):
        os.makedirs(os.path.dirname(fasta_file), exist_ok=True)
    with open(fasta_file, "w") as f:
        for header, seq in zip(headers, seqs):
            f.write(f">{header}\n{seq}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="Input FASTA of long sequences")
    ap.add_argument("--output", default="shrunk_sequences.fasta")
    ap.add_argument(
        "--ckpt",
        default=os.path.join(SCRIPT_DIR, "weights", "SCISOR_U90_S.ckpt"),
        help="Local SCISOR checkpoint path",
    )
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--shrink-pct", type=float, help="Delete this %% of each seq")
    group.add_argument(
        "--target-length", type=int, help="Delete down to this absolute length"
    )
    ap.add_argument("--temperature", type=float, default=1.0,
                    help="0 = deterministic top-k; >0 = stochastic multinomial")
    ap.add_argument("--num-samples", type=int, default=1,
                    help="Independent stochastic shrinks per input sequence")
    ap.add_argument("--batch-size", type=int, default=0, help="0 = model default")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    if args.seed is not None:
        torch.manual_seed(args.seed)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Loading checkpoint: {args.ckpt}")
    model = ShorteningSCUD.load_from_checkpoint(args.ckpt, map_location=device)
    model.to(device)
    model.eval()
    model.p0 = torch.load(os.path.join(SCRIPT_DIR, "p0.pt"), map_location=device)
    # The checkpoint pickles an EsmTokenizer built under an older transformers that
    # is missing attributes (e.g. _unk_token) on newer transformers. Replace it with
    # a fresh ESM2 tokenizer (identical standard vocab used during training).
    model.tokenizer = AutoTokenizer.from_pretrained(ESM_TOKENIZER)
    rate = 1 / 1.1
    model.alpha = lambda t: (1 - t) ** rate
    model.beta = lambda t: rate / (1 - t)

    base = read_fasta_to_df(args.input)
    # Replicate each input sequence num_samples times.
    df = base.loc[base.index.repeat(args.num_samples)].reset_index(drop=True)
    df["Sample"] = list(range(args.num_samples)) * len(base)

    seq_lengths = torch.tensor(df.Length.values, device=model.device)
    if args.shrink_pct is not None:
        num_deletions = torch.ceil(seq_lengths * args.shrink_pct / 100).int()
        mode = f"{args.shrink_pct}pct"
    else:
        num_deletions = torch.clamp(seq_lengths - args.target_length, min=0).int()
        mode = f"L{args.target_length}"
    print(
        f"Shrinking {len(base)} target(s) x {args.num_samples} sample(s) "
        f"= {len(df)} seqs (mode {mode}, T={args.temperature}) on {device}"
    )
    for h, L in zip(base.Header, base.Length):
        if args.target_length is not None:
            nd = max(L - args.target_length, 0)
            print(f"  {h[:48]:48s} L={L:5d} delete={nd:5d} -> {L - nd}")
        else:
            nd = -(-L * args.shrink_pct // 100)  # ceil
            print(f"  {h[:48]:48s} L={L:5d} delete={int(nd):5d} -> {L - int(nd)}")

    input_ids = [model.tokenizer(s).input_ids for s in df.Sequence]
    max_len = max(len(x) for x in input_ids)
    x = torch.vstack(
        [
            F.pad(
                torch.tensor(t, device=model.device),
                (0, max_len - len(t)),
                value=model.tokenizer.pad_token_id,
            )
            for t in input_ids
        ]
    )

    batch_size = args.batch_size or model.hparams.batch_size
    sampled_sequences, deleted_indices = [], []
    for i in tqdm(range(0, len(x), batch_size)):
        sequences, preserved_indices = model.shrink_sequence(
            x[i : i + batch_size],
            num_deletions[i : i + batch_size],
            temperature=args.temperature,
        )
        del_idx = [
            list(set(range(len(s))) - set([j - 1 for j in p]))
            for s, p in zip(df.Sequence[i : i + batch_size], preserved_indices)
        ]
        sampled_sequences.extend(untokenize(s, model.tokenizer) for s in sequences)
        deleted_indices.extend(del_idx)

    assert all(
        "".join(c for i, c in enumerate(s) if i not in d) == n
        for s, d, n in zip(df.Sequence, deleted_indices, sampled_sequences)
    ), "Reconstruction check failed"

    del_str = [
        ",".join(f"{c}{i}" for i, c in enumerate(s) if i in d)
        for s, d in zip(df.Sequence, deleted_indices)
    ]
    new_headers = (
        df.Header
        + "|sample " + df.Sample.astype(str)
        + "|mode " + mode
        + "|deletions " + pd.Series(del_str, index=df.index)
    )
    save_sequences_to_fasta(args.output, sampled_sequences, new_headers)
    print(f"Saved {len(sampled_sequences)} shrunk sequences to {args.output}")


if __name__ == "__main__":
    main()
