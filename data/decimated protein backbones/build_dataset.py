#!/usr/bin/env python3
"""Build the decimated-curve dataset: tidy CSV of coordinates plus a threshold table.

The decimator is the one used throughout: uniform stride from index 0, the last
vertex always appended if the stride missed it, and the SAME index set applied to
both chains of a pair (lockstep), because delta3 is index-wise and the index
correspondence is the correspondence.

Two coordinate frames are emitted for every chain:
  deposited  the coordinates exactly as in the PDB entry, no superposition
  aligned    the same chains after a rigid Kabsch fit onto 1o7j.a

The decimation is purely index-based, so the retained source rows are identical
in both frames; only the coordinates differ.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cps_harness import feasibility_floor
from frechet import continuous_frechet, discrete_frechet

FRAMES = {"deposited": "/home/claude/proteins/curves",
          "aligned": "/home/claude/datasets/proteins_aligned"}
OUT = "/home/claude/datasets/decimated"
REF = "1o7j.a"
PARTNERS = ["1hfj.c", "4eca.b", "4eca.c", "4eca.d", "1qd1.b", "1toh.a", "1d9q.d"]
CHAINS = [REF] + PARTNERS
WS = (1, 2, 4, 8, 16)
ALPHAS = (1.0, 2.0, 4.0)
GAP_FACTOR = 1.3          # an edge above this multiple of the median is a gap


def keep_indices(n: int, w: int) -> np.ndarray:
    """The retained source rows. Endpoints are always in."""
    if w == 1:
        return np.arange(n)
    idx = list(range(0, n, w))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return np.array(idx)


def edge_stats(P):
    e = np.linalg.norm(np.diff(P, axis=0), axis=1)
    return dict(edge_mean=float(e.mean()), edge_min=float(e.min()),
                edge_max=float(e.max()), edge_last=float(e[-1]),
                arc_length=float(e.sum()))


def gaps(P, factor=GAP_FACTOR):
    e = np.linalg.norm(np.diff(P, axis=0), axis=1)
    hit = np.where(e > factor * np.median(e))[0]
    return [(int(i), round(float(e[i]), 2)) for i in hit]


def main():
    os.makedirs(OUT, exist_ok=True)

    raw = {f: {c: np.loadtxt(f"{d}/{c}.txt") for c in CHAINS} for f, d in FRAMES.items()}
    n = len(raw["aligned"][REF])
    assert all(len(P) == n for f in raw for P in raw[f].values()), "chains differ in length"

    # ---- coordinates -------------------------------------------------------
    rows = []
    for w in WS:
        idx = keep_indices(n, w)
        for frame in FRAMES:
            for c in CHAINS:
                P = raw[frame][c][idx]
                for t, (src, xyz) in enumerate(zip(idx, P)):
                    rows.append((c, frame, w, t, int(src), *np.round(xyz, 3)))
    curves = pd.DataFrame(rows, columns=["chain", "frame", "w", "t", "source_row",
                                         "x", "y", "z"])
    curves.to_csv(f"{OUT}/decimated_curves.csv", index=False)

    # ---- per-chain summary -------------------------------------------------
    crows = []
    for w in WS:
        idx = keep_indices(n, w)
        stride_last = int(idx[-1] - idx[-2])
        for frame in FRAMES:
            for c in CHAINS:
                P = raw[frame][c][idx]
                g = gaps(P)
                crows.append(dict(chain=c, frame=frame, w=w, n_vertices=len(idx),
                                  n_before=n, keep_frac=round(len(idx) / n, 4),
                                  last_stride=stride_last,
                                  appended_stub=bool(stride_last != w and w > 1),
                                  **{k: round(v, 3) for k, v in edge_stats(P).items()},
                                  n_gaps=len(g), gaps=json.dumps(g)))
    chains = pd.DataFrame(crows)
    base = chains[chains.w == 1].set_index(["chain", "frame"]).edge_mean
    chains["edge_growth_factor"] = [
        round(r.edge_mean / base[(r.chain, r.frame)], 3) for r in chains.itertuples()]
    chains["coarse_trace"] = chains.edge_mean > 10.0
    chains.to_csv(f"{OUT}/decimated_chains.csv", index=False)

    # ---- per-pair thresholds ----------------------------------------------
    prows = []
    for w in WS:
        idx = keep_indices(n, w)
        for frame in FRAMES:
            A = raw[frame][REF][idx]
            for nb in PARTNERS:
                B = raw[frame][nb][idx]
                ddf = discrete_frechet(A, B)
                dfc = continuous_frechet(A, B)
                fl = feasibility_floor(A, B)
                d3 = float(np.ceil(ddf))
                e = 0.5 * (edge_stats(A)["edge_mean"] + edge_stats(B)["edge_mean"])
                r = dict(pair=f"{REF}/{nb}", chain_A=REF, chain_B=nb, frame=frame, w=w,
                         superimposed=(frame == "aligned"),
                         n_A=len(idx), n_B=len(idx), n_before=n,
                         edge_mean=round(e, 3),
                         dF_discrete=round(ddf, 3), dF_continuous=round(dfc, 3),
                         feasibility_floor=round(fl, 3),
                         delta3=d3, feasible=bool(d3 >= fl),
                         endpoint_floor=round(float(max(
                             np.linalg.norm(A[0] - B[0]),
                             np.linalg.norm(A[-1] - B[-1]))), 3))
                for a in ALPHAS:
                    r[f"delta1_alpha{a:g}"] = round(a * e, 3)
                prows.append(r)
    pairs = pd.DataFrame(prows)
    pairs.to_csv(f"{OUT}/decimated_pairs.csv", index=False)

    print(f"decimated_curves.csv  {len(curves):6d} rows")
    print(f"decimated_chains.csv  {len(chains):6d} rows")
    print(f"decimated_pairs.csv   {len(pairs):6d} rows")
    return curves, chains, pairs


if __name__ == "__main__":
    main()
