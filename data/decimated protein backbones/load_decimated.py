#!/usr/bin/env python3
"""Loader for the decimated protein backbone dataset.

    from load_decimated import load, load_pair, keep_indices, CHAINS, WS

    A = load("1o7j.a", w=4)              # (82, 3) float array, aligned frame
    A, B = load_pair("4eca.c", w=4)      # both chains of the pair, same indices
    rows = keep_indices(325, 4)          # which source rows survived

Everything is index based, so `keep_indices` alone reproduces the decimation
from any 325-row chain without reading the CSV.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CHAINS = ["1o7j.a", "1hfj.c", "4eca.b", "4eca.c", "4eca.d",
          "1qd1.b", "1toh.a", "1d9q.d"]
REF = "1o7j.a"
WS = (1, 2, 4, 8, 16)
FRAMES = ("aligned", "deposited")

_cache = {}


def keep_indices(n: int, w: int) -> np.ndarray:
    """The exact decimator. Uniform stride from 0, last vertex always kept."""
    if w == 1:
        return np.arange(n)
    idx = list(range(0, n, w))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return np.array(idx)


def _table(name):
    if name not in _cache:
        _cache[name] = pd.read_csv(os.path.join(HERE, name))
    return _cache[name]


def curves():
    return _table("decimated_curves.csv")


def chains():
    return _table("decimated_chains.csv")


def pairs():
    return _table("decimated_pairs.csv")


def load(chain: str, w: int = 1, frame: str = "aligned") -> np.ndarray:
    """One chain as an (n, 3) array, vertices in order."""
    t = curves()
    g = t[(t.chain == chain) & (t.w == w) & (t.frame == frame)].sort_values("t")
    if g.empty:
        raise KeyError(f"no rows for chain={chain!r} w={w} frame={frame!r}")
    return g[["x", "y", "z"]].to_numpy(float)


def source_rows(chain: str, w: int = 1, frame: str = "aligned") -> np.ndarray:
    """Which rows of the original 325-vertex chain these vertices came from."""
    t = curves()
    g = t[(t.chain == chain) & (t.w == w) & (t.frame == frame)].sort_values("t")
    return g.source_row.to_numpy(int)


def load_pair(partner: str, w: int = 1, frame: str = "aligned"):
    """(A, B) for the pair 1o7j.a / partner, decimated at the same indices."""
    return load(REF, w, frame), load(partner, w, frame)


def thresholds(partner: str, w: int = 1, frame: str = "aligned", alpha: float = 1.0):
    """delta1, delta2, delta3 for one instance, as the pipeline defines them."""
    t = pairs()
    r = t[(t.pair == f"{REF}/{partner}") & (t.w == w) & (t.frame == frame)]
    if r.empty:
        raise KeyError(f"no row for partner={partner!r} w={w} frame={frame!r}")
    r = r.iloc[0]
    d1 = float(r[f"delta1_alpha{alpha:g}"])
    return d1, d1, float(r.delta3)


if __name__ == "__main__":
    print(f"{'pair':18} {'w':>3} {'n':>4} {'d_dF':>7} {'delta3':>7} {'floor':>7} {'ok':>5}")
    for p in CHAINS[1:]:
        for w in WS:
            A, B = load_pair(p, w)
            d1, d2, d3 = thresholds(p, w)
            r = pairs()
            r = r[(r.pair == f"{REF}/{p}") & (r.w == w) & (r.frame == "aligned")].iloc[0]
            print(f"{REF}/{p:10} {w:>3} {len(A):>4} {r.dF_discrete:>7.2f} "
                  f"{r.delta3:>7.0f} {r.feasibility_floor:>7.2f} {str(bool(r.feasible)):>5}")
