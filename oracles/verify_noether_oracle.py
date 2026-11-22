#!/usr/bin/env python3
"""
PHY v4.0.0 — Independent Exact-Rational Level-2 Verification Oracle (Gate 5).
Independently re-verifies all 15 discovered physics/algebra modes (entangled, yangbox,
joint, q4, yb8, dual4, ame8, fracton, biquandle8, tetra4, ame4-nogo, simplex4, chiral6,
qscar, fracton-comp) using Python arbitrary-precision ints, fractions.Fraction
Gaussian elimination, multi-length Floquet / L=32 circuit simulations, exact Z[i] quantum
circuit traces, HSF L=6,7,8 scaling collapse, .fir bit-exact parity checks, live
Lean 4 + leanchecker replay with 0 sorryAx, and QUA/ELE cross-ecosystem projection
(when sibling builds are present).
"""

from __future__ import annotations

import concurrent.futures
import json
import math
import multiprocessing
import os
import random
import shutil
import subprocess
import sys
from collections import deque
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def _find_elan_tool(name: str) -> Path:
    env_override = os.environ.get(name.upper())
    if env_override:
        return Path(env_override)
    which_p = shutil.which(name)
    if which_p:
        return Path(which_p)
    home_p = Path.home() / ".elan" / "bin" / name
    if home_p.exists():
        return home_p
    return Path("/Users/kuangtalin/.elan/bin") / name


ROOT = Path(__file__).resolve().parent.parent
BUILD = ROOT / "build"
FLOW_DIR = Path(os.environ.get("FLOW_DIR", str(ROOT.parent / "FLOW"))).resolve()
PHY_BIN = BUILD / "phy"
LEAN = _find_elan_tool("lean")
LEANCHECKER = _find_elan_tool("leanchecker")


def pack2(l: int, r: int) -> int:
    return ((l & 3) << 2) | (r & 3)


def hi2(u: int) -> int:
    return (u >> 2) & 3


def lo2(u: int) -> int:
    return u & 3


def step3(p: List[int], a: int, b: int, c: int) -> int:
    pa, pb, pc = p[a], p[b], p[c]
    odd_l = p[pack2(lo2(pa), hi2(pb))]
    odd_r = p[pack2(lo2(pb), hi2(pc))]
    return pack2(lo2(odd_l), hi2(odd_r))


def anf_degree(p: List[int]) -> int:
    max_d = 0
    for bit in range(4):
        f = [(p[x] >> bit) & 1 for x in range(16)]
        for v in range(4):
            for mask in range(16):
                if mask & (1 << v):
                    f[mask] ^= f[mask ^ (1 << v)]
        for mask in range(16):
            if f[mask]:
                max_d = max(max_d, bin(mask).count("1"))
    return max_d


def is_factorizable(p: List[int]) -> bool:
    for a in range(4):
        h0 = hi2(p[pack2(a, 0)])
        for b in range(1, 4):
            if hi2(p[pack2(a, b)]) != h0:
                return False
    for b in range(4):
        l0 = lo2(p[pack2(0, b)])
        for a in range(1, 4):
            if lo2(p[pack2(a, b)]) != l0:
                return False
    return True


def is_lyubashenko(p: List[int]) -> bool:
    for a in range(4):
        for b in range(4):
            if hi2(p[pack2(a, b)]) != hi2(p[pack2(0, b)]):
                return False
            if lo2(p[pack2(a, b)]) != lo2(p[pack2(a, 0)]):
                return False
    return True


def yb_violations(p: List[int]) -> int:
    viol = 0
    for a in range(4):
        for b in range(4):
            for c in range(4):
                ab = p[pack2(a, b)]
                a1, b1 = hi2(ab), lo2(ab)
                bc = p[pack2(b1, c)]
                b2, c2 = hi2(bc), lo2(bc)
                ab2 = p[pack2(a1, b2)]
                lhs = (hi2(ab2), lo2(ab2), c2)

                bc0 = p[pack2(b, c)]
                b0, c0 = hi2(bc0), lo2(bc0)
                ab0 = p[pack2(a, b0)]
                a0, b00 = hi2(ab0), lo2(ab0)
                bc1 = p[pack2(b00, c0)]
                rhs = (a0, hi2(bc1), lo2(bc1))

                if lhs != rhs:
                    viol += 1
    return viol


def fraction_rank(rows: List[List[Fraction]], n_cols: int) -> int:
    pivs: Dict[int, List[Fraction]] = {}
    for row in rows:
        v = list(row)
        for col in sorted(pivs.keys()):
            if v[col] != 0:
                fac = v[col]
                pr = pivs[col]
                for j in range(col, n_cols):
                    v[j] -= fac * pr[j]
        lead = next((j for j, x in enumerate(v) if x != 0), None)
        if lead is not None:
            inv = Fraction(1) / v[lead]
            pivs[lead] = [x * inv for x in v]
    return len(pivs)


def verify_nonseparable_fraction(q: List[int]) -> bool:
    sep_rows = []
    for h in (1, 2, 3):
        sep_rows.append([Fraction(1 if hi2(u) == h else 0) for u in range(16)])
    for l in (1, 2, 3):
        sep_rows.append([Fraction(1 if lo2(u) == l else 0) for u in range(16)])
    const_row = [Fraction(1) for _ in range(16)]
    base = sep_rows + [const_row]
    r0 = fraction_rank(base, 16)
    r1 = fraction_rank(base + [[Fraction(x) for x in q]], 16)
    return r1 > r0


def verify_continuity_exact(p: List[int], q: List[int]) -> Tuple[int, int]:
    f00 = [step3(p, 0, 0, x) for x in range(16)]
    f0 = [[step3(p, 0, x, y) for y in range(16)] for x in range(16)]
    j = [[(q[f00[b]] - q[0]) + (q[f0[b][c]] - q[b]) for c in range(16)] for b in range(16)]
    j_nz = sum(1 for b in range(16) for c in range(16) if j[b][c] != 0)
    for a in range(16):
        for b in range(16):
            for c in range(16):
                lhs = q[step3(p, a, b, c)] - q[b]
                rhs = j[b][c] - j[a][b]
                if lhs != rhs:
                    raise AssertionError(
                        f"Continuity mismatch at ({a},{b},{c}): lhs={lhs} != rhs={rhs}"
                    )
    return j_nz, sum(abs(x) for x in q)


def floquet_step_list(p: List[int], u: List[int]) -> List[int]:
    n = len(u)
    return [step3(p, u[(i - 1) % n], u[i], u[(i + 1) % n]) for i in range(n)]


def floquet_inv_step_list(pinv: List[int], u: List[int]) -> List[int]:
    """Exact inverse Floquet step using odd(P^-1) followed by even(P^-1) on 2n sites."""
    n = len(u)
    sites = []
    for cell in u:
        sites.append(hi2(cell))
        sites.append(lo2(cell))
    m = 2 * n
    # Inverse of odd layer on (2i+1, 2i+2)
    mid = sites[:]
    for i in range(n):
        idx0, idx1 = (2 * i + 1) % m, (2 * i + 2) % m
        val = pinv[pack2(sites[idx0], sites[idx1])]
        mid[idx0], mid[idx1] = hi2(val), lo2(val)
    # Inverse of even layer on (2i, 2i+1)
    out = [0] * n
    for i in range(n):
        out[i] = pinv[pack2(mid[2 * i], mid[2 * i + 1])]
    return out


def verify_multi_length_floquet(p: List[int], q: List[int], seed: int) -> None:
    assert len(set(p)) == 16 and p[0] == 0
    pinv = [0] * 16
    for u, v in enumerate(p):
        pinv[v] = u
    rng = random.Random(seed)
    for n in (3, 5, 11, 29, 64):
        for _ in range(10):
            u0 = [rng.randrange(16) for _ in range(n)]
            q0 = sum(q[x] for x in u0)
            u = u0[:]
            for _ in range(20):
                u = floquet_step_list(p, u)
                qt = sum(q[x] for x in u)
                assert qt == q0, f"Charge drift at N={n}: {qt} != {q0}"
            for _ in range(20):
                u = floquet_inv_step_list(pinv, u)
            assert u == u0, f"Microscopic reversibility mismatch at N={n}"


def probe_4site_conservation(p: List[int], trials: int = 200) -> int:
    """Quick Fraction nullspace dimension check over 2-block (4-site) densities q(u, v) on N=4 rings."""
    rng = random.Random(99)
    rows: List[List[Fraction]] = []
    for _ in range(trials):
        u0 = [rng.randrange(16) for _ in range(4)]
        u1 = floquet_step_list(p, u0)
        row = [Fraction(0)] * 255
        for i in range(4):
            pair0 = (u0[i] << 4) | u0[(i + 1) % 4]
            pair1 = (u1[i] << 4) | u1[(i + 1) % 4]
            if pair1 > 0:
                row[pair1 - 1] += 1
            if pair0 > 0:
                row[pair0 - 1] -= 1
        if any(row):
            rows.append(row)
    return 255 - fraction_rank(rows, 255)


def verify_fir_eval(fir_path: Path, expected: List[int], subcmd: str = "eval-fir") -> None:
    out = subprocess.check_output([str(PHY_BIN), subcmd, str(fir_path)], text=True).strip()
    got = [int(tok) for tok in out.split()]
    assert got == expected, f".fir ({subcmd}) mismatch: got len={len(got)} vs expected len={len(expected)}"


def verify_q4_continuity_exact(p: List[int], q4: List[int]) -> int:
    f3 = [[[step3(p, a, b, c) for c in range(16)] for b in range(16)] for a in range(16)]
    f000 = f3[0][0][0]
    j4 = [
        [
            [
                (q4[(f000 << 4) | f3[0][0][a]] - q4[0])
                + (q4[(f3[0][0][a] << 4) | f3[0][a][b]] - q4[a])
                + (q4[(f3[0][a][b] << 4) | f3[a][b][c]] - q4[(a << 4) | b])
                for c in range(16)
            ]
            for b in range(16)
        ]
        for a in range(16)
    ]
    j4_nz = sum(1 for a in range(16) for b in range(16) for c in range(16) if j4[a][b][c] != 0)
    for a in range(16):
        for b in range(16):
            for c in range(16):
                u1 = f3[a][b][c]
                jabc = j4[a][b][c]
                for d in range(16):
                    u2 = f3[b][c][d]
                    lhs = q4[(u1 << 4) | u2] - q4[(b << 4) | c]
                    rhs = j4[b][c][d] - jabc
                    if lhs != rhs:
                        raise AssertionError(
                            f"Q4 continuity mismatch at ({a},{b},{c},{d}): {lhs} != {rhs}"
                        )
    return j4_nz


def verify_q4_nontrivial_fraction(q4: List[int]) -> bool:
    base: List[List[Fraction]] = []
    for k in range(1, 16):
        row = [Fraction(0)] * 256
        for u in range(16):
            for v in range(16):
                row[(u << 4) | v] = Fraction((1 if v == k else 0) - (1 if u == k else 0))
        base.append(row)
    for k in range(1, 16):
        row = [Fraction(0)] * 256
        for v in range(16):
            row[(k << 4) | v] = Fraction(1)
        base.append(row)
    for h1 in (1, 2, 3):
        for h2 in (1, 2, 3):
            row = [Fraction(0)] * 256
            for u in range(16):
                if hi2(u) != h1:
                    continue
                for v in range(16):
                    if hi2(v) == h2:
                        row[(u << 4) | v] = Fraction(1)
            base.append(row)
    for l1 in (1, 2, 3):
        for l2 in (1, 2, 3):
            row = [Fraction(0)] * 256
            for u in range(16):
                if lo2(u) != l1:
                    continue
                for v in range(16):
                    if lo2(v) == l2:
                        row[(u << 4) | v] = Fraction(1)
            base.append(row)
    base.append([Fraction(1)] * 256)
    r0 = fraction_rank(base, 256)
    r1 = fraction_rank(base + [[Fraction(x) for x in q4]], 256)
    return r1 > r0


def verify_multi_length_floquet_q4(p: List[int], q4: List[int], seed: int) -> None:
    assert len(set(p)) == 16 and p[0] == 0
    pinv = [0] * 16
    for u, v in enumerate(p):
        pinv[v] = u
    rng = random.Random(seed)
    for n in (3, 5, 11, 29, 64):
        for _ in range(10):
            u0 = [rng.randrange(16) for _ in range(n)]
            q0 = sum(q4[(u0[i] << 4) | u0[(i + 1) % n]] for i in range(n))
            u = u0[:]
            for _ in range(20):
                u = floquet_step_list(p, u)
                qt = sum(q4[(u[i] << 4) | u[(i + 1) % n]] for i in range(n))
                assert qt == q0, f"Q4 charge drift at N={n}: {qt} != {q0}"
            for _ in range(20):
                u = floquet_inv_step_list(pinv, u)
            assert u == u0, f"Q4 microscopic reversibility mismatch at N={n}"


def pack3(l: int, r: int) -> int:
    return ((l & 7) << 3) | (r & 7)


def hi3(u: int) -> int:
    return (u >> 3) & 7


def lo3(u: int) -> int:
    return u & 7


def verify_yb8_exact(r: List[int]) -> Tuple[int, int]:
    assert len(r) == 64 and len(set(r)) == 64, "R8 is not a bijection on Fin 64"
    for a in range(8):
        for b in range(8):
            ab = r[pack3(a, b)]
            a1, b1 = hi3(ab), lo3(ab)
            for c in range(8):
                bc = r[pack3(b1, c)]
                b2, c2 = hi3(bc), lo3(bc)
                ab2 = r[pack3(a1, b2)]
                lhs = (hi3(ab2), lo3(ab2), c2)

                bc0 = r[pack3(b, c)]
                b0, c0 = hi3(bc0), lo3(bc0)
                ab0 = r[pack3(a, b0)]
                a0, b00 = hi3(ab0), lo3(ab0)
                bc1 = r[pack3(b00, c0)]
                rhs = (a0, hi3(bc1), lo3(bc1))
                assert lhs == rhs, f"YB8 failure at ({a},{b},{c}): {lhs} != {rhs}"

    # Non-degenerate
    for a in range(8):
        assert {hi3(r[pack3(a, b)]) for b in range(8)} == set(range(8))
    for b in range(8):
        assert {lo3(r[pack3(a, b)]) for a in range(8)} == set(range(8))

    # Non-involutive
    assert any(r[r[u]] != u for u in range(64)), "R8 is involutive!"

    # Indecomposable (undirected orbit of 0)
    reached = {0}
    changed = True
    while changed:
        changed = False
        cur = list(reached)
        for u in cur:
            for x in range(8):
                for nxt in (hi3(r[pack3(x, u)]), lo3(r[pack3(u, x)])):
                    if nxt not in reached:
                        reached.add(nxt)
                        changed = True
                for v in range(8):
                    if hi3(r[pack3(x, v)]) == u or lo3(r[pack3(v, x)]) == u:
                        if v not in reached:
                            reached.add(v)
                            changed = True
    assert len(reached) == 8, f"R8 is decomposable (orbit={reached})"

    # Non-Lyubashenko
    is_lyub = all(
        hi3(r[pack3(a, b)]) == hi3(r[pack3(0, b)]) and lo3(r[pack3(a, b)]) == lo3(r[pack3(a, 0)])
        for a in range(8)
        for b in range(8)
    )
    assert not is_lyub, "R8 is Lyubashenko!"

    return anf_degree6(r), perm_order(r)


def anf_degree6(tbl: List[int]) -> int:
    max_d = 0
    for bit in range(6):
        f = [(tbl[x] >> bit) & 1 for x in range(64)]
        for v in range(6):
            for mask in range(64):
                if mask & (1 << v):
                    f[mask] ^= f[mask ^ (1 << v)]
        for mask in range(64):
            if f[mask]:
                max_d = max(max_d, bin(mask).count("1"))
    return max_d


def perm_order(tbl: List[int]) -> int:
    n = len(tbl)
    vis = [False] * n
    order = 1
    for i in range(n):
        if not vis[i]:
            clen = 0
            cur = i
            while not vis[cur]:
                vis[cur] = True
                cur = tbl[cur]
                clen += 1
            order = math.lcm(order, clen)
    return order


def _mat4_mul(a: List[List[int]], b: List[List[int]]) -> List[List[int]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(4)) for j in range(4)]
        for i in range(4)
    ]


def verify_dual4_exact(p: List[int], d4: dict) -> None:
    assert len(p) == 16 and len(set(p)) == 16 and p[0] == 0
    assert anf_degree(p) == d4["anf_deg"] >= 2
    assert len({(a, hi2(p[pack2(a, b)])) for a in range(4) for b in range(4)}) == 16
    assert len({(lo2(p[pack2(a, b)]), b) for a in range(4) for b in range(4)}) == 16

    m_plus = [
        [sum(1 for b in range(4) if lo2(p[pack2(a, b)]) == d) for d in range(4)]
        for a in range(4)
    ]
    m_minus = [
        [sum(1 for a in range(4) if hi2(p[pack2(a, b)]) == c) for c in range(4)]
        for b in range(4)
    ]
    assert m_plus == d4["m_plus"] and m_minus == d4["m_minus"]
    mp2 = _mat4_mul(m_plus, m_plus)
    mp4 = _mat4_mul(mp2, mp2)
    mm2 = _mat4_mul(m_minus, m_minus)
    mm4 = _mat4_mul(mm2, mm2)
    assert all(mp4[i][j] > 0 and mm4[i][j] > 0 for i in range(4) for j in range(4))

    # Exact Fraction proof: mixing Dual-Unitary gate is ergodic => 0-dim nontrivial Q2 nullspace
    f00 = [step3(p, 0, 0, x) for x in range(16)]
    f0 = [[step3(p, 0, x, y) for y in range(16)] for x in range(16)]
    pivs: Dict[int, List[Fraction]] = {}
    for a in range(16):
        for b in range(16):
            for c in range(16):
                lhs_s = step3(p, a, b, c)
                v = [Fraction(0)] * 15
                for idx_p, sgn in (
                    (lhs_s, 1),
                    (f00[b], -1),
                    (f0[b][c], -1),
                    (f00[a], 1),
                    (f0[a][b], 1),
                    (a, -1),
                ):
                    if idx_p > 0:
                        v[idx_p - 1] += sgn
                if not any(v):
                    continue
                for col in sorted(pivs.keys()):
                    if v[col] != 0:
                        fac = v[col]
                        pr = pivs[col]
                        for j in range(col, 15):
                            v[j] -= fac * pr[j]
                lead = next((j for j, x in enumerate(v) if x != 0), None)
                if lead is not None:
                    inv = Fraction(1) / v[lead]
                    pivs[lead] = [x * inv for x in v]
    assert 15 - len(pivs) == 0, f"Expected 0-dim conserved Q2 nullspace for mixing dual-unitary gate, got {15 - len(pivs)}"
    assert d4["rank_ent"] == 0


def verify_ame8_exact(r: List[int], a8: dict) -> int:
    assert len(r) == 64 and len(set(r)) == 64
    deg = anf_degree6(r)
    assert deg == a8["anf_deg"] >= 2
    assert perm_order(r) == a8["order"]
    for a in range(8):
        assert len({hi3(r[pack3(a, b)]) for b in range(8)}) == 8
        assert len({lo3(r[pack3(a, b)]) for b in range(8)}) == 8
    for b in range(8):
        assert len({lo3(r[pack3(a, b)]) for a in range(8)}) == 8
        assert len({hi3(r[pack3(a, b)]) for a in range(8)}) == 8
    return deg


def _frac_apply(t: List[int], st: int, pos: int) -> int:
    sh = pos * 2
    a = (st >> sh) & 3
    b = (st >> (sh + 2)) & 3
    c = (st >> (sh + 4)) & 3
    o = t[(a << 4) | (b << 2) | c]
    oa, ob, oc = (o >> 4) & 3, (o >> 2) & 3, o & 3
    nw = oa | (ob << 2) | (oc << 4)
    return (st & ~(0x3F << sh)) | (nw << sh)


def _frac_qd(st: int, n_sites: int) -> Tuple[int, int]:
    q, d = 0, 0
    for i in range(n_sites):
        v = (st >> (i * 2)) & 3
        q += v
        d += i * v
    return q, d


def verify_fracton_exact(t: List[int], fr: dict) -> int:
    assert len(t) == 64 and len(set(t)) == 64
    tinv = [0] * 64
    for u, v in enumerate(t):
        tinv[v] = u
    deg = anf_degree6(t)
    assert deg == fr["anf_deg"] >= 3
    act = sum(1 for u in range(64) if t[u] != u)
    assert act == fr["active_entries"] > 0

    for u in range(64):
        a, b, c = (u >> 4) & 3, (u >> 2) & 3, u & 3
        o = t[u]
        oa, ob, oc = (o >> 4) & 3, (o >> 2) & 3, o & 3
        phi = oa - a
        assert (ob - b) == -(2 * phi) and (oc - c) == phi

    # Exhaustive L=6 (4096 states) Q & D conservation + largest-sector BFS
    bq, bd = fr["hsf_target_q"], fr["hsf_target_d"]
    sec_states: List[int] = []
    for st in range(4096):
        q0, d0 = _frac_qd(st, 6)
        if q0 == bq and d0 == bd and len(sec_states) < 128:
            sec_states.append(st)
        for pos in range(4):
            q1, d1 = _frac_qd(_frac_apply(t, st, pos), 6)
            assert (q1, d1) == (q0, d0)
    assert len(sec_states) == fr["hsf_sector_size"]

    vis = set()
    n_comps = 0
    max_comp = 0
    frozen = 0
    for root in sec_states:
        if root in vis:
            continue
        q = deque([root])
        vis.add(root)
        csz = 0
        while q:
            cur = q.popleft()
            csz += 1
            for pos in range(4):
                for tbl in (t, tinv):
                    nxt = _frac_apply(tbl, cur, pos)
                    if nxt not in vis:
                        vis.add(nxt)
                        q.append(nxt)
        n_comps += 1
        max_comp = max(max_comp, csz)
        if csz == 1:
            frozen += 1
    assert n_comps == fr["hsf_num_comps"]
    assert max_comp == fr["hsf_largest_comp"]
    assert frozen == fr["hsf_frozen"] and n_comps > 1

    # L=32 isolated single charge immobility vs dipole probe
    iso_init = 1 << (16 * 2)
    dip_init = 2 << (16 * 2)
    iso_q0, iso_d0 = _frac_qd(iso_init, 32)
    dip_q0, dip_d0 = _frac_qd(dip_init, 32)
    s_iso, s_dip = iso_init, dip_init
    iso_span, dip_span = 0, 0
    for _ in range(32):
        for layer in range(3):
            for pos in range(layer, 30, 3):
                s_iso = _frac_apply(t, s_iso, pos)
                s_dip = _frac_apply(t, s_dip, pos)
        assert _frac_qd(s_iso, 32) == (iso_q0, iso_d0)
        assert _frac_qd(s_dip, 32) == (dip_q0, dip_d0)
        for site in range(32):
            dist = abs(site - 16)
            if ((s_iso >> (site * 2)) & 3) != ((iso_init >> (site * 2)) & 3):
                iso_span = max(iso_span, dist)
            if ((s_dip >> (site * 2)) & 3) != ((dip_init >> (site * 2)) & 3):
                dip_span = max(dip_span, dist)
    assert iso_span == 0 == fr["iso_span"] and dip_span == fr["dip_span"] > 0
    return deg


def verify_biquandle8_exact(r: List[int], bq: dict) -> int:
    deg, order = verify_yb8_exact(r)
    assert deg == bq["anf_deg"] >= 2 and order == bq["order"]
    assert all(r[pack3(x, x)] == pack3(x, x) for x in range(8))

    sigs = [
        (tuple(hi3(r[pack3(x, k)]) for k in range(8)), tuple(lo3(r[pack3(k, x)]) for k in range(8)))
        for x in range(8)
    ]
    assert len(set(sigs)) == 8 == bq["retract_classes"]

    rinv = [0] * 64
    for u, v in enumerate(r):
        rinv[v] = u

    def fix_pow(k: int) -> int:
        cnt = 0
        for u in range(64):
            x = u
            for _ in range(k):
                x = r[x]
            if x == u:
                cnt += 1
        return cnt

    col_fig8 = 0
    for a in range(8):
        for b in range(8):
            for c in range(8):
                x, y, z = a, b, c
                for _ in range(2):
                    ab = r[pack3(x, y)]
                    x, y = hi3(ab), lo3(ab)
                    bc = rinv[pack3(y, z)]
                    y, z = hi3(bc), lo3(bc)
                if (x, y, z) == (a, b, c):
                    col_fig8 += 1

    assert fix_pow(1) == bq["col_unknot"]
    assert fix_pow(2) == bq["col_hopf"]
    assert fix_pow(3) == bq["col_trefoil"]
    assert fix_pow(5) == bq["col_torus5"]
    assert fix_pow(7) == bq["col_torus7"]
    assert col_fig8 == bq["col_fig8"]
    return deg


def _tetra4_apply(t: List[int], st: int, i: int, j: int, k: int) -> int:
    si, sj, sk = 2 * (5 - i), 2 * (5 - j), 2 * (5 - k)
    o = t[(((st >> si) & 3) << 4) | (((st >> sj) & 3) << 2) | ((st >> sk) & 3)]
    clr = ~((3 << si) | (3 << sj) | (3 << sk)) & 0xFFF
    return (st & clr) | (((o >> 4) & 3) << si) | (((o >> 2) & 3) << sj) | ((o & 3) << sk)


def verify_tetra4_exact(t: List[int], t4: dict) -> int:
    assert len(t) == 64 and len(set(t)) == 64
    deg = anf_degree6(t)
    assert deg == t4["anf_deg"] >= 3
    assert perm_order(t) == t4["order"]
    assert any(t[t[u]] != u for u in range(64))

    # 3-channel middle coupling witnesses
    left_c = any(
        ((t[(a << 4) | (b << 2) | c] >> 2) & 3) != ((t[(b << 2) | c] >> 2) & 3)
        for a in range(1, 4)
        for b in range(4)
        for c in range(4)
    )
    right_c = any(
        ((t[(a << 4) | (b << 2) | c] >> 2) & 3) != ((t[(a << 4) | (b << 2)] >> 2) & 3)
        for a in range(4)
        for b in range(4)
        for c in range(1, 4)
    )
    assert left_c and right_c

    for s in range(4096):
        lhs = _tetra4_apply(
            t, _tetra4_apply(t, _tetra4_apply(t, _tetra4_apply(t, s, 2, 4, 5), 1, 3, 5), 0, 3, 4), 0, 1, 2
        )
        rhs = _tetra4_apply(
            t, _tetra4_apply(t, _tetra4_apply(t, _tetra4_apply(t, s, 0, 1, 2), 0, 3, 4), 1, 3, 5), 2, 4, 5
        )
        assert lhs == rhs, f"Tetra4 failure at state {s}"
    return deg


def _anf_degree_bits(table: List[int], n_bits: int) -> int:
    n = 1 << n_bits
    max_deg = 0
    for bit in range(n_bits):
        anf = [(table[u] >> bit) & 1 for u in range(n)]
        for i in range(n_bits):
            step = 1 << i
            for mask in range(n):
                if mask & step:
                    anf[mask] ^= anf[mask ^ step]
        for mask in range(n):
            if anf[mask]:
                d = bin(mask).count("1")
                if d > max_deg:
                    max_deg = d
    return max_deg


def verify_ame4_nogo_exact(ng: dict) -> None:
    cls_out = subprocess.check_output([str(PHY_BIN), "classify-ame4"], text=True).strip()
    assert (
        "latin4=576" in cls_out
        and "total=6912" in cls_out
        and "affine=6912" in cls_out
        and "nonlinear=0" in cls_out
        and "bernoulli=6912" in cls_out
    ), cls_out

    # Exhaustive pure-Python enumeration of all 4x4 Latin squares -> OLS -> GF(2)-affinity
    latins: List[Tuple[int, ...]] = []
    grid = [0] * 16
    row_m = [0] * 4
    col_m = [0] * 4

    def bt(pos: int) -> None:
        if pos == 16:
            latins.append(tuple(grid))
            return
        r, c = pos >> 2, pos & 3
        used = row_m[r] | col_m[c]
        for v in range(4):
            bit = 1 << v
            if not (used & bit):
                grid[pos] = v
                row_m[r] |= bit
                col_m[c] |= bit
                bt(pos + 1)
                row_m[r] ^= bit
                col_m[c] ^= bit

    bt(0)
    assert len(latins) == 576 == ng["latin4_count"]

    # Map each value 0..3 of R to its 4 cell indices to test orthogonality in O(1)
    occ_list = []
    for sq in latins:
        occ = ([], [], [], [])
        for idx, v in enumerate(sq):
            occ[v].append(idx)
        occ_list.append(occ)

    total_ols = 0
    affine_ols = 0
    nonlin_ols = 0
    bernoulli_ols = 0
    for r_occ, r_sq in zip(occ_list, latins):
        for c_sq in latins:
            ok = True
            for v in range(4):
                i0, i1, i2, i3 = r_occ[v]
                if ((1 << c_sq[i0]) | (1 << c_sq[i1]) | (1 << c_sq[i2]) | (1 << c_sq[i3])) != 0xF:
                    ok = False
                    break
            if not ok:
                continue
            total_ols += 1
            p0 = (r_sq[0] << 2) | c_sq[0]
            b1 = ((r_sq[1] << 2) | c_sq[1]) ^ p0
            b2 = ((r_sq[2] << 2) | c_sq[2]) ^ p0
            b4 = ((r_sq[4] << 2) | c_sq[4]) ^ p0
            b8 = ((r_sq[8] << 2) | c_sq[8]) ^ p0
            is_aff = True
            for x in range(16):
                px = (r_sq[x] << 2) | c_sq[x]
                lin = (b1 if (x & 1) else 0) ^ (b2 if (x & 2) else 0) ^ (b4 if (x & 4) else 0) ^ (b8 if (x & 8) else 0)
                if (px ^ p0) != lin:
                    is_aff = False
                    break
            if is_aff:
                affine_ols += 1
                bernoulli_ols += 1
            else:
                nonlin_ols += 1

    assert (total_ols, affine_ols, nonlin_ols, bernoulli_ols) == (
        ng["total_ols"],
        ng["affine_ols"],
        ng["nonlin_ols"],
        ng["bernoulli_ols"],
    )
    assert total_ols == 6912 and affine_ols == 6912 and nonlin_ols == 0 and bernoulli_ols == 6912

    sp = ng["sample_p"]
    assert len(sp) == 16 and len(set(sp)) == 16
    assert anf_degree(sp) == 1
    assert len({(a, sp[(a << 2) | b] >> 2) for a in range(4) for b in range(4)}) == 16
    assert len({(a, sp[(a << 2) | b] & 3) for a in range(4) for b in range(4)}) == 16
    assert len({(sp[(b << 2) | a] >> 2, a) for a in range(4) for b in range(4)}) == 16
    assert len({(sp[(b << 2) | a] & 3, a) for a in range(4) for b in range(4)}) == 16
    packed = 0
    for i in range(16):
        packed |= (sp[i] & 0xF) << (i * 4)
    assert f"0x{packed:016x}" == ng["sample_lut"]


def _s4_apply2(t: List[int], st: int, i: int, j: int, k: int, l: int) -> int:
    si, sj, sk, sl = 9 - i, 9 - j, 9 - k, 9 - l
    u = (
        (((st >> si) & 1) << 3)
        | (((st >> sj) & 1) << 2)
        | (((st >> sk) & 1) << 1)
        | ((st >> sl) & 1)
    )
    o = t[u]
    clr = ~((1 << si) | (1 << sj) | (1 << sk) | (1 << sl)) & 0x3FF
    return (
        (st & clr)
        | (((o >> 3) & 1) << si)
        | (((o >> 2) & 1) << sj)
        | (((o >> 1) & 1) << sk)
        | ((o & 1) << sl)
    )


_S4_TUPLES = ((6, 7, 8, 9), (3, 4, 5, 9), (1, 2, 5, 8), (0, 2, 4, 7), (0, 1, 3, 6))
_S4_SHIFTS = tuple(
    (2 * (9 - i), 2 * (9 - j), 2 * (9 - k), 2 * (9 - m)) for i, j, k, m in _S4_TUPLES
)
_S4_KEEPS = tuple(
    0xFFFFF ^ ((3 << s0) | (3 << s1) | (3 << s2) | (3 << s3)) for s0, s1, s2, s3 in _S4_SHIFTS
)


def _s4_viol4(t: List[int]) -> int:
    (
        (s0a, s0b, s0c, s0d),
        (s1a, s1b, s1c, s1d),
        (s2a, s2b, s2c, s2d),
        (s3a, s3b, s3c, s3d),
        (s4a, s4b, s4c, s4d),
    ) = _S4_SHIFTS
    k0, k1, k2, k3, k4 = _S4_KEEPS
    viol = 0
    for st in range(1 << 20):
        o = t[
            (((st >> s0a) & 3) << 6)
            | (((st >> s0b) & 3) << 4)
            | (((st >> s0c) & 3) << 2)
            | ((st >> s0d) & 3)
        ]
        l = (
            (st & k0)
            | (((o >> 6) & 3) << s0a)
            | (((o >> 4) & 3) << s0b)
            | (((o >> 2) & 3) << s0c)
            | ((o & 3) << s0d)
        )
        o = t[
            (((l >> s1a) & 3) << 6)
            | (((l >> s1b) & 3) << 4)
            | (((l >> s1c) & 3) << 2)
            | ((l >> s1d) & 3)
        ]
        l = (
            (l & k1)
            | (((o >> 6) & 3) << s1a)
            | (((o >> 4) & 3) << s1b)
            | (((o >> 2) & 3) << s1c)
            | ((o & 3) << s1d)
        )
        o = t[
            (((l >> s2a) & 3) << 6)
            | (((l >> s2b) & 3) << 4)
            | (((l >> s2c) & 3) << 2)
            | ((l >> s2d) & 3)
        ]
        l = (
            (l & k2)
            | (((o >> 6) & 3) << s2a)
            | (((o >> 4) & 3) << s2b)
            | (((o >> 2) & 3) << s2c)
            | ((o & 3) << s2d)
        )
        o = t[
            (((l >> s3a) & 3) << 6)
            | (((l >> s3b) & 3) << 4)
            | (((l >> s3c) & 3) << 2)
            | ((l >> s3d) & 3)
        ]
        l = (
            (l & k3)
            | (((o >> 6) & 3) << s3a)
            | (((o >> 4) & 3) << s3b)
            | (((o >> 2) & 3) << s3c)
            | ((o & 3) << s3d)
        )
        o = t[
            (((l >> s4a) & 3) << 6)
            | (((l >> s4b) & 3) << 4)
            | (((l >> s4c) & 3) << 2)
            | ((l >> s4d) & 3)
        ]
        l = (
            (l & k4)
            | (((o >> 6) & 3) << s4a)
            | (((o >> 4) & 3) << s4b)
            | (((o >> 2) & 3) << s4c)
            | ((o & 3) << s4d)
        )

        o = t[
            (((st >> s4a) & 3) << 6)
            | (((st >> s4b) & 3) << 4)
            | (((st >> s4c) & 3) << 2)
            | ((st >> s4d) & 3)
        ]
        r = (
            (st & k4)
            | (((o >> 6) & 3) << s4a)
            | (((o >> 4) & 3) << s4b)
            | (((o >> 2) & 3) << s4c)
            | ((o & 3) << s4d)
        )
        o = t[
            (((r >> s3a) & 3) << 6)
            | (((r >> s3b) & 3) << 4)
            | (((r >> s3c) & 3) << 2)
            | ((r >> s3d) & 3)
        ]
        r = (
            (r & k3)
            | (((o >> 6) & 3) << s3a)
            | (((o >> 4) & 3) << s3b)
            | (((o >> 2) & 3) << s3c)
            | ((o & 3) << s3d)
        )
        o = t[
            (((r >> s2a) & 3) << 6)
            | (((r >> s2b) & 3) << 4)
            | (((r >> s2c) & 3) << 2)
            | ((r >> s2d) & 3)
        ]
        r = (
            (r & k2)
            | (((o >> 6) & 3) << s2a)
            | (((o >> 4) & 3) << s2b)
            | (((o >> 2) & 3) << s2c)
            | ((o & 3) << s2d)
        )
        o = t[
            (((r >> s1a) & 3) << 6)
            | (((r >> s1b) & 3) << 4)
            | (((r >> s1c) & 3) << 2)
            | ((r >> s1d) & 3)
        ]
        r = (
            (r & k1)
            | (((o >> 6) & 3) << s1a)
            | (((o >> 4) & 3) << s1b)
            | (((o >> 2) & 3) << s1c)
            | ((o & 3) << s1d)
        )
        o = t[
            (((r >> s0a) & 3) << 6)
            | (((r >> s0b) & 3) << 4)
            | (((r >> s0c) & 3) << 2)
            | ((r >> s0d) & 3)
        ]
        r = (
            (r & k0)
            | (((o >> 6) & 3) << s0a)
            | (((o >> 4) & 3) << s0b)
            | (((o >> 2) & 3) << s0c)
            | ((o & 3) << s0d)
        )
        if l != r:
            viol += 1
    return viol


def verify_simplex4_exact(s4: dict) -> int:
    ccx02 = s4["ccx02"]
    assert len(ccx02) == 16 and len(set(ccx02)) == 16
    assert anf_degree(ccx02) == 2
    ccx13 = [u ^ ((((u >> 2) & 1) & (u & 1)) << 1) for u in range(16)]
    cccx1 = [u ^ ((((u >> 3) & 1) & ((u >> 1) & 1) & (u & 1)) << 2) for u in range(16)]
    assert anf_degree(cccx1) == 3

    def viol2(tbl: List[int]) -> int:
        cnt = 0
        for st in range(1024):
            l = st
            for qu in ((6, 7, 8, 9), (3, 4, 5, 9), (1, 2, 5, 8), (0, 2, 4, 7), (0, 1, 3, 6)):
                l = _s4_apply2(tbl, l, *qu)
            r = st
            for qu in ((0, 1, 3, 6), (0, 2, 4, 7), (1, 2, 5, 8), (3, 4, 5, 9), (6, 7, 8, 9)):
                r = _s4_apply2(tbl, r, *qu)
            if l != r:
                cnt += 1
        return cnt

    assert viol2(ccx02) == 0 == s4["viol2_ccx"]
    assert viol2(ccx13) == 0
    assert viol2(cccx1) == s4["viol2_cccx"] > 0

    s = s4["s"]
    sc = s4["s_coupled"]
    assert len(s) == 256 and len(set(s)) == 256
    assert len(sc) == 256 and len(set(sc)) == 256
    deg = _anf_degree_bits(s, 8)
    assert deg == s4["anf_deg"] >= 3
    assert perm_order(s) == s4["order4"] > 2
    assert any(s[s[u]] != u for u in range(256)) and s4["is_involutive"] == 0

    c0 = any(
        ((sc[(a << 6) | (b << 4) | (c << 2) | d] >> 4) & 3)
        != ((sc[(b << 4) | (c << 2) | d] >> 4) & 3)
        for a in range(1, 4)
        for b in range(4)
        for c in range(4)
        for d in range(4)
    )
    c1 = any(
        ((sc[(a << 6) | (b << 4) | (c << 2) | d] >> 2) & 3)
        != ((sc[(a << 6) | (c << 2) | d] >> 2) & 3)
        for a in range(4)
        for b in range(1, 4)
        for c in range(4)
        for d in range(4)
    )
    c2 = any(
        ((sc[(a << 6) | (b << 4) | (c << 2) | d] >> 4) & 3)
        != ((sc[(a << 6) | (b << 4) | d] >> 4) & 3)
        for a in range(4)
        for b in range(4)
        for c in range(1, 4)
        for d in range(4)
    )
    c3 = any(
        ((sc[(a << 6) | (b << 4) | (c << 2) | d] >> 2) & 3)
        != ((sc[(a << 6) | (b << 4) | (c << 2)] >> 2) & 3)
        for a in range(4)
        for b in range(4)
        for c in range(4)
        for d in range(1, 4)
    )
    assert c0 and c1 and c2 and c3 and s4["is_4coupled"] == 1

    assert _s4_viol4(s) == 0 == s4["viol4"]
    assert _s4_viol4(sc) == 0 == s4["viol4_coupled"]
    return deg


def verify_chiral6_exact(ch: dict) -> None:
    r = ch["r"]
    phi = ch["phi"]
    assert len(r) == 36 and len(set(r)) == 36
    assert len(phi) == 36 and all(0 <= v < 4 for v in phi)
    rinv = [0] * 36
    for u, v in enumerate(r):
        rinv[v] = u
    pa = [u // 6 for u in range(36)]
    pb = [u % 6 for u in range(36)]

    # R1 idempotence & 216-triple Yang-Baxter Equation
    for x in range(6):
        assert r[x * 6 + x] == x * 6 + x
        assert (phi[x * 6 + x] & 3) == 0
    for a in range(6):
        for b in range(6):
            ab_i = a * 6 + b
            ab = r[ab_i]
            a1, b1 = pa[ab], pb[ab]
            bc0_i_base = b * 6
            for c in range(6):
                bc_i = b1 * 6 + c
                bc = r[bc_i]
                b2, c2 = pa[bc], pb[bc]
                ab2_i = a1 * 6 + b2
                ab2 = r[ab2_i]
                a2, b3 = pa[ab2], pb[ab2]

                bc0_i = bc0_i_base + c
                bc0 = r[bc0_i]
                b0, c0 = pa[bc0], pb[bc0]
                ab0_i = a * 6 + b0
                ab0 = r[ab0_i]
                a0, b00 = pa[ab0], pb[ab0]
                bc1_i = b00 * 6 + c0
                bc1 = r[bc1_i]
                b01, c01 = pa[bc1], pb[bc1]
                assert (a2, b3, c2) == (a0, b01, c01)

                lhs = (phi[ab_i] + phi[bc_i] + phi[ab2_i]) % 4
                rhs = (phi[bc0_i] + phi[ab0_i] + phi[bc1_i]) % 4
                assert lhs == rhs
    assert ch["cocycle_viol"] == 0

    # Exhaustively verify phi is NOT any Z_4 coboundary over all 4^6 = 4096 1-cochains
    for mask in range(4**6):
        f = [(mask >> (2 * i)) & 3 for i in range(6)]
        if all(((f[pa[r[u]]] + f[pb[r[u]]] - f[pa[u]] - f[pb[u]]) % 4) == phi[u] for u in range(36)):
            raise AssertionError("phi unexpectedly turned out to be a coboundary!")

    unk = [0, 0, 0, 0]
    hopf = [0, 0, 0, 0]
    tref_l = [0, 0, 0, 0]
    tref_r = [0, 0, 0, 0]
    fig8 = [0, 0, 0, 0]
    for u in range(36):
        x, s = u, 0
        for step in (1, 2, 3):
            s = (s + phi[x]) % 4
            x = r[x]
            if x == u:
                if step == 1:
                    unk[s] += 1
                elif step == 2:
                    hopf[s] += 1
                else:
                    tref_l[s] += 1
        x, s = u, 0
        for _ in range(3):
            y = rinv[x]
            s = (s - phi[y]) % 4
            x = y
        if x == u:
            tref_r[s] += 1
    for a in range(6):
        for b in range(6):
            for c in range(6):
                x, y, z, s = a, b, c, 0
                for _ in range(2):
                    u1 = x * 6 + y
                    s = (s + phi[u1]) % 4
                    o1 = r[u1]
                    x, y = pa[o1], pb[o1]
                    u2 = y * 6 + z
                    o2 = rinv[u2]
                    s = (s - phi[o2]) % 4
                    y, z = pa[o2], pb[o2]
                if (x, y, z) == (a, b, c):
                    fig8[s] += 1

    assert unk == ch["unknot_w"] == [6, 0, 0, 0]
    assert hopf == ch["hopf_w"] == [6, 0, 6, 0]
    assert tref_l == ch["tref_l_w"] == [6, 24, 0, 0]
    assert tref_r == ch["tref_r_w"] == [6, 0, 0, 24]
    assert fig8 == ch["fig8_w"] == [6, 0, 0, 0]
    assert tref_l != tref_r and ch["is_chiral"] == 1 and ch["is_non_coboundary"] == 1


def verify_qscar_exact(qs: dict) -> None:
    wre = qs["w_re"]
    wim = qs["w_im"]
    w = [[complex(wre[i][j], wim[i][j]) for j in range(4)] for i in range(4)]

    # Exact integer unitarity (W† W = 4I), Dual-Unitarity (W~† W~ = 4I), full superposition
    for i in range(4):
        for j in range(4):
            assert wre[i][j] != 0 or wim[i][j] != 0
    for j in range(4):
        for k in range(4):
            sr = sum(wre[i][j] * wre[i][k] + wim[i][j] * wim[i][k] for i in range(4))
            si = sum(wre[i][j] * wim[i][k] - wim[i][j] * wre[i][k] for i in range(4))
            assert sr == (4 if j == k else 0) and si == 0
    for b1 in range(2):
        for d1 in range(2):
            for b2 in range(2):
                for d2 in range(2):
                    col1 = (b1 << 1) | d1
                    col2 = (b2 << 1) | d2
                    sr, si = 0, 0
                    for a in range(2):
                        for c in range(2):
                            r1, c1 = (c << 1) | d1, (a << 1) | b1
                            r2, c2 = (c << 1) | d2, (a << 1) | b2
                            sr += wre[r1][c1] * wre[r2][c2] + wim[r1][c1] * wim[r2][c2]
                            si += wre[r1][c1] * wim[r2][c2] - wim[r1][c1] * wre[r2][c2]
                    assert sr == (4 if col1 == col2 else 0) and si == 0
    assert qs["du_unitary"] == 1 and qs["du_dual"] == 1 and qs["du_full_sup"] == 1

    def apply_2q(st_r: List[int], st_i: List[int], q0: int, q1: int) -> Tuple[List[int], List[int]]:
        nr = [0] * 64
        ni = [0] * 64
        for s in range(64):
            ar, ai = st_r[s], st_i[s]
            if ar == 0 and ai == 0:
                continue
            iin = (((s >> q0) & 1) << 1) | ((s >> q1) & 1)
            base = s & ~((1 << q0) | (1 << q1))
            for o in range(4):
                wr, wi = wre[o][iin], wim[o][iin]
                ns = base | (((o >> 1) & 1) << q0) | ((o & 1) << q1)
                nr[ns] += wr * ar - wi * ai
                ni[ns] += wr * ai + wi * ar
        return nr, ni

    def apply_scar4(st_r: List[int], st_i: List[int], q0: int, q1: int, q2: int, q3: int) -> Tuple[List[int], List[int]]:
        nr = [0] * 64
        ni = [0] * 64
        for s in range(64):
            ar, ai = st_r[s], st_i[s]
            if ar == 0 and ai == 0:
                continue
            if ((s >> q0) & 1) != ((s >> q3) & 1):
                nr[s] += 2 * ar
                ni[s] += 2 * ai
            else:
                iin = (((s >> q1) & 1) << 1) | ((s >> q2) & 1)
                base = s & ~((1 << q1) | (1 << q2))
                for o in range(4):
                    wr, wi = wre[o][iin], wim[o][iin]
                    ns = base | (((o >> 1) & 1) << q1) | ((o & 1) << q2)
                    nr[ns] += wr * ar - wi * ai
                    ni[ns] += wr * ai + wi * ar
        return nr, ni

    def metrics(st_r: List[int], st_i: List[int], ref_s0: int) -> Tuple[Fraction, Fraction, int]:
        norm2 = sum(st_r[s] * st_r[s] + st_i[s] * st_i[s] for s in range(64))
        sup = sum(1 for s in range(64) if st_r[s] != 0 or st_i[s] != 0)
        ret = Fraction(st_r[ref_s0] * st_r[ref_s0] + st_i[ref_s0] * st_i[ref_s0], norm2)
        tr_rho2 = 0
        for a1 in range(8):
            for a2 in range(8):
                mr = sum(st_r[(b << 3) | a1] * st_r[(b << 3) | a2] + st_i[(b << 3) | a1] * st_i[(b << 3) | a2] for b in range(8))
                mi = sum(st_i[(b << 3) | a1] * st_r[(b << 3) | a2] - st_r[(b << 3) | a1] * st_i[(b << 3) | a2] for b in range(8))
                tr_rho2 += mr * mr + mi * mi
        return Fraction(tr_rho2, norm2 * norm2), ret, sup

    # 1. Dual-Unitary brick-wall on L=6
    sr, si = [0] * 64, [0] * 64
    sr[0] = 1
    bw_pairs = ((0, 1), (2, 3), (4, 5), (1, 2), (3, 4), (5, 0))
    du_pur = []
    du_sup = []
    for _ in range(4):
        for q0, q1 in bw_pairs:
            sr, si = apply_2q(sr, si, q0, q1)
        p, _, sup = metrics(sr, si, 0)
        du_pur.append(p)
        du_sup.append(sup)
    assert du_sup == qs["du_sup"]
    assert du_pur[1] == Fraction(qs["du_purity_t2_num"], qs["du_purity_t2_den"]) == Fraction(1, 8)

    # 2. Scar + X^6 DTC on Neel |010101> (0x15)
    sr, si = [0] * 64, [0] * 64
    sr[0x15] = 1
    neel_sup, neel_ret = [], []
    for _ in range(4):
        for i in range(6):
            sr, si = apply_scar4(sr, si, i, (i + 1) % 6, (i + 2) % 6, (i + 3) % 6)
        sr = [sr[s ^ 63] for s in range(64)]
        si = [si[s ^ 63] for s in range(64)]
        _, ret, sup = metrics(sr, si, 0x15)
        neel_sup.append(sup)
        neel_ret.append(ret)
    assert neel_sup == qs["neel_sup"] == [1, 1, 1, 1]
    assert neel_ret == [Fraction(qs["neel_ret_num"][t], qs["neel_ret_den"][t]) for t in range(4)]
    assert neel_ret == [Fraction(0), Fraction(1), Fraction(0), Fraction(1)]
    assert qs["scar_period"] == 2

    # 3. Single-defect thermalization from |010111> (0x17)
    sr, si = [0] * 64, [0] * 64
    sr[0x17] = 1
    th_pur, th_sup = [], []
    for _ in range(4):
        for i in range(6):
            sr, si = apply_scar4(sr, si, i, (i + 1) % 6, (i + 2) % 6, (i + 3) % 6)
        sr = [sr[s ^ 63] for s in range(64)]
        si = [si[s ^ 63] for s in range(64)]
        p, _, sup = metrics(sr, si, 0x17)
        th_pur.append(p)
        th_sup.append(sup)
    assert th_sup == qs["therm_sup"] and th_sup[3] >= 32
    exp_th = Fraction(qs["therm_purity_t4_num"], qs["therm_purity_t4_den"])
    assert th_pur[3] == exp_th < Fraction(1, 3)


def verify_fracton_comp_exact(fc: dict) -> None:
    t_s = fc["t_sample"]
    t_tf = fc["t_toff"]
    assert len(t_s) == 64 and len(set(t_s)) == 64
    assert len(t_tf) == 64 and len(set(t_tf)) == 64
    for u in range(64):
        assert t_tf[t_tf[u]] == u
        a, b, c = (u >> 4) & 3, (u >> 2) & 3, u & 3
        oa, ob, oc = (t_tf[u] >> 4) & 3, (t_tf[u] >> 2) & 3, t_tf[u] & 3
        phi = oa - a
        assert (ob - b) == -(2 * phi) and (oc - c) == phi
    assert fc["is_involutive"] == 1 and fc["laplacian_ok"] == 1

    for ca in (0, 1):
        for cc in (0, 1):
            for tgt in (0, 1):
                in_u = 34 if (ca and cc and tgt) else ((ca << 4) | ((0 if tgt else 2) << 2) | cc)
                out_u = t_tf[in_u]
                got_tgt = 1 if ((out_u >> 2) & 3) == 0 else 0
                assert got_tgt == (tgt ^ (ca & cc))
    assert fc["ccx_ok"] == 1

    t_s_inv = [0] * 64
    for u, v in enumerate(t_s):
        t_s_inv[v] = u

    sec_l, comp_l, larg_l, froz_l = [], [], [], []
    for n_sites in (6, 7, 8):
        total = 1 << (2 * n_sites)
        sectors: Dict[Tuple[int, int], List[int]] = {}
        for st in range(total):
            qd = _frac_qd(st, n_sites)
            sectors.setdefault(qd, []).append(st)
        best_qd = max(sectors, key=lambda k: len(sectors[k]))
        sec_states = sectors[best_qd]
        assert len(sec_states) == fc["sec"][n_sites - 6]
        vis = bytearray(total)
        n_comps = 0
        max_comp = 0
        frozen = 0
        for root in sec_states:
            if vis[root]:
                continue
            q = deque([root])
            vis[root] = 1
            csz = 0
            while q:
                cur = q.popleft()
                csz += 1
                for pos in range(n_sites - 2):
                    for tbl in (t_s, t_s_inv):
                        nxt = _frac_apply(tbl, cur, pos)
                        if not vis[nxt]:
                            vis[nxt] = 1
                            q.append(nxt)
            n_comps += 1
            if csz > max_comp:
                max_comp = csz
            if csz == 1:
                frozen += 1
        sec_l.append(len(sec_states))
        comp_l.append(n_comps)
        larg_l.append(max_comp)
        froz_l.append(frozen)
    assert sec_l == fc["sec"]
    assert comp_l == fc["comps"]
    assert larg_l == fc["largest"]
    assert froz_l == fc["frozen"]
    assert larg_l[1] * sec_l[0] < larg_l[0] * sec_l[1]
    assert larg_l[2] * sec_l[1] < larg_l[1] * sec_l[2]
    assert fc["scaling_strict"] == 1


def verify_cross_ecosystem_bridge(fir_paths: List[Tuple[str, Path]]) -> None:
    qcc = ROOT.parent / "QUA" / "qcc" / "build" / "qcc"
    ele = ROOT.parent / "ELE" / "build" / "ele_bootstrap"
    pdk = ROOT.parent / "ELE" / "build" / "pdk_sky130.flow"
def _run_qcc_one(qcc: Path, tag: str, fp: Path) -> None:
    kmb = BUILD / f"qcc_{tag}.kmb"
    subprocess.run([str(qcc), "opt", str(fp), "-o", str(kmb)], capture_output=True, text=True, check=True)
    vres = subprocess.run(
        [str(qcc), "verify", str(kmb), "--ref", str(fp), "--states", "4096"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "PASS" in vres.stdout or "val_err=0" in vres.stdout, f"qcc verify failed on {tag}: {vres.stdout}"


def _run_ele_one(ele: Path, pdk: Path, tag: str, fp: Path) -> None:
    sv = BUILD / f"ele_{tag}.sv"
    gds = BUILD / f"ele_{tag}.gds"
    eres = subprocess.run(
        [str(ele), "--fir-synth", str(fp), "-o", str(sv), "--pdk", str(pdk), "--gds-out", str(gds), "--verify-fir"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "DRC=0" in eres.stdout and "LVS=0" in eres.stdout, f"ELE synth failed on {tag}: {eres.stdout}"


def verify_cross_ecosystem_bridge(
    fir_paths: List[Tuple[str, Path]], pool: Optional[concurrent.futures.ThreadPoolExecutor] = None
) -> List[str]:
    qcc = ROOT.parent / "QUA" / "qcc" / "build" / "qcc"
    ele = ROOT.parent / "ELE" / "build" / "ele_bootstrap"
    pdk = ROOT.parent / "ELE" / "build" / "pdk_sky130.flow"
    lines: List[str] = []
    if pool is None:
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 4)) as local_pool:
            return verify_cross_ecosystem_bridge(fir_paths, pool=local_pool)

    qcc_futs = [pool.submit(_run_qcc_one, qcc, tag, fp) for tag, fp in fir_paths] if qcc.exists() else []
    ele_futs = (
        [pool.submit(_run_ele_one, ele, pdk, tag, fp) for tag, fp in fir_paths]
        if (ele.exists() and pdk.exists())
        else []
    )
    for fut in qcc_futs:
        fut.result()
    if qcc_futs:
        lines.append(
            f"  [Oracle QUA-qcc]   PASS ({len(fir_paths)}/{len(fir_paths)} .fir -> qcc .kmb 4096-state verify val_err=0)"
        )
    for fut in ele_futs:
        fut.result()
    if ele_futs:
        lines.append(
            f"  [Oracle ELE-ASIC]  PASS ({len(fir_paths)}/{len(fir_paths)} .fir -> SKY130 SV+GDSII DRC=0 LVS=0)"
        )
    return lines


def compile_and_check_lean_cert(lean_path: Path, mod_name: str) -> None:
    olean = BUILD / "lean" / f"{mod_name}.olean"
    env = os.environ.copy()
    env["LEAN_PATH"] = f"{BUILD / 'lean'}:{FLOW_DIR / 'build' / 'lean'}"

    subprocess.run(
        [
            str(LEAN),
            "-DwarningAsError=true",
            f"--root={BUILD}",
            "-o",
            str(olean),
            str(lean_path),
        ],
        check=True,
        env=env,
    )
    subprocess.run([str(LEANCHECKER), mod_name], check=True, env=env)

    src_text = lean_path.read_text(encoding="utf-8")
    assert "sorry" not in src_text.lower(), f"Found sorry in {lean_path}"


def batch_audit_lean_axioms(certs: List[Tuple[str, List[str]]]) -> None:
    env = os.environ.copy()
    env["LEAN_PATH"] = f"{BUILD / 'lean'}:{FLOW_DIR / 'build' / 'lean'}"
    imports = "\n".join(f"import {mod_name}" for mod_name, _ in certs)
    axioms = "\n".join(f"#print axioms {t}" for _, targets in certs for t in targets)
    axiom_script = f"{imports}\n{axioms}\n"
    proc = subprocess.run(
        [str(LEAN), "--stdin"],
        input=axiom_script,
        text=True,
        capture_output=True,
        check=True,
        env=env,
    )
    assert "sorryAx" not in proc.stdout, f"sorryAx leaked in batch axiom audit:\n{proc.stdout}"
    for _, targets in certs:
        for t in targets:
            assert t in proc.stdout, f"Missing #print axioms output for {t}"


def run_hunt(
    mode: str, seed: int, trials: int, lean_out: Path, fir_out: Optional[Path], json_out: Path
) -> Tuple[dict, str]:
    cmd = [
        str(PHY_BIN),
        "hunt",
        mode,
        "--seed",
        str(seed),
        "--trials",
        str(trials),
        "--out-lean",
        str(lean_out),
        "--out-json",
        str(json_out),
    ]
    if fir_out is not None:
        cmd.extend(["--out-fir", str(fir_out)])
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    cli_line = f"  [CLI {mode}] {res.stdout.strip()}"
    return json.loads(json_out.read_text(encoding="utf-8")), cli_line


ModeResult = Tuple[int, str, str, str, List[str], Optional[Tuple[str, Path]]]


def _mode_01_entangled() -> ModeResult:
    ent_lean = BUILD / "PhyDiscoveredCert.lean"
    ent_fir = BUILD / "discovered_noether.fir"
    ent_json = BUILD / "discovered_noether.json"
    ent, cli_line = run_hunt("entangled", 2026, 2000, ent_lean, ent_fir, ent_json)
    p_ent, q_ent = ent["p"], ent["q"]
    assert len(set(p_ent)) == 16 and p_ent[0] == 0
    assert anf_degree(p_ent) == ent["anf_deg"] >= 2
    assert not is_factorizable(p_ent)
    assert verify_nonseparable_fraction(q_ent), "Entangled charge lies in 6-dim separable span!"
    j_nz, _ = verify_continuity_exact(p_ent, q_ent)
    assert j_nz == ent["j_nonzero"] and j_nz > 0
    verify_multi_length_floquet(p_ent, q_ent, seed=2026)
    dim4 = probe_4site_conservation(p_ent, trials=240)
    verify_fir_eval(ent_fir, p_ent)
    compile_and_check_lean_cert(ent_lean, "PhyDiscoveredCert")
    oracle_line = (
        f"  [Oracle ENTANGLED] PASS: deg={ent['anf_deg']} ent={ent['rank_ent']} j_nz={j_nz}/256 "
        f"N=(3,5,11,29,64) reversible+exact, 4site_null_ub={dim4}, Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredNoether.global_conservation",
        "Phy.DiscoveredNoether.charge_nonseparable",
        "Phy.DiscoveredNoether.current_active",
        "Phy.DiscoveredNoether.ruleP_nonlinear",
        "Phy.DiscoveredNoether.ruleP_nonfactorizable",
    ]
    return (1, cli_line, oracle_line, "PhyDiscoveredCert", targets, None)


def _mode_02_yangbox() -> ModeResult:
    yb_lean = BUILD / "PhyDiscoveredYB.lean"
    yb_fir = BUILD / "discovered_yb.fir"
    yb_json = BUILD / "discovered_yb.json"
    yb, cli_line = run_hunt("yangbox", 2026, 6000, yb_lean, yb_fir, yb_json)
    p_yb = yb["p"]
    assert len(set(p_yb)) == 16 and p_yb[0] == 0
    assert anf_degree(p_yb) == yb["anf_deg"] >= 2
    assert not is_factorizable(p_yb)
    assert not is_lyubashenko(p_yb)
    assert yb_violations(p_yb) == 0
    verify_fir_eval(yb_fir, p_yb)
    compile_and_check_lean_cert(yb_lean, "PhyDiscoveredYB")
    oracle_line = (
        f"  [Oracle YANGBOX]   PASS: deg={yb['anf_deg']} yb_viol=0/64 non-Lyubashenko "
        f"FIR_AND={yb['fir_and_count']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredYB.ruleYB_bijective",
        "Phy.DiscoveredYB.ruleYB_satisfies_ybe",
        "Phy.DiscoveredYB.ruleYB_nonlinear",
        "Phy.DiscoveredYB.ruleYB_nonfactorizable",
        "Phy.DiscoveredYB.ruleYB_ne_swap",
    ]
    return (2, cli_line, oracle_line, "PhyDiscoveredYB", targets, None)


def _mode_03_joint() -> ModeResult:
    jnt_lean = BUILD / "PhyDiscoveredJoint.lean"
    jnt_fir = BUILD / "discovered_joint.fir"
    jnt_json = BUILD / "discovered_joint.json"
    jnt, cli_line = run_hunt("joint", 2026, 20000, jnt_lean, jnt_fir, jnt_json)
    p_jnt, q_jnt = jnt["p"], jnt["q"]
    assert len(set(p_jnt)) == 16 and p_jnt[0] == 0
    assert anf_degree(p_jnt) == jnt["anf_deg"] >= 2
    assert not is_factorizable(p_jnt)
    assert not is_lyubashenko(p_jnt)
    assert yb_violations(p_jnt) == 0
    assert verify_nonseparable_fraction(q_jnt), "Joint charge lies in 6-dim separable span!"
    jnt_nz, _ = verify_continuity_exact(p_jnt, q_jnt)
    verify_multi_length_floquet(p_jnt, q_jnt, seed=1337)
    verify_fir_eval(jnt_fir, p_jnt)
    compile_and_check_lean_cert(jnt_lean, "PhyDiscoveredJoint")
    oracle_line = (
        f"  [Oracle JOINT]     PASS: deg={jnt['anf_deg']} yb_viol=0/64 tot={jnt['rank_tot']} "
        f"ent={jnt['rank_ent']} j_nz={jnt_nz}/256 FIR_AND={jnt['fir_and_count']} "
        f"Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredJoint.global_conservation",
        "Phy.DiscoveredJoint.ruleP_satisfies_ybe",
        "Phy.DiscoveredJoint.charge_nonseparable",
        "Phy.DiscoveredJoint.ruleP_nonlinear",
        "Phy.DiscoveredJoint.ruleP_nonfactorizable",
    ]
    return (3, cli_line, oracle_line, "PhyDiscoveredJoint", targets, None)


def _mode_04_q4() -> ModeResult:
    q4_lean = BUILD / "PhyDiscoveredQ4.lean"
    q4_fir = BUILD / "discovered_q4.fir"
    q4_json = BUILD / "discovered_q4.json"
    q4_res, cli_line = run_hunt("q4", 2026, 12000, q4_lean, q4_fir, q4_json)
    p_q4, q4_tbl = q4_res["p"], q4_res["q4"]
    assert len(set(p_q4)) == 16 and p_q4[0] == 0
    assert anf_degree(p_q4) >= 2 and not is_factorizable(p_q4)
    assert verify_q4_nontrivial_fraction(q4_tbl), "Q4 density lies in 48-dim trivial span!"
    j4_nz = verify_q4_continuity_exact(p_q4, q4_tbl)
    assert j4_nz == q4_res["j4_nonzero"] and j4_nz > 0
    verify_multi_length_floquet_q4(p_q4, q4_tbl, seed=2026)
    verify_fir_eval(q4_fir, p_q4, "eval-fir")
    compile_and_check_lean_cert(q4_lean, "PhyDiscoveredQ4")
    oracle_line = (
        f"  [Oracle Q4]        PASS: tot2={q4_res['tot2']} tot4={q4_res['tot4']} "
        f"genuine_q4={q4_res['genuine_q4']} j4_nz={j4_nz}/4096 bf_span={q4_res['bf_max_span']} "
        f"Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredQ4.global_conservation4",
        "Phy.DiscoveredQ4.q4_nontrivial",
        "Phy.DiscoveredQ4.j4_active",
        "Phy.DiscoveredQ4.ruleP_bijective",
    ]
    return (4, cli_line, oracle_line, "PhyDiscoveredQ4", targets, ("q4", q4_fir))


def _mode_05_yb8() -> ModeResult:
    yb8_lean = BUILD / "PhyDiscoveredYB8.lean"
    yb8_fir = BUILD / "discovered_yb8.fir"
    yb8_json = BUILD / "discovered_yb8.json"
    yb8_res, cli_line = run_hunt("yb8", 2026, 4000, yb8_lean, yb8_fir, yb8_json)
    r_yb8 = yb8_res["r"]
    deg8, order8 = verify_yb8_exact(r_yb8)
    assert deg8 == yb8_res["anf_deg"] >= 2
    assert order8 == yb8_res["order"] >= 4
    verify_fir_eval(yb8_fir, r_yb8, "eval-fir6")
    compile_and_check_lean_cert(yb8_lean, "PhyDiscoveredYB8")
    oracle_line = (
        f"  [Oracle YB8]       PASS: order={order8} deg={deg8} yb8_viol=0/512 nondeg=1 "
        f"noninvol=1 indecomp=1 FIR6_AND={yb8_res['fir_and_count']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredYB8.ruleR8_bijective",
        "Phy.DiscoveredYB8.ruleR8_satisfies_ybe",
        "Phy.DiscoveredYB8.ruleR8_nondegenerate",
        "Phy.DiscoveredYB8.ruleR8_noninvolutive",
        "Phy.DiscoveredYB8.ruleR8_indecomposable",
        "Phy.DiscoveredYB8.ruleR8_nonlyubashenko",
        "Phy.DiscoveredYB8.ruleR8_nonlinear",
    ]
    return (5, cli_line, oracle_line, "PhyDiscoveredYB8", targets, None)


def _mode_06_dual4() -> ModeResult:
    du4_lean = BUILD / "PhyDiscoveredDual4.lean"
    du4_fir = BUILD / "discovered_dual4.fir"
    du4_json = BUILD / "discovered_dual4.json"
    du4_res, cli_line = run_hunt("dual4", 2026, 4000, du4_lean, du4_fir, du4_json)
    p_du4 = du4_res["p"]
    verify_dual4_exact(p_du4, du4_res)
    verify_fir_eval(du4_fir, p_du4, "eval-fir")
    compile_and_check_lean_cert(du4_lean, "PhyDiscoveredDual4")
    oracle_line = (
        f"  [Oracle DUAL4]     PASS: deg={du4_res['anf_deg']} mix=1 rank_ent={du4_res['rank_ent']} "
        f"FIR_AND={du4_res['fir_and_count']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredDual4.ruleDU4_bijective",
        "Phy.DiscoveredDual4.ruleDU4_nondegenerate",
        "Phy.DiscoveredDual4.ruleDU4_dual_unitary",
        "Phy.DiscoveredDual4.ruleDU4_q2_nogo",
        "Phy.DiscoveredDual4.ruleDU4_mixing",
        "Phy.DiscoveredDual4.ruleDU4_nonlinear",
    ]
    return (6, cli_line, oracle_line, "PhyDiscoveredDual4", targets, ("dual4", du4_fir))


def _mode_07_ame8() -> ModeResult:
    ame_lean = BUILD / "PhyDiscoveredAME8.lean"
    ame_fir = BUILD / "discovered_ame8.fir"
    ame_json = BUILD / "discovered_ame8.json"
    ame_res, cli_line = run_hunt("ame8", 2026, 4000, ame_lean, ame_fir, ame_json)
    r_ame = ame_res["r"]
    deg_ame = verify_ame8_exact(r_ame, ame_res)
    verify_fir_eval(ame_fir, r_ame, "eval-fir6")
    compile_and_check_lean_cert(ame_lean, "PhyDiscoveredAME8")
    oracle_line = (
        f"  [Oracle AME8]      PASS: order={ame_res['order']} deg={deg_ame} 2unitary=1 "
        f"FIR6_AND={ame_res['fir_and_count']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredAME8.ruleAME8_bijective",
        "Phy.DiscoveredAME8.ruleAME8_nondegenerate",
        "Phy.DiscoveredAME8.ruleAME8_crossbijective",
        "Phy.DiscoveredAME8.ruleAME8_is_2unitary",
        "Phy.DiscoveredAME8.ruleAME8_nonlinear",
    ]
    return (7, cli_line, oracle_line, "PhyDiscoveredAME8", targets, ("ame8", ame_fir))


def _mode_08_fracton() -> ModeResult:
    fr_lean = BUILD / "PhyDiscoveredFracton.lean"
    fr_fir = BUILD / "discovered_fracton.fir"
    fr_json = BUILD / "discovered_fracton.json"
    fr_res, cli_line = run_hunt("fracton", 2026, 4000, fr_lean, fr_fir, fr_json)
    t_fr = fr_res["t"]
    deg_fr = verify_fracton_exact(t_fr, fr_res)
    verify_fir_eval(fr_fir, t_fr, "eval-fir6")
    compile_and_check_lean_cert(fr_lean, "PhyDiscoveredFracton")
    oracle_line = (
        f"  [Oracle FRACTON]   PASS: deg={deg_fr} comps={fr_res['hsf_num_comps']}/{fr_res['hsf_sector_size']} "
        f"largest={fr_res['hsf_largest_comp']} iso=0 dip={fr_res['dip_span']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredFracton.ruleFrac_bijective",
        "Phy.DiscoveredFracton.ruleFrac_dipole_laplacian",
        "Phy.DiscoveredFracton.ruleFrac_global_Q_conservation",
        "Phy.DiscoveredFracton.ruleFrac_global_D_conservation",
        "Phy.DiscoveredFracton.ruleFrac_single_charge_frozen",
        "Phy.DiscoveredFracton.ruleFrac_quadrupole_active",
        "Phy.DiscoveredFracton.ruleFrac_nonlinear",
    ]
    return (8, cli_line, oracle_line, "PhyDiscoveredFracton", targets, ("fracton", fr_fir))


def _mode_09_biquandle8() -> ModeResult:
    bq_lean = BUILD / "PhyDiscoveredBiquandle8.lean"
    bq_fir = BUILD / "discovered_biquandle8.fir"
    bq_json = BUILD / "discovered_biquandle8.json"
    bq_res, cli_line = run_hunt("biquandle8", 2026, 4000, bq_lean, bq_fir, bq_json)
    r_bq = bq_res["r"]
    deg_bq = verify_biquandle8_exact(r_bq, bq_res)
    verify_fir_eval(bq_fir, r_bq, "eval-fir6")
    compile_and_check_lean_cert(bq_lean, "PhyDiscoveredBiquandle8")
    oracle_line = (
        f"  [Oracle BIQ8]      PASS: order={bq_res['order']} deg={deg_bq} ret=8 R1=1 "
        f"T7={bq_res['col_torus7']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredBiquandle8.ruleBQ8_bijective",
        "Phy.DiscoveredBiquandle8.ruleBQ8_satisfies_ybe",
        "Phy.DiscoveredBiquandle8.ruleBQ8_nondegenerate",
        "Phy.DiscoveredBiquandle8.ruleBQ8_r1_biquandle",
        "Phy.DiscoveredBiquandle8.ruleBQ8_noninvolutive",
        "Phy.DiscoveredBiquandle8.ruleBQ8_indecomposable",
        "Phy.DiscoveredBiquandle8.ruleBQ8_irretractable",
        "Phy.DiscoveredBiquandle8.ruleBQ8_nonlyubashenko",
        "Phy.DiscoveredBiquandle8.ruleBQ8_nonlinear",
        "Phy.DiscoveredBiquandle8.ruleBQ8_knot_separation",
    ]
    return (9, cli_line, oracle_line, "PhyDiscoveredBiquandle8", targets, ("biquandle8", bq_fir))


def _mode_10_tetra4() -> ModeResult:
    cls_out = subprocess.check_output([str(PHY_BIN), "classify-tetra2"], text=True).strip()
    assert "total=26" in cls_out and "nonlinear=2" in cls_out and "noninvol_nl=0" in cls_out, cls_out
    t4_lean = BUILD / "PhyDiscoveredTetra4.lean"
    t4_fir = BUILD / "discovered_tetra4.fir"
    t4_json = BUILD / "discovered_tetra4.json"
    t4_res, cli_line = run_hunt("tetra4", 2026, 4000, t4_lean, t4_fir, t4_json)
    t_t4 = t4_res["t"]
    deg_t4 = verify_tetra4_exact(t_t4, t4_res)
    verify_fir_eval(t4_fir, t_t4, "eval-fir6")
    compile_and_check_lean_cert(t4_lean, "PhyDiscoveredTetra4")
    oracle_line = (
        f"  [Oracle TETRA4]    PASS: S8=(26,2,0) order={t4_res['order']} deg={deg_t4} viol=0/4096 "
        f"Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredTetra4.ruleT4_bijective",
        "Phy.DiscoveredTetra4.ruleT4_satisfies_tetra4",
        "Phy.DiscoveredTetra4.ruleT4_noninvolutive",
        "Phy.DiscoveredTetra4.ruleT4_3coupled",
        "Phy.DiscoveredTetra4.ruleT4_nonlinear",
    ]
    return (10, cli_line, oracle_line, "PhyDiscoveredTetra4", targets, ("tetra4", t4_fir))


def _mode_11_ame4_nogo() -> ModeResult:
    ame4_lean = BUILD / "PhyDiscoveredAME4NoGo.lean"
    ame4_json = BUILD / "discovered_ame4_nogo.json"
    ame4_res, cli_line = run_hunt("ame4-nogo", 2026, 1, ame4_lean, None, ame4_json)
    verify_ame4_nogo_exact(ame4_res)
    compile_and_check_lean_cert(ame4_lean, "PhyDiscoveredAME4NoGo")
    oracle_line = (
        f"  [Oracle AME4-NOGO] PASS: ls={ame4_res['latin4_count']} ols={ame4_res['total_ols']} "
        f"aff={ame4_res['affine_ols']} nl={ame4_res['nonlin_ols']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredAME4NoGo.universal_ame4_affine_rigidity",
        "Phy.DiscoveredAME4NoGo.sampleAME4_is_2unitary",
        "Phy.DiscoveredAME4NoGo.sampleAME4_forced_affine",
    ]
    return (11, cli_line, oracle_line, "PhyDiscoveredAME4NoGo", targets, None)


def _mode_12_simplex4() -> ModeResult:
    s4_lean = BUILD / "PhyDiscoveredSimplex4.lean"
    s4_fir = BUILD / "discovered_simplex4_ccx.fir"
    s4_json = BUILD / "discovered_simplex4.json"
    s4_res, cli_line = run_hunt("simplex4", 2026, 16, s4_lean, s4_fir, s4_json)
    deg_s4 = verify_simplex4_exact(s4_res)
    verify_fir_eval(s4_fir, s4_res["ccx02"], "eval-fir")
    compile_and_check_lean_cert(s4_lean, "PhyDiscoveredSimplex4")
    oracle_line = (
        f"  [Oracle SIMPLEX4]  PASS: viol4=0/1048576 order={s4_res['order4']} deg={deg_s4} "
        f"ccx=0 cccx={s4_res['viol2_cccx']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredSimplex4.ccx02_4simplex",
        "Phy.DiscoveredSimplex4.ccx13_4simplex",
        "Phy.DiscoveredSimplex4.cccx1_not_4simplex",
        "Phy.DiscoveredSimplex4.ruleS4_bijective",
        "Phy.DiscoveredSimplex4.ruleS4_satisfies_4simplex",
        "Phy.DiscoveredSimplex4.ruleS4_noninvolutive",
        "Phy.DiscoveredSimplex4.coupled4ch_satisfies_4simplex",
        "Phy.DiscoveredSimplex4.coupled4ch_4coupled",
    ]
    return (12, cli_line, oracle_line, "PhyDiscoveredSimplex4", targets, ("simplex4", s4_fir))


def _mode_13_chiral6() -> ModeResult:
    ch_lean = BUILD / "PhyDiscoveredChiral6.lean"
    ch_json = BUILD / "discovered_chiral6.json"
    ch_res, cli_line = run_hunt("chiral6", 2026, 1, ch_lean, None, ch_json)
    verify_chiral6_exact(ch_res)
    compile_and_check_lean_cert(ch_lean, "PhyDiscoveredChiral6")
    oracle_line = (
        f"  [Oracle CHIRAL6]   PASS: rank={ch_res['cocycle_rank']} L3_1={ch_res['tref_l_w']} "
        f"R3_1={ch_res['tref_r_w']} Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredChiral6.q6_quandle_ybe",
        "Phy.DiscoveredChiral6.q6_quandle_r1",
        "Phy.DiscoveredChiral6.discPhi_is_2cocycle",
        "Phy.DiscoveredChiral6.discPhi_not_coboundary",
        "Phy.DiscoveredChiral6.discPhi_chiral_separation",
    ]
    return (13, cli_line, oracle_line, "PhyDiscoveredChiral6", targets, None)


def _mode_14_qscar() -> ModeResult:
    qs_lean = BUILD / "PhyDiscoveredQScar.lean"
    qs_json = BUILD / "discovered_qscar.json"
    qs_res, cli_line = run_hunt("qscar", 2026, 1, qs_lean, None, qs_json)
    verify_qscar_exact(qs_res)
    compile_and_check_lean_cert(qs_lean, "PhyDiscoveredQScar")
    oracle_line = (
        f"  [Oracle QSCAR]     PASS: du_t2={qs_res['du_purity_t2_num']}/{qs_res['du_purity_t2_den']} "
        f"dtc={qs_res['scar_period']} therm_t4={qs_res['therm_purity_t4_num']}/{qs_res['therm_purity_t4_den']} "
        f"Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredQScar.discW_unitary",
        "Phy.DiscoveredQScar.discW_dual_unitary",
        "Phy.DiscoveredQScar.discW_full_superposition",
        "Phy.DiscoveredQScar.qscar_local_gate_unitary",
        "Phy.DiscoveredQScar.qscar_dtc_neel_protected",
        "Phy.DiscoveredQScar.qscar_dtc_exact_period2",
    ]
    return (14, cli_line, oracle_line, "PhyDiscoveredQScar", targets, None)


def _mode_15_fracton_comp() -> ModeResult:
    fc_lean = BUILD / "PhyDiscoveredFracComp.lean"
    fc_fir = BUILD / "discovered_fracton_comp.fir"
    fc_json = BUILD / "discovered_fracton_comp.json"
    fc_res, cli_line = run_hunt("fracton-comp", 2026, 4000, fc_lean, fc_fir, fc_json)
    verify_fracton_comp_exact(fc_res)
    verify_fir_eval(fc_fir, fc_res["t_toff"], "eval-fir6")
    compile_and_check_lean_cert(fc_lean, "PhyDiscoveredFracComp")
    oracle_line = (
        f"  [Oracle FRACCOMP]  PASS: comps={fc_res['comps']} largest={fc_res['largest']} "
        f"ccx=1 Lean+leanchecker PASS (0 sorryAx)"
    )
    targets = [
        "Phy.DiscoveredFracComp.ruleFracToff_bijective",
        "Phy.DiscoveredFracComp.ruleFracToff_laplacian",
        "Phy.DiscoveredFracComp.ruleFracToff_universal_involution",
        "Phy.DiscoveredFracComp.ruleFracToff_static_control_nogo",
        "Phy.DiscoveredFracComp.ruleFracToff_implements_ccx",
        "Phy.DiscoveredFracComp.ruleFracToff_single_charge_frozen",
    ]
    return (15, cli_line, oracle_line, "PhyDiscoveredFracComp", targets, ("fracton_comp", fc_fir))


_MODE_WORKERS = (
    _mode_12_simplex4,
    _mode_01_entangled,
    _mode_02_yangbox,
    _mode_03_joint,
    _mode_04_q4,
    _mode_05_yb8,
    _mode_06_dual4,
    _mode_07_ame8,
    _mode_08_fracton,
    _mode_09_biquandle8,
    _mode_10_tetra4,
    _mode_11_ame4_nogo,
    _mode_13_chiral6,
    _mode_14_qscar,
    _mode_15_fracton_comp,
)


def main() -> int:
    if not PHY_BIN.exists():
        print(f"ERROR: {PHY_BIN} not found", file=sys.stderr)
        return 1

    (BUILD / "lean").mkdir(parents=True, exist_ok=True)
    print("[Gate 5] Running PHY discovery hunts + independent L2 exact-rational verification...")

    workers = min(8, os.cpu_count() or 4)
    ctx = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as proc_pool:
        futures = [proc_pool.submit(fn) for fn in _MODE_WORKERS]
        results = [fut.result() for fut in futures]

    results.sort(key=lambda r: r[0])
    certs: List[Tuple[str, List[str]]] = []
    bridge_items: List[Tuple[str, Path]] = []
    for _, cli_line, oracle_line, mod_name, targets, bridge_item in results:
        print(cli_line)
        print(oracle_line)
        certs.append((mod_name, targets))
        if bridge_item is not None:
            bridge_items.append(bridge_item)

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as io_pool:
        ax_fut = io_pool.submit(batch_audit_lean_axioms, certs)
        bridge_lines = verify_cross_ecosystem_bridge(bridge_items, pool=io_pool)
        ax_fut.result()

    for line in bridge_lines:
        print(line)

    print("GATE 5 ORACLE VERIFICATION: ALL 15 MODES PASSED (0 sorryAx)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
