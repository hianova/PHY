#!/usr/bin/env python3
"""
PHY v3.0.0 — Independent Exact-Rational Level-2 Verification Oracle (Gate 5).
Independently re-verifies all 10 discovered physics/algebra modes (entangled, yangbox,
joint, q4, yb8, dual4, ame8, fracton, biquandle8, tetra4) using Python arbitrary-precision
ints, fractions.Fraction Gaussian elimination, multi-length Floquet / L=32 circuit
simulations, .fir bit-exact parity checks, live Lean 4 + leanchecker replay with
0 sorryAx, and QUA/ELE cross-ecosystem projection (when sibling builds are present).
"""

from __future__ import annotations

import json
import math
import os
import random
import shutil
import subprocess
import sys
from collections import deque
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple


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


def verify_cross_ecosystem_bridge(fir_paths: List[Tuple[str, Path]]) -> None:
    qcc = ROOT.parent / "QUA" / "qcc" / "build" / "qcc"
    ele = ROOT.parent / "ELE" / "build" / "ele_bootstrap"
    pdk = ROOT.parent / "ELE" / "build" / "pdk_sky130.flow"
    if qcc.exists():
        for tag, fp in fir_paths:
            kmb = BUILD / f"qcc_{tag}.kmb"
            subprocess.run([str(qcc), "opt", str(fp), "-o", str(kmb)], capture_output=True, text=True, check=True)
            vres = subprocess.run(
                [str(qcc), "verify", str(kmb), "--ref", str(fp), "--states", "4096"],
                capture_output=True,
                text=True,
                check=True,
            )
            assert "PASS" in vres.stdout or "val_err=0" in vres.stdout, f"qcc verify failed on {tag}: {vres.stdout}"
        print(f"  [Oracle QUA-qcc]   PASS ({len(fir_paths)}/{len(fir_paths)} .fir -> qcc .kmb 4096-state verify val_err=0)")
    if ele.exists() and pdk.exists():
        for tag, fp in fir_paths:
            sv = BUILD / f"ele_{tag}.sv"
            gds = BUILD / f"ele_{tag}.gds"
            eres = subprocess.run(
                [str(ele), "--fir-synth", str(fp), "-o", str(sv), "--pdk", str(pdk), "--gds-out", str(gds), "--verify-fir"],
                capture_output=True,
                text=True,
                check=True,
            )
            assert "DRC=0" in eres.stdout and "LVS=0" in eres.stdout, f"ELE synth failed on {tag}: {eres.stdout}"
        print(f"  [Oracle ELE-ASIC]  PASS ({len(fir_paths)}/{len(fir_paths)} .fir -> SKY130 SV+GDSII DRC=0 LVS=0)")


def verify_lean_cert(lean_path: Path, mod_name: str, axiom_targets: List[str]) -> None:
    olean = BUILD / "lean" / f"{mod_name}.olean"
    env = os.environ.copy()
    env["LEAN_PATH"] = str(BUILD / "lean")

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

    axiom_script = f"import {mod_name}\n" + "\n".join(f"#print axioms {t}" for t in axiom_targets) + "\n"
    proc = subprocess.run(
        [str(LEAN), "--stdin"],
        input=axiom_script,
        text=True,
        capture_output=True,
        check=True,
        env=env,
    )
    assert "sorryAx" not in proc.stdout, f"sorryAx leaked in {mod_name}:\n{proc.stdout}"


def run_hunt(mode: str, seed: int, trials: int, lean_out: Path, fir_out: Path, json_out: Path) -> dict:
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
        "--out-fir",
        str(fir_out),
        "--out-json",
        str(json_out),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    print(f"  [CLI {mode}] {res.stdout.strip()}")
    return json.loads(json_out.read_text(encoding="utf-8"))


def main() -> int:
    if not PHY_BIN.exists():
        print(f"ERROR: {PHY_BIN} not found", file=sys.stderr)
        return 1

    print("[Gate 5] Running PHY discovery hunts + independent L2 exact-rational verification...")

    # 1. Entangled Noether Universe
    ent_lean = BUILD / "PhyDiscoveredCert.lean"
    ent_fir = BUILD / "discovered_noether.fir"
    ent_json = BUILD / "discovered_noether.json"
    ent = run_hunt("entangled", 2026, 2000, ent_lean, ent_fir, ent_json)
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
    verify_lean_cert(
        ent_lean,
        "PhyDiscoveredCert",
        [
            "Phy.DiscoveredNoether.global_conservation",
            "Phy.DiscoveredNoether.charge_nonseparable",
            "Phy.DiscoveredNoether.current_active",
            "Phy.DiscoveredNoether.ruleP_nonlinear",
            "Phy.DiscoveredNoether.ruleP_nonfactorizable",
        ],
    )
    print(
        f"  [Oracle ENTANGLED] PASS: deg={ent['anf_deg']} ent={ent['rank_ent']} j_nz={j_nz}/256 "
        f"N=(3,5,11,29,64) reversible+exact, 4site_null_ub={dim4}, Lean+leanchecker PASS (0 sorryAx)"
    )

    # 2. Nonlinear Non-Lyubashenko Yang-Baxter Universe
    yb_lean = BUILD / "PhyDiscoveredYB.lean"
    yb_fir = BUILD / "discovered_yb.fir"
    yb_json = BUILD / "discovered_yb.json"
    yb = run_hunt("yangbox", 2026, 6000, yb_lean, yb_fir, yb_json)
    p_yb = yb["p"]
    assert len(set(p_yb)) == 16 and p_yb[0] == 0
    assert anf_degree(p_yb) == yb["anf_deg"] >= 2
    assert not is_factorizable(p_yb)
    assert not is_lyubashenko(p_yb)
    assert yb_violations(p_yb) == 0
    verify_fir_eval(yb_fir, p_yb)
    verify_lean_cert(
        yb_lean,
        "PhyDiscoveredYB",
        [
            "Phy.DiscoveredYB.ruleYB_bijective",
            "Phy.DiscoveredYB.ruleYB_satisfies_ybe",
            "Phy.DiscoveredYB.ruleYB_nonlinear",
            "Phy.DiscoveredYB.ruleYB_nonfactorizable",
            "Phy.DiscoveredYB.ruleYB_ne_swap",
        ],
    )
    print(
        f"  [Oracle YANGBOX]   PASS: deg={yb['anf_deg']} yb_viol=0/64 non-Lyubashenko "
        f"FIR_AND={yb['fir_and_count']} Lean+leanchecker PASS (0 sorryAx)"
    )

    # 3. Joint Yang-Baxter + Entangled Noether Universe
    jnt_lean = BUILD / "PhyDiscoveredJoint.lean"
    jnt_fir = BUILD / "discovered_joint.fir"
    jnt_json = BUILD / "discovered_joint.json"
    jnt = run_hunt("joint", 2026, 20000, jnt_lean, jnt_fir, jnt_json)
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
    verify_lean_cert(
        jnt_lean,
        "PhyDiscoveredJoint",
        [
            "Phy.DiscoveredJoint.global_conservation",
            "Phy.DiscoveredJoint.ruleP_satisfies_ybe",
            "Phy.DiscoveredJoint.charge_nonseparable",
            "Phy.DiscoveredJoint.ruleP_nonlinear",
            "Phy.DiscoveredJoint.ruleP_nonfactorizable",
        ],
    )
    print(
        f"  [Oracle JOINT]     PASS: deg={jnt['anf_deg']} yb_viol=0/64 tot={jnt['rank_tot']} "
        f"ent={jnt['rank_ent']} j_nz={jnt_nz}/256 FIR_AND={jnt['fir_and_count']} "
        f"Lean+leanchecker PASS (0 sorryAx)"
    )

    # 4. Higher-Order 4-Site (2-Block) Noether Hierarchy (Q2 < Q4)
    q4_lean = BUILD / "PhyDiscoveredQ4.lean"
    q4_fir = BUILD / "discovered_q4.fir"
    q4_json = BUILD / "discovered_q4.json"
    q4_res = run_hunt("q4", 2026, 12000, q4_lean, q4_fir, q4_json)
    p_q4, q4_tbl = q4_res["p"], q4_res["q4"]
    assert len(set(p_q4)) == 16 and p_q4[0] == 0
    assert anf_degree(p_q4) >= 2 and not is_factorizable(p_q4)
    assert verify_q4_nontrivial_fraction(q4_tbl), "Q4 density lies in 48-dim trivial span!"
    j4_nz = verify_q4_continuity_exact(p_q4, q4_tbl)
    assert j4_nz == q4_res["j4_nonzero"] and j4_nz > 0
    verify_multi_length_floquet_q4(p_q4, q4_tbl, seed=2026)
    verify_fir_eval(q4_fir, p_q4, "eval-fir")
    verify_lean_cert(
        q4_lean,
        "PhyDiscoveredQ4",
        [
            "Phy.DiscoveredQ4.global_conservation4",
            "Phy.DiscoveredQ4.q4_nontrivial",
            "Phy.DiscoveredQ4.j4_active",
            "Phy.DiscoveredQ4.ruleP_bijective",
        ],
    )
    print(
        f"  [Oracle Q4]        PASS: tot2={q4_res['tot2']} tot4={q4_res['tot4']} "
        f"genuine_q4={q4_res['genuine_q4']} j4_nz={j4_nz}/4096 bf_span={q4_res['bf_max_span']} "
        f"Lean+leanchecker PASS (0 sorryAx)"
    )

    # 5. |X|=8 Non-Involutive Non-Degenerate Indecomposable Set-Theoretic YBE
    yb8_lean = BUILD / "PhyDiscoveredYB8.lean"
    yb8_fir = BUILD / "discovered_yb8.fir"
    yb8_json = BUILD / "discovered_yb8.json"
    yb8_res = run_hunt("yb8", 2026, 4000, yb8_lean, yb8_fir, yb8_json)
    r_yb8 = yb8_res["r"]
    deg8, order8 = verify_yb8_exact(r_yb8)
    assert deg8 == yb8_res["anf_deg"] >= 2
    assert order8 == yb8_res["order"] >= 4
    verify_fir_eval(yb8_fir, r_yb8, "eval-fir6")
    verify_lean_cert(
        yb8_lean,
        "PhyDiscoveredYB8",
        [
            "Phy.DiscoveredYB8.ruleR8_bijective",
            "Phy.DiscoveredYB8.ruleR8_satisfies_ybe",
            "Phy.DiscoveredYB8.ruleR8_nondegenerate",
            "Phy.DiscoveredYB8.ruleR8_noninvolutive",
            "Phy.DiscoveredYB8.ruleR8_indecomposable",
            "Phy.DiscoveredYB8.ruleR8_nonlyubashenko",
            "Phy.DiscoveredYB8.ruleR8_nonlinear",
        ],
    )
    print(
        f"  [Oracle YB8]       PASS: order={order8} deg={deg8} yb8_viol=0/512 nondeg=1 "
        f"noninvol=1 indecomp=1 FIR6_AND={yb8_res['fir_and_count']} Lean+leanchecker PASS (0 sorryAx)"
    )

    # 6. Dual-Unitary Mixing & Universal Q2 Separability No-Go
    du4_lean = BUILD / "PhyDiscoveredDual4.lean"
    du4_fir = BUILD / "discovered_dual4.fir"
    du4_json = BUILD / "discovered_dual4.json"
    du4_res = run_hunt("dual4", 2026, 4000, du4_lean, du4_fir, du4_json)
    p_du4 = du4_res["p"]
    verify_dual4_exact(p_du4, du4_res)
    verify_fir_eval(du4_fir, p_du4, "eval-fir")
    verify_lean_cert(
        du4_lean,
        "PhyDiscoveredDual4",
        [
            "Phy.DiscoveredDual4.ruleDU4_bijective",
            "Phy.DiscoveredDual4.ruleDU4_nondegenerate",
            "Phy.DiscoveredDual4.ruleDU4_dual_unitary",
            "Phy.DiscoveredDual4.ruleDU4_q2_nogo",
            "Phy.DiscoveredDual4.ruleDU4_mixing",
            "Phy.DiscoveredDual4.ruleDU4_nonlinear",
        ],
    )
    print(
        f"  [Oracle DUAL4]     PASS: deg={du4_res['anf_deg']} mix=1 rank_ent={du4_res['rank_ent']} "
        f"FIR_AND={du4_res['fir_and_count']} Lean+leanchecker PASS (0 sorryAx)"
    )

    # 7. 2-Unitary Perfect Tensor AME(4,8) on BitVec 6
    ame_lean = BUILD / "PhyDiscoveredAME8.lean"
    ame_fir = BUILD / "discovered_ame8.fir"
    ame_json = BUILD / "discovered_ame8.json"
    ame_res = run_hunt("ame8", 2026, 4000, ame_lean, ame_fir, ame_json)
    r_ame = ame_res["r"]
    deg_ame = verify_ame8_exact(r_ame, ame_res)
    verify_fir_eval(ame_fir, r_ame, "eval-fir6")
    verify_lean_cert(
        ame_lean,
        "PhyDiscoveredAME8",
        [
            "Phy.DiscoveredAME8.ruleAME8_bijective",
            "Phy.DiscoveredAME8.ruleAME8_nondegenerate",
            "Phy.DiscoveredAME8.ruleAME8_crossbijective",
            "Phy.DiscoveredAME8.ruleAME8_is_2unitary",
            "Phy.DiscoveredAME8.ruleAME8_nonlinear",
        ],
    )
    print(
        f"  [Oracle AME8]      PASS: order={ame_res['order']} deg={deg_ame} 2unitary=1 "
        f"FIR6_AND={ame_res['fir_and_count']} Lean+leanchecker PASS (0 sorryAx)"
    )

    # 8. 3-Site Fracton Dipole-Laplacian & Strong HSF
    fr_lean = BUILD / "PhyDiscoveredFracton.lean"
    fr_fir = BUILD / "discovered_fracton.fir"
    fr_json = BUILD / "discovered_fracton.json"
    fr_res = run_hunt("fracton", 2026, 4000, fr_lean, fr_fir, fr_json)
    t_fr = fr_res["t"]
    deg_fr = verify_fracton_exact(t_fr, fr_res)
    verify_fir_eval(fr_fir, t_fr, "eval-fir6")
    verify_lean_cert(
        fr_lean,
        "PhyDiscoveredFracton",
        [
            "Phy.DiscoveredFracton.ruleFrac_bijective",
            "Phy.DiscoveredFracton.ruleFrac_dipole_laplacian",
            "Phy.DiscoveredFracton.ruleFrac_global_Q_conservation",
            "Phy.DiscoveredFracton.ruleFrac_global_D_conservation",
            "Phy.DiscoveredFracton.ruleFrac_single_charge_frozen",
            "Phy.DiscoveredFracton.ruleFrac_quadrupole_active",
            "Phy.DiscoveredFracton.ruleFrac_nonlinear",
        ],
    )
    print(
        f"  [Oracle FRACTON]   PASS: deg={deg_fr} comps={fr_res['hsf_num_comps']}/{fr_res['hsf_sector_size']} "
        f"largest={fr_res['hsf_largest_comp']} iso=0 dip={fr_res['dip_span']} Lean+leanchecker PASS (0 sorryAx)"
    )

    # 9. Irretractable |X|=8 R1+R2+R3 Biquandle & Knot Separation
    bq_lean = BUILD / "PhyDiscoveredBiquandle8.lean"
    bq_fir = BUILD / "discovered_biquandle8.fir"
    bq_json = BUILD / "discovered_biquandle8.json"
    bq_res = run_hunt("biquandle8", 2026, 4000, bq_lean, bq_fir, bq_json)
    r_bq = bq_res["r"]
    deg_bq = verify_biquandle8_exact(r_bq, bq_res)
    verify_fir_eval(bq_fir, r_bq, "eval-fir6")
    verify_lean_cert(
        bq_lean,
        "PhyDiscoveredBiquandle8",
        [
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
        ],
    )
    print(
        f"  [Oracle BIQ8]      PASS: order={bq_res['order']} deg={deg_bq} ret=8 R1=1 "
        f"T7={bq_res['col_torus7']} Lean+leanchecker PASS (0 sorryAx)"
    )

    # 10. 2+1D Zamolodchikov Tetrahedron Equation (|X|=2 S_8 + |X|=4)
    cls_out = subprocess.check_output([str(PHY_BIN), "classify-tetra2"], text=True).strip()
    assert "total=26" in cls_out and "nonlinear=2" in cls_out and "noninvol_nl=0" in cls_out, cls_out
    t4_lean = BUILD / "PhyDiscoveredTetra4.lean"
    t4_fir = BUILD / "discovered_tetra4.fir"
    t4_json = BUILD / "discovered_tetra4.json"
    t4_res = run_hunt("tetra4", 2026, 4000, t4_lean, t4_fir, t4_json)
    t_t4 = t4_res["t"]
    deg_t4 = verify_tetra4_exact(t_t4, t4_res)
    verify_fir_eval(t4_fir, t_t4, "eval-fir6")
    verify_lean_cert(
        t4_lean,
        "PhyDiscoveredTetra4",
        [
            "Phy.DiscoveredTetra4.ruleT4_bijective",
            "Phy.DiscoveredTetra4.ruleT4_satisfies_tetra4",
            "Phy.DiscoveredTetra4.ruleT4_noninvolutive",
            "Phy.DiscoveredTetra4.ruleT4_3coupled",
            "Phy.DiscoveredTetra4.ruleT4_nonlinear",
        ],
    )
    print(
        f"  [Oracle TETRA4]    PASS: S8=(26,2,0) order={t4_res['order']} deg={deg_t4} viol=0/4096 "
        f"Lean+leanchecker PASS (0 sorryAx)"
    )

    verify_cross_ecosystem_bridge(
        [
            ("q4", q4_fir),
            ("dual4", du4_fir),
            ("ame8", ame_fir),
            ("fracton", fr_fir),
            ("biquandle8", bq_fir),
            ("tetra4", t4_fir),
        ]
    )

    print("GATE 5 ORACLE VERIFICATION: ALL 10 MODES PASSED (0 sorryAx)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
