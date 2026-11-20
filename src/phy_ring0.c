#include "phy.h"

/* --- Ring-0 Non-Elidable Volatile Zero/Copy & Pure-GPR Helpers --- */

static void phy_vzero_u64(uint64_t *dst, size_t n) {
    volatile uint64_t *vp = (volatile uint64_t *)dst;
    for (size_t i = 0U; i < n; ++i) {
        vp[i] = UINT64_C(0);
    }
}

static void phy_vzero_i64(int64_t *dst, size_t n) {
    volatile int64_t *vp = (volatile int64_t *)dst;
    for (size_t i = 0U; i < n; ++i) {
        vp[i] = INT64_C(0);
    }
}

static void phy_vzero_i32(int32_t *dst, size_t n) {
    volatile int32_t *vp = (volatile int32_t *)dst;
    for (size_t i = 0U; i < n; ++i) {
        vp[i] = INT32_C(0);
    }
}

static void phy_vzero_u8(uint8_t *dst, size_t n) {
    volatile uint8_t *vp = (volatile uint8_t *)dst;
    for (size_t i = 0U; i < n; ++i) {
        vp[i] = UINT8_C(0);
    }
}

static void phy_vcopy_u8(uint8_t *dst, const uint8_t *src, size_t n) {
    volatile uint8_t *vd = (volatile uint8_t *)dst;
    for (size_t i = 0U; i < n; ++i) {
        vd[i] = src[i];
    }
}

static void phy_vcopy_i32(int32_t *dst, const int32_t *src, size_t n) {
    volatile int32_t *vd = (volatile int32_t *)dst;
    for (size_t i = 0U; i < n; ++i) {
        vd[i] = src[i];
    }
}

static void phy_vcopy_u64(uint64_t *dst, const uint64_t *src, size_t n) {
    volatile uint64_t *vd = (volatile uint64_t *)dst;
    for (size_t i = 0U; i < n; ++i) {
        vd[i] = src[i];
    }
}

static uint64_t phy_abs_i64(int64_t x) {
    uint64_t ux = (uint64_t)x;
    return (x < 0) ? (UINT64_C(0) - ux) : ux;
}

static uint64_t phy_abs_i128(flow_i128 x) {
    flow_u128 ux = (flow_u128)x;
    return (x < 0) ? (uint64_t)((flow_u128)0U - ux) : (uint64_t)ux;
}

/* --- FNV-1a & Division-Free Mersenne-61 + Binary Stein GCD --- */

uint64_t phy_fnv1a64_u64(uint64_t h, uint64_t v) {
    for (uint32_t b = 0U; b < 8U; ++b) {
        h ^= (v >> (b * 8U)) & UINT64_C(0xFF);
        h *= UINT64_C(0x100000001B3);
    }
    return h;
}

uint64_t phy_m61_add(uint64_t a, uint64_t b) {
    uint64_t s = a + b;
    uint64_t mask = flow_mask_bool(s >= PHY_M61);
    return s - (PHY_M61 & mask);
}

uint64_t phy_m61_sub(uint64_t a, uint64_t b) {
    uint64_t diff = a - b;
    uint64_t mask = flow_mask_bool(a < b);
    return diff + (PHY_M61 & mask);
}

static uint64_t phy_m61_reduce128(flow_u128 x) {
    uint64_t lo = (uint64_t)x & PHY_M61;
    uint64_t mid = (uint64_t)(x >> 61U) & PHY_M61;
    uint64_t hi = (uint64_t)(x >> 122U);
    uint64_t s = lo + mid + hi;
    s = (s & PHY_M61) + (s >> 61U);
    uint64_t mask = flow_mask_bool(s >= PHY_M61);
    return s - (PHY_M61 & mask);
}

uint64_t phy_m61_mul(uint64_t a, uint64_t b) {
    return phy_m61_reduce128((flow_u128)a * (flow_u128)b);
}

uint64_t phy_m61_inv(uint64_t a) {
    if (a == UINT64_C(0)) {
        return UINT64_C(0);
    }
    uint64_t exp = PHY_M61 - UINT64_C(2);
    uint64_t res = UINT64_C(1);
    uint64_t base = a;
    while (exp != UINT64_C(0)) {
        if ((exp & UINT64_C(1)) != UINT64_C(0)) {
            res = phy_m61_mul(res, base);
        }
        base = phy_m61_mul(base, base);
        exp >>= 1U;
    }
    return res;
}

uint64_t phy_stein_gcd64(uint64_t a, uint64_t b) {
    if (a == UINT64_C(0)) {
        return b;
    }
    if (b == UINT64_C(0)) {
        return a;
    }
    uint32_t az = searu_ctz64(a);
    uint32_t bz = searu_ctz64(b);
    uint32_t shift = (az < bz) ? az : bz;
    a >>= az;
    b >>= bz;
    while (a != b) {
        if (a > b) {
            a -= b;
            a >>= searu_ctz64(a);
        } else {
            b -= a;
            b >>= searu_ctz64(b);
        }
    }
    return a << shift;
}

uint64_t phy_div_u64(uint64_t num, uint64_t den, uint64_t *rem) {
    if (den == UINT64_C(0)) {
        if (rem != NULL) {
            *rem = UINT64_C(0);
        }
        return UINT64_C(0);
    }
    if (den == UINT64_C(1)) {
        if (rem != NULL) {
            *rem = UINT64_C(0);
        }
        return num;
    }
    uint64_t q = UINT64_C(0);
    uint64_t r = UINT64_C(0);
    for (uint32_t idx = 64U; idx > 0U; --idx) {
        uint32_t b = idx - 1U;
        r = (r << 1U) | ((num >> b) & UINT64_C(1));
        if (r >= den) {
            r -= den;
            q |= (UINT64_C(1) << b);
        }
    }
    if (rem != NULL) {
        *rem = r;
    }
    return q;
}

static int64_t phy_div_s64(int64_t num, uint64_t den) {
    if (den == UINT64_C(0)) {
        return INT64_C(0);
    }
    uint64_t mag = phy_div_u64(phy_abs_i64(num), den, NULL);
    return (num < 0) ? -(int64_t)mag : (int64_t)mag;
}

static int64_t phy_div_i128_by_u64(flow_i128 num, uint64_t den) {
    if (den == UINT64_C(0)) {
        return INT64_C(0);
    }
    flow_u128 an = (num < 0) ? ((flow_u128)0U - (flow_u128)num) : (flow_u128)num;
    flow_u128 q = (flow_u128)0U;
    flow_u128 r = (flow_u128)0U;
    for (uint32_t idx = 128U; idx > 0U; --idx) {
        uint32_t b = idx - 1U;
        r = (r << 1U) | ((an >> b) & (flow_u128)1U);
        if (r >= (flow_u128)den) {
            r -= (flow_u128)den;
            q |= ((flow_u128)1U << b);
        }
    }
    uint64_t q64 = (uint64_t)q;
    return (num < 0) ? -(int64_t)q64 : (int64_t)q64;
}

/* --- 2-Bit Pair Rule Operations, ANF, Non-Factorizability, Yang-Baxter --- */

static uint8_t phy_hi2(uint8_t u) {
    return (uint8_t)(((uint32_t)u >> 2U) & UINT32_C(3));
}

static uint8_t phy_lo2(uint8_t u) {
    return (uint8_t)((uint32_t)u & UINT32_C(3));
}

static uint8_t phy_pack2(uint8_t l, uint8_t r) {
    return (uint8_t)((((uint32_t)l & UINT32_C(3)) << 2U) | ((uint32_t)r & UINT32_C(3)));
}

uint64_t phy_pack_lut(const uint8_t p[16]) {
    uint64_t acc = UINT64_C(0);
    for (uint32_t i = 0U; i < 16U; ++i) {
        acc |= ((uint64_t)(p[i] & 0xFU)) << (i * 4U);
    }
    return acc;
}

void phy_unpack_lut(uint64_t packed, uint8_t p[16]) {
    for (uint32_t i = 0U; i < 16U; ++i) {
        p[i] = (uint8_t)((packed >> (i * 4U)) & 0xFU);
    }
}

uint8_t phy_step3(const uint8_t p[16], uint8_t a, uint8_t b, uint8_t c) {
    uint8_t pa = p[a & 0xFU];
    uint8_t pb = p[b & 0xFU];
    uint8_t pc = p[c & 0xFU];
    uint8_t oddL = p[phy_pack2(phy_lo2(pa), phy_hi2(pb))];
    uint8_t oddR = p[phy_pack2(phy_lo2(pb), phy_hi2(pc))];
    return phy_pack2(phy_lo2(oddL), phy_hi2(oddR));
}

uint8_t phy_anf_degree4(const uint8_t p[16], uint8_t *out_bit, uint8_t *out_a, uint8_t *out_b) {
    uint8_t max_deg = 0U;
    uint8_t best_bit = 0U;
    uint8_t best_a = 0U;
    uint8_t best_b = 0U;
    bool found_nl = false;
    for (uint8_t bit = 0U; bit < 4U; ++bit) {
        uint8_t anf[16];
        for (uint32_t x = 0U; x < 16U; ++x) {
            anf[x] = (uint8_t)(((uint32_t)p[x] >> (uint32_t)bit) & UINT32_C(1));
        }
        for (uint32_t var = 0U; var < 4U; ++var) {
            for (uint32_t mask = 0U; mask < 16U; ++mask) {
                if ((mask & (UINT32_C(1) << var)) != 0U) {
                    anf[mask] ^= anf[mask ^ (UINT32_C(1) << var)];
                }
            }
        }
        for (uint32_t mask = 0U; mask < 16U; ++mask) {
            if (anf[mask] != 0U) {
                uint8_t d = (uint8_t)searu_popcnt64((uint64_t)mask);
                if (d > max_deg) {
                    max_deg = d;
                }
            }
        }
        if (!found_nl) {
            for (uint8_t a = 1U; a < 16U && !found_nl; ++a) {
                for (uint8_t b = 1U; b < 16U && !found_nl; ++b) {
                    uint32_t xab = (uint32_t)p[a ^ b] >> (uint32_t)bit;
                    uint32_t xa = (uint32_t)p[a] >> (uint32_t)bit;
                    uint32_t xb = (uint32_t)p[b] >> (uint32_t)bit;
                    uint32_t x0 = (uint32_t)p[0] >> (uint32_t)bit;
                    uint8_t d2 = (uint8_t)((xab ^ xa ^ xb ^ x0) & UINT32_C(1));
                    if (d2 != 0U) {
                        found_nl = true;
                        best_bit = bit;
                        best_a = a;
                        best_b = b;
                    }
                }
            }
        }
    }
    if (out_bit != NULL) {
        *out_bit = best_bit;
    }
    if (out_a != NULL) {
        *out_a = best_a;
    }
    if (out_b != NULL) {
        *out_b = best_b;
    }
    return max_deg;
}

uint8_t phy_yb_check(const uint8_t p[16]) {
    uint8_t viol = 0U;
    for (uint8_t a = 0U; a < 4U; ++a) {
        for (uint8_t b = 0U; b < 4U; ++b) {
            for (uint8_t c = 0U; c < 4U; ++c) {
                /* R12 * R23 * R12 */
                uint8_t ab = p[phy_pack2(a, b)];
                uint8_t a1 = phy_hi2(ab);
                uint8_t b1 = phy_lo2(ab);
                uint8_t bc = p[phy_pack2(b1, c)];
                uint8_t b2 = phy_hi2(bc);
                uint8_t c2 = phy_lo2(bc);
                uint8_t ab2 = p[phy_pack2(a1, b2)];
                uint16_t lhs = (uint16_t)(((uint32_t)ab2 << 2U) | (uint32_t)c2);

                /* R23 * R12 * R23 */
                uint8_t bc0 = p[phy_pack2(b, c)];
                uint8_t b0 = phy_hi2(bc0);
                uint8_t c0 = phy_lo2(bc0);
                uint8_t ab0 = p[phy_pack2(a, b0)];
                uint8_t a0 = phy_hi2(ab0);
                uint8_t b00 = phy_lo2(ab0);
                uint8_t bc1 = p[phy_pack2(b00, c0)];
                uint16_t rhs = (uint16_t)(((uint32_t)a0 << 4U) | (uint32_t)bc1);

                if (lhs != rhs) {
                    ++viol;
                }
            }
        }
    }
    return viol;
}

static bool phy_is_lyubashenko(const uint8_t p[16]) {
    for (uint8_t a = 0U; a < 4U; ++a) {
        for (uint8_t b = 0U; b < 4U; ++b) {
            if (phy_hi2(p[phy_pack2(a, b)]) != phy_hi2(p[phy_pack2(0U, b)])) {
                return false;
            }
            if (phy_lo2(p[phy_pack2(a, b)]) != phy_lo2(p[phy_pack2(a, 0U)])) {
                return false;
            }
        }
    }
    return true;
}

uint32_t phy_perm_order(const uint8_t p[16]) {
    uint8_t vis[16];
    phy_vzero_u8(vis, 16U);
    uint64_t lcm = UINT64_C(1);
    for (uint8_t i = 0U; i < 16U; ++i) {
        if (vis[i] == 0U) {
            uint64_t clen = UINT64_C(0);
            uint8_t cur = i;
            while (vis[cur] == 0U && clen < 17U) {
                vis[cur] = 1U;
                cur = p[cur] & 0xFU;
                ++clen;
            }
            if (clen > UINT64_C(0)) {
                uint64_t g = phy_stein_gcd64(lcm, clen);
                lcm = phy_div_u64(lcm, g, NULL) * clen;
            }
        }
    }
    return (uint32_t)lcm;
}

/* --- Streaming Online M61 Echelon & Exact Integer Hattori-Takesue Solver --- */

static uint64_t phy_m61_from_signed(int64_t v) {
    if (v >= 0) {
        return (uint64_t)v;
    }
    return PHY_M61 - (uint64_t)(-v);
}

static void phy_continuity_row_terms(const uint8_t f00[16], const uint8_t f0[16][16],
                                     const uint8_t p[16], uint8_t a, uint8_t b, uint8_t c,
                                     int64_t row[15]) {
    phy_vzero_i64(row, 15U);
    uint8_t u0 = phy_step3(p, a, b, c);
    uint8_t u1 = f0[b][c];
    uint8_t u2 = f00[b];
    uint8_t u3 = f00[a];
    uint8_t u4 = f0[a][b];
    uint8_t u5 = a;
    if (u0 > 0U) {
        row[u0 - 1U] += INT64_C(1);
    }
    if (u1 > 0U) {
        row[u1 - 1U] -= INT64_C(1);
    }
    if (u2 > 0U) {
        row[u2 - 1U] -= INT64_C(1);
    }
    if (u3 > 0U) {
        row[u3 - 1U] += INT64_C(1);
    }
    if (u4 > 0U) {
        row[u4 - 1U] += INT64_C(1);
    }
    if (u5 > 0U) {
        row[u5 - 1U] -= INT64_C(1);
    }
}

static uint32_t phy_stream_m61_pivots(const uint8_t p[16], const uint8_t f00[16],
                                      const uint8_t f0[16][16], uint16_t piv_triples[15],
                                      uint64_t rref[15 * 15], uint8_t piv_col_present[15]) {
    phy_vzero_u64(rref, 225U);
    phy_vzero_u8(piv_col_present, 15U);
    uint32_t rank = 0U;

    for (uint8_t a = 0U; a < 16U && rank < 15U; ++a) {
        for (uint8_t b = 0U; b < 16U && rank < 15U; ++b) {
            for (uint8_t c = 0U; c < 16U && rank < 15U; ++c) {
                int64_t srow[15];
                phy_continuity_row_terms(f00, f0, p, a, b, c, srow);
                uint64_t mrow[15];
                for (uint32_t k = 0U; k < 15U; ++k) {
                    mrow[k] = phy_m61_from_signed(srow[k]);
                }
                for (uint32_t col = 0U; col < 15U; ++col) {
                    uint64_t val = mrow[col];
                    if (val != UINT64_C(0) && piv_col_present[col] != 0U) {
                        const uint64_t *prow = &rref[(size_t)col * 15U];
                        for (uint32_t j = col; j < 15U; ++j) {
                            uint64_t sub = phy_m61_mul(val, prow[j]);
                            mrow[j] = phy_m61_sub(mrow[j], sub);
                        }
                    }
                }
                uint32_t lead = 15U;
                for (uint32_t col = 0U; col < 15U; ++col) {
                    if (mrow[col] != UINT64_C(0)) {
                        lead = col;
                        break;
                    }
                }
                if (lead < 15U) {
                    uint64_t inv = phy_m61_inv(mrow[lead]);
                    uint64_t *dst = &rref[(size_t)lead * 15U];
                    for (uint32_t j = lead; j < 15U; ++j) {
                        dst[j] = phy_m61_mul(mrow[j], inv);
                    }
                    piv_col_present[lead] = 1U;
                    uint32_t enc = (((uint32_t)a) << 8U) | (((uint32_t)b) << 4U) | (uint32_t)c;
                    piv_triples[rank] = (uint16_t)enc;
                    ++rank;
                    if (rank == 15U) {
                        return 15U;
                    }
                }
            }
        }
    }
    return rank;
}

static uint32_t phy_extract_m61_null_basis(uint64_t rref[15 * 15],
                                           const uint8_t piv_col_present[15],
                                           uint64_t null_basis[15 * 15]) {
    /* Back-substitute strictly to reduced row echelon form over M61 */
    for (uint32_t idx = 15U; idx > 0U; --idx) {
        uint32_t col = idx - 1U;
        if (piv_col_present[col] != 0U) {
            const uint64_t *prow = &rref[(size_t)col * 15U];
            for (uint32_t up = 0U; up < col; ++up) {
                if (piv_col_present[up] != 0U) {
                    uint64_t fac = rref[(size_t)up * 15U + col];
                    if (fac != UINT64_C(0)) {
                        for (uint32_t j = col; j < 15U; ++j) {
                            uint64_t sub = phy_m61_mul(fac, prow[j]);
                            rref[(size_t)up * 15U + j] =
                                phy_m61_sub(rref[(size_t)up * 15U + j], sub);
                        }
                    }
                }
            }
        }
    }

    phy_vzero_u64(null_basis, 225U);
    uint32_t n_null = 0U;
    for (uint32_t f = 0U; f < 15U; ++f) {
        if (piv_col_present[f] == 0U) {
            uint64_t *vec = &null_basis[(size_t)n_null * 15U];
            vec[f] = UINT64_C(1);
            for (uint32_t col = 0U; col < f; ++col) {
                if (piv_col_present[col] != 0U) {
                    uint64_t v = rref[(size_t)col * 15U + f];
                    vec[col] = (v == UINT64_C(0)) ? UINT64_C(0) : (PHY_M61 - v);
                }
            }
            ++n_null;
        }
    }
    return n_null;
}

static uint8_t phy_count_entangled_dim_m61(const uint64_t null_basis[15 * 15], uint32_t n_null) {
    uint64_t ech[15 * 15];
    uint8_t has[15];
    phy_vzero_u64(ech, 225U);
    phy_vzero_u8(has, 15U);

    /* Insert 6 exact chiral separable basis vectors first */
    for (uint32_t s = 0U; s < 6U; ++s) {
        uint64_t vec[15];
        phy_vzero_u64(vec, 15U);
        uint8_t target = (uint8_t)((s < 3U) ? (s + 1U) : (s - 2U));
        for (uint8_t u = 1U; u < 16U; ++u) {
            uint8_t val = (s < 3U) ? phy_hi2(u) : phy_lo2(u);
            if (val == target) {
                vec[u - 1U] = UINT64_C(1);
            }
        }
        for (uint32_t col = 0U; col < 15U; ++col) {
            uint64_t v = vec[col];
            if (v != UINT64_C(0)) {
                if (has[col] != 0U) {
                    const uint64_t *prow = &ech[(size_t)col * 15U];
                    for (uint32_t j = col; j < 15U; ++j) {
                        vec[j] = phy_m61_sub(vec[j], phy_m61_mul(v, prow[j]));
                    }
                } else {
                    uint64_t inv = phy_m61_inv(v);
                    uint64_t *dst = &ech[(size_t)col * 15U];
                    for (uint32_t j = col; j < 15U; ++j) {
                        dst[j] = phy_m61_mul(vec[j], inv);
                    }
                    has[col] = 1U;
                    break;
                }
            }
        }
    }

    uint8_t ent_dim = 0U;
    for (uint32_t k = 0U; k < n_null; ++k) {
        uint64_t vec[15];
        phy_vcopy_u64(vec, &null_basis[(size_t)k * 15U], 15U);
        for (uint32_t col = 0U; col < 15U; ++col) {
            uint64_t v = vec[col];
            if (v != UINT64_C(0)) {
                if (has[col] != 0U) {
                    const uint64_t *prow = &ech[(size_t)col * 15U];
                    for (uint32_t j = col; j < 15U; ++j) {
                        vec[j] = phy_m61_sub(vec[j], phy_m61_mul(v, prow[j]));
                    }
                } else {
                    uint64_t inv = phy_m61_inv(v);
                    uint64_t *dst = &ech[(size_t)col * 15U];
                    for (uint32_t j = col; j < 15U; ++j) {
                        dst[j] = phy_m61_mul(vec[j], inv);
                    }
                    has[col] = 1U;
                    ++ent_dim;
                    break;
                }
            }
        }
    }
    return ent_dim;
}

static bool phy_find_nonsep_witness(const int32_t q[16], uint8_t *l0, uint8_t *l1, uint8_t *r0,
                                    uint8_t *r1) {
    for (uint8_t a0 = 0U; a0 < 4U; ++a0) {
        for (uint8_t a1 = (uint8_t)(a0 + 1U); a1 < 4U; ++a1) {
            for (uint8_t b0 = 0U; b0 < 4U; ++b0) {
                for (uint8_t b1 = (uint8_t)(b0 + 1U); b1 < 4U; ++b1) {
                    int64_t d = (int64_t)q[phy_pack2(a1, b1)] - (int64_t)q[phy_pack2(a1, b0)] -
                                (int64_t)q[phy_pack2(a0, b1)] + (int64_t)q[phy_pack2(a0, b0)];
                    if (d != INT64_C(0)) {
                        *l0 = a0;
                        *l1 = a1;
                        *r0 = b0;
                        *r1 = b1;
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

static uint16_t phy_fill_canon_current(const uint8_t f00[16], const uint8_t f0[16][16],
                                       const int32_t q[16], int32_t out_j[16][16]) {
    uint16_t j_nz = 0U;
    for (uint8_t b = 0U; b < 16U; ++b) {
        int32_t hb = q[f00[b]] - q[0];
        for (uint8_t c = 0U; c < 16U; ++c) {
            int32_t val = hb + (q[f0[b][c]] - q[b]);
            out_j[b][c] = val;
            if (val != 0) {
                ++j_nz;
            }
        }
    }
    return j_nz;
}

static int32_t phy_verify_continuity_i64(const uint8_t p[16], const uint8_t f00[16],
                                         const uint8_t f0[16][16], const int32_t q[16]) {
    int32_t fails = 0;
    for (uint8_t a = 0U; a < 16U; ++a) {
        int64_t ha = (int64_t)q[f00[a]] - (int64_t)q[0];
        for (uint8_t b = 0U; b < 16U; ++b) {
            int64_t jab = ha + ((int64_t)q[f0[a][b]] - (int64_t)q[a]);
            int64_t hb = (int64_t)q[f00[b]] - (int64_t)q[0];
            for (uint8_t c = 0U; c < 16U; ++c) {
                int64_t jbc = hb + ((int64_t)q[f0[b][c]] - (int64_t)q[b]);
                int64_t lhs = (int64_t)q[phy_step3(p, a, b, c)] - (int64_t)q[b];
                int64_t rhs = jbc - jab;
                if (lhs != rhs) {
                    ++fails;
                }
            }
        }
    }
    return fails;
}

static void phy_normalize_int64_row(int64_t row[15]) {
    uint64_t g = UINT64_C(0);
    for (uint32_t j = 0U; j < 15U; ++j) {
        if (row[j] != INT64_C(0)) {
            g = phy_stein_gcd64(g, phy_abs_i64(row[j]));
        }
    }
    if (g > UINT64_C(1)) {
        for (uint32_t j = 0U; j < 15U; ++j) {
            row[j] = phy_div_s64(row[j], g);
        }
    }
}

static bool phy_select_best_int_charge(const uint8_t p[16], const uint8_t f00[16],
                                       const uint8_t f0[16][16], const uint16_t piv_triples[15],
                                       uint32_t rank, phy_rule_analysis_t *out) {
    if (rank >= 15U) {
        return false;
    }
    int64_t imat[15 * 15];
    phy_vzero_i64(imat, 225U);
    for (uint32_t r = 0U; r < rank; ++r) {
        uint32_t enc = (uint32_t)piv_triples[r];
        uint8_t a = (uint8_t)((enc >> 8U) & UINT32_C(0xF));
        uint8_t b = (uint8_t)((enc >> 4U) & UINT32_C(0xF));
        uint8_t c = (uint8_t)(enc & UINT32_C(0xF));
        phy_continuity_row_terms(f00, f0, p, a, b, c, &imat[(size_t)r * 15U]);
    }

    uint8_t is_free[15];
    uint8_t piv_c[15];
    for (uint32_t c = 0U; c < 15U; ++c) {
        is_free[c] = 1U;
        piv_c[c] = 0U;
    }

    uint32_t cur_r = 0U;
    for (uint32_t col = 0U; col < 15U && cur_r < rank; ++col) {
        uint32_t found = rank;
        for (uint32_t r = cur_r; r < rank; ++r) {
            if (imat[(size_t)r * 15U + col] != INT64_C(0)) {
                found = r;
                break;
            }
        }
        if (found == rank) {
            continue;
        }
        if (found != cur_r) {
            for (uint32_t j = 0U; j < 15U; ++j) {
                int64_t tmp = imat[(size_t)cur_r * 15U + j];
                imat[(size_t)cur_r * 15U + j] = imat[(size_t)found * 15U + j];
                imat[(size_t)found * 15U + j] = tmp;
            }
        }
        if (imat[(size_t)cur_r * 15U + col] < 0) {
            for (uint32_t j = 0U; j < 15U; ++j) {
                imat[(size_t)cur_r * 15U + j] = -imat[(size_t)cur_r * 15U + j];
            }
        }
        phy_normalize_int64_row(&imat[(size_t)cur_r * 15U]);
        is_free[col] = 0U;
        piv_c[cur_r] = (uint8_t)col;

        int64_t pval = imat[cur_r * 15U + col];
        for (uint32_t k = 0U; k < rank; ++k) {
            if (k != cur_r && imat[k * 15U + col] != INT64_C(0)) {
                int64_t kval = imat[k * 15U + col];
                uint64_t g = phy_stein_gcd64(phy_abs_i64(pval), phy_abs_i64(kval));
                int64_t mk = phy_div_s64(pval, g);
                int64_t mp = phy_div_s64(kval, g);
                flow_i128 tmp_row[15];
                uint64_t rg = UINT64_C(0);
                for (uint32_t j = 0U; j < 15U; ++j) {
                    tmp_row[j] = (flow_i128)mk * (flow_i128)imat[k * 15U + j] -
                                 (flow_i128)mp * (flow_i128)imat[cur_r * 15U + j];
                    if (tmp_row[j] != 0) {
                        rg = phy_stein_gcd64(rg, phy_abs_i128(tmp_row[j]));
                    }
                }
                if (rg == UINT64_C(0)) {
                    rg = UINT64_C(1);
                }
                for (uint32_t j = 0U; j < 15U; ++j) {
                    imat[k * 15U + j] = phy_div_i128_by_u64(tmp_row[j], rg);
                }
            }
        }
        ++cur_r;
    }

    if (cur_r >= 15U) {
        return false;
    }

    uint64_t diag_lcm = UINT64_C(1);
    for (uint32_t r = 0U; r < cur_r; ++r) {
        uint64_t dv = phy_abs_i64(imat[r * 15U + piv_c[r]]);
        if (dv > UINT64_C(1)) {
            uint64_t g = phy_stein_gcd64(diag_lcm, dv);
            diag_lcm = phy_div_u64(diag_lcm, g, NULL) * dv;
        }
    }

    int32_t best_q[16];
    int32_t act_q[16];
    phy_vzero_i32(best_q, 16U);
    phy_vzero_i32(act_q, 16U);
    bool has_best = false;
    bool has_act = false;
    uint16_t best_jnz = 0U;
    uint16_t best_act_jnz = 0U;
    uint64_t best_max_abs = UINT64_MAX;

    for (uint32_t f = 0U; f < 15U; ++f) {
        if (is_free[f] == 0U) {
            continue;
        }
        int64_t vec64[16];
        phy_vzero_i64(vec64, 16U);
        vec64[f + 1U] = (int64_t)diag_lcm;
        for (uint32_t r = 0U; r < cur_r; ++r) {
            int64_t dv = imat[(size_t)r * 15U + piv_c[r]];
            int64_t num = -(int64_t)diag_lcm * imat[(size_t)r * 15U + f];
            vec64[piv_c[r] + 1U] = phy_div_s64(num, phy_abs_i64(dv));
            if (dv < 0) {
                vec64[piv_c[r] + 1U] = -vec64[piv_c[r] + 1U];
            }
        }
        uint64_t vg = UINT64_C(0);
        for (uint32_t u = 1U; u < 16U; ++u) {
            if (vec64[u] != INT64_C(0)) {
                vg = phy_stein_gcd64(vg, phy_abs_i64(vec64[u]));
            }
        }
        int32_t cand_q[16];
        phy_vzero_i32(cand_q, 16U);
        bool in_range = true;
        uint64_t max_abs = UINT64_C(0);
        for (uint32_t u = 1U; u < 16U; ++u) {
            int64_t nv = (vg > UINT64_C(1)) ? phy_div_s64(vec64[u], vg) : vec64[u];
            uint64_t av = phy_abs_i64(nv);
            if (av > UINT64_C(100000000)) {
                in_range = false;
            }
            if (av > max_abs) {
                max_abs = av;
            }
            cand_q[u] = (int32_t)nv;
        }
        if (!in_range || max_abs == UINT64_C(0)) {
            continue;
        }

        int32_t tmp_j[16][16];
        uint16_t jnz = phy_fill_canon_current(f00, f0, cand_q, tmp_j);
        if (jnz > best_act_jnz) {
            phy_vcopy_i32(act_q, cand_q, 16U);
            best_act_jnz = jnz;
            has_act = true;
        }
        uint8_t sl0 = 0U;
        uint8_t sl1 = 0U;
        uint8_t sr0 = 0U;
        uint8_t sr1 = 0U;
        if (phy_find_nonsep_witness(cand_q, &sl0, &sl1, &sr0, &sr1)) {
            if (phy_verify_continuity_i64(p, f00, f0, cand_q) == 0) {
                bool better = !has_best;
                if (!better) {
                    if (jnz > best_jnz || (jnz == best_jnz && max_abs < best_max_abs)) {
                        better = true;
                    }
                }
                if (better) {
                    has_best = true;
                    best_jnz = jnz;
                    best_max_abs = max_abs;
                    phy_vcopy_i32(best_q, cand_q, 16U);
                }
            }
        }
    }

    /* If best non-separable basis vector happens to have zero current while another conserved
       vector has non-zero current, combine them and normalize by Stein GCD */
    if (has_best && best_jnz == 0U && has_act) {
        int64_t comb64[16];
        uint64_t cg = UINT64_C(0);
        for (uint32_t u = 0U; u < 16U; ++u) {
            comb64[u] = (int64_t)best_q[u] + (int64_t)act_q[u];
            if (comb64[u] != INT64_C(0)) {
                cg = phy_stein_gcd64(cg, phy_abs_i64(comb64[u]));
            }
        }
        int32_t comb_q[16];
        bool ok = (cg > UINT64_C(0));
        for (uint32_t u = 0U; u < 16U && ok; ++u) {
            int64_t nv = (cg > UINT64_C(1)) ? phy_div_s64(comb64[u], cg) : comb64[u];
            if (phy_abs_i64(nv) > UINT64_C(100000000)) {
                ok = false;
            }
            comb_q[u] = (int32_t)nv;
        }
        if (ok && phy_verify_continuity_i64(p, f00, f0, comb_q) == 0) {
            uint8_t sl0 = 0U;
            uint8_t sl1 = 0U;
            uint8_t sr0 = 0U;
            uint8_t sr1 = 0U;
            if (phy_find_nonsep_witness(comb_q, &sl0, &sl1, &sr0, &sr1)) {
                phy_vcopy_i32(best_q, comb_q, 16U);
                (void)phy_fill_canon_current(f00, f0, best_q, out->j);
            }
        }
    }

    if (!has_best) {
        return false;
    }

    phy_vcopy_i32(out->q, best_q, 16U);
    out->j_nonzero = phy_fill_canon_current(f00, f0, out->q, out->j);
    out->continuity_fail = phy_verify_continuity_i64(p, f00, f0, out->q);
    (void)phy_find_nonsep_witness(out->q, &out->sep_l0, &out->sep_l1, &out->sep_r0, &out->sep_r1);
    return out->continuity_fail == 0;
}

static bool phy_check_factorizable(const uint8_t p[16], uint8_t *out_axis, uint8_t *out_l,
                                   uint8_t *out_r0, uint8_t *out_r1) {
    for (uint8_t a = 0U; a < 4U; ++a) {
        uint8_t h0 = phy_hi2(p[phy_pack2(a, 0U)]);
        for (uint8_t b = 1U; b < 4U; ++b) {
            if (phy_hi2(p[phy_pack2(a, b)]) != h0) {
                *out_axis = 0U;
                *out_l = a;
                *out_r0 = 0U;
                *out_r1 = b;
                return false;
            }
        }
    }
    for (uint8_t b = 0U; b < 4U; ++b) {
        uint8_t l0 = phy_lo2(p[phy_pack2(0U, b)]);
        for (uint8_t a = 1U; a < 4U; ++a) {
            if (phy_lo2(p[phy_pack2(a, b)]) != l0) {
                *out_axis = 1U;
                *out_l = b;
                *out_r0 = 0U;
                *out_r1 = a;
                return false;
            }
        }
    }
    return true;
}

static void phy_zero_rule_analysis(phy_rule_analysis_t *out) {
    phy_vcopy_u8(out->p, out->p, 0U);
    phy_vzero_u8(out->p, 16U);
    phy_vzero_u8(out->pinv, 16U);
    phy_vzero_i32(out->q, 16U);
    for (uint32_t r = 0U; r < 16U; ++r) {
        phy_vzero_i32(out->j[r], 16U);
    }
    out->lut_packed = UINT64_C(0);
    out->pinv_packed = UINT64_C(0);
    out->merkle = UINT64_C(0);
    out->seed = UINT64_C(0);
    out->flags = UINT32_C(0);
    out->cycle_lcm = UINT32_C(1);
    out->j_nonzero = UINT16_C(0);
    out->anf_deg = UINT8_C(0);
    out->yb_violations = UINT8_C(64);
    out->rank_tot = UINT8_C(0);
    out->rank_ent = UINT8_C(0);
    out->is_bijective = UINT8_C(0);
    out->is_swap = UINT8_C(0);
    out->is_id = UINT8_C(0);
    out->is_factorizable = UINT8_C(0);
    out->nl_bit = UINT8_C(0);
    out->nl_a = UINT8_C(0);
    out->nl_b = UINT8_C(0);
    out->sep_l0 = UINT8_C(0);
    out->sep_l1 = UINT8_C(0);
    out->sep_r0 = UINT8_C(0);
    out->sep_r1 = UINT8_C(0);
    out->nf_axis = UINT8_C(0);
    out->nf_l = UINT8_C(0);
    out->nf_r0 = UINT8_C(0);
    out->nf_r1 = UINT8_C(0);
    out->yb_neq_swap_u = UINT8_C(0);
    out->continuity_fail = 0;
}

void phy_analyze_rule(const uint8_t p[16], uint64_t seed, phy_rule_analysis_t *out) {
    if (out == NULL) {
        return;
    }
    phy_zero_rule_analysis(out);
    phy_vcopy_u8(out->p, p, 16U);
    out->seed = seed;
    out->lut_packed = phy_pack_lut(out->p);

    uint32_t seen_mask = UINT32_C(0);
    bool swap_match = true;
    bool id_match = true;
    uint8_t neq_swap = 0U;
    for (uint8_t u = 0U; u < 16U; ++u) {
        uint8_t v = out->p[u] & 0xFU;
        out->p[u] = v;
        seen_mask |= (UINT32_C(1) << v);
        out->pinv[v] = u;
        uint8_t sw = phy_pack2(phy_lo2(u), phy_hi2(u));
        if (v != sw) {
            if (swap_match) {
                neq_swap = u;
            }
            swap_match = false;
        }
        if (v != u) {
            id_match = false;
        }
    }
    out->is_bijective = (seen_mask == UINT32_C(0xFFFF)) ? UINT8_C(1) : UINT8_C(0);
    out->is_swap = swap_match ? UINT8_C(1) : UINT8_C(0);
    out->is_id = id_match ? UINT8_C(1) : UINT8_C(0);
    out->yb_neq_swap_u = neq_swap;
    out->pinv_packed = (out->is_bijective != 0U) ? phy_pack_lut(out->pinv) : UINT64_C(0);
    out->cycle_lcm = (out->is_bijective != 0U) ? phy_perm_order(out->p) : UINT32_C(0);

    out->anf_deg = phy_anf_degree4(out->p, &out->nl_bit, &out->nl_a, &out->nl_b);
    out->is_factorizable =
        phy_check_factorizable(out->p, &out->nf_axis, &out->nf_l, &out->nf_r0, &out->nf_r1)
            ? UINT8_C(1)
            : UINT8_C(0);
    out->yb_violations = phy_yb_check(out->p);

    uint8_t f00[16];
    uint8_t f0[16][16];
    for (uint8_t a = 0U; a < 16U; ++a) {
        f00[a] = phy_step3(out->p, 0U, 0U, a);
        for (uint8_t b = 0U; b < 16U; ++b) {
            f0[a][b] = phy_step3(out->p, 0U, a, b);
        }
    }

    uint16_t piv_triples[15];
    uint64_t rref[15 * 15];
    uint8_t has_piv[15];
    phy_vzero_u8((uint8_t *)piv_triples, sizeof(piv_triples));
    uint32_t rank = phy_stream_m61_pivots(out->p, f00, f0, piv_triples, rref, has_piv);
    out->rank_tot = (uint8_t)(15U - rank);

    if (out->rank_tot > 0U) {
        uint64_t null_basis[15 * 15];
        uint32_t n_null = phy_extract_m61_null_basis(rref, has_piv, null_basis);
        out->rank_ent = phy_count_entangled_dim_m61(null_basis, n_null);
        if (out->rank_ent > 0U) {
            (void)phy_select_best_int_charge(out->p, f00, f0, piv_triples, rank, out);
        }
    }

    uint32_t fl = PHY_RULE_FLAG_VALID;
    if (out->is_bijective != 0U) {
        fl |= PHY_RULE_FLAG_BIJECTIVE;
    }
    if (out->rank_ent > 0U && out->continuity_fail == 0 &&
        (out->sep_l0 != out->sep_l1 || out->sep_r0 != out->sep_r1)) {
        fl |= PHY_RULE_FLAG_ENTANGLED;
    }
    if (out->j_nonzero > 0U) {
        fl |= PHY_RULE_FLAG_ACTIVE_J;
    }
    if (out->anf_deg >= 2U) {
        fl |= PHY_RULE_FLAG_NONLINEAR;
    }
    if (out->is_factorizable == 0U && out->is_swap == 0U) {
        fl |= PHY_RULE_FLAG_NONFACT;
    }
    if (out->yb_violations == 0U && out->is_swap == 0U && out->is_id == 0U &&
        !phy_is_lyubashenko(out->p)) {
        fl |= PHY_RULE_FLAG_YANGBOX;
    }
    if ((fl & (PHY_RULE_FLAG_YANGBOX | PHY_RULE_FLAG_ENTANGLED | PHY_RULE_FLAG_NONLINEAR |
               PHY_RULE_FLAG_NONFACT)) == (PHY_RULE_FLAG_YANGBOX | PHY_RULE_FLAG_ENTANGLED |
                                           PHY_RULE_FLAG_NONLINEAR | PHY_RULE_FLAG_NONFACT)) {
        fl |= PHY_RULE_FLAG_JOINT;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    mh = phy_fnv1a64_u64(mh, out->lut_packed);
    mh = phy_fnv1a64_u64(mh, (uint64_t)fl);
    mh = phy_fnv1a64_u64(mh, ((uint64_t)out->rank_ent << 16U) | (uint64_t)out->j_nonzero);
    for (uint32_t u = 0U; u < 16U; ++u) {
        mh = phy_fnv1a64_u64(mh, (uint64_t)(uint32_t)out->q[u]);
    }
    out->merkle = mh;
}

/* --- 64-Site Periodic Floquet Lattice + FLOW 2-Jet Front Tracker --- */

static uint8_t phy_get_site64(uint64_t lo, uint64_t hi, uint32_t site) {
    uint32_t s = site & 63U;
    uint64_t w = (s < 32U) ? lo : hi;
    uint32_t sh = (s & 31U) * 2U;
    return (uint8_t)((w >> sh) & UINT64_C(3));
}

static void phy_set_site64(uint64_t *lo, uint64_t *hi, uint32_t site, uint8_t val) {
    uint32_t s = site & 63U;
    uint32_t sh = (s & 31U) * 2U;
    uint64_t clr = ~(UINT64_C(3) << sh);
    uint64_t ins = ((uint64_t)(val & 3U)) << sh;
    if (s < 32U) {
        *lo = (*lo & clr) | ins;
    } else {
        *hi = (*hi & clr) | ins;
    }
}

void phy_step_64ring(const uint8_t p[16], uint64_t lo, uint64_t hi, uint64_t *out_lo,
                     uint64_t *out_hi) {
    /* Even half-layer on blocks (2m, 2m+1) */
    uint64_t elo = UINT64_C(0);
    uint64_t ehi = UINT64_C(0);
    for (uint32_t m = 0U; m < 32U; ++m) {
        uint8_t l = phy_get_site64(lo, hi, 2U * m);
        uint8_t r = phy_get_site64(lo, hi, 2U * m + 1U);
        uint8_t y = p[phy_pack2(l, r)];
        phy_set_site64(&elo, &ehi, 2U * m, phy_hi2(y));
        phy_set_site64(&elo, &ehi, 2U * m + 1U, phy_lo2(y));
    }
    /* Odd half-layer on staggered blocks (2m+1, (2m+2) mod 64) */
    uint64_t olo = UINT64_C(0);
    uint64_t ohi = UINT64_C(0);
    for (uint32_t m = 0U; m < 32U; ++m) {
        uint32_t s0 = (2U * m + 1U) & 63U;
        uint32_t s1 = (2U * m + 2U) & 63U;
        uint8_t l = phy_get_site64(elo, ehi, s0);
        uint8_t r = phy_get_site64(elo, ehi, s1);
        uint8_t z = p[phy_pack2(l, r)];
        phy_set_site64(&olo, &ohi, s0, phy_hi2(z));
        phy_set_site64(&olo, &ohi, s1, phy_lo2(z));
    }
    *out_lo = olo;
    *out_hi = ohi;
}

int64_t phy_ring_charge_64(const int32_t q[16], uint64_t lo, uint64_t hi) {
    int64_t total = INT64_C(0);
    for (uint32_t m = 0U; m < 32U; ++m) {
        uint8_t l = phy_get_site64(lo, hi, 2U * m);
        uint8_t r = phy_get_site64(lo, hi, 2U * m + 1U);
        total += (int64_t)q[phy_pack2(l, r)];
    }
    return total;
}

void phy_sim_soliton_64(phy_arena_u *arena, const phy_rule_analysis_t *rule, uint64_t init_lo,
                        uint64_t init_hi, uint32_t steps, uint8_t defect_site, uint8_t defect_xor,
                        phy_sim_trace_t *out_trace) {
    if (rule == NULL || out_trace == NULL) {
        return;
    }
    phy_fiber_t fib;
    fib.state_lo = init_lo;
    fib.state_hi = init_hi;
    fib.pert_lo = init_lo;
    fib.pert_hi = init_hi;
    uint32_t ds = (uint32_t)(defect_site & 63U);
    uint8_t dv = (uint8_t)(phy_get_site64(init_lo, init_hi, ds) ^ (defect_xor & 3U));
    if (dv == phy_get_site64(init_lo, init_hi, ds)) {
        dv = (uint8_t)(dv ^ 1U);
    }
    phy_set_site64(&fib.pert_lo, &fib.pert_hi, ds, dv);
    fib.merkle = rule->merkle;
    fib.total_charge = phy_ring_charge_64(rule->q, fib.state_lo, fib.state_hi);
    fib.step_idx = 0U;
    fib.damage_count = 0U;
    fib.front_left = (uint8_t)ds;
    fib.front_right = (uint8_t)ds;

    phy_vzero_u64(fib.flow_state, 8U);
    fib.flow_state[1] = UINT64_C(65536); /* J1 = 1.0 Q16 for lossless unit-gain tracking */

    int64_t q0 = fib.total_charge;
    int64_t max_drift = INT64_C(0);
    uint32_t max_span = 0U;
    uint64_t vel_sum = UINT64_C(0);
    uint32_t frames = 0U;

    for (uint32_t t = 0U; t <= steps; ++t) {
        uint64_t diff_mask = UINT64_C(0);
        uint32_t dmg = 0U;
        uint32_t span = 0U;
        uint32_t fl = ds;
        uint32_t fr = ds;
        uint32_t best_dl = 0U;
        uint32_t best_dr = 0U;
        for (uint32_t s = 0U; s < 64U; ++s) {
            if (phy_get_site64(fib.state_lo, fib.state_hi, s) !=
                phy_get_site64(fib.pert_lo, fib.pert_hi, s)) {
                diff_mask |= (UINT64_C(1) << s);
                ++dmg;
                uint32_t d_right = (s - ds) & 63U;
                uint32_t d_left = (ds - s) & 63U;
                uint32_t cdist = (d_right < d_left) ? d_right : d_left;
                if (cdist > span) {
                    span = cdist;
                }
                if (d_right <= 32U && d_right >= best_dr) {
                    best_dr = d_right;
                    fr = s;
                }
                if (d_left <= 32U && d_left >= best_dl) {
                    best_dl = d_left;
                    fl = s;
                }
            }
        }
        int64_t cur_q = phy_ring_charge_64(rule->q, fib.state_lo, fib.state_hi);
        int64_t drift = cur_q - q0;
        uint64_t adrift = phy_abs_i64(drift);
        if ((int64_t)adrift > max_drift) {
            max_drift = (int64_t)adrift;
        }
        if (span > max_span) {
            max_span = span;
        }

        uint64_t bound[8];
        phy_vzero_u64(bound, 8U);
        int64_t span_q16 = (int64_t)((uint64_t)span << 16U);
        bound[0] = (uint64_t)span_q16;
        bound[3] =
            ((uint64_t)(UINT32_C(32) << 16U) << 32U) | (uint64_t)(uint32_t)(-(INT32_C(32) << 16));
        bound[4] = ((uint64_t)(UINT32_C(64) << 16U) << 32U);
        fib.flow_state[0] = fib.flow_state[3];
        fib.flow_state[1] = UINT64_C(65536);

        uint64_t flow_next[8];
        flow_step(flow_next, fib.flow_state, bound);
        phy_vcopy_u64(fib.flow_state, flow_next, 8U);
        if (t > 0U) {
            vel_sum += fib.flow_state[4];
        }

        fib.merkle = phy_fnv1a64_u64(fib.merkle, fib.state_lo ^ fib.state_hi ^ diff_mask);
        if (arena != NULL && t < PHY_COUNT_SIM_FRAMES) {
            uint32_t slot = PHY_SLOT_SIM_FRAMES + t;
            phy_vzero_u64(arena->cells[slot].u64, 8U);
            arena->map.sim_frames[t].state_lo = fib.state_lo;
            arena->map.sim_frames[t].state_hi = fib.state_hi;
            arena->map.sim_frames[t].pert_lo = fib.pert_lo;
            arena->map.sim_frames[t].pert_hi = fib.pert_hi;
            arena->map.sim_frames[t].diff_mask = diff_mask;
            arena->map.sim_frames[t].total_charge = cur_q;
            arena->map.sim_frames[t].step_idx = t;
            arena->map.sim_frames[t].damage_count = (uint16_t)dmg;
            arena->map.sim_frames[t].front_left = (uint8_t)fl;
            arena->map.sim_frames[t].front_right = (uint8_t)fr;
            arena->map.sim_frames[t].jet_width_lo = (uint32_t)fib.flow_state[3];
            arena->map.sim_frames[t].frame_merkle = (uint32_t)fib.merkle;
            ++frames;
        }

        if (t < steps) {
            uint64_t nlo = UINT64_C(0);
            uint64_t nhi = UINT64_C(0);
            phy_step_64ring(rule->p, fib.state_lo, fib.state_hi, &nlo, &nhi);
            fib.state_lo = nlo;
            fib.state_hi = nhi;

            uint64_t plo = UINT64_C(0);
            uint64_t phi = UINT64_C(0);
            phy_step_64ring(rule->p, fib.pert_lo, fib.pert_hi, &plo, &phi);
            fib.pert_lo = plo;
            fib.pert_hi = phi;
        }
    }

    out_trace->steps_run = steps;
    out_trace->max_span = max_span;
    out_trace->vel_num_q16 =
        (steps > 0U) ? (int32_t)phy_div_u64(vel_sum, (uint64_t)steps, NULL) : 0;
    out_trace->initial_charge = q0;
    out_trace->final_charge = phy_ring_charge_64(rule->q, fib.state_lo, fib.state_hi);
    out_trace->max_charge_drift = max_drift;
    out_trace->final_merkle = fib.merkle;
    out_trace->frames_recorded = frames;
}

/* --- Arena Initialization & SEARU Cayley Graph Discovery Engine --- */

void phy_arena_init(phy_arena_u *arena, uint64_t seed) {
    if (arena == NULL) {
        return;
    }
    for (uint32_t slot = 0U; slot < PHY_TOTAL_CELLS; ++slot) {
        phy_vzero_u64(arena->cells[slot].u64, 8U);
    }
    arena->map.img_hdr.magic = PHY_MAGIC;
    arena->map.img_hdr.version = PHY_VERSION;
    arena->map.img_hdr.total_cells = PHY_TOTAL_CELLS;
    arena->map.img_hdr.active_rules = 0U;
    arena->map.img_hdr.merkle_root = PHY_SCHEMA_HASH;
    arena->map.img_hdr.seed = seed;
}

static bool phy_rule_matches_mode(uint32_t mode, const phy_rule_analysis_t *cand) {
    uint32_t req_ent = PHY_RULE_FLAG_BIJECTIVE | PHY_RULE_FLAG_ENTANGLED | PHY_RULE_FLAG_ACTIVE_J |
                       PHY_RULE_FLAG_NONLINEAR | PHY_RULE_FLAG_NONFACT;
    if (mode == PHY_HUNT_ENTANGLED) {
        return (cand->flags & req_ent) == req_ent;
    }
    if (mode == PHY_HUNT_YANGBOX) {
        uint32_t req_yb = PHY_RULE_FLAG_BIJECTIVE | PHY_RULE_FLAG_YANGBOX |
                          PHY_RULE_FLAG_NONLINEAR | PHY_RULE_FLAG_NONFACT;
        return (cand->flags & req_yb) == req_yb;
    }
    return (cand->flags & PHY_RULE_FLAG_JOINT) == PHY_RULE_FLAG_JOINT;
}

static uint64_t phy_score_rule(uint32_t mode, const phy_rule_analysis_t *cand) {
    if (mode == PHY_HUNT_ENTANGLED) {
        uint64_t sc = ((uint64_t)cand->rank_ent << 40U) | ((uint64_t)cand->anf_deg << 32U) |
                      ((uint64_t)cand->j_nonzero << 16U) | (uint64_t)cand->rank_tot;
        return sc;
    }
    if (mode == PHY_HUNT_YANGBOX) {
        uint64_t non_triv = (uint64_t)(cand->is_swap == 0U && cand->is_id == 0U) * UINT64_C(1);
        uint64_t non_ly = (uint64_t)(!phy_is_lyubashenko(cand->p)) * UINT64_C(1);
        return (non_triv << 48U) | (non_ly << 44U) | ((uint64_t)cand->anf_deg << 32U) |
               ((uint64_t)cand->rank_ent << 24U) | (uint64_t)cand->cycle_lcm;
    }
    uint64_t active = (uint64_t)(cand->j_nonzero > 0U) * UINT64_C(1);
    return (active << 48U) | ((uint64_t)cand->rank_ent << 40U) | ((uint64_t)cand->anf_deg << 32U) |
           ((uint64_t)cand->j_nonzero << 16U) | (uint64_t)cand->rank_tot;
}

static void phy_store_arena_rule(phy_arena_u *arena, uint32_t idx, const phy_rule_analysis_t *rule,
                                 const phy_sim_trace_t *sim) {
    if (arena == NULL || idx >= PHY_COUNT_RULES) {
        return;
    }
    uint32_t rslot = PHY_SLOT_RULES + idx;
    phy_vzero_u64(arena->cells[rslot].u64, 8U);
    arena->map.rules[idx].lut_packed = rule->lut_packed;
    arena->map.rules[idx].pinv_packed = rule->pinv_packed;
    arena->map.rules[idx].merkle = rule->merkle;
    arena->map.rules[idx].seed = rule->seed;
    arena->map.rules[idx].flags = rule->flags;
    arena->map.rules[idx].anf_deg = rule->anf_deg;
    arena->map.rules[idx].yb_violations = rule->yb_violations;
    arena->map.rules[idx].rank_tot = rule->rank_tot;
    arena->map.rules[idx].rank_ent = rule->rank_ent;
    arena->map.rules[idx].j_nonzero = rule->j_nonzero;
    arena->map.rules[idx].mc_and_count = 0U;
    arena->map.rules[idx].cycle_lcm = rule->cycle_lcm;
    arena->map.rules[idx].soliton_vel_num = (sim != NULL) ? sim->vel_num_q16 : 0;
    arena->map.rules[idx].soliton_max_span = (sim != NULL) ? sim->max_span : 0U;

    uint32_t cslot = PHY_SLOT_CHARGES + idx;
    phy_vzero_u64(arena->cells[cslot].u64, 8U);
    phy_vcopy_i32(arena->map.charges[idx].q, rule->q, 16U);

    if (idx < 8U) {
        uint32_t j_base = idx * 16U;
        for (uint32_t r = 0U; r < 16U; ++r) {
            uint32_t slot = PHY_SLOT_CURRENTS + j_base + r;
            phy_vzero_u64(arena->cells[slot].u64, 8U);
            phy_vcopy_i32(arena->map.currents[j_base + r].j, rule->j[r], 16U);
        }
    }
}

static void phy_build_defect_sector_p(uint64_t *rng, uint8_t out_p[16]) {
    uint8_t d = (uint8_t)(1U + searu_Bounded_u32(searu_rng_next(rng), 3U));
    uint8_t others[3];
    uint32_t oi = 0U;
    for (uint8_t c = 0U; c < 4U; ++c) {
        if (c != d && oi < 3U) {
            others[oi++] = c;
        }
    }
    uint8_t inv_idx[4];
    phy_vzero_u8(inv_idx, 4U);
    for (uint32_t k = 0U; k < 3U; ++k) {
        inv_idx[others[k]] = (uint8_t)k;
    }

    uint8_t sig_l[3];
    uint8_t sig_r[3];
    for (uint32_t k = 0U; k < 3U; ++k) {
        sig_l[k] = others[k];
        sig_r[k] = others[k];
    }
    for (uint32_t k = 2U; k > 0U; --k) {
        uint32_t jl = searu_Bounded_u32(searu_rng_next(rng), k + 1U);
        uint8_t tl = sig_l[k];
        sig_l[k] = sig_l[jl];
        sig_l[jl] = tl;
        uint32_t jr = searu_Bounded_u32(searu_rng_next(rng), k + 1U);
        uint8_t tr = sig_r[k];
        sig_r[k] = sig_r[jr];
        sig_r[jr] = tr;
    }

    uint8_t b_l[3];
    uint8_t b_r[3];
    for (uint32_t k = 0U; k < 3U; ++k) {
        b_l[k] = others[k];
        b_r[k] = others[k];
    }
    if ((searu_rng_next(rng) & UINT64_C(1)) != UINT64_C(0)) {
        uint8_t t = b_l[1];
        b_l[1] = b_l[2];
        b_l[2] = t;
    }
    if ((searu_rng_next(rng) & UINT64_C(1)) != UINT64_C(0)) {
        uint8_t t = b_r[1];
        b_r[1] = b_r[2];
        b_r[2] = t;
    }

    for (uint8_t a = 0U; a < 4U; ++a) {
        for (uint8_t b = 0U; b < 4U; ++b) {
            uint8_t u = phy_pack2(a, b);
            if (a == d && b == d) {
                out_p[u] = phy_pack2(d, d);
            } else if (a == d) {
                out_p[u] = phy_pack2(d, sig_l[inv_idx[b]]);
            } else if (b == d) {
                out_p[u] = phy_pack2(sig_r[inv_idx[a]], d);
            } else {
                out_p[u] = phy_pack2(b_l[inv_idx[b]], b_r[inv_idx[a]]);
            }
        }
    }
}

static void phy_pick_block_swap(uint64_t *rng, bool hi_blocks, uint8_t p[16]) {
    uint32_t blk = searu_Bounded_u32(searu_rng_next(rng), 4U);
    uint32_t m1 = (blk >> 1U) & UINT32_C(1);
    uint32_t m2 = blk & UINT32_C(1);
    uint8_t elems[4];
    uint32_t cnt = 0U;
    for (uint32_t x = 0U; x < 2U; ++x) {
        for (uint32_t y = 0U; y < 2U; ++y) {
            uint32_t a = hi_blocks ? ((m1 << 1U) | x) : ((x << 1U) | m1);
            uint32_t b = hi_blocks ? ((m2 << 1U) | y) : ((y << 1U) | m2);
            uint8_t u = phy_pack2((uint8_t)a, (uint8_t)b);
            if (u != 0U && cnt < 4U) {
                elems[cnt++] = u;
            }
        }
    }
    if (cnt < 2U) {
        return;
    }
    uint32_t i = searu_Bounded_u32(searu_rng_next(rng), cnt);
    uint32_t j = searu_Bounded_u32(searu_rng_next(rng), cnt);
    if (i == j) {
        j = (j + 1U == cnt) ? 0U : (j + 1U);
    }
    uint8_t tmp = p[elems[i]];
    p[elems[i]] = p[elems[j]];
    p[elems[j]] = tmp;
}

uint32_t phy_hunt(phy_arena_u *arena, uint32_t mode, uint64_t seed, uint64_t max_trials,
                  phy_rule_analysis_t *out_best) {
    if (out_best == NULL) {
        return 0U;
    }
    phy_zero_rule_analysis(out_best);
    if (arena != NULL) {
        phy_arena_init(arena, seed);
    }

    uint8_t base_swap[16];
    for (uint8_t u = 0U; u < 16U; ++u) {
        base_swap[u] = phy_pack2(phy_lo2(u), phy_hi2(u));
    }

    uint64_t rng = seed ^ UINT64_C(0x9E3779B97F4A7C15);
    if (rng == UINT64_C(0)) {
        rng = UINT64_C(0xA5A5A5A55A5A5A5A);
    }

    uint64_t seen_luts[16];
    phy_vzero_u64(seen_luts, 16U);
    uint32_t found_count = 0U;
    uint64_t best_score = UINT64_C(0);
    bool have_best = false;

    phy_rule_analysis_t cand;
    uint8_t cur_p[16];
    phy_vcopy_u8(cur_p, base_swap, 16U);

    uint64_t trials = 0U;
    uint64_t bijective_tested = 0U;
    uint64_t entangled_found = 0U;
    uint64_t yb_found = 0U;
    uint64_t joint_found = 0U;
    uint64_t rej_sep = 0U;
    uint64_t rej_lin = 0U;

    if (mode == PHY_HUNT_ENTANGLED) {
        for (uint64_t t = 0U; t < max_trials; ++t) {
            ++trials;
            uint8_t trial_p[16];
            if ((t & UINT64_C(1)) == UINT64_C(0)) {
                phy_build_defect_sector_p(&rng, trial_p);
            } else {
                if ((t & UINT64_C(63)) == UINT64_C(1)) {
                    phy_vcopy_u8(cur_p, base_swap, 16U);
                }
                phy_vcopy_u8(trial_p, cur_p, 16U);
                uint64_t bswan = UINT64_C(0);
                uint32_t jv = (uint32_t)(searu_zipf_jump(&rng, 48U, &bswan) & UINT64_C(3));
                uint32_t n_swaps = (jv == 3U) ? 2U : (1U + jv);
                for (uint32_t k = 0U; k < n_swaps; ++k) {
                    uint32_t i = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
                    uint32_t j = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
                    if (i == j) {
                        j = (j == 15U) ? 1U : (j + 1U);
                    }
                    uint8_t tmp = trial_p[i];
                    trial_p[i] = trial_p[j];
                    trial_p[j] = tmp;
                }
            }
            ++bijective_tested;
            phy_analyze_rule(trial_p, seed + t, &cand);
            if (cand.anf_deg < 2U) {
                ++rej_lin;
            } else if (cand.rank_tot > 0U && cand.rank_ent == 0U) {
                ++rej_sep;
            }
            if ((cand.flags & (PHY_RULE_FLAG_ENTANGLED | PHY_RULE_FLAG_ACTIVE_J)) ==
                (PHY_RULE_FLAG_ENTANGLED | PHY_RULE_FLAG_ACTIVE_J)) {
                ++entangled_found;
            }
            if (phy_rule_matches_mode(mode, &cand)) {
                bool dup = false;
                for (uint32_t s = 0U; s < found_count; ++s) {
                    if (seen_luts[s] == cand.lut_packed) {
                        dup = true;
                        break;
                    }
                }
                if (!dup && found_count < PHY_COUNT_RULES) {
                    seen_luts[found_count] = cand.lut_packed;
                    phy_sim_trace_t tr;
                    phy_sim_soliton_64((found_count == 0U) ? arena : NULL, &cand,
                                       UINT64_C(0x123456789ABCDEF0), UINT64_C(0x0FEDCBA987654321),
                                       32U, 32U, 3U, &tr);
                    phy_store_arena_rule(arena, found_count, &cand, &tr);
                    ++found_count;
                }
                uint64_t sc = phy_score_rule(mode, &cand);
                if (!have_best || sc > best_score) {
                    have_best = true;
                    best_score = sc;
                    phy_vcopy_u8(cur_p, trial_p, 16U);
                    phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
                }
            }
        }
    } else {
        /* YANGBOX / JOINT mode: multi-seed SA across ID, SWAP, PSWAP, and PSWAP_HI with
           bit-channel block-fibered moves in JOINT mode */
        uint8_t base_pswap[16];
        uint8_t base_pswap_hi[16];
        for (uint8_t u = 0U; u < 16U; ++u) {
            uint32_t h = (uint32_t)phy_hi2(u);
            uint32_t l = (uint32_t)phy_lo2(u);
            base_pswap[u] = phy_pack2((uint8_t)((h & UINT32_C(2)) | (l & UINT32_C(1))),
                                      (uint8_t)((l & UINT32_C(2)) | (h & UINT32_C(1))));
            base_pswap_hi[u] = phy_pack2((uint8_t)((l & UINT32_C(2)) | (h & UINT32_C(1))),
                                         (uint8_t)((h & UINT32_C(2)) | (l & UINT32_C(1))));
        }
        uint64_t restarts = max_trials >> 6U;
        if (restarts < 32U) {
            restarts = 32U;
        }
        uint64_t budget = max_trials;
        uint32_t rs_mod3 = 0U;
        for (uint64_t rs = 0U; rs < restarts && trials < budget; ++rs) {
            bool fibered_joint = (mode == PHY_HUNT_JOINT) && ((rs & UINT64_C(1)) == UINT64_C(0));
            bool use_hi_blk = ((rs & UINT64_C(2)) == UINT64_C(0));
            if (fibered_joint) {
                phy_vcopy_u8(cur_p, use_hi_blk ? base_pswap : base_pswap_hi, 16U);
                uint32_t kick = 1U + searu_Bounded_u32(searu_rng_next(&rng), 3U);
                for (uint32_t k = 0U; k < kick; ++k) {
                    phy_pick_block_swap(&rng, use_hi_blk, cur_p);
                }
            } else {
                if (rs_mod3 == 0U) {
                    for (uint8_t u = 0U; u < 16U; ++u) {
                        cur_p[u] = u;
                    }
                } else if (rs_mod3 == 1U) {
                    phy_vcopy_u8(cur_p, base_swap, 16U);
                } else {
                    phy_vcopy_u8(cur_p, base_pswap, 16U);
                }
                rs_mod3 = (rs_mod3 == 2U) ? 0U : (rs_mod3 + 1U);
                uint32_t kick = 1U + searu_Bounded_u32(searu_rng_next(&rng), 4U);
                for (uint32_t k = 0U; k < kick; ++k) {
                    uint32_t i = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
                    uint32_t j = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
                    if (i == j) {
                        j = (j == 15U) ? 1U : (j + 1U);
                    }
                    uint8_t tmp = cur_p[i];
                    cur_p[i] = cur_p[j];
                    cur_p[j] = tmp;
                }
            }
            uint8_t cur_yb = phy_yb_check(cur_p);
            uint32_t step_lim = fibered_joint ? 100U : 180U;
            for (uint32_t step = 0U; step < step_lim && trials < budget; ++step) {
                ++trials;
                if (cur_yb == 0U) {
                    ++bijective_tested;
                    phy_analyze_rule(cur_p, seed + rs * UINT64_C(1000) + step, &cand);
                    if ((cand.flags & PHY_RULE_FLAG_YANGBOX) != 0U) {
                        ++yb_found;
                    }
                    if ((cand.flags & PHY_RULE_FLAG_JOINT) != 0U) {
                        ++joint_found;
                    }
                    if (phy_rule_matches_mode(mode, &cand)) {
                        bool dup = false;
                        for (uint32_t s = 0U; s < found_count; ++s) {
                            if (seen_luts[s] == cand.lut_packed) {
                                dup = true;
                                break;
                            }
                        }
                        if (!dup && found_count < PHY_COUNT_RULES) {
                            seen_luts[found_count] = cand.lut_packed;
                            phy_sim_trace_t tr;
                            phy_sim_soliton_64((found_count == 0U) ? arena : NULL, &cand,
                                               UINT64_C(0x123456789ABCDEF0),
                                               UINT64_C(0x0FEDCBA987654321), 32U, 32U, 3U, &tr);
                            phy_store_arena_rule(arena, found_count, &cand, &tr);
                            ++found_count;
                        }
                        uint64_t sc = phy_score_rule(mode, &cand);
                        if (!have_best || sc > best_score) {
                            have_best = true;
                            best_score = sc;
                            phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
                        }
                    }
                    break;
                }
                uint8_t saved_p[16];
                if (fibered_joint) {
                    phy_vcopy_u8(saved_p, cur_p, 16U);
                    phy_pick_block_swap(&rng, use_hi_blk, cur_p);
                } else {
                    uint32_t i = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
                    uint32_t j = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
                    if (i == j) {
                        j = (j == 15U) ? 1U : (j + 1U);
                    }
                    phy_vcopy_u8(saved_p, cur_p, 16U);
                    uint8_t tmp = cur_p[i];
                    cur_p[i] = cur_p[j];
                    cur_p[j] = tmp;
                }
                uint8_t nxt_yb = phy_yb_check(cur_p);
                uint32_t accept_r = searu_Bounded_u32(searu_rng_next(&rng), 100U);
                if (nxt_yb <= cur_yb || accept_r < 8U) {
                    cur_yb = nxt_yb;
                } else {
                    phy_vcopy_u8(cur_p, saved_p, 16U);
                }
            }
        }
    }

    if (have_best && arena != NULL) {
        phy_sim_trace_t tr;
        phy_sim_soliton_64(arena, out_best, UINT64_C(0x123456789ABCDEF0),
                           UINT64_C(0x0FEDCBA987654321), 32U, 32U, 3U, &tr);
        arena->map.img_hdr.active_rules = found_count;
        arena->map.img_hdr.merkle_root = out_best->merkle;
        arena->map.hunt_stats.trials_done = trials;
        arena->map.hunt_stats.bijective_tested = bijective_tested;
        arena->map.hunt_stats.entangled_found = entangled_found;
        arena->map.hunt_stats.yangbaxter_found = yb_found;
        arena->map.hunt_stats.joint_found = joint_found;
        arena->map.hunt_stats.rejected_separable = rej_sep;
        arena->map.hunt_stats.rejected_linear = rej_lin;
        arena->map.hunt_stats.best_rule_merkle = out_best->merkle;
    } else if (arena != NULL) {
        arena->map.hunt_stats.trials_done = trials;
        arena->map.hunt_stats.bijective_tested = bijective_tested;
        arena->map.hunt_stats.entangled_found = entangled_found;
        arena->map.hunt_stats.yangbaxter_found = yb_found;
        arena->map.hunt_stats.joint_found = joint_found;
        arena->map.hunt_stats.rejected_separable = rej_sep;
        arena->map.hunt_stats.rejected_linear = rej_lin;
    }

    return found_count;
}

/* ============================================================================
 * PHY v2.0: Multi-Point Butterfly Velocity, 4-Site Q4 Hierarchy & |X|=8 YBE
 * ========================================================================== */

void phy_measure_butterfly_64(const phy_rule_analysis_t *rule, uint64_t init_lo, uint64_t init_hi,
                              uint32_t steps, phy_butterfly_t *out) {
    if (rule == NULL || out == NULL) {
        return;
    }
    out->max_span = 0U;
    out->mean_span_q16 = 0U;
    out->max_vel_q16 = 0;
    out->mean_vel_q16 = 0;
    out->active_probes = 0U;
    out->max_charge_drift = 0;

    uint64_t span_sum = UINT64_C(0);
    int64_t vel_sum = INT64_C(0);
    for (uint32_t k = 0U; k < 8U; ++k) {
        uint8_t site = (uint8_t)(k * 8U);
        uint32_t xm = k & UINT32_C(3);
        uint8_t dxor = (xm == UINT32_C(0)) ? UINT8_C(1) : (uint8_t)xm;
        phy_sim_trace_t tr;
        phy_sim_soliton_64(NULL, rule, init_lo, init_hi, steps, site, dxor, &tr);
        if (tr.max_span > 0U) {
            out->active_probes += 1U;
        }
        if (tr.max_span > out->max_span) {
            out->max_span = tr.max_span;
        }
        if (tr.vel_num_q16 > out->max_vel_q16) {
            out->max_vel_q16 = tr.vel_num_q16;
        }
        if (tr.max_charge_drift > out->max_charge_drift) {
            out->max_charge_drift = tr.max_charge_drift;
        }
        span_sum += (uint64_t)tr.max_span;
        vel_sum += (int64_t)tr.vel_num_q16;
    }
    out->mean_span_q16 = (uint32_t)phy_div_u64(span_sum << 16U, UINT64_C(8), NULL);
    out->mean_vel_q16 = (int32_t)phy_div_s64(vel_sum, UINT64_C(8));
}

static uint64_t phy_q4_canon_j(const uint8_t f3[16][16][16], const int32_t q4[256], uint8_t a,
                               uint8_t b, uint8_t c) {
    uint8_t f000 = f3[0][0][0];
    uint8_t f00a = f3[0][0][a];
    uint8_t f0ab = f3[0][a][b];
    uint8_t fabc = f3[a][b][c];
    int64_t t1 = (int64_t)q4[((uint32_t)f000 << 4U) | (uint32_t)f00a] - (int64_t)q4[0];
    int64_t t2 = (int64_t)q4[((uint32_t)f00a << 4U) | (uint32_t)f0ab] - (int64_t)q4[a];
    int64_t t3 =
        (int64_t)q4[((uint32_t)f0ab << 4U) | (uint32_t)fabc] - (int64_t)q4[((uint32_t)a << 4U) | b];
    return (uint64_t)(t1 + t2 + t3);
}

static int32_t phy_q4_verify_continuity_i64(const uint8_t f3[16][16][16], const int32_t q4[256],
                                            uint16_t *out_j4_nz, uint8_t *wa, uint8_t *wb,
                                            uint8_t *wc) {
    int32_t fails = 0;
    uint16_t j_nz = 0U;
    bool got_j_wit = false;
    for (uint8_t a = 0U; a < 16U; ++a) {
        for (uint8_t b = 0U; b < 16U; ++b) {
            int64_t jab = (int64_t)phy_q4_canon_j(f3, q4, 0U, a, b);
            (void)jab;
            for (uint8_t c = 0U; c < 16U; ++c) {
                int64_t jabc = (int64_t)phy_q4_canon_j(f3, q4, a, b, c);
                if (jabc != INT64_C(0)) {
                    ++j_nz;
                    if (!got_j_wit) {
                        got_j_wit = true;
                        *wa = a;
                        *wb = b;
                        *wc = c;
                    }
                }
                for (uint8_t d = 0U; d < 16U; ++d) {
                    uint8_t u1 = f3[a][b][c];
                    uint8_t u2 = f3[b][c][d];
                    int64_t lhs = (int64_t)q4[((uint32_t)u1 << 4U) | (uint32_t)u2] -
                                  (int64_t)q4[((uint32_t)b << 4U) | (uint32_t)c];
                    int64_t jbcd = (int64_t)phy_q4_canon_j(f3, q4, b, c, d);
                    int64_t rhs = jbcd - jabc;
                    if (lhs != rhs) {
                        ++fails;
                    }
                }
            }
        }
    }
    *out_j4_nz = j_nz;
    return fails;
}

static bool phy_q4_find_nontriv_witness(const int32_t q4[256], uint8_t *axis, uint8_t *f1,
                                        uint8_t *f2, uint8_t *x0, uint8_t *x1, uint8_t *y0,
                                        uint8_t *y1) {
    /* Inner cross-difference with fixed outer (aL, bR) detects genuine 4-site coupling */
    for (uint8_t aL = 0U; aL < 4U; ++aL) {
        for (uint8_t bR = 0U; bR < 4U; ++bR) {
            for (uint8_t r0 = 0U; r0 < 4U; ++r0) {
                for (uint8_t r1 = (uint8_t)(r0 + 1U); r1 < 4U; ++r1) {
                    for (uint8_t l0 = 0U; l0 < 4U; ++l0) {
                        for (uint8_t l1 = (uint8_t)(l0 + 1U); l1 < 4U; ++l1) {
                            uint32_t u0 = phy_pack2(aL, r0);
                            uint32_t u1 = phy_pack2(aL, r1);
                            uint32_t v0 = phy_pack2(l0, bR);
                            uint32_t v1 = phy_pack2(l1, bR);
                            int64_t cd =
                                (int64_t)q4[(u1 << 4U) | v1] - (int64_t)q4[(u1 << 4U) | v0] -
                                (int64_t)q4[(u0 << 4U) | v1] + (int64_t)q4[(u0 << 4U) | v0];
                            if (cd != INT64_C(0)) {
                                *axis = 0U;
                                *f1 = aL;
                                *f2 = bR;
                                *x0 = r0;
                                *x1 = r1;
                                *y0 = l0;
                                *y1 = l1;
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    /* Outer cross-difference with fixed inner (aR, bL) */
    for (uint8_t aR = 0U; aR < 4U; ++aR) {
        for (uint8_t bL = 0U; bL < 4U; ++bL) {
            for (uint8_t l0 = 0U; l0 < 4U; ++l0) {
                for (uint8_t l1 = (uint8_t)(l0 + 1U); l1 < 4U; ++l1) {
                    for (uint8_t r0 = 0U; r0 < 4U; ++r0) {
                        for (uint8_t r1 = (uint8_t)(r0 + 1U); r1 < 4U; ++r1) {
                            uint32_t u0 = phy_pack2(l0, aR);
                            uint32_t u1 = phy_pack2(l1, aR);
                            uint32_t v0 = phy_pack2(bL, r0);
                            uint32_t v1 = phy_pack2(bL, r1);
                            int64_t cd =
                                (int64_t)q4[(u1 << 4U) | v1] - (int64_t)q4[(u1 << 4U) | v0] -
                                (int64_t)q4[(u0 << 4U) | v1] + (int64_t)q4[(u0 << 4U) | v0];
                            if (cd != INT64_C(0)) {
                                *axis = 1U;
                                *f1 = aR;
                                *f2 = bL;
                                *x0 = l0;
                                *x1 = l1;
                                *y0 = r0;
                                *y1 = r1;
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

static bool phy_q4_ech_insert(uint64_t *ech, uint8_t *has, uint64_t row[255]) {
    for (uint32_t col = 0U; col < 255U; ++col) {
        uint64_t val = row[col];
        if (val != UINT64_C(0)) {
            if (has[col] != 0U) {
                const uint64_t *prow = &ech[(size_t)col * 255U];
                for (uint32_t j = col; j < 255U; ++j) {
                    uint64_t pj = prow[j];
                    if (pj != UINT64_C(0)) {
                        row[j] = phy_m61_sub(row[j], phy_m61_mul(val, pj));
                    }
                }
            } else {
                uint64_t inv = phy_m61_inv(val);
                uint64_t *dst = &ech[(size_t)col * 255U];
                for (uint32_t j = col; j < 255U; ++j) {
                    uint64_t rj = row[j];
                    dst[j] = (rj != UINT64_C(0)) ? phy_m61_mul(rj, inv) : UINT64_C(0);
                }
                has[col] = 1U;
                return true;
            }
        }
    }
    return false;
}

static void phy_q4_add_term(uint64_t row[255], uint32_t idx, int64_t s) {
    if (idx == 0U) {
        return;
    }
    uint32_t col = idx - 1U;
    if (s > 0) {
        row[col] = phy_m61_add(row[col], (uint64_t)s);
    } else if (s < 0) {
        row[col] = phy_m61_sub(row[col], (uint64_t)(-s));
    }
}

static void phy_q4_build_trivial_basis(phy_q4_ws_t *ws) {
    phy_vzero_u64(ws->triv_ech, (size_t)255U * 255U);
    phy_vzero_u8(ws->triv_has, 255U);

    /* 1. Shift-gauge: q(u,w) = delta_k(w) - delta_k(u) for k=1..15 */
    for (uint8_t k = 1U; k < 16U; ++k) {
        phy_vzero_u64(ws->work_row, 255U);
        for (uint8_t u = 0U; u < 16U; ++u) {
            for (uint8_t w = 0U; w < 16U; ++w) {
                int64_t s =
                    ((w == k) ? INT64_C(1) : INT64_C(0)) - ((u == k) ? INT64_C(1) : INT64_C(0));
                if (s != INT64_C(0)) {
                    phy_q4_add_term(ws->work_row, ((uint32_t)u << 4U) | (uint32_t)w, s);
                }
            }
        }
        (void)phy_q4_ech_insert(ws->triv_ech, ws->triv_has, ws->work_row);
    }

    /* 2. All 2-site densities q2(u) */
    for (uint8_t k = 1U; k < 16U; ++k) {
        phy_vzero_u64(ws->work_row, 255U);
        for (uint8_t w = 0U; w < 16U; ++w) {
            phy_q4_add_term(ws->work_row, ((uint32_t)k << 4U) | (uint32_t)w, INT64_C(1));
        }
        (void)phy_q4_ech_insert(ws->triv_ech, ws->triv_has, ws->work_row);
    }

    /* 3. Right-chiral 2-site pair correlations h(hi2(u), hi2(w)) */
    for (uint8_t h1 = 1U; h1 < 4U; ++h1) {
        for (uint8_t h2 = 1U; h2 < 4U; ++h2) {
            phy_vzero_u64(ws->work_row, 255U);
            for (uint8_t u = 0U; u < 16U; ++u) {
                if (phy_hi2(u) != h1) {
                    continue;
                }
                for (uint8_t w = 0U; w < 16U; ++w) {
                    if (phy_hi2(w) == h2) {
                        phy_q4_add_term(ws->work_row, ((uint32_t)u << 4U) | (uint32_t)w,
                                        INT64_C(1));
                    }
                }
            }
            (void)phy_q4_ech_insert(ws->triv_ech, ws->triv_has, ws->work_row);
        }
    }

    /* 4. Left-chiral 2-site pair correlations k(lo2(u), lo2(w)) */
    for (uint8_t l1 = 1U; l1 < 4U; ++l1) {
        for (uint8_t l2 = 1U; l2 < 4U; ++l2) {
            phy_vzero_u64(ws->work_row, 255U);
            for (uint8_t u = 0U; u < 16U; ++u) {
                if (phy_lo2(u) != l1) {
                    continue;
                }
                for (uint8_t w = 0U; w < 16U; ++w) {
                    if (phy_lo2(w) == l2) {
                        phy_q4_add_term(ws->work_row, ((uint32_t)u << 4U) | (uint32_t)w,
                                        INT64_C(1));
                    }
                }
            }
            (void)phy_q4_ech_insert(ws->triv_ech, ws->triv_has, ws->work_row);
        }
    }
}

static bool phy_q4_m61_to_small_int(const uint64_t mvec[255], int32_t out_q4[256],
                                    uint16_t *out_nnz) {
    uint64_t half_m = PHY_M61 >> 1U;
    for (uint64_t scale = UINT64_C(1); scale <= UINT64_C(64); ++scale) {
        int64_t sv[256];
        sv[0] = INT64_C(0);
        uint64_t g = UINT64_C(0);
        bool ok = true;
        for (uint32_t c = 0U; c < 255U; ++c) {
            uint64_t val = (scale == UINT64_C(1)) ? mvec[c] : phy_m61_mul(mvec[c], scale);
            int64_t s = (val <= half_m) ? (int64_t)val : -(int64_t)(PHY_M61 - val);
            uint64_t av = phy_abs_i64(s);
            if (av > UINT64_C(10000000)) {
                ok = false;
                break;
            }
            sv[c + 1U] = s;
            if (av != UINT64_C(0)) {
                g = phy_stein_gcd64(g, av);
            }
        }
        if (ok && g > UINT64_C(0)) {
            uint16_t nnz = 0U;
            out_q4[0] = 0;
            for (uint32_t i = 1U; i < 256U; ++i) {
                int64_t nv = (g > UINT64_C(1)) ? phy_div_s64(sv[i], g) : sv[i];
                out_q4[i] = (int32_t)nv;
                if (nv != INT64_C(0)) {
                    ++nnz;
                }
            }
            *out_nnz = nnz;
            return nnz > 0U;
        }
    }
    return false;
}

void phy_analyze_q4(const uint8_t p[16], phy_q4_ws_t *ws, phy_q4_analysis_t *out) {
    if (p == NULL || ws == NULL || out == NULL) {
        return;
    }
    phy_vcopy_u8(out->p, p, 16U);
    phy_vzero_i32(out->q4, 256U);
    out->lut_packed = phy_pack_lut(out->p);
    out->q4_merkle = UINT64_C(0);
    out->tot2 = 0U;
    out->tot4 = 0U;
    out->irred4 = 0U;
    out->genuine_q4 = 0U;
    out->j4_nonzero = 0U;
    out->q4_nnz = 0U;
    out->continuity4_fail = 0;
    out->cross_axis = 0U;
    out->cross_f1 = 0U;
    out->cross_f2 = 0U;
    out->cross_x0 = 0U;
    out->cross_x1 = 0U;
    out->cross_y0 = 0U;
    out->cross_y1 = 0U;
    out->j4_wit_a = 0U;
    out->j4_wit_b = 0U;
    out->j4_wit_c = 0U;

    /* Precompute 3-block transition table f3[a][b][c] */
    for (uint8_t a = 0U; a < 16U; ++a) {
        for (uint8_t b = 0U; b < 16U; ++b) {
            for (uint8_t c = 0U; c < 16U; ++c) {
                ws->f3[a][b][c] = phy_step3(out->p, a, b, c);
            }
        }
    }

    /* Count 2-site tot2 exactly */
    uint8_t f00[16];
    uint8_t f0[16][16];
    for (uint8_t a = 0U; a < 16U; ++a) {
        f00[a] = ws->f3[0][0][a];
        for (uint8_t b = 0U; b < 16U; ++b) {
            f0[a][b] = ws->f3[0][a][b];
        }
    }
    uint16_t piv2[15];
    uint64_t rref2[15 * 15];
    uint8_t has2[15];
    phy_vzero_u8((uint8_t *)piv2, sizeof(piv2));
    uint32_t rk2 = phy_stream_m61_pivots(out->p, f00, f0, piv2, rref2, has2);
    out->tot2 = (uint16_t)(15U - rk2);

    /* Stream all 65,536 4-tuples (a, b, c, d) into 255-col M61 RREF */
    phy_vzero_u64(ws->rref, (size_t)255U * 255U);
    phy_vzero_u8(ws->has_piv, 255U);
    uint8_t f000 = ws->f3[0][0][0];
    uint32_t rank4 = 0U;

    for (uint8_t a = 0U; a < 16U; ++a) {
        uint32_t ta1 = ((uint32_t)f000 << 4U) | (uint32_t)ws->f3[0][0][a];
        for (uint8_t b = 0U; b < 16U; ++b) {
            uint32_t tb1 = ((uint32_t)f000 << 4U) | (uint32_t)ws->f3[0][0][b];
            uint32_t ta2 = ((uint32_t)ws->f3[0][0][a] << 4U) | (uint32_t)ws->f3[0][a][b];
            uint32_t ab = ((uint32_t)a << 4U) | (uint32_t)b;
            for (uint8_t c = 0U; c < 16U; ++c) {
                uint32_t tb2 = ((uint32_t)ws->f3[0][0][b] << 4U) | (uint32_t)ws->f3[0][b][c];
                uint32_t ta3 = ((uint32_t)ws->f3[0][a][b] << 4U) | (uint32_t)ws->f3[a][b][c];
                uint32_t fabc = (uint32_t)ws->f3[a][b][c];
                uint32_t f0bc = (uint32_t)ws->f3[0][b][c];
                for (uint8_t d = 0U; d < 16U; ++d) {
                    uint32_t fbcd = (uint32_t)ws->f3[b][c][d];
                    phy_vzero_u64(ws->work_row, 255U);
                    phy_q4_add_term(ws->work_row, (fabc << 4U) | fbcd, INT64_C(1));
                    phy_q4_add_term(ws->work_row, tb1, INT64_C(-1));
                    phy_q4_add_term(ws->work_row, tb2, INT64_C(-1));
                    phy_q4_add_term(ws->work_row, (uint32_t)b, INT64_C(1));
                    phy_q4_add_term(ws->work_row, (f0bc << 4U) | fbcd, INT64_C(-1));
                    phy_q4_add_term(ws->work_row, ta1, INT64_C(1));
                    phy_q4_add_term(ws->work_row, ta2, INT64_C(1));
                    phy_q4_add_term(ws->work_row, (uint32_t)a, INT64_C(-1));
                    phy_q4_add_term(ws->work_row, ta3, INT64_C(1));
                    phy_q4_add_term(ws->work_row, ab, INT64_C(-1));
                    if (phy_q4_ech_insert(ws->rref, ws->has_piv, ws->work_row)) {
                        ++rank4;
                    }
                }
            }
        }
    }

    out->tot4 = (uint16_t)(255U - rank4);
    out->irred4 =
        (out->tot4 >= (uint16_t)(15U + out->tot2)) ? (uint16_t)(out->tot4 - 15U - out->tot2) : 0U;

    /* Back-substitute 255-col RREF strictly */
    for (uint32_t idx = 255U; idx > 0U; --idx) {
        uint32_t col = idx - 1U;
        if (ws->has_piv[col] != 0U) {
            const uint64_t *prow = &ws->rref[(size_t)col * 255U];
            for (uint32_t up = 0U; up < col; ++up) {
                if (ws->has_piv[up] != 0U) {
                    uint64_t fac = ws->rref[(size_t)up * 255U + col];
                    if (fac != UINT64_C(0)) {
                        for (uint32_t j = col; j < 255U; ++j) {
                            uint64_t pj = prow[j];
                            if (pj != UINT64_C(0)) {
                                ws->rref[(size_t)up * 255U + j] = phy_m61_sub(
                                    ws->rref[(size_t)up * 255U + j], phy_m61_mul(fac, pj));
                            }
                        }
                    }
                }
            }
        }
    }

    phy_q4_build_trivial_basis(ws);

    bool have_best = false;
    uint16_t best_jnz = 0U;
    uint16_t best_nnz = 256U;
    int32_t best_q4[256];
    phy_vzero_i32(best_q4, 256U);

    for (uint32_t fcol = 0U; fcol < 255U; ++fcol) {
        if (ws->has_piv[fcol] != 0U) {
            continue;
        }
        phy_vzero_u64(ws->null_vec, 255U);
        ws->null_vec[fcol] = UINT64_C(1);
        for (uint32_t pc = 0U; pc < fcol; ++pc) {
            if (ws->has_piv[pc] != 0U) {
                uint64_t v = ws->rref[(size_t)pc * 255U + fcol];
                ws->null_vec[pc] = (v == UINT64_C(0)) ? UINT64_C(0) : (PHY_M61 - v);
            }
        }
        phy_vcopy_u64(ws->work_row, ws->null_vec, 255U);
        if (phy_q4_ech_insert(ws->triv_ech, ws->triv_has, ws->work_row)) {
            out->genuine_q4 += 1U;
            int32_t cand_q4[256];
            uint16_t cand_nnz = 0U;
            if (phy_q4_m61_to_small_int(ws->null_vec, cand_q4, &cand_nnz)) {
                uint8_t ax = 0U;
                uint8_t f1 = 0U;
                uint8_t f2 = 0U;
                uint8_t x0 = 0U;
                uint8_t x1 = 0U;
                uint8_t y0 = 0U;
                uint8_t y1 = 0U;
                if (phy_q4_find_nontriv_witness(cand_q4, &ax, &f1, &f2, &x0, &x1, &y0, &y1)) {
                    uint16_t jnz = 0U;
                    uint8_t wa = 0U;
                    uint8_t wb = 0U;
                    uint8_t wc = 0U;
                    if (phy_q4_verify_continuity_i64(ws->f3, cand_q4, &jnz, &wa, &wb, &wc) == 0) {
                        bool better = !have_best;
                        if (!better) {
                            bool old_act = (best_jnz > 0U);
                            bool new_act = (jnz > 0U);
                            if ((new_act && !old_act) ||
                                (new_act == old_act && cand_nnz < best_nnz) ||
                                (new_act == old_act && cand_nnz == best_nnz && jnz > best_jnz)) {
                                better = true;
                            }
                        }
                        if (better) {
                            have_best = true;
                            best_jnz = jnz;
                            best_nnz = cand_nnz;
                            phy_vcopy_i32(best_q4, cand_q4, 256U);
                        }
                    }
                }
            }
        }
    }

    if (have_best) {
        phy_vcopy_i32(out->q4, best_q4, 256U);
        out->q4_nnz = best_nnz;
        out->continuity4_fail = phy_q4_verify_continuity_i64(
            ws->f3, out->q4, &out->j4_nonzero, &out->j4_wit_a, &out->j4_wit_b, &out->j4_wit_c);
        (void)phy_q4_find_nontriv_witness(out->q4, &out->cross_axis, &out->cross_f1, &out->cross_f2,
                                          &out->cross_x0, &out->cross_x1, &out->cross_y0,
                                          &out->cross_y1);
    }

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    mh = phy_fnv1a64_u64(mh, out->lut_packed);
    mh = phy_fnv1a64_u64(mh, ((uint64_t)out->genuine_q4 << 16U) | (uint64_t)out->j4_nonzero);
    for (uint32_t u = 0U; u < 256U; ++u) {
        mh = phy_fnv1a64_u64(mh, (uint64_t)(uint32_t)out->q4[u]);
    }
    out->q4_merkle = mh;
}

int64_t phy_ring_charge4_64(const int32_t q4[256], uint64_t lo, uint64_t hi) {
    int64_t total = INT64_C(0);
    for (uint32_t m = 0U; m < 32U; ++m) {
        uint32_t mn = (m + 1U) & 31U;
        uint8_t u0 = phy_pack2(phy_get_site64(lo, hi, 2U * m), phy_get_site64(lo, hi, 2U * m + 1U));
        uint8_t u1 =
            phy_pack2(phy_get_site64(lo, hi, 2U * mn), phy_get_site64(lo, hi, 2U * mn + 1U));
        total += (int64_t)q4[((uint32_t)u0 << 4U) | (uint32_t)u1];
    }
    return total;
}

static void phy_store_arena_q4(phy_arena_u *arena, uint32_t idx, const phy_q4_analysis_t *q4) {
    if (arena == NULL || idx >= PHY_COUNT_Q4_META) {
        return;
    }
    uint32_t mslot = PHY_SLOT_Q4_META + idx;
    phy_vzero_u64(arena->cells[mslot].u64, 8U);
    arena->map.q4_meta[idx].lut_packed = q4->lut_packed;
    arena->map.q4_meta[idx].q4_merkle = q4->q4_merkle;
    arena->map.q4_meta[idx].tot2 = q4->tot2;
    arena->map.q4_meta[idx].tot4 = q4->tot4;
    arena->map.q4_meta[idx].irred4 = q4->irred4;
    arena->map.q4_meta[idx].genuine_q4 = q4->genuine_q4;
    arena->map.q4_meta[idx].j4_nonzero = q4->j4_nonzero;
    arena->map.q4_meta[idx].q4_nnz = q4->q4_nnz;
    arena->map.q4_meta[idx].continuity4_fail = q4->continuity4_fail;

    if (idx < 4U) {
        uint32_t base = idx * 16U;
        for (uint32_t r = 0U; r < 16U; ++r) {
            uint32_t dslot = PHY_SLOT_Q4_DENSITIES + base + r;
            phy_vzero_u64(arena->cells[dslot].u64, 8U);
            phy_vcopy_i32(arena->map.q4_densities[base + r].q, &q4->q4[(size_t)r * 16U], 16U);
        }
    }
}

uint32_t phy_hunt_q4(phy_arena_u *arena, phy_q4_ws_t *ws, uint64_t seed, uint64_t max_trials,
                     phy_rule_analysis_t *out_rule, phy_q4_analysis_t *out_q4) {
    if (ws == NULL || out_rule == NULL || out_q4 == NULL) {
        return 0U;
    }
    phy_zero_rule_analysis(out_rule);
    if (arena != NULL) {
        phy_arena_init(arena, seed);
    }

    uint8_t base_swap[16];
    for (uint8_t u = 0U; u < 16U; ++u) {
        base_swap[u] = phy_pack2(phy_lo2(u), phy_hi2(u));
    }

    uint64_t rng = seed ^ UINT64_C(0xD1B54A32D192ED03);
    if (rng == UINT64_C(0)) {
        rng = UINT64_C(0x5A5A5A5AA5A5A5A5);
    }

    uint64_t seen_luts[16];
    phy_vzero_u64(seen_luts, 16U);
    uint32_t found = 0U;
    uint64_t best_sc = UINT64_C(0);
    bool have_best = false;

    uint64_t restarts = max_trials >> 6U;
    if (restarts < 24U) {
        restarts = 24U;
    }
    uint64_t trials = 0U;
    uint32_t rs_mod2 = 0U;
    phy_rule_analysis_t ra;
    phy_q4_analysis_t qa;

    for (uint64_t rs = 0U; rs < restarts && trials < max_trials; ++rs) {
        uint8_t cur_p[16];
        if (rs_mod2 == 0U) {
            phy_vcopy_u8(cur_p, base_swap, 16U);
        } else {
            for (uint8_t u = 0U; u < 16U; ++u) {
                cur_p[u] = u;
            }
        }
        rs_mod2 ^= 1U;
        uint32_t kick = 1U + searu_Bounded_u32(searu_rng_next(&rng), 4U);
        for (uint32_t k = 0U; k < kick; ++k) {
            uint32_t i = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
            uint32_t j = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
            if (i == j) {
                j = (j == 15U) ? 1U : (j + 1U);
            }
            uint8_t tmp = cur_p[i];
            cur_p[i] = cur_p[j];
            cur_p[j] = tmp;
        }
        uint8_t cur_yb = phy_yb_check(cur_p);
        for (uint32_t step = 0U; step < 220U && trials < max_trials; ++step) {
            ++trials;
            if (cur_yb == 0U) {
                phy_analyze_rule(cur_p, seed + rs * UINT64_C(1000) + step, &ra);
                uint32_t req = PHY_RULE_FLAG_BIJECTIVE | PHY_RULE_FLAG_YANGBOX |
                               PHY_RULE_FLAG_NONLINEAR | PHY_RULE_FLAG_NONFACT;
                if ((ra.flags & req) == req) {
                    bool dup = false;
                    for (uint32_t s = 0U; s < found; ++s) {
                        if (seen_luts[s] == ra.lut_packed) {
                            dup = true;
                            break;
                        }
                    }
                    if (!dup) {
                        phy_analyze_q4(cur_p, ws, &qa);
                        if (qa.genuine_q4 > 0U && qa.continuity4_fail == 0 && qa.j4_nonzero > 0U) {
                            if (found < 16U) {
                                seen_luts[found] = ra.lut_packed;
                                phy_store_arena_q4(arena, found, &qa);
                                ++found;
                            }
                            uint64_t sc = ((uint64_t)(ra.rank_ent == 0U ? 1U : 0U) << 48U) |
                                          ((uint64_t)ra.anf_deg << 40U) |
                                          ((uint64_t)qa.genuine_q4 << 24U) |
                                          (uint64_t)qa.j4_nonzero;
                            if (!have_best || sc > best_sc) {
                                have_best = true;
                                best_sc = sc;
                                phy_vcopy_u8((uint8_t *)out_rule, (const uint8_t *)&ra, sizeof(ra));
                                phy_vcopy_u8((uint8_t *)out_q4, (const uint8_t *)&qa, sizeof(qa));
                            }
                        }
                    }
                }
                break;
            }
            uint32_t i = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
            uint32_t j = 1U + searu_Bounded_u32(searu_rng_next(&rng), 15U);
            if (i == j) {
                j = (j == 15U) ? 1U : (j + 1U);
            }
            uint8_t tmp = cur_p[i];
            cur_p[i] = cur_p[j];
            cur_p[j] = tmp;
            uint8_t nxt_yb = phy_yb_check(cur_p);
            if (nxt_yb <= cur_yb || searu_Bounded_u32(searu_rng_next(&rng), 100U) < 6U) {
                cur_yb = nxt_yb;
            } else {
                tmp = cur_p[i];
                cur_p[i] = cur_p[j];
                cur_p[j] = tmp;
            }
        }
    }

    if (have_best && arena != NULL) {
        arena->map.img_hdr.active_rules = found;
        arena->map.img_hdr.merkle_root = out_q4->q4_merkle;
    }
    return found;
}

/* --- |X|=8 Non-Involutive Non-Degenerate Indecomposable Set-Theoretic YBE --- */

static uint8_t phy_hi3(uint8_t u) {
    return (uint8_t)(((uint32_t)u >> 3U) & UINT32_C(7));
}

static uint8_t phy_lo3(uint8_t u) {
    return (uint8_t)((uint32_t)u & UINT32_C(7));
}

static uint8_t phy_pack3(uint8_t l, uint8_t r) {
    return (uint8_t)((((uint32_t)l & UINT32_C(7)) << 3U) | ((uint32_t)r & UINT32_C(7)));
}

uint16_t phy_yb8_check(const uint8_t r[64]) {
    uint16_t yb8_viol = 0U;
    for (uint8_t a = 0U; a < 8U; ++a) {
        for (uint8_t b = 0U; b < 8U; ++b) {
            uint8_t ab = r[phy_pack3(a, b)];
            uint8_t a1 = phy_hi3(ab);
            uint8_t b1 = phy_lo3(ab);
            for (uint8_t c = 0U; c < 8U; ++c) {
                uint8_t bc = r[phy_pack3(b1, c)];
                uint8_t b2 = phy_hi3(bc);
                uint8_t c2 = phy_lo3(bc);
                uint8_t ab2 = r[phy_pack3(a1, b2)];
                uint16_t lhs8 = (uint16_t)(((uint32_t)ab2 << 3U) | (uint32_t)c2);

                uint8_t bc0 = r[phy_pack3(b, c)];
                uint8_t b0 = phy_hi3(bc0);
                uint8_t c0 = phy_lo3(bc0);
                uint8_t ab0 = r[phy_pack3(a, b0)];
                uint8_t a0 = phy_hi3(ab0);
                uint8_t b00 = phy_lo3(ab0);
                uint8_t bc1 = r[phy_pack3(b00, c0)];
                uint16_t rhs8 = (uint16_t)(((uint32_t)a0 << 6U) | (uint32_t)bc1);

                if (lhs8 != rhs8) {
                    ++yb8_viol;
                }
            }
        }
    }
    return yb8_viol;
}

static uint8_t phy_anf_degree6(const uint8_t r[64], uint8_t *out_bit, uint8_t *out_a,
                               uint8_t *out_b) {
    uint8_t max_deg6 = 0U;
    uint8_t best_bit = 0U;
    uint8_t best_a = 0U;
    uint8_t best_b = 0U;
    bool found_nl = false;
    for (uint8_t bit = 0U; bit < 6U; ++bit) {
        uint8_t anf[64];
        for (uint32_t x = 0U; x < 64U; ++x) {
            anf[x] = (uint8_t)(((uint32_t)r[x] >> (uint32_t)bit) & UINT32_C(1));
        }
        for (uint32_t var = 0U; var < 6U; ++var) {
            for (uint32_t mask = 0U; mask < 64U; ++mask) {
                if ((mask & (UINT32_C(1) << var)) != 0U) {
                    anf[mask] ^= anf[mask ^ (UINT32_C(1) << var)];
                }
            }
        }
        for (uint32_t mask = 0U; mask < 64U; ++mask) {
            if (anf[mask] != 0U) {
                uint8_t d = (uint8_t)searu_popcnt64((uint64_t)mask);
                if (d > max_deg6) {
                    max_deg6 = d;
                }
            }
        }
        if (!found_nl) {
            for (uint8_t a = 1U; a < 64U && !found_nl; ++a) {
                for (uint8_t b = 1U; b < 64U && !found_nl; ++b) {
                    uint32_t xab = (uint32_t)r[a ^ b] >> (uint32_t)bit;
                    uint32_t xa = (uint32_t)r[a] >> (uint32_t)bit;
                    uint32_t xb = (uint32_t)r[b] >> (uint32_t)bit;
                    uint32_t x0 = (uint32_t)r[0] >> (uint32_t)bit;
                    if (((xab ^ xa ^ xb ^ x0) & UINT32_C(1)) != 0U) {
                        found_nl = true;
                        best_bit = bit;
                        best_a = a;
                        best_b = b;
                    }
                }
            }
        }
    }
    *out_bit = best_bit;
    *out_a = best_a;
    *out_b = best_b;
    return max_deg6;
}

static uint32_t phy_perm64_order(const uint8_t r[64]) {
    uint8_t vis[64];
    phy_vzero_u8(vis, 64U);
    uint64_t lcm = UINT64_C(1);
    for (uint8_t i = 0U; i < 64U; ++i) {
        if (vis[i] == 0U) {
            uint64_t clen = UINT64_C(0);
            uint8_t cur = i;
            while (vis[cur] == 0U && clen < 65U) {
                vis[cur] = 1U;
                cur = (uint8_t)(r[cur] & 63U);
                ++clen;
            }
            if (clen > UINT64_C(0)) {
                uint64_t g = phy_stein_gcd64(lcm, clen);
                lcm = phy_div_u64(lcm, g, NULL) * clen;
            }
        }
    }
    return (uint32_t)lcm;
}

void phy_analyze_yb8(const uint8_t r[64], phy_yb8_analysis_t *out) {
    if (r == NULL || out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    uint64_t seen = UINT64_C(0);
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint8_t v = (uint8_t)(r[u] & 63U);
        out->r[u] = v;
        out->rinv[v] = u;
        seen |= (UINT64_C(1) << v);
    }
    out->is_bijective = (seen == UINT64_C(0xFFFFFFFFFFFFFFFF)) ? UINT8_C(1) : UINT8_C(0);
    out->yb8_violations = phy_yb8_check(out->r);
    out->anf_deg = phy_anf_degree6(out->r, &out->nl_bit, &out->nl_a, &out->nl_b);
    out->order = (out->is_bijective != 0U) ? phy_perm64_order(out->r) : UINT32_C(0);

    /* Non-involutivity check */
    out->is_involutive = 1U;
    for (uint8_t u = 0U; u < 64U; ++u) {
        if (out->r[out->r[u]] != u) {
            out->is_involutive = 0U;
            out->noninvol_u = u;
            break;
        }
    }

    /* Non-degeneracy: sigma_a(b) and tau_b(a) bijections */
    bool nondeg = true;
    for (uint8_t a = 0U; a < 8U && nondeg; ++a) {
        uint32_t mask = UINT32_C(0);
        for (uint8_t b = 0U; b < 8U; ++b) {
            uint8_t s = phy_hi3(out->r[phy_pack3(a, b)]);
            mask |= (UINT32_C(1) << s);
            out->sig_inv[a][s] = b;
        }
        if (mask != UINT32_C(0xFF)) {
            nondeg = false;
        }
    }
    for (uint8_t b = 0U; b < 8U && nondeg; ++b) {
        uint32_t mask = UINT32_C(0);
        for (uint8_t a = 0U; a < 8U; ++a) {
            uint8_t t = phy_lo3(out->r[phy_pack3(a, b)]);
            mask |= (UINT32_C(1) << t);
            out->tau_inv[b][t] = a;
        }
        if (mask != UINT32_C(0xFF)) {
            nondeg = false;
        }
    }
    out->is_nondeg = nondeg ? UINT8_C(1) : UINT8_C(0);

    /* Indecomposability via 7-step BFS reachability from 0 */
    uint8_t reached = 1U;
    for (uint32_t k = 0U; k < 7U; ++k) {
        uint8_t nxt = reached;
        for (uint8_t u = 0U; u < 8U; ++u) {
            if ((((uint32_t)reached >> u) & 1U) == 0U) {
                continue;
            }
            for (uint8_t x = 0U; x < 8U; ++x) {
                uint8_t s1 = phy_hi3(out->r[phy_pack3(x, u)]);
                uint8_t t1 = phy_lo3(out->r[phy_pack3(u, x)]);
                nxt |= (uint8_t)(UINT32_C(1) << s1);
                nxt |= (uint8_t)(UINT32_C(1) << t1);
                for (uint8_t v = 0U; v < 8U; ++v) {
                    if (phy_hi3(out->r[phy_pack3(x, v)]) == u ||
                        phy_lo3(out->r[phy_pack3(v, x)]) == u) {
                        nxt |= (uint8_t)(UINT32_C(1) << v);
                    }
                }
            }
        }
        reached = nxt;
    }
    out->is_indecomposable = (reached == UINT8_C(0xFF)) ? UINT8_C(1) : UINT8_C(0);

    /* Non-Lyubashenko cross witness */
    bool lyub = true;
    for (uint8_t a = 1U; a < 8U && lyub; ++a) {
        for (uint8_t b = 0U; b < 8U; ++b) {
            if (phy_hi3(out->r[phy_pack3(a, b)]) != phy_hi3(out->r[phy_pack3(0U, b)])) {
                lyub = false;
                out->nly_axis = 0U;
                out->nly_a = 0U;
                out->nly_b0 = a;
                out->nly_b1 = b;
                break;
            }
        }
    }
    for (uint8_t a = 0U; a < 8U && lyub; ++a) {
        for (uint8_t b = 1U; b < 8U; ++b) {
            if (phy_lo3(out->r[phy_pack3(a, b)]) != phy_lo3(out->r[phy_pack3(a, 0U)])) {
                lyub = false;
                out->nly_axis = 1U;
                out->nly_a = a;
                out->nly_b0 = 0U;
                out->nly_b1 = b;
                break;
            }
        }
    }
    out->is_lyubashenko = lyub ? UINT8_C(1) : UINT8_C(0);

    uint32_t fl = UINT32_C(0);
    if (out->is_bijective != 0U) {
        fl |= PHY_YB8_FLAG_BIJECTIVE;
    }
    if (out->yb8_violations == 0U) {
        fl |= PHY_YB8_FLAG_YBE;
    }
    if (out->is_nondeg != 0U) {
        fl |= PHY_YB8_FLAG_NONDEG;
    }
    if (out->is_involutive == 0U) {
        fl |= PHY_YB8_FLAG_NONINVOL;
    }
    if (out->is_indecomposable != 0U) {
        fl |= PHY_YB8_FLAG_INDECOMP;
    }
    if (out->is_lyubashenko == 0U) {
        fl |= PHY_YB8_FLAG_NONLYUB;
    }
    if (out->anf_deg >= 2U) {
        fl |= PHY_YB8_FLAG_NONLINEAR;
    }
    out->flags = fl;

    uint64_t yb8_mh = UINT64_C(0xCBF29CE484222325);
    for (uint32_t w = 0U; w < 8U; ++w) {
        uint64_t blk = UINT64_C(0);
        for (uint32_t k = 0U; k < 8U; ++k) {
            blk |= ((uint64_t)out->r[w * 8U + k]) << (k * 8U);
        }
        yb8_mh = phy_fnv1a64_u64(yb8_mh, blk);
    }
    yb8_mh = phy_fnv1a64_u64(yb8_mh, ((uint64_t)fl << 32U) | (uint64_t)out->order);
    out->merkle = yb8_mh;
}

static void phy_store_arena_yb8(phy_arena_u *arena, uint32_t idx, const phy_yb8_analysis_t *yb8) {
    if (arena == NULL || idx >= PHY_COUNT_YB8_RULES) {
        return;
    }
    uint32_t rslot = PHY_SLOT_YB8_RULES + idx;
    phy_vzero_u64(arena->cells[rslot].u64, 8U);
    phy_vcopy_u8(arena->map.yb8_rules[idx].r, yb8->r, 64U);

    uint32_t mslot = PHY_SLOT_YB8_META + idx;
    phy_vzero_u64(arena->cells[mslot].u64, 8U);
    arena->map.yb8_meta[idx].merkle = yb8->merkle;
    arena->map.yb8_meta[idx].flags = yb8->flags;
    arena->map.yb8_meta[idx].order = yb8->order;
    arena->map.yb8_meta[idx].yb8_violations = yb8->yb8_violations;
    arena->map.yb8_meta[idx].anf_deg = yb8->anf_deg;
    arena->map.yb8_meta[idx].is_bijective = yb8->is_bijective;
    arena->map.yb8_meta[idx].is_nondeg = yb8->is_nondeg;
    arena->map.yb8_meta[idx].is_involutive = yb8->is_involutive;
    arena->map.yb8_meta[idx].is_indecomposable = yb8->is_indecomposable;
    arena->map.yb8_meta[idx].is_lyubashenko = yb8->is_lyubashenko;
    arena->map.yb8_meta[idx].noninvol_u = yb8->noninvol_u;
}

uint32_t phy_hunt_yb8(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                      phy_yb8_analysis_t *out_best) {
    if (out_best == NULL) {
        return 0U;
    }
    phy_vzero_u8((uint8_t *)out_best, sizeof(*out_best));
    if (arena != NULL) {
        phy_arena_init(arena, seed);
    }

    uint64_t rng = seed ^ UINT64_C(0x94D049BB133111EB);
    if (rng == UINT64_C(0)) {
        rng = UINT64_C(0x13572468ACE09753);
    }

    uint64_t seen_merkles[16];
    phy_vzero_u64(seen_merkles, 16U);
    uint32_t found = 0U;
    uint64_t best_sc = UINT64_C(0);
    bool have_best = false;

    uint32_t req = PHY_YB8_FLAG_BIJECTIVE | PHY_YB8_FLAG_YBE | PHY_YB8_FLAG_NONDEG |
                   PHY_YB8_FLAG_NONINVOL | PHY_YB8_FLAG_INDECOMP | PHY_YB8_FLAG_NONLYUB |
                   PHY_YB8_FLAG_NONLINEAR;

    phy_yb8_analysis_t cand;
    for (uint64_t t = 0U; t < max_trials; ++t) {
        uint8_t s = (uint8_t)(1U + 2U * searu_Bounded_u32(searu_rng_next(&rng), 4U));
        uint8_t tt = (uint8_t)(1U + 2U * searu_Bounded_u32(searu_rng_next(&rng), 4U));
        uint8_t v = (uint8_t)searu_Bounded_u32(searu_rng_next(&rng), 8U);
        uint8_t w = (uint8_t)searu_Bounded_u32(searu_rng_next(&rng), 8U);
        uint8_t c1 = (uint8_t)searu_Bounded_u32(searu_rng_next(&rng), 8U);
        uint8_t c2 = (uint8_t)searu_Bounded_u32(searu_rng_next(&rng), 8U);

        uint8_t pi[8];
        uint8_t pi_inv[8];
        for (uint8_t x = 0U; x < 8U; ++x) {
            pi[x] = x;
        }
        if ((t & UINT64_C(3)) != UINT64_C(0)) {
            uint32_t nsw = 1U + searu_Bounded_u32(searu_rng_next(&rng), 2U);
            for (uint32_t k = 0U; k < nsw; ++k) {
                uint32_t i = searu_Bounded_u32(searu_rng_next(&rng), 8U);
                uint32_t j = searu_Bounded_u32(searu_rng_next(&rng), 8U);
                uint8_t tmp = pi[i];
                pi[i] = pi[j];
                pi[j] = tmp;
            }
        }
        for (uint8_t x = 0U; x < 8U; ++x) {
            pi_inv[pi[x]] = x;
        }

        uint8_t r[64];
        for (uint8_t a = 0U; a < 8U; ++a) {
            uint8_t pa = pi[a];
            for (uint8_t b = 0U; b < 8U; ++b) {
                uint8_t pb = pi[b];
                uint8_t oa = (uint8_t)(((uint32_t)s * pb + (uint32_t)v * pa + c1) & UINT32_C(7));
                uint8_t ob = (uint8_t)(((uint32_t)tt * pa + (uint32_t)w * pb + c2) & UINT32_C(7));
                r[phy_pack3(a, b)] = phy_pack3(pi_inv[oa], pi_inv[ob]);
            }
        }

        if (phy_yb8_check(r) != 0U) {
            continue;
        }
        phy_analyze_yb8(r, &cand);
        if ((cand.flags & req) == req && cand.order >= 4U) {
            bool dup = false;
            for (uint32_t k = 0U; k < found; ++k) {
                if (seen_merkles[k] == cand.merkle) {
                    dup = true;
                    break;
                }
            }
            if (!dup && found < PHY_COUNT_YB8_RULES) {
                seen_merkles[found] = cand.merkle;
                phy_store_arena_yb8(arena, found, &cand);
                ++found;
            }
            /* Favor high-order & clean deg=2..3 rules compatible with low-AND .fir synthesis */
            uint64_t low_mc_bonus = (cand.anf_deg <= 3U) ? UINT64_C(1) : UINT64_C(0);
            uint64_t sc = (low_mc_bonus << 40U) | ((uint64_t)(cand.order >= 8U ? 1U : 0U) << 36U) |
                          ((uint64_t)cand.anf_deg << 24U) | (uint64_t)cand.order;
            if (!have_best || sc > best_sc) {
                have_best = true;
                best_sc = sc;
                phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
            }
        }
    }

    if (have_best && arena != NULL) {
        arena->map.img_hdr.active_rules = found;
        arena->map.img_hdr.merkle_root = out_best->merkle;
    }
    return found;
}

/* ============================================================================
 * PHY v3.0: Dual-Unitary, AME(4,8), Fracton HSF, Biquandle8 Knots, Tetrahedron
 * ========================================================================== */

static uint8_t phy_gf8_mul(uint8_t a, uint8_t b) {
    uint32_t aa = (uint32_t)(a & 7U);
    uint32_t bb = (uint32_t)(b & 7U);
    uint32_t res = UINT32_C(0);
    for (uint32_t k = 0U; k < 3U; ++k) {
        if ((bb & UINT32_C(1)) != UINT32_C(0)) {
            res ^= aa;
        }
        uint32_t hi = aa & UINT32_C(4);
        aa = (aa << 1U) & UINT32_C(7);
        if (hi != UINT32_C(0)) {
            aa ^= UINT32_C(3);
        }
        bb >>= 1U;
    }
    return (uint8_t)(res & UINT32_C(7));
}

static uint64_t phy_pack_mat4(const uint8_t m[4][4]) {
    uint64_t acc = UINT64_C(0);
    for (uint32_t r = 0U; r < 4U; ++r) {
        for (uint32_t c = 0U; c < 4U; ++c) {
            uint32_t idx = (r << 2U) | c;
            acc |= ((uint64_t)(m[r][c] & 0xFU)) << (idx * 4U);
        }
    }
    return acc;
}

static bool phy_mat4_pow4_all_pos(const uint8_t m[4][4]) {
    uint32_t m2[4][4];
    for (uint32_t i = 0U; i < 4U; ++i) {
        for (uint32_t j = 0U; j < 4U; ++j) {
            uint32_t acc = UINT32_C(0);
            for (uint32_t k = 0U; k < 4U; ++k) {
                acc += (uint32_t)m[i][k] * (uint32_t)m[k][j];
            }
            m2[i][j] = acc;
        }
    }
    for (uint32_t i = 0U; i < 4U; ++i) {
        for (uint32_t j = 0U; j < 4U; ++j) {
            uint32_t acc = UINT32_C(0);
            for (uint32_t k = 0U; k < 4U; ++k) {
                acc += m2[i][k] * m2[k][j];
            }
            if (acc == UINT32_C(0)) {
                return false;
            }
        }
    }
    return true;
}

void phy_analyze_dual4(const uint8_t p[16], phy_dual4_analysis_t *out) {
    if (p == NULL || out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));

    phy_rule_analysis_t ra;
    phy_analyze_rule(p, UINT64_C(0), &ra);
    phy_vcopy_u8(out->p, ra.p, 16U);
    phy_vcopy_u8(out->pinv, ra.pinv, 16U);
    out->lut_packed = ra.lut_packed;
    out->anf_deg = ra.anf_deg;
    out->is_bijective = ra.is_bijective;
    out->yb_violations = ra.yb_violations;
    out->rank_ent = ra.rank_ent;
    out->nl_bit = ra.nl_bit;
    out->nl_a = ra.nl_a;
    out->nl_b = ra.nl_b;

    bool nondeg = (out->is_bijective != 0U);
    for (uint8_t a = 0U; a < 4U && nondeg; ++a) {
        uint32_t mask = UINT32_C(0);
        for (uint8_t b = 0U; b < 4U; ++b) {
            uint8_t c = phy_hi2(out->p[phy_pack2(a, b)]);
            mask |= (UINT32_C(1) << c);
            out->sig_inv[a][c] = b;
        }
        if (mask != UINT32_C(0xF)) {
            nondeg = false;
        }
    }
    for (uint8_t b = 0U; b < 4U && nondeg; ++b) {
        uint32_t mask = UINT32_C(0);
        for (uint8_t a = 0U; a < 4U; ++a) {
            uint8_t d = phy_lo2(out->p[phy_pack2(a, b)]);
            mask |= (UINT32_C(1) << d);
            out->tau_inv[b][d] = a;
        }
        if (mask != UINT32_C(0xF)) {
            nondeg = false;
        }
    }
    out->is_nondeg = nondeg ? UINT8_C(1) : UINT8_C(0);
    out->is_dual_unitary = (out->is_bijective != 0U && nondeg) ? UINT8_C(1) : UINT8_C(0);

    bool non_bern = false;
    for (uint8_t a = 0U; a < 4U; ++a) {
        for (uint8_t b = 0U; b < 4U; ++b) {
            uint8_t cd = out->p[phy_pack2(a, b)];
            uint8_t c = phy_hi2(cd);
            uint8_t d = phy_lo2(cd);
            out->m_plus[a][d] += 1U;
            out->m_minus[b][c] += 1U;
        }
    }
    for (uint32_t r = 0U; r < 4U; ++r) {
        for (uint32_t c = 0U; c < 4U; ++c) {
            if (out->m_plus[r][c] != 1U || out->m_minus[r][c] != 1U) {
                non_bern = true;
            }
        }
    }
    out->mplus_packed = phy_pack_mat4(out->m_plus);
    out->mminus_packed = phy_pack_mat4(out->m_minus);
    out->is_non_bernoulli = non_bern ? UINT8_C(1) : UINT8_C(0);
    out->is_mixing = (out->is_dual_unitary != 0U && phy_mat4_pow4_all_pos(out->m_plus) &&
                      phy_mat4_pow4_all_pos(out->m_minus))
                         ? UINT8_C(1)
                         : UINT8_C(0);

    uint32_t fl = UINT32_C(0);
    if (out->is_bijective != 0U) {
        fl |= PHY_DUAL4_FLAG_BIJECTIVE;
    }
    if (out->is_nondeg != 0U) {
        fl |= PHY_DUAL4_FLAG_NONDEG;
    }
    if (out->is_dual_unitary != 0U) {
        fl |= PHY_DUAL4_FLAG_DUAL;
    }
    if (out->is_mixing != 0U) {
        fl |= PHY_DUAL4_FLAG_MIXING;
    }
    if (out->is_non_bernoulli != 0U) {
        fl |= PHY_DUAL4_FLAG_NONBERN;
    }
    if (out->anf_deg >= 2U) {
        fl |= PHY_DUAL4_FLAG_NONLINEAR;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    mh = phy_fnv1a64_u64(mh, out->lut_packed);
    mh = phy_fnv1a64_u64(mh, out->mplus_packed ^ (out->mminus_packed << 1U));
    mh = phy_fnv1a64_u64(mh, (uint64_t)fl);
    out->merkle = mh;
}

static void phy_store_arena_dual4(phy_arena_u *arena, uint32_t idx,
                                  const phy_dual4_analysis_t *d4) {
    if (arena == NULL || idx >= PHY_COUNT_DUAL4_META) {
        return;
    }
    uint32_t slot = PHY_SLOT_DUAL4_META + idx;
    phy_vzero_u64(arena->cells[slot].u64, 8U);
    arena->map.dual4_meta[idx].lut_packed = d4->lut_packed;
    arena->map.dual4_meta[idx].merkle = d4->merkle;
    arena->map.dual4_meta[idx].mplus_packed = d4->mplus_packed;
    arena->map.dual4_meta[idx].mminus_packed = d4->mminus_packed;
    arena->map.dual4_meta[idx].flags = d4->flags;
    arena->map.dual4_meta[idx].anf_deg = d4->anf_deg;
    arena->map.dual4_meta[idx].is_dual_unitary = d4->is_dual_unitary;
    arena->map.dual4_meta[idx].is_mixing = d4->is_mixing;
    arena->map.dual4_meta[idx].is_non_bernoulli = d4->is_non_bernoulli;
    arena->map.dual4_meta[idx].yb_violations = d4->yb_violations;
    arena->map.dual4_meta[idx].rank_ent = d4->rank_ent;
}

static void phy_gen_s4_perms(uint8_t s4[24][4], uint8_t fix0_idx[6]) {
    uint32_t n = 0U;
    uint32_t f = 0U;
    for (uint8_t a = 0U; a < 4U; ++a) {
        for (uint8_t b = 0U; b < 4U; ++b) {
            if (b == a) {
                continue;
            }
            for (uint8_t c = 0U; c < 4U; ++c) {
                if (c == a || c == b) {
                    continue;
                }
                uint8_t d = (uint8_t)(6U - a - b - c);
                s4[n][0] = a;
                s4[n][1] = b;
                s4[n][2] = c;
                s4[n][3] = d;
                if (a == 0U && f < 6U) {
                    fix0_idx[f++] = (uint8_t)n;
                }
                ++n;
            }
        }
    }
}

uint32_t phy_hunt_dual4(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                        phy_dual4_analysis_t *out_best) {
    if (out_best == NULL) {
        return 0U;
    }
    phy_vzero_u8((uint8_t *)out_best, sizeof(*out_best));
    if (arena != NULL) {
        phy_arena_init(arena, seed);
    }

    uint8_t s4[24][4];
    uint8_t fix0_idx[6];
    phy_gen_s4_perms(s4, fix0_idx);

    uint64_t rng = seed ^ UINT64_C(0xA4093822299F31D0);
    if (rng == UINT64_C(0)) {
        rng = UINT64_C(0x13198A2E03707344);
    }

    uint64_t seen[PHY_COUNT_DUAL4_META];
    phy_vzero_u64(seen, PHY_COUNT_DUAL4_META);
    uint32_t found = 0U;
    uint64_t best_sc = UINT64_C(0);
    bool have_best = false;

    uint32_t req = PHY_DUAL4_FLAG_BIJECTIVE | PHY_DUAL4_FLAG_NONDEG | PHY_DUAL4_FLAG_DUAL |
                   PHY_DUAL4_FLAG_MIXING | PHY_DUAL4_FLAG_NONBERN | PHY_DUAL4_FLAG_NONLINEAR;

    phy_dual4_analysis_t cand;
    for (uint64_t t = 0U; t < max_trials; ++t) {
        uint8_t c_row[4];
        c_row[0] = fix0_idx[searu_Bounded_u32(searu_rng_next(&rng), 6U)];
        c_row[1] = (uint8_t)searu_Bounded_u32(searu_rng_next(&rng), 24U);
        c_row[2] = (uint8_t)searu_Bounded_u32(searu_rng_next(&rng), 24U);
        c_row[3] = (uint8_t)searu_Bounded_u32(searu_rng_next(&rng), 24U);

        uint8_t r_off[4];
        uint8_t try_k[4];
        uint8_t d_col[4];
        uint16_t mask_stk[5];
        for (uint32_t b = 0U; b < 4U; ++b) {
            r_off[b] = (uint8_t)searu_Bounded_u32(searu_rng_next(&rng), 24U);
            try_k[b] = 0U;
            d_col[b] = 0U;
        }
        mask_stk[0] = 0U;

        /* Iterative non-recursive backtracking over the 4 column permutations */
        int32_t col = 0;
        bool solved = false;
        while (col >= 0) {
            if (try_k[col] >= 24U) {
                try_k[col] = 0U;
                --col;
                continue;
            }
            uint8_t pid = (uint8_t)(r_off[col] + try_k[col]);
            if (pid >= 24U) {
                pid = (uint8_t)(pid - 24U);
            }
            ++try_k[col];
            if (col == 0 && s4[pid][0] != 0U) {
                continue;
            }
            uint16_t cur_mask = mask_stk[col];
            bool ok = true;
            for (uint8_t a = 0U; a < 4U; ++a) {
                uint8_t pair = phy_pack2(s4[c_row[a]][col], s4[pid][a]);
                uint16_t bit = (uint16_t)(UINT32_C(1) << pair);
                if ((cur_mask & bit) != 0U) {
                    ok = false;
                    break;
                }
                cur_mask |= bit;
            }
            if (ok) {
                d_col[col] = pid;
                mask_stk[col + 1] = cur_mask;
                if (col == 3) {
                    solved = true;
                    break;
                }
                ++col;
                try_k[col] = 0U;
            }
        }
        if (!solved) {
            continue;
        }

        uint8_t p_cand[16];
        for (uint8_t a = 0U; a < 4U; ++a) {
            for (uint8_t b = 0U; b < 4U; ++b) {
                p_cand[phy_pack2(a, b)] = phy_pack2(s4[c_row[a]][b], s4[d_col[b]][a]);
            }
        }

        phy_analyze_dual4(p_cand, &cand);
        if ((cand.flags & req) == req && cand.anf_deg >= 3U && cand.rank_ent == 0U) {
            bool dup = false;
            for (uint32_t k = 0U; k < found; ++k) {
                if (seen[k] == cand.merkle) {
                    dup = true;
                    break;
                }
            }
            if (!dup && found < PHY_COUNT_DUAL4_META) {
                seen[found] = cand.merkle;
                phy_store_arena_dual4(arena, found, &cand);
                ++found;
            }
            uint32_t zeros_p = 0U;
            for (uint32_t r = 0U; r < 4U; ++r) {
                for (uint32_t c = 0U; c < 4U; ++c) {
                    if (cand.m_plus[r][c] == 0U) {
                        ++zeros_p;
                    }
                    if (cand.m_minus[r][c] == 0U) {
                        ++zeros_p;
                    }
                }
            }
            uint64_t sc = ((uint64_t)cand.anf_deg << 32U) | ((uint64_t)zeros_p << 16U) |
                          ((uint64_t)cand.yb_violations);
            if (!have_best || sc > best_sc) {
                have_best = true;
                best_sc = sc;
                phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
            }
        }
    }

    if (have_best && arena != NULL) {
        arena->map.img_hdr.active_rules = found;
        arena->map.img_hdr.merkle_root = out_best->merkle;
    }
    return found;
}

static void phy_store_arena_v3_rule6(phy_arena_u *arena, uint32_t idx, const uint8_t r[64]) {
    if (arena == NULL || idx >= PHY_COUNT_V3_RULES6) {
        return;
    }
    uint32_t slot = PHY_SLOT_V3_RULES6 + idx;
    phy_vzero_u64(arena->cells[slot].u64, 8U);
    phy_vcopy_u8(arena->map.v3_rules6[idx].r, r, 64U);
}

void phy_analyze_ame8(const uint8_t r[64], phy_ame8_analysis_t *out) {
    if (r == NULL || out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    uint64_t seen = UINT64_C(0);
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint8_t v = (uint8_t)(r[u] & 63U);
        out->r[u] = v;
        out->rinv[v] = u;
        seen |= (UINT64_C(1) << v);
    }
    out->is_bijective = (seen == UINT64_C(0xFFFFFFFFFFFFFFFF)) ? UINT8_C(1) : UINT8_C(0);
    out->anf_deg = phy_anf_degree6(out->r, &out->nl_bit, &out->nl_a, &out->nl_b);
    out->order = (out->is_bijective != 0U) ? phy_perm64_order(out->r) : UINT32_C(0);

    bool nondeg = (out->is_bijective != 0U);
    for (uint8_t a = 0U; a < 8U && nondeg; ++a) {
        uint32_t mask = UINT32_C(0);
        for (uint8_t b = 0U; b < 8U; ++b) {
            uint8_t c = phy_hi3(out->r[phy_pack3(a, b)]);
            mask |= (UINT32_C(1) << c);
            out->sig_inv[a][c] = b;
        }
        if (mask != UINT32_C(0xFF)) {
            nondeg = false;
        }
    }
    for (uint8_t b = 0U; b < 8U && nondeg; ++b) {
        uint32_t mask = UINT32_C(0);
        for (uint8_t a = 0U; a < 8U; ++a) {
            uint8_t d = phy_lo3(out->r[phy_pack3(a, b)]);
            mask |= (UINT32_C(1) << d);
            out->tau_inv[b][d] = a;
        }
        if (mask != UINT32_C(0xFF)) {
            nondeg = false;
        }
    }
    out->is_nondeg = nondeg ? UINT8_C(1) : UINT8_C(0);

    bool cross = (out->is_bijective != 0U);
    for (uint8_t a = 0U; a < 8U && cross; ++a) {
        uint32_t mask = UINT32_C(0);
        for (uint8_t b = 0U; b < 8U; ++b) {
            uint8_t d = phy_lo3(out->r[phy_pack3(a, b)]);
            mask |= (UINT32_C(1) << d);
            out->cross_l_inv[a][d] = b;
        }
        if (mask != UINT32_C(0xFF)) {
            cross = false;
        }
    }
    for (uint8_t b = 0U; b < 8U && cross; ++b) {
        uint32_t mask = UINT32_C(0);
        for (uint8_t a = 0U; a < 8U; ++a) {
            uint8_t c = phy_hi3(out->r[phy_pack3(a, b)]);
            mask |= (UINT32_C(1) << c);
            out->cross_r_inv[b][c] = a;
        }
        if (mask != UINT32_C(0xFF)) {
            cross = false;
        }
    }
    out->is_crossbij = cross ? UINT8_C(1) : UINT8_C(0);
    out->is_2unitary = (out->is_bijective != 0U && nondeg && cross) ? UINT8_C(1) : UINT8_C(0);

    uint32_t fl = UINT32_C(0);
    if (out->is_bijective != 0U) {
        fl |= PHY_AME8_FLAG_BIJECTIVE;
    }
    if (out->is_nondeg != 0U) {
        fl |= PHY_AME8_FLAG_NONDEG;
    }
    if (out->is_crossbij != 0U) {
        fl |= PHY_AME8_FLAG_CROSSBIJ;
    }
    if (out->is_2unitary != 0U) {
        fl |= PHY_AME8_FLAG_2UNITARY;
    }
    if (out->anf_deg >= 2U) {
        fl |= PHY_AME8_FLAG_NONLINEAR;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    for (uint32_t w = 0U; w < 8U; ++w) {
        uint64_t blk = UINT64_C(0);
        for (uint32_t k = 0U; k < 8U; ++k) {
            blk |= ((uint64_t)out->r[w * 8U + k]) << (k * 8U);
        }
        mh = phy_fnv1a64_u64(mh, blk);
    }
    mh = phy_fnv1a64_u64(mh, ((uint64_t)fl << 32U) | (uint64_t)out->order);
    out->merkle = mh;
}

static void phy_store_arena_ame8(phy_arena_u *arena, uint32_t idx, const phy_ame8_analysis_t *a8) {
    if (arena == NULL || idx >= PHY_COUNT_AME8_META) {
        return;
    }
    uint32_t slot = PHY_SLOT_AME8_META + idx;
    phy_vzero_u64(arena->cells[slot].u64, 8U);
    arena->map.ame8_meta[idx].merkle = a8->merkle;
    arena->map.ame8_meta[idx].flags = a8->flags;
    arena->map.ame8_meta[idx].order = a8->order;
    arena->map.ame8_meta[idx].anf_deg = a8->anf_deg;
    arena->map.ame8_meta[idx].is_bijective = a8->is_bijective;
    arena->map.ame8_meta[idx].is_nondeg = a8->is_nondeg;
    arena->map.ame8_meta[idx].is_crossbij = a8->is_crossbij;
    arena->map.ame8_meta[idx].is_2unitary = a8->is_2unitary;
    arena->map.ame8_meta[idx].nl_bit = a8->nl_bit;
    arena->map.ame8_meta[idx].nl_a = a8->nl_a;
    arena->map.ame8_meta[idx].nl_b = a8->nl_b;
    phy_store_arena_v3_rule6(arena, idx, a8->r);
}

uint32_t phy_hunt_ame8(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                       phy_ame8_analysis_t *out_best) {
    if (out_best == NULL) {
        return 0U;
    }
    phy_vzero_u8((uint8_t *)out_best, sizeof(*out_best));
    if (arena != NULL) {
        phy_arena_init(arena, seed);
    }

    uint64_t rng = seed ^ UINT64_C(0x517CC1B727220A95);
    if (rng == UINT64_C(0)) {
        rng = UINT64_C(0x243F6A8885A308D3);
    }

    uint64_t seen[PHY_COUNT_AME8_META];
    phy_vzero_u64(seen, PHY_COUNT_AME8_META);
    uint32_t found = 0U;
    uint64_t best_sc = UINT64_C(0);
    bool have_best = false;

    uint32_t req = PHY_AME8_FLAG_BIJECTIVE | PHY_AME8_FLAG_NONDEG | PHY_AME8_FLAG_CROSSBIJ |
                   PHY_AME8_FLAG_2UNITARY | PHY_AME8_FLAG_NONLINEAR;

    phy_ame8_analysis_t cand;
    for (uint64_t t = 0U; t < max_trials; ++t) {
        uint8_t s1 = (uint8_t)(1U + searu_Bounded_u32(searu_rng_next(&rng), 7U));
        uint8_t s2 = (uint8_t)(1U + searu_Bounded_u32(searu_rng_next(&rng), 7U));
        uint8_t t1 = (uint8_t)(1U + searu_Bounded_u32(searu_rng_next(&rng), 7U));
        uint8_t t2 = (uint8_t)(1U + searu_Bounded_u32(searu_rng_next(&rng), 7U));
        if (phy_gf8_mul(s1, t2) == phy_gf8_mul(s2, t1)) {
            continue;
        }

        uint8_t pic[8];
        uint8_t pid[8];
        for (uint8_t x = 0U; x < 8U; ++x) {
            pic[x] = x;
            pid[x] = x;
        }
        if (t == UINT64_C(0)) {
            pic[6] = 7U;
            pic[7] = 6U;
            pid[6] = 7U;
            pid[7] = 6U;
        } else {
            uint32_t nc = 1U + searu_Bounded_u32(searu_rng_next(&rng), 2U);
            for (uint32_t k = 0U; k < nc; ++k) {
                uint32_t i = 1U + searu_Bounded_u32(searu_rng_next(&rng), 7U);
                uint32_t j = 1U + searu_Bounded_u32(searu_rng_next(&rng), 7U);
                if (i == j) {
                    j = (j == 7U) ? 1U : (j + 1U);
                }
                uint8_t tmp = pic[i];
                pic[i] = pic[j];
                pic[j] = tmp;
            }
            if ((t & UINT64_C(1)) != UINT64_C(0)) {
                uint32_t i = 1U + searu_Bounded_u32(searu_rng_next(&rng), 7U);
                uint32_t j = 1U + searu_Bounded_u32(searu_rng_next(&rng), 7U);
                if (i == j) {
                    j = (j == 7U) ? 1U : (j + 1U);
                }
                uint8_t tmp = pid[i];
                pid[i] = pid[j];
                pid[j] = tmp;
            }
        }

        uint8_t r[64];
        for (uint8_t a = 0U; a < 8U; ++a) {
            for (uint8_t b = 0U; b < 8U; ++b) {
                uint8_t c = (uint8_t)(phy_gf8_mul(s1, a) ^ phy_gf8_mul(s2, b));
                uint8_t d = (uint8_t)(phy_gf8_mul(t1, a) ^ phy_gf8_mul(t2, b));
                r[phy_pack3(a, b)] = phy_pack3(pic[c], pid[d]);
            }
        }

        phy_analyze_ame8(r, &cand);
        if ((cand.flags & req) == req && cand.r[0] == 0U && cand.anf_deg >= 2U) {
            bool dup = false;
            for (uint32_t k = 0U; k < found; ++k) {
                if (seen[k] == cand.merkle) {
                    dup = true;
                    break;
                }
            }
            if (!dup && found < PHY_COUNT_AME8_META) {
                seen[found] = cand.merkle;
                phy_store_arena_ame8(arena, found, &cand);
                ++found;
            }
            /* Favor clean deg=2..3 rules for compact QTT linear .fir circuits */
            uint64_t low_deg_bonus = (cand.anf_deg <= 3U) ? UINT64_C(1) : UINT64_C(0);
            uint64_t sc = (low_deg_bonus << 40U) |
                          ((uint64_t)(cand.anf_deg == 2U ? 1U : 0U) << 36U) | (uint64_t)cand.order;
            if (!have_best || sc > best_sc) {
                have_best = true;
                best_sc = sc;
                phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
            }
        }
    }

    if (have_best && arena != NULL) {
        arena->map.img_hdr.active_rules = found;
        arena->map.img_hdr.merkle_root = out_best->merkle;
    }
    return found;
}

static uint64_t phy_frac32_apply_gate(const uint8_t t[64], uint64_t st, uint32_t pos) {
    uint32_t sh = pos * 2U;
    uint32_t w = (uint32_t)((st >> sh) & UINT64_C(0x3F));
    uint8_t a = (uint8_t)(w & 3U);
    uint8_t b = (uint8_t)((w >> 2U) & 3U);
    uint8_t c = (uint8_t)((w >> 4U) & 3U);
    uint8_t u = (uint8_t)(((uint32_t)a << 4U) | ((uint32_t)b << 2U) | (uint32_t)c);
    uint32_t o = (uint32_t)t[u];
    uint8_t oa = (uint8_t)((o >> 4U) & UINT32_C(3));
    uint8_t ob = (uint8_t)((o >> 2U) & UINT32_C(3));
    uint8_t oc = (uint8_t)(o & UINT32_C(3));
    uint64_t nw = (uint64_t)oa | ((uint64_t)ob << 2U) | ((uint64_t)oc << 4U);
    uint64_t clr = ~(UINT64_C(0x3F) << sh);
    return (st & clr) | (nw << sh);
}

static uint64_t phy_frac32_step(const uint8_t t[64], uint64_t st) {
    for (uint32_t layer = 0U; layer < 3U; ++layer) {
        for (uint32_t pos = layer; pos + 2U < 32U; pos += 3U) {
            st = phy_frac32_apply_gate(t, st, pos);
        }
    }
    return st;
}

static void phy_frac_qd(uint64_t st, uint32_t n_sites, int64_t *out_q, int64_t *out_d) {
    int64_t q = INT64_C(0);
    int64_t d = INT64_C(0);
    for (uint32_t i = 0U; i < n_sites; ++i) {
        int64_t v = (int64_t)((st >> (i * 2U)) & UINT64_C(3));
        q += v;
        d += (int64_t)i * v;
    }
    *out_q = q;
    *out_d = d;
}

static uint16_t phy_frac_apply6(const uint8_t t[64], uint16_t st, uint32_t pos) {
    uint64_t next = phy_frac32_apply_gate(t, (uint64_t)st, pos);
    return (uint16_t)(next & UINT64_C(0xFFF));
}

void phy_analyze_fracton(const uint8_t t[64], phy_fracton_hsf_ws_t *ws,
                         phy_fracton_analysis_t *out) {
    if (t == NULL || ws == NULL || out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    uint64_t seen = UINT64_C(0);
    uint16_t act = 0U;
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint8_t v = (uint8_t)(t[u] & 63U);
        out->t[u] = v;
        out->tinv[v] = u;
        seen |= (UINT64_C(1) << v);
        if (v != u) {
            ++act;
        }
    }
    out->is_bijective = (seen == UINT64_C(0xFFFFFFFFFFFFFFFF)) ? UINT8_C(1) : UINT8_C(0);
    out->active_entries = act;
    out->anf_deg = phy_anf_degree6(out->t, &out->nl_bit, &out->nl_a, &out->nl_b);
    out->order = (out->is_bijective != 0U) ? phy_perm64_order(out->t) : UINT32_C(0);

    bool lap_ok = true;
    bool got_wit = false;
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint32_t uu = (uint32_t)u;
        int32_t a = (int32_t)((uu >> 4U) & UINT32_C(3));
        int32_t b = (int32_t)((uu >> 2U) & UINT32_C(3));
        int32_t c = (int32_t)(uu & UINT32_C(3));
        uint32_t o = (uint32_t)out->t[u];
        int32_t oa = (int32_t)((o >> 4U) & UINT32_C(3));
        int32_t ob = (int32_t)((o >> 2U) & UINT32_C(3));
        int32_t oc = (int32_t)(o & UINT32_C(3));
        int32_t phi = oa - a;
        if ((ob - b) != -(2 * phi) || (oc - c) != phi) {
            lap_ok = false;
        } else if (phi != 0 && !got_wit) {
            got_wit = true;
            out->wit_u = u;
        }
    }
    out->laplacian_ok = lap_ok ? UINT8_C(1) : UINT8_C(0);

    /* Exhaustive L=6 (4096-state) (Q,D) check & largest-sector HSF BFS */
    uint16_t sec_cnt[19][46];
    for (uint32_t q = 0U; q < 19U; ++q) {
        for (uint32_t d = 0U; d < 46U; ++d) {
            sec_cnt[q][d] = 0U;
        }
    }
    bool qd_all_ok = (out->is_bijective != 0U);
    uint8_t best_q = 0U;
    uint8_t best_d = 0U;
    uint16_t best_sz = 0U;
    for (uint32_t st = 0U; st < 4096U; ++st) {
        int64_t q0 = INT64_C(0);
        int64_t d0 = INT64_C(0);
        phy_frac_qd((uint64_t)st, 6U, &q0, &d0);
        if (q0 >= 0 && q0 < 19 && d0 >= 0 && d0 < 46) {
            sec_cnt[q0][d0] += 1U;
            if (sec_cnt[q0][d0] > best_sz) {
                best_sz = sec_cnt[q0][d0];
                best_q = (uint8_t)q0;
                best_d = (uint8_t)d0;
            }
        }
        for (uint32_t pos = 0U; pos < 4U && qd_all_ok; ++pos) {
            uint16_t nst = phy_frac_apply6(out->t, (uint16_t)st, pos);
            int64_t q1 = INT64_C(0);
            int64_t d1 = INT64_C(0);
            phy_frac_qd((uint64_t)nst, 6U, &q1, &d1);
            if (q1 != q0 || d1 != d0) {
                qd_all_ok = false;
            }
        }
    }
    out->cons_qd_4096_ok = qd_all_ok ? UINT8_C(1) : UINT8_C(0);
    out->hsf_target_q = best_q;
    out->hsf_target_d = best_d;
    out->hsf_sector_size = best_sz;

    uint32_t n_sec = 0U;
    for (uint32_t st = 0U; st < 4096U; ++st) {
        int64_t q0 = INT64_C(0);
        int64_t d0 = INT64_C(0);
        phy_frac_qd((uint64_t)st, 6U, &q0, &d0);
        if (q0 == (int64_t)best_q && d0 == (int64_t)best_d && n_sec < 128U) {
            ws->sec_states[n_sec++] = (uint16_t)st;
        }
    }
    phy_vzero_u8(ws->vis, 4096U);
    uint16_t n_comps = 0U;
    uint16_t max_comp = 0U;
    uint16_t frozen_cnt = 0U;
    for (uint32_t i = 0U; i < n_sec; ++i) {
        uint16_t root = ws->sec_states[i];
        if (ws->vis[root] != 0U) {
            continue;
        }
        uint32_t head = 0U;
        uint32_t tail = 0U;
        ws->vis[root] = 1U;
        ws->queue[tail++] = root;
        while (head < tail) {
            uint16_t cur = ws->queue[head++];
            for (uint32_t pos = 0U; pos < 4U; ++pos) {
                uint16_t n1 = phy_frac_apply6(out->t, cur, pos);
                if (ws->vis[n1] == 0U) {
                    ws->vis[n1] = 1U;
                    ws->queue[tail++] = n1;
                }
                uint16_t n2 = phy_frac_apply6(out->tinv, cur, pos);
                if (ws->vis[n2] == 0U) {
                    ws->vis[n2] = 1U;
                    ws->queue[tail++] = n2;
                }
            }
        }
        uint16_t csz = (uint16_t)tail;
        ++n_comps;
        if (csz > max_comp) {
            max_comp = csz;
        }
        if (csz == 1U) {
            ++frozen_cnt;
        }
    }
    out->hsf_num_comps = n_comps;
    out->hsf_largest_comp = max_comp;
    out->hsf_frozen = frozen_cnt;

    /* L=32 staggered circuit simulation: isolated single charge vs dipole probe */
    uint64_t iso_init = UINT64_C(1) << (16U * 2U);
    uint64_t dip_init = UINT64_C(2) << (16U * 2U);
    int64_t iso_q0 = INT64_C(0);
    int64_t iso_d0 = INT64_C(0);
    int64_t dip_q0 = INT64_C(0);
    int64_t dip_d0 = INT64_C(0);
    phy_frac_qd(iso_init, 32U, &iso_q0, &iso_d0);
    phy_frac_qd(dip_init, 32U, &dip_q0, &dip_d0);

    uint64_t s_iso = iso_init;
    uint64_t s_dip = dip_init;
    uint16_t max_iso_span = 0U;
    uint16_t max_dip_span = 0U;
    int64_t max_qd = INT64_C(0);
    int64_t max_dd = INT64_C(0);

    for (uint32_t step = 0U; step < 32U; ++step) {
        s_iso = phy_frac32_step(out->t, s_iso);
        s_dip = phy_frac32_step(out->t, s_dip);

        int64_t iq = INT64_C(0);
        int64_t id = INT64_C(0);
        int64_t dq = INT64_C(0);
        int64_t dd = INT64_C(0);
        phy_frac_qd(s_iso, 32U, &iq, &id);
        phy_frac_qd(s_dip, 32U, &dq, &dd);

        uint64_t e1 = phy_abs_i64(iq - iso_q0);
        uint64_t e2 = phy_abs_i64(dq - dip_q0);
        uint64_t e3 = phy_abs_i64(id - iso_d0);
        uint64_t e4 = phy_abs_i64(dd - dip_d0);
        if ((int64_t)e1 > max_qd) {
            max_qd = (int64_t)e1;
        }
        if ((int64_t)e2 > max_qd) {
            max_qd = (int64_t)e2;
        }
        if ((int64_t)e3 > max_dd) {
            max_dd = (int64_t)e3;
        }
        if ((int64_t)e4 > max_dd) {
            max_dd = (int64_t)e4;
        }

        for (uint32_t site = 0U; site < 32U; ++site) {
            uint16_t dist = (site >= 16U) ? (uint16_t)(site - 16U) : (uint16_t)(16U - site);
            if (((s_iso >> (site * 2U)) & UINT64_C(3)) !=
                ((iso_init >> (site * 2U)) & UINT64_C(3))) {
                if (dist > max_iso_span) {
                    max_iso_span = dist;
                }
            }
            if (((s_dip >> (site * 2U)) & UINT64_C(3)) !=
                ((dip_init >> (site * 2U)) & UINT64_C(3))) {
                if (dist > max_dip_span) {
                    max_dip_span = dist;
                }
            }
        }
    }
    out->iso_span = max_iso_span;
    out->dip_span = max_dip_span;
    out->max_q_drift = max_qd;
    out->max_d_drift = max_dd;

    uint32_t fl = UINT32_C(0);
    if (out->is_bijective != 0U) {
        fl |= PHY_FRACTON_FLAG_BIJECTIVE;
    }
    if (out->laplacian_ok != 0U && out->cons_qd_4096_ok != 0U && max_qd == 0 && max_dd == 0) {
        fl |= PHY_FRACTON_FLAG_LAPLACIAN;
    }
    if (got_wit && out->active_entries > 0U) {
        fl |= PHY_FRACTON_FLAG_ACTIVE;
    }
    if (out->hsf_num_comps >= 8U && (uint32_t)out->hsf_largest_comp * 2U < out->hsf_sector_size &&
        out->hsf_frozen >= 4U) {
        fl |= PHY_FRACTON_FLAG_STRONG_HSF;
    }
    if (out->iso_span == 0U) {
        fl |= PHY_FRACTON_FLAG_ISO_FROZEN;
    }
    if (out->dip_span > 0U) {
        fl |= PHY_FRACTON_FLAG_DIP_MOBILE;
    }
    if (out->anf_deg >= 3U) {
        fl |= PHY_FRACTON_FLAG_NONLINEAR;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    for (uint32_t w = 0U; w < 8U; ++w) {
        uint64_t blk = UINT64_C(0);
        for (uint32_t k = 0U; k < 8U; ++k) {
            blk |= ((uint64_t)out->t[w * 8U + k]) << (k * 8U);
        }
        mh = phy_fnv1a64_u64(mh, blk);
    }
    mh = phy_fnv1a64_u64(mh, ((uint64_t)fl << 32U) | (uint64_t)out->hsf_num_comps);
    out->merkle = mh;
}

static void phy_store_arena_fracton(phy_arena_u *arena, uint32_t idx,
                                    const phy_fracton_analysis_t *fr) {
    if (arena == NULL || idx >= PHY_COUNT_FRACTON_META) {
        return;
    }
    uint32_t slot = PHY_SLOT_FRACTON_META + idx;
    phy_vzero_u64(arena->cells[slot].u64, 8U);
    arena->map.fracton_meta[idx].merkle = fr->merkle;
    arena->map.fracton_meta[idx].flags = fr->flags;
    arena->map.fracton_meta[idx].active_entries = fr->active_entries;
    arena->map.fracton_meta[idx].hsf_sector_size = fr->hsf_sector_size;
    arena->map.fracton_meta[idx].hsf_num_comps = fr->hsf_num_comps;
    arena->map.fracton_meta[idx].hsf_largest_comp = fr->hsf_largest_comp;
    arena->map.fracton_meta[idx].hsf_frozen = fr->hsf_frozen;
    arena->map.fracton_meta[idx].iso_span = fr->iso_span;
    arena->map.fracton_meta[idx].dip_span = fr->dip_span;
    arena->map.fracton_meta[idx].anf_deg = fr->anf_deg;
    arena->map.fracton_meta[idx].is_bijective = fr->is_bijective;
    arena->map.fracton_meta[idx].laplacian_ok = fr->laplacian_ok;
    arena->map.fracton_meta[idx].wit_u = fr->wit_u;
    arena->map.fracton_meta[idx].hsf_target_q = fr->hsf_target_q;
    arena->map.fracton_meta[idx].hsf_target_d = fr->hsf_target_d;
    phy_store_arena_v3_rule6(arena, idx, fr->t);
}

uint32_t phy_hunt_fracton(phy_arena_u *arena, phy_fracton_hsf_ws_t *ws, uint64_t seed,
                          uint64_t max_trials, phy_fracton_analysis_t *out_best) {
    if (ws == NULL || out_best == NULL) {
        return 0U;
    }
    phy_vzero_u8((uint8_t *)out_best, sizeof(*out_best));
    if (arena != NULL) {
        phy_arena_init(arena, seed);
    }

    uint8_t pair_u[18];
    uint8_t pair_v[18];
    uint8_t pair_is_s1[18];
    uint8_t pair_even_q[18];
    uint32_t np = 0U;
    uint32_t core_idx = 0U;
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint32_t uu = (uint32_t)u;
        uint32_t a = (uu >> 4U) & UINT32_C(3);
        uint32_t b = (uu >> 2U) & UINT32_C(3);
        uint32_t c = uu & UINT32_C(3);
        if (a + 1U <= 3U && b >= 2U && c + 1U <= 3U) {
            uint8_t v = (uint8_t)(((a + 1U) << 4U) | ((b - 2U) << 2U) | (c + 1U));
            if (np < 18U) {
                pair_u[np] = u;
                pair_v[np] = v;
                pair_is_s1[np] =
                    (a <= 2U && b <= 2U && c <= 2U && ((a + 1U) <= 2U) && ((c + 1U) <= 2U)) ? 1U
                                                                                            : 0U;
                pair_even_q[np] = (((a + b + c) & UINT32_C(1)) == UINT32_C(0)) ? 1U : 0U;
                if (u == 8U && v == 17U) {
                    core_idx = np;
                }
                ++np;
            }
        }
    }

    uint64_t rng = seed ^ UINT64_C(0x3C6EF372FE94F82B);
    if (rng == UINT64_C(0)) {
        rng = UINT64_C(0x9E3779B97F4A7C15);
    }

    uint64_t seen[PHY_COUNT_FRACTON_META];
    phy_vzero_u64(seen, PHY_COUNT_FRACTON_META);
    uint32_t found = 0U;
    uint64_t best_sc = UINT64_C(0);
    bool have_best = false;

    uint32_t req = PHY_FRACTON_FLAG_BIJECTIVE | PHY_FRACTON_FLAG_LAPLACIAN |
                   PHY_FRACTON_FLAG_ACTIVE | PHY_FRACTON_FLAG_STRONG_HSF |
                   PHY_FRACTON_FLAG_ISO_FROZEN | PHY_FRACTON_FLAG_DIP_MOBILE |
                   PHY_FRACTON_FLAG_NONLINEAR;

    phy_fracton_analysis_t cand;
    for (uint64_t trial = 0U; trial < max_trials; ++trial) {
        uint8_t t[64];
        for (uint8_t u = 0U; u < 64U; ++u) {
            t[u] = u;
        }
        uint32_t active_pairs = 0U;
        if (trial == UINT64_C(0)) {
            for (uint32_t i = 0U; i < np; ++i) {
                if (pair_is_s1[i] != 0U) {
                    t[pair_u[i]] = pair_v[i];
                    t[pair_v[i]] = pair_u[i];
                    ++active_pairs;
                }
            }
        } else if (trial == UINT64_C(1)) {
            for (uint32_t i = 0U; i < np; ++i) {
                if (pair_even_q[i] != 0U) {
                    t[pair_u[i]] = pair_v[i];
                    t[pair_v[i]] = pair_u[i];
                    ++active_pairs;
                }
            }
        } else {
            uint32_t mask = (uint32_t)(searu_rng_next(&rng) & UINT64_C(0x3FFFF));
            mask |= (UINT32_C(1) << core_idx);
            for (uint32_t i = 0U; i < np; ++i) {
                if ((mask & (UINT32_C(1) << i)) != UINT32_C(0)) {
                    t[pair_u[i]] = pair_v[i];
                    t[pair_v[i]] = pair_u[i];
                    ++active_pairs;
                }
            }
        }
        if (active_pairs < 4U) {
            continue;
        }

        phy_analyze_fracton(t, ws, &cand);
        if ((cand.flags & req) == req && cand.t[0] == 0U) {
            bool dup = false;
            for (uint32_t k = 0U; k < found; ++k) {
                if (seen[k] == cand.merkle) {
                    dup = true;
                    break;
                }
            }
            if (!dup && found < PHY_COUNT_FRACTON_META) {
                seen[found] = cand.merkle;
                phy_store_arena_fracton(arena, found, &cand);
                ++found;
            }
            uint64_t comp_bonus = ((uint64_t)cand.hsf_num_comps << 24U) |
                                  ((uint64_t)cand.hsf_frozen << 12U) |
                                  (uint64_t)(256U - cand.active_entries);
            uint64_t sc = ((uint64_t)cand.anf_deg << 40U) | comp_bonus;
            if (!have_best || sc > best_sc) {
                have_best = true;
                best_sc = sc;
                phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
            }
        }
    }

    if (have_best && arena != NULL) {
        arena->map.img_hdr.active_rules = found;
        arena->map.img_hdr.merkle_root = out_best->merkle;
    }
    return found;
}

static uint16_t phy_bq8_fix_pow(const uint8_t r[64], uint32_t k) {
    uint16_t cnt = 0U;
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint8_t x = u;
        for (uint32_t s = 0U; s < k; ++s) {
            x = r[x];
        }
        if (x == u) {
            ++cnt;
        }
    }
    return cnt;
}

static uint16_t phy_bq8_fig8(const uint8_t r[64], const uint8_t rinv[64]) {
    uint16_t cnt = 0U;
    for (uint8_t a = 0U; a < 8U; ++a) {
        for (uint8_t b = 0U; b < 8U; ++b) {
            for (uint8_t c = 0U; c < 8U; ++c) {
                uint8_t x = a;
                uint8_t y = b;
                uint8_t z = c;
                for (uint32_t s = 0U; s < 2U; ++s) {
                    uint8_t ab = r[phy_pack3(x, y)];
                    x = phy_hi3(ab);
                    y = phy_lo3(ab);
                    uint8_t bc = rinv[phy_pack3(y, z)];
                    y = phy_hi3(bc);
                    z = phy_lo3(bc);
                }
                if (x == a && y == b && z == c) {
                    ++cnt;
                }
            }
        }
    }
    return cnt;
}

void phy_analyze_biquandle8(const uint8_t r[64], phy_biquandle8_analysis_t *out) {
    if (r == NULL || out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    phy_analyze_yb8(r, &out->base);

    bool r1 = (out->base.is_bijective != 0U);
    for (uint8_t x = 0U; x < 8U && r1; ++x) {
        if (out->base.r[phy_pack3(x, x)] != phy_pack3(x, x)) {
            r1 = false;
        }
    }
    out->is_r1 = r1 ? UINT8_C(1) : UINT8_C(0);

    uint8_t ret_cls = 0U;
    for (uint8_t x = 0U; x < 8U; ++x) {
        bool dup = false;
        for (uint8_t y = 0U; y < x; ++y) {
            bool same = true;
            for (uint8_t k = 0U; k < 8U; ++k) {
                if (phy_hi3(out->base.r[phy_pack3(x, k)]) !=
                        phy_hi3(out->base.r[phy_pack3(y, k)]) ||
                    phy_lo3(out->base.r[phy_pack3(k, x)]) !=
                        phy_lo3(out->base.r[phy_pack3(k, y)])) {
                    same = false;
                    break;
                }
            }
            if (same) {
                dup = true;
                break;
            }
        }
        if (!dup) {
            ++ret_cls;
        }
    }
    out->retract_classes = ret_cls;
    out->is_irretract = (ret_cls == 8U) ? UINT8_C(1) : UINT8_C(0);

    out->col_unknot = phy_bq8_fix_pow(out->base.r, 1U);
    out->col_hopf = phy_bq8_fix_pow(out->base.r, 2U);
    out->col_trefoil = phy_bq8_fix_pow(out->base.r, 3U);
    out->col_torus5 = phy_bq8_fix_pow(out->base.r, 5U);
    out->col_torus7 = phy_bq8_fix_pow(out->base.r, 7U);
    out->col_fig8 = phy_bq8_fig8(out->base.r, out->base.rinv);

    uint32_t fl = out->base.flags;
    if (out->is_r1 != 0U) {
        fl |= PHY_BQ8_FLAG_R1;
    }
    if (out->is_irretract != 0U) {
        fl |= PHY_BQ8_FLAG_IRRETRACT;
    }
    if (out->col_unknot != out->col_torus7 || out->col_unknot != out->col_trefoil ||
        out->col_unknot != out->col_fig8) {
        fl |= PHY_BQ8_FLAG_KNOT_SEP;
    }
    out->base.flags = fl;
    out->base.merkle = phy_fnv1a64_u64(out->base.merkle, ((uint64_t)fl << 32U) |
                                                             ((uint64_t)out->col_torus7 << 16U) |
                                                             (uint64_t)out->col_trefoil);
}

static void phy_store_arena_bq8(phy_arena_u *arena, uint32_t idx,
                                const phy_biquandle8_analysis_t *bq) {
    if (arena == NULL || idx >= PHY_COUNT_BQ8_META) {
        return;
    }
    uint32_t slot = PHY_SLOT_BQ8_META + idx;
    phy_vzero_u64(arena->cells[slot].u64, 8U);
    arena->map.bq8_meta[idx].merkle = bq->base.merkle;
    arena->map.bq8_meta[idx].flags = bq->base.flags;
    arena->map.bq8_meta[idx].order = bq->base.order;
    arena->map.bq8_meta[idx].col_unknot = bq->col_unknot;
    arena->map.bq8_meta[idx].col_hopf = bq->col_hopf;
    arena->map.bq8_meta[idx].col_trefoil = bq->col_trefoil;
    arena->map.bq8_meta[idx].col_fig8 = bq->col_fig8;
    arena->map.bq8_meta[idx].col_torus5 = bq->col_torus5;
    arena->map.bq8_meta[idx].col_torus7 = bq->col_torus7;
    arena->map.bq8_meta[idx].anf_deg = bq->base.anf_deg;
    arena->map.bq8_meta[idx].retract_classes = bq->retract_classes;
    arena->map.bq8_meta[idx].is_r1 = bq->is_r1;
    arena->map.bq8_meta[idx].is_irretract = bq->is_irretract;
    phy_store_arena_v3_rule6(arena, idx, bq->base.r);
}

uint32_t phy_hunt_biquandle8(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                             phy_biquandle8_analysis_t *out_best) {
    if (out_best == NULL) {
        return 0U;
    }
    phy_vzero_u8((uint8_t *)out_best, sizeof(*out_best));
    if (arena != NULL) {
        phy_arena_init(arena, seed);
    }

    uint64_t rng = seed ^ UINT64_C(0x85EBCA6B49D76B91);
    if (rng == UINT64_C(0)) {
        rng = UINT64_C(0xC2B2AE3D27D4EB4F);
    }

    uint64_t seen[PHY_COUNT_BQ8_META];
    phy_vzero_u64(seen, PHY_COUNT_BQ8_META);
    uint32_t found = 0U;
    uint64_t best_sc = UINT64_C(0);
    bool have_best = false;

    uint32_t req = PHY_YB8_FLAG_BIJECTIVE | PHY_YB8_FLAG_YBE | PHY_YB8_FLAG_NONDEG |
                   PHY_YB8_FLAG_NONINVOL | PHY_YB8_FLAG_INDECOMP | PHY_YB8_FLAG_NONLYUB |
                   PHY_YB8_FLAG_NONLINEAR | PHY_BQ8_FLAG_R1 | PHY_BQ8_FLAG_IRRETRACT |
                   PHY_BQ8_FLAG_KNOT_SEP;

    phy_biquandle8_analysis_t cand;
    for (uint64_t trial = 0U; trial < max_trials; ++trial) {
        uint8_t s = (uint8_t)(1U + searu_Bounded_u32(searu_rng_next(&rng), 7U));
        uint8_t tt = (uint8_t)(1U + searu_Bounded_u32(searu_rng_next(&rng), 7U));
        uint8_t v = (uint8_t)(1U ^ s);
        uint8_t w = (uint8_t)(1U ^ tt);
        if (v == 0U && w == 0U) {
            s = 1U;
            tt = 2U;
            v = 0U;
            w = 3U;
        }

        uint8_t pi[8];
        uint8_t pi_inv[8];
        for (uint8_t x = 0U; x < 8U; ++x) {
            pi[x] = x;
        }
        if (trial == UINT64_C(0)) {
            pi[6] = 7U;
            pi[7] = 6U;
        } else {
            uint32_t nsw = 1U + searu_Bounded_u32(searu_rng_next(&rng), 2U);
            for (uint32_t k = 0U; k < nsw; ++k) {
                uint32_t i = searu_Bounded_u32(searu_rng_next(&rng), 8U);
                uint32_t j = searu_Bounded_u32(searu_rng_next(&rng), 8U);
                if (i == j) {
                    j = (j + 1U) & 7U;
                }
                uint8_t tmp = pi[i];
                pi[i] = pi[j];
                pi[j] = tmp;
            }
        }
        for (uint8_t x = 0U; x < 8U; ++x) {
            pi_inv[pi[x]] = x;
        }

        uint8_t r[64];
        for (uint8_t a = 0U; a < 8U; ++a) {
            uint8_t pa = pi[a];
            for (uint8_t b = 0U; b < 8U; ++b) {
                uint8_t pb = pi[b];
                uint8_t oa = (uint8_t)(phy_gf8_mul(s, pb) ^ phy_gf8_mul(v, pa));
                uint8_t ob = (uint8_t)(phy_gf8_mul(tt, pa) ^ phy_gf8_mul(w, pb));
                r[phy_pack3(a, b)] = phy_pack3(pi_inv[oa], pi_inv[ob]);
            }
        }

        phy_analyze_biquandle8(r, &cand);
        if ((cand.base.flags & req) == req && cand.col_unknot == 8U && cand.col_torus7 == 64U) {
            bool dup = false;
            for (uint32_t k = 0U; k < found; ++k) {
                if (seen[k] == cand.base.merkle) {
                    dup = true;
                    break;
                }
            }
            if (!dup && found < PHY_COUNT_BQ8_META) {
                seen[found] = cand.base.merkle;
                phy_store_arena_bq8(arena, found, &cand);
                ++found;
            }
            uint64_t deg_pref = (cand.base.anf_deg == 3U) ? UINT64_C(1) : UINT64_C(0);
            uint64_t sc = (deg_pref << 40U) | ((uint64_t)cand.base.anf_deg << 32U) |
                          ((uint64_t)cand.col_torus7 << 16U) | (uint64_t)cand.base.order;
            if (!have_best || sc > best_sc) {
                have_best = true;
                best_sc = sc;
                phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
            }
        }
    }

    if (have_best && arena != NULL) {
        arena->map.img_hdr.active_rules = found;
        arena->map.img_hdr.merkle_root = out_best->base.merkle;
    }
    return found;
}

static uint8_t phy_tetra2_apply(const uint8_t p[8], uint8_t st, uint32_t i, uint32_t j,
                                uint32_t k) {
    uint32_t si = 5U - i;
    uint32_t sj = 5U - j;
    uint32_t sk = 5U - k;
    uint32_t ust = (uint32_t)st;
    uint32_t b1 = (ust >> si) & UINT32_C(1);
    uint32_t b2 = (ust >> sj) & UINT32_C(1);
    uint32_t b3 = (ust >> sk) & UINT32_C(1);
    uint32_t o = (uint32_t)p[(b1 << 2U) | (b2 << 1U) | b3];
    uint32_t o1 = (o >> 2U) & UINT32_C(1);
    uint32_t o2 = (o >> 1U) & UINT32_C(1);
    uint32_t o3 = o & UINT32_C(1);
    uint32_t clr =
        ~((UINT32_C(1) << si) | (UINT32_C(1) << sj) | (UINT32_C(1) << sk)) & UINT32_C(0x3F);
    return (uint8_t)((ust & clr) | (o1 << si) | (o2 << sj) | (o3 << sk));
}

static void phy_tetra2_inspect(const uint8_t p[8], uint16_t *tot, uint16_t *nl, uint16_t *ni_nl) {
    for (uint8_t s = 0U; s < 64U; ++s) {
        uint8_t l = phy_tetra2_apply(
            p,
            phy_tetra2_apply(p, phy_tetra2_apply(p, phy_tetra2_apply(p, s, 2U, 4U, 5U), 1U, 3U, 5U),
                             0U, 3U, 4U),
            0U, 1U, 2U);
        uint8_t r = phy_tetra2_apply(
            p,
            phy_tetra2_apply(p, phy_tetra2_apply(p, phy_tetra2_apply(p, s, 0U, 1U, 2U), 0U, 3U, 4U),
                             1U, 3U, 5U),
            2U, 4U, 5U);
        if (l != r) {
            return;
        }
    }
    *tot += 1U;

    uint8_t md = 0U;
    for (uint8_t bit = 0U; bit < 3U; ++bit) {
        uint8_t anf[8];
        for (uint8_t x = 0U; x < 8U; ++x) {
            anf[x] = (uint8_t)(((uint32_t)p[x] >> bit) & 1U);
        }
        for (uint32_t v = 0U; v < 3U; ++v) {
            for (uint32_t m = 0U; m < 8U; ++m) {
                if ((m & (UINT32_C(1) << v)) != 0U) {
                    anf[m] ^= anf[m ^ (UINT32_C(1) << v)];
                }
            }
        }
        for (uint32_t m = 0U; m < 8U; ++m) {
            if (anf[m] != 0U) {
                uint8_t d = (uint8_t)searu_popcnt64((uint64_t)m);
                if (d > md) {
                    md = d;
                }
            }
        }
    }
    if (md >= 2U) {
        *nl += 1U;
        for (uint8_t x = 0U; x < 8U; ++x) {
            if (p[p[x]] != x) {
                *ni_nl += 1U;
                break;
            }
        }
    }
}

void phy_classify_tetra2(uint16_t *out_total, uint16_t *out_nl, uint16_t *out_noninvol_nl) {
    uint16_t tot = 0U;
    uint16_t nl = 0U;
    uint16_t ni_nl = 0U;

    uint8_t p[8];
    uint8_t c[8];
    for (uint8_t i = 0U; i < 8U; ++i) {
        p[i] = i;
        c[i] = 0U;
    }
    phy_tetra2_inspect(p, &tot, &nl, &ni_nl);

    uint32_t idx = 0U;
    while (idx < 8U) {
        if ((uint32_t)c[idx] < idx) {
            uint32_t sw = ((idx & 1U) == 0U) ? 0U : (uint32_t)c[idx];
            uint8_t tmp = p[sw];
            p[sw] = p[idx];
            p[idx] = tmp;
            phy_tetra2_inspect(p, &tot, &nl, &ni_nl);
            c[idx] += 1U;
            idx = 0U;
        } else {
            c[idx] = 0U;
            ++idx;
        }
    }
    if (out_total != NULL) {
        *out_total = tot;
    }
    if (out_nl != NULL) {
        *out_nl = nl;
    }
    if (out_noninvol_nl != NULL) {
        *out_noninvol_nl = ni_nl;
    }
}

static uint16_t phy_tetra4_apply(const uint8_t t[64], uint16_t st, uint32_t i, uint32_t j,
                                 uint32_t k) {
    uint32_t si = 2U * (5U - i);
    uint32_t sj = 2U * (5U - j);
    uint32_t sk = 2U * (5U - k);
    uint32_t ust = (uint32_t)st;
    uint32_t x1 = (ust >> si) & UINT32_C(3);
    uint32_t x2 = (ust >> sj) & UINT32_C(3);
    uint32_t x3 = (ust >> sk) & UINT32_C(3);
    uint32_t o = (uint32_t)t[(x1 << 4U) | (x2 << 2U) | x3];
    uint32_t o1 = (o >> 4U) & UINT32_C(3);
    uint32_t o2 = (o >> 2U) & UINT32_C(3);
    uint32_t o3 = o & UINT32_C(3);
    uint32_t clr =
        ~((UINT32_C(3) << si) | (UINT32_C(3) << sj) | (UINT32_C(3) << sk)) & UINT32_C(0xFFF);
    return (uint16_t)((ust & clr) | (o1 << si) | (o2 << sj) | (o3 << sk));
}

uint16_t phy_tetra4_check(const uint8_t t[64]) {
    uint16_t viol = 0U;
    for (uint16_t s = 0U; s < 4096U; ++s) {
        uint16_t l = phy_tetra4_apply(
            t,
            phy_tetra4_apply(t, phy_tetra4_apply(t, phy_tetra4_apply(t, s, 2U, 4U, 5U), 1U, 3U, 5U),
                             0U, 3U, 4U),
            0U, 1U, 2U);
        uint16_t r = phy_tetra4_apply(
            t,
            phy_tetra4_apply(t, phy_tetra4_apply(t, phy_tetra4_apply(t, s, 0U, 1U, 2U), 0U, 3U, 4U),
                             1U, 3U, 5U),
            2U, 4U, 5U);
        if (l != r) {
            ++viol;
        }
    }
    return viol;
}

void phy_analyze_tetra4(const uint8_t t[64], phy_tetra4_analysis_t *out) {
    if (t == NULL || out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    uint64_t seen = UINT64_C(0);
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint8_t v = (uint8_t)(t[u] & 63U);
        out->t[u] = v;
        out->tinv[v] = u;
        seen |= (UINT64_C(1) << v);
    }
    out->is_bijective = (seen == UINT64_C(0xFFFFFFFFFFFFFFFF)) ? UINT8_C(1) : UINT8_C(0);
    out->tetra4_viol = phy_tetra4_check(out->t);
    out->anf_deg = phy_anf_degree6(out->t, &out->nl_bit, &out->nl_a, &out->nl_b);
    out->order = (out->is_bijective != 0U) ? phy_perm64_order(out->t) : UINT32_C(0);

    out->is_involutive = 1U;
    for (uint8_t u = 0U; u < 64U; ++u) {
        if (out->t[out->t[u]] != u) {
            out->is_involutive = 0U;
            out->noninvol_u = u;
            break;
        }
    }

    /* 3-channel coupling check on middle output channel (site1) vs channel 0 and channel 2 */
    bool c_left = false;
    bool c_right = false;
    for (uint8_t b = 0U; b < 4U && !c_left; ++b) {
        for (uint8_t c = 0U; c < 4U && !c_left; ++c) {
            uint32_t idx0 = ((uint32_t)b << 2U) | (uint32_t)c;
            uint8_t m0 = (uint8_t)(((uint32_t)out->t[idx0] >> 2U) & UINT32_C(3));
            for (uint8_t a = 1U; a < 4U; ++a) {
                uint32_t idx1 = ((uint32_t)a << 4U) | idx0;
                uint8_t m1 = (uint8_t)(((uint32_t)out->t[idx1] >> 2U) & UINT32_C(3));
                if (m1 != m0) {
                    c_left = true;
                    out->c0_a0 = 0U;
                    out->c0_a1 = a;
                    out->c0_b = b;
                    out->c0_c = c;
                    break;
                }
            }
        }
    }
    for (uint8_t a = 0U; a < 4U && !c_right; ++a) {
        for (uint8_t b = 0U; b < 4U && !c_right; ++b) {
            uint32_t base = ((uint32_t)a << 4U) | ((uint32_t)b << 2U);
            uint8_t m0 = (uint8_t)(((uint32_t)out->t[base] >> 2U) & UINT32_C(3));
            for (uint8_t c = 1U; c < 4U; ++c) {
                uint8_t m1 = (uint8_t)(((uint32_t)out->t[base | (uint32_t)c] >> 2U) & UINT32_C(3));
                if (m1 != m0) {
                    c_right = true;
                    out->c2_a = a;
                    out->c2_b = b;
                    out->c2_c0 = 0U;
                    out->c2_c1 = c;
                    break;
                }
            }
        }
    }
    out->is_3coupled = (c_left && c_right) ? UINT8_C(1) : UINT8_C(0);

    uint32_t fl = UINT32_C(0);
    if (out->is_bijective != 0U) {
        fl |= PHY_TETRA4_FLAG_BIJECTIVE;
    }
    if (out->tetra4_viol == 0U) {
        fl |= PHY_TETRA4_FLAG_TETRA;
    }
    if (out->is_involutive == 0U) {
        fl |= PHY_TETRA4_FLAG_NONINVOL;
    }
    if (out->is_3coupled != 0U) {
        fl |= PHY_TETRA4_FLAG_3COUPLED;
    }
    if (out->anf_deg >= 2U) {
        fl |= PHY_TETRA4_FLAG_NONLINEAR;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    for (uint32_t w = 0U; w < 8U; ++w) {
        uint64_t blk = UINT64_C(0);
        for (uint32_t k = 0U; k < 8U; ++k) {
            blk |= ((uint64_t)out->t[w * 8U + k]) << (k * 8U);
        }
        mh = phy_fnv1a64_u64(mh, blk);
    }
    mh = phy_fnv1a64_u64(mh, ((uint64_t)fl << 32U) | (uint64_t)out->order);
    out->merkle = mh;
}

static void phy_store_arena_tetra4(phy_arena_u *arena, uint32_t idx,
                                   const phy_tetra4_analysis_t *t4) {
    if (arena == NULL || idx >= PHY_COUNT_TETRA4_META) {
        return;
    }
    uint32_t slot = PHY_SLOT_TETRA4_META + idx;
    phy_vzero_u64(arena->cells[slot].u64, 8U);
    arena->map.tetra4_meta[idx].merkle = t4->merkle;
    arena->map.tetra4_meta[idx].flags = t4->flags;
    arena->map.tetra4_meta[idx].order = t4->order;
    arena->map.tetra4_meta[idx].tetra4_viol = t4->tetra4_viol;
    arena->map.tetra4_meta[idx].tetra2_total = t4->tetra2_total;
    arena->map.tetra4_meta[idx].tetra2_nl = t4->tetra2_nl;
    arena->map.tetra4_meta[idx].tetra2_noninvol_nl = t4->tetra2_noninvol_nl;
    arena->map.tetra4_meta[idx].anf_deg = t4->anf_deg;
    arena->map.tetra4_meta[idx].is_bijective = t4->is_bijective;
    arena->map.tetra4_meta[idx].is_involutive = t4->is_involutive;
    arena->map.tetra4_meta[idx].is_3coupled = t4->is_3coupled;
    arena->map.tetra4_meta[idx].noninvol_u = t4->noninvol_u;
    phy_store_arena_v3_rule6(arena, idx, t4->t);
}

uint32_t phy_hunt_tetra4(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                         phy_tetra4_analysis_t *out_best) {
    if (out_best == NULL) {
        return 0U;
    }
    phy_vzero_u8((uint8_t *)out_best, sizeof(*out_best));
    if (arena != NULL) {
        phy_arena_init(arena, seed);
    }

    uint16_t t2_tot = 0U;
    uint16_t t2_nl = 0U;
    uint16_t t2_ni = 0U;
    phy_classify_tetra2(&t2_tot, &t2_nl, &t2_ni);

    uint8_t s4[24][4];
    uint8_t fix0_idx[6];
    phy_gen_s4_perms(s4, fix0_idx);

    uint64_t rng = seed ^ UINT64_C(0xBB67AE8584CAA73B);
    if (rng == UINT64_C(0)) {
        rng = UINT64_C(0x3C6EF372FE94F82A);
    }

    uint64_t seen[PHY_COUNT_TETRA4_META];
    phy_vzero_u64(seen, PHY_COUNT_TETRA4_META);
    uint32_t found = 0U;
    uint64_t best_sc = UINT64_C(0);
    bool have_best = false;

    uint32_t req = PHY_TETRA4_FLAG_BIJECTIVE | PHY_TETRA4_FLAG_TETRA | PHY_TETRA4_FLAG_NONINVOL |
                   PHY_TETRA4_FLAG_3COUPLED | PHY_TETRA4_FLAG_NONLINEAR;

    phy_tetra4_analysis_t cand;
    for (uint64_t trial = 0U; trial < max_trials; ++trial) {
        uint8_t k = (uint8_t)(((searu_rng_next(&rng) & 1U) != 0U) ? 3U : 1U);
        uint8_t c_off = 0U;
        uint8_t p_id = 0U;
        if (trial == UINT64_C(0)) {
            k = 1U;
            p_id = fix0_idx[0];
        } else {
            p_id = fix0_idx[searu_Bounded_u32(searu_rng_next(&rng), 6U)];
        }
        uint8_t pi[4];
        uint8_t pi_inv[4];
        for (uint8_t x = 0U; x < 4U; ++x) {
            pi[x] = s4[p_id][x];
            pi_inv[pi[x]] = x;
        }

        uint8_t t[64];
        for (uint8_t a = 0U; a < 4U; ++a) {
            for (uint8_t b = 0U; b < 4U; ++b) {
                for (uint8_t c = 0U; c < 4U; ++c) {
                    uint32_t ob =
                        ((uint32_t)pi[b] + (uint32_t)k * pi[a] * pi[c] + (uint32_t)c_off) &
                        UINT32_C(3);
                    uint32_t tidx = ((uint32_t)a << 4U) | ((uint32_t)b << 2U) | (uint32_t)c;
                    t[tidx] =
                        (uint8_t)(((uint32_t)a << 4U) | ((uint32_t)pi_inv[ob] << 2U) | (uint32_t)c);
                }
            }
        }

        if (phy_tetra4_check(t) != 0U) {
            continue;
        }
        phy_analyze_tetra4(t, &cand);
        cand.tetra2_total = t2_tot;
        cand.tetra2_nl = t2_nl;
        cand.tetra2_noninvol_nl = t2_ni;

        if ((cand.flags & req) == req && cand.anf_deg >= 3U && cand.t[0] == 0U) {
            bool dup = false;
            for (uint32_t k_i = 0U; k_i < found; ++k_i) {
                if (seen[k_i] == cand.merkle) {
                    dup = true;
                    break;
                }
            }
            if (!dup && found < PHY_COUNT_TETRA4_META) {
                seen[found] = cand.merkle;
                phy_store_arena_tetra4(arena, found, &cand);
                ++found;
            }
            /* Favor canonical Z4 identity conjugate (fast bv_decide + minimal ANDs) */
            uint64_t clean_bonus = (trial == UINT64_C(0)) ? UINT64_C(1) : UINT64_C(0);
            uint64_t sc =
                (clean_bonus << 40U) | ((uint64_t)cand.anf_deg << 24U) | (uint64_t)cand.order;
            if (!have_best || sc > best_sc) {
                have_best = true;
                best_sc = sc;
                phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
            }
        }
    }

    if (have_best && arena != NULL) {
        arena->map.img_hdr.active_rules = found;
        arena->map.img_hdr.merkle_root = out_best->merkle;
    }
    return found;
}

/* --- PHY v4.0: Exact Rigidity, 3+1D 4-Simplex, Chiral Cohomology, Z[i] Scar, Fracton Univ --- */

static uint32_t phy_ctz_u64(uint64_t x) {
    if (x == UINT64_C(0)) {
        return UINT32_C(64);
    }
    uint32_t tz = UINT32_C(0);
    while (((x >> tz) & UINT64_C(1)) == UINT64_C(0) && tz < UINT32_C(63)) {
        ++tz;
    }
    return tz;
}

static uint8_t phy_anf_degree8(const uint8_t s[256]) {
    uint8_t max_d = 0U;
    uint8_t f[256];
    for (uint32_t bit = 0U; bit < 8U; ++bit) {
        for (uint32_t x = 0U; x < 256U; ++x) {
            f[x] = (uint8_t)(((uint32_t)s[x] >> bit) & UINT32_C(1));
        }
        for (uint32_t v = 0U; v < 8U; ++v) {
            for (uint32_t mask = 0U; mask < 256U; ++mask) {
                if ((mask & (UINT32_C(1) << v)) != UINT32_C(0)) {
                    f[mask] = (uint8_t)((uint32_t)f[mask] ^ (uint32_t)f[mask ^ (UINT32_C(1) << v)]);
                }
            }
        }
        for (uint32_t mask = 0U; mask < 256U; ++mask) {
            if (f[mask] != 0U) {
                uint8_t d = 0U;
                uint32_t m = mask;
                while (m != UINT32_C(0)) {
                    d += (uint8_t)(m & UINT32_C(1));
                    m >>= UINT32_C(1);
                }
                if (d > max_d) {
                    max_d = d;
                }
            }
        }
    }
    return max_d;
}

static uint32_t phy_perm256_order(const uint8_t s[256]) {
    uint8_t vis[256];
    phy_vzero_u8(vis, 256U);
    uint64_t order = UINT64_C(1);
    for (uint32_t i = 0U; i < 256U; ++i) {
        if (vis[i] != 0U) {
            continue;
        }
        uint64_t clen = UINT64_C(0);
        uint32_t cur = i;
        while (vis[cur] == 0U && clen < UINT64_C(256)) {
            vis[cur] = 1U;
            cur = s[cur];
            ++clen;
        }
        uint64_t g = phy_stein_gcd64(order, clen);
        uint64_t q = phy_div_u64(order, (g == UINT64_C(0)) ? UINT64_C(1) : g, NULL);
        order = q * clen;
    }
    return (uint32_t)order;
}

bool phy_is_2unitary4(const uint8_t p[16]) {
    if (p == NULL) {
        return false;
    }
    uint32_t seen = UINT32_C(0);
    for (uint32_t i = 0U; i < 16U; ++i) {
        uint32_t v = (uint32_t)p[i] & UINT32_C(0xF);
        seen |= (UINT32_C(1) << v);
    }
    if (seen != UINT32_C(0xFFFF)) {
        return false;
    }
    for (uint32_t a = 0U; a < 4U; ++a) {
        uint32_t rc_m = UINT32_C(0);
        uint32_t rd_m = UINT32_C(0);
        uint32_t cc_m = UINT32_C(0);
        uint32_t cd_m = UINT32_C(0);
        for (uint32_t b = 0U; b < 4U; ++b) {
            uint32_t rv = (uint32_t)p[(a << 2U) | b] & UINT32_C(0xF);
            uint32_t cv = (uint32_t)p[(b << 2U) | a] & UINT32_C(0xF);
            rc_m |= (UINT32_C(1) << ((rv >> 2U) & UINT32_C(3)));
            rd_m |= (UINT32_C(1) << (rv & UINT32_C(3)));
            cc_m |= (UINT32_C(1) << ((cv >> 2U) & UINT32_C(3)));
            cd_m |= (UINT32_C(1) << (cv & UINT32_C(3)));
        }
        if (rc_m != UINT32_C(0xF) || rd_m != UINT32_C(0xF) || cc_m != UINT32_C(0xF) ||
            cd_m != UINT32_C(0xF)) {
            return false;
        }
    }
    return true;
}

void phy_classify_ame4(phy_arena_u *arena, phy_ame4_nogo_t *out) {
    if (out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));

    uint32_t ls[576];
    uint32_t ls_cnt = UINT32_C(0);
    uint8_t grid[16];
    uint8_t row_m[4];
    uint8_t col_m[4];
    int32_t next_v[16];
    phy_vzero_u8(grid, 16U);
    phy_vzero_u8(row_m, 4U);
    phy_vzero_u8(col_m, 4U);
    for (uint32_t i = 0U; i < 16U; ++i) {
        next_v[i] = 0;
    }

    int32_t pos = 0;
    while (pos >= 0) {
        if (pos == 16) {
            if (ls_cnt < UINT32_C(576)) {
                uint32_t packed = UINT32_C(0);
                for (uint32_t i = 0U; i < 16U; ++i) {
                    packed |= ((uint32_t)grid[i] & UINT32_C(3)) << (i * UINT32_C(2));
                }
                ls[ls_cnt++] = packed;
            }
            --pos;
            uint32_t rp = (uint32_t)pos >> 2U;
            uint32_t cp = (uint32_t)pos & UINT32_C(3);
            uint32_t bit = UINT32_C(1) << ((uint32_t)grid[pos] & UINT32_C(3));
            row_m[rp] = (uint8_t)((uint32_t)row_m[rp] ^ bit);
            col_m[cp] = (uint8_t)((uint32_t)col_m[cp] ^ bit);
            continue;
        }

        uint32_t r = (uint32_t)pos >> 2U;
        uint32_t c = (uint32_t)pos & UINT32_C(3);
        uint32_t used = (uint32_t)row_m[r] | (uint32_t)col_m[c];
        int32_t v = next_v[pos];
        while (v < 4 && (used & (UINT32_C(1) << (uint32_t)v)) != UINT32_C(0)) {
            ++v;
        }
        if (v < 4) {
            grid[pos] = (uint8_t)v;
            uint32_t bit = UINT32_C(1) << (uint32_t)v;
            row_m[r] = (uint8_t)((uint32_t)row_m[r] | bit);
            col_m[c] = (uint8_t)((uint32_t)col_m[c] | bit);
            next_v[pos] = v + 1;
            ++pos;
            if (pos < 16) {
                next_v[pos] = 0;
            }
        } else {
            next_v[pos] = 0;
            --pos;
            if (pos >= 0) {
                uint32_t rp = (uint32_t)pos >> 2U;
                uint32_t cp = (uint32_t)pos & UINT32_C(3);
                uint32_t bit = UINT32_C(1) << ((uint32_t)grid[pos] & UINT32_C(3));
                row_m[rp] = (uint8_t)((uint32_t)row_m[rp] ^ bit);
                col_m[cp] = (uint8_t)((uint32_t)col_m[cp] ^ bit);
            }
        }
    }
    out->latin4_count = ls_cnt;

    uint32_t total_ols = UINT32_C(0);
    uint32_t affine_ols = UINT32_C(0);
    uint32_t nonlin_ols = UINT32_C(0);
    uint32_t bern_ols = UINT32_C(0);
    bool have_sample = false;

    for (uint32_t i = 0U; i < ls_cnt; ++i) {
        uint32_t ci = ls[i];
        for (uint32_t j = 0U; j < ls_cnt; ++j) {
            uint32_t dj = ls[j];
            uint32_t seen_mask = UINT32_C(0);
            uint8_t p[16];
            for (uint32_t k = 0U; k < 16U; ++k) {
                uint32_t ck = (ci >> (k * UINT32_C(2))) & UINT32_C(3);
                uint32_t dk = (dj >> (k * UINT32_C(2))) & UINT32_C(3);
                uint32_t val = (ck << UINT32_C(2)) | dk;
                p[k] = (uint8_t)val;
                seen_mask |= (UINT32_C(1) << val);
            }
            if (seen_mask != UINT32_C(0xFFFF)) {
                continue;
            }
            ++total_ols;
            uint8_t deg = phy_anf_degree4(p, NULL, NULL, NULL);
            if (deg == 1U) {
                ++affine_ols;
            } else {
                ++nonlin_ols;
            }
            ++bern_ols;
            if (!have_sample && p[0] == 0U && phy_is_2unitary4(p)) {
                have_sample = true;
                phy_vcopy_u8(out->sample_p, p, 16U);
                out->sample_lut = phy_pack_lut(p);
            }
        }
    }

    out->total_ols = total_ols;
    out->affine_ols = affine_ols;
    out->nonlin_ols = nonlin_ols;
    out->bernoulli_ols = bern_ols;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    mh = phy_fnv1a64_u64(mh, out->sample_lut);
    mh = phy_fnv1a64_u64(mh, ((uint64_t)total_ols << 32U) | (uint64_t)affine_ols);
    out->merkle = mh;

    if (arena != NULL) {
        phy_vzero_u64(arena->cells[PHY_SLOT_AME4_NOGO].u64, 8U);
        arena->map.ame4_nogo.sample_lut = out->sample_lut;
        arena->map.ame4_nogo.merkle = out->merkle;
        arena->map.ame4_nogo.latin4_count = out->latin4_count;
        arena->map.ame4_nogo.total_ols = out->total_ols;
        arena->map.ame4_nogo.affine_ols = out->affine_ols;
        arena->map.ame4_nogo.nonlin_ols = out->nonlin_ols;
        arena->map.ame4_nogo.bernoulli_ols = out->bernoulli_ols;
    }
}

static uint16_t phy_simplex4_apply2(const uint8_t s[16], uint16_t st, uint32_t i, uint32_t j,
                                    uint32_t k, uint32_t l) {
    uint32_t si = UINT32_C(9) - i;
    uint32_t sj = UINT32_C(9) - j;
    uint32_t sk = UINT32_C(9) - k;
    uint32_t sl = UINT32_C(9) - l;
    uint32_t ust = (uint32_t)st;
    uint32_t u = (((ust >> si) & UINT32_C(1)) << UINT32_C(3)) |
                 (((ust >> sj) & UINT32_C(1)) << UINT32_C(2)) |
                 (((ust >> sk) & UINT32_C(1)) << UINT32_C(1)) | ((ust >> sl) & UINT32_C(1));
    uint32_t o = (uint32_t)s[u] & UINT32_C(0xF);
    uint32_t clr =
        ~((UINT32_C(1) << si) | (UINT32_C(1) << sj) | (UINT32_C(1) << sk) | (UINT32_C(1) << sl)) &
        UINT32_C(0x3FF);
    return (uint16_t)((ust & clr) | (((o >> 3U) & UINT32_C(1)) << si) |
                      (((o >> 2U) & UINT32_C(1)) << sj) | (((o >> 1U) & UINT32_C(1)) << sk) |
                      ((o & UINT32_C(1)) << sl));
}

uint16_t phy_simplex4_check2(const uint8_t s[16], uint16_t *out_first_bad) {
    if (s == NULL) {
        return UINT16_C(0xFFFF);
    }
    uint16_t viol = UINT16_C(0);
    bool got = false;
    for (uint16_t st = UINT16_C(0); st < UINT16_C(1024); ++st) {
        uint16_t l = st;
        l = phy_simplex4_apply2(s, l, 6U, 7U, 8U, 9U);
        l = phy_simplex4_apply2(s, l, 3U, 4U, 5U, 9U);
        l = phy_simplex4_apply2(s, l, 1U, 2U, 5U, 8U);
        l = phy_simplex4_apply2(s, l, 0U, 2U, 4U, 7U);
        l = phy_simplex4_apply2(s, l, 0U, 1U, 3U, 6U);

        uint16_t r = st;
        r = phy_simplex4_apply2(s, r, 0U, 1U, 3U, 6U);
        r = phy_simplex4_apply2(s, r, 0U, 2U, 4U, 7U);
        r = phy_simplex4_apply2(s, r, 1U, 2U, 5U, 8U);
        r = phy_simplex4_apply2(s, r, 3U, 4U, 5U, 9U);
        r = phy_simplex4_apply2(s, r, 6U, 7U, 8U, 9U);
        if (l != r) {
            if (!got && out_first_bad != NULL) {
                got = true;
                *out_first_bad = st;
            }
            ++viol;
        }
    }
    return viol;
}

void phy_analyze_simplex4_2(phy_simplex4_2_t *out) {
    if (out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    for (uint8_t u = 0U; u < 16U; ++u) {
        uint32_t a = ((uint32_t)u >> 3U) & UINT32_C(1);
        uint32_t b = ((uint32_t)u >> 2U) & UINT32_C(1);
        uint32_t c = ((uint32_t)u >> 1U) & UINT32_C(1);
        uint32_t d = (uint32_t)u & UINT32_C(1);
        out->ccx02[u] = (uint8_t)((uint32_t)u ^ ((a & c) << 2U));
        out->ccx13[u] = (uint8_t)((uint32_t)u ^ ((b & d) << 1U));
        out->cccx1[u] = (uint8_t)((uint32_t)u ^ ((a & c & d) << 2U));
    }
    out->viol2_ccx02 = phy_simplex4_check2(out->ccx02, NULL);
    out->viol2_ccx13 = phy_simplex4_check2(out->ccx13, NULL);
    out->viol2_cccx1 = phy_simplex4_check2(out->cccx1, &out->cccx1_wit_s);

    uint8_t nl = 0U;
    for (uint32_t mask_b = 0U; mask_b < 256U; ++mask_b) {
        uint8_t s1[16];
        for (uint32_t u = 0U; u < 16U; ++u) {
            uint32_t idx = (((u >> 3U) & UINT32_C(1)) << 2U) | (((u >> 1U) & UINT32_C(1)) << 1U) |
                           (u & UINT32_C(1));
            s1[u] = (uint8_t)(u ^ (((mask_b >> idx) & UINT32_C(1)) << 2U));
        }
        if (phy_anf_degree4(s1, NULL, NULL, NULL) >= 2U && phy_simplex4_check2(s1, NULL) == 0U) {
            ++nl;
        }
    }
    out->nl_count_b_ctrl = nl;
}

static uint32_t phy_simplex4_apply4(const uint8_t s[256], uint32_t st, uint32_t i, uint32_t j,
                                    uint32_t k, uint32_t l) {
    uint32_t si = UINT32_C(2) * (UINT32_C(9) - i);
    uint32_t sj = UINT32_C(2) * (UINT32_C(9) - j);
    uint32_t sk = UINT32_C(2) * (UINT32_C(9) - k);
    uint32_t sl = UINT32_C(2) * (UINT32_C(9) - l);
    uint32_t u = (((st >> si) & UINT32_C(3)) << 6U) | (((st >> sj) & UINT32_C(3)) << 4U) |
                 (((st >> sk) & UINT32_C(3)) << 2U) | ((st >> sl) & UINT32_C(3));
    uint32_t o = (uint32_t)s[u];
    uint32_t clr =
        ~((UINT32_C(3) << si) | (UINT32_C(3) << sj) | (UINT32_C(3) << sk) | (UINT32_C(3) << sl)) &
        UINT32_C(0xFFFFF);
    return (st & clr) | (((o >> 6U) & UINT32_C(3)) << si) | (((o >> 4U) & UINT32_C(3)) << sj) |
           (((o >> 2U) & UINT32_C(3)) << sk) | ((o & UINT32_C(3)) << sl);
}

uint32_t phy_simplex4_check4(const uint8_t s[256]) {
    if (s == NULL) {
        return UINT32_C(0xFFFFFFFF);
    }
    uint32_t viol = UINT32_C(0);
    for (uint32_t st = UINT32_C(0); st < UINT32_C(1048576); ++st) {
        uint32_t l = st;
        l = phy_simplex4_apply4(s, l, 6U, 7U, 8U, 9U);
        l = phy_simplex4_apply4(s, l, 3U, 4U, 5U, 9U);
        l = phy_simplex4_apply4(s, l, 1U, 2U, 5U, 8U);
        l = phy_simplex4_apply4(s, l, 0U, 2U, 4U, 7U);
        l = phy_simplex4_apply4(s, l, 0U, 1U, 3U, 6U);

        uint32_t r = st;
        r = phy_simplex4_apply4(s, r, 0U, 1U, 3U, 6U);
        r = phy_simplex4_apply4(s, r, 0U, 2U, 4U, 7U);
        r = phy_simplex4_apply4(s, r, 1U, 2U, 5U, 8U);
        r = phy_simplex4_apply4(s, r, 3U, 4U, 5U, 9U);
        r = phy_simplex4_apply4(s, r, 6U, 7U, 8U, 9U);
        if (l != r) {
            ++viol;
        }
    }
    return viol;
}

static uint8_t phy_check_4ch_coupled4(const uint8_t s[256]) {
    bool c0 = false;
    bool c1 = false;
    bool c2 = false;
    bool c3 = false;
    for (uint32_t a = 0U; a < 4U; ++a) {
        for (uint32_t b = 0U; b < 4U; ++b) {
            for (uint32_t c = 0U; c < 4U; ++c) {
                for (uint32_t d = 0U; d < 4U; ++d) {
                    uint32_t base = (a << 6U) | (b << 4U) | (c << 2U) | d;
                    uint32_t o = s[base];
                    if (a > 0U) {
                        uint32_t o0 = s[(b << 4U) | (c << 2U) | d];
                        if (((o >> 4U) & 3U) != ((o0 >> 4U) & 3U)) {
                            c0 = true;
                        }
                    }
                    if (b > 0U) {
                        uint32_t o0 = s[(a << 6U) | (c << 2U) | d];
                        if (((o >> 2U) & 3U) != ((o0 >> 2U) & 3U)) {
                            c1 = true;
                        }
                    }
                    if (c > 0U) {
                        uint32_t o0 = s[(a << 6U) | (b << 4U) | d];
                        if (((o >> 4U) & 3U) != ((o0 >> 4U) & 3U)) {
                            c2 = true;
                        }
                    }
                    if (d > 0U) {
                        uint32_t o0 = s[(a << 6U) | (b << 4U) | (c << 2U)];
                        if (((o >> 2U) & 3U) != ((o0 >> 2U) & 3U)) {
                            c3 = true;
                        }
                    }
                }
            }
        }
    }
    return (c0 && c1 && c2 && c3) ? UINT8_C(1) : UINT8_C(0);
}

void phy_analyze_simplex4_4(const uint8_t s[256], phy_simplex4_4_t *out) {
    if (s == NULL || out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    uint64_t seen[4];
    phy_vzero_u64(seen, 4U);
    for (uint32_t u = 0U; u < 256U; ++u) {
        uint8_t v = s[u];
        out->s[u] = v;
        out->sinv[v] = (uint8_t)u;
        seen[(uint32_t)v >> 6U] |= (UINT64_C(1) << ((uint32_t)v & UINT32_C(63)));
    }
    out->is_bijective =
        (seen[0] == UINT64_C(0xFFFFFFFFFFFFFFFF) && seen[1] == UINT64_C(0xFFFFFFFFFFFFFFFF) &&
         seen[2] == UINT64_C(0xFFFFFFFFFFFFFFFF) && seen[3] == UINT64_C(0xFFFFFFFFFFFFFFFF))
            ? UINT8_C(1)
            : UINT8_C(0);
    out->viol4 = phy_simplex4_check4(out->s);
    out->anf_deg = phy_anf_degree8(out->s);
    out->order4 = (out->is_bijective != 0U) ? phy_perm256_order(out->s) : UINT32_C(0);

    out->is_involutive = 1U;
    for (uint32_t u = 0U; u < 256U; ++u) {
        if (out->s[out->s[u]] != u) {
            out->is_involutive = 0U;
            out->noninvol_u = (uint8_t)u;
            break;
        }
    }

    for (uint32_t a = 0U; a < 4U; ++a) {
        for (uint32_t b = 0U; b < 4U; ++b) {
            for (uint32_t c = 0U; c < 4U; ++c) {
                for (uint32_t d = 0U; d < 4U; ++d) {
                    uint32_t nb = b ^ (a & c & UINT32_C(1));
                    uint32_t nc = c ^ (b & d & UINT32_C(2));
                    out->s_coupled[(a << 6U) | (b << 4U) | (c << 2U) | d] =
                        (uint8_t)((a << 6U) | (nb << 4U) | (nc << 2U) | d);
                }
            }
        }
    }
    out->viol4_coupled = phy_simplex4_check4(out->s_coupled);
    out->is_4coupled = phy_check_4ch_coupled4(out->s_coupled);

    phy_simplex4_2_t s2;
    phy_analyze_simplex4_2(&s2);
    phy_vcopy_u8(out->ccx02, s2.ccx02, 16U);
    out->viol2_ccx = s2.viol2_ccx02;
    out->viol2_cccx = s2.viol2_cccx1;

    uint32_t fl = UINT32_C(0);
    if (out->is_bijective != 0U) {
        fl |= PHY_S4_FLAG_BIJECTIVE;
    }
    if (out->viol4 == UINT32_C(0) && out->viol4_coupled == UINT32_C(0) && out->viol2_ccx == 0U &&
        out->viol2_cccx > 0U) {
        fl |= PHY_S4_FLAG_4SIMPLEX;
    }
    if (out->is_involutive == 0U) {
        fl |= PHY_S4_FLAG_NONINVOL;
    }
    if (out->is_4coupled != 0U) {
        fl |= PHY_S4_FLAG_4COUPLED;
    }
    if (out->anf_deg >= 3U) {
        fl |= PHY_S4_FLAG_NONLINEAR;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    for (uint32_t u = 0U; u < 256U; ++u) {
        mh = phy_fnv1a64_u64(mh, ((uint64_t)out->s[u] << 8U) | (uint64_t)out->s_coupled[u]);
    }
    mh = phy_fnv1a64_u64(mh, ((uint64_t)fl << 32U) | (uint64_t)out->order4);
    out->merkle = mh;
}

static void phy_store_arena_simplex4(phy_arena_u *arena, uint32_t idx, const phy_simplex4_4_t *s4) {
    if (arena == NULL || idx >= PHY_COUNT_SIMPLEX4_META) {
        return;
    }
    uint32_t slot = PHY_SLOT_SIMPLEX4_META + idx;
    phy_vzero_u64(arena->cells[slot].u64, 8U);
    arena->map.simplex4_meta[idx].merkle = s4->merkle;
    arena->map.simplex4_meta[idx].flags = s4->flags;
    arena->map.simplex4_meta[idx].viol4 = s4->viol4;
    arena->map.simplex4_meta[idx].order4 = s4->order4;
    arena->map.simplex4_meta[idx].viol2_ccx = s4->viol2_ccx;
    arena->map.simplex4_meta[idx].viol2_cccx = s4->viol2_cccx;
    arena->map.simplex4_meta[idx].anf_deg = s4->anf_deg;
    arena->map.simplex4_meta[idx].is_bijective = s4->is_bijective;
    arena->map.simplex4_meta[idx].is_involutive = s4->is_involutive;
    arena->map.simplex4_meta[idx].is_4coupled = s4->is_4coupled;
}

uint32_t phy_hunt_simplex4(phy_arena_u *arena, uint64_t seed, uint64_t max_trials,
                           phy_simplex4_4_t *out_best) {
    if (out_best == NULL || max_trials == UINT64_C(0)) {
        return UINT32_C(0);
    }
    uint64_t rng = (seed == UINT64_C(0)) ? UINT64_C(2026) : seed;
    uint8_t s4p[24][4];
    uint8_t fix0[6];
    phy_gen_s4_perms(s4p, fix0);

    uint64_t seen[PHY_COUNT_SIMPLEX4_META];
    phy_vzero_u64(seen, PHY_COUNT_SIMPLEX4_META);
    uint32_t found = UINT32_C(0);
    bool have_best = false;
    uint64_t best_sc = UINT64_C(0);

    uint32_t req = PHY_S4_FLAG_BIJECTIVE | PHY_S4_FLAG_4SIMPLEX | PHY_S4_FLAG_NONINVOL |
                   PHY_S4_FLAG_4COUPLED | PHY_S4_FLAG_NONLINEAR;

    phy_simplex4_4_t cand;
    for (uint64_t trial = 0U; trial < max_trials; ++trial) {
        uint8_t p_id =
            (trial == UINT64_C(0)) ? fix0[0] : fix0[searu_Bounded_u32(searu_rng_next(&rng), 6U)];
        uint8_t k =
            (trial == UINT64_C(0)) ? 1U : (uint8_t)(((searu_rng_next(&rng) & 1U) != 0U) ? 3U : 1U);
        uint8_t pi[4];
        uint8_t pi_inv[4];
        for (uint8_t x = 0U; x < 4U; ++x) {
            pi[x] = s4p[p_id][x];
            pi_inv[pi[x]] = x;
        }
        uint8_t s[256];
        for (uint32_t a = 0U; a < 4U; ++a) {
            for (uint32_t b = 0U; b < 4U; ++b) {
                for (uint32_t c = 0U; c < 4U; ++c) {
                    for (uint32_t d = 0U; d < 4U; ++d) {
                        uint32_t ob = ((uint32_t)pi[b] + (uint32_t)k * pi[a] * pi[c]) & UINT32_C(3);
                        s[(a << 6U) | (b << 4U) | (c << 2U) | d] =
                            (uint8_t)((a << 6U) | ((uint32_t)pi_inv[ob] << 4U) | (c << 2U) | d);
                    }
                }
            }
        }
        phy_analyze_simplex4_4(s, &cand);
        if ((cand.flags & req) == req && cand.s[0] == 0U) {
            bool dup = false;
            for (uint32_t i = 0U; i < found; ++i) {
                if (seen[i] == cand.merkle) {
                    dup = true;
                    break;
                }
            }
            if (!dup && found < PHY_COUNT_SIMPLEX4_META) {
                seen[found] = cand.merkle;
                phy_store_arena_simplex4(arena, found, &cand);
                ++found;
            }
            uint64_t clean_bonus = (trial == UINT64_C(0)) ? UINT64_C(1) : UINT64_C(0);
            uint64_t sc =
                (clean_bonus << 40U) | ((uint64_t)cand.anf_deg << 24U) | (uint64_t)cand.order4;
            if (!have_best || sc > best_sc) {
                have_best = true;
                best_sc = sc;
                phy_vcopy_u8((uint8_t *)out_best, (const uint8_t *)&cand, sizeof(cand));
            }
        }
    }
    if (have_best && arena != NULL) {
        arena->map.img_hdr.active_rules = found;
        arena->map.img_hdr.merkle_root = out_best->merkle;
    }
    return found;
}

uint16_t phy_verify_q6_cocycle(const uint8_t r[36], const uint8_t phi[36]) {
    if (r == NULL || phi == NULL) {
        return UINT16_C(0xFFFF);
    }
    uint8_t pa[36];
    uint8_t pb[36];
    for (uint8_t a = 0U; a < 6U; ++a) {
        for (uint8_t b = 0U; b < 6U; ++b) {
            uint32_t idx = (uint32_t)a * UINT32_C(6) + (uint32_t)b;
            uint32_t rv = (uint32_t)r[idx];
            if (rv >= UINT32_C(36)) {
                return UINT16_C(0xFFFF);
            }
        }
    }
    for (uint8_t a = 0U; a < 6U; ++a) {
        for (uint8_t b = 0U; b < 6U; ++b) {
            uint32_t idx = (uint32_t)a * UINT32_C(6) + (uint32_t)b;
            pa[idx] = a;
            pb[idx] = b;
        }
    }
    uint16_t viol = UINT16_C(0);
    for (uint8_t x = 0U; x < 6U; ++x) {
        if ((phi[(uint32_t)x * UINT32_C(6) + (uint32_t)x] & 3U) != 0U) {
            ++viol;
        }
    }
    for (uint8_t a = 0U; a < 6U; ++a) {
        for (uint8_t b = 0U; b < 6U; ++b) {
            uint32_t ab_i = (uint32_t)a * UINT32_C(6) + (uint32_t)b;
            uint32_t ab = r[ab_i];
            uint32_t a1 = pa[ab];
            uint32_t b1 = pb[ab];
            for (uint8_t c = 0U; c < 6U; ++c) {
                uint32_t bc_i = b1 * UINT32_C(6) + (uint32_t)c;
                uint32_t bc = r[bc_i];
                uint32_t b2 = pa[bc];
                uint32_t ab2_i = a1 * UINT32_C(6) + b2;

                uint32_t bc0_i = (uint32_t)b * UINT32_C(6) + (uint32_t)c;
                uint32_t bc0 = r[bc0_i];
                uint32_t b0 = pa[bc0];
                uint32_t c0 = pb[bc0];
                uint32_t ab0_i = (uint32_t)a * UINT32_C(6) + b0;
                uint32_t ab0 = r[ab0_i];
                uint32_t b00 = pb[ab0];
                uint32_t bc1_i = b00 * UINT32_C(6) + c0;

                uint32_t lhs = ((uint32_t)phi[ab_i] + (uint32_t)phi[bc_i] + (uint32_t)phi[ab2_i]) &
                               UINT32_C(3);
                uint32_t rhs =
                    ((uint32_t)phi[bc0_i] + (uint32_t)phi[ab0_i] + (uint32_t)phi[bc1_i]) &
                    UINT32_C(3);
                if (lhs != rhs) {
                    ++viol;
                }
            }
        }
    }
    return viol;
}

static void phy_q6_insert_z4_row(uint8_t piv[36][36], uint8_t has[36], uint8_t row[36]) {
    for (uint32_t col = 0U; col < 36U; ++col) {
        uint32_t rc = (uint32_t)row[col] & UINT32_C(3);
        if (rc == UINT32_C(0)) {
            continue;
        }
        if (has[col] == 0U) {
            if ((rc & UINT32_C(1)) != UINT32_C(0)) {
                for (uint32_t j = col; j < 36U; ++j) {
                    piv[col][j] = (uint8_t)(((uint32_t)row[j] * rc) & UINT32_C(3));
                }
                has[col] = 1U;
            } else {
                for (uint32_t j = col; j < 36U; ++j) {
                    piv[col][j] = (uint8_t)((uint32_t)row[j] & UINT32_C(3));
                }
                has[col] = 2U;
            }
            return;
        }
        if (has[col] == 1U) {
            for (uint32_t j = col; j < 36U; ++j) {
                uint32_t sub = (rc * (uint32_t)piv[col][j]) & UINT32_C(3);
                row[j] = (uint8_t)(((uint32_t)row[j] + UINT32_C(4) - sub) & UINT32_C(3));
            }
        } else {
            if ((rc & UINT32_C(1)) != UINT32_C(0)) {
                uint8_t tmp[36];
                for (uint32_t j = col; j < 36U; ++j) {
                    tmp[j] = (uint8_t)(((uint32_t)row[j] * rc) & UINT32_C(3));
                }
                for (uint32_t j = col; j < 36U; ++j) {
                    uint32_t sub = (UINT32_C(2) * (uint32_t)tmp[j]) & UINT32_C(3);
                    row[j] = (uint8_t)(((uint32_t)piv[col][j] + UINT32_C(4) - sub) & UINT32_C(3));
                    piv[col][j] = tmp[j];
                }
                has[col] = 1U;
            } else {
                for (uint32_t j = col; j < 36U; ++j) {
                    row[j] = (uint8_t)(((uint32_t)row[j] + UINT32_C(4) - (uint32_t)piv[col][j]) &
                                       UINT32_C(3));
                }
            }
        }
    }
}

static bool phy_q6_backsub_z4(const uint8_t piv[36][36], const uint8_t has[36], uint64_t *rng,
                              uint8_t out_phi[36]) {
    uint8_t vec[36];
    phy_vzero_u8(vec, 36U);
    for (int32_t col = 35; col >= 0; --col) {
        if (has[col] == 0U) {
            vec[col] = (uint8_t)(searu_rng_next(rng) & UINT64_C(3));
        } else {
            uint32_t sum = UINT32_C(0);
            for (uint32_t j = (uint32_t)col + UINT32_C(1); j < 36U; ++j) {
                sum += (uint32_t)piv[col][j] * (uint32_t)vec[j];
            }
            uint32_t rhs = (UINT32_C(4) - (sum & UINT32_C(3))) & UINT32_C(3);
            if (has[col] == 1U) {
                vec[col] = (uint8_t)rhs;
            } else {
                if ((rhs & UINT32_C(1)) != UINT32_C(0)) {
                    return false;
                }
                uint32_t half = (rhs >> 1U) & UINT32_C(1);
                uint32_t shift = ((searu_rng_next(rng) & UINT64_C(1)) != UINT64_C(0)) ? 2U : 0U;
                vec[col] = (uint8_t)((half + shift) & UINT32_C(3));
            }
        }
    }
    phy_vcopy_u8(out_phi, vec, 36U);
    return true;
}

static void phy_q6_eval_knots(const uint8_t r[36], const uint8_t rinv[36], const uint8_t phi[36],
                              const uint8_t pa[36], const uint8_t pb[36],
                              phy_chiral6_analysis_t *out) {
    for (uint32_t w = 0U; w < 4U; ++w) {
        out->unknot_w[w] = 0U;
        out->hopf_w[w] = 0U;
        out->tref_l_w[w] = 0U;
        out->tref_r_w[w] = 0U;
        out->fig8_w[w] = 0U;
    }
    for (uint32_t u = 0U; u < 36U; ++u) {
        /* k = 1, 2, 3 positive braids */
        uint32_t x = u;
        uint32_t s = UINT32_C(0);
        for (uint32_t step = 1U; step <= 3U; ++step) {
            s = (s + (uint32_t)phi[x]) & UINT32_C(3);
            x = r[x];
            if (x == u) {
                if (step == 1U) {
                    out->unknot_w[s] += 1U;
                } else if (step == 2U) {
                    out->hopf_w[s] += 1U;
                } else {
                    out->tref_l_w[s] += 1U;
                }
            }
        }
        /* Negative braid T(2,-3) */
        x = u;
        s = UINT32_C(0);
        for (uint32_t step = 0U; step < 3U; ++step) {
            uint32_t y = rinv[x];
            s = (s + (UINT32_C(4) - ((uint32_t)phi[y] & UINT32_C(3)))) & UINT32_C(3);
            x = y;
        }
        if (x == u) {
            out->tref_r_w[s] += 1U;
        }
    }

    for (uint8_t a = 0U; a < 6U; ++a) {
        for (uint8_t b = 0U; b < 6U; ++b) {
            for (uint8_t c = 0U; c < 6U; ++c) {
                uint32_t x = a;
                uint32_t y = b;
                uint32_t z = c;
                uint32_t s = UINT32_C(0);
                for (uint32_t rep = 0U; rep < 2U; ++rep) {
                    uint32_t u1 = x * UINT32_C(6) + y;
                    s = (s + (uint32_t)phi[u1]) & UINT32_C(3);
                    uint32_t o1 = r[u1];
                    x = pa[o1];
                    y = pb[o1];

                    uint32_t u2 = y * UINT32_C(6) + z;
                    uint32_t o2 = rinv[u2];
                    s = (s + (UINT32_C(4) - ((uint32_t)phi[o2] & UINT32_C(3)))) & UINT32_C(3);
                    y = pa[o2];
                    z = pb[o2];
                }
                if (x == (uint32_t)a && y == (uint32_t)b && z == (uint32_t)c) {
                    out->fig8_w[s] += 1U;
                }
            }
        }
    }
}

void phy_analyze_chiral6_custom(const uint8_t act[6][6], const uint8_t phi_in[36],
                                phy_chiral6_analysis_t *out) {
    if (act == NULL || out == NULL) {
        return;
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    uint8_t pa[36];
    uint8_t pb[36];
    for (uint8_t a = 0U; a < 6U; ++a) {
        for (uint8_t b = 0U; b < 6U; ++b) {
            out->act[a][b] = act[a][b];
            uint32_t idx = (uint32_t)a * UINT32_C(6) + (uint32_t)b;
            pa[idx] = a;
            pb[idx] = b;
            uint8_t res = act[b][a];
            if (res >= 6U) {
                return;
            }
            uint32_t rv = (uint32_t)b * UINT32_C(6) + (uint32_t)res;
            out->r[idx] = (uint8_t)rv;
            out->rinv[rv] = (uint8_t)idx;
        }
    }

    bool r1_ok = true;
    for (uint8_t x = 0U; x < 6U; ++x) {
        uint32_t d = (uint32_t)x * UINT32_C(6) + (uint32_t)x;
        if (out->r[d] != d) {
            r1_ok = false;
        }
    }
    uint16_t yb_viol = 0U;
    for (uint8_t a = 0U; a < 6U; ++a) {
        for (uint8_t b = 0U; b < 6U; ++b) {
            uint32_t ab = out->r[(uint32_t)a * UINT32_C(6) + (uint32_t)b];
            uint32_t a1 = pa[ab];
            uint32_t b1 = pb[ab];
            for (uint8_t c = 0U; c < 6U; ++c) {
                uint32_t bc = out->r[b1 * UINT32_C(6) + (uint32_t)c];
                uint32_t b2 = pa[bc];
                uint32_t c2 = pb[bc];
                uint32_t ab2 = out->r[a1 * UINT32_C(6) + b2];

                uint32_t bc0 = out->r[(uint32_t)b * UINT32_C(6) + (uint32_t)c];
                uint32_t b0 = pa[bc0];
                uint32_t c0 = pb[bc0];
                uint32_t ab0 = out->r[(uint32_t)a * UINT32_C(6) + b0];
                uint32_t a0 = pa[ab0];
                uint32_t b00 = pb[ab0];
                uint32_t bc1 = out->r[b00 * UINT32_C(6) + c0];
                if (pa[ab2] != a0 || pb[ab2] != pa[bc1] || c2 != pb[bc1]) {
                    ++yb_viol;
                }
            }
        }
    }

    uint8_t piv[36][36];
    uint8_t has[36];
    phy_vzero_u8(&piv[0][0], sizeof(piv));
    phy_vzero_u8(has, 36U);
    uint8_t row[36];
    for (uint8_t x = 0U; x < 6U; ++x) {
        phy_vzero_u8(row, 36U);
        row[(uint32_t)x * UINT32_C(6) + (uint32_t)x] = 1U;
        phy_q6_insert_z4_row(piv, has, row);
    }
    for (uint8_t a = 0U; a < 6U; ++a) {
        for (uint8_t b = 0U; b < 6U; ++b) {
            uint32_t ab_i = (uint32_t)a * UINT32_C(6) + (uint32_t)b;
            uint32_t ab = out->r[ab_i];
            uint32_t a1 = pa[ab];
            uint32_t b1 = pb[ab];
            for (uint8_t c = 0U; c < 6U; ++c) {
                uint32_t bc_i = b1 * UINT32_C(6) + (uint32_t)c;
                uint32_t bc = out->r[bc_i];
                uint32_t ab2_i = a1 * UINT32_C(6) + (uint32_t)pa[bc];
                uint32_t bc0_i = (uint32_t)b * UINT32_C(6) + (uint32_t)c;
                uint32_t bc0 = out->r[bc0_i];
                uint32_t ab0_i = (uint32_t)a * UINT32_C(6) + (uint32_t)pa[bc0];
                uint32_t bc1_i = (uint32_t)pb[out->r[ab0_i]] * UINT32_C(6) + (uint32_t)pb[bc0];

                phy_vzero_u8(row, 36U);
                row[ab_i] = (uint8_t)((row[ab_i] + 1U) & 3U);
                row[bc_i] = (uint8_t)((row[bc_i] + 1U) & 3U);
                row[ab2_i] = (uint8_t)((row[ab2_i] + 1U) & 3U);
                row[bc0_i] = (uint8_t)((row[bc0_i] + 3U) & 3U);
                row[ab0_i] = (uint8_t)((row[ab0_i] + 3U) & 3U);
                row[bc1_i] = (uint8_t)((row[bc1_i] + 3U) & 3U);
                phy_q6_insert_z4_row(piv, has, row);
            }
        }
    }
    uint16_t rk = 0U;
    for (uint32_t c = 0U; c < 36U; ++c) {
        if (has[c] != 0U) {
            ++rk;
        }
    }
    out->cocycle_rank = rk;

    if (phi_in != NULL) {
        for (uint32_t i = 0U; i < 36U; ++i) {
            out->phi[i] = (uint8_t)(phi_in[i] & 3U);
        }
        out->cocycle_viol = phy_verify_q6_cocycle(out->r, out->phi);
        phy_q6_eval_knots(out->r, out->rinv, out->phi, pa, pb, out);
    } else {
        uint64_t rng = UINT64_C(42);
        uint8_t cand_phi[36];
        for (uint32_t tr = 0U; tr < 128U; ++tr) {
            if (!phy_q6_backsub_z4(piv, has, &rng, cand_phi)) {
                continue;
            }
            if (phy_verify_q6_cocycle(out->r, cand_phi) != 0U) {
                continue;
            }
            phy_q6_eval_knots(out->r, out->rinv, cand_phi, pa, pb, out);
            bool diff = false;
            for (uint32_t w = 0U; w < 4U; ++w) {
                if (out->tref_l_w[w] != out->tref_r_w[w]) {
                    diff = true;
                }
            }
            if (diff) {
                if (out->tref_l_w[3] > out->tref_l_w[1]) {
                    for (uint32_t i = 0U; i < 36U; ++i) {
                        cand_phi[i] =
                            (uint8_t)((UINT32_C(4) - ((uint32_t)cand_phi[i] & UINT32_C(3))) &
                                      UINT32_C(3));
                    }
                    phy_q6_eval_knots(out->r, out->rinv, cand_phi, pa, pb, out);
                }
                phy_vcopy_u8(out->phi, cand_phi, 36U);
                out->cocycle_viol = 0U;
                break;
            }
        }
    }

    bool chiral = false;
    for (uint32_t w = 0U; w < 4U; ++w) {
        if (out->tref_l_w[w] != out->tref_r_w[w]) {
            chiral = true;
        }
    }
    out->is_chiral = chiral ? UINT8_C(1) : UINT8_C(0);

    bool non_cob = false;
    for (uint32_t u = 0U; u < 36U && !non_cob; ++u) {
        uint32_t x = u;
        uint32_t s = UINT32_C(0);
        for (uint32_t st = 0U; st < 3U; ++st) {
            s = (s + (uint32_t)out->phi[x]) & UINT32_C(3);
            x = out->r[x];
        }
        if (x == u && s != UINT32_C(0)) {
            non_cob = true;
            out->wit_a = pa[u];
            out->wit_b = pb[u];
        }
    }
    out->is_non_coboundary = (out->cocycle_viol == 0U && non_cob) ? UINT8_C(1) : UINT8_C(0);

    uint32_t fl = UINT32_C(0);
    if (yb_viol == 0U) {
        fl |= PHY_CHIRAL6_FLAG_YBE;
    }
    if (r1_ok) {
        fl |= PHY_CHIRAL6_FLAG_R1;
    }
    if (out->cocycle_viol == 0U) {
        fl |= PHY_CHIRAL6_FLAG_COCYCLE;
    }
    if (out->is_chiral != 0U) {
        fl |= PHY_CHIRAL6_FLAG_CHIRAL;
    }
    if (out->is_non_coboundary != 0U) {
        fl |= PHY_CHIRAL6_FLAG_NONCOBOUND;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    for (uint32_t i = 0U; i < 36U; ++i) {
        mh = phy_fnv1a64_u64(mh, ((uint64_t)out->r[i] << 8U) | (uint64_t)out->phi[i]);
    }
    mh = phy_fnv1a64_u64(mh, (uint64_t)fl);
    out->merkle = mh;
}

uint32_t phy_solve_chiral6(phy_arena_u *arena, phy_chiral6_analysis_t *out) {
    if (out == NULL) {
        return UINT32_C(0);
    }
    static const uint8_t q6_act[6][6] = {
        {0U, 1U, 3U, 4U, 5U, 2U}, {0U, 1U, 5U, 2U, 3U, 4U}, {5U, 3U, 2U, 0U, 4U, 1U},
        {2U, 4U, 1U, 3U, 0U, 5U}, {3U, 5U, 2U, 1U, 4U, 0U}, {4U, 2U, 0U, 3U, 1U, 5U},
    };
    phy_analyze_chiral6_custom(q6_act, NULL, out);
    uint32_t req = PHY_CHIRAL6_FLAG_YBE | PHY_CHIRAL6_FLAG_R1 | PHY_CHIRAL6_FLAG_COCYCLE |
                   PHY_CHIRAL6_FLAG_CHIRAL | PHY_CHIRAL6_FLAG_NONCOBOUND;
    if ((out->flags & req) != req) {
        return UINT32_C(0);
    }
    if (arena != NULL) {
        phy_vzero_u64(arena->cells[PHY_SLOT_CHIRAL6_META].u64, 8U);
        arena->map.chiral6_meta[0].merkle = out->merkle;
        arena->map.chiral6_meta[0].flags = out->flags;
        arena->map.chiral6_meta[0].cocycle_rank = out->cocycle_rank;
        arena->map.chiral6_meta[0].unknot_w0 = out->unknot_w[0];
        arena->map.chiral6_meta[0].hopf_w0 = out->hopf_w[0];
        arena->map.chiral6_meta[0].hopf_w2 = out->hopf_w[2];
        arena->map.chiral6_meta[0].tref_l_w1 = out->tref_l_w[1];
        arena->map.chiral6_meta[0].tref_r_w3 = out->tref_r_w[3];
        arena->map.chiral6_meta[0].fig8_w0 = out->fig8_w[0];
        arena->map.chiral6_meta[0].is_chiral = out->is_chiral;
        arena->map.chiral6_meta[0].is_non_coboundary = out->is_non_coboundary;
        arena->map.img_hdr.active_rules = UINT32_C(1);
        arena->map.img_hdr.merkle_root = out->merkle;
    }
    return UINT32_C(1);
}

void phy_verify_zi_du4(const int8_t re[4][4], const int8_t im[4][4], uint8_t *out_unitary,
                       uint8_t *out_dual, uint8_t *out_full_sup) {
    if (re == NULL || im == NULL) {
        return;
    }
    bool uni = true;
    bool full = true;
    for (uint32_t i = 0U; i < 4U; ++i) {
        for (uint32_t j = 0U; j < 4U; ++j) {
            if (re[i][j] == 0 && im[i][j] == 0) {
                full = false;
            }
        }
    }
    for (uint32_t j = 0U; j < 4U && uni; ++j) {
        for (uint32_t k = 0U; k < 4U && uni; ++k) {
            int64_t sr = INT64_C(0);
            int64_t si = INT64_C(0);
            for (uint32_t i = 0U; i < 4U; ++i) {
                int64_t ar = (int64_t)re[i][j];
                int64_t ai = (int64_t)im[i][j];
                int64_t br = (int64_t)re[i][k];
                int64_t bi = (int64_t)im[i][k];
                sr += ar * br + ai * bi;
                si += ar * bi - ai * br;
            }
            int64_t exp_r = (j == k) ? INT64_C(4) : INT64_C(0);
            if (sr != exp_r || si != INT64_C(0)) {
                uni = false;
            }
        }
    }
    bool dual = true;
    for (uint32_t b = 0U; b < 2U && dual; ++b) {
        for (uint32_t d = 0U; d < 2U && dual; ++d) {
            for (uint32_t b2 = 0U; b2 < 2U && dual; ++b2) {
                for (uint32_t d2 = 0U; d2 < 2U && dual; ++d2) {
                    uint32_t col1 = (b << 1U) | d;
                    uint32_t col2 = (b2 << 1U) | d2;
                    int64_t sr = INT64_C(0);
                    int64_t si = INT64_C(0);
                    for (uint32_t a = 0U; a < 2U; ++a) {
                        for (uint32_t c = 0U; c < 2U; ++c) {
                            uint32_t row1 = (c << 1U) | ((col1 & 1U));
                            uint32_t row2 = (c << 1U) | ((col2 & 1U));
                            uint32_t c1 = (a << 1U) | b;
                            uint32_t c2 = (a << 1U) | b2;
                            int64_t ar = (int64_t)re[row1][c1];
                            int64_t ai = (int64_t)im[row1][c1];
                            int64_t br = (int64_t)re[row2][c2];
                            int64_t bi = (int64_t)im[row2][c2];
                            sr += ar * br + ai * bi;
                            si += ar * bi - ai * br;
                        }
                    }
                    int64_t exp_r = (col1 == col2) ? INT64_C(4) : INT64_C(0);
                    if (sr != exp_r || si != INT64_C(0)) {
                        dual = false;
                    }
                }
            }
        }
    }
    if (out_unitary != NULL) {
        *out_unitary = uni ? UINT8_C(1) : UINT8_C(0);
    }
    if (out_dual != NULL) {
        *out_dual = dual ? UINT8_C(1) : UINT8_C(0);
    }
    if (out_full_sup != NULL) {
        *out_full_sup = full ? UINT8_C(1) : UINT8_C(0);
    }
}

static void phy_qscar_reduce_pow2(int64_t re[64], int64_t im[64]) {
    uint64_t mask = UINT64_C(0);
    for (uint32_t s = 0U; s < 64U; ++s) {
        mask |= phy_abs_i64(re[s]) | phy_abs_i64(im[s]);
    }
    if (mask == UINT64_C(0)) {
        return;
    }
    uint32_t tz = phy_ctz_u64(mask);
    if (tz > 0U && tz < 63U) {
        for (uint32_t s = 0U; s < 64U; ++s) {
            uint64_t rmag = phy_abs_i64(re[s]) >> tz;
            uint64_t imag = phy_abs_i64(im[s]) >> tz;
            re[s] = (re[s] < 0) ? -(int64_t)rmag : (int64_t)rmag;
            im[s] = (im[s] < 0) ? -(int64_t)imag : (int64_t)imag;
        }
    }
}

static void phy_qscar_apply_2q(int64_t re[64], int64_t im[64], const int8_t w_re[4][4],
                               const int8_t w_im[4][4], uint32_t q0, uint32_t q1) {
    int64_t nre[64];
    int64_t nim[64];
    phy_vzero_i64(nre, 64U);
    phy_vzero_i64(nim, 64U);
    for (uint32_t s = 0U; s < 64U; ++s) {
        int64_t ar = re[s];
        int64_t ai = im[s];
        if (ar == 0 && ai == 0) {
            continue;
        }
        uint32_t idx_in = (((s >> q0) & UINT32_C(1)) << 1U) | ((s >> q1) & UINT32_C(1));
        uint32_t base = s & ~((UINT32_C(1) << q0) | (UINT32_C(1) << q1));
        for (uint32_t o = 0U; o < 4U; ++o) {
            int64_t wr = (int64_t)w_re[o][idx_in];
            int64_t wi = (int64_t)w_im[o][idx_in];
            uint32_t ns = base | (((o >> 1U) & UINT32_C(1)) << q0) | ((o & UINT32_C(1)) << q1);
            nre[ns] += wr * ar - wi * ai;
            nim[ns] += wr * ai + wi * ar;
        }
    }
    for (uint32_t s = 0U; s < 64U; ++s) {
        re[s] = nre[s];
        im[s] = nim[s];
    }
}

static void phy_qscar_apply_scar4(int64_t re[64], int64_t im[64], const int8_t w_re[4][4],
                                  const int8_t w_im[4][4], uint32_t q0, uint32_t q1, uint32_t q2,
                                  uint32_t q3, bool inv_ctrl) {
    int64_t nre[64];
    int64_t nim[64];
    phy_vzero_i64(nre, 64U);
    phy_vzero_i64(nim, 64U);
    for (uint32_t s = 0U; s < 64U; ++s) {
        int64_t ar = re[s];
        int64_t ai = im[s];
        if (ar == 0 && ai == 0) {
            continue;
        }
        uint32_t b0 = (s >> q0) & UINT32_C(1);
        uint32_t b3 = (s >> q3) & UINT32_C(1);
        bool keep = inv_ctrl ? (b0 == b3) : (b0 != b3);
        if (keep) {
            nre[s] += INT64_C(2) * ar;
            nim[s] += INT64_C(2) * ai;
        } else {
            uint32_t idx_in = (((s >> q1) & UINT32_C(1)) << 1U) | ((s >> q2) & UINT32_C(1));
            uint32_t base = s & ~((UINT32_C(1) << q1) | (UINT32_C(1) << q2));
            for (uint32_t o = 0U; o < 4U; ++o) {
                int64_t wr = (int64_t)w_re[o][idx_in];
                int64_t wi = (int64_t)w_im[o][idx_in];
                uint32_t ns = base | (((o >> 1U) & UINT32_C(1)) << q1) | ((o & UINT32_C(1)) << q2);
                nre[ns] += wr * ar - wi * ai;
                nim[ns] += wr * ai + wi * ar;
            }
        }
    }
    for (uint32_t s = 0U; s < 64U; ++s) {
        re[s] = nre[s];
        im[s] = nim[s];
    }
}

static void phy_qscar_half_purity(const int64_t re[64], const int64_t im[64], uint32_t ref_s0,
                                  uint64_t *out_pur_num, uint64_t *out_pur_den,
                                  uint64_t *out_ret_num, uint64_t *out_ret_den, uint16_t *out_sup) {
    uint64_t norm2 = UINT64_C(0);
    uint16_t sup = 0U;
    for (uint32_t s = 0U; s < 64U; ++s) {
        if (re[s] != 0 || im[s] != 0) {
            ++sup;
        }
        uint64_t r = phy_abs_i64(re[s]);
        uint64_t i = phy_abs_i64(im[s]);
        norm2 += r * r + i * i;
    }
    uint64_t r0 = phy_abs_i64(re[ref_s0]);
    uint64_t i0 = phy_abs_i64(im[ref_s0]);
    uint64_t ret_n = r0 * r0 + i0 * i0;
    uint64_t ret_g = phy_stein_gcd64(ret_n, norm2);
    if (out_ret_num != NULL) {
        *out_ret_num = (ret_g == UINT64_C(0)) ? UINT64_C(0) : phy_div_u64(ret_n, ret_g, NULL);
    }
    if (out_ret_den != NULL) {
        *out_ret_den = (ret_g == UINT64_C(0)) ? UINT64_C(1) : phy_div_u64(norm2, ret_g, NULL);
    }
    if (out_sup != NULL) {
        *out_sup = sup;
    }

    uint64_t tr_rho2 = UINT64_C(0);
    for (uint32_t a1 = 0U; a1 < 8U; ++a1) {
        for (uint32_t a2 = 0U; a2 < 8U; ++a2) {
            int64_t mr = INT64_C(0);
            int64_t mi = INT64_C(0);
            for (uint32_t b = 0U; b < 8U; ++b) {
                uint32_t idx1 = (b << 3U) | a1;
                uint32_t idx2 = (b << 3U) | a2;
                mr += re[idx1] * re[idx2] + im[idx1] * im[idx2];
                mi += im[idx1] * re[idx2] - re[idx1] * im[idx2];
            }
            tr_rho2 += phy_abs_i64(mr) * phy_abs_i64(mr) + phy_abs_i64(mi) * phy_abs_i64(mi);
        }
    }
    uint64_t den = norm2 * norm2;
    uint64_t g = phy_stein_gcd64(tr_rho2, den);
    if (out_pur_num != NULL) {
        *out_pur_num = (g == UINT64_C(0)) ? UINT64_C(0) : phy_div_u64(tr_rho2, g, NULL);
    }
    if (out_pur_den != NULL) {
        *out_pur_den = (g == UINT64_C(0)) ? UINT64_C(1) : phy_div_u64(den, g, NULL);
    }
}

uint32_t phy_analyze_qscar_custom(phy_arena_u *arena, const int8_t re[4][4], const int8_t im[4][4],
                                  bool invert_scar_ctrl, phy_qscar_analysis_t *out) {
    if (re == NULL || im == NULL || out == NULL) {
        return UINT32_C(0);
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    for (uint32_t i = 0U; i < 4U; ++i) {
        for (uint32_t j = 0U; j < 4U; ++j) {
            out->w_re[i][j] = re[i][j];
            out->w_im[i][j] = im[i][j];
        }
    }
    phy_verify_zi_du4(out->w_re, out->w_im, &out->du_unitary, &out->du_dual, &out->du_full_sup);
    out->scar_unitary = out->du_unitary;

    static const uint32_t bw_pairs[6][2] = {{0U, 1U}, {2U, 3U}, {4U, 5U},
                                            {1U, 2U}, {3U, 4U}, {5U, 0U}};

    /* 1. Dual-Unitary brick-wall on L=6 from |000000> */
    int64_t st_r[64];
    int64_t st_i[64];
    phy_vzero_i64(st_r, 64U);
    phy_vzero_i64(st_i, 64U);
    st_r[0] = INT64_C(1);
    for (uint32_t t = 0U; t < 4U; ++t) {
        for (uint32_t k = 0U; k < 6U; ++k) {
            phy_qscar_apply_2q(st_r, st_i, out->w_re, out->w_im, bw_pairs[k][0], bw_pairs[k][1]);
        }
        phy_qscar_reduce_pow2(st_r, st_i);
        phy_qscar_half_purity(st_r, st_i, 0U, &out->du_purity_num[t], &out->du_purity_den[t], NULL,
                              NULL, &out->du_sup[t]);
    }

    /* 2. Local Scar + X^{\otimes 6} Period-2 DTC on Neel |010101> (0x15) */
    phy_vzero_i64(st_r, 64U);
    phy_vzero_i64(st_i, 64U);
    st_r[0x15] = INT64_C(1);
    uint16_t scar_period = 0U;
    for (uint32_t t = 0U; t < 4U; ++t) {
        for (uint32_t i = 0U; i < 6U; ++i) {
            phy_qscar_apply_scar4(st_r, st_i, out->w_re, out->w_im, i, (i + 1U) % 6U, (i + 2U) % 6U,
                                  (i + 3U) % 6U, invert_scar_ctrl);
        }
        for (uint32_t s = 0U; s < 32U; ++s) {
            uint32_t ms = s ^ UINT32_C(63);
            int64_t tr = st_r[s];
            int64_t ti = st_i[s];
            st_r[s] = st_r[ms];
            st_i[s] = st_i[ms];
            st_r[ms] = tr;
            st_i[ms] = ti;
        }
        phy_qscar_reduce_pow2(st_r, st_i);
        phy_qscar_half_purity(st_r, st_i, 0x15U, NULL, NULL, &out->neel_ret_num[t],
                              &out->neel_ret_den[t], &out->neel_sup[t]);
        if (scar_period == 0U && out->neel_sup[t] == 1U && out->neel_ret_num[t] == UINT64_C(1) &&
            out->neel_ret_den[t] == UINT64_C(1)) {
            scar_period = (uint16_t)(t + 1U);
        }
    }
    out->scar_period = scar_period;

    /* 3. Single-defect thermal trajectory from |010111> (0x17) */
    phy_vzero_i64(st_r, 64U);
    phy_vzero_i64(st_i, 64U);
    st_r[0x17] = INT64_C(1);
    for (uint32_t t = 0U; t < 4U; ++t) {
        for (uint32_t i = 0U; i < 6U; ++i) {
            phy_qscar_apply_scar4(st_r, st_i, out->w_re, out->w_im, i, (i + 1U) % 6U, (i + 2U) % 6U,
                                  (i + 3U) % 6U, invert_scar_ctrl);
        }
        for (uint32_t s = 0U; s < 32U; ++s) {
            uint32_t ms = s ^ UINT32_C(63);
            int64_t tr = st_r[s];
            int64_t ti = st_i[s];
            st_r[s] = st_r[ms];
            st_i[s] = st_i[ms];
            st_r[ms] = tr;
            st_i[ms] = ti;
        }
        phy_qscar_reduce_pow2(st_r, st_i);
        phy_qscar_half_purity(st_r, st_i, 0x17U, &out->therm_purity_num[t],
                              &out->therm_purity_den[t], NULL, NULL, &out->therm_sup[t]);
    }

    uint32_t fl = UINT32_C(0);
    if (out->du_unitary != 0U) {
        fl |= PHY_QSCAR_FLAG_UNITARY;
    }
    if (out->du_dual != 0U) {
        fl |= PHY_QSCAR_FLAG_DUAL;
    }
    if (out->du_full_sup != 0U) {
        fl |= PHY_QSCAR_FLAG_FULL_SUP;
    }
    if (out->scar_period == 2U && out->neel_sup[0] == 1U && out->neel_ret_num[0] == UINT64_C(0) &&
        out->neel_ret_num[1] == UINT64_C(1)) {
        fl |= PHY_QSCAR_FLAG_DTC_P2;
    }
    if (out->therm_sup[3] >= 32U && out->therm_purity_den[3] > UINT64_C(0) &&
        out->therm_purity_num[3] * UINT64_C(3) < out->therm_purity_den[3]) {
        fl |= PHY_QSCAR_FLAG_THERMAL;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    mh = phy_fnv1a64_u64(mh, out->du_purity_num[1] ^ (out->du_purity_den[1] << 8U));
    mh = phy_fnv1a64_u64(mh, out->therm_purity_num[3] ^ out->therm_purity_den[3]);
    mh = phy_fnv1a64_u64(mh, (uint64_t)fl);
    out->merkle = mh;

    if (arena != NULL && !invert_scar_ctrl) {
        phy_vzero_u64(arena->cells[PHY_SLOT_QSCAR_META].u64, 8U);
        arena->map.qscar_meta[0].merkle = out->merkle;
        arena->map.qscar_meta[0].du_purity_t2_num = out->du_purity_num[1];
        arena->map.qscar_meta[0].du_purity_t2_den = out->du_purity_den[1];
        arena->map.qscar_meta[0].therm_purity_t4_num = out->therm_purity_num[3];
        arena->map.qscar_meta[0].therm_purity_t4_den = out->therm_purity_den[3];
        arena->map.qscar_meta[0].flags = out->flags;
        arena->map.qscar_meta[0].scar_period = out->scar_period;
        arena->map.qscar_meta[0].therm_max_sup = out->therm_sup[3];
        arena->map.qscar_meta[0].du_unitary = out->du_unitary;
        arena->map.qscar_meta[0].du_dual = out->du_dual;
        arena->map.qscar_meta[0].du_full_sup = out->du_full_sup;
        arena->map.img_hdr.active_rules = UINT32_C(1);
        arena->map.img_hdr.merkle_root = out->merkle;
    }
    uint32_t req = PHY_QSCAR_FLAG_UNITARY | PHY_QSCAR_FLAG_DUAL | PHY_QSCAR_FLAG_FULL_SUP |
                   PHY_QSCAR_FLAG_DTC_P2 | PHY_QSCAR_FLAG_THERMAL;
    return ((out->flags & req) == req) ? UINT32_C(1) : UINT32_C(0);
}

uint32_t phy_analyze_qscar(phy_arena_u *arena, phy_qscar_analysis_t *out) {
    static const int8_t w_re[4][4] = {
        {1, 0, 1, 0},
        {1, 0, -1, 0},
        {1, 0, 1, 0},
        {1, 0, -1, 0},
    };
    static const int8_t w_im[4][4] = {
        {0, 1, 0, -1},
        {0, 1, 0, 1},
        {0, -1, 0, 1},
        {0, -1, 0, -1},
    };
    return phy_analyze_qscar_custom(arena, w_re, w_im, false, out);
}

void phy_build_dipole_toffoli(uint8_t t[64]) {
    if (t == NULL) {
        return;
    }
    for (uint8_t u = 0U; u < 64U; ++u) {
        t[u] = u;
    }
    /* (1,2,1)=25 <-> (2,0,2)=34 */
    t[25] = 34U;
    t[34] = 25U;
}

static void phy_frac_hsf_sector_stats(const uint8_t t[64], const uint8_t tinv[64], uint32_t n_sites,
                                      phy_fracton_scale_ws_t *ws, uint16_t *out_sec,
                                      uint16_t *out_comps, uint16_t *out_largest,
                                      uint16_t *out_frozen) {
    uint16_t sec_cnt[25][85];
    for (uint32_t q = 0U; q < 25U; ++q) {
        for (uint32_t d = 0U; d < 85U; ++d) {
            sec_cnt[q][d] = 0U;
        }
    }
    uint32_t total_st = UINT32_C(1) << (n_sites * UINT32_C(2));
    uint32_t best_q = 0U;
    uint32_t best_d = 0U;
    uint16_t best_sz = 0U;
    for (uint32_t st = 0U; st < total_st; ++st) {
        int64_t q0 = INT64_C(0);
        int64_t d0 = INT64_C(0);
        phy_frac_qd((uint64_t)st, n_sites, &q0, &d0);
        if (q0 >= 0 && q0 < 25 && d0 >= 0 && d0 < 85) {
            sec_cnt[q0][d0] += 1U;
            if (sec_cnt[q0][d0] > best_sz) {
                best_sz = sec_cnt[q0][d0];
                best_q = (uint32_t)q0;
                best_d = (uint32_t)d0;
            }
        }
    }
    uint32_t n_sec = 0U;
    for (uint32_t st = 0U; st < total_st; ++st) {
        int64_t q0 = INT64_C(0);
        int64_t d0 = INT64_C(0);
        phy_frac_qd((uint64_t)st, n_sites, &q0, &d0);
        if (q0 == (int64_t)best_q && d0 == (int64_t)best_d && n_sec < UINT32_C(2048)) {
            ws->sec_states[n_sec++] = st;
        }
    }
    phy_vzero_u8(ws->vis, total_st);
    uint16_t n_comps = 0U;
    uint16_t max_comp = 0U;
    uint16_t frozen = 0U;
    uint32_t st_mask = total_st - UINT32_C(1);

    for (uint32_t i = 0U; i < n_sec; ++i) {
        uint32_t root = ws->sec_states[i];
        if (ws->vis[root] != 0U) {
            continue;
        }
        uint32_t head = 0U;
        uint32_t tail = 0U;
        ws->vis[root] = 1U;
        ws->queue[tail++] = root;
        while (head < tail) {
            uint32_t cur = ws->queue[head++];
            for (uint32_t pos = 0U; pos + 2U < n_sites; ++pos) {
                uint32_t n1 = (uint32_t)(phy_frac32_apply_gate(t, (uint64_t)cur, pos) & st_mask);
                if (ws->vis[n1] == 0U && tail < UINT32_C(2048)) {
                    ws->vis[n1] = 1U;
                    ws->queue[tail++] = n1;
                }
                uint32_t n2 = (uint32_t)(phy_frac32_apply_gate(tinv, (uint64_t)cur, pos) & st_mask);
                if (ws->vis[n2] == 0U && tail < UINT32_C(2048)) {
                    ws->vis[n2] = 1U;
                    ws->queue[tail++] = n2;
                }
            }
        }
        uint16_t csz = (uint16_t)tail;
        ++n_comps;
        if (csz > max_comp) {
            max_comp = csz;
        }
        if (csz == 1U) {
            ++frozen;
        }
    }
    *out_sec = (uint16_t)n_sec;
    *out_comps = n_comps;
    *out_largest = max_comp;
    *out_frozen = frozen;
}

uint32_t phy_analyze_fracton_comp(phy_arena_u *arena, const uint8_t t_sample[64],
                                  const uint8_t t_toff[64], phy_fracton_scale_ws_t *ws,
                                  phy_fracton_comp_t *out) {
    if (t_sample == NULL || t_toff == NULL || ws == NULL || out == NULL) {
        return UINT32_C(0);
    }
    phy_vzero_u8((uint8_t *)out, sizeof(*out));
    uint8_t ts_inv[64];
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint8_t v = (uint8_t)(t_sample[u] & 63U);
        out->t_sample[u] = v;
        ts_inv[v] = u;
    }

    uint64_t seen = UINT64_C(0);
    bool invol = true;
    bool lap_ok = true;
    for (uint8_t u = 0U; u < 64U; ++u) {
        uint8_t v = (uint8_t)(t_toff[u] & 63U);
        out->t_toff[u] = v;
        out->t_toff_inv[v] = u;
        seen |= (UINT64_C(1) << v);
    }
    bool bij = (seen == UINT64_C(0xFFFFFFFFFFFFFFFF));
    for (uint8_t u = 0U; u < 64U; ++u) {
        if (out->t_toff[out->t_toff[u]] != u) {
            invol = false;
        }
        uint32_t uu = (uint32_t)u;
        int32_t a = (int32_t)((uu >> 4U) & UINT32_C(3));
        int32_t b = (int32_t)((uu >> 2U) & UINT32_C(3));
        int32_t c = (int32_t)(uu & UINT32_C(3));
        uint32_t o = (uint32_t)out->t_toff[u];
        int32_t oa = (int32_t)((o >> 4U) & UINT32_C(3));
        int32_t ob = (int32_t)((o >> 2U) & UINT32_C(3));
        int32_t oc = (int32_t)(o & UINT32_C(3));
        int32_t phi = oa - a;
        if ((ob - b) != -(2 * phi) || (oc - c) != phi) {
            lap_ok = false;
        }
    }
    out->is_involutive = (bij && invol) ? UINT8_C(1) : UINT8_C(0);
    out->laplacian_ok = lap_ok ? UINT8_C(1) : UINT8_C(0);

    bool ccx_ok = true;
    for (uint32_t ca = 0U; ca < 2U; ++ca) {
        for (uint32_t cc = 0U; cc < 2U; ++cc) {
            for (uint32_t tgt = 0U; tgt < 2U; ++tgt) {
                uint32_t in_u =
                    (ca != 0U && cc != 0U && tgt != 0U)
                        ? UINT32_C(34)
                        : ((ca << 4U) | ((tgt != 0U ? UINT32_C(0) : UINT32_C(2)) << 2U) | cc);
                uint32_t o_u = (uint32_t)out->t_toff[in_u];
                uint32_t got_tgt =
                    ((((o_u >> 2U) & UINT32_C(3)) == UINT32_C(0)) ? UINT32_C(1) : UINT32_C(0));
                uint32_t exp_tgt = tgt ^ (ca & cc);
                if (got_tgt != exp_tgt) {
                    ccx_ok = false;
                }
            }
        }
    }
    out->ccx_ok = ccx_ok ? UINT8_C(1) : UINT8_C(0);

    for (uint32_t idx = 0U; idx < 3U; ++idx) {
        phy_frac_hsf_sector_stats(out->t_sample, ts_inv, 6U + idx, ws, &out->sec[idx],
                                  &out->comps[idx], &out->largest[idx], &out->frozen[idx]);
    }
    bool strict = ((uint32_t)out->largest[1] * (uint32_t)out->sec[0] <
                   (uint32_t)out->largest[0] * (uint32_t)out->sec[1]) &&
                  ((uint32_t)out->largest[2] * (uint32_t)out->sec[1] <
                   (uint32_t)out->largest[1] * (uint32_t)out->sec[2]);
    out->scaling_strict = strict ? UINT8_C(1) : UINT8_C(0);

    uint32_t fl = UINT32_C(0);
    if (bij) {
        fl |= PHY_FCOMP_FLAG_BIJECTIVE;
    }
    if (out->is_involutive != 0U) {
        fl |= PHY_FCOMP_FLAG_INVOLUTIVE;
    }
    if (out->laplacian_ok != 0U) {
        fl |= PHY_FCOMP_FLAG_LAPLACIAN;
    }
    if (out->ccx_ok != 0U) {
        fl |= PHY_FCOMP_FLAG_CCX;
    }
    if (out->scaling_strict != 0U) {
        fl |= PHY_FCOMP_FLAG_SCALING;
    }
    out->flags = fl;

    uint64_t mh = UINT64_C(0xCBF29CE484222325);
    for (uint32_t u = 0U; u < 64U; ++u) {
        mh = phy_fnv1a64_u64(mh, ((uint64_t)out->t_toff[u] << 8U) | (uint64_t)out->t_sample[u]);
    }
    mh = phy_fnv1a64_u64(mh, ((uint64_t)fl << 32U) | (uint64_t)out->sec[2]);
    out->merkle = mh;

    if (arena != NULL) {
        phy_vzero_u64(arena->cells[PHY_SLOT_FRAC_COMP_META].u64, 8U);
        arena->map.frac_comp_meta[0].merkle = out->merkle;
        arena->map.frac_comp_meta[0].flags = out->flags;
        arena->map.frac_comp_meta[0].sec6 = out->sec[0];
        arena->map.frac_comp_meta[0].sec7 = out->sec[1];
        arena->map.frac_comp_meta[0].sec8 = out->sec[2];
        arena->map.frac_comp_meta[0].comp6 = out->comps[0];
        arena->map.frac_comp_meta[0].comp7 = out->comps[1];
        arena->map.frac_comp_meta[0].comp8 = out->comps[2];
        arena->map.frac_comp_meta[0].largest6 = out->largest[0];
        arena->map.frac_comp_meta[0].largest7 = out->largest[1];
        arena->map.frac_comp_meta[0].largest8 = out->largest[2];
        arena->map.frac_comp_meta[0].frozen6 = out->frozen[0];
        arena->map.frac_comp_meta[0].frozen7 = out->frozen[1];
        arena->map.frac_comp_meta[0].frozen8 = out->frozen[2];
        arena->map.frac_comp_meta[0].ccx_ok = out->ccx_ok;
        arena->map.frac_comp_meta[0].laplacian_ok = out->laplacian_ok;
        arena->map.frac_comp_meta[0].is_involutive = out->is_involutive;
        arena->map.frac_comp_meta[0].scaling_strict = out->scaling_strict;
        arena->map.img_hdr.active_rules = UINT32_C(1);
        arena->map.img_hdr.merkle_root = out->merkle;
    }
    uint32_t req = PHY_FCOMP_FLAG_BIJECTIVE | PHY_FCOMP_FLAG_INVOLUTIVE | PHY_FCOMP_FLAG_LAPLACIAN |
                   PHY_FCOMP_FLAG_CCX | PHY_FCOMP_FLAG_SCALING;
    return ((out->flags & req) == req) ? UINT32_C(1) : UINT32_C(0);
}
