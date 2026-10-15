#include <stddef.h>
#include <stdint.h>
#include <torsion_constants.h>
#if 0
#elif 8*DIGIT_LEN == 16
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7d;
const uint64_t TORSION_ODD_PRIMES[20] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f, 0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_ODD_POWERS[20] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[7] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f};
const size_t TORSION_PLUS_ODD_POWERS[7] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_MINUS_ODD_PRIMES[13] = {0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_MINUS_ODD_POWERS[13] = {0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_C1_PRIMES[7] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f};
const size_t TORSION_C1_POWERS[7] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_C2_PRIMES[13] = {0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_C2_POWERS[13] = {0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xdfff,0x825b,0xc160,0x88ed,0x869c,0xcb33,0x5ff4,0xfccb,0x37f7}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x507,0x7c81,0x2209,0x64bf,0x3960,0x8543,0x783d,0x2343,0x9fb2,0x5145,0x5546,0xfce2,0x243d,0x7363,0x4691,0x4a}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[20] = {{{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x61}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x107}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1df}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3d1}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5a7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x61f}}}, {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xe821,0x291f,0xb452,0xa8b8}}}, {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0xf5e1,0x5}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x79}}}, {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0xaa5d,0x5}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xe9}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xef}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x119}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x199}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1d3}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x22d}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3cb}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x427}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x761}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xee47,0xb0bb,0x36dc,0x2a1}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 12, ._mp_d = (mp_limb_t[]) {0xd341,0x26f3,0xacc,0xcab0,0xbd0b,0x22c3,0xc549,0xf47f,0xb7f3,0x31d5,0x9450,0x1c3e}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x2000}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x2000}}};
const ibz_t TORSION_C1 = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xee47,0xb0bb,0x36dc,0x2a1}}};
const ibz_t TORSION_C2 = {{._mp_alloc = 0, ._mp_size = 12, ._mp_d = (mp_limb_t[]) {0xd341,0x26f3,0xacc,0xcab0,0xbd0b,0x22c3,0xc549,0xf47f,0xb7f3,0x31d5,0x9450,0x1c3e}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0xaa9,0x34b3,0x22b8,0x437c,0xaa}}};
const ibz_t TORSION_PLUS_2_C1 = {{._mp_alloc = 0, ._mp_size = 12, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0xe000,0x7dc8,0x9617,0x26db,0x54}}};
const ibz_t TORSION_PRIMES_PRODUCT = {{._mp_alloc = 0, ._mp_size = 14, ._mp_d = (mp_limb_t[]) {0x53c6,0x119,0x829,0xb51b,0xb79f,0xfec0,0x5586,0xa8ef,0x85c9,0xb34c,0x90b,0x21c4,0x2ff1,0x4eb2}}};
#elif 8*DIGIT_LEN == 32
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7d;
const uint64_t TORSION_ODD_PRIMES[20] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f, 0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_ODD_POWERS[20] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[7] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f};
const size_t TORSION_PLUS_ODD_POWERS[7] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_MINUS_ODD_PRIMES[13] = {0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_MINUS_ODD_POWERS[13] = {0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_C1_PRIMES[7] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f};
const size_t TORSION_C1_POWERS[7] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_C2_PRIMES[13] = {0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_C2_POWERS[13] = {0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0xffffffff,0xffffffff,0xffffffff,0xdfffffff,0xc160825b,0x869c88ed,0x5ff4cb33,0x37f7fccb}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x7c810507,0x64bf2209,0x85433960,0x2343783d,0x51459fb2,0xfce25546,0x7363243d,0x4a4691}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[20] = {{{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x61}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x107}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1df}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3d1}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5a7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x61f}}}, {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x291fe821,0xa8b8b452}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5f5e1}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x79}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5aa5d}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xe9}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xef}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x119}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x199}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1d3}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x22d}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3cb}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x427}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x761}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0xb0bbee47,0x2a136dc}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 6, ._mp_d = (mp_limb_t[]) {0x26f3d341,0xcab00acc,0x22c3bd0b,0xf47fc549,0x31d5b7f3,0x1c3e9450}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x20000000}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x20000000}}};
const ibz_t TORSION_C1 = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0xb0bbee47,0x2a136dc}}};
const ibz_t TORSION_C2 = {{._mp_alloc = 0, ._mp_size = 6, ._mp_d = (mp_limb_t[]) {0x26f3d341,0xcab00acc,0x22c3bd0b,0xf47fc549,0x31d5b7f3,0x1c3e9450}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x34b30aa9,0x437c22b8,0xaa}}};
const ibz_t TORSION_PLUS_2_C1 = {{._mp_alloc = 0, ._mp_size = 6, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0xe0000000,0x96177dc8,0x5426db}}};
const ibz_t TORSION_PRIMES_PRODUCT = {{._mp_alloc = 0, ._mp_size = 7, ._mp_d = (mp_limb_t[]) {0x11953c6,0xb51b0829,0xfec0b79f,0xa8ef5586,0xb34c85c9,0x21c4090b,0x4eb22ff1}}};
#elif 8*DIGIT_LEN == 64
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7d;
const uint64_t TORSION_ODD_PRIMES[20] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f, 0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_ODD_POWERS[20] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[7] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f};
const size_t TORSION_PLUS_ODD_POWERS[7] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_MINUS_ODD_PRIMES[13] = {0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_MINUS_ODD_POWERS[13] = {0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_C1_PRIMES[7] = {0x7, 0x61, 0x107, 0x1df, 0x3d1, 0x5a7, 0x61f};
const size_t TORSION_C1_POWERS[7] = {0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const uint64_t TORSION_C2_PRIMES[13] = {0x3, 0x5, 0xb, 0xd, 0xe9, 0xef, 0x119, 0x199, 0x1d3, 0x22d, 0x3cb, 0x427, 0x761};
const size_t TORSION_C2_POWERS[13] = {0x28, 0x8, 0x2, 0x5, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xffffffffffffffff,0xdfffffffffffffff,0x869c88edc160825b,0x37f7fccb5ff4cb33}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x64bf22097c810507,0x2343783d85433960,0xfce2554651459fb2,0x4a46917363243d}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[20] = {{{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x61}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x107}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1df}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3d1}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5a7}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x61f}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xa8b8b452291fe821}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5f5e1}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x79}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5aa5d}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xe9}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xef}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x119}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x199}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1d3}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x22d}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3cb}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x427}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x761}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x2a136dcb0bbee47}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xcab00acc26f3d341,0xf47fc54922c3bd0b,0x1c3e945031d5b7f3}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x0,0x2000000000000000}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x0,0x2000000000000000}}};
const ibz_t TORSION_C1 = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x2a136dcb0bbee47}}};
const ibz_t TORSION_C2 = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xcab00acc26f3d341,0xf47fc54922c3bd0b,0x1c3e945031d5b7f3}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x437c22b834b30aa9,0xaa}}};
const ibz_t TORSION_PLUS_2_C1 = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x0,0xe000000000000000,0x5426db96177dc8}}};
const ibz_t TORSION_PRIMES_PRODUCT = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xb51b0829011953c6,0xa8ef5586fec0b79f,0x21c4090bb34c85c9,0x4eb22ff1}}};
#endif
