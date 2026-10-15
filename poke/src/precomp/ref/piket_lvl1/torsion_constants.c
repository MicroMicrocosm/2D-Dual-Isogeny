#include <stddef.h>
#include <stdint.h>
#include <torsion_constants.h>
#if 0
#elif 8*DIGIT_LEN == 16
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7e;
const uint64_t TORSION_ODD_PRIMES[1] = {0x3};
const uint64_t TORSION_ODD_POWERS[1] = {0x9e};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x9e};
const size_t DEGREE_COMMITMENT_POWERS[2] = {0x9e, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 24, ._mp_d = (mp_limb_t[]) {0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xbfff,0xd59,0x8d84,0xace8,0x6547,0x68b4,0xd623,0x433e,0x475e,0x8eb9,0xb91f,0xce38,0xf981,0x3f59,0x7556,0x20c1,0x6fdf}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xa2b9,0x2c28,0xfcd8,0xeb5c,0xce96,0x9fdf,0x3b5f,0x4cff,0xbc02,0xdf4d,0x3a7c,0xb07b,0x2d32,0x756d,0xe9c5,0x55d}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0x3c77,0xb241,0xf954,0xdd59,0xf008,0xe2d1,0xc3ab,0x6c26,0xa}}};
const ibz_t TORSION_D_COFACTOR = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x4572,0x5851,0xf9b0,0xd6b9,0x9d2d,0x3fbf,0x76bf,0x99fe,0x7804,0xbe9b,0x74f9,0x60f6,0x5a65,0xeada,0xd38a,0xabb}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[2] = {{{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xa2b9,0x2c28,0xfcd8,0xeb5c,0xce96,0x9fdf,0x3b5f,0x4cff,0xbc02,0xdf4d,0x3a7c,0xb07b,0x2d32,0x756d,0xe9c5,0x55d}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xa2b9,0x2c28,0xfcd8,0xeb5c,0xce96,0x9fdf,0x3b5f,0x4cff,0xbc02,0xdf4d,0x3a7c,0xb07b,0x2d32,0x756d,0xe9c5,0x55d}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x4000}}};
const ibz_t TORSION_MINUS_3POWER = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xa2b9,0x2c28,0xfcd8,0xeb5c,0xce96,0x9fdf,0x3b5f,0x4cff,0xbc02,0xdf4d,0x3a7c,0xb07b,0x2d32,0x756d,0xe9c5,0x55d}}};
const ibz_t TORSION_PLUS_CPOWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 24, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x4000,0x28ae,0xb0a,0x3f36,0xbad7,0xf3a5,0xe7f7,0xced7,0x933f,0x6f00,0x37d3,0xce9f,0xac1e,0x4b4c,0x5d5b,0x7a71,0x157}}};
const ibz_t TORSION_PLUS_3CPOWER = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xe82b,0x8479,0xf688,0xc216,0x6bc4,0xdf9f,0xb21e,0xe6fd,0x3406,0x9de9,0xaf76,0x1171,0x8798,0x6047,0xbd50,0x1019}}};
const ibz_t TORSION_PLUS_2CPOWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0xc000}}};
const ibz_t TORSION_PLUS_23CPOWER = {{._mp_alloc = 0, ._mp_size = 24, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0xc000,0x7a0a,0x211e,0xbda2,0x3085,0xdaf1,0xb7e7,0x6c87,0xb9bf,0x4d01,0xa77a,0x6bdd,0x45c,0xe1e6,0x1811,0x6f54,0x406}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xa2b9,0x2c28,0xfcd8,0xeb5c,0xce96,0x9fdf,0x3b5f,0x4cff,0xbc02,0xdf4d,0x3a7c,0xb07b,0x2d32,0x756d,0xe9c5,0x55d}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xa2b9,0x2c28,0xfcd8,0xeb5c,0xce96,0x9fdf,0x3b5f,0x4cff,0xbc02,0xdf4d,0x3a7c,0xb07b,0x2d32,0x756d,0xe9c5,0x55d}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x4000}}};
#elif 8*DIGIT_LEN == 32
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7e;
const uint64_t TORSION_ODD_PRIMES[1] = {0x3};
const uint64_t TORSION_ODD_POWERS[1] = {0x9e};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x9e};
const size_t DEGREE_COMMITMENT_POWERS[2] = {0x9e, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 12, ._mp_d = (mp_limb_t[]) {0xffffffff,0xffffffff,0xffffffff,0xbfffffff,0x8d840d59,0x6547ace8,0xd62368b4,0x475e433e,0xb91f8eb9,0xf981ce38,0x75563f59,0x6fdf20c1}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x2c28a2b9,0xeb5cfcd8,0x9fdfce96,0x4cff3b5f,0xdf4dbc02,0xb07b3a7c,0x756d2d32,0x55de9c5}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0xb2413c77,0xdd59f954,0xe2d1f008,0x6c26c3ab,0xa}}};
const ibz_t TORSION_D_COFACTOR = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x58514572,0xd6b9f9b0,0x3fbf9d2d,0x99fe76bf,0xbe9b7804,0x60f674f9,0xeada5a65,0xabbd38a}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[2] = {{{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x2c28a2b9,0xeb5cfcd8,0x9fdfce96,0x4cff3b5f,0xdf4dbc02,0xb07b3a7c,0x756d2d32,0x55de9c5}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x2c28a2b9,0xeb5cfcd8,0x9fdfce96,0x4cff3b5f,0xdf4dbc02,0xb07b3a7c,0x756d2d32,0x55de9c5}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x40000000}}};
const ibz_t TORSION_MINUS_3POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x2c28a2b9,0xeb5cfcd8,0x9fdfce96,0x4cff3b5f,0xdf4dbc02,0xb07b3a7c,0x756d2d32,0x55de9c5}}};
const ibz_t TORSION_PLUS_CPOWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 12, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x40000000,0xb0a28ae,0xbad73f36,0xe7f7f3a5,0x933fced7,0x37d36f00,0xac1ece9f,0x5d5b4b4c,0x1577a71}}};
const ibz_t TORSION_PLUS_3CPOWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x8479e82b,0xc216f688,0xdf9f6bc4,0xe6fdb21e,0x9de93406,0x1171af76,0x60478798,0x1019bd50}}};
const ibz_t TORSION_PLUS_2CPOWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0xc0000000}}};
const ibz_t TORSION_PLUS_23CPOWER = {{._mp_alloc = 0, ._mp_size = 12, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0xc0000000,0x211e7a0a,0x3085bda2,0xb7e7daf1,0xb9bf6c87,0xa77a4d01,0x45c6bdd,0x1811e1e6,0x4066f54}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x2c28a2b9,0xeb5cfcd8,0x9fdfce96,0x4cff3b5f,0xdf4dbc02,0xb07b3a7c,0x756d2d32,0x55de9c5}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x2c28a2b9,0xeb5cfcd8,0x9fdfce96,0x4cff3b5f,0xdf4dbc02,0xb07b3a7c,0x756d2d32,0x55de9c5}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x40000000}}};
#elif 8*DIGIT_LEN == 64
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7e;
const uint64_t TORSION_ODD_PRIMES[1] = {0x3};
const uint64_t TORSION_ODD_POWERS[1] = {0x9e};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x9e};
const size_t DEGREE_COMMITMENT_POWERS[2] = {0x9e, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 6, ._mp_d = (mp_limb_t[]) {0xffffffffffffffff,0xbfffffffffffffff,0x6547ace88d840d59,0x475e433ed62368b4,0xf981ce38b91f8eb9,0x6fdf20c175563f59}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xeb5cfcd82c28a2b9,0x4cff3b5f9fdfce96,0xb07b3a7cdf4dbc02,0x55de9c5756d2d32}}};
const ibz_t TORSION_D = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0xdd59f954b2413c77,0x6c26c3abe2d1f008,0xa}}};
const ibz_t TORSION_D_COFACTOR = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xd6b9f9b058514572,0x99fe76bf3fbf9d2d,0x60f674f9be9b7804,0xabbd38aeada5a65}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[2] = {{{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xeb5cfcd82c28a2b9,0x4cff3b5f9fdfce96,0xb07b3a7cdf4dbc02,0x55de9c5756d2d32}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xeb5cfcd82c28a2b9,0x4cff3b5f9fdfce96,0xb07b3a7cdf4dbc02,0x55de9c5756d2d32}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x0,0x4000000000000000}}};
const ibz_t TORSION_MINUS_3POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xeb5cfcd82c28a2b9,0x4cff3b5f9fdfce96,0xb07b3a7cdf4dbc02,0x55de9c5756d2d32}}};
const ibz_t TORSION_PLUS_CPOWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x3}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 6, ._mp_d = (mp_limb_t[]) {0x0,0x4000000000000000,0xbad73f360b0a28ae,0x933fced7e7f7f3a5,0xac1ece9f37d36f00,0x1577a715d5b4b4c}}};
const ibz_t TORSION_PLUS_3CPOWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xc216f6888479e82b,0xe6fdb21edf9f6bc4,0x1171af769de93406,0x1019bd5060478798}}};
const ibz_t TORSION_PLUS_2CPOWER = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x0,0xc000000000000000}}};
const ibz_t TORSION_PLUS_23CPOWER = {{._mp_alloc = 0, ._mp_size = 6, ._mp_d = (mp_limb_t[]) {0x0,0xc000000000000000,0x3085bda2211e7a0a,0xb9bf6c87b7e7daf1,0x45c6bdda77a4d01,0x4066f541811e1e6}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xeb5cfcd82c28a2b9,0x4cff3b5f9fdfce96,0xb07b3a7cdf4dbc02,0x55de9c5756d2d32}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xeb5cfcd82c28a2b9,0x4cff3b5f9fdfce96,0xb07b3a7cdf4dbc02,0x55de9c5756d2d32}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x0,0x4000000000000000}}};
#endif
