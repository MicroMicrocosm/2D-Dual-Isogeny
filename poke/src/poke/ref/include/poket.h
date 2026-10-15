#include <ec.h>

/** @brief POKE secret key
 *
 * @typedef poket_sk_t
 *
 * @struct poket_sk_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct poket_sk_t {
    digit_t deg[NWORDS_ORDER];
    digit_t alpha[NWORDS_ORDER];
    digit_t beta[NWORDS_ORDER];
    digit_t delta[NWORDS_ORDER];
} poket_sk_t;

/** @brief POKE public key
 *
 * @typedef poket_pk_t
 *
 * @struct poket_pk_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct poket_pk_t {
    ec_curve_t EA;
    ec_basis_t PQ2;
    ec_basis_t PQ3;
    ec_point_t imX;
} poket_pk_t;

/** @brief POKE ciphertext
 *
 * @typedef poket_ct_t
 *
 * @struct poket_ct_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct poket_ct_t {
    ec_curve_t EB;
    ec_curve_t EAB;
    ec_basis_t PQ2_B;
    ec_basis_t PQ2_AB;
    ec_point_t imX_B;
    uint8_t ct[32];
} poket_ct_t;

int decaps(unsigned char *key, poket_ct_t *ct, const poket_pk_t *pk, const poket_sk_t *sk, unsigned char *dummy_m);
int encaps(unsigned char *key, poket_ct_t *ct, const poket_pk_t *pk);
int ct_encode(unsigned char *encoded_ct, poket_ct_t *ct);
int ct_decode(poket_ct_t *ct, const unsigned char *encoded_ct);
int keygen(poket_sk_t *sk, poket_pk_t *pk);
int encrypt(poket_ct_t *ct, const poket_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len);
int decrypt(unsigned char *m, size_t *m_len, const poket_ct_t *ct, const poket_sk_t *sk);