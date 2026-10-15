#include <ec.h>

/** @brief POKE secret key
 *
 * @typedef pokeb_sk_t
 *
 * @struct pokeb_sk_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct pokeb_sk_t {
    digit_t deg[NWORDS_ORDER];
    digit_t alpha[NWORDS_ORDER];
    digit_t beta[NWORDS_ORDER];
    digit_t delta[NWORDS_ORDER];
} pokeb_sk_t;

/** @brief POKE public key
 *
 * @typedef pokeb_pk_t
 *
 * @struct pokeb_pk_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct pokeb_pk_t {
    ec_curve_t EA;
    ec_basis_t PQ2;
    ec_basis_t PQC1;
    ec_basis_t PQC2;
    ec_point_t imX;
} pokeb_pk_t;

/** @brief POKE ciphertext
 *
 * @typedef pokeb_ct_t
 *
 * @struct pokeb_ct_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct pokeb_ct_t {
    ec_curve_t EB;
    ec_curve_t EAB;
    ec_basis_t PQ2_B;
    ec_basis_t PQ2_AB;
    ec_point_t imX_B;
    uint8_t ct[32];
} pokeb_ct_t;

int decaps(unsigned char *key, pokeb_ct_t *ct, const pokeb_pk_t *pk, const pokeb_sk_t *sk, unsigned char *dummy_m);
int encaps(unsigned char *key, pokeb_ct_t *ct, const pokeb_pk_t *pk);
int ct_encode(unsigned char *encoded_ct, pokeb_ct_t *ct);
int ct_decode(pokeb_ct_t *ct, const unsigned char *encoded_ct);
int keygen(pokeb_sk_t *sk, pokeb_pk_t *pk);
int encrypt(pokeb_ct_t *ct, const pokeb_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len);
int decrypt(unsigned char *m, size_t *m_len, const pokeb_ct_t *ct, const pokeb_sk_t *sk);
