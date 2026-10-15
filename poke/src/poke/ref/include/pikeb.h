#include "ec.h"

/** @brief PIKE secret key
 *
 * @typedef pikeb_sk_t
 * @struct pikeb_sk_t
 */
typedef struct pikeb_sk_t {
    digit_t deg[NWORDS_ORDER];      // q
    digit_t alpha[NWORDS_ORDER];    // mask P in E[2^a]
    digit_t beta[NWORDS_ORDER];     // mask Q in E[2^a]
    digit_t delta[NWORDS_ORDER];    // mask P in E[D]
} pikeb_sk_t;

/** @brief PIKE public key
 *
 * @typedef pikeb_pk_t
 * @struct pikeb_pk_t
 */
typedef struct pikeb_pk_t {
    ec_curve_t EA;
    ec_basis_t PQ2;
    ec_basis_t PQC1;
    ec_basis_t PQC2;
    ec_point_t imX;
} pikeb_pk_t;

/** @brief PIKE ciphertext
 *
 * @typedef pikeb_ct_t
 * @struct pikeb_ct_t
 */
typedef struct pikeb_ct_t {
    ec_curve_t EB;
    ec_curve_t EAB;
    ec_basis_t PQ2_B;
    ec_basis_t PQ2_AB;
    ec_point_t imY_B;
    ec_point_t imX_AB;
    uint8_t ct[32];
} pikeb_ct_t;

int decaps(unsigned char *key, pikeb_ct_t *ct, const pikeb_pk_t *pk, const pikeb_sk_t *sk, unsigned char *dummy_m);
int encaps(unsigned char *key, pikeb_ct_t *ct, const pikeb_pk_t *pk);
int ct_encode(unsigned char *encoded_ct, pikeb_ct_t *ct);
int ct_decode(pikeb_ct_t *ct, const unsigned char *encoded_ct);
int keygen(pikeb_sk_t *sk, pikeb_pk_t *pk);
int encrypt(pikeb_ct_t *ct, const pikeb_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len);
int decrypt(unsigned char *m, size_t *m_len, const pikeb_ct_t *ct, const pikeb_sk_t *sk);