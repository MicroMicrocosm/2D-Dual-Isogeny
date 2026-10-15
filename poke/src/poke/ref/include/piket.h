#include <ec.h>

/** @brief PIKE secret key
 *
 * @typedef piket_sk_t
 *
 * @struct piket_sk_t
 *
 */
typedef struct piket_sk_t {
    digit_t deg[NWORDS_ORDER];
    digit_t alpha[NWORDS_ORDER];
    digit_t beta[NWORDS_ORDER];
    digit_t delta[NWORDS_ORDER];
} piket_sk_t;

/** @brief PIKE public key
 *
 * @typedef piket_pk_t
 *
 * @struct piket_pk_t
 *
 */
typedef struct piket_pk_t {
    ec_curve_t EA;
    ec_basis_t PQ2;
    ec_basis_t PQ3;
    ec_point_t imX;
} piket_pk_t;

/** @brief PIKE ciphertext
 *
 * @typedef piket_ct_t
 *
 * @struct piket_ct_t
 *
 */
typedef struct piket_ct_t {
    ec_curve_t EB;
    ec_curve_t EAB;
    ec_basis_t PQ2_B;
    ec_basis_t PQ2_AB;
    ec_point_t imY_B;
    ec_point_t imX_AB;
    uint8_t ct[32];
} piket_ct_t;

int decaps(unsigned char *key, piket_ct_t *ct, const piket_pk_t *pk, const piket_sk_t *sk, unsigned char *dummy_m);
int encaps(unsigned char *key, piket_ct_t *ct, const piket_pk_t *pk);
int ct_encode(unsigned char *encoded_ct, piket_ct_t *ct);
int ct_decode(piket_ct_t *ct, const unsigned char *encoded_ct);
int keygen(piket_sk_t *sk, piket_pk_t *pk);
int encrypt(piket_ct_t *ct, const piket_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len);
int decrypt(unsigned char *m, size_t *m_len, const piket_ct_t *ct, const piket_sk_t *sk);
