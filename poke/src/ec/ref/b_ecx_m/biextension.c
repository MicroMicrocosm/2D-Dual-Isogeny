#include "curve_extras.h"
#include <ec_params.h>
#include <assert.h>
#include "const_PIKE.h"

/*
 * We implement the biextension arithmetic by using the cubical torsor representation. For now only
 * implement the 2^e-ladder.
 *
 * Warning: both cubicalDBL and cubicalADD are off by a factor x4 with respect
 * to the cubical arithmetic.
 * Since this factor is the same, this means that the biextension
 * arithmetic is correct, so the pairings are ok (they only rely on the
 * biextension arithmetic).
 * In the special case where Q=P (self pairings), we use the cubical ladder
 * rather than the biextension ladder because this is faster. In that case,
 * when we do a ladder we are off by a factor 4^m, m the number of bits.
 * This factor thus disappear in the Weil pairing since we take a quotient,
 * and also in the Tate pairing due to the final exponentiation; so
 * everything is ok too.
 * (Note that when the curves are supersingular as in our case, the Tate
 * self pairing is always trivial anyway because the Galois structure of the
 * isogeneous curves are all the same, so the étale torsor representing the
 * Tate pairing has to be trivial).
 */

/* return the *normalised* point (A+2C)/4C */
void
A24_from_AC(ec_point_t *A24, ec_point_t const *AC)
{
    fp2_add(&A24->z, &AC->z, &AC->z);
    fp2_add(&A24->x, &AC->x, &A24->z);
    fp2_add(&A24->z, &A24->z, &A24->z); //(A+2C: 4C)
    ec_normalize_point(A24);
}

// this is exactly like xDBLv2, except we use the fact that P is normalised
// to gain a multiplication
// Warning: for now we need to assume that A24 is normalised, ie C24=1.
// (maybe add an assert?)
void
cubicalDBL(ec_point_t *Q, ec_point_t const *P, ec_point_t const *A24)
{
    // A24 = (A+2C:4C)
    fp2_t t0, t1, t2;

    assert(fp2_is_one(&A24->z));
    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    // fp2_mul(&t1, &t1, &A24->z);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_mul(&t0, &t2, &A24->x);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

// this would be exactly like xADD if PQ was 'antinormalised' as (1,z)
void
cubicalADD(ec_point_t *R, ec_point_t const *P, ec_point_t const *Q, fp2_t const *ixPQ)
{
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_add(&t2, &Q->x, &Q->z);
    fp2_sub(&t3, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t3);
    fp2_mul(&t1, &t1, &t2);
    fp2_add(&t2, &t0, &t1);
    fp2_sub(&t3, &t0, &t1);
    fp2_sqr(&R->z, &t3);
    fp2_sqr(&t2, &t2);
    fp2_mul(&R->x, ixPQ, &t2);
}

// given cubical reps of P+Q, Q, P, return P+2Q, 2Q
void
biextDBL(ec_point_t *PQQ,
         ec_point_t *QQ,
         ec_point_t const *PQ,
         ec_point_t const *Q,
         fp2_t const *ixP,
         ec_point_t const *A24)
{
    cubicalADD(PQQ, PQ, Q, ixP);
    cubicalDBL(QQ, Q, A24);
}

// iterative biextension doubling
void
biext_ladder_2e(uint64_t e,
                ec_point_t *PnQ,
                ec_point_t *nQ,
                ec_point_t const *PQ,
                ec_point_t const *Q,
                fp2_t const *ixP,
                ec_point_t const *A24)
{
    copy_point(PnQ, PQ);
    copy_point(nQ, Q);
    for (uint64_t i = 0; i < e; i++) {
        biextDBL(PnQ, nQ, PnQ, nQ, ixP, A24);
    }
}

// iterative biextension
void
biext_ladder(const digit_t *n,
            const int n_bitlen,
            ec_point_t *QnP,
            ec_point_t *nP,
            ec_point_t const *P,
            ec_point_t const *Q,
            ec_point_t const *PQ,
            fp2_t const *ixP,
            fp2_t const *ixQ,
            fp2_t const *ixPQ,
            ec_point_t const *A24)
{
    ec_point_t S1, R;
    copy_point(QnP, Q);
    ec_point_init(nP);
    ec_point_init(&R);
    copy_point(&S1, P);
    for (int i = n_bitlen / RADIX; i >= 0; i--) {
        digit_t t = 1ULL << 63;
        for (int j = RADIX - 1; j >= 0; j--) {
            if(n_bitlen >=  i * RADIX + j + 1)
            {
                cubicalADD(&R, nP, &S1, ixP);
                if((t & n[i]) == 0)
                {
                    biextDBL(QnP, nP, QnP, nP, ixQ, A24);
                    copy_point(&S1, &R);
                }
                else
                {
                    biextDBL(QnP, &S1, QnP, &S1, ixPQ, A24);
                    copy_point(nP, &R);
                }
            }
            t >>= 1;
        };
    };
    
}

// compute the monodromy ratio of cubical points [(P+nQ)/P] / [(nQ)/0]
void
ratio(fp2_t *r, ec_point_t const *PnQ, ec_point_t const *nQ, ec_point_t const *P)
{
    // Sanity tests
    assert(ec_is_zero(nQ));
    assert(is_point_equal(PnQ, P));

    fp2_mul(r, &nQ->x, &P->z);
    fp2_inv(r);
    fp2_mul(r, r, &PnQ->z);
}

// Compute the ratio X/Z above as a (X:Z) point to avoid a division
void
point_ratio(ec_point_t *R, ec_point_t const *PnQ, ec_point_t const *nQ, ec_point_t const *P)
{
    // Sanity tests
    assert(ec_is_zero(nQ));
    assert(is_point_equal(PnQ, P));

    fp2_mul(&R->x, &nQ->x, &P->x);
    fp2_copy(&R->z, &PnQ->x);
}

// (X(P):Z(P))->x(P)
void
x_coord(fp2_t *r, ec_point_t const *P)
{
    fp2_copy(r, &P->z);
    fp2_inv(r);
    fp2_mul(r, r, &P->x);
}

// compute the cubical translation of P by a point of 2-torsion T
void
translate(ec_point_t *P, ec_point_t const *T)
{
    fp2_t t0, t1, t2;
    if (fp2_is_zero(&T->z)) {
        // do nothing
    } else if (fp2_is_zero(&T->x)) {
        fp2_copy(&t0, &P->x);
        fp2_copy(&P->x, &P->z);
        fp2_copy(&P->z, &t0);
    } else {
        fp2_mul(&t0, &T->x, &P->x);
        fp2_mul(&t1, &T->z, &P->z);
        fp2_sub(&t2, &t0, &t1);
        fp2_mul(&t0, &T->z, &P->x);
        fp2_mul(&t1, &T->x, &P->z);
        fp2_sub(&P->z, &t0, &t1);
        fp2_copy(&P->x, &t2);
    }
}

// Compute the monodromy P+2^e Q (in level 1)
// The suffix _i means that we are given 1/x(P) as parameter.
// Warning: to get meaningful result when using the monodromy to compute
// pairings, we need P, Q, PQ, A24 to be normalised
// (this is not strictly necessary, but care need to be taken when they are not normalised. Only
// handle the normalised case for now)
void
monodromy_i(ec_point_t *r,
            uint64_t e,
            ec_point_t const *PQ,
            ec_point_t const *Q,
            ec_point_t const *P,
            fp2_t const *ixP,
            ec_point_t const *A24)
{
    ec_point_t PnQ, nQ;
    biext_ladder_2e(e - 1, &PnQ, &nQ, PQ, Q, ixP, A24);
    translate(&PnQ, &nQ);
    translate(&nQ, &nQ);
    point_ratio(r, &PnQ, &nQ, P);
}

// void
// monodromy(ec_point_t *r,
//           uint64_t e,
//           ec_point_t const *PQ,
//           ec_point_t const *Q,
//           ec_point_t const *P,
//           ec_point_t const *A24)
// {
//     fp2_t ixP;
//     fp2_copy(&ixP, &P->x);
//     fp2_inv(&ixP);
//     monodromy_i(r, e, PQ, Q, P, &ixP, A24);
// }

// This version computes the monodromy with respect to the biextension
// associated to 2(0_E), so the square of the monodromy above
// (not used)
void
monodromy2(fp2_t *r,
           uint64_t e,
           ec_point_t const *PQ,
           ec_point_t const *Q,
           ec_point_t const *P,
           ec_point_t const *A24)
{
    fp2_t ixP;
    ec_point_t PnQ, nQ;
    fp2_copy(&ixP, &P->x);
    fp2_inv(&ixP);
    biext_ladder_2e(e, &PnQ, &nQ, PQ, Q, &ixP, A24);
    ratio(r, &PnQ, &nQ, P);
}

void
monodromy(fp2_t *r,
           const digit_t *n,
           const int n_bitlen,
           ec_point_t const *PQ,
           ec_point_t const *Q,
           ec_point_t const *P,
           fp2_t const *ixP,
           fp2_t const *ixQ,
           fp2_t const *ixPQ,
           ec_point_t const *A24)
{
    ec_point_t PnQ, nQ;
    biext_ladder(n, n_bitlen, &PnQ, &nQ, Q, P, PQ, ixQ, ixP, ixPQ, A24);
    ratio(r, &PnQ, &nQ, P);
}

// TODO: use only one inversion
// And normalize A24 at the same time (if needed), to save another inversion
void
to_cubical(ec_point_t *Q, ec_point_t *P)
{
    // ec_normalize_point(A24);
    ec_normalize_point(P);
    ec_normalize_point(Q);
    // ec_normalize_point(PQ);
}

// Normalize the points and also store 1/x(P), 1/x(Q)
void
to_cubical_i(ec_point_t *P, ec_point_t *Q, fp2_t *ixP, fp2_t *ixQ)
{
    /*
    //ec_normalize_point(A24);
    ec_normalize_point(P);
    ec_normalize_point(Q);
    //ec_normalize_point(PQ);
    fp2_copy(ixP, &P->x);
    fp2_inv(ixP);
    fp2_copy(ixQ, &Q->x);
    fp2_inv(ixQ);
    */
    fp2_t t[4];
    fp2_copy(&t[0], &P->x);
    fp2_copy(&t[1], &P->z);
    fp2_copy(&t[2], &Q->x);
    fp2_copy(&t[3], &Q->z);
    fp2_batched_inv(t, 4);
    fp2_mul(ixP, &P->z, &t[0]);
    fp2_mul(ixQ, &Q->z, &t[2]);
    fp2_mul(&P->x, &P->x, &t[1]);
    fp2_mul(&Q->x, &Q->x, &t[3]);
    fp2_set_one(&P->z);
    fp2_set_one(&Q->z);
}

void
to_cubical_odd_i(ec_point_t *P, ec_point_t *Q, ec_point_t *PQ, fp2_t *ixP, fp2_t *ixQ, fp2_t *ixPQ)
{
    fp2_t t[6];
    fp2_copy(&t[0], &P->x);
    fp2_copy(&t[1], &P->z);
    fp2_copy(&t[2], &Q->x);
    fp2_copy(&t[3], &Q->z);
    fp2_copy(&t[4], &PQ->x);
    fp2_copy(&t[5], &PQ->z);
    fp2_batched_inv(t, 6);
    fp2_mul(ixP, &P->z, &t[0]);
    fp2_mul(ixQ, &Q->z, &t[2]);
    fp2_mul(ixPQ, &PQ->z, &t[4]);
    fp2_mul(&P->x, &P->x, &t[1]);
    fp2_mul(&Q->x, &Q->x, &t[3]);
    fp2_mul(&PQ->x, &PQ->x, &t[5]);
    fp2_set_one(&P->z);
    fp2_set_one(&Q->z);
    fp2_set_one(&PQ->z);
}

/* (Do we need this?)
void to_cubical_c(ec_point_t* P, ec_point_t* A24, ec_point_t const* P_, ec_point_t const* A24_) {
    copy_point(P, P_);
    copy_point(A24, A24_);
    inline_to_cubical(P, A24);
}
*/

// non reduced Tate pairing, PQ should be P+Q in (X:Z) coordinates
// Assume the cubical points are normalised, and that we have 1/x(P)
// The _n suffix means we assume the points are normalised
void
non_reduced_tate_n(fp2_t *r,
                   uint64_t e,
                   ec_point_t *P,
                   ec_point_t *Q,
                   ec_point_t *PQ,
                   fp2_t const *ixP,
                   ec_point_t *A24)
{
    ec_point_t R;
    monodromy_i(&R, e, PQ, Q, P, ixP, A24);
    x_coord(r, &R);
}

// Same as above, but first normalise the points
void
non_reduced_tate(fp2_t *r,
                 uint64_t e,
                 ec_point_t *P,
                 ec_point_t *Q,
                 ec_point_t *PQ,
                 ec_point_t *A24)
{
    // to_cubical(Q, P);
    fp2_t ixP, ixQ;
    to_cubical_i(P, Q, &ixP, &ixQ); // TODO: ixQ not used
    non_reduced_tate_n(r, e, PQ, Q, P, &ixP, A24);
}

// P, Q, PQ are being normalized in the function.
// return the Weil pairing e_n(P,Q)^2
void weil_odd(fp2_t *r, const digit_t* n, const int n_bitlen, ec_point_t *P, ec_point_t *Q, ec_point_t *PQ, const ec_point_t *A24)
{
    fp2_t ixP, ixQ, ixPQ;
    fp2_t r1, r2;
    to_cubical_odd_i(P, Q, PQ, &ixP, &ixQ, &ixPQ);
    monodromy(&r1, n, n_bitlen, PQ, P, Q, &ixQ, &ixP, &ixPQ, A24);
    to_cubical_odd_i(P, Q, PQ, &ixP, &ixQ, &ixPQ);
    monodromy(&r2, n, n_bitlen, PQ, Q, P, &ixP, &ixQ, &ixPQ, A24);
    fp2_inv(&r2);
    fp2_mul(r, &r1, &r2);
}

// Used in PIKEB and PIKET  //D divides p-1
void tate_odd_m(fp_t *rr, const digit_t* n, const int n_bitlen, ec_point_t *P, ec_point_t *Q, ec_point_t *PQ, ec_point_t *A24)
{
    fp2_t ixP, ixQ, ixPQ;
    fp2_t r;
    to_cubical_odd_i(P, Q, PQ, &ixP, &ixQ, &ixPQ);
    monodromy(&r, n, n_bitlen, PQ, Q, P, &ixP, &ixQ, &ixPQ, A24);
    fp_mul(&r.re, &r.re, &MAGIC_FACTOR);
    fp_mul(&r.im, &r.im, &MAGIC_FACTOR);
    fp2_t tmp;
    fp2_copy(&tmp, &r);
    fp_neg(&tmp.im, &tmp.im);       //^(p)
    fp2_mul(&r, &r, &tmp);            //^(p+1)
    fp_copy(rr, &r.re);
    digit_t exp[NWORDS_ORDER] = {0};
    ibz_to_digits(exp, &TORSION_D_COFACTOR);
    fp_pow_vartime(rr, rr, exp, TORSION_D_COFACTOR->_mp_size);   //^C2 = ^((p-1)/D)
}

// Weil pairing, PQ should be P+Q in (X:Z) coordinates
// We assume the points are normalised correctly
// Do we need a weil_c version?
void
weil_n(fp2_t *r,
       uint64_t e,
       ec_point_t const *P,
       ec_point_t const *Q,
       ec_point_t const *PQ,
       fp2_t const *ixP,
       fp2_t const *ixQ,
       ec_point_t const *A24)
{
    ec_point_t R0, R1;
    monodromy_i(&R0, e, PQ, Q, P, ixP, A24);
    monodromy_i(&R1, e, PQ, P, Q, ixQ, A24);
    // TODO: check if that's the Weil pairing or its inverse
    fp2_mul(r, &R0.x, &R1.z);
    fp2_inv(r);
    fp2_mul(r, r, &R0.z);
    fp2_mul(r, r, &R1.x);
}

// Weil pairing, PQ should be P+Q in (X:Z) coordinates
// Normalise the points and call the code above
// The code will crash (division by 0) if either P or Q is (0:1)
void
weil(fp2_t *r, uint64_t e, ec_point_t *P, ec_point_t *Q, ec_point_t *PQ, ec_point_t *A24)
{
    fp2_t ixP, ixQ;
    to_cubical_i(P, Q, &ixP, &ixQ);
    weil_n(r, e, P, Q, PQ, &ixP, &ixQ, A24);
}

// a <- b^3
void fp2_pow3(fp2_t *a, const fp2_t *b) {
    fp2_t tmp;
    fp2_sqr(&tmp, b);
    fp2_mul(a, &tmp, b);
}

// Normalize a "basis" (P, Q), (P1, P2) and store their inverse
void
to_cubical_basis_i(ec_point_t *P,
                   ec_point_t *Q,
                   ec_point_t *P1,
                   ec_point_t *P2,
                   fp2_t *ixP,
                   fp2_t *ixQ,
                   fp2_t *ixP1,
                   fp2_t *ixP2)
{
    fp2_t t[8];
    fp2_copy(&t[0], &P->x);
    fp2_copy(&t[1], &P->z);
    fp2_copy(&t[2], &Q->x);
    fp2_copy(&t[3], &Q->z);
    fp2_copy(&t[4], &P1->x);
    fp2_copy(&t[5], &P1->z);
    fp2_copy(&t[6], &P2->x);
    fp2_copy(&t[7], &P2->z);
    fp2_batched_inv(t, 8);
    fp2_mul(ixP, &P->z, &t[0]);
    fp2_mul(&P->x, &P->x, &t[1]);
    fp2_set_one(&P->z);
    fp2_mul(ixQ, &Q->z, &t[2]);
    fp2_mul(&Q->x, &Q->x, &t[3]);
    fp2_set_one(&Q->z);
    fp2_mul(ixP1, &P1->z, &t[4]);
    fp2_mul(&P1->x, &P1->x, &t[5]);
    fp2_set_one(&P1->z);
    fp2_mul(ixP2, &P2->z, &t[6]);
    fp2_mul(&P2->x, &P2->x, &t[7]);
    fp2_set_one(&P2->z);
}
