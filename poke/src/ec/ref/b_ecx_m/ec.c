#include "curve_extras.h"
#include <ec_params.h>
#include <assert.h>
#include <stdio.h>
#include "const_PIKE.h"

void
ec_normalize_point(ec_point_t *P)
{
    fp2_inv(&P->z);
    fp2_mul(&P->x, &P->x, &P->z);
    fp2_set_one(&(P->z));
}

void
ec_normalize_curve(ec_curve_t *E)
{
    fp2_inv(&E->C);
    fp2_mul(&E->A, &E->A, &E->C);
    fp2_set_one(&E->C);
}

void
ec_curve_normalize_A24(ec_curve_t *E)
{
    if (E->is_A24_computed_and_normalized == 0) {
        AC_to_A24(&E->A24, E);
        ec_normalize_point(&E->A24);
        E->is_A24_computed_and_normalized = 1;
    }
}

bool
ec_is_zero(ec_point_t const *P)
{
    return fp2_is_zero(&P->z);
}

uint32_t
ec_has_zero_coordinate(const ec_point_t *P)
{
    return fp2_is_zero(&P->x) | fp2_is_zero(&P->z);
}

void
ec_point_init(ec_point_t *P)
{ // Initialize point as identity element (1:0)
    fp2_set_one(&(P->x));
    fp2_set_zero(&(P->z));
}

// Initialise the curve struct
void
ec_curve_init(ec_curve_t *E)
{
    // Initialise the constants
    fp2_set_zero(&(E->A));
    fp2_set_one(&(E->C));

    // Initalise the point (A+2 : 4C)
    ec_point_init(&(E->A24));

    // Set the bool to be false by default
    E->is_A24_computed_and_normalized = 0;
}

// TODO: this is simply ec_point_init, we don't need both?
void
ec_set_zero(ec_point_t *P)
{
    fp2_set_one(&(P->x));
    fp2_set_zero(&(P->z));
    return;
}

void
copy_curve(ec_curve_t *E1, ec_curve_t const *E2)
{
    fp2_copy(&(E1->A), &(E2->A));
    fp2_copy(&(E1->C), &(E2->C));
    E1->is_A24_computed_and_normalized = E2->is_A24_computed_and_normalized;
    copy_point(&E1->A24, &E2->A24);
}

void
xDBL_E0(ec_point_t *Q, const ec_point_t *P)
{ // Doubling of a Montgomery point in projective coordinates (X:Z) on the curve E0 with (A:C) = (0:1).
  // Input: projective Montgomery x-coordinates P = (XP:ZP), where xP=XP/ZP, and Montgomery curve constants (A:C) = (0:1). 
  // Output: projective Montgomery x-coordinates Q <- 2*P = (XQ:ZQ) such that x(2P)=XQ/ZQ.
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_add(&t1, &t1, &t1);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_add(&Q->z, &t1, &t2);
    fp2_mul(&Q->z, &Q->z, &t2);
}

void
xDBL(ec_point_t *Q, ec_point_t const *P, ec_point_t const *AC)
{
    // This version computes the coefficient values A+2C and 4C on-the-fly
    // The curve coefficients are passed via AC = (A:C)
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_add(&t3, &AC->z, &AC->z);
    fp2_mul(&t1, &t1, &t3);
    fp2_add(&t1, &t1, &t1);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_add(&t0, &t3, &AC->x);
    fp2_mul(&t0, &t0, &t2);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void
xDBL_A24(ec_point_t *Q, ec_point_t const *P, ec_point_t const *A24)
{
    // This version receives the coefficient value A24 = (A+2C:4C)
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_mul(&t1, &t1, &A24->z);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_mul(&t0, &t2, &A24->x);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void
xDBL_A24_normalized(ec_point_t *Q, ec_point_t const *P, ec_point_t const *A24)
{
    // This version receives the coefficient value A24 = (A+2C/4C:1)
    fp2_t t0, t1, t2;

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

void
xADD(ec_point_t *R, ec_point_t const *P, ec_point_t const *Q, ec_point_t const *PQ)
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
    fp2_sqr(&t2, &t2);
    fp2_sqr(&t3, &t3);
    fp2_mul(&t2, &PQ->z, &t2);
    fp2_mul(&R->z, &PQ->x, &t3);
    fp2_copy(&R->x, &t2);
}

void
xDBLADD(ec_point_t *R,
        ec_point_t *S,
        ec_point_t const *P,
        ec_point_t const *Q,
        ec_point_t const *PQ,
        ec_point_t const *A24)
{
    // Requires precomputation of A24 = (A+2C:4C)
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&R->x, &t0);
    fp2_sub(&t2, &Q->x, &Q->z);
    fp2_add(&S->x, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t2);
    fp2_sqr(&R->z, &t1);
    fp2_mul(&t1, &t1, &S->x);
    fp2_sub(&t2, &R->x, &R->z);
    fp2_mul(&R->z, &R->z, &A24->z);
    fp2_mul(&R->x, &R->x, &R->z);
    fp2_mul(&S->x, &A24->x, &t2);
    fp2_sub(&S->z, &t0, &t1);
    fp2_add(&R->z, &R->z, &S->x);
    fp2_add(&S->x, &t0, &t1);
    fp2_mul(&R->z, &R->z, &t2);
    fp2_sqr(&S->z, &S->z);
    fp2_sqr(&S->x, &S->x);
    fp2_mul(&S->z, &S->z, &PQ->x);
    fp2_mul(&S->x, &S->x, &PQ->z);
}

// assume that A24 has been normalized
void
xDBLADD_normalized(ec_point_t *R,
                   ec_point_t *S,
                   ec_point_t const *P,
                   ec_point_t const *Q,
                   ec_point_t const *PQ,
                   ec_point_t const *A24)
{
    // Requires precomputation of A24 = (A+2C:4C)
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&R->x, &t0);
    fp2_sub(&t2, &Q->x, &Q->z);
    fp2_add(&S->x, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t2);
    fp2_sqr(&R->z, &t1);
    fp2_mul(&t1, &t1, &S->x);
    fp2_sub(&t2, &R->x, &R->z);
    // fp2_mul(&R->z, &R->z, &A24->z);
    fp2_mul(&R->x, &R->x, &R->z);
    fp2_mul(&S->x, &A24->x, &t2);
    fp2_sub(&S->z, &t0, &t1);
    fp2_add(&R->z, &R->z, &S->x);
    fp2_add(&S->x, &t0, &t1);
    fp2_mul(&R->z, &R->z, &t2);
    fp2_sqr(&S->z, &S->z);
    fp2_sqr(&S->x, &S->x);
    fp2_mul(&S->z, &S->z, &PQ->x);
    fp2_mul(&S->x, &S->x, &PQ->z);
}

bool
is_point_equal(const ec_point_t *P, const ec_point_t *Q)
{ // Evaluate if two points in Montgomery coordinates (X:Z) are equal
  // Returns 1 (true) if P=Q, 0 (false) otherwise
    fp2_t t0, t1;

    if ((fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) || (fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z))) {
        return fp2_is_zero(&P->x) && fp2_is_zero(&P->z) && fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z);
    }

    fp2_mul(&t0, &P->x, &Q->z);
    fp2_mul(&t1, &Q->x, &P->z);
    fp2_sub(&t0, &t0, &t1);

    return fp2_is_zero(&t0);
}

uint32_t
ec_is_two_torsion(const ec_point_t *P, const ec_curve_t *E)
{
    if (ec_is_zero(P))
        return 0;

    uint32_t x_is_zero, tmp_is_zero;
    fp2_t t0, t1, t2;
    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_add(&t1, &t0, &t1);
    fp2_mul(&t2, &t2, &E->A);
    fp2_mul(&t1, &t1, &E->C);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t0, &t1, &t2); // 4 (CX^2+CZ^2+AXZ)

    x_is_zero = fp2_is_zero(&P->x);
    tmp_is_zero = fp2_is_zero(&t0);

    // two torsion if x or x^2 + Ax + 1 is zero
    return x_is_zero | tmp_is_zero;
}

void
select_point(ec_point_t *Q, ec_point_t const *P1, ec_point_t const *P2, const digit_t option)
{ // Select points
  // If option = 0 then Q <- P1, else if option = 0xFF...FF then Q <- P2
    fp2_select(&(Q->x), &(P1->x), &(P2->x), option);
    fp2_select(&(Q->z), &(P1->z), &(P2->z), option);
}

void
cswap_points(ec_point_t *P, ec_point_t *Q, const digit_t option)
{ // Swap points in constant time
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;

    for (int i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->x.re[i] ^ Q->x.re[i]);
        P->x.re[i] = temp ^ P->x.re[i];
        Q->x.re[i] = temp ^ Q->x.re[i];
        temp = option & (P->x.im[i] ^ Q->x.im[i]);
        P->x.im[i] = temp ^ P->x.im[i];
        Q->x.im[i] = temp ^ Q->x.im[i];
        temp = option & (P->z.re[i] ^ Q->z.re[i]);
        P->z.re[i] = temp ^ P->z.re[i];
        Q->z.re[i] = temp ^ Q->z.re[i];
        temp = option & (P->z.im[i] ^ Q->z.im[i]);
        P->z.im[i] = temp ^ P->z.im[i];
        Q->z.im[i] = temp ^ Q->z.im[i];
    }
}

void
swap_points(ec_point_t *P, ec_point_t *Q, const digit_t option)
{ // Swap points
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    fp2_cswap(&(P->x), &(Q->x), option);
    fp2_cswap(&(P->z), &(Q->z), option);
}

void
copy_point(ec_point_t *P, ec_point_t const *Q)
{
    fp2_copy(&(P->x), &(Q->x));
    fp2_copy(&(P->z), &(Q->z));
}

void
ec_neg(ec_point_t *res, const ec_point_t *P)
{
    // DOES NOTHING
    copy_point(res, P);
}

void
xMUL(ec_point_t *Q, ec_point_t const *P, digit_t const *k, ec_curve_t const *curve)
{
    ec_point_t R0, R1, A24;
    digit_t mask;
    unsigned int bit = 0, prevbit = 0, swap;

    fp2_add(&A24.x, &curve->C, &curve->C); // Precomputation of A24=(A+2C:4C)
    fp2_add(&A24.z, &A24.x, &A24.x);
    fp2_add(&A24.x, &A24.x, &curve->A);

    // R0 <- (1:0), R1 <- P
    ec_point_init(&R0);
    fp2_copy(&R1.x, &P->x);
    fp2_copy(&R1.z, &P->z);

    // Main loop
    for (int i = BITS - 1; i >= 0; i--) {
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX - 1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, &R0, &R1, P, &A24);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(&R0, &R1, mask);

    fp2_copy(&Q->x, &R0.x);
    fp2_copy(&Q->z, &R0.z);
}

void
xMULv2(ec_point_t *Q, ec_point_t const *P, digit_t const *k, const int kbits, ec_point_t const *A24)
{
    // This version receives the coefficient value A24 = (A+2C:4C)
    ec_point_t R0, R1;
    digit_t mask;
    unsigned int bit = 0, prevbit = 0, swap;

    // R0 <- (1:0), R1 <- P
    ec_point_init(&R0);
    fp2_copy(&R1.x, &P->x);
    fp2_copy(&R1.z, &P->z);

    // Main loop
    for (int i = kbits - 1; i >= 0; i--) {
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX - 1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, &R0, &R1, P, A24);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(&R0, &R1, mask);

    fp2_copy(&Q->x, &R0.x);
    fp2_copy(&Q->z, &R0.z);
}

// Compute S = k*P + l*Q, with PQ = P+Q
void
xDBLMUL(ec_point_t *S,
        ec_point_t const *P,
        digit_t const *k,
        ec_point_t const *Q,
        digit_t const *l,
        ec_point_t const *PQ,
        const ec_curve_t *curve)
{

    int i;
    digit_t evens, mevens, bitk0, bitl0, maskk, maskl, temp, bs1_ip1, bs2_ip1, bs1_i, bs2_i, h;
    digit_t sigma[2] = { 0 }, pre_sigma = 0;
    digit_t k_t[NWORDS_ORDER], l_t[NWORDS_ORDER], one[NWORDS_ORDER] = { 0 }, r[2 * BITS] = { 0 };
    ec_point_t A24, DIFF1a, DIFF1b, DIFF2a, DIFF2b, R[3] = { 0 }, T[3];

    // Derive sigma according to parity
    bitk0 = (k[0] & 1);
    bitl0 = (l[0] & 1);
    maskk = 0 - bitk0; // Parity masks: 0 if even, otherwise 1...1
    maskl = 0 - bitl0;
    sigma[0] = (bitk0 ^ 1);
    sigma[1] = (bitl0 ^ 1);
    evens = sigma[0] + sigma[1]; // Count number of even scalars
    mevens =
        0 - (evens & 1); // Mask mevens <- 0 if # even scalars = 0 or 2, otherwise mevens = 1...1

    // If k and l are both even or both odd, pick sigma = (0,1)
    sigma[0] = (sigma[0] & mevens);
    sigma[1] = (sigma[1] & mevens) | (1 & ~mevens);

    // Convert even scalars to odd
    one[0] = 1;
    mp_sub(k_t, k, one, NWORDS_ORDER);
    mp_sub(l_t, l, one, NWORDS_ORDER);
    select_ct(k_t, k_t, k, maskk, NWORDS_ORDER);
    select_ct(l_t, l_t, l, maskl, NWORDS_ORDER);

    // Scalar recoding
    for (i = 0; i < BITS; i++) {
        // If sigma[0] = 1 swap k_t and l_t
        maskk = 0 - (sigma[0] ^ pre_sigma);
        swap_ct(k_t, l_t, maskk, NWORDS_ORDER);

        if (i == BITS - 1) {
            bs1_ip1 = 0;
            bs2_ip1 = 0;
        } else {
            bs1_ip1 = mp_shiftr(k_t, 1, NWORDS_ORDER);
            bs2_ip1 = mp_shiftr(l_t, 1, NWORDS_ORDER);
        }
        bs1_i = k_t[0] & 1;
        bs2_i = l_t[0] & 1;

        r[2 * i] = bs1_i ^ bs1_ip1;
        r[2 * i + 1] = bs2_i ^ bs2_ip1;

        // Revert sigma if second bit, r_(2i+1), is 1
        pre_sigma = sigma[0];
        maskk = 0 - r[2 * i + 1];
        select_ct(&temp, &sigma[0], &sigma[1], maskk, 1);
        select_ct(&sigma[1], &sigma[1], &sigma[0], maskk, 1);
        sigma[0] = temp;
    }

    // Point initialization
    ec_point_init(&R[0]);
    maskk = 0 - sigma[0];
    select_point(&R[1], P, Q, maskk);
    select_point(&R[2], Q, P, maskk);

    fp2_copy(&DIFF1a.x, &R[1].x);
    fp2_copy(&DIFF1a.z, &R[1].z);
    fp2_copy(&DIFF1b.x, &R[2].x);
    fp2_copy(&DIFF1b.z, &R[2].z);

    // Initialize DIFF2a <- P+Q, DIFF2b <- P-Q
    xADD(&R[2], &R[1], &R[2], PQ);
    fp2_copy(&DIFF2a.x, &R[2].x);
    fp2_copy(&DIFF2a.z, &R[2].z);
    fp2_copy(&DIFF2b.x, &PQ->x);
    fp2_copy(&DIFF2b.z, &PQ->z);

    // fp2_add(&A24.x, &curve->C, &curve->C);    // Precomputation of A24=(A+2C:4C)
    // fp2_add(&A24.z, &A24.x, &A24.x);
    // fp2_add(&A24.x, &A24.x, &curve->A);

    // normalizing
    ec_curve_t E;
    copy_curve(&E, curve);
    ec_curve_normalize_A24(&E);
    copy_point(&A24, &E.A24);

    // Main loop
    for (i = BITS - 1; i >= 0; i--) {
        h = r[2 * i] + r[2 * i + 1]; // in {0, 1, 2}
        maskk = 0 - (h & 1);
        select_point(&T[0], &R[0], &R[1], maskk);
        maskk = 0 - (h >> 1);
        select_point(&T[0], &T[0], &R[2], maskk);
        xDBL_A24_normalized(&T[0], &T[0], &A24);

        maskk = 0 - r[2 * i + 1]; // in {0, 1}
        select_point(&T[1], &R[0], &R[1], maskk);
        select_point(&T[2], &R[1], &R[2], maskk);

        swap_points(&DIFF1a, &DIFF1b, maskk);
        xADD(&T[1], &T[1], &T[2], &DIFF1a);
        xADD(&T[2], &R[0], &R[2], &DIFF2a);

        // If hw (mod 2) = 1 then swap DIFF2a and DIFF2b
        maskk = 0 - (h & 1);
        swap_points(&DIFF2a, &DIFF2b, maskk);

        // R <- T
        copy_point(&R[0], &T[0]);
        copy_point(&R[1], &T[1]);
        copy_point(&R[2], &T[2]);
    }

    // Output R[evens]
    select_point(S, &R[0], &R[1], mevens);

    maskk = 0 - (bitk0 & bitl0);
    select_point(S, S, &R[2], maskk);
}

int xDBLMUL_bounded(ec_point_t *S,
        const ec_point_t *P,
        const digit_t *k,
        const ec_point_t *Q,
        const digit_t *l,
        const ec_point_t *PQ,
        const ec_curve_t *curve,
        const int kbits)
{ // The Montgomery biladder
  // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, scalars k and l of
  //         bitlength kbits, the difference PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants (A:C).
  // Output: projective Montgomery point S <- k*P + l*Q = (XS:ZS) such that x(k*P + l*Q)=XS/ZS.

    int i, A_is_zero;
    digit_t evens, mevens, bitk0, bitl0, maskk, maskl, temp, bs1_ip1, bs2_ip1, bs1_i, bs2_i, h;
    digit_t sigma[2] = { 0 }, pre_sigma = 0;
    digit_t k_t[NWORDS_ORDER], l_t[NWORDS_ORDER], one[NWORDS_ORDER] = { 0 }, r[2 * BITS] = { 0 };
    ec_point_t DIFF1a, DIFF1b, DIFF2a, DIFF2b, R[3] = { 0 }, T[3];

    // differential additions formulas are invalid in this case
    if (ec_has_zero_coordinate(P) | ec_has_zero_coordinate(Q) | ec_has_zero_coordinate(PQ))
        return 0;

    // Derive sigma according to parity
    bitk0 = (k[0] & 1);
    bitl0 = (l[0] & 1);
    maskk = 0 - bitk0; // Parity masks: 0 if even, otherwise 1...1
    maskl = 0 - bitl0;
    sigma[0] = (bitk0 ^ 1);
    sigma[1] = (bitl0 ^ 1);
    evens = sigma[0] + sigma[1]; // Count number of even scalars
    mevens = 0 - (evens & 1);    // Mask mevens <- 0 if # even of scalars = 0 or 2, otherwise mevens = 1...1

    // If k and l are both even or both odd, pick sigma = (0,1)
    sigma[0] = (sigma[0] & mevens);
    sigma[1] = (sigma[1] & mevens) | (1 & ~mevens);

    // Convert even scalars to odd
    one[0] = 1;
    mp_sub(k_t, k, one, NWORDS_ORDER);
    mp_sub(l_t, l, one, NWORDS_ORDER);
    select_ct(k_t, k_t, k, maskk, NWORDS_ORDER);
    select_ct(l_t, l_t, l, maskl, NWORDS_ORDER);

    // Scalar recoding
    for (i = 0; i < kbits; i++) {
        // If sigma[0] = 1 swap k_t and l_t
        maskk = 0 - (sigma[0] ^ pre_sigma);
        swap_ct(k_t, l_t, maskk, NWORDS_ORDER);

        if (i == kbits - 1) {
            bs1_ip1 = 0;
            bs2_ip1 = 0;
        } else {
            bs1_ip1 = mp_shiftr(k_t, 1, NWORDS_ORDER);
            bs2_ip1 = mp_shiftr(l_t, 1, NWORDS_ORDER);
        }
        bs1_i = k_t[0] & 1;
        bs2_i = l_t[0] & 1;

        r[2 * i] = bs1_i ^ bs1_ip1;
        r[2 * i + 1] = bs2_i ^ bs2_ip1;

        // Revert sigma if second bit, r_(2i+1), is 1
        pre_sigma = sigma[0];
        maskk = 0 - r[2 * i + 1];
        select_ct(&temp, &sigma[0], &sigma[1], maskk, 1);
        select_ct(&sigma[1], &sigma[1], &sigma[0], maskk, 1);
        sigma[0] = temp;
    }

    // Point initialization
    ec_point_init(&R[0]);
    maskk = 0 - sigma[0];
    select_point(&R[1], P, Q, maskk);
    select_point(&R[2], Q, P, maskk);

    fp2_copy(&DIFF1a.x, &R[1].x);
    fp2_copy(&DIFF1a.z, &R[1].z);
    fp2_copy(&DIFF1b.x, &R[2].x);
    fp2_copy(&DIFF1b.z, &R[2].z);

    // Initialize DIFF2a <- P+Q, DIFF2b <- P-Q
    xADD(&R[2], &R[1], &R[2], PQ);
    if (ec_has_zero_coordinate(&R[2]))
        return 0; // non valid formulas

    fp2_copy(&DIFF2a.x, &R[2].x);
    fp2_copy(&DIFF2a.z, &R[2].z);
    fp2_copy(&DIFF2b.x, &PQ->x);
    fp2_copy(&DIFF2b.z, &PQ->z);

    A_is_zero = fp2_is_zero(&curve->A);

    // Main loop
    for (i = kbits - 1; i >= 0; i--) {
        h = r[2 * i] + r[2 * i + 1]; // in {0, 1, 2}
        maskk = 0 - (h & 1);
        select_point(&T[0], &R[0], &R[1], maskk);
        maskk = 0 - (h >> 1);
        select_point(&T[0], &T[0], &R[2], maskk);
        if (A_is_zero) {
            xDBL_E0(&T[0], &T[0]);
        } else {
            assert(fp2_is_one(&curve->A24.z));
            xDBL_A24(&T[0], &T[0], &curve->A24);
        }

        maskk = 0 - r[2 * i + 1]; // in {0, 1}
        select_point(&T[1], &R[0], &R[1], maskk);
        select_point(&T[2], &R[1], &R[2], maskk);

        cswap_points(&DIFF1a, &DIFF1b, maskk);
        xADD(&T[1], &T[1], &T[2], &DIFF1a);
        xADD(&T[2], &R[0], &R[2], &DIFF2a);

        // If hw (mod 2) = 1 then swap DIFF2a and DIFF2b
        maskk = 0 - (h & 1);
        cswap_points(&DIFF2a, &DIFF2b, maskk);

        // R <- T
        copy_point(&R[0], &T[0]);
        copy_point(&R[1], &T[1]);
        copy_point(&R[2], &T[2]);
    }

    // Output R[evens]
    select_point(S, &R[0], &R[1], mevens);

    maskk = 0 - (bitk0 & bitl0);
    select_point(S, S, &R[2], maskk);
    return 1;
}

void
ec_ladder3pt(ec_point_t *R,
             digit_t *const m,
             ec_point_t const *P,
             ec_point_t const *Q,
             ec_point_t const *PQ,
             ec_curve_t const *A)
{

    // Curve constant in the form A24=(A+2C:4C)
    // ec_point_t A24;
    // fp2_add(&A24.z, &A->C, &A->C);
    // fp2_add(&A24.x, &A->A, &A24.z);
    // fp2_add(&A24.z, &A24.z, &A24.z);

    ec_point_t X0, X1, X2;
    copy_point(&X0, Q);
    copy_point(&X1, P);
    copy_point(&X2, PQ);

    int i, j;
    uint64_t t;
    for (i = 0; i < NWORDS_FIELD; i++) {
        t = 1;
        for (j = 0; j < 64; j++) {
            swap_points(&X1, &X2, -((t & m[i]) == 0));
            xDBLADD_normalized(&X0, &X1, &X0, &X1, &X2, &A->A24);
            swap_points(&X1, &X2, -((t & m[i]) == 0));
            t <<= 1;
        };
    };
    copy_point(R, &X1);
}

void
ec_j_inv(fp2_t *j_inv, const ec_curve_t *curve)
{
    /* j-invariant computation for montgommery coefficient A2=(A+2C:4C) */
    fp2_t t0, t1;

    fp2_sqr(&t1, &curve->C);
    fp2_sqr(j_inv, &curve->A);
    fp2_add(&t0, &t1, &t1);
    fp2_sub(&t0, j_inv, &t0);
    fp2_sub(&t0, &t0, &t1);
    fp2_sub(j_inv, &t0, &t1);
    fp2_sqr(&t1, &t1);
    fp2_mul(j_inv, j_inv, &t1);
    fp2_add(&t0, &t0, &t0);
    fp2_add(&t0, &t0, &t0);
    fp2_sqr(&t1, &t0);
    fp2_mul(&t0, &t0, &t1);
    fp2_add(&t0, &t0, &t0);
    fp2_add(&t0, &t0, &t0);
    fp2_inv(j_inv);
    fp2_mul(j_inv, &t0, j_inv);
}

static void
jac_init(jac_point_t *P)
{ // Initialize Montgomery in Jacobian coordinates as identity element (0:1:0)
    fp2_set_zero(&(P->x));
    fp2_set_one(&(P->y));
    fp2_set_zero(&(P->z));
}

bool
is_jac_equal(const jac_point_t *P, const jac_point_t *Q)
{ // Evaluate if two points in Jacobian coordinates (X:Y:Z) are equal
  // Returns 1 (true) if P=Q, 0 (false) otherwise
    fp2_t t0, t1, t2, t3;

    fp2_sqr(&t0, &Q->z);
    fp2_mul(&t2, &P->x, &t0); // x1*z2^2
    fp2_sqr(&t1, &P->z);
    fp2_mul(&t3, &Q->x, &t1); // x2*z1^2
    fp2_sub(&t2, &t2, &t3);

    fp2_mul(&t0, &t0, &Q->z);
    fp2_mul(&t0, &P->y, &t0); // y1*z2^3
    fp2_mul(&t1, &t1, &P->z);
    fp2_mul(&t1, &Q->y, &t1); // y2*z1^3
    fp2_sub(&t0, &t0, &t1);

    return fp2_is_zero(&t0) && fp2_is_zero(&t2);
}

void
jac_to_xz(ec_point_t *P, const jac_point_t *xyP)
{
    fp2_copy(&P->x, &xyP->x);
    fp2_copy(&P->z, &xyP->z);
    fp2_sqr(&P->z, &P->z);
}

void jac_to_proj(proj_point_t *P, const jac_point_t *xyP)
{
    fp2_mul(&P->x, &xyP->x, &xyP->z);
    fp2_copy(&P->y, &xyP->y);
    fp2_sqr(&P->z, &xyP->z);
    fp2_mul(&P->z, &P->z, &xyP->z);
}

bool
is_jac_xz_equal(const jac_point_t *P, const ec_point_t *Q)
{ // Evaluate if point P in Jacobian coordinates is equal to Q in homogeneous projective coordinates
  // (X:Z) Comparison is up to sign (only compares X and Z coordinates) Returns 1 (true) if P=Q, 0
  // (false) otherwise
    fp2_t t0, t1;

    fp2_mul(&t0, &P->x, &Q->z); // x1*z2
    fp2_sqr(&t1, &P->z);
    fp2_mul(&t1, &Q->x, &t1); // x2*z1^2
    fp2_sub(&t0, &t0, &t1);

    return fp2_is_zero(&t0);
}

void
copy_jac_point(jac_point_t *P, jac_point_t const *Q)
{
    fp2_copy(&(P->x), &(Q->x));
    fp2_copy(&(P->y), &(Q->y));
    fp2_copy(&(P->z), &(Q->z));
}

void
jac_neg(jac_point_t *Q, jac_point_t const *P)
{
    fp2_copy(&Q->x, &P->x);
    fp2_neg(&Q->y, &P->y);
    fp2_copy(&Q->z, &P->z);
}

void
DBL(jac_point_t *Q, jac_point_t const *P, ec_curve_t const *AC)
{ // Doubling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to
  // (X/Z^2,Y/Z^3) This version receives the coefficient value A
    fp2_t t0, t1, t2, t3;

    if (fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) {
        jac_init(Q);
        return;
    }

    fp2_sqr(&t0, &P->x); // t0 = x1^2
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t0, &t0, &t1); // t0 = 3x1^2
    fp2_sqr(&t1, &P->z);    // t1 = z1^2
    fp2_mul(&t2, &P->x, &AC->A);
    fp2_add(&t2, &t2, &t2); // t2 = 2Ax1
    fp2_add(&t2, &t1, &t2); // t2 = 2Ax1+z1^2
    fp2_mul(&t2, &t1, &t2); // t2 = z1^2(2Ax1+z1^2)
    fp2_add(&t2, &t0, &t2); // t2 = alpha = 3x1^2 + z1^2(2Ax1+z1^2)
    fp2_mul(&Q->z, &P->y, &P->z);
    fp2_add(&Q->z, &Q->z, &Q->z); // z2 = 2y1z1
    fp2_sqr(&t0, &Q->z);
    fp2_mul(&t0, &t0, &AC->A); // t0 = 4Ay1^2z1^2
    fp2_sqr(&t1, &P->y);
    fp2_add(&t1, &t1, &t1);     // t1 = 2y1^2
    fp2_add(&t3, &P->x, &P->x); // t3 = 2x1
    fp2_mul(&t3, &t1, &t3);     // t3 = 4x1y1^2
    fp2_sqr(&Q->x, &t2);        // x2 = alpha^2
    fp2_sub(&Q->x, &Q->x, &t0); // x2 = alpha^2 - 4Ay1^2z1^2
    fp2_sub(&Q->x, &Q->x, &t3);
    fp2_sub(&Q->x, &Q->x, &t3); // x2 = alpha^2 - 4Ay1^2z1^2 - 8x1y1^2
    fp2_sub(&Q->y, &t3, &Q->x); // y2 = 4x1y1^2 - x2
    fp2_mul(&Q->y, &Q->y, &t2); // y2 = alpha(4x1y1^2 - x2)
    fp2_sqr(&t1, &t1);          // t1 = 4y1^4
    fp2_sub(&Q->y, &Q->y, &t1);
    fp2_sub(&Q->y, &Q->y, &t1); // y2 = alpha(4x1y1^2 - x2) - 8y1^4
}

void
ADD(jac_point_t *R, jac_point_t const *P, jac_point_t const *Q, ec_curve_t const *AC)
{ // Addition on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to
  // (X/Z^2,Y/Z^3) This version receives the coefficient value A
    fp2_t t0, t1, t2, t3, t4, t5, t6;
    jac_point_t T;

    if (is_jac_equal(P, Q)) {
        DBL(R, P, AC);
        return;
    }
    jac_neg(&T, P);
    if (is_jac_equal(&T, Q)) {
        jac_init(R);
        return;
    }
    if (fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) {
        copy_jac_point(R, Q);
        return;
    } else if (fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z)) {
        copy_jac_point(R, P);
        return;
    }

    fp2_sqr(&t0, &P->z);        // t0 = z1^2
    fp2_mul(&t1, &t0, &P->z);   // t1 = z1^3
    fp2_sqr(&t2, &Q->z);        // t2 = z2^2
    fp2_mul(&t3, &t2, &Q->z);   // t3 = z2^3
    fp2_mul(&t1, &t1, &Q->y);   // t1 = y2z1^3
    fp2_mul(&t3, &t3, &P->y);   // t3 = y1z2^3
    fp2_sub(&t1, &t1, &t3);     // t1 = lambda1 = y2z1^3 - y1z2^3
    fp2_mul(&t0, &t0, &Q->x);   // t0 = x2z1^2
    fp2_mul(&t2, &t2, &P->x);   // t2 = x1z2^2
    fp2_sub(&t4, &t0, &t2);     // t4 = lambda3 = x2z1^2 - x1z2^2
    fp2_add(&t0, &t0, &t2);     // t0 = lambda2 = x2z1^2 + x1z2^2
    fp2_mul(&t5, &P->z, &Q->z); // t5 = z1z2
    fp2_mul(&R->z, &t4, &t5);   // z3 = z1z2*lambda3
    fp2_sqr(&t5, &t5);          // t5 = z1^2z2^2
    fp2_mul(&t5, &AC->A, &t5);  // t5 = Az1^2z2^2
    fp2_add(&t0, &t0, &t5);     // t0 = Az1^2z2^2 + lambda2
    fp2_sqr(&t6, &t4);          // t6 = lambda3^2
    fp2_mul(&t5, &t0, &t6);     // t5 = lambda3^2(Az1^2z2^2 + lambda2)
    fp2_sqr(&R->x, &t1);        // x3 = lambda1^2
    fp2_sub(&R->x, &R->x, &t5); // x3 = lambda1^2 - lambda3^2(Az1^2z2^2 + lambda2)
    fp2_mul(&t3, &t3, &t4);     // t3 = y1z2^3*lambda3
    fp2_mul(&t3, &t3, &t6);     // t3 = y1z2^3*lambda3^3
    fp2_mul(&t2, &t2, &t6);     // t2 = x1z2^2*lambda3^2
    fp2_sub(&R->y, &t2, &R->x); // y3 = x1z2^2*lambda3^2 - x3
    fp2_mul(&R->y, &R->y, &t1); // y3 = lambda1(x1z2^2*lambda3^2 - x3)
    fp2_sub(&R->y, &R->y, &t3); // y3 = lambda1(x1z2^2*lambda3^2 - x3) - y1z2^3*lambda3^3
}

void
TPL(jac_point_t *Q, jac_point_t const *P, ec_curve_t const *AC)
{ // Naive tripling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z)
  // corresponding to (X/Z^2,Y/Z^3) This version receives the coefficient value A
    jac_point_t R;

    DBL(&R, P, AC);
    ADD(Q, &R, P, AC);
}

void
MUL_FIVE(jac_point_t *Q, jac_point_t const *P, ec_curve_t const *AC)
{ // Naive tripling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z)
  // corresponding to (X/Z^2,Y/Z^3) This version receives the coefficient value A
    jac_point_t R;

    DBL(&R, P, AC);
    DBL(&R, &R, AC);
    ADD(Q, P, &R, AC);
}

void
recover_y(fp2_t *y, fp2_t const *Px, ec_curve_t const *curve)
{ // Recover y-coordinate of a point on the Montgomery curve y^2 = x^3 + Ax^2 + x
    fp2_t t0;

    fp2_sqr(&t0, Px);
    fp2_mul(y, &t0, &curve->A); // Ax^2
    fp2_add(y, y, Px);          // Ax^2 + x
    fp2_mul(&t0, &t0, Px);
    fp2_add(y, y, &t0); // x^3 + Ax^2 + x
    fp2_sqrt(y);
}

// scalar multiplication R <- k*P
void MUL_generic(jac_point_t *R,
                const jac_point_t *P,
                const digit_t *k,
                const int nwords,
                const ec_curve_t *curve
            )
{
    digit_t k_t;
    jac_point_t R_t;

    jac_init(&R_t);
    for(int i = nwords; i > 0; i--) {
        for (int j = 0; j < RADIX; j++) {
            k_t = k[i - 1] >> (RADIX - 1 - j);
            k_t &= 0x1;
            DBL(&R_t, &R_t, curve);
            if (k_t == 1) {
                ADD(&R_t, &R_t, P, curve);
            }
        }
    }
    copy_jac_point(R, &R_t);
}

// Double-scalar multiplication R <- k*P + l*Q, fixed for (64*size)-bit scalars
void
DBLMUL_generic(jac_point_t *R,
               const jac_point_t *P,
               const digit_t *k,
               const jac_point_t *Q,
               const digit_t *l,
               const ec_curve_t *curve,
               int size)
{

    digit_t k_t, l_t;
    jac_point_t PQ;

    ADD(&PQ, P, Q, curve);
    jac_init(R);
    for (int j = size; j > 0; j--) {
        for (int i = 0; i < 64; i++) {
            k_t = k[j - 1] >> (63 - i);
            k_t &= 0x01;
            l_t = l[j - 1] >> (63 - i);
            l_t &= 0x01;
            DBL(R, R, curve);
            if (k_t == 1 && l_t == 1) {
                ADD(R, R, &PQ, curve);
            } else if (k_t == 1) {
                ADD(R, R, P, curve);
            } else if (l_t == 1) {
                ADD(R, R, Q, curve);
            }
        }
    }
}

void
DBLMUL(jac_point_t *R,
       const jac_point_t *P,
       const digit_t k,
       const jac_point_t *Q,
       const digit_t l,
       const ec_curve_t *curve)
{ // Double-scalar multiplication R <- k*P + l*Q, fixed for 64-bit scalars
    digit_t k_t, l_t;
    jac_point_t PQ;

    ADD(&PQ, P, Q, curve);
    jac_init(R);

    for (int i = 0; i < 64; i++) {
        k_t = k >> (63 - i);
        k_t &= 0x01;
        l_t = l >> (63 - i);
        l_t &= 0x01;
        DBL(R, R, curve);
        if (k_t == 1 && l_t == 1) {
            ADD(R, R, &PQ, curve);
        } else if (k_t == 1) {
            ADD(R, R, P, curve);
        } else if (l_t == 1) {
            ADD(R, R, Q, curve);
        }
    }
}

#define SCALAR_BITS 128 //// TODO: this could be defined by level

void
DBLMUL2(jac_point_t *R,
        const jac_point_t *P,
        const digit_t *k,
        const jac_point_t *Q,
        const digit_t *l,
        const ec_curve_t *curve)
{ // Double-scalar multiplication R <- k*P + l*Q, fixed for 128-bit scalars
    digit_t k_t[2], l_t[2];
    jac_point_t PQ;

    ADD(&PQ, P, Q, curve);
    jac_init(R);

    for (int i = 0; i < SCALAR_BITS - 64; i++) {
        k_t[1] = k[1] >> (SCALAR_BITS - 64 - 1 - i);
        k_t[1] &= 0x01;
        l_t[1] = l[1] >> (SCALAR_BITS - 64 - 1 - i);
        l_t[1] &= 0x01;
        DBL(R, R, curve);
        if (k_t[1] == 1 && l_t[1] == 1) {
            ADD(R, R, &PQ, curve);
        } else if (k_t[1] == 1) {
            ADD(R, R, P, curve);
        } else if (l_t[1] == 1) {
            ADD(R, R, Q, curve);
        }
    }

    for (int i = 0; i < 64; i++) {
        k_t[0] = k[0] >> (64 - 1 - i);
        k_t[0] &= 0x01;
        l_t[0] = l[0] >> (64 - 1 - i);
        l_t[0] &= 0x01;
        DBL(R, R, curve);
        if (k_t[0] == 1 && l_t[0] == 1) {
            ADD(R, R, &PQ, curve);
        } else if (k_t[0] == 1) {
            ADD(R, R, P, curve);
        } else if (l_t[0] == 1) {
            ADD(R, R, Q, curve);
        }
    }
}

void
lift_point(jac_point_t *P, ec_point_t *Q, ec_curve_t *E)
{
    fp2_t inverses[2];
    if (fp2_is_zero(&Q->z)) {
        fp2_set_one(&P->x);
        fp2_set_one(&P->y);
        fp2_set_zero(&P->z);
    } else {
        fp2_copy(&inverses[0], &Q->z);
        fp2_copy(&inverses[1], &E->C);
        fp2_batched_inv(inverses, 2);

        fp2_set_one(&Q->z);
        fp2_set_one(&E->C);

        fp2_mul(&Q->x, &Q->x, &inverses[0]);
        fp2_mul(&E->A, &E->A, &inverses[1]);

        fp2_copy(&P->x, &Q->x);
        fp2_set_one(&P->z);
        recover_y(&P->y, &P->x, E);
    }
}

void
lift_basis(jac_point_t *P, jac_point_t *Q, ec_basis_t *B, ec_curve_t *E)
{
    fp2_t inverses[2];
    fp2_copy(&inverses[0], &B->P.z);
    // fp2_copy(&inverses[1],&B->Q.z);
    fp2_copy(&inverses[1], &E->C);

    fp2_batched_inv(inverses, 2);
    fp2_set_one(&B->P.z);
    fp2_set_one(&E->C);

    fp2_mul(&B->P.x, &B->P.x, &inverses[0]);
    fp2_mul(&E->A, &E->A, &inverses[1]);

    fp2_copy(&P->x, &B->P.x);
    fp2_copy(&Q->x, &B->Q.x);
    fp2_copy(&Q->z, &B->Q.z);
    fp2_set_one(&P->z);
    recover_y(&P->y, &P->x, E);

    // Algorithm of Okeya-Sakurai to recover y.Q in the montgomery model
    fp2_t v1, v2, v3, v4;
    fp2_mul(&v1, &P->x, &Q->z);
    fp2_add(&v2, &Q->x, &v1);
    fp2_sub(&v3, &Q->x, &v1);
    fp2_sqr(&v3, &v3);
    fp2_mul(&v3, &v3, &B->PmQ.x);
    fp2_add(&v1, &E->A, &E->A);
    fp2_mul(&v1, &v1, &Q->z);
    fp2_add(&v2, &v2, &v1);
    fp2_mul(&v4, &P->x, &Q->x);
    fp2_add(&v4, &v4, &Q->z);
    fp2_mul(&v2, &v2, &v4);
    fp2_mul(&v1, &v1, &Q->z);
    fp2_sub(&v2, &v2, &v1);
    fp2_mul(&v2, &v2, &B->PmQ.z);
    fp2_sub(&Q->y, &v3, &v2);
    fp2_add(&v1, &P->y, &P->y);
    fp2_mul(&v1, &v1, &Q->z);
    fp2_mul(&v1, &v1, &B->PmQ.z);
    fp2_mul(&Q->x, &Q->x, &v1);
    fp2_mul(&Q->z, &Q->z, &v1);

    // Transforming to a jacobian coordinate
    fp2_sqr(&v1, &Q->z);
    fp2_mul(&Q->y, &Q->y, &v1);
    fp2_mul(&Q->x, &Q->x, &Q->z);
}

// WRAPPERS

void
ec_dbl(ec_point_t *res, const ec_curve_t *curve, const ec_point_t *P)
{
    xDBL(res, P, (ec_point_t const *)curve);
}

void
ec_dbl_iter(ec_point_t *res, int n, ec_curve_t *curve, const ec_point_t *P)
{
    if (n > 0) {

        // then it's probably worth it to normalize
        if (n > 50) {
            ec_curve_normalize_A24(curve);
            xDBL_A24(res, P, &curve->A24);
            for (int i = 0; i < n - 1; i++) {
                xDBL_A24(res, res, &curve->A24);
            }
        } else {

            // TODO might still be worth it to compute the A24 and use xDBLv2
            ec_dbl(res, curve, P);
            for (int i = 0; i < n - 1; i++) {
                ec_dbl(res, curve, res);
            }
        }
    }
}

void
ec_mul(ec_point_t *res, const ec_curve_t *curve, const digit_t *scalar, const ec_point_t *P)
{
    xMUL(res, P, scalar, curve);
}

void
ec_biscalar_mul(ec_point_t *res,
                const ec_curve_t *curve,
                const digit_t *scalarP,
                const digit_t *scalarQ,
                const ec_basis_t *PQ)
{
    xDBLMUL(res, &PQ->P, scalarP, &PQ->Q, scalarQ, &PQ->PmQ, curve);
}

void
xMULv2_normalized(ec_point_t *Q, ec_point_t const *P, digit_t const *k, const int kbits, ec_point_t const *A24)
{
    // This version assumes the coefficient value A24 = (A+2C/4C:1) has been normalized
    ec_point_t R0, R1;
    digit_t mask;
    unsigned int bit = 0, prevbit = 0, swap;

    // R0 <- (1:0), R1 <- P
    ec_point_init(&R0);
    fp2_copy(&R1.x, &P->x);
    fp2_copy(&R1.z, &P->z);

    // Main loop
    for (int i = kbits - 1; i >= 0; i--) {
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX - 1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(&R0, &R1, mask);
        
        // Call the optimized normalized addition/doubling function
        xDBLADD_normalized(&R0, &R1, &R0, &R1, P, A24);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(&R0, &R1, mask);

    fp2_copy(&Q->x, &R0.x);
    fp2_copy(&Q->z, &R0.z);
}

void
ec_mul_normalized_A24(ec_point_t *res, const ec_curve_t *curve, const digit_t *scalar, const ec_point_t *P, const int kbits)
{
    ec_curve_t E;
    
    // Copy the curve to avoid mutating the original constant curve pointer
    copy_curve(&E, curve);
    
    // Normalize the curve to compute A24 and ensure its Z coordinate is 1
    ec_curve_normalize_A24(&E);
    
    // Execute the scalar multiplication using the normalized A24 parameter
    xMULv2_normalized(res, P, scalar, kbits, &E.A24);
}


int xDBLMUL_bounded_normalized(ec_point_t *S,
        const ec_point_t *P,
        const digit_t *k,
        const ec_point_t *Q,
        const digit_t *l,
        const ec_point_t *PQ,
        const ec_curve_t *curve,
        const int kbits)
{ // The Montgomery biladder
  // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, scalars k and l of
  //         bitlength kbits, the difference PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants (A:C).
  // Output: projective Montgomery point S <- k*P + l*Q = (XS:ZS) such that x(k*P + l*Q)=XS/ZS.

    int i, A_is_zero;
    digit_t evens, mevens, bitk0, bitl0, maskk, maskl, temp, bs1_ip1, bs2_ip1, bs1_i, bs2_i, h;
    digit_t sigma[2] = { 0 }, pre_sigma = 0;
    digit_t k_t[NWORDS_ORDER], l_t[NWORDS_ORDER], one[NWORDS_ORDER] = { 0 }, r[2 * BITS] = { 0 };
    ec_point_t DIFF1a, DIFF1b, DIFF2a, DIFF2b, R[3] = { 0 }, T[3];

    // differential additions formulas are invalid in this case
    if (ec_has_zero_coordinate(P) | ec_has_zero_coordinate(Q) | ec_has_zero_coordinate(PQ))
        return 0;

    // Derive sigma according to parity
    bitk0 = (k[0] & 1);
    bitl0 = (l[0] & 1);
    maskk = 0 - bitk0; // Parity masks: 0 if even, otherwise 1...1
    maskl = 0 - bitl0;
    sigma[0] = (bitk0 ^ 1);
    sigma[1] = (bitl0 ^ 1);
    evens = sigma[0] + sigma[1]; // Count number of even scalars
    mevens = 0 - (evens & 1);    // Mask mevens <- 0 if # even of scalars = 0 or 2, otherwise mevens = 1...1

    // If k and l are both even or both odd, pick sigma = (0,1)
    sigma[0] = (sigma[0] & mevens);
    sigma[1] = (sigma[1] & mevens) | (1 & ~mevens);

    // Convert even scalars to odd
    one[0] = 1;
    mp_sub(k_t, k, one, NWORDS_ORDER);
    mp_sub(l_t, l, one, NWORDS_ORDER);
    select_ct(k_t, k_t, k, maskk, NWORDS_ORDER);
    select_ct(l_t, l_t, l, maskl, NWORDS_ORDER);

    // Scalar recoding
    for (i = 0; i < kbits; i++) {
        // If sigma[0] = 1 swap k_t and l_t
        maskk = 0 - (sigma[0] ^ pre_sigma);
        swap_ct(k_t, l_t, maskk, NWORDS_ORDER);

        if (i == kbits - 1) {
            bs1_ip1 = 0;
            bs2_ip1 = 0;
        } else {
            bs1_ip1 = mp_shiftr(k_t, 1, NWORDS_ORDER);
            bs2_ip1 = mp_shiftr(l_t, 1, NWORDS_ORDER);
        }
        bs1_i = k_t[0] & 1;
        bs2_i = l_t[0] & 1;

        r[2 * i] = bs1_i ^ bs1_ip1;
        r[2 * i + 1] = bs2_i ^ bs2_ip1;

        // Revert sigma if second bit, r_(2i+1), is 1
        pre_sigma = sigma[0];
        maskk = 0 - r[2 * i + 1];
        select_ct(&temp, &sigma[0], &sigma[1], maskk, 1);
        select_ct(&sigma[1], &sigma[1], &sigma[0], maskk, 1);
        sigma[0] = temp;
    }

    // Point initialization
    ec_point_init(&R[0]);
    maskk = 0 - sigma[0];
    select_point(&R[1], P, Q, maskk);
    select_point(&R[2], Q, P, maskk);

    fp2_copy(&DIFF1a.x, &R[1].x);
    fp2_copy(&DIFF1a.z, &R[1].z);
    fp2_copy(&DIFF1b.x, &R[2].x);
    fp2_copy(&DIFF1b.z, &R[2].z);

    // Initialize DIFF2a <- P+Q, DIFF2b <- P-Q
    xADD(&R[2], &R[1], &R[2], PQ);
    if (ec_has_zero_coordinate(&R[2]))
        return 0; // non valid formulas

    fp2_copy(&DIFF2a.x, &R[2].x);
    fp2_copy(&DIFF2a.z, &R[2].z);
    fp2_copy(&DIFF2b.x, &PQ->x);
    fp2_copy(&DIFF2b.z, &PQ->z);

    A_is_zero = fp2_is_zero(&curve->A);

    // Main loop
    for (i = kbits - 1; i >= 0; i--) {
        h = r[2 * i] + r[2 * i + 1]; // in {0, 1, 2}
        maskk = 0 - (h & 1);
        select_point(&T[0], &R[0], &R[1], maskk);
        maskk = 0 - (h >> 1);
        select_point(&T[0], &T[0], &R[2], maskk);
        if (A_is_zero) {
            xDBL_E0(&T[0], &T[0]);
        } else {
            assert(fp2_is_one(&curve->A24.z));
            xDBL_A24_normalized(&T[0], &T[0], &curve->A24);
        }

        maskk = 0 - r[2 * i + 1]; // in {0, 1}
        select_point(&T[1], &R[0], &R[1], maskk);
        select_point(&T[2], &R[1], &R[2], maskk);

        cswap_points(&DIFF1a, &DIFF1b, maskk);
        xADD(&T[1], &T[1], &T[2], &DIFF1a);
        xADD(&T[2], &R[0], &R[2], &DIFF2a);

        // If hw (mod 2) = 1 then swap DIFF2a and DIFF2b
        maskk = 0 - (h & 1);
        cswap_points(&DIFF2a, &DIFF2b, maskk);

        // R <- T
        copy_point(&R[0], &T[0]);
        copy_point(&R[1], &T[1]);
        copy_point(&R[2], &T[2]);
    }

    // Output R[evens]
    select_point(S, &R[0], &R[1], mevens);

    maskk = 0 - (bitk0 & bitl0);
    select_point(S, S, &R[2], maskk);
    return 1;
}

int
ec_biscalar_mul_bounded(ec_point_t *res,
                        const ec_curve_t *curve,
                        const digit_t *scalarP,
                        const digit_t *scalarQ,
                        const ec_basis_t *PQ,
                        int kbits)
{
    if (fp2_is_zero(&PQ->PmQ.z))
        return 0;

    /* Differential additions behave badly when PmQ = (0:1), so we need to
     * treat this case specifically. Since we assume P, Q are a basis, this
     * can happen only if kbits==1 */
    if (kbits == 1) {
        // Sanity check: our basis should be given by 2-torsion points
        if (!ec_is_two_torsion(&PQ->P, curve) || !ec_is_two_torsion(&PQ->Q, curve) ||
            !ec_is_two_torsion(&PQ->PmQ, curve))
            return 0;
        digit_t bP, bQ;
        bP = (scalarP[0] & 1);
        bQ = (scalarQ[0] & 1);
        if (bP == 0 && bQ == 0)
            ec_point_init(res); //(1: 0)
        else if (bP == 1 && bQ == 0)
            copy_point(res, &PQ->P);
        else if (bP == 0 && bQ == 1)
            copy_point(res, &PQ->Q);
        else if (bP == 1 && bQ == 1)
            copy_point(res, &PQ->PmQ);
        else // should never happen
            assert(0);
        return 1;
    } else {
        ec_curve_t E;
        copy_curve(&E, curve);

        if (!fp2_is_zero(&curve->A)) { // If A is not zero normalize
            ec_curve_normalize_A24(&E);
        }
        return xDBLMUL_bounded_normalized(res, &PQ->P, scalarP, &PQ->Q, scalarQ, &PQ->PmQ, (const ec_curve_t *)&E, kbits);
    }
}

void
recover_y_twist(fp2_t *y, fp2_t const *Px, ec_curve_t const *curve)
{ // Recover y-coordinate of a point on the Montgomery curve y^2 = x^3 + Ax^2 + x
    fp2_t t0;

    fp2_sqr(&t0, Px);
    fp2_mul(y, &t0, &curve->A); // Ax^2
    fp2_add(y, y, Px);          // Ax^2 + x
    fp2_mul(&t0, &t0, Px);
    fp2_add(y, y, &t0); // x^3 + Ax^2 + x

    fp2_mul(y, y , &B_qnr_inv);
    fp2_sqrt(y);
}

void
lift_point_twist(jac_point_t *P, ec_point_t *Q, ec_curve_t *E)
{
    fp2_t inverses[2];
    if (fp2_is_zero(&Q->z)) {
        fp2_set_one(&P->x);
        fp2_set_one(&P->y);
        fp2_set_zero(&P->z);
    } else {
        fp2_copy(&inverses[0], &Q->z);
        fp2_copy(&inverses[1], &E->C);
        fp2_batched_inv(inverses, 2);

        fp2_set_one(&Q->z);
        fp2_set_one(&E->C);

        fp2_mul(&Q->x, &Q->x, &inverses[0]);
        fp2_mul(&E->A, &E->A, &inverses[1]);

        fp2_copy(&P->x, &Q->x);
        fp2_set_one(&P->z);
        recover_y_twist(&P->y, &P->x, E);
    }
}

void
DBL_twist(jac_point_t *Q, jac_point_t const *P, ec_curve_t const *AC)
{ // Doubling on a twisted Montgomery curve By^2 = x^3 + Ax^2 + x, representation in Jacobian coordinates (X:Y:Z) corresponding to
  // (X/Z^2,Y/Z^3) This version receives the coefficient value A and twisting factor B
    fp2_t t0, t1, t2, t3;

    if (fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) {
        jac_init(Q);
        return;
    }

    fp2_sqr(&t0, &P->x); // t0 = x1^2
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t0, &t0, &t1); // t0 = 3x1^2
    fp2_sqr(&t1, &P->z);    // t1 = z1^2
    fp2_mul(&t2, &P->x, &AC->A);
    fp2_add(&t2, &t2, &t2); // t2 = 2Ax1
    fp2_add(&t2, &t1, &t2); // t2 = 2Ax1+z1^2
    fp2_mul(&t2, &t1, &t2); // t2 = z1^2(2Ax1+z1^2)
    fp2_add(&t2, &t0, &t2); // t2 = alpha = 3x1^2 + z1^2(2Ax1+z1^2)

    fp2_mul(&Q->z, &P->y, &P->z);
    fp2_add(&Q->z, &Q->z, &Q->z); // z2 = 2y1z1
    fp2_mul(&Q->z, &Q->z, &B_qnr); // z2 = 2 * B * y1 * z1 (NEW)

    fp2_sqr(&t0, &Q->z);
    fp2_mul(&t0, &t0, &AC->A); // t0 = A * Z2^2 = 4 * A * B^2 * y1^2 * z1^2

    fp2_sqr(&t1, &P->y);
    fp2_add(&t1, &t1, &t1);     // t1 = 2y1^2
    fp2_mul(&t1, &t1, &B_qnr);   // t1 = 2 * B * y1^2 (NEW)

    fp2_add(&t3, &P->x, &P->x); // t3 = 2x1
    fp2_mul(&t3, &t1, &t3);     // t3 = 4 * B * x1 * y1^2
    fp2_mul(&t3, &t3, &B_qnr);   // t3 = 4 * B^2 * x1 * y1^2 (NEW)

    fp2_sqr(&Q->x, &t2);          // x2 = alpha^2
    fp2_mul(&Q->x, &Q->x, &B_qnr); // x2 = B * alpha^2 (NEW)
    fp2_sub(&Q->x, &Q->x, &t0);   // x2 = B*alpha^2 - A*Z2^2
    fp2_sub(&Q->x, &Q->x, &t3);
    fp2_sub(&Q->x, &Q->x, &t3);   // x2 = B*alpha^2 - A*Z2^2 - 8 * B^2 * x1 * y1^2

    fp2_sub(&Q->y, &t3, &Q->x); // y2 = 4*B^2*x1*y1^2 - x2
    fp2_mul(&Q->y, &Q->y, &t2); // y2 = alpha(4*B^2*x1*y1^2 - x2)
    
    fp2_sqr(&t1, &t1);          // t1 = 4 * B^2 * y1^4
    fp2_mul(&t1, &t1, &B_qnr);   // t1 = 4 * B^3 * y1^4 (NEW)
    
    fp2_sub(&Q->y, &Q->y, &t1);
    fp2_sub(&Q->y, &Q->y, &t1); // y2 = alpha(4*B^2*x1*y1^2 - x2) - 8 * B^3 * y1^4
}

void
DBL_twist_inv(jac_point_t *Q, jac_point_t const *P, ec_curve_t const *AC)
{ // Doubling on a twisted Montgomery curve By^2 = x^3 + Ax^2 + x
  // Uses explicitly calculated/precomputed B^{-1} to avoid cascading Z scaling
    fp2_t t0, t1, t2, t3;

    if (fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) {
        jac_init(Q);
        return;
    }

    fp2_sqr(&t0, &P->x); // t0 = x1^2
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t0, &t0, &t1); // t0 = 3x1^2
    fp2_sqr(&t1, &P->z);    // t1 = z1^2
    fp2_mul(&t2, &P->x, &AC->A);
    fp2_add(&t2, &t2, &t2); // t2 = 2Ax1
    fp2_add(&t2, &t1, &t2); // t2 = 2Ax1+z1^2
    fp2_mul(&t2, &t1, &t2); // t2 = z1^2(2Ax1+z1^2)
    fp2_add(&t2, &t0, &t2); // t2 = alpha = 3x1^2 + z1^2(2Ax1+z1^2)

    fp2_mul(&Q->z, &P->y, &P->z);
    fp2_add(&Q->z, &Q->z, &Q->z); // z2 = 2y1z1

    fp2_sqr(&t0, &Q->z);
    fp2_mul(&t0, &t0, &AC->A); // t0 = 4Ay1^2z1^2
    fp2_sqr(&t1, &P->y);
    fp2_add(&t1, &t1, &t1);     // t1 = 2y1^2
    fp2_add(&t3, &P->x, &P->x); // t3 = 2x1
    fp2_mul(&t3, &t1, &t3);     // t3 = 4x1y1^2

    fp2_sqr(&Q->x, &t2);               // x2 = alpha^2
    fp2_mul(&Q->x, &Q->x, &B_qnr_inv); // x2 = B^{-1} * alpha^2 (NEW)
    fp2_sub(&Q->x, &Q->x, &t0);        // x2 = B^{-1}*alpha^2 - 4Ay1^2z1^2
    fp2_sub(&Q->x, &Q->x, &t3);
    fp2_sub(&Q->x, &Q->x, &t3);        // x2 = B^{-1}*alpha^2 - 4Ay1^2z1^2 - 8x1y1^2

    fp2_sub(&Q->y, &t3, &Q->x);        // y2 = 4x1y1^2 - x2
    fp2_mul(&Q->y, &Q->y, &t2);        // y2 = alpha(4x1y1^2 - x2)
    fp2_mul(&Q->y, &Q->y, &B_qnr_inv); // y2 = B^{-1} * alpha(4x1y1^2 - x2) (NEW)
    fp2_sqr(&t1, &t1);                 // t1 = 4y1^4
    fp2_sub(&Q->y, &Q->y, &t1);
    fp2_sub(&Q->y, &Q->y, &t1);        // y2 = B^{-1}*alpha(4x1y1^2 - x2) - 8y1^4
}

void
ADD_twist(jac_point_t *R, jac_point_t const *P, jac_point_t const *Q, ec_curve_t const *AC)
{ // Addition on a Twisted Montgomery curve By^2 = x^3 + Ax^2 + x
    fp2_t t0, t1, t2, t3, t4, t5, t6;
    jac_point_t T;

    if (is_jac_equal(P, Q)) {
        DBL_twist(R, P, AC);
        return;
    }
    jac_neg(&T, P);
    if (is_jac_equal(&T, Q)) {
        jac_init(R);
        return;
    }
    if (fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) {
        copy_jac_point(R, Q);
        return;
    } else if (fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z)) {
        copy_jac_point(R, P);
        return;
    }

    fp2_sqr(&t0, &P->z);        // t0 = z1^2
    fp2_mul(&t1, &t0, &P->z);   // t1 = z1^3
    fp2_sqr(&t2, &Q->z);        // t2 = z2^2
    fp2_mul(&t3, &t2, &Q->z);   // t3 = z2^3
    fp2_mul(&t1, &t1, &Q->y);   // t1 = y2z1^3
    fp2_mul(&t3, &t3, &P->y);   // t3 = y1z2^3
    fp2_sub(&t1, &t1, &t3);     // t1 = lambda1 = y2z1^3 - y1z2^3
    fp2_mul(&t0, &t0, &Q->x);   // t0 = x2z1^2
    fp2_mul(&t2, &t2, &P->x);   // t2 = x1z2^2
    fp2_sub(&t4, &t0, &t2);     // t4 = lambda3 = x2z1^2 - x1z2^2
    fp2_add(&t0, &t0, &t2);     // t0 = lambda2 = x2z1^2 + x1z2^2
    fp2_mul(&t5, &P->z, &Q->z); // t5 = z1z2
    fp2_mul(&R->z, &t4, &t5);   // z3 = z1z2*lambda3
    fp2_sqr(&t5, &t5);          // t5 = z1^2z2^2
    fp2_mul(&t5, &AC->A, &t5);  // t5 = Az1^2z2^2
    fp2_add(&t0, &t0, &t5);     // t0 = Az1^2z2^2 + lambda2
    fp2_sqr(&t6, &t4);          // t6 = lambda3^2
    fp2_mul(&t5, &t0, &t6);     // t5 = lambda3^2(Az1^2z2^2 + lambda2)
    
    fp2_sqr(&R->x, &t1);          // x3 = lambda1^2
    fp2_mul(&R->x, &R->x, &B_qnr); // x3 = B * lambda1^2 
    fp2_sub(&R->x, &R->x, &t5);   // x3 = B*lambda1^2 - lambda3^2(Az1^2z2^2 + lambda2)

    fp2_mul(&t3, &t3, &t4);     // t3 = y1z2^3*lambda3
    fp2_mul(&t3, &t3, &t6);     // t3 = y1z2^3*lambda3^3
    fp2_mul(&t2, &t2, &t6);     // t2 = x1z2^2*lambda3^2
    fp2_sub(&R->y, &t2, &R->x); // y3 = x1z2^2*lambda3^2 - x3
    fp2_mul(&R->y, &R->y, &t1); // y3 = lambda1(x1z2^2*lambda3^2 - x3)
    fp2_sub(&R->y, &R->y, &t3); // y3 = lambda1(x1z2^2*lambda3^2 - x3) - y1z2^3*lambda3^3
}

// scalar multiplication R <- k*P
void MUL_generic_twist(jac_point_t *R,
                const jac_point_t *P,
                const digit_t *k,
                const int nwords,
                const ec_curve_t *curve
            )
{
    digit_t k_t;
    jac_point_t R_t;

    jac_init(&R_t);
    for(int i = nwords; i > 0; i--) {
        for (int j = 0; j < RADIX; j++) {
            k_t = k[i - 1] >> (RADIX - 1 - j);
            k_t &= 0x1;
            DBL_twist_inv(&R_t, &R_t, curve);
            if (k_t == 1) {
                ADD_twist(&R_t, &R_t, P, curve);
            }
        }
    }
    copy_jac_point(R, &R_t);
}