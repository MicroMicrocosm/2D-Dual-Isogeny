#include "isog.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

static inline void A24_to_AC(ec_curve_t *curve, const ec_point_t *A24)
{
    // (A:C) = ((A+2C)*2-4C : 4C)
    fp2_add(&curve->A, &A24->x, &A24->x);
    fp2_sub(&curve->A, &curve->A, &A24->z);
    fp2_add(&curve->A, &curve->A, &curve->A);
    fp2_copy(&curve->C, &A24->z);
}

void ec_dual_chain_init(ec_dual_chain_t *dual_chain)
{
    memset(dual_chain, 0, sizeof(*dual_chain));
}

static void cache_push(ec_dual_chain_t *dual_chain, const ec_dual_params_t *step)
{
    assert(dual_chain->length < MAX_DUAL_STEPS);
    dual_chain->steps[dual_chain->length] = *step;
    dual_chain->length += 1;
}

static void xisog_2_params(ec_point_t *B, const ec_point_t *P,
                            ec_dual_params_t *step)
{
    fp2_sqr(&B->x, &P->x);
    fp2_sqr(&B->z, &P->z);
    fp2_sub(&B->x, &B->z, &B->x);
    fp2_add(&K[0].x, &P->x, &P->z);
    fp2_sub(&K[0].z, &P->x, &P->z);

    memset(step, 0, sizeof(*step));
    step->kind = DUAL_2_GENERIC;
    fp2_copy(&step->a, &P->x);
    fp2_copy(&step->b, &P->z);
}

static void xisog_4_params(ec_point_t *B, const ec_point_t *P,
                            ec_dual_params_t *step)
{
    fp2_t x2, z2;

    fp2_sqr(&x2, &P->x);
    fp2_sqr(&z2, &P->z);
    fp2_copy(&K[0].x, &x2);
    fp2_copy(&K[0].z, &z2);

    fp2_sqr(&B->z, &z2);
    fp2_sqr(&B->x, &x2);
    fp2_sub(&B->x, &B->z, &B->x);

    fp2_add(&K[2].x, &P->x, &P->z);
    fp2_sub(&K[1].x, &P->x, &P->z);
    fp2_add(&K[0].x, &K[0].z, &K[0].z);
    fp2_add(&K[0].x, &K[0].x, &K[0].x);

    memset(step, 0, sizeof(*step));
    step->kind = DUAL_4_GENERIC;
    fp2_copy(&step->a, &B->x);
    fp2_copy(&step->b, &B->z);
    fp2_sub(&step->c, &x2, &z2);
    fp2_add(&step->d, &x2, &z2);
    fp2_mul(&step->e, &P->x, &P->z);
    fp2_add(&step->e, &step->e, &step->e);
}

static void xisog_2_singular_params(ec_point_t *B24,
                                    const ec_point_t *A24,
                                    ec_dual_params_t *step)
{
    fp2_t sX, sZ;

    fp2_copy(&sX, &A24->x);
    fp2_sqrt(&sX);
    fp2_copy(&sZ, &A24->z);
    fp2_sqrt(&sZ);

    fp2_copy(&K[0].x, &sX);
    fp2_copy(&K[0].z, &sZ);

    fp2_add(&B24->x, &sX, &sZ);
    fp2_sqr(&B24->x, &B24->x);

    fp2_mul(&B24->z, &sX, &sZ);
    fp2_add(&B24->z, &B24->z, &B24->z);
    fp2_add(&B24->z, &B24->z, &B24->z);

    memset(step, 0, sizeof(*step));
    step->kind = DUAL_2_SINGULAR;
    fp2_copy(&step->a, &sX);
    fp2_copy(&step->b, &sZ);
}

static void xisog_4_singular_params(ec_point_t *B24, const ec_point_t *P,
                                     const ec_point_t *A24,
                                     ec_dual_params_t *step)
{
    int sign_positive = 0;

    fp2_copy(&K[0].z, &A24->z);
    if (fp2_is_equal(&P->x, &P->z)) {
        sign_positive = 1;
        fp2_copy(&K[0].x, &A24->x);
        fp2_sub(&K[1].x, &A24->x, &A24->z);
        fp2_neg(&B24->z, &K[1].x);
    } else {
        fp2_copy(&K[1].x, &A24->x);
        fp2_sub(&K[0].x, &A24->x, &A24->z);
        fp2_neg(&B24->z, &K[0].x);
        fp2_copy(&B24->z, &K[1].x);
    }
    fp2_copy(&B24->x, &K[0].z);

    memset(step, 0, sizeof(*step));
    step->kind = DUAL_4_SINGULAR;
    step->sign_positive = (uint8_t)sign_positive;
    fp2_copy(&step->a, &B24->x);
    fp2_copy(&step->b, &B24->z);
}

static void cache_dual4_params(ec_dual_chain_t *dual_chain, ec_point_t *A24,
                               const ec_point_t *kernel, const bool singular)
{
    ec_dual_params_t step;

    if (singular)
        xisog_4_singular_params(A24, kernel, A24, &step);
    else
        xisog_4_params(A24, kernel, &step);

    cache_push(dual_chain, &step);
}

static void cache_dual2_params(ec_dual_chain_t *dual_chain, ec_point_t *A24,
                               const ec_point_t *kernel, const bool singular)
{
    ec_dual_params_t step;

    if (singular) {
        xisog_2_singular_params(A24, A24, &step);
    } else {
        xisog_2_params(A24, kernel, &step);
    }

    cache_push(dual_chain, &step);
}

static void ec_eval_even_params_sequential(ec_curve_t *image,
                                           ec_dual_chain_t *dual_chain,
                                           ec_point_t *points,
                                           unsigned short points_len,
                                           ec_point_t *A24,
                                           const ec_point_t *kernel,
                                           int isog_len)
{
    int i, j;
    int rem = isog_len;
    ec_point_t Q, K;

    copy_point(&Q, kernel);
    if (rem & 1) {
        copy_point(&K, &Q);
        for (i = 0; i < rem - 1; i++)
            xDBLv2(&K, &K, A24);
        bool singular = fp2_is_zero(&K.x);
        cache_dual2_params(dual_chain, A24, &K, singular);
        if (singular) {
            xeval_2_singular(&Q, &Q, 1);
            xeval_2_singular(points, points, points_len);
        } else {
            xeval_2(&Q, &Q, 1);
            xeval_2(points, points, points_len);
        }
        rem -= 1;
    }

    for (j = 0; j < rem / 2; j++) {
        copy_point(&K, &Q);
        for (i = 0; i < rem - 2 * j - 2; i++)
            xDBLv2(&K, &K, A24);
        cache_dual4_params(dual_chain, A24, &K, 0);
        if (j < rem / 2 - 1)
            xeval_4(&Q, &Q, 1);
        xeval_4(points, points, points_len);
    }

    A24_to_AC(image, A24);
    copy_point(&image->A24, A24);
}

static void ec_eval_even_params_strategy(ec_curve_t *image,
                                         ec_dual_chain_t *dual_chain,
                                         ec_point_t *points,
                                         unsigned short points_len,
                                         ec_point_t *A24,
                                         const ec_point_t *kernel,
                                         int isog_len)
{
    if (isog_len >= 0 && isog_len != POWER_OF_2 - 2) {
        ec_eval_even_params_sequential(image, dual_chain, points, points_len, A24,
                                       kernel, isog_len);
        return;
    }

    assert(isog_len == POWER_OF_2 - 2);
    uint8_t log2_of_e, tmp;
    uint32_t e_half = (uint32_t)(isog_len >> 1);
    for (tmp = (uint8_t)e_half, log2_of_e = 0; tmp > 0;
         tmp >>= 1, ++log2_of_e) {
    }
    log2_of_e *= 2;

    ec_point_t SPLITTING_POINTS[log2_of_e];
    int XDBLs[log2_of_e];
    copy_point(&SPLITTING_POINTS[0], kernel);

    int strategy = 0;
    int BLOCK = 0;
    int current = 0;
    int i, j;

    if (isog_len & 1) {
        copy_point(&SPLITTING_POINTS[1], &SPLITTING_POINTS[0]);
        for (i = 0; i < isog_len - 1; i++)
            xDBLv2(&SPLITTING_POINTS[1], &SPLITTING_POINTS[1], A24);
        assert(!fp2_is_zero(&SPLITTING_POINTS[1].x));
        cache_dual2_params(dual_chain, A24, &SPLITTING_POINTS[1], 0);
        xeval_2(SPLITTING_POINTS, SPLITTING_POINTS, 1);
        xeval_2(points, points, points_len);
    }

    for (j = 0; j < (int)e_half - 1; j++) {
        while (BLOCK != (int)e_half - 1 - j) {
            current += 1;
            copy_point(&SPLITTING_POINTS[current],
                       &SPLITTING_POINTS[current - 1]);
            for (i = 0; i < 2 * STRATEGY4[strategy]; i++)
                xDBLv2(&SPLITTING_POINTS[current],
                       &SPLITTING_POINTS[current], A24);
            XDBLs[current] = STRATEGY4[strategy];
            BLOCK += STRATEGY4[strategy];
            strategy += 1;
        }

        cache_dual4_params(dual_chain, A24, &SPLITTING_POINTS[current], 0);
        xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current);
        xeval_4(points, points, points_len);

        BLOCK -= XDBLs[current];
        XDBLs[current] = 0;
        current -= 1;
    }

    cache_dual4_params(dual_chain, A24, &SPLITTING_POINTS[current], 0);
    xeval_4(points, points, points_len);

    A24_to_AC(image, A24);
    copy_point(&image->A24, A24);
}

void ec_eval_even_params(ec_curve_t *image, ec_dual_chain_t *dual_chain,
                         const ec_isog_even_t *phi, ec_point_t *points,
                         unsigned short points_len)
{
    ec_point_t Q, A24;
    assert(dual_chain != NULL);
    assert(phi->length <= POWER_OF_2);
    dual_chain->length = 0;

    AC_to_A24(&A24, &phi->curve);
    copy_point(&Q, &phi->kernel);

    if (phi->length < 2) {
        ec_eval_even_params_strategy(image, dual_chain, points, points_len, &A24,
                                     &Q, (int)phi->length);
        return;
    }

    ec_point_t Q4, K2;
    copy_point(&Q4, &Q);
    for (int i = 0; i < (int)phi->length - 2; i++)
        xDBLv2(&Q4, &Q4, &A24);
    xDBLv2(&K2, &Q4, &A24);
    bool singular = fp2_is_zero(&K2.x);
    cache_dual4_params(dual_chain, &A24, &Q4, singular);
    if (singular) {
        xeval_4_singular(points, points, points_len, Q4);
        xeval_4_singular(&Q, &phi->kernel, 1, Q4);
    } else {
        xeval_4(points, points, points_len);
        xeval_4(&Q, &phi->kernel, 1);
    }

    ec_eval_even_params_strategy(image, dual_chain, points, points_len, &A24, &Q,
                                 (int)phi->length - 2);
}

void xeval_2_dual(ec_point_t *out, const ec_dual_params_t *step,
                  const ec_point_t *in)
{
    fp2_t t0, t1;

    ec_point_t tmp;
    copy_point(&tmp, in);

    if (step->kind == DUAL_2_SINGULAR) {
        /* g-hat(X:Z) = (X*(sX*X+sZ*Z) : Z*(sZ*X+sX*Z)). */
        fp2_mul(&t0, &step->a, &tmp.x);
        fp2_mul(&t1, &step->b, &tmp.z);
        fp2_add(&t0, &t0, &t1);
        fp2_mul(&out->x, &tmp.x, &t0);

        fp2_mul(&t0, &step->b, &tmp.x);
        fp2_mul(&t1, &step->a, &tmp.z);
        fp2_add(&t0, &t0, &t1);
        fp2_mul(&out->z, &tmp.z, &t0);
    } else {
        fp2_add(&t0, &tmp.x, &tmp.z);
        fp2_sqr(&t0, &t0);
        fp2_mul(&out->x, &step->b, &t0);

        fp2_mul(&t1, &tmp.x, &tmp.z);
        fp2_mul(&out->z, &step->a, &t1);
        fp2_add(&out->z, &out->z, &out->z);
        fp2_add(&out->z, &out->z, &out->z);
    }
}

void xeval_4_dual(ec_point_t *out, const ec_dual_params_t *step,
                  const ec_point_t *in)
{
    fp2_t t0, t1, t2, t3, t4, t5;
    fp2_t X;

    ec_point_t tmp;
    copy_point(&tmp, in);

    if (step->kind != DUAL_4_GENERIC) {
        /* x([2]K)=0, so x(K)=±1.
         * For x(K)=1 the dual 4-isogeny is
         * φ̂4(x) = (x+1)²((x-1)²+(A'+2)x)/((2-A')x(x-1)²),
         * with A' the codomain parameter.
         * In A24 coordinates A'=4d-2 and the projectivized form is:
         *   X = (X+Z)²(step.b·(X-Z)² + 4·step.a·XZ)
         *   Z = 4(step.b-step.a)·XZ·(X-Z)² */
        if (step->sign_positive)
            fp2_neg(&X, &tmp.x);              /* x=1: iota on input */
        else
            fp2_copy(&X, &tmp.x);             /* x=-1: no input iota */

        fp2_add(&t0, &X, &tmp.z);
        fp2_sub(&t1, &X, &tmp.z);
        fp2_sqr(&t0, &t0);
        fp2_sqr(&t1, &t1);
        fp2_mul(&t2, &X, &tmp.z);

        if (step->sign_positive) {
            /* x=1: X=(X+Z)^2*(step.b*(X-Z)^2+4*(step.b-step.a)*XZ),
             *      Z=4*step.a*XZ*(X-Z)^2 */
            fp2_sub(&t3, &step->b, &step->a);
            fp2_mul(&t4, &t3, &t2);
            fp2_add(&t4, &t4, &t4);
            fp2_add(&t4, &t4, &t4);
            fp2_mul(&t3, &step->b, &t1);
            fp2_add(&t3, &t3, &t4);
            fp2_mul(&out->x, &t0, &t3);

            fp2_mul(&t4, &step->a, &t2);
            fp2_add(&t4, &t4, &t4);
            fp2_add(&t4, &t4, &t4);
            fp2_mul(&out->z, &t4, &t1);
        } else {
            /* x=-1: dual = C's forward formula with A24 replaced by codomain.
             * X = (X-Z)^2*(step.b*(X-Z)^2 + 4*step.a*XZ)
             * Z = 4*step.a*XZ*(X+Z)^2 */
            fp2_mul(&t3, &step->b, &t1);
            fp2_mul(&t4, &step->a, &t2);
            fp2_add(&t4, &t4, &t4);
            fp2_add(&t4, &t4, &t4);
            fp2_add(&t3, &t3, &t4);
            fp2_mul(&out->x, &t1, &t3);

            fp2_mul(&t4, &step->a, &t2);
            fp2_add(&t4, &t4, &t4);
            fp2_add(&t4, &t4, &t4);
            fp2_mul(&out->z, &t4, &t0);
        }
        return;
    }

    /* C's xisog_4 emits codomain A24=(Z^4-X^4:Z^4), while the paper's φ2
     * uses A24=(X^4:Z^4). They differ by the isomorphism ι (A -> -A), which
     * on the Kummer line is x -> -x; apply it to the input x-coordinate. */
    fp2_neg(&X, &tmp.x);

    /* Input (X:Z), constants:
     *   C=step.b (codomain A24.z), T3=4XZ*step.a, where step.a equals the
     *   paper codomain's ĉ24-â24. K0=step.c, K1=step.d, K2=step.e.
     * T0=(X+Z)^2, T1=(X-Z)^2, T2=T0-T1=4XZ, T3=T2*step.a.
     * X := T0*(C*T0 - T3), Z := T3*T1.
     * out.x = K0*X + K1*Z, out.z = K2*Z. */
    fp2_add(&t0, &X, &tmp.z);             /* X+Z */
    fp2_sub(&t1, &X, &tmp.z);             /* X-Z */
    fp2_sqr(&t0, &t0);                    /* (X+Z)^2 */
    fp2_sqr(&t1, &t1);                    /* (X-Z)^2 */
    fp2_sub(&t2, &t0, &t1);               /* 4XZ */
    fp2_mul(&t3, &t2, &step->a);          /* T3 = 4XZ*(Z^4-X^4) */

    fp2_mul(&t4, &step->b, &t0);          /* C*(X+Z)^2 */
    fp2_sub(&t4, &t4, &t3);               /* C*(X+Z)^2 - T3 */
    fp2_mul(&t4, &t0, &t4);               /* X */

    fp2_mul(&t5, &t3, &t1);               /* Z = T3*(X-Z)^2 */

    fp2_mul(&t0, &step->c, &t4);          /* K0*X */
    fp2_mul(&t1, &step->d, &t5);          /* K1*Z */
    fp2_add(&out->x, &t0, &t1);
    fp2_mul(&out->z, &step->e, &t5);      /* K2*Z */
}

static void xeval_dual(ec_point_t *out, const ec_dual_chain_t *dual_chain, const ec_point_t *K_prime)
{
    ec_point_t cur, next;

    copy_point(&cur, K_prime);
    for (int i = (int)dual_chain->length - 1; i >= 0; i--) {
        const ec_dual_params_t *step = &dual_chain->steps[i];

        if (step->kind == DUAL_2_GENERIC || step->kind == DUAL_2_SINGULAR) {
            xeval_2_dual(&next, step, &cur);
        } else {
            xeval_4_dual(&next, step, &cur);
        }
        copy_point(&cur, &next);
    }
    copy_point(out, &cur);
}

void ec_eval_even_dual(const ec_dual_chain_t *dual_chain,
                       ec_point_t *points,
                       unsigned short points_len)
{
    for (int i = 0; i < points_len; i++) {
        xeval_dual(&points[i], dual_chain, &points[i]);
    }
}

// void ec_eval_even_dual(const ec_dual_chain_t *dual_chain,
//                        ec_point_t *points,
//                        unsigned short points_len)
// {
//     for (int i = (int)dual_chain->length - 1; i >= 0; i--) {
//         const ec_dual_params_t *step = &dual_chain->steps[i];
// 
//         for (int j = 0; j < points_len; j++) {
//             if (step->kind == DUAL_2_GENERIC || step->kind == DUAL_2_SINGULAR) {
//                 xeval_2_dual(&points[j], step, &points[j]);
//             } else {
//                 xeval_4_dual(&points[j], step, &points[j]);
//             }
//         }
//     }
// }
