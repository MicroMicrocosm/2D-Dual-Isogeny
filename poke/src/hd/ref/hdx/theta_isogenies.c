#include "theta_isogenies.h"
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <tools.h>

#ifndef NDEBUG
#define DEBUG_PRINT 1
#else
#define DEBUG_PRINT 0
#endif

#ifndef __FILE_NAME__
#define __FILE_NAME__ "NA"
#endif

#ifndef __LINE__
#define __LINE__ 0
#endif

#ifndef __func__
#define __func__ "NA"
#endif

//debug print
#define debug_print(fmt)                                                                           \
    do {                                                                                           \
        if (DEBUG_PRINT)                                                                           \
            printf("warning: %s, file %s, line %d, function %s().\n",                              \
                   fmt,                                                                            \
                   __FILE_NAME__,                                                                  \
                   __LINE__,                                                                       \
                   __func__);                                                                      \
    } while (0)

void
theta_print(char *name, theta_point_t P)
{
    // fp2_t a;
    // assert(!fp2_is_zero(&P.y));
    // fp2_copy(&a, &P.y);
    // fp2_inv(&a);
    // fp2_mul(&a, &a, &P.t);
    // fp2_print(name, &a);
    printf("%s\n", name);
    fp2_print("x = ", &P.x);
    fp2_print("y = ", &P.y);
    fp2_print("z = ", &P.z);
    fp2_print("t = ", &P.t);
    printf("\n");
}

void
fp2_sett(fp2_t *a, int t)
{
    if (t == 1) {
        fp2_set_one(a);
    } else if (t == 0) {
        fp2_set_zero(a);
    } else if (t == -1) {
        fp2_set_one(a);
        fp2_neg(a, a);
    } else {
        fp2_set_small(a, t);
    }
}

void
choose_index_theta_point(fp2_t *res, int ind, const theta_point_t *T)
{
    int t = ind % 4;
    if (t == 0) {
        fp2_copy(res, &T->x);
    } else if (t == 1) {
        fp2_copy(res, &T->y);
    } else if (t == 2) {
        fp2_copy(res, &T->z);
    } else if (t == 3) {
        fp2_copy(res, &T->t);
    } else {
        assert(0);
    }
}

void
set_index_theta_point(theta_point_t *res, int ind, const fp2_t *val)
{
    int t = ind % 4;
    if (t == 0) {
        fp2_copy(&res->x, val);
    } else if (t == 1) {
        fp2_copy(&res->y, val);
    } else if (t == 2) {
        fp2_copy(&res->z, val);
    } else if (t == 3) {
        fp2_copy(&res->t, val);
    } else {
        assert(0);
    }
}

void
get_matrix(fp2_t *a00, fp2_t *a01, fp2_t *a10, fp2_t *a11, const ec_point_t *P, const ec_curve_t *E)
{
    ec_point_t P2;
    fp2_t det, u, invz, temp;

    // temp = [2]P
    ec_dbl(&P2, E, P);

    // invz = 1/P.z
    invz = P->z;
    fp2_inv(&invz);

    // P = (x,z) P2 = (u,w)
    // det = x w - u z
    fp2_mul(&det, &P->x, &P2.z);
    fp2_mul(&temp, &P->z, &P2.x);
    fp2_sub(&det, &det, &temp);

    // det = 1/det
    fp2_inv(&det);

    // a10 = ux /det - x/z
    fp2_mul(&temp, &P->x, &invz);
    fp2_mul(a10, &P->x, &P2.x);
    fp2_mul(a10, a10, &det);
    fp2_sub(a10, a10, &temp);

    // a11 = u z * det
    fp2_mul(a11, &P2.x, &det);
    fp2_mul(a11, a11, &P->z);

    // a00 = -a11
    fp2_neg(a00, a11);

    // a01 = - w z det
    fp2_mul(a01, &P2.z, &det);
    fp2_mul(a01, a01, &P->z);
    fp2_neg(a01, a01);
}

// compute the theta_point corresponding to the couple of point T on an elliptic product
void
base_change(theta_point_t *out, const theta_gluing_t *phi, const theta_couple_point_t *T)
{
    fp2_t a, b, c, d, x1, x2;
    if (fp2_is_zero(&T->P1.z) && fp2_is_zero(&T->P1.x)) {
        fp2_sett(&x1, 1);
    } else {
        x1 = T->P1.x;
    }
    if (fp2_is_zero(&T->P2.z) && fp2_is_zero(&T->P2.x)) {
        fp2_sett(&x2, 1);
    } else {
        x2 = T->P2.x;
    }

    // a = P1.x P2.x, b = P1.x P2.z, c =P1.z P2.x, d = P1.z P2.z
    fp2_mul(&a, &x1, &x2);
    fp2_mul(&b, &x1, &T->P2.z);
    fp2_mul(&c, &x2, &T->P1.z);
    fp2_mul(&d, &T->P1.z, &T->P2.z);

    // Apply the matrix
    fp2_mul(&out->x, &a, &phi->M00);
    fp2_mul(&x1, &b, &phi->M01);
    fp2_add(&out->x, &out->x, &x1);
    fp2_mul(&x1, &c, &phi->M02);
    fp2_add(&out->x, &out->x, &x1);
    fp2_mul(&x1, &d, &phi->M03);
    fp2_add(&out->x, &out->x, &x1);

    fp2_mul(&out->y, &a, &phi->M10);
    fp2_mul(&x1, &b, &phi->M11);
    fp2_add(&out->y, &out->y, &x1);
    fp2_mul(&x1, &c, &phi->M12);
    fp2_add(&out->y, &out->y, &x1);
    fp2_mul(&x1, &d, &phi->M13);
    fp2_add(&out->y, &out->y, &x1);

    fp2_mul(&out->z, &a, &phi->M20);
    fp2_mul(&x1, &b, &phi->M21);
    fp2_add(&out->z, &out->z, &x1);
    fp2_mul(&x1, &c, &phi->M22);
    fp2_add(&out->z, &out->z, &x1);
    fp2_mul(&x1, &d, &phi->M23);
    fp2_add(&out->z, &out->z, &x1);

    fp2_mul(&out->t, &a, &phi->M30);
    fp2_mul(&x1, &b, &phi->M31);
    fp2_add(&out->t, &out->t, &x1);
    fp2_mul(&x1, &c, &phi->M32);
    fp2_add(&out->t, &out->t, &x1);
    fp2_mul(&x1, &d, &phi->M33);
    fp2_add(&out->t, &out->t, &x1);
}

static int
get_symmetric_element(
    coherent_symmetric_element_t *out,
    const ec_point_t *P1,
    const ec_point_t *P2)
{
    // Assume zeta = sqrt(-1)
    fp2_t x1, z1, x2, z2;

    fp2_copy(&x1, &P1->x);
    fp2_copy(&x2, &P2->x);
    fp2_copy(&z1, &P1->z);
    fp2_copy(&z2, &P2->z);

    fp2_t xz_p1, xz_p2, xz_m1, xz_m2, xz_pm12, xz_mp12, tmp, xz_sqr, xz_dbl;
    fp2_add(&xz_p1, &x1, &z1);
    fp2_add(&xz_p2, &x2, &z2);
    fp2_sub(&xz_m1, &x1, &z1);
    fp2_sub(&xz_m2, &x2, &z2);
    fp2_mul(&xz_pm12, &xz_p1, &xz_m2);
    fp2_mul(&xz_mp12, &xz_m1, &xz_p2);
    fp2_neg(&tmp, &xz_mp12);

    // xz_pm12   = ( a,  b) = a + bi
    // xz_mp12   = ( c,  d)
    // tmp       = (-c, -d)
    // i*xz_pm12 = (-b,  a) = ai - b

    fp2_t zeta;
    fp2_set_zero(&zeta);
    fp_set_one(&zeta.im);

    // Determine the position and compute for each case
    // TODO Constant time, no if-else statement
    if (fp2_is_zero(&xz_m1)) {
        out->pos = 0; // x1 = z1
        fp2_mul(&out->lambda, &xz_p2, &xz_m2);

        fp2_set_zero(&out->gT1[0][0]);
        fp2_copy(&out->gT1[0][1], &out->lambda);
        fp2_copy(&out->gT1[1][0], &out->lambda);
        fp2_set_zero(&out->gT1[1][1]);

        fp2_sqr(&xz_sqr, &x2);
        fp2_sqr(&xz_dbl, &z2);
        fp2_add(&xz_sqr, &xz_sqr, &xz_dbl);
        fp2_mul(&xz_dbl, &x2, &z2);
        fp2_add(&xz_dbl, &xz_dbl, &xz_dbl);
        fp2_copy(&out->gT2[0][0], &xz_sqr);
        fp2_neg(&out->gT2[0][1], &xz_dbl);
        fp2_copy(&out->gT2[1][0], &xz_dbl);
        fp2_neg(&out->gT2[1][1], &xz_sqr);

        fp2_mul(&xz_sqr, &xz_sqr, &zeta);
        fp2_mul(&xz_dbl, &xz_dbl, &zeta);
        fp2_copy(&out->gT3[0][0], &xz_dbl);
        fp2_neg(&out->gT3[0][1], &xz_sqr);
        fp2_copy(&out->gT3[1][0], &xz_sqr);
        fp2_neg(&out->gT3[1][1], &xz_dbl);
    } else if (fp2_is_zero(&xz_p1)) {
        out->pos = 1; // x1 = -z1
        fp2_mul(&out->lambda, &xz_p2, &xz_m2);

        fp2_set_zero(&out->gT1[0][0]);
        fp2_neg(&out->gT1[0][1], &out->lambda);
        fp2_neg(&out->gT1[1][0], &out->lambda);
        fp2_set_zero(&out->gT1[1][1]);

        fp2_sqr(&xz_sqr, &x2);
        fp2_sqr(&xz_dbl, &z2);
        fp2_add(&xz_sqr, &xz_sqr, &xz_dbl);
        fp2_mul(&xz_dbl, &x2, &z2);
        fp2_add(&xz_dbl, &xz_dbl, &xz_dbl);
        fp2_copy(&out->gT2[0][0], &xz_sqr);
        fp2_neg(&out->gT2[0][1], &xz_dbl);
        fp2_copy(&out->gT2[1][0], &xz_dbl);
        fp2_neg(&out->gT2[1][1], &xz_sqr);

        fp2_mul(&xz_sqr, &xz_sqr, &zeta);
        fp2_mul(&xz_dbl, &xz_dbl, &zeta);
        fp2_neg(&out->gT3[0][0], &xz_dbl);
        fp2_copy(&out->gT3[0][1], &xz_sqr);
        fp2_neg(&out->gT3[1][0], &xz_sqr);
        fp2_copy(&out->gT3[1][1], &xz_dbl);
    } else if (fp2_is_zero(&xz_m2)) {
        out->pos = 2; // x2 = z2
        fp2_mul(&out->lambda, &xz_p1, &xz_m1);

        fp2_sqr(&xz_sqr, &x1);
        fp2_sqr(&xz_dbl, &z1);
        fp2_add(&xz_sqr, &xz_sqr, &xz_dbl);
        fp2_mul(&xz_dbl, &x1, &z1);
        fp2_add(&xz_dbl, &xz_dbl, &xz_dbl);
        fp2_copy(&out->gT1[0][0], &xz_sqr);
        fp2_neg(&out->gT1[0][1], &xz_dbl);
        fp2_copy(&out->gT1[1][0], &xz_dbl);
        fp2_neg(&out->gT1[1][1], &xz_sqr);

        fp2_set_zero(&out->gT2[0][0]);
        fp2_copy(&out->gT2[0][1], &out->lambda);
        fp2_copy(&out->gT2[1][0], &out->lambda);
        fp2_set_zero(&out->gT2[1][1]);

        fp2_mul(&xz_sqr, &xz_sqr, &zeta);
        fp2_mul(&xz_dbl, &xz_dbl, &zeta);
        fp2_neg(&out->gT3[0][0], &xz_dbl);
        fp2_copy(&out->gT3[0][1], &xz_sqr);
        fp2_neg(&out->gT3[1][0], &xz_sqr);
        fp2_copy(&out->gT3[1][1], &xz_dbl);
    } else if (fp2_is_zero(&xz_p2)) {
        out->pos = 3; // x2 = -z2
        fp2_mul(&out->lambda, &xz_p1, &xz_m1);

        fp2_sqr(&xz_sqr, &x1);
        fp2_sqr(&xz_dbl, &z1);
        fp2_add(&xz_sqr, &xz_sqr, &xz_dbl);
        fp2_mul(&xz_dbl, &x1, &z1);
        fp2_add(&xz_dbl, &xz_dbl, &xz_dbl);
        fp2_copy(&out->gT1[0][0], &xz_sqr);
        fp2_neg(&out->gT1[0][1], &xz_dbl);
        fp2_copy(&out->gT1[1][0], &xz_dbl);
        fp2_neg(&out->gT1[1][1], &xz_sqr);

        fp2_set_zero(&out->gT2[0][0]);
        fp2_neg(&out->gT2[0][1], &out->lambda);
        fp2_neg(&out->gT2[1][0], &out->lambda);
        fp2_set_zero(&out->gT2[1][1]);

        fp2_mul(&xz_sqr, &xz_sqr, &zeta);
        fp2_mul(&xz_dbl, &xz_dbl, &zeta);
        fp2_copy(&out->gT3[0][0], &xz_dbl);
        fp2_neg(&out->gT3[0][1], &xz_sqr);
        fp2_copy(&out->gT3[1][0], &xz_sqr);
        fp2_neg(&out->gT3[1][1], &xz_dbl);
    } else if (fp_is_equal(&xz_pm12.im, &xz_mp12.re) & fp_is_equal(&xz_pm12.re, &tmp.im)) {
        out->pos = 4; // i*xz_pm12 = -xz_mp12, b = c, a = -d
        fp2_mul(&out->lambda, &xz_p1, &xz_m1);

        fp2_sqr(&xz_sqr, &x1);
        fp2_sqr(&xz_dbl, &z1);
        fp2_add(&xz_sqr, &xz_sqr, &xz_dbl);
        fp2_mul(&xz_dbl, &x1, &z1);
        fp2_add(&xz_dbl, &xz_dbl, &xz_dbl);
        fp2_copy(&out->gT1[0][0], &xz_sqr);
        fp2_neg(&out->gT1[0][1], &xz_dbl);
        fp2_copy(&out->gT1[1][0], &xz_dbl);
        fp2_neg(&out->gT1[1][1], &xz_sqr);

        fp2_mul(&xz_sqr, &xz_sqr, &zeta);
        fp2_mul(&xz_dbl, &xz_dbl, &zeta);
        fp2_neg(&out->gT2[0][0], &xz_dbl);
        fp2_copy(&out->gT2[0][1], &xz_sqr);
        fp2_neg(&out->gT2[1][0], &xz_sqr);
        fp2_copy(&out->gT2[1][1], &xz_dbl);

        fp2_set_zero(&out->gT3[0][0]);
        fp2_neg(&out->gT3[0][1], &out->lambda);
        fp2_neg(&out->gT3[1][0], &out->lambda);
        fp2_set_zero(&out->gT3[1][1]);
    } else if (fp_is_equal(&xz_pm12.im, &tmp.re) & fp_is_equal(&xz_pm12.re, &xz_mp12.im)) {
        out->pos = 5; // i*xz_pm12 = xz_mp12, b = -c, a = d
        fp2_mul(&out->lambda, &xz_p1, &xz_m1);

        fp2_sqr(&xz_sqr, &x1);
        fp2_sqr(&xz_dbl, &z1);
        fp2_add(&xz_sqr, &xz_sqr, &xz_dbl);
        fp2_mul(&xz_dbl, &x1, &z1);
        fp2_add(&xz_dbl, &xz_dbl, &xz_dbl);
        fp2_copy(&out->gT1[0][0], &xz_sqr);
        fp2_neg(&out->gT1[0][1], &xz_dbl);
        fp2_copy(&out->gT1[1][0], &xz_dbl);
        fp2_neg(&out->gT1[1][1], &xz_sqr);

        fp2_mul(&xz_sqr, &xz_sqr, &zeta);
        fp2_mul(&xz_dbl, &xz_dbl, &zeta);
        fp2_copy(&out->gT2[0][0], &xz_dbl);
        fp2_neg(&out->gT2[0][1], &xz_sqr);
        fp2_copy(&out->gT2[1][0], &xz_sqr);
        fp2_neg(&out->gT2[1][1], &xz_dbl);

        fp2_set_zero(&out->gT3[0][0]);
        fp2_copy(&out->gT3[0][1], &out->lambda);
        fp2_copy(&out->gT3[1][0], &out->lambda);
        fp2_set_zero(&out->gT3[1][1]);
    }
    else
    {
        debug_print("no matching case found.");
        return 0;
    }

    return 1;
}

void
get_superglue_coefficients(
    theta_gluing_t *out,
    const theta_couple_point_t *K1_4,
    const theta_couple_point_t *K2_4)
{
    coherent_symmetric_element_t symm_E1, symm_E2;

    get_symmetric_element(&symm_E1, &K1_4->P1, &K2_4->P1);
    get_symmetric_element(&symm_E2, &K1_4->P2, &K2_4->P2);

    out->coeff.pos = 3 * (symm_E1.pos >> 1) + (symm_E2.pos >> 1);

    bool sgl1 = (out->coeff.pos == 0);
    bool sgl2 = (out->coeff.pos == 4);
    bool sgl3 = (out->coeff.pos == 8);
    bool nsgl = !(sgl1 | sgl2 | sgl3);

    fp2_t lambda;
    fp2_mul(&lambda, &symm_E1.lambda, &symm_E2.lambda);

    fp2_t M00, M01, M02, M03;
    fp2_t M10, M11, M12, M13;
    fp2_t M20, M21, M22, M23;
    fp2_t M30, M31, M32, M33;
    // fp2_t Mi1;
    fp2_t tp0, tp1, tp2;

    // The trace with id
    // fp2_copy(&M00, &lambda);
    // fp2_set_zero(&M01);
    // fp2_set_zero(&M02);
    // fp2_set_zero(&M03);

    // gT1
    fp2_mul(&tp0, &symm_E1.gT1[0][0], &symm_E2.gT1[0][0]);
    fp2_add(&M00, &lambda, &tp0);
    fp2_mul(&M01, &symm_E1.gT1[0][0], &symm_E2.gT1[0][1]);
    fp2_mul(&M02, &symm_E1.gT1[0][1], &symm_E2.gT1[0][0]);
    fp2_mul(&M03, &symm_E1.gT1[0][1], &symm_E2.gT1[0][1]);

    // gT2
    fp2_mul(&tp0, &symm_E1.gT2[0][0], &symm_E2.gT2[0][0]);
    fp2_add(&M00, &M00, &tp0);
    fp2_mul(&tp0, &symm_E1.gT2[0][0], &symm_E2.gT2[0][1]);
    fp2_add(&M01, &M01, &tp0);
    fp2_mul(&tp0, &symm_E1.gT2[0][1], &symm_E2.gT2[0][0]);
    fp2_add(&M02, &M02, &tp0);
    fp2_mul(&tp0, &symm_E1.gT2[0][1], &symm_E2.gT2[0][1]);
    fp2_add(&M03, &M03, &tp0);

    // gT1 gT2
    fp2_mul(&tp0, &symm_E1.gT3[0][0], &symm_E2.gT3[0][0]);
    fp2_sub(&M00, &M00, &tp0);
    fp2_mul(&tp0, &symm_E1.gT3[0][0], &symm_E2.gT3[0][1]);
    fp2_sub(&M01, &M01, &tp0);
    fp2_mul(&tp0, &symm_E1.gT3[0][1], &symm_E2.gT3[0][0]);
    fp2_sub(&M02, &M02, &tp0);
    fp2_mul(&tp0, &symm_E1.gT3[0][1], &symm_E2.gT3[0][1]);
    fp2_sub(&M03, &M03, &tp0);

    // Act by S1' = (0, Q2')
    fp2_mul(&M10, &M00, &symm_E2.gT2[0][0]);
    fp2_mul(&tp0, &M01, &symm_E2.gT2[1][0]);
    fp2_add(&M10, &M10, &tp0);
    fp2_mul(&M11, &M00, &symm_E2.gT2[0][1]);
    fp2_mul(&tp0, &M01, &symm_E2.gT2[1][1]);
    fp2_add(&M11, &M11, &tp0);
    fp2_mul(&M12, &M02, &symm_E2.gT2[0][0]);
    fp2_mul(&tp0, &M03, &symm_E2.gT2[1][0]);
    fp2_add(&M12, &M12, &tp0);
    fp2_mul(&M13, &M02, &symm_E2.gT2[0][1]);
    fp2_mul(&tp0, &M03, &symm_E2.gT2[1][1]);
    fp2_add(&M13, &M13, &tp0);

    // Act by S2' = (P1', 0)
    fp2_mul(&M20, &M00, &symm_E1.gT1[0][0]);
    fp2_mul(&tp0, &M02, &symm_E1.gT1[1][0]);
    fp2_add(&M20, &M20, &tp0);
    fp2_mul(&M21, &M01, &symm_E1.gT1[0][0]);
    fp2_mul(&tp0, &M03, &symm_E1.gT1[1][0]);
    fp2_add(&M21, &M21, &tp0);
    fp2_mul(&M22, &M00, &symm_E1.gT1[0][1]);
    fp2_mul(&tp0, &M02, &symm_E1.gT1[1][1]);
    fp2_add(&M22, &M22, &tp0);
    fp2_mul(&M23, &M01, &symm_E1.gT1[0][1]);
    fp2_mul(&tp0, &M03, &symm_E1.gT1[1][1]);
    fp2_add(&M23, &M23, &tp0);

    // Act by S1' + S2' = (P1', Q2')
    fp2_mul(&M30, &M10, &symm_E1.gT1[0][0]);
    fp2_mul(&tp0, &M12, &symm_E1.gT1[1][0]);
    fp2_add(&M30, &M30, &tp0);
    fp2_mul(&M31, &M11, &symm_E1.gT1[0][0]);
    fp2_mul(&tp0, &M13, &symm_E1.gT1[1][0]);
    fp2_add(&M31, &M31, &tp0);
    fp2_mul(&M32, &M10, &symm_E1.gT1[0][1]);
    fp2_mul(&tp0, &M12, &symm_E1.gT1[1][1]);
    fp2_add(&M32, &M32, &tp0);
    fp2_mul(&M33, &M11, &symm_E1.gT1[0][1]);
    fp2_mul(&tp0, &M13, &symm_E1.gT1[1][1]);
    fp2_add(&M33, &M33, &tp0);

    fp2_mul(&M00, &M00, &lambda);
    fp2_mul(&M01, &M01, &lambda);
    fp2_mul(&M02, &M02, &lambda);
    fp2_mul(&M03, &M03, &lambda);

    fp2_mul(&M10, &M10, &symm_E1.lambda);
    fp2_mul(&M11, &M11, &symm_E1.lambda);
    fp2_mul(&M12, &M12, &symm_E1.lambda);
    fp2_mul(&M13, &M13, &symm_E1.lambda);

    fp2_mul(&M20, &M20, &symm_E2.lambda);
    fp2_mul(&M21, &M21, &symm_E2.lambda);
    fp2_mul(&M22, &M22, &symm_E2.lambda);
    fp2_mul(&M23, &M23, &symm_E2.lambda);

    fp2_copy(&out->M00, &M00);
    fp2_copy(&out->M01, &M01);
    fp2_copy(&out->M02, &M02);
    fp2_copy(&out->M03, &M03);
    fp2_copy(&out->M10, &M10);
    fp2_copy(&out->M11, &M11);
    fp2_copy(&out->M12, &M12);
    fp2_copy(&out->M13, &M13);
    fp2_copy(&out->M20, &M20);
    fp2_copy(&out->M21, &M21);
    fp2_copy(&out->M22, &M22);
    fp2_copy(&out->M23, &M23);
    fp2_copy(&out->M30, &M30);
    fp2_copy(&out->M31, &M31);
    fp2_copy(&out->M32, &M32);
    fp2_copy(&out->M33, &M33);
    
    // if (sgl1 || sgl3) {
    //     fp2_copy(&Mi1, &M11);
    //     // fp2_mul(&Mi1, &M00, &symm_E2.gT2[0][1]);
    //     // fp2_mul(&tp0, &M01, &symm_E2.gT2[1][1]);
    //     // fp2_add(&Mi1, &Mi1, &tp0);
    //     // fp2_mul(&Mi1, &Mi1, &symm_E1.lambda);
    // }
    // if (sgl2) {
    //     fp2_copy(&Mi1, &M21);
    //     // fp2_mul(&Mi1, &M01, &symm_E1.gT1[0][0]);
    //     // fp2_mul(&tp0, &M03, &symm_E1.gT1[1][0]);
    //     // fp2_add(&Mi1, &Mi1, &tp0);
    //     // fp2_mul(&Mi1, &Mi1, &symm_E2.lambda);
    // }
    // fp2_select(&Mi1, &M11, &M21, -(uint32_t)sgl2);

    // g00
    fp2_sqr(&tp0, &M00);
    fp2_sqr(&tp1, &M10);
    fp2_sqr(&tp2, &M20);
    fp2_add(&out->coeff.g00[0], &tp1, &tp2);
    fp2_sub(&out->coeff.g00[1], &tp0, &tp1);
    fp2_sub(&out->coeff.g00[2], &tp0, &tp2);

    out->coeff.sign_flip = !fp2_is_equal(&M01, &M02);
    fp2_mul(&tp0, &M00, &M01);
    switch (out->coeff.pos)
    {
    case 0: // Position 00
    case 8: // Position 22
        // g01
        // fp2_mul(&tp1, &M10, &Mi1);
        fp2_mul(&tp1, &M10, &M11);
        fp2_add(&out->coeff.g01[0], &tp0, &tp1);
        fp2_sub(&out->coeff.g01[1], &tp0, &tp1);
        // Dummy computation
        // fp2_add(&out->coeff.g02[0], &tp0, &tp1);
        // fp2_sub(&out->coeff.g02[1], &tp0, &tp1);
        // fp2_add(&out->coeff.g03[0], &tp0, &tp1);
        // fp2_sub(&out->coeff.g03[1], &tp0, &tp1);
        break;
    case 4: // Position 11
        // g01
        // fp2_mul(&tp1, &M20, &Mi1);
        fp2_mul(&tp1, &M20, &M21);
        fp2_add(&out->coeff.g01[0], &tp0, &tp1);
        fp2_sub(&out->coeff.g01[1], &tp0, &tp1);
        // Dummy computation
        // fp2_add(&out->coeff.g02[0], &tp0, &tp1);
        // fp2_sub(&out->coeff.g02[1], &tp0, &tp1);
        // fp2_add(&out->coeff.g03[0], &tp0, &tp1);
        // fp2_sub(&out->coeff.g03[1], &tp0, &tp1);
        break;
    default: // Other position
        out->coeff.sign_flip = false;
        // g01
        fp2_mul(&tp1, &M02, &M03);
        fp2_add(&out->coeff.g01[0], &tp0, &tp1);
        fp2_sub(&out->coeff.g01[1], &tp0, &tp1);
        // g02
        fp2_mul(&tp0, &M00, &M02);
        fp2_mul(&tp1, &M01, &M03);
        fp2_add(&out->coeff.g02[0], &tp0, &tp1);
        fp2_sub(&out->coeff.g02[1], &tp0, &tp1);
        // g03
        fp2_mul(&tp0, &M00, &M03);
        fp2_mul(&tp1, &M01, &M02);
        fp2_add(&out->coeff.g03[0], &tp0, &tp1);
        fp2_sub(&out->coeff.g03[1], &tp0, &tp1);
        break;
    }
}

void superglue_evaluation(
    theta_point_t *out,
    const superglue_coefficient_t *coeff,
    const add_components_t *bary1,
    const add_components_t *bary2,
    const fp2_t *twist_P,
    const fp2_t *twist_Q)
{
    fp2_t uuvv_p_ww1, uuvv_p_ww2;
    fp2_t uuvv_m_ww1, uuvv_m_ww2;
    fp2_t uw1, uw2, vw12;
    fp2_t pmw1, tmp1, tmp2, tmp3, tmp4;

    fp2_t delta;
    fp2_mul(&delta, twist_P, twist_Q);

    fp2_copy(&pmw1, &bary1->w);
    if (coeff->sign_flip) {
        fp2_neg(&pmw1, &pmw1);
    }

    fp2_sqr(&tmp1, &bary1->u);
    fp2_sqr(&tmp2, &bary1->v);
    fp2_mul(&tmp2, &tmp2, &delta);
    fp2_sub(&tmp1, &tmp1, &tmp2);
    fp2_sqr(&tmp2, &bary1->w);
    fp2_sub(&uuvv_m_ww1, &tmp1, &tmp2);
    fp2_add(&uuvv_p_ww1, &tmp1, &tmp2);

    fp2_sqr(&tmp1, &bary2->u);
    fp2_sqr(&tmp2, &bary2->v);
    fp2_mul(&tmp2, &tmp2, &delta);
    fp2_sub(&tmp1, &tmp1, &tmp2);
    fp2_sqr(&tmp2, &bary2->w);
    fp2_sub(&uuvv_m_ww2, &tmp1, &tmp2);
    fp2_add(&uuvv_p_ww2, &tmp1, &tmp2);

    fp2_mul(&uw1, &bary1->u, &pmw1);
    fp2_add(&uw1, &uw1, &uw1);

    fp2_mul(&uw2, &bary2->u, &bary2->w);
    fp2_add(&uw2, &uw2, &uw2);

    fp2_mul(&vw12, &bary1->v, &pmw1);
    fp2_mul(&tmp1, &bary2->v, &bary2->w);
    fp2_mul(&vw12, &vw12, &tmp1);
    fp2_add(&vw12, &vw12, &vw12);
    fp2_add(&vw12, &vw12, &vw12);
    fp2_mul(&vw12, &vw12, &delta);

    switch (coeff->pos)
    {
    case 1:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->y, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->z, &tmp3, &uuvv_m_ww1);

        fp2_mul(&out->t, &coeff->g03[1], &vw12);
        break;
    case 2:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->x, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->z, &tmp3, &uuvv_m_ww1);

        fp2_mul(&out->t, &coeff->g03[1], &vw12);
        break;
    case 3:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->z, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->y, &tmp3, &uuvv_m_ww1);

        fp2_mul(&out->t, &coeff->g03[1], &vw12);
        break;
    case 5:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->x, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->y, &tmp3, &uuvv_m_ww1);

        fp2_mul(&out->t, &coeff->g03[1], &vw12);
        break;
    case 6:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->z, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->x, &tmp3, &uuvv_m_ww1);

        fp2_mul(&out->t, &coeff->g03[1], &vw12);
        break;
    case 7:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->y, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->x, &tmp3, &uuvv_m_ww1);

        fp2_mul(&out->t, &coeff->g03[1], &vw12);
        break;
    case 0:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uuvv_m_ww2, &uuvv_m_ww1);
        fp2_mul(&out->z, &tmp1, &coeff->g00[2]);

        fp2_mul(&out->t, &coeff->g00[2], &vw12);
        break;
    case 4:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[2]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uuvv_m_ww2, &uuvv_m_ww1);
        fp2_mul(&out->y, &tmp1, &coeff->g00[1]);

        fp2_mul(&out->t, &coeff->g00[1], &vw12);
        break;
    case 8:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[2]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uuvv_m_ww2, &uuvv_m_ww1);
        fp2_mul(&out->x, &tmp1, &coeff->g00[0]);

        fp2_mul(&out->t, &coeff->g00[0], &vw12);
        break;
    default:
        debug_print("position should be between 0 and 8");
        break;
    }
}

// bary1 = (u1, v1^2, w1)
// bary2 = (u2, v2^2, w2)
// v1 == 0 or v2 == 0
void superglue_evaluation_half_special(
    theta_point_t *out,
    const superglue_coefficient_t *coeff,
    const add_components_t *bary1,
    const add_components_t *bary2)
{
    assert(fp2_is_zero(&bary1->v) || fp2_is_zero(&bary2->v));
    fp2_t uuvv_p_ww1, uuvv_p_ww2;
    fp2_t uuvv_m_ww1, uuvv_m_ww2;
    fp2_t uw1, uw2;
    fp2_t pmw1, tmp1, tmp2, tmp3, tmp4;

    fp2_copy(&pmw1, &bary1->w);
    if (coeff->sign_flip) {
        fp2_neg(&pmw1, &pmw1);
    }

    fp2_sqr(&tmp1, &bary1->u);
    fp2_sub(&tmp1, &tmp1, &bary1->v);
    fp2_sqr(&tmp2, &bary1->w);
    fp2_sub(&uuvv_m_ww1, &tmp1, &tmp2);
    fp2_add(&uuvv_p_ww1, &tmp1, &tmp2);

    fp2_sqr(&tmp1, &bary2->u);
    fp2_sub(&tmp1, &tmp1, &bary2->v);
    fp2_sqr(&tmp2, &bary2->w);
    fp2_sub(&uuvv_m_ww2, &tmp1, &tmp2);
    fp2_add(&uuvv_p_ww2, &tmp1, &tmp2);

    fp2_mul(&uw1, &bary1->u, &pmw1);
    fp2_add(&uw1, &uw1, &uw1);

    fp2_mul(&uw2, &bary2->u, &bary2->w);
    fp2_add(&uw2, &uw2, &uw2);

    fp2_set_zero(&out->t);

    switch (coeff->pos)
    {
    case 1:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->y, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->z, &tmp3, &uuvv_m_ww1);

        break;
    case 2:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->x, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->z, &tmp3, &uuvv_m_ww1);

        break;
    case 3:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->z, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->y, &tmp3, &uuvv_m_ww1);

        break;
    case 5:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->x, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->y, &tmp3, &uuvv_m_ww1);

        break;
    case 6:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->z, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->x, &tmp3, &uuvv_m_ww1);

        break;
    case 7:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww1, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->y, &tmp1, &uuvv_m_ww2);

        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->x, &tmp3, &uuvv_m_ww1);

        break;
    case 0:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uuvv_m_ww2, &uuvv_m_ww1);
        fp2_mul(&out->z, &tmp1, &coeff->g00[2]);

        break;
    case 4:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[2]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uuvv_m_ww2, &uuvv_m_ww1);
        fp2_mul(&out->y, &tmp1, &coeff->g00[1]);

        break;
    case 8:
        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[2]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uuvv_p_ww2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uuvv_p_ww1);
        fp2_mul(&tmp3, &uuvv_p_ww2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uuvv_m_ww2, &uuvv_m_ww1);
        fp2_mul(&out->x, &tmp1, &coeff->g00[0]);

        break;
    default:
        debug_print("position should be between 0 and 8");
        break;
    }
}

// Assume P = K1_8 or K2_8 and Q = (0, 0), i.e. v1 = v2 = 0
void superglue_evaluation_Ki8(
    theta_point_t *out,
    const superglue_coefficient_t *coeff,
    const theta_couple_point_t *Ki_8)
{
    fp2_t FP2_ZERO;
    fp2_set_zero(&FP2_ZERO);
    add_components_t bary1 = {Ki_8->P1.x, FP2_ZERO, Ki_8->P1.z};
    add_components_t bary2 = {Ki_8->P2.x, FP2_ZERO, Ki_8->P2.z};
    // if (fp2_is_zero(&bary1.w)) {
    //     fp2_set_one(&bary1.u);
    // }
    // if (fp2_is_zero(&bary2.w)) {
    //     fp2_set_one(&bary2.u);
    // }
    fp2_set_zero(&out->t);

    fp2_t uwp1, uwp2, uwm1, uwm2;
    fp2_t uw1, uw2;
    fp2_t pmw1, tmp1, tmp2, tmp3, tmp4;

    fp2_copy(&pmw1, &bary1.w);
    if (coeff->sign_flip) {
        fp2_neg(&pmw1, &pmw1);
    }

    fp2_sqr(&tmp1, &bary1.u);
    fp2_sqr(&tmp2, &bary1.w);
    fp2_add(&uwp1, &tmp1, &tmp2);
    fp2_sub(&uwm1, &tmp1, &tmp2);

    fp2_sqr(&tmp3, &bary2.u);
    fp2_sqr(&tmp4, &bary2.w);
    fp2_add(&uwp2, &tmp3, &tmp4);
    fp2_sub(&uwm2, &tmp3, &tmp4);

    fp2_mul(&uw1, &bary1.u, &pmw1);
    fp2_add(&uw1, &uw1, &uw1);

    fp2_mul(&uw2, &bary2.u, &bary2.w);
    fp2_add(&uw2, &uw2, &uw2);

    switch (coeff->pos)
    {
    case 1:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp1, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->y, &tmp1, &uwm2);

        fp2_mul(&tmp3, &uwp2, &coeff->g00[2]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->z, &tmp3, &uwm1);

        break;
    case 2:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp1, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->x, &tmp1, &uwm2);

        fp2_mul(&tmp3, &uwp2, &coeff->g00[2]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->z, &tmp3, &uwm1);

        break;
    case 3:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp1, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->z, &tmp1, &uwm2);

        fp2_mul(&tmp3, &uwp2, &coeff->g00[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->y, &tmp3, &uwm1);

        break;
    case 5:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp1, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->x, &tmp1, &uwm2);

        fp2_mul(&tmp3, &uwp2, &coeff->g00[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->y, &tmp3, &uwm1);

        break;
    case 6:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp1, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->z, &tmp1, &uwm2);

        fp2_mul(&tmp3, &uwp2, &coeff->g00[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->x, &tmp3, &uwm1);

        break;
    case 7:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g02[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g03[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp1, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw1, &coeff->g02[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&out->y, &tmp1, &uwm2);

        fp2_mul(&tmp3, &uwp2, &coeff->g00[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g01[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&out->x, &tmp3, &uwm1);

        break;
    case 0:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uwm2, &uwm1);
        fp2_mul(&out->z, &tmp1, &coeff->g00[2]);

        break;
    case 4:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[0]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[2]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->x, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[0]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uwm2, &uwm1);
        fp2_mul(&out->y, &tmp1, &coeff->g00[1]);

        break;
    case 8:
        fp2_mul(&tmp1, &uwp2, &coeff->g00[1]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[1]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g01[0]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[2]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->y, &tmp1, &tmp3);

        fp2_mul(&tmp1, &uwp2, &coeff->g00[2]);
        fp2_mul(&tmp2, &uw2, &coeff->g01[0]);
        fp2_add(&tmp1, &tmp1, &tmp2);
        fp2_mul(&tmp1, &tmp1, &uwp1);
        fp2_mul(&tmp3, &uwp2, &coeff->g01[1]);
        fp2_mul(&tmp4, &uw2, &coeff->g00[1]);
        fp2_add(&tmp3, &tmp3, &tmp4);
        fp2_mul(&tmp3, &tmp3, &uw1);
        fp2_add(&out->z, &tmp1, &tmp2);

        fp2_mul(&tmp1, &uwm2, &uwm1);
        fp2_mul(&out->x, &tmp1, &coeff->g00[0]);

        break;
    default:
        debug_print("position should be between 0 and 8");
        break;
    }
}

/**
 * @brief Compute the gluing isogeny from an elliptic product
 *
 * @param out Output: the theta_gluing
 * @param K1_8 a couple point
 * @param E12 an elliptic curve product
 * @param K2_8 a point in E2[8]
 *
 * out : E1xE2 -> A of kernel [4](K1_8,K2_8)
 * if the kernel supplied has the incorrect order, or gluing seems malformed,
 * returns 0, otherwise returns 1.
 */
int superglue_gluing_compute(
    theta_gluing_t *out,
    const theta_couple_curve_t *E12,
    const theta_couple_jac_point_t *xyK1_8,
    const theta_couple_jac_point_t *xyK2_8,
    bool verify)
{
    // Ensure that we have been given the eight torsion
#ifndef NDEBUG
    {
        int check = test_jac_order_twof(&xyK1_8->P1, &E12->E1, 3);
        if (!check)
            debug_print("xyK1_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK2_8->P1, &E12->E1, 3);
        if (!check)
            debug_print("xyK2_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK1_8->P2, &E12->E2, 3);
        if (!check)
            debug_print("xyK2_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK2_8->P2, &E12->E2, 3);
        if (!check)
            debug_print("xyK2_8->P2 does not have order 8");
    }
#endif

    out->domain = *E12;

    // Given points in E[8] x E[8] we need the four torsion below
    theta_couple_jac_point_t xyK1_4, xyK2_4;

    double_couple_jac_point_iter(&xyK1_4, 1, E12, xyK1_8);
    double_couple_jac_point_iter(&xyK2_4, 1, E12, xyK2_8);

    // Convert from (X:Y:Z) coordinates to (X:Z)
    theta_couple_point_t K1_8, K2_8;
    theta_couple_point_t K1_4, K2_4;

    couple_jac_to_xz(&K1_8, xyK1_8);
    couple_jac_to_xz(&K2_8, xyK2_8);
    couple_jac_to_xz(&K1_4, &xyK1_4);
    couple_jac_to_xz(&K2_4, &xyK2_4);

    copy_point(&out->K1_8.P1, &K1_8.P1);
    copy_point(&out->K1_8.P2, &K1_8.P2);

    copy_jac_point(&out->jac_K1_8.P1, &xyK1_8->P1);
    copy_jac_point(&out->jac_K1_8.P2, &xyK1_8->P2);

    fp2_mul(&out->proj_K1_8.P1.x, &xyK1_8->P1.x, &xyK1_8->P1.z);
    fp2_copy(&out->proj_K1_8.P1.y, &xyK1_8->P1.y);
    fp2_mul(&out->proj_K1_8.P1.z, &K1_8.P1.z, &xyK1_8->P1.z);

    fp2_mul(&out->proj_K1_8.P2.x, &xyK1_8->P2.x, &xyK1_8->P2.z);
    fp2_copy(&out->proj_K1_8.P2.y, &xyK1_8->P2.y);
    fp2_mul(&out->proj_K1_8.P2.z, &K1_8.P2.z, &xyK1_8->P2.z);

    // Check whether we have been given a valid K[8] for computation
    //    theta_couple_jac_point_t xyK1_2, xyK2_2;
    //    theta_couple_point_t K1_2, K2_2;
    //    double_couple_jac_point_iter(&xyK1_2, 1, E12, &xyK1_4);
    //    double_couple_jac_point_iter(&xyK2_2, 1, E12, &xyK2_4);
    //    couple_jac_to_xz(&K1_2, &xyK1_2);
    //    couple_jac_to_xz(&K2_2, &xyK2_2);
    //    if (!verify_two_torsion(&K1_2, &K2_2, E12)) {
    //        debug_print("gluing failed as kernel does not have correct order");
    //        return 0;
    //    }

    // Get the coefficients for superglue_evaluation
    get_superglue_coefficients(out, &K1_4, &K2_4);

    // apply the base change to the kernel
    theta_point_t TT1, TT2;

    // add_components_t bary1 = {K1_8.P1.x, {ZERO, ZERO}, K1_8.P1.z};
    // add_components_t bary2 = {K1_8.P2.x, {ZERO, ZERO}, K1_8.P2.z};
    // superglue_evaluation(&TT1, &out->coeff, &bary1, &bary2);
    // add_components_t bary3 = {K2_8.P1.x, {ZERO, ZERO}, K2_8.P1.z};
    // add_components_t bary4 = {K2_8.P2.x, {ZERO, ZERO}, K2_8.P2.z};
    // superglue_evaluation(&TT2, &out->coeff, &bary3, &bary4);
    superglue_evaluation_Ki8(&TT1, &out->coeff, &K1_8); // TT1 = [xA : xB : yC : yD]
    superglue_evaluation_Ki8(&TT2, &out->coeff, &K2_8); // TT2 = [zA : tB : zC : tD]

    // If the kernel is well formed then TT1.t and TT2.t are zero
    // if they are not, we exit early as the signature we are validating
    // is probably malformed
    if (!(fp2_is_zero(&TT1.t) & fp2_is_zero(&TT2.t))) {
        debug_print("gluing failed TT1.t or TT2.t is not zero");
        return 0;
    }
    out->zero_idx = 3;

    // Test our projective factors are non zero
    if (fp2_is_zero(&TT1.x) | fp2_is_zero(&TT2.x) | fp2_is_zero(&TT1.y) | fp2_is_zero(&TT2.z) | fp2_is_zero(&TT1.z))
        return 0; // invalid input

    // Projective factor: Ax
    fp2_mul(&out->codomain.x, &TT1.x, &TT2.x);
    fp2_mul(&out->codomain.y, &TT1.y, &TT2.x);
    fp2_mul(&out->codomain.z, &TT1.x, &TT2.z);
    fp2_set_zero(&out->codomain.t);
    // Projective factor: ABCxz
    fp2_mul(&out->precomputation.x, &TT1.y, &TT2.z);
    fp2_copy(&out->precomputation.y, &out->codomain.z);
    fp2_copy(&out->precomputation.z, &out->codomain.y);
    fp2_set_zero(&out->precomputation.t);

    // Compute the two components of phi(K1_8) = (x:x:y:y).
    fp2_mul(&out->imageK1_8.x, &TT1.x, &out->precomputation.x);
    fp2_mul(&out->imageK1_8.y, &TT1.z, &out->precomputation.z);

    // If K1_8 and K2_8 are our 8-torsion points, this ensures that the
    // 4-torsion points [2]K1_8 and [2]K2_8 are isotropic.
    if (verify) {
        fp2_t t1, t2;
        fp2_mul(&t1, &TT1.y, &out->precomputation.y);
        if (!fp2_is_equal(&out->imageK1_8.x, &t1))
            return 0;
        fp2_mul(&t1, &TT2.x, &out->precomputation.x);
        fp2_mul(&t2, &TT2.z, &out->precomputation.z);
        if (!fp2_is_equal(&t2, &t1))
            return 0;
    }

    // compute the final codomain
    hadamard(&out->codomain, &out->codomain);
    return 1;
}

int superglue_gluing_eval_point_jac(theta_point_t *image,
    const theta_couple_jac_point_t *P,
    const theta_gluing_t *phi,
    bool special_case,
    const fp2_t *twist_P)
{
    add_components_t add_comp1, add_comp2;
    fp2_t FP2_ONE;
    fp2_set_one(&FP2_ONE);

    // Compute the cross addition components of P1+Q1 and P2+Q2
    if (fp2_is_zero(&P->P1.x) && fp2_is_zero(&P->P1.z)) {
        fp2_copy(&add_comp1.u, &phi->K1_8.P1.x);
        fp2_set_zero(&add_comp1.v);
        fp2_copy(&add_comp1.w, &phi->K1_8.P1.z);
    } else {
        jac_to_xz_add_components(&add_comp1, &P->P1, &phi->jac_K1_8.P1, &phi->domain.E1, twist_P, &FP2_ONE);
    }

    if (fp2_is_zero(&P->P2.x) && fp2_is_zero(&P->P2.z)) {
        fp2_copy(&add_comp2.u, &phi->K1_8.P2.x);
        fp2_set_zero(&add_comp2.v);
        fp2_copy(&add_comp2.w, &phi->K1_8.P2.z);
    } else {
        jac_to_xz_add_components(&add_comp2, &P->P2, &phi->jac_K1_8.P2, &phi->domain.E2, twist_P, &FP2_ONE);
    }

    superglue_evaluation(image, &phi->coeff, &add_comp1, &add_comp2, twist_P, &FP2_ONE);

    if (special_case && !fp2_is_zero(&image->t))
        return 0;

    // Compute (x, y, z, t)
    // As imageK1_8 = (x:x:y:y), its inverse is (y:y:x:x).
    fp2_mul(&image->x, &image->x, &phi->imageK1_8.y);
    fp2_mul(&image->y, &image->y, &phi->imageK1_8.y);
    fp2_mul(&image->z, &image->z, &phi->imageK1_8.x);
    fp2_mul(&image->t, &image->t, &phi->imageK1_8.x);

    hadamard(image, image);
    return 1;
}

int superglue_gluing_eval_point_proj(theta_point_t *image,
    const theta_couple_proj_point_t *P,
    const theta_gluing_t *phi,
    bool special_case,
    const fp2_t *twist_P)
{
    add_components_t add_comp1, add_comp2;
    fp2_t FP2_ONE;
    fp2_set_one(&FP2_ONE);

    // Compute the cross addition components of P1+Q1 and P2+Q2
    if (fp2_is_zero(&P->P1.x) && fp2_is_zero(&P->P1.z)) {
        fp2_copy(&add_comp1.u, &phi->K1_8.P1.x);
        fp2_set_zero(&add_comp1.v);
        fp2_copy(&add_comp1.w, &phi->K1_8.P1.z);
    } else {
        proj_to_xz_add_components(&add_comp1, &P->P1, &phi->proj_K1_8.P1, &phi->domain.E1, twist_P, &FP2_ONE);
    }

    if (fp2_is_zero(&P->P2.x) && fp2_is_zero(&P->P2.z)) {
        fp2_copy(&add_comp2.u, &phi->K1_8.P2.x);
        fp2_set_zero(&add_comp2.v);
        fp2_copy(&add_comp2.w, &phi->K1_8.P2.z);
    } else {
        proj_to_xz_add_components(&add_comp2, &P->P2, &phi->proj_K1_8.P2, &phi->domain.E2, twist_P, &FP2_ONE);
    }

    superglue_evaluation(image, &phi->coeff, &add_comp1, &add_comp2, twist_P, &FP2_ONE);

    if (special_case && !fp2_is_zero(&image->t))
        return 0;

    // Compute (x, y, z, t)
    // As imageK1_8 = (x:x:y:y), its inverse is (y:y:x:x).
    fp2_mul(&image->x, &image->x, &phi->imageK1_8.y);
    fp2_mul(&image->y, &image->y, &phi->imageK1_8.y);
    fp2_mul(&image->z, &image->z, &phi->imageK1_8.x);
    fp2_mul(&image->t, &image->t, &phi->imageK1_8.x);

    hadamard(image, image);
    return 1;
}

int superglue_gluing_eval_point_ec(
    theta_point_t *image,
    const theta_couple_point_t *P,
    const theta_gluing_t *phi)
{
    add_components_t add_comp1, add_comp2;

    // Compute the cross addition components of P1+Q1 and P2+Q2
    if (fp2_is_zero(&P->P1.x) && fp2_is_zero(&P->P1.z)) {
        fp2_copy(&add_comp1.u, &phi->K1_8.P1.x);
        fp2_set_zero(&add_comp1.v);
        fp2_copy(&add_comp1.w, &phi->K1_8.P1.z);
    } else {
        ec_to_xz_add_components(&add_comp1, &P->P1, &phi->K1_8.P1, &phi->domain.E1);
    }

    if (fp2_is_zero(&P->P2.x) && fp2_is_zero(&P->P2.z)) {
        fp2_copy(&add_comp2.u, &phi->K1_8.P2.x);
        fp2_set_zero(&add_comp2.v);
        fp2_copy(&add_comp2.w, &phi->K1_8.P2.z);
    } else {
        ec_to_xz_add_components(&add_comp2, &P->P2, &phi->K1_8.P2, &phi->domain.E2);
    }

    superglue_evaluation_half_special(image, &phi->coeff, &add_comp1, &add_comp2);

    if (!fp2_is_zero(&image->t))
        return 0;

    // Compute (x, y, z, t)
    // As imageK1_8 = (x:x:y:y), its inverse is (y:y:x:x).
    fp2_mul(&image->x, &image->x, &phi->imageK1_8.y);
    fp2_mul(&image->y, &image->y, &phi->imageK1_8.y);
    fp2_mul(&image->z, &image->z, &phi->imageK1_8.x);
    fp2_mul(&image->t, &image->t, &phi->imageK1_8.x);

    hadamard(image, image);
    return 1;
}

/**
 * @brief Evaluate a theta isogeny
 *
 * @param out Output: the evaluating point
 * @param phi a theta isogeny
 * @param P a point in the codomain of phi
 *
 * out = \tilde{phi}(P)
 *
 */
void dual_theta_isogeny_eval(theta_point_t *out, theta_isogeny_t *phi, const theta_point_t *P)
{
    if (phi->bool2) {
        hadamard(out, P);
        to_squared_theta(out, out);
    } else {
        to_squared_theta(out, P);
    }

    if (!phi->domain.precomputation) {
        theta_point_t domain_null_point;
        if (phi->bool1) {
            hadamard(&domain_null_point, &phi->domain.null_point);
        } else {
            fp2_copy(&domain_null_point.x, &phi->domain.null_point.x);
            fp2_copy(&domain_null_point.y, &phi->domain.null_point.y);
            fp2_copy(&domain_null_point.z, &phi->domain.null_point.z);
            fp2_copy(&domain_null_point.t, &phi->domain.null_point.t);
        }
        assert(!fp2_is_zero(&domain_null_point.x));
        assert(!fp2_is_zero(&domain_null_point.y));
        assert(!fp2_is_zero(&domain_null_point.z));
        if (fp2_is_zero(&domain_null_point.t)) {
            fp2_set_zero(&phi->domain.xyz0);
            fp2_mul(&phi->domain.xyt0, &domain_null_point.x, &domain_null_point.y);
            fp2_mul(&phi->domain.yzt0, &domain_null_point.y, &domain_null_point.z);
            fp2_mul(&phi->domain.xzt0, &domain_null_point.z, &domain_null_point.x);
        } else {
            fp2_t xy, zt;
            fp2_mul(&xy, &domain_null_point.x, &domain_null_point.y);
            fp2_mul(&zt, &domain_null_point.z, &domain_null_point.t);
            fp2_mul(&phi->domain.xyz0, &xy, &domain_null_point.z);
            fp2_mul(&phi->domain.xyt0, &xy, &domain_null_point.t);
            fp2_mul(&phi->domain.yzt0, &zt, &domain_null_point.y);
            fp2_mul(&phi->domain.xzt0, &zt, &domain_null_point.x);
        }
    }

    fp2_mul(&out->x, &out->x, &phi->domain.yzt0);
    fp2_mul(&out->y, &out->y, &phi->domain.xzt0);
    fp2_mul(&out->z, &out->z, &phi->domain.xyt0);
    fp2_mul(&out->t, &out->t, &phi->domain.xyz0);

    if (phi->bool1) {
        hadamard(out, out);
    }
}

void dual_apply_isomorphism(theta_point_t *res, theta_splitting_t *phi, const theta_point_t *P)
{
    // assert(phi->split_type != 0);

    theta_point_t temp;

    fp2_t x1;
    fp2_t Mk0, Mk1, Mk2, Mk3;

    // k = 0
    fp2_conj(&Mk0, &phi->M00);
    fp2_conj(&Mk1, &phi->M10);
    fp2_conj(&Mk2, &phi->M20);
    fp2_conj(&Mk3, &phi->M30);
    fp2_mul(&temp.x, &P->x, &Mk0);
    fp2_mul(&x1, &P->y, &Mk1);
    fp2_add(&temp.x, &temp.x, &x1);
    fp2_mul(&x1, &P->z, &Mk2);
    fp2_add(&temp.x, &temp.x, &x1);
    fp2_mul(&x1, &P->t, &Mk3);
    fp2_add(&temp.x, &temp.x, &x1);

    // k = 1
    fp2_conj(&Mk0, &phi->M01);
    fp2_conj(&Mk1, &phi->M11);
    fp2_conj(&Mk2, &phi->M21);
    fp2_conj(&Mk3, &phi->M31);
    fp2_mul(&temp.y, &P->x, &Mk0);
    fp2_mul(&x1, &P->y, &Mk1);
    fp2_add(&temp.y, &temp.y, &x1);
    fp2_mul(&x1, &P->z, &Mk2);
    fp2_add(&temp.y, &temp.y, &x1);
    fp2_mul(&x1, &P->t, &Mk3);
    fp2_add(&temp.y, &temp.y, &x1);

    // k = 2
    fp2_conj(&Mk0, &phi->M02);
    fp2_conj(&Mk1, &phi->M12);
    fp2_conj(&Mk2, &phi->M22);
    fp2_conj(&Mk3, &phi->M32);
    fp2_mul(&temp.z, &P->x, &Mk0);
    fp2_mul(&x1, &P->y, &Mk1);
    fp2_add(&temp.z, &temp.z, &x1);
    fp2_mul(&x1, &P->z, &Mk2);
    fp2_add(&temp.z, &temp.z, &x1);
    fp2_mul(&x1, &P->t, &Mk3);
    fp2_add(&temp.z, &temp.z, &x1);

    // k = 3
    fp2_conj(&Mk0, &phi->M03);
    fp2_conj(&Mk1, &phi->M13);
    fp2_conj(&Mk2, &phi->M23);
    fp2_conj(&Mk3, &phi->M33);
    fp2_mul(&temp.t, &P->x, &Mk0);
    fp2_mul(&x1, &P->y, &Mk1);
    fp2_add(&temp.t, &temp.t, &x1);
    fp2_mul(&x1, &P->z, &Mk2);
    fp2_add(&temp.t, &temp.t, &x1);
    fp2_mul(&x1, &P->t, &Mk3);
    fp2_add(&temp.t, &temp.t, &x1);

    fp2_copy(&res->x, &temp.x);
    fp2_copy(&res->y, &temp.y);
    fp2_copy(&res->z, &temp.z);
    fp2_copy(&res->t, &temp.t);
}

void montgomery_point_to_theta_point(
    theta_point_t *P,
    const theta_couple_point_t *P12,
    const theta_structure_t *A)
{
    fp2_t xpz, xmz;
    fp2_t a, b;
    fp2_t r1, s1, r2, s2;

    // theta coordinates of P1
    fp2_copy(&a, &A->null_point.x);
    fp2_copy(&b, &A->null_point.y);
    fp2_sub(&xmz, &P12->P1.x, &P12->P1.z);
    fp2_add(&xpz, &P12->P1.x, &P12->P1.z);
    fp2_mul(&r1, &a, &xmz);
    fp2_mul(&s1, &b, &xpz);

    // theta coordinates of P2
    fp2_copy(&a, &A->null_point.y);
    fp2_copy(&b, &A->null_point.t);
    fp2_sub(&xmz, &P12->P2.x, &P12->P2.z);
    fp2_add(&xpz, &P12->P2.x, &P12->P2.z);
    fp2_mul(&r2, &a, &xmz);
    fp2_mul(&s2, &b, &xpz);

    // P = (r1r2, s1r2, r1s2, s1s2)
    fp2_mul(&P->x, &r1, &r2);
    fp2_mul(&P->y, &s1, &r2);
    fp2_mul(&P->z, &r1, &s2);
    fp2_mul(&P->t, &s1, &s2);
}

void apply_isomorphism(theta_point_t *res, const theta_splitting_t *out, const theta_point_t *P);
void theta_point_to_montgomery_point(theta_couple_point_t *P12, const theta_point_t *P, const theta_structure_t *A);
void dual_theta_chain_eval_special_case(
    theta_couple_point_t *out,
    theta_chain_t *phi,
    const theta_couple_point_t *P12)
{
    theta_point_t tmp;

    montgomery_point_to_theta_point(&tmp, P12, &phi->last_step.B);
    dual_apply_isomorphism(&tmp, &phi->last_step, &tmp);
#ifndef NDEBUG
    // Check the splitting isogeny and its dual
    theta_point_t tmp_Q;
    apply_isomorphism(&tmp_Q, &phi->last_step, &tmp);
    theta_couple_point_t tmp_P12;
    theta_point_to_montgomery_point(&tmp_P12, &tmp_Q, &phi->last_step.B);
    assert(ec_is_equal(&P12->P1, &tmp_P12.P1));
    assert(ec_is_equal(&P12->P2, &tmp_P12.P2));
    dual_apply_isomorphism(&tmp_Q, &phi->last_step, &tmp_Q);
    assert(!fp2_is_zero(&tmp_Q.x));
    assert(!fp2_is_zero(&tmp_Q.y));
    assert(!fp2_is_zero(&tmp_Q.z));
    assert(!fp2_is_zero(&tmp_Q.t));
    fp2_t a, b;
    fp2_mul(&a, &tmp.x, &tmp_Q.y);
    fp2_mul(&b, &tmp.y, &tmp_Q.x);
    assert(fp2_is_equal(&a, &b));
    fp2_mul(&a, &tmp.x, &tmp_Q.z);
    fp2_mul(&b, &tmp.z, &tmp_Q.x);
    assert(fp2_is_equal(&a, &b));
    fp2_mul(&a, &tmp.x, &tmp_Q.t);
    fp2_mul(&b, &tmp.t, &tmp_Q.x);
    assert(fp2_is_equal(&a, &b));
#endif

    int n = phi->length;
    for (int i = 0; i < n - 1; i++) {
#ifndef NDEBUG // [ ]
        theta_point_t tmp_backup_2;
        {
            fp2_copy(&tmp_backup_2.x, &tmp.x);
            fp2_copy(&tmp_backup_2.y, &tmp.y);
            fp2_copy(&tmp_backup_2.z, &tmp.z);
            fp2_copy(&tmp_backup_2.t, &tmp.t);
            if (!phi->steps[n - 2 - i].bool2) {
                theta_point_t null_point, sq_dual_point;
                hadamard(&null_point, &phi->steps[n - 2 - i].codomain.null_point);
                to_squared_theta(&sq_dual_point, &null_point);
                fp2_t xyz, xyt, xzt, yzt;
                fp2_t XYZ, XYT, XZT, YZT;
                fp2_mul(&XYZ, &sq_dual_point.x, &sq_dual_point.y);
                fp2_mul(&XYT, &XYZ, &sq_dual_point.t);
                fp2_mul(&XYZ, &XYZ, &sq_dual_point.z);
                fp2_mul(&XZT, &sq_dual_point.z, &sq_dual_point.t);
                fp2_mul(&YZT, &XZT, &sq_dual_point.y);
                fp2_mul(&XZT, &XZT, &sq_dual_point.x);
                assert(!fp2_is_zero(&null_point.x));
                assert(!fp2_is_zero(&null_point.y));
                assert(!fp2_is_zero(&null_point.z));
                if (i == 1) {
                    assert(fp2_is_zero(&null_point.t));
                    fp2_mul(&xyt, &null_point.x, &null_point.y);
                    fp2_mul(&xzt, &null_point.x, &null_point.z);
                    fp2_mul(&yzt, &null_point.y, &null_point.z);
                    fp2_set_zero(&xyz);
                } else {
                    assert(!fp2_is_zero(&null_point.t));
                    fp2_mul(&xyz, &null_point.x, &null_point.y);
                    fp2_mul(&xyt, &xyz, &null_point.t);
                    fp2_mul(&xyz, &xyz, &null_point.z);
                    fp2_mul(&xzt, &null_point.z, &null_point.t);
                    fp2_mul(&yzt, &xzt, &null_point.y);
                    fp2_mul(&xzt, &xzt, &null_point.x);
                }

                hadamard(&tmp_backup_2, &tmp_backup_2);
                to_squared_theta(&tmp_backup_2, &tmp_backup_2);
                fp2_sqr(&tmp_backup_2.x, &tmp_backup_2.x);
                fp2_sqr(&tmp_backup_2.y, &tmp_backup_2.y);
                fp2_sqr(&tmp_backup_2.z, &tmp_backup_2.z);
                fp2_sqr(&tmp_backup_2.t, &tmp_backup_2.t);
                fp2_mul(&tmp_backup_2.x, &tmp_backup_2.x, &YZT);
                fp2_mul(&tmp_backup_2.y, &tmp_backup_2.y, &XZT);
                fp2_mul(&tmp_backup_2.z, &tmp_backup_2.z, &XYT);
                fp2_mul(&tmp_backup_2.t, &tmp_backup_2.t, &XYZ);
                hadamard(&tmp_backup_2, &tmp_backup_2);
                fp2_mul(&tmp_backup_2.x, &tmp_backup_2.x, &yzt);
                fp2_mul(&tmp_backup_2.y, &tmp_backup_2.y, &xzt);
                fp2_mul(&tmp_backup_2.z, &tmp_backup_2.z, &xyt);
                fp2_mul(&tmp_backup_2.t, &tmp_backup_2.t, &xyz);
            } else {
                double_point(&tmp_backup_2, &phi->steps[n - 2 - i].codomain, &tmp_backup_2);
            }
        }
#endif

        // Dual isogeny evaluation
        dual_theta_isogeny_eval(&tmp, &phi->steps[n - 2 - i], &tmp);

#ifndef NDEBUG // [ ]
        {
            theta_point_t tmp_2P, tmp_ffP;
            if (phi->steps[n - 2 - i].bool1) {
                theta_point_t domain_null_point, sq_codomain_dual_null_point;
                hadamard(&domain_null_point, &phi->steps[n - 2 - i].domain.null_point);
                assert(!fp2_is_zero(&domain_null_point.x));
                assert(!fp2_is_zero(&domain_null_point.y));
                assert(!fp2_is_zero(&domain_null_point.z));
                assert(fp2_is_zero(&domain_null_point.t));
                to_squared_theta(&sq_codomain_dual_null_point, &domain_null_point);
                fp2_t xy, yz, zx, XYZ, XYT, XZT, YZT;
                fp2_mul(&XYZ, &sq_codomain_dual_null_point.x, &sq_codomain_dual_null_point.y);
                fp2_mul(&XYT, &XYZ, &sq_codomain_dual_null_point.t);
                fp2_mul(&XYZ, &XYZ, &sq_codomain_dual_null_point.z);
                fp2_mul(&XZT, &sq_codomain_dual_null_point.z, &sq_codomain_dual_null_point.t);
                fp2_mul(&YZT, &XZT, &sq_codomain_dual_null_point.y);
                fp2_mul(&XZT, &XZT, &sq_codomain_dual_null_point.x);
                fp2_mul(&xy, &domain_null_point.x, &domain_null_point.y);
                fp2_mul(&yz, &domain_null_point.y, &domain_null_point.z);
                fp2_mul(&zx, &domain_null_point.z, &domain_null_point.x);
                hadamard(&tmp_2P, &tmp);
                to_squared_theta(&tmp_2P, &tmp_2P);
                fp2_sqr(&tmp_2P.x, &tmp_2P.x);
                fp2_sqr(&tmp_2P.y, &tmp_2P.y);
                fp2_sqr(&tmp_2P.z, &tmp_2P.z);
                fp2_sqr(&tmp_2P.t, &tmp_2P.t);
                fp2_mul(&tmp_2P.x, &tmp_2P.x, &YZT);
                fp2_mul(&tmp_2P.y, &tmp_2P.y, &XZT);
                fp2_mul(&tmp_2P.z, &tmp_2P.z, &XYT);
                fp2_mul(&tmp_2P.t, &tmp_2P.t, &XYZ);
                hadamard(&tmp_2P, &tmp_2P);
                assert(!fp2_is_zero(&tmp_2P.x));
                assert(!fp2_is_zero(&tmp_2P.y));
                assert(!fp2_is_zero(&tmp_2P.z));
                assert(fp2_is_zero(&tmp_2P.t));
                fp2_mul(&tmp_2P.x, &tmp_2P.x, &yz);
                fp2_mul(&tmp_2P.y, &tmp_2P.y, &zx);
                fp2_mul(&tmp_2P.z, &tmp_2P.z, &xy);
            } else {
                assert(!fp2_is_zero(&phi->steps[n - 2 - i].domain.null_point.x));
                assert(!fp2_is_zero(&phi->steps[n - 2 - i].domain.null_point.y));
                assert(!fp2_is_zero(&phi->steps[n - 2 - i].domain.null_point.z));
                assert(!fp2_is_zero(&phi->steps[n - 2 - i].domain.null_point.t));
                double_point(&tmp_2P, &phi->steps[n - 2 - i].domain, &tmp);
            }
            theta_isogeny_eval(&tmp_ffP, &phi->steps[n - 2 - i], &tmp);
            {
                theta_point_t tmp_ffP_hadamard;
                if (!phi->steps[n - 2 - i].bool2) {
                    hadamard(&tmp_ffP_hadamard, &tmp_ffP);
                } else {
                    fp2_copy(&tmp_ffP_hadamard.x, &tmp_ffP.x);
                    fp2_copy(&tmp_ffP_hadamard.y, &tmp_ffP.y);
                    fp2_copy(&tmp_ffP_hadamard.z, &tmp_ffP.z);
                    fp2_copy(&tmp_ffP_hadamard.t, &tmp_ffP.t);
                }
                fp2_t a, b;
                fp2_mul(&a, &tmp_backup_2.x, &tmp_ffP_hadamard.y);
                fp2_mul(&b, &tmp_backup_2.y, &tmp_ffP_hadamard.x);
                assert(!fp2_is_zero(&a));
                assert(fp2_is_equal(&a, &b));
                fp2_mul(&a, &tmp_backup_2.x, &tmp_ffP_hadamard.z);
                fp2_mul(&b, &tmp_backup_2.z, &tmp_ffP_hadamard.x);
                assert(!fp2_is_zero(&a));
                assert(fp2_is_equal(&a, &b));
                fp2_mul(&a, &tmp_backup_2.x, &tmp_ffP_hadamard.t);
                fp2_mul(&b, &tmp_backup_2.t, &tmp_ffP_hadamard.x);
                // assert(!fp2_is_zero(&tmp_ffP_hadamard.t));
                // assert(!fp2_is_zero(&tmp_backup_2.t));
                // assert(!fp2_is_zero(&a));
                assert(fp2_is_equal(&a, &b));
            }
            dual_theta_isogeny_eval(&tmp_ffP, &phi->steps[n - 2 - i], &tmp_ffP);
            if (phi->steps[n - 2 - i].bool1) {
                hadamard(&tmp_ffP, &tmp_ffP);
            }
            assert(!fp2_is_zero(&tmp_ffP.x));
            assert(!fp2_is_zero(&tmp_ffP.y));
            assert(!fp2_is_zero(&tmp_ffP.z));
            if (i > 0) {
                assert(!fp2_is_zero(&tmp_ffP.t));
            } else {
                assert(fp2_is_zero(&tmp_ffP.t));
            }
            fp2_t a, b;
            fp2_mul(&a, &tmp_2P.x, &tmp_ffP.y);
            fp2_mul(&b, &tmp_2P.y, &tmp_ffP.x);
            assert(!fp2_is_zero(&a));
            assert(fp2_is_equal(&a, &b));
            fp2_mul(&a, &tmp_2P.x, &tmp_ffP.z);
            fp2_mul(&b, &tmp_2P.z, &tmp_ffP.x);
            assert(!fp2_is_zero(&a));
            assert(fp2_is_equal(&a, &b));
            fp2_mul(&a, &tmp_2P.x, &tmp_ffP.t);
            fp2_mul(&b, &tmp_2P.t, &tmp_ffP.x);
            // assert(!fp2_is_zero(&a));
            assert(fp2_is_equal(&a, &b));
        }
#endif
    }

    hadamard(&tmp, &tmp);
    to_squared_theta(&tmp, &tmp);

    fp2_t x0, y0, z0, t0;
    fp2_mul(&x0, &phi->first_step.M20, &phi->first_step.M30);
    fp2_mul(&y0, &x0, &phi->first_step.M00);
    fp2_mul(&x0, &x0, &phi->first_step.M10);
    fp2_mul(&z0, &phi->first_step.M00, &phi->first_step.M10);
    fp2_mul(&t0, &z0, &phi->first_step.M20);
    fp2_mul(&z0, &z0, &phi->first_step.M30);

    fp2_mul(&tmp.x, &tmp.x, &x0);
    fp2_mul(&tmp.y, &tmp.y, &y0);
    fp2_mul(&tmp.z, &tmp.z, &z0);
    fp2_mul(&tmp.t, &tmp.t, &t0);

    // Compute the adjugate matrix of the change of basis matrix for phi->first_step
    fp2_t adjM00, adjM01, adjM02, adjM03;
    fp2_t adjM10, adjM11, adjM12, adjM13;
    fp2_t adjM20, adjM21, adjM22, adjM23;
    fp2_t adjM30, adjM31, adjM32, adjM33;
    fp2_t det_01_01, det_01_02, det_01_03, det_01_12, det_01_13, det_01_23;
    fp2_t det_23_01, det_23_02, det_23_03, det_23_12, det_23_13, det_23_23;
    
    fp2_mul(&x0, &phi->first_step.M00, &phi->first_step.M11);
    fp2_mul(&y0, &phi->first_step.M10, &phi->first_step.M01);
    fp2_sub(&det_01_01, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M00, &phi->first_step.M12);
    fp2_mul(&y0, &phi->first_step.M10, &phi->first_step.M02);
    fp2_sub(&det_01_02, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M00, &phi->first_step.M13);
    fp2_mul(&y0, &phi->first_step.M10, &phi->first_step.M03);
    fp2_sub(&det_01_03, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M01, &phi->first_step.M12);
    fp2_mul(&y0, &phi->first_step.M11, &phi->first_step.M02);
    fp2_sub(&det_01_12, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M01, &phi->first_step.M13);
    fp2_mul(&y0, &phi->first_step.M11, &phi->first_step.M03);
    fp2_sub(&det_01_13, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M02, &phi->first_step.M13);
    fp2_mul(&y0, &phi->first_step.M12, &phi->first_step.M03);
    fp2_sub(&det_01_23, &x0, &y0);
    
    fp2_mul(&x0, &phi->first_step.M20, &phi->first_step.M31);
    fp2_mul(&y0, &phi->first_step.M30, &phi->first_step.M21);
    fp2_sub(&det_23_01, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M20, &phi->first_step.M32);
    fp2_mul(&y0, &phi->first_step.M30, &phi->first_step.M22);
    fp2_sub(&det_23_02, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M20, &phi->first_step.M33);
    fp2_mul(&y0, &phi->first_step.M30, &phi->first_step.M23);
    fp2_sub(&det_23_03, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M21, &phi->first_step.M32);
    fp2_mul(&y0, &phi->first_step.M31, &phi->first_step.M22);
    fp2_sub(&det_23_12, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M21, &phi->first_step.M33);
    fp2_mul(&y0, &phi->first_step.M31, &phi->first_step.M23);
    fp2_sub(&det_23_13, &x0, &y0);
    fp2_mul(&x0, &phi->first_step.M22, &phi->first_step.M33);
    fp2_mul(&y0, &phi->first_step.M32, &phi->first_step.M23);
    fp2_sub(&det_23_23, &x0, &y0);

    fp2_mul(&x0, &phi->first_step.M11, &det_23_23);
    fp2_mul(&y0, &phi->first_step.M12, &det_23_13);
    fp2_mul(&z0, &phi->first_step.M13, &det_23_12);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM00, &t0, &y0);
    fp2_mul(&x0, &phi->first_step.M10, &det_23_23);
    fp2_mul(&y0, &phi->first_step.M12, &det_23_03);
    fp2_mul(&z0, &phi->first_step.M13, &det_23_02);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM10, &y0, &t0);
    fp2_mul(&x0, &phi->first_step.M10, &det_23_13);
    fp2_mul(&y0, &phi->first_step.M11, &det_23_03);
    fp2_mul(&z0, &phi->first_step.M13, &det_23_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM20, &t0, &y0);
    fp2_mul(&x0, &phi->first_step.M10, &det_23_12);
    fp2_mul(&y0, &phi->first_step.M11, &det_23_02);
    fp2_mul(&z0, &phi->first_step.M12, &det_23_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM30, &y0, &t0);

    fp2_mul(&x0, &phi->first_step.M01, &det_23_23);
    fp2_mul(&y0, &phi->first_step.M02, &det_23_13);
    fp2_mul(&z0, &phi->first_step.M03, &det_23_12);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM01, &y0, &t0);
    fp2_mul(&x0, &phi->first_step.M00, &det_23_23);
    fp2_mul(&y0, &phi->first_step.M02, &det_23_03);
    fp2_mul(&z0, &phi->first_step.M03, &det_23_02);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM11, &t0, &y0);
    fp2_mul(&x0, &phi->first_step.M00, &det_23_13);
    fp2_mul(&y0, &phi->first_step.M01, &det_23_03);
    fp2_mul(&z0, &phi->first_step.M03, &det_23_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM21, &y0, &t0);
    fp2_mul(&x0, &phi->first_step.M00, &det_23_12);
    fp2_mul(&y0, &phi->first_step.M01, &det_23_02);
    fp2_mul(&z0, &phi->first_step.M02, &det_23_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM31, &t0, &y0);

    fp2_mul(&x0, &phi->first_step.M31, &det_01_23);
    fp2_mul(&y0, &phi->first_step.M32, &det_01_13);
    fp2_mul(&z0, &phi->first_step.M33, &det_01_12);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM02, &t0, &y0);
    fp2_mul(&x0, &phi->first_step.M30, &det_01_23);
    fp2_mul(&y0, &phi->first_step.M32, &det_01_03);
    fp2_mul(&z0, &phi->first_step.M33, &det_01_02);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM12, &y0, &t0);
    fp2_mul(&x0, &phi->first_step.M30, &det_01_13);
    fp2_mul(&y0, &phi->first_step.M31, &det_01_03);
    fp2_mul(&z0, &phi->first_step.M33, &det_01_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM22, &t0, &y0);
    fp2_mul(&x0, &phi->first_step.M30, &det_01_12);
    fp2_mul(&y0, &phi->first_step.M31, &det_01_02);
    fp2_mul(&z0, &phi->first_step.M32, &det_01_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM32, &y0, &t0);

    fp2_mul(&x0, &phi->first_step.M21, &det_01_23);
    fp2_mul(&y0, &phi->first_step.M22, &det_01_13);
    fp2_mul(&z0, &phi->first_step.M23, &det_01_12);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM03, &y0, &t0);
    fp2_mul(&x0, &phi->first_step.M20, &det_01_23);
    fp2_mul(&y0, &phi->first_step.M22, &det_01_03);
    fp2_mul(&z0, &phi->first_step.M23, &det_01_02);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM13, &t0, &y0);
    fp2_mul(&x0, &phi->first_step.M20, &det_01_13);
    fp2_mul(&y0, &phi->first_step.M21, &det_01_03);
    fp2_mul(&z0, &phi->first_step.M23, &det_01_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM23, &y0, &t0);
    fp2_mul(&x0, &phi->first_step.M20, &det_01_12);
    fp2_mul(&y0, &phi->first_step.M21, &det_01_02);
    fp2_mul(&z0, &phi->first_step.M22, &det_01_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&adjM33, &t0, &y0);

    // Apply the adjugate matrxi for image point
    theta_point_t out_theta;
    
    fp2_mul(&out_theta.x, &adjM00, &tmp.x);
    fp2_mul(&x0, &adjM01, &tmp.y);
    fp2_add(&out_theta.x, &out_theta.x, &x0);
    fp2_mul(&x0, &adjM02, &tmp.z);
    fp2_add(&out_theta.x, &out_theta.x, &x0);
    fp2_mul(&x0, &adjM03, &tmp.t);
    fp2_add(&out_theta.x, &out_theta.x, &x0);

    fp2_mul(&out_theta.y, &adjM10, &tmp.x);
    fp2_mul(&x0, &adjM11, &tmp.y);
    fp2_add(&out_theta.y, &out_theta.y, &x0);
    fp2_mul(&x0, &adjM12, &tmp.z);
    fp2_add(&out_theta.y, &out_theta.y, &x0);
    fp2_mul(&x0, &adjM13, &tmp.t);
    fp2_add(&out_theta.y, &out_theta.y, &x0);

    fp2_mul(&out_theta.z, &adjM20, &tmp.x);
    fp2_mul(&x0, &adjM21, &tmp.y);
    fp2_add(&out_theta.z, &out_theta.z, &x0);
    fp2_mul(&x0, &adjM22, &tmp.z);
    fp2_add(&out_theta.z, &out_theta.z, &x0);
    fp2_mul(&x0, &adjM23, &tmp.t);
    fp2_add(&out_theta.z, &out_theta.z, &x0);

    fp2_mul(&out_theta.t, &adjM30, &tmp.x);
    fp2_mul(&x0, &adjM31, &tmp.y);
    fp2_add(&out_theta.t, &out_theta.t, &x0);
    fp2_mul(&x0, &adjM32, &tmp.z);
    fp2_add(&out_theta.t, &out_theta.t, &x0);
    fp2_mul(&x0, &adjM33, &tmp.t);
    fp2_add(&out_theta.t, &out_theta.t, &x0);

#ifndef NDEBUG
    // Check : out_theta == (xx, xz, zx, zz)
    fp2_mul(&x0, &out_theta.x, &out_theta.t);
    fp2_mul(&y0, &out_theta.y, &out_theta.z);
    assert(fp2_is_equal(&x0, &y0));
#endif

    // P1.X =   null_point.y * P.x + null_point.x * P.y
    // P1.Z = - null_point.y * P.x + null_point.x * P.y
    // P2.X =   null_point.t * P.y + null_point.y * P.t
    // P2.Z = - null_point.t * P.y + null_point.y * P.t
    // fp2_mul(&x0, &phi->first_step.M10, &out_theta.x);
    // fp2_mul(&y0, &phi->first_step.M00, &out_theta.y);
    // fp2_add(&out->P1.x, &y0, &x0);
    // fp2_sub(&out->P1.z, &y0, &x0);
    // fp2_mul(&x0, &phi->first_step.M30, &out_theta.y);
    // fp2_mul(&y0, &phi->first_step.M10, &out_theta.t);
    // fp2_add(&out->P1.x, &y0, &x0);
    // fp2_sub(&out->P1.z, &y0, &x0);
    fp2_copy(&out->P1.x, &out_theta.x);
    fp2_copy(&out->P1.z, &out_theta.y);
    fp2_copy(&out->P2.x, &out_theta.x);
    fp2_copy(&out->P2.z, &out_theta.z);
}

/**
 * @brief Compute the gluing isogeny from an elliptic product
 *
 * @param out Output: the theta_gluing
 * @param K1_8 a couple point
 * @param E12 an elliptic curve product
 * @param K2_8 a point in E2[8]
 *
 * out : E1xE2 -> A of kernel [4](K1_8,K2_8)
 *
 */
void
gluing_comput(theta_gluing_t *out,
              theta_couple_curve_t *E12,
              const theta_couple_jac_point_t *xyK1_8,
              const theta_couple_jac_point_t *xyK2_8
              // const theta_couple_point_t *K1_8,const theta_couple_point_t *K2_8, const
              // theta_couple_point_t *K1m2_8
)
{

    // var init
    fp2_t M1100;
    fp2_t M1101;
    fp2_t M1110;
    fp2_t M1111;
    fp2_t M1200;
    fp2_t M1201;
    fp2_t M1210;
    fp2_t M1211;
    fp2_t M2100;
    fp2_t M2101;
    fp2_t M2110;
    fp2_t M2111;
    fp2_t M2200;
    fp2_t M2201;
    fp2_t M2210;
    fp2_t M2211;
    fp2_t t001, t101, t002, t102, temp;

    theta_point_t TT1, TT2;
    theta_couple_point_t K1m2_4;
    ec_basis_t B;
    theta_couple_point_t K1_8, K2_8;

    double_couple_jac_point_iter(&out->xyK1_4, 1, E12, xyK1_8);
    double_couple_jac_point_iter(&out->xyK2_4, 1, E12, xyK2_8);

    couple_jac_to_xz(&K1_8, xyK1_8);
    couple_jac_to_xz(&K2_8, xyK2_8);
    couple_jac_to_xz(&out->K1_4, &out->xyK1_4);
    couple_jac_to_xz(&out->K2_4, &out->xyK2_4);
    // assert(test_point_order_twof(&K1_8.P1,&E12->E1,8));

    // computing the base change matrix
    ec_point_t P11, P12, P21, P22;
    fp2_t t[8];
    fp2_t temp_fp;

    ec_dbl(&P11, &E12->E1, &out->K1_4.P1);
    ec_dbl(&P12, &E12->E2, &out->K1_4.P2);
    ec_dbl(&P21, &E12->E1, &out->K2_4.P1);
    ec_dbl(&P22, &E12->E2, &out->K2_4.P2);

    fp2_copy(&t[0], &out->K1_4.P1.z);
    fp2_copy(&t[1], &out->K1_4.P2.z);
    fp2_copy(&t[2], &out->K2_4.P1.z);
    fp2_copy(&t[3], &out->K2_4.P2.z);

    // Ki_4.Pj = (xij,zij)
    // Pij = (uij,wij)
    // detij = xij wij - uij zij
    fp2_mul(&t[4], &out->K1_4.P1.x, &P11.z);
    fp2_mul(&temp_fp, &out->K1_4.P1.z, &P11.x);
    fp2_sub(&t[4], &t[4], &temp_fp);
    fp2_mul(&t[5], &out->K1_4.P2.x, &P12.z);
    fp2_mul(&temp_fp, &out->K1_4.P2.z, &P12.x);
    fp2_sub(&t[5], &t[5], &temp_fp);
    fp2_mul(&t[6], &out->K2_4.P1.x, &P21.z);
    fp2_mul(&temp_fp, &out->K2_4.P1.z, &P21.x);
    fp2_sub(&t[6], &t[6], &temp_fp);
    fp2_mul(&t[7], &out->K2_4.P2.x, &P22.z);
    fp2_mul(&temp_fp, &out->K2_4.P2.z, &P22.x);
    fp2_sub(&t[7], &t[7], &temp_fp);

    fp2_batched_inv(t, 8);

    // Mij10 = uij xij /detij - xij/zij
    fp2_mul(&temp_fp, &out->K1_4.P1.x, &t[0]);
    fp2_mul(&M1110, &out->K1_4.P1.x, &P11.x);
    fp2_mul(&M1110, &M1110, &t[4]);
    fp2_sub(&M1110, &M1110, &temp_fp);
    fp2_mul(&temp_fp, &out->K1_4.P2.x, &t[1]);
    fp2_mul(&M1210, &out->K1_4.P2.x, &P12.x);
    fp2_mul(&M1210, &M1210, &t[5]);
    fp2_sub(&M1210, &M1210, &temp_fp);
    fp2_mul(&temp_fp, &out->K2_4.P1.x, &t[2]);
    fp2_mul(&M2110, &out->K2_4.P1.x, &P21.x);
    fp2_mul(&M2110, &M2110, &t[6]);
    fp2_sub(&M2110, &M2110, &temp_fp);
    fp2_mul(&temp_fp, &out->K2_4.P2.x, &t[3]);
    fp2_mul(&M2210, &out->K2_4.P2.x, &P22.x);
    fp2_mul(&M2210, &M2210, &t[7]);
    fp2_sub(&M2210, &M2210, &temp_fp);

    // Mij11 = uij zij * detij
    fp2_mul(&M1111, &P11.x, &t[4]);
    fp2_mul(&M1111, &M1111, &out->K1_4.P1.z);
    fp2_mul(&M1211, &P12.x, &t[5]);
    fp2_mul(&M1211, &M1211, &out->K1_4.P2.z);
    fp2_mul(&M2111, &P21.x, &t[6]);
    fp2_mul(&M2111, &M2111, &out->K2_4.P1.z);
    fp2_mul(&M2211, &P22.x, &t[7]);
    fp2_mul(&M2211, &M2211, &out->K2_4.P2.z);

    // Mij00 = -Mij11
    fp2_neg(&M1100, &M1111);
    fp2_neg(&M1200, &M1211);
    fp2_neg(&M2100, &M2111);
    fp2_neg(&M2200, &M2211);

    // Mij01 = - wij zij detij
    fp2_mul(&M1101, &P11.z, &t[4]);
    fp2_mul(&M1101, &M1101, &out->K1_4.P1.z);
    fp2_neg(&M1101, &M1101);
    fp2_mul(&M1201, &P12.z, &t[5]);
    fp2_mul(&M1201, &M1201, &out->K1_4.P2.z);
    fp2_neg(&M1201, &M1201);
    fp2_mul(&M2101, &P21.z, &t[6]);
    fp2_mul(&M2101, &M2101, &out->K2_4.P1.z);
    fp2_neg(&M2101, &M2101);
    fp2_mul(&M2201, &P22.z, &t[7]);
    fp2_mul(&M2201, &M2201, &out->K2_4.P2.z);
    fp2_neg(&M2201, &M2201);

    // multiplication of the matrices
    // t001,t101 (resp t002,t102) first column of M11 * M21 (resp M12 * M22)
    fp2_mul(&t001, &M1100, &M2100);
    fp2_mul(&temp, &M1101, &M2110);
    fp2_add(&t001, &t001, &temp);

    fp2_mul(&t101, &M1110, &M2100);
    fp2_mul(&temp, &M1111, &M2110);
    fp2_add(&t101, &t101, &temp);

    fp2_mul(&t002, &M1200, &M2200);
    fp2_mul(&temp, &M1201, &M2210);
    fp2_add(&t002, &t002, &temp);

    fp2_mul(&t102, &M1210, &M2200);
    fp2_mul(&temp, &M1211, &M2210);
    fp2_add(&t102, &t102, &temp);

    // trace for the first row
    fp2_sett(&out->M00, 1);
    fp2_mul(&temp, &t001, &t002);
    fp2_add(&out->M00, &out->M00, &temp);
    fp2_mul(&temp, &M2100, &M2200);
    fp2_add(&out->M00, &out->M00, &temp);
    fp2_mul(&temp, &M1100, &M1200);
    fp2_add(&out->M00, &out->M00, &temp);

    fp2_mul(&out->M01, &t001, &t102);
    fp2_mul(&temp, &M2100, &M2210);
    fp2_add(&out->M01, &out->M01, &temp);
    fp2_mul(&temp, &M1100, &M1210);
    fp2_add(&out->M01, &out->M01, &temp);

    fp2_mul(&out->M02, &t101, &t002);
    fp2_mul(&temp, &M2110, &M2200);
    fp2_add(&out->M02, &out->M02, &temp);
    fp2_mul(&temp, &M1110, &M1200);
    fp2_add(&out->M02, &out->M02, &temp);

    fp2_mul(&out->M03, &t101, &t102);
    fp2_mul(&temp, &M2110, &M2210);
    fp2_add(&out->M03, &out->M03, &temp);
    fp2_mul(&temp, &M1110, &M1210);
    fp2_add(&out->M03, &out->M03, &temp);

    // Compute the action of (0,out.K2_4.P2) for the second row
    fp2_mul(&temp, &M2201, &out->M01);
    fp2_mul(&out->M10, &M2200, &out->M00);
    fp2_add(&out->M10, &out->M10, &temp);

    fp2_mul(&temp, &M2211, &out->M01);
    fp2_mul(&out->M11, &M2210, &out->M00);
    fp2_add(&out->M11, &out->M11, &temp);

    fp2_mul(&temp, &M2201, &out->M03);
    fp2_mul(&out->M12, &M2200, &out->M02);
    fp2_add(&out->M12, &out->M12, &temp);

    fp2_mul(&temp, &M2211, &out->M03);
    fp2_mul(&out->M13, &M2210, &out->M02);
    fp2_add(&out->M13, &out->M13, &temp);

    // compute the action of (K1_4.P1,0) for the third row
    fp2_mul(&temp, &M1101, &out->M02);
    fp2_mul(&out->M20, &M1100, &out->M00);
    fp2_add(&out->M20, &out->M20, &temp);

    fp2_mul(&temp, &M1101, &out->M03);
    fp2_mul(&out->M21, &M1100, &out->M01);
    fp2_add(&out->M21, &out->M21, &temp);

    fp2_mul(&temp, &M1111, &out->M02);
    fp2_mul(&out->M22, &M1110, &out->M00);
    fp2_add(&out->M22, &out->M22, &temp);

    fp2_mul(&temp, &M1111, &out->M03);
    fp2_mul(&out->M23, &M1110, &out->M01);
    fp2_add(&out->M23, &out->M23, &temp);

    // compute the action of (K1_4.P1,K2_4.P2) for the final row
    fp2_mul(&temp, &M1101, &out->M12);
    fp2_mul(&out->M30, &M1100, &out->M10);
    fp2_add(&out->M30, &out->M30, &temp);

    fp2_mul(&temp, &M1101, &out->M13);
    fp2_mul(&out->M31, &M1100, &out->M11);
    fp2_add(&out->M31, &out->M31, &temp);

    fp2_mul(&temp, &M1111, &out->M12);
    fp2_mul(&out->M32, &M1110, &out->M10);
    fp2_add(&out->M32, &out->M32, &temp);

    fp2_mul(&temp, &M1111, &out->M13);
    fp2_mul(&out->M33, &M1110, &out->M11);
    fp2_add(&out->M33, &out->M33, &temp);

    // apply the base change
    base_change(&out->T1_8, out, &K1_8);
    base_change(&out->T2_8, out, &K2_8);

    // computing the codomain
    // computation of the zero index
    to_squared_theta(&TT1, &out->T1_8);
    to_squared_theta(&TT2, &out->T2_8);

    if (fp2_is_zero(&TT1.x)) {
        out->zero_idx = 0;
    } else if (fp2_is_zero(&TT1.y)) {
        out->zero_idx = 1;
    } else if (fp2_is_zero(&TT1.z)) {
        out->zero_idx = 2;
    } else {
        out->zero_idx = 3;
    }

#ifndef NDEBUG
    fp2_t a1, a2;
    choose_index_theta_point(&a1, out->zero_idx, &TT1);
    choose_index_theta_point(&a2, out->zero_idx, &TT2);
    assert(fp2_is_zero(&a1) && fp2_is_zero(&a2));
#endif

    choose_index_theta_point(&t001, 1 ^ out->zero_idx, &TT2);
    choose_index_theta_point(&t002, 2 ^ out->zero_idx, &TT1);
    choose_index_theta_point(&t101, 3 ^ out->zero_idx, &TT2);
    choose_index_theta_point(&t102, 3 ^ out->zero_idx, &TT1);

    fp2_t t1, t2, t3, t4;
    // t1 = t001;t2=t002;t3=t101;t4=t102;
    fp2_copy(&t1, &t001);
    fp2_copy(&t2, &t002);
    fp2_copy(&t3, &t101);
    fp2_copy(&t4, &t102);

    fp2_t t_inv[4];
    fp2_copy(&t_inv[0], &t1);
    fp2_copy(&t_inv[1], &t2);
    fp2_copy(&t_inv[2], &t3);
    fp2_copy(&t_inv[3], &t4);
    fp2_batched_inv(t_inv, 4);
    fp2_copy(&t001, &t_inv[0]);
    fp2_copy(&t002, &t_inv[1]);
    fp2_copy(&t101, &t_inv[2]);
    fp2_copy(&t102, &t_inv[3]);

    // Compute A,B,C,D
    fp2_sett(&temp, 0);
    set_index_theta_point(&out->codomain, 0 ^ out->zero_idx, &temp);
    fp2_mul(&temp, &t101, &t1);
    set_index_theta_point(&out->codomain, 1 ^ out->zero_idx, &temp);
    fp2_mul(&temp, &t2, &t102);
    set_index_theta_point(&out->codomain, 2 ^ out->zero_idx, &temp);
    fp2_sett(&temp, 1);
    set_index_theta_point(&out->codomain, 3 ^ out->zero_idx, &temp);
    choose_index_theta_point(&temp, 2 ^ out->zero_idx, &out->codomain);

    // compute precomp
    fp2_sett(&temp, 0);
    set_index_theta_point(&out->precomputation, out->zero_idx, &temp);
    fp2_mul(&temp, &t001, &t3);
    set_index_theta_point(&out->precomputation, 1 ^ out->zero_idx, &temp);
    fp2_mul(&temp, &t4, &t002);
    set_index_theta_point(&out->precomputation, 2 ^ out->zero_idx, &temp);
    fp2_sett(&temp, 1);
    set_index_theta_point(&out->precomputation, 3 ^ out->zero_idx, &temp);
    choose_index_theta_point(&temp, 2 ^ out->zero_idx, &out->precomputation);

    // compute the final codomain
    hadamard(&out->codomain, &out->codomain);
}

// sub routine of the gluing eval
void
gluing_eval_point(theta_point_t *image1,
                  const theta_couple_point_t *P,
                  const theta_couple_point_t *Pt,
                  const theta_gluing_t *phi)
{

    theta_point_t T, Tt;
    fp2_t x, y, z, t, temp;

    // apply the basis change
    base_change(&T, phi, P);
    base_change(&Tt, phi, Pt);

    // apply the to_squared_theta transform
    to_squared_theta(&T, &T);
    to_squared_theta(&Tt, &Tt);

    // compute y,z,t
    choose_index_theta_point(&y, 1 ^ phi->zero_idx, &T);
    choose_index_theta_point(&temp, 1 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&y, &y, &temp);
    choose_index_theta_point(&z, 2 ^ phi->zero_idx, &T);
    choose_index_theta_point(&temp, 2 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&z, &z, &temp);
    choose_index_theta_point(&t, 3 ^ phi->zero_idx, &T);

    //  normalize
    if (!fp2_is_zero(&z)) {
        fp2_copy(&x, &z);
        choose_index_theta_point(&temp, 3 ^ phi->zero_idx, &Tt);

        fp2_mul(&y, &y, &temp);
        fp2_mul(&z, &z, &temp);
        fp2_mul(&t, &t, &temp);
    } else {
        choose_index_theta_point(&temp, 2 ^ phi->zero_idx, &Tt);
        choose_index_theta_point(&x, 2 ^ phi->zero_idx, &phi->precomputation);
        fp2_mul(&temp, &temp, &x);
        fp2_copy(&x, &t);

        fp2_mul(&y, &y, &temp);
        fp2_mul(&z, &z, &temp);
        fp2_mul(&t, &t, &temp);
    }

    // recover x
    choose_index_theta_point(&temp, 1 ^ phi->zero_idx, &Tt);
    fp2_mul(&x, &x, &temp);
    choose_index_theta_point(&temp, 1 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&x, &x, &temp);

    // fill the image coordinates
    set_index_theta_point(image1, 0 ^ phi->zero_idx, &x);
    set_index_theta_point(image1, 1 ^ phi->zero_idx, &y);
    set_index_theta_point(image1, 2 ^ phi->zero_idx, &z);
    set_index_theta_point(image1, 3 ^ phi->zero_idx, &t);

    // hadamard
    hadamard(image1, image1);
}

// same as gluing_eval_point but in the very special case where we already know that the point will
// have a zero coordinate at the place where the zero coordinate of the dual_theta_nullpoint would
// have made the computation difficult
void
gluing_eval_point_special_case(theta_point_t *image,
                               const theta_couple_point_t *P,
                               const theta_gluing_t *phi)
{
    theta_point_t T;
    fp2_t x, y, z, t, temp;

    // apply the basis change
    base_change(&T, phi, P);

    // apply the to_squared_theta transform
    to_squared_theta(&T, &T);

    // compute y,z,t
    choose_index_theta_point(&y, 1 ^ phi->zero_idx, &T);
    choose_index_theta_point(&temp, 1 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&y, &y, &temp);
    choose_index_theta_point(&z, 2 ^ phi->zero_idx, &T);
    choose_index_theta_point(&temp, 2 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&z, &z, &temp);
    choose_index_theta_point(&t, 3 ^ phi->zero_idx, &T);

    // this is always 0 in the special case
    fp2_set_zero(&x);

    // fill the image coordinates
    set_index_theta_point(image, 0 ^ phi->zero_idx, &x);
    set_index_theta_point(image, 1 ^ phi->zero_idx, &y);
    set_index_theta_point(image, 2 ^ phi->zero_idx, &z);
    set_index_theta_point(image, 3 ^ phi->zero_idx, &t);

    // hadamard
    hadamard(image, image);
}

// sub routine of the gluing eval
void
gluing_eval_point_no_help(theta_point_t *image1,
                          const theta_couple_jac_point_t *P,
                          const theta_couple_curve_t *E12,
                          const theta_gluing_t *phi)
{

    theta_couple_point_t Pt;
    theta_couple_jac_point_t tmp;
    theta_point_t T, Tt;
    fp2_t x, y, z, t, temp;

#ifndef NDEBUG
    jac_to_xz(&Pt.P2, &phi->xyK1_4.P2);
    assert(test_point_order_twof(&phi->K1_4.P2, &E12->E2, 2));
    assert(test_point_order_twof(&Pt.P2, &E12->E2, 2));
#endif

    // we compute Pt
    ADD(&tmp.P1, &P->P1, &phi->xyK1_4.P1, &E12->E1);
    ADD(&tmp.P2, &P->P2, &phi->xyK1_4.P2, &E12->E2);

    couple_jac_to_xz(&Pt, &tmp);

    // apply the basis change
    // on the translated point
    base_change(&Tt, phi, &Pt);
    // on the point to evaluate
    couple_jac_to_xz(&Pt, P);
    base_change(&T, phi, &Pt);

    // apply the to_squared_theta transform
    to_squared_theta(&T, &T);
    to_squared_theta(&Tt, &Tt);

    // compute y,z,t
    choose_index_theta_point(&y, 1 ^ phi->zero_idx, &T);
    choose_index_theta_point(&temp, 1 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&y, &y, &temp);
    choose_index_theta_point(&z, 2 ^ phi->zero_idx, &T);
    choose_index_theta_point(&temp, 2 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&z, &z, &temp);
    choose_index_theta_point(&t, 3 ^ phi->zero_idx, &T);

    //  normalize
    if (!fp2_is_zero(&z)) {
        fp2_copy(&x, &z);
        choose_index_theta_point(&temp, 3 ^ phi->zero_idx, &Tt);

        fp2_mul(&y, &y, &temp);
        fp2_mul(&z, &z, &temp);
        fp2_mul(&t, &t, &temp);
    } else {
        choose_index_theta_point(&temp, 2 ^ phi->zero_idx, &Tt);
        choose_index_theta_point(&x, 2 ^ phi->zero_idx, &phi->precomputation);
        fp2_mul(&temp, &temp, &x);
        fp2_copy(&x, &t);

        fp2_mul(&y, &y, &temp);
        fp2_mul(&z, &z, &temp);
        fp2_mul(&t, &t, &temp);
    }

    // recover x
    choose_index_theta_point(&temp, 1 ^ phi->zero_idx, &Tt);
    fp2_mul(&x, &x, &temp);
    choose_index_theta_point(&temp, 1 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&x, &x, &temp);

    // fill the image coordinates
    set_index_theta_point(image1, 0 ^ phi->zero_idx, &x);
    set_index_theta_point(image1, 1 ^ phi->zero_idx, &y);
    set_index_theta_point(image1, 2 ^ phi->zero_idx, &z);
    set_index_theta_point(image1, 3 ^ phi->zero_idx, &t);

    // hadamard
    hadamard(image1, image1);
}

/**
 * @brief Evaluate a gluing isogeny from an elliptic product on a basis
 *
 * @param image1 Output: the theta_point of the image of the first couple of points
 * @param image2 Output : the theta point of the image of the second couple of points
 * @param P : a couple point in E12
 * @param Q : a couple point in E12
 * @param PmQ : a couple point in E12 corresponding to P-Q
 * @param a : ibz
 * @param b : ibz
 * @param E12 : an elliptic product
 * @param phi : a gluing isogeny E1 x E2 -> A
 *
 * The integers a,b are such that the phi.K1_4 = a P + b Q
 * out : phi( P ), Phi (Q)
 *
 */
void
gluing_eval_basis(theta_point_t *image1,
                  theta_point_t *image2,
                  // const theta_couple_point_t *P,const theta_couple_point_t *Q, const
                  // theta_couple_point_t *PmQ,const ibz_t *a, const ibz_t *b,
                  const theta_couple_jac_point_t *xyT1,
                  const theta_couple_jac_point_t *xyT2,
                  theta_couple_curve_t *E12,
                  const theta_gluing_t *phi)
{

    theta_couple_point_t P, Pt;
    theta_couple_jac_point_t T;

    // add the point to push with xyK1_4
    ADD(&T.P1, &xyT1->P1, &phi->xyK1_4.P1, &E12->E1);
    ADD(&T.P2, &xyT1->P2, &phi->xyK1_4.P2, &E12->E2);

    couple_jac_to_xz(&Pt, &T);
    couple_jac_to_xz(&P, xyT1);

    // then we evaluate the gluing
    gluing_eval_point(image1, &P, &Pt, phi);

    // then we do the same on the second point
    ADD(&T.P1, &xyT2->P1, &phi->xyK1_4.P1, &E12->E1);
    ADD(&T.P2, &xyT2->P2, &phi->xyK1_4.P2, &E12->E2);

    couple_jac_to_xz(&Pt, &T);
    couple_jac_to_xz(&P, xyT2);
    // then we evaluate the gluing
    gluing_eval_point(image2, &P, &Pt, phi);
}

/**
 * @brief Compute  a (2,2) isogeny in dimension 2 in the theta_model
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_8 a point in A[8]
 * @param T2_8 a point in A[8]
 * @param bool1 a boolean
 * @param boo2 a boolean
 *
 * out : A -> B of kernel [4](T1_8,T2_8)
 * bool1 controls if the domain is in standard or dual coordinates
 * bool2 controls if the codomain is in standard or dual coordinates
 *
 */
void
theta_isogeny_comput(theta_isogeny_t *out,
                     const theta_structure_t *A,
                     const theta_point_t *T1_8,
                     const theta_point_t *T2_8,
                     int bool1,
                     int bool2)
{
    out->bool1 = bool1;
    out->bool2 = bool2;
    out->domain = *A;
    out->T1_8 = *T1_8;
    out->T2_8 = *T2_8;
    out->codomain.precomputation = 0;

    theta_point_t TT1, TT2;

    // fp2_t xA_inv,zA_inv,tB_inv;

    if (bool1) {
        hadamard(&TT1, T1_8);
        to_squared_theta(&TT1, &TT1);
        hadamard(&TT2, T2_8);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT1, T1_8);
        to_squared_theta(&TT2, T2_8);
    }

    if (!bool1 && A->precomputation) {
        // xA_inv = TT1.x;
        // zA_inv = TT2.x;
        // tB_inv = TT2.y;
        // fp2_t invs[3];
        // fp2_copy(&invs[0],&xA_inv);
        // fp2_copy(&invs[1],&zA_inv);
        // fp2_copy(&invs[2],&tB_inv);
        // fp2_batched_inv(invs,3);
        // fp2_copy(&xA_inv,&invs[0]);
        // fp2_copy(&zA_inv,&invs[1]);
        // fp2_copy(&tB_inv,&invs[2]);

        // computation of B,C,D for the codomain
        // fp2_sett(&out->codomain.null_point.x,1);
        // fp2_mul(&out->codomain.null_point.y,&xA_inv,&TT1.y);
        // fp2_mul(&out->codomain.null_point.z,&zA_inv,&TT2.z);
        // fp2_mul(&out->codomain.null_point.t,&tB_inv,&TT2.t);
        // fp2_mul(&out->codomain.null_point.t,&out->codomain.null_point.t,&out->codomain.null_point.y);

        fp2_mul(&out->codomain.null_point.x, &TT2.x, &TT2.y);
        fp2_mul(&out->codomain.null_point.y, &out->codomain.null_point.x, &TT1.y);
        fp2_mul(&out->codomain.null_point.x, &out->codomain.null_point.x, &TT1.x);
        fp2_mul(&out->codomain.null_point.z, &TT2.z, &TT1.x);
        fp2_mul(&out->codomain.null_point.z, &out->codomain.null_point.z, &TT2.y);
        fp2_mul(&out->codomain.null_point.t, &TT2.t, &TT2.x);
        fp2_mul(&out->codomain.null_point.t, &out->codomain.null_point.t, &TT1.y);

        // computation B_inv,C_inv,D_inv for the precomputation
        // fp2_sett(&out->precomputation.x,1);
        //     fp2_copy(&out->precomputation.x,&out->codomain.null_point.x);
        //     fp2_mul(&out->precomputation.y,&out->codomain.null_point.y,&out->domain.Y0);
        //     fp2_mul(&out->precomputation.z,&out->codomain.null_point.z,&out->domain.Z0);
        //     fp2_mul(&out->precomputation.t,&out->codomain.null_point.t,&out->domain.T0);

        fp2_mul(&out->precomputation.x, &out->codomain.null_point.x, &out->domain.YZT0);
        fp2_mul(&out->precomputation.y, &out->codomain.null_point.y, &out->domain.XZT0);
        fp2_mul(&out->precomputation.z, &out->codomain.null_point.z, &out->domain.XYT0);
        fp2_mul(&out->precomputation.t, &out->codomain.null_point.t, &out->domain.XYZ0);

        // fp2_mul(&out->precomputation.x,&out->precomputation.x,&out->domain.ZT0);
        // fp2_mul(&out->precomputation.y,&out->precomputation.y,&out->domain.ZT0);
        // fp2_mul(&out->precomputation.z,&out->precomputation.z,&out->domain.XY0);
        // fp2_mul(&out->precomputation.t,&out->precomputation.t,&out->domain.XY0);

    } else {
        // fp2_t xB_inv,zC_inv,tD_inv;
        // xA_inv = TT1.x;
        // zA_inv = TT2.x;
        // tB_inv = TT2.y;
        // xB_inv=TT1.y;
        // zC_inv=TT2.z;
        // tD_inv=TT2.t;
        // fp2_t invs[3];

        // fp2_copy(&invs[0],&xA_inv);
        // fp2_copy(&invs[1],&zA_inv);
        // fp2_copy(&invs[2],&tB_inv);
        // // fp2_copy(&invs[3],&xB_inv);
        // // fp2_copy(&invs[4],&zC_inv);
        // // fp2_copy(&invs[5],&tD_inv);
        // fp2_batched_inv(invs,3);
        // fp2_copy(&xA_inv,&invs[0]);
        // fp2_copy(&zA_inv,&invs[1]);
        // fp2_copy(&tB_inv,&invs[2]);
        // fp2_copy(&xB_inv,&invs[3]);
        // fp2_copy(&zC_inv,&invs[4]);
        // fp2_copy(&tD_inv,&invs[5]);

        fp2_mul(&out->codomain.null_point.x, &TT2.x, &TT2.y);
        fp2_mul(&out->codomain.null_point.y, &out->codomain.null_point.x, &TT1.y);
        fp2_mul(&out->codomain.null_point.x, &out->codomain.null_point.x, &TT1.x);
        fp2_mul(&out->codomain.null_point.z, &TT2.z, &TT1.x);
        fp2_mul(&out->codomain.null_point.z, &out->codomain.null_point.z, &TT2.y);
        fp2_mul(&out->codomain.null_point.t, &TT2.t, &TT2.x);
        fp2_mul(&out->codomain.null_point.t, &out->codomain.null_point.t, &TT1.y);

        // fp2_sett(&out->codomain.null_point.x,1);

        // // computation of B,C,D for the codomain
        // fp2_mul(&out->codomain.null_point.y,&xA_inv,&TT1.y);
        // fp2_mul(&out->codomain.null_point.z,&zA_inv,&TT2.z);
        // fp2_mul(&out->codomain.null_point.t,&tB_inv,&TT2.t);
        // fp2_mul(&out->codomain.null_point.t,&out->codomain.null_point.t,&out->codomain.null_point.y);

        // computation of B_inv,C_inv,D_inv for the precomputation
        // fp2_sett(&out->precomputation.x,1);
        // fp2_mul(&out->precomputation.y,&xB_inv,&TT1.x);
        // fp2_mul(&out->precomputation.z,&zC_inv,&TT2.x);
        // fp2_mul(&out->precomputation.t,&tD_inv,&TT2.y);
        // fp2_mul(&out->precomputation.t,&out->precomputation.t,&out->precomputation.y);

        fp2_mul(&out->precomputation.x, &TT2.t, &TT2.z);
        fp2_mul(&out->precomputation.y, &out->precomputation.x, &TT1.x);
        fp2_mul(&out->precomputation.x, &out->precomputation.x, &TT1.y);
        fp2_mul(&out->precomputation.z, &TT2.x, &TT1.y);
        fp2_mul(&out->precomputation.z, &out->precomputation.z, &TT2.t);
        fp2_mul(&out->precomputation.t, &TT2.z, &TT2.y);
        fp2_mul(&out->precomputation.t, &out->precomputation.t, &TT1.x);
    }
    if (bool2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
}

/**
 * @brief Compute a (2,2) isogeny when only the 4 torsion above the kernel is known and not the 8
 * torsion
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_4 a point in A[4]
 * @param T2_4 a point in A[4]
 * @param bool1 a boolean
 * @param boo2 a boolean
 *
 * out : A -> B of kernel [2](T1_4,T2_4)
 * bool1 controls if the domain is in standard or dual coordinates
 * bool2 controls if the codomain is in standard or dual coordinates
 *
 */
void
theta_isogeny_comput4(theta_isogeny_t *out,
                      const theta_structure_t *A,
                      const theta_point_t *T1_4,
                      const theta_point_t *T2_4,
                      int bool1,
                      int bool2)
{

    out->bool1 = bool1;
    out->bool2 = bool2;
    out->domain = *A;
    out->T1_8 = *T1_4;
    out->T2_8 = *T2_4;
    out->codomain.precomputation = 0;

    theta_point_t TT1, TT2;
    // we will compute:
    // TT1 = (xAB, _ , xCD, _)
    // TT2 = (AA,BB,CC,DD)

    // fp2_t xA_inv,zA_inv,tB_inv;

    if (bool1) {

        hadamard(&TT1, T1_4);
        to_squared_theta(&TT1, &TT1);

        hadamard(&TT2, &A->null_point);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT1, T1_4);
        to_squared_theta(&TT2, &A->null_point);
    }

    fp2_t sqaabb, sqaacc;
    fp2_mul(&sqaabb, &TT2.x, &TT2.y);
    fp2_mul(&sqaacc, &TT2.x, &TT2.z);
    // sqaabb = sqrt(AA*BB)
    fp2_sqrt(&sqaabb);
    // sqaacc = sqrt(AA*CC)
    fp2_sqrt(&sqaacc);

    // we compute out->codomain.null_point = (xAB * sqaacc * AA, xAB *sqaabb *sqaacc, xCD*sqaabb *
    // AA) out->precomputation = (xAB * BB * CC *DD , sqaabb * CC * DD * xAB , sqaacc * BB* DD * xAB
    // , xCD * sqaabb *sqaacc * BB)

    fp2_mul(&out->codomain.null_point.y, &sqaabb, &sqaacc);
    fp2_mul(&out->precomputation.t, &out->codomain.null_point.y, &TT1.z);
    fp2_mul(&out->codomain.null_point.y,
            &out->codomain.null_point.y,
            &TT1.x); // done for out->codomain.null_point.y

    fp2_mul(&out->codomain.null_point.t, &TT1.z, &sqaabb);
    fp2_mul(&out->codomain.null_point.t,
            &out->codomain.null_point.t,
            &TT2.x); // done for out->codomain.null_point.t

    fp2_mul(&out->codomain.null_point.x, &TT1.x, &TT2.x);
    fp2_mul(&out->codomain.null_point.z,
            &out->codomain.null_point.x,
            &TT2.z); // done for out->codomain.null_point.z
    fp2_mul(&out->codomain.null_point.x,
            &out->codomain.null_point.x,
            &sqaacc); // done for out->codomain.null_point.x

    fp2_mul(&out->precomputation.x, &TT1.x, &TT2.t);
    fp2_mul(&out->precomputation.z, &out->precomputation.x, &TT2.y);
    fp2_mul(&out->precomputation.x, &out->precomputation.x, &TT2.z);
    fp2_mul(
        &out->precomputation.y, &out->precomputation.x, &sqaabb); // done for out->precomputation.y
    fp2_mul(
        &out->precomputation.x, &out->precomputation.x, &TT2.y); // done for out->precomputation.x
    fp2_mul(
        &out->precomputation.z, &out->precomputation.z, &sqaacc); // done for out->precomputation.z
    fp2_mul(
        &out->precomputation.t, &out->precomputation.t, &TT2.y); // done for out->precomputation.t

    if (bool2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
}

/**
 * @brief Compute a (2,2) isogeny when only the kernel is known and not the 8 or 4 torsion above
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_2 a point in A[2]
 * @param T2_2 a point in A[2]
 * @param bool1 a boolean
 * @param boo2 a boolean
 *
 * out : A -> B of kernel (T1_2,T2_2)
 * bool1 controls if the domain is in standard or dual coordinates
 * bool2 controls if the codomain is in standard or dual coordinates
 *
 */
void
theta_isogeny_comput2(theta_isogeny_t *out,
                      const theta_structure_t *A,
                      const theta_point_t *T1_2,
                      const theta_point_t *T2_2,
                      int bool1,
                      int bool2)
{

    out->bool1 = bool1;
    out->bool2 = bool2;
    out->domain = *A;
    out->T1_8 = *T1_2;
    out->T2_8 = *T2_2;
    out->codomain.precomputation = 0;

    theta_point_t TT2;
    // we will compute:
    // TT2 = (AA,BB,CC,DD)

    if (bool1) {
        hadamard(&TT2, &A->null_point);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT2, &A->null_point);
    }

    // we compute out->codomain.null_point = (AA,sqaabb, sqaacc, sqaadd)
    // out->precomputation = (  BB * CC *DD , sqaabb * CC * DD , sqaacc * BB* DD , sqaadd * BB * CC)
    fp2_copy(&out->codomain.null_point.x, &TT2.x);
    fp2_mul(&out->codomain.null_point.y, &TT2.x, &TT2.y);
    fp2_mul(&out->codomain.null_point.z, &TT2.x, &TT2.z);
    fp2_mul(&out->codomain.null_point.t, &TT2.x, &TT2.t);
    fp2_sqrt(&out->codomain.null_point.y);
    fp2_sqrt(&out->codomain.null_point.z);
    fp2_sqrt(&out->codomain.null_point.t);
    // fp2_neg(&out->codomain.null_point.t,&out->codomain.null_point.t);

    fp2_mul(&out->precomputation.x, &TT2.z, &TT2.t);
    fp2_mul(&out->precomputation.y,
            &out->precomputation.x,
            &out->codomain.null_point.y); // done for out->precomputation.y
    fp2_mul(
        &out->precomputation.x, &out->precomputation.x, &TT2.y); // done for out->precomputation.x
    fp2_mul(&out->precomputation.z, &TT2.t, &out->codomain.null_point.z);
    fp2_mul(
        &out->precomputation.z, &out->precomputation.z, &TT2.y); // done for out->precomputation.z
    fp2_mul(&out->precomputation.t, &TT2.z, &out->codomain.null_point.t);
    fp2_mul(
        &out->precomputation.t, &out->precomputation.t, &TT2.y); // done for out->precomputation.t

    if (bool2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
}

/**
 * @brief Evaluate a theta isogeny
 *
 * @param out Output: the evaluating point
 * @param phi a theta isogeny
 * @param P a point in the domain of phi
 *
 * out = phi(P)
 *
 */
void
theta_isogeny_eval(theta_point_t *out, const theta_isogeny_t *phi, const theta_point_t *P)
{

    if (phi->bool1) {
        hadamard(out, P);
        to_squared_theta(out, out);
    } else {
        to_squared_theta(out, P);
    }
    fp2_mul(&out->x, &out->x, &phi->precomputation.x);
    fp2_mul(&out->y, &out->y, &phi->precomputation.y);
    fp2_mul(&out->z, &out->z, &phi->precomputation.z);
    fp2_mul(&out->t, &out->t, &phi->precomputation.t);

    if (phi->bool2) {
        hadamard(out, out);
    }
}

void
apply_isomorphism(theta_point_t *res, const theta_splitting_t *out, const theta_point_t *P)
{

    theta_point_t temp;

    fp2_t x1;
    fp2_mul(&temp.x, &P->x, &out->M00);
    fp2_mul(&x1, &P->y, &out->M01);
    fp2_add(&temp.x, &temp.x, &x1);
    fp2_mul(&x1, &P->z, &out->M02);
    fp2_add(&temp.x, &temp.x, &x1);
    fp2_mul(&x1, &P->t, &out->M03);
    fp2_add(&temp.x, &temp.x, &x1);

    fp2_mul(&temp.y, &P->x, &out->M10);
    fp2_mul(&x1, &P->y, &out->M11);
    fp2_add(&temp.y, &temp.y, &x1);
    fp2_mul(&x1, &P->z, &out->M12);
    fp2_add(&temp.y, &temp.y, &x1);
    fp2_mul(&x1, &P->t, &out->M13);
    fp2_add(&temp.y, &temp.y, &x1);

    fp2_mul(&temp.z, &P->x, &out->M20);
    fp2_mul(&x1, &P->y, &out->M21);
    fp2_add(&temp.z, &temp.z, &x1);
    fp2_mul(&x1, &P->z, &out->M22);
    fp2_add(&temp.z, &temp.z, &x1);
    fp2_mul(&x1, &P->t, &out->M23);
    fp2_add(&temp.z, &temp.z, &x1);

    fp2_mul(&temp.t, &P->x, &out->M30);
    fp2_mul(&x1, &P->y, &out->M31);
    fp2_add(&temp.t, &temp.t, &x1);
    fp2_mul(&x1, &P->z, &out->M32);
    fp2_add(&temp.t, &temp.t, &x1);
    fp2_mul(&x1, &P->t, &out->M33);
    fp2_add(&temp.t, &temp.t, &x1);

    fp2_copy(&res->x, &temp.x);
    fp2_copy(&res->y, &temp.y);
    fp2_copy(&res->z, &temp.z);
    fp2_copy(&res->t, &temp.t);
}

/**
 * @brief Compute the splitting isomorphism from a theta structure to the product theta structure,
 * returns false if the given theta structure is not isomorphic to an elliptic product
 *
 * @return a boolean indicating if A is isomorphic to an elliptic product
 * @param out: the splitting isomorphism
 * @param A : the theta_structure in consideration
 *
 * out : A -> B where B is theta product associated to ExF an elliptic product
 */
int
splitting_comput(theta_splitting_t *out, const theta_structure_t *A)
{

    // init
    int good = 0;
    int even_index[10][2] = {
        { 0, 0 },
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 1, 0 },
        { 1, 2 },
        { 2, 0 },
        { 2, 1 },
        { 3, 0 },
        { 3, 3 }
    };
    int chi_eval[4][4] = {
        { 1,  1,  1,  1 },
        { 1, -1,  1, -1 },
        { 1,  1, -1, -1 },
        { 1, -1, -1,  1 }
    };

    fp2_t U_cst, temp, temp2;

    // enumerate through all indices
    for (int i = 0; i < 10; i++) {
        fp2_sett(&U_cst, 0);
        for (int t = 0; t < 4; t++) {
            choose_index_theta_point(&temp2, t, &A->null_point);
            choose_index_theta_point(&temp, t ^ even_index[i][1], &A->null_point);
            fp2_mul(&temp, &temp, &temp2);
            fp2_sett(&temp2, chi_eval[even_index[i][0]][t]);
            fp2_mul(&temp, &temp, &temp2);
            fp2_add(&U_cst, &U_cst, &temp);
        }
        if (fp2_is_zero(&U_cst)) {
            good = 1 + i;
            break;
        }
    }

    // temp = sqrt{-1}
    // TODO precompute this?
    fp2_sett(&temp, -1);
    fp2_sqrt(&temp);

    // compute the matrix
    if (good == 1) {
        // (0, 0): [1, i, 1, i,
        //          1, -i, -1, i,
        //          1, i, -1, -i,
        //          -1, i, -1, i],
        fp2_sett(&out->M00, 1);
        fp2_copy(&out->M01, &temp);
        fp2_sett(&out->M02, 1);
        fp2_copy(&out->M03, &temp);
        fp2_sett(&out->M10, 1);
        fp2_neg(&out->M11, &temp);
        fp2_sett(&out->M12, -1);
        fp2_copy(&out->M13, &temp);
        fp2_sett(&out->M20, 1);
        fp2_copy(&out->M21, &temp);
        fp2_neg(&out->M23, &temp);
        fp2_sett(&out->M22, -1);
        fp2_sett(&out->M30, -1);
        fp2_copy(&out->M31, &temp);
        fp2_sett(&out->M32, -1);
        fp2_copy(&out->M33, &temp);

    } else if (good == 2) {
        // (0, 1): [1, 0, 0, 0,
        //          0, 0, 0, 1,
        //          0, 0, 1, 0,
        //          0, -1, 0, 0],
        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 0);
        fp2_sett(&out->M02, 0);
        fp2_sett(&out->M03, 0);
        fp2_sett(&out->M10, 0);
        fp2_sett(&out->M11, 0);
        fp2_sett(&out->M12, 0);
        fp2_sett(&out->M13, 1);
        fp2_sett(&out->M20, 0);
        fp2_sett(&out->M21, 0);
        fp2_sett(&out->M22, 1);
        fp2_sett(&out->M23, 0);
        fp2_sett(&out->M30, 0);
        fp2_sett(&out->M31, -1);
        fp2_sett(&out->M32, 0);
        fp2_sett(&out->M33, 0);
    } else if (good == 3) {
        // (0, 2): [1, 0, 0, 0,
        //          0, 1, 0, 0,
        //          0, 0, 0, 1,
        //      0, 0, -1, 0],
        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 0);
        fp2_sett(&out->M02, 0);
        fp2_sett(&out->M03, 0);
        fp2_sett(&out->M10, 0);
        fp2_sett(&out->M11, 1);
        fp2_sett(&out->M12, 0);
        fp2_sett(&out->M13, 0);
        fp2_sett(&out->M20, 0);
        fp2_sett(&out->M21, 0);
        fp2_sett(&out->M22, 0);
        fp2_sett(&out->M23, 1);
        fp2_sett(&out->M30, 0);
        fp2_sett(&out->M31, 0);
        fp2_sett(&out->M32, -1);
        fp2_sett(&out->M33, 0);
    } else if (good == 4) {
        // (0, 3): [1, 0, 0, 0,
        //          0, 1, 0, 0,
        //          0, 0, 1, 0,
        //          0, 0, 0, -1],
        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 0);
        fp2_sett(&out->M02, 0);
        fp2_sett(&out->M03, 0);
        fp2_sett(&out->M10, 0);
        fp2_sett(&out->M11, 1);
        fp2_sett(&out->M12, 0);
        fp2_sett(&out->M13, 0);
        fp2_sett(&out->M20, 0);
        fp2_sett(&out->M21, 0);
        fp2_sett(&out->M22, 1);
        fp2_sett(&out->M23, 0);
        fp2_sett(&out->M30, 0);
        fp2_sett(&out->M31, 0);
        fp2_sett(&out->M32, 0);
        fp2_sett(&out->M33, -1);

    } else if (good == 7) {
        // (2, 0): [1, 1, 1, 1,
        //          1, -1, 1, -1,
        //          1, -1, -1, 1,
        //          -1, -1, 1, 1],
        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 1);
        fp2_sett(&out->M02, 1);
        fp2_sett(&out->M03, 1);
        fp2_sett(&out->M10, 1);
        fp2_sett(&out->M11, -1);
        fp2_sett(&out->M12, 1);
        fp2_sett(&out->M13, -1);
        fp2_sett(&out->M20, 1);
        fp2_sett(&out->M21, -1);
        fp2_sett(&out->M22, -1);
        fp2_sett(&out->M23, 1);
        fp2_sett(&out->M30, -1);
        fp2_sett(&out->M31, -1);
        fp2_sett(&out->M32, 1);
        fp2_sett(&out->M33, 1);
    } else if (good == 8) {
        //(2, 1): [1, 1, 1, 1,
        //          1, -1, 1, -1,
        //          1, -1, -1, 1,
        //          1, 1, -1, -1],

        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 1);
        fp2_sett(&out->M02, 1);
        fp2_sett(&out->M03, 1);
        fp2_sett(&out->M10, 1);
        fp2_sett(&out->M11, -1);
        fp2_sett(&out->M12, 1);
        fp2_sett(&out->M13, -1);
        fp2_sett(&out->M20, 1);
        fp2_sett(&out->M21, -1);
        fp2_sett(&out->M22, -1);
        fp2_sett(&out->M23, 1);
        fp2_sett(&out->M30, 1);
        fp2_sett(&out->M31, 1);
        fp2_sett(&out->M32, -1);
        fp2_sett(&out->M33, -1);
    } else if (good == 5) {
        // (1, 0): [1, 1, 1, 1,
        //          1, -1, -1, 1,
        //          1, 1, -1, -1,
        //          -1, 1, -1, 1],
        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 1);
        fp2_sett(&out->M02, 1);
        fp2_sett(&out->M03, 1);
        fp2_sett(&out->M10, 1);
        fp2_sett(&out->M11, -1);
        fp2_sett(&out->M12, -1);
        fp2_sett(&out->M13, 1);
        fp2_sett(&out->M20, 1);
        fp2_sett(&out->M21, 1);
        fp2_sett(&out->M22, -1);
        fp2_sett(&out->M23, -1);
        fp2_sett(&out->M30, -1);
        fp2_sett(&out->M31, 1);
        fp2_sett(&out->M32, -1);
        fp2_sett(&out->M33, 1);

    } else if (good == 6) {
        //(1, 2): [1, 0, 0, 0,
        //          0, 1, 0, 0,
        //          0, 0, 0, 1,
        //          0, 0, 1, 0],
        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 0);
        fp2_sett(&out->M02, 0);
        fp2_sett(&out->M03, 0);
        fp2_sett(&out->M10, 0);
        fp2_sett(&out->M11, 1);
        fp2_sett(&out->M12, 0);
        fp2_sett(&out->M13, 0);
        fp2_sett(&out->M20, 0);
        fp2_sett(&out->M21, 0);
        fp2_sett(&out->M22, 0);
        fp2_sett(&out->M23, 1);
        fp2_sett(&out->M30, 0);
        fp2_sett(&out->M31, 0);
        fp2_sett(&out->M32, 1);
        fp2_sett(&out->M33, 0);

    } else if (good == 9) {
        // (3, 0): [1, 1, 1, 1,
        //          1, -1, 1, -1,
        //          1, 1, -1, -1,
        //          -1, 1, 1, -1],
        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 1);
        fp2_sett(&out->M02, 1);
        fp2_sett(&out->M03, 1);
        fp2_sett(&out->M10, 1);
        fp2_sett(&out->M11, -1);
        fp2_sett(&out->M12, 1);
        fp2_sett(&out->M13, -1);
        fp2_sett(&out->M20, 1);
        fp2_sett(&out->M21, 1);
        fp2_sett(&out->M22, -1);
        fp2_sett(&out->M23, -1);
        fp2_sett(&out->M30, -1);
        fp2_sett(&out->M31, 1);
        fp2_sett(&out->M32, 1);
        fp2_sett(&out->M33, -1);
    } else if (good == 10) {
        // (3, 3): [1, 0, 0, 0,
        //          0, 1, 0, 0,
        //          0, 0, 1, 0,
        //          0, 0, 0, 1],
        fp2_sett(&out->M00, 1);
        fp2_sett(&out->M01, 0);
        fp2_sett(&out->M02, 0);
        fp2_sett(&out->M03, 0);
        fp2_sett(&out->M10, 0);
        fp2_sett(&out->M11, 1);
        fp2_sett(&out->M12, 0);
        fp2_sett(&out->M13, 0);
        fp2_sett(&out->M20, 0);
        fp2_sett(&out->M21, 0);
        fp2_sett(&out->M22, 1);
        fp2_sett(&out->M23, 0);
        fp2_sett(&out->M30, 0);
        fp2_sett(&out->M31, 0);
        fp2_sett(&out->M32, 0);
        fp2_sett(&out->M33, 1);
    }

    // now we apply the isomorphism if it was computed
    if (good) {
        apply_isomorphism(&out->B.null_point, out, &A->null_point);
    }
    return good;
}

void
theta_product_structure_to_elliptic_product(theta_couple_curve_t *E12, theta_structure_t *A)
{
    fp2_t xx, yy, temp1, temp2;

    ec_curve_init(&(E12->E1));
    ec_curve_init(&(E12->E2));

    // xx = x², yy = y²
    fp2_sqr(&xx, &A->null_point.x);
    fp2_sqr(&yy, &A->null_point.y);

    // A1 = ( (xx² + yy²)² + (xx² - yy²)² )  / (xx² - yy²)
    fp2_add(&temp1, &xx, &yy);
    fp2_sub(&temp2, &xx, &yy);
    fp2_mul(&E12->E1.C, &temp1, &temp2);
    fp2_sqr(&temp1, &temp1);
    fp2_sqr(&temp2, &temp2);
    fp2_add(&temp1, &temp1, &temp2);
    // TODO fix this!
    fp2_sett(&E12->E1.A, 1);
    fp2_mul(&E12->E1.A, &E12->E1.A, &temp1);
    fp2_neg(&E12->E1.A, &E12->E1.A);
    // fp2_sett(&E12->E1.C,1);

    // same with y,t
    fp2_sqr(&xx, &A->null_point.y);
    fp2_sqr(&yy, &A->null_point.t);

    // A2 = ( (xx² + yy²)² + (xx² - yy²)² )  / (xx² - yy²)
    fp2_add(&temp1, &xx, &yy);
    fp2_sub(&temp2, &xx, &yy);
    fp2_mul(&E12->E2.C, &temp1, &temp2);
    fp2_sqr(&temp1, &temp1);
    fp2_sqr(&temp2, &temp2);
    fp2_add(&temp1, &temp1, &temp2);
    fp2_sett(&E12->E2.A, 1);
    fp2_mul(&E12->E2.A, &E12->E2.A, &temp1);
    fp2_neg(&E12->E2.A, &E12->E2.A);
}

/**
 * @brief Compute  a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the theta_chain
 * @param n : the length of the isogeny chain
 * @param E12 an elliptic curve product
 * @param T1 a couple point on E12[2^(n+2)]
 * @param T2 a couple point on E12[2^(n+2)]
 * @param T1m2 a couple point on E12[2^(n+2)] equal to T1-T2
 *
 * out : E1xE2 -> E3xE4 of kernel [4](T1,T2)
 *
 */
void
theta_chain_comput_naive(theta_chain_t *out,
                         int n,
                         theta_couple_curve_t *E12,
                         const theta_couple_point_t *T1,
                         const theta_couple_point_t *T2,
                         const theta_couple_point_t *T1m2)
{

    theta_couple_point_t P1, P2, P1m2;
    theta_point_t Q1, Q2, R1, R2;
    theta_isogeny_t steps[n - 1];
    theta_structure_t codomain;

    ibz_t a, b;
    ibz_init(&a);
    ibz_init(&b);

    // init of the isogeny chain
    out->domain = *E12;
    out->length = n;
    out->T1 = *T1;
    out->T2 = *T2;
    out->steps = malloc((n - 1) * sizeof(theta_isogeny_t));

    // First, we compute the first step
    // multiply by 2^n-1
    double_couple_point_iter(&P1, n - 1, E12, T1);
    double_couple_point_iter(&P2, n - 1, E12, T2);
    double_couple_point_iter(&P1m2, n - 1, E12, T1m2);

#ifndef NDEBUG
    // checking that the points have order 8
    ec_point_t test1, test2;
    test1 = P1.P1;
    test2 = P1.P2;
    ec_dbl(&test1, &E12->E1, &test1);
    ec_dbl(&test1, &E12->E1, &test1);
    ec_dbl(&test2, &E12->E2, &test2);
    ec_dbl(&test2, &E12->E2, &test2);
    assert(!fp2_is_zero(&test1.z));
    assert(!fp2_is_zero(&test2.z));
    ec_dbl(&test1, &E12->E1, &test1);
    ec_dbl(&test2, &E12->E2, &test2);
    assert(fp2_is_zero(&test1.z));
    assert(fp2_is_zero(&test2.z));
    test1 = P2.P1;
    test2 = P2.P2;
    ec_dbl(&test1, &E12->E1, &test1);
    ec_dbl(&test1, &E12->E1, &test1);
    ec_dbl(&test2, &E12->E2, &test2);
    ec_dbl(&test2, &E12->E2, &test2);
    assert(!fp2_is_zero(&test1.z));
    assert(!fp2_is_zero(&test2.z));
    ec_dbl(&test1, &E12->E1, &test1);
    ec_dbl(&test2, &E12->E2, &test2);
    assert(fp2_is_zero(&test1.z));
    assert(fp2_is_zero(&test2.z));
#endif

    // TODO : lift T1,T2,T1m2 to xy coordinates
    theta_couple_jac_point_t xyT1, xyT2;
    assert(0);

    // compute the gluing isogeny
    gluing_comput(&out->first_step, E12, &xyT1, &xyT2);

    // set-up the theta_structure for the first codomain
    codomain.null_point = out->first_step.codomain;
    codomain.precomputation = 0;
    theta_precomputation(&codomain);

    // push the kernel through the gluing isogeny
    // need to setup the input before
    ibz_pow(&a, &ibz_const_two, n);
    ibz_set(&b, 0);
    gluing_eval_basis(&Q1,
                      &Q2,
                      // T1,T2,T1m2,&a,&b,
                      &xyT1,
                      &xyT2,
                      E12,
                      &out->first_step);

    for (int i = 0; i < n - 1; i++) {

        // computing the kernel of the next step
        double_iter(&R1, &codomain, &Q1, n - i - 2);
        double_iter(&R2, &codomain, &Q2, n - i - 2);

        // computing the next step
        if (i == n - 3) {
            theta_isogeny_comput(&steps[i], &codomain, &R1, &R2, 0, 0);
        } else if (i == n - 2) {
            theta_isogeny_comput(&steps[i], &codomain, &R1, &R2, 1, 0);
        } else {
            theta_isogeny_comput(&steps[i], &codomain, &R1, &R2, 0, 1);
        }

        // updating the codomain
        codomain = steps[i].codomain;

        // pushing the kernel
        if (i < n - 2) {
            theta_isogeny_eval(&Q1, &steps[i], &Q1);
            theta_isogeny_eval(&Q2, &steps[i], &Q2);
        }
    }

    // copying the steps
    out->steps = steps;

    // final splitting step
    int is_split = splitting_comput(&out->last_step, &steps[n - 2].codomain);
    assert(is_split);

    // computing the curves of the codomain
    theta_product_structure_to_elliptic_product(&out->codomain, &out->last_step.B);

    ibz_finalize(&a);
    ibz_finalize(&b);
}

void
theta_chain_comput_rec(theta_chain_t *out,
                       theta_structure_t *codomain,
                       theta_point_t *R1,
                       theta_point_t *R2,
                       long len,
                       long index,
                       bool advance,
                       theta_point_t *P1,
                       theta_point_t *P2,
                       long stacklen,
                       long total_length)
{

    if (len == 0) {
        return;
    }
    if (len == 1) {

        // first we compute the isogeny and update the codomain

        if (index == total_length - 2) {
            theta_isogeny_comput(&out->steps[index], codomain, R1, R2, 0, 0);
        } else if (index == total_length - 1) {
            theta_isogeny_comput(&out->steps[index], codomain, R1, R2, 1, 0);
        } else {
            theta_isogeny_comput(&out->steps[index], codomain, R1, R2, 0, 1);
        }

        *codomain = out->steps[index].codomain;
        // push points
        for (int i = 0; i < stacklen; i++) {
            theta_isogeny_eval(P1 + i, &out->steps[index], P1 + i);
            theta_isogeny_eval(P2 + i, &out->steps[index], P2 + i);
        }

    } else {
        long right = 2 * len / 3;
        long left = len - right;
        P1[stacklen] = *R1;
        P2[stacklen] = *R2;
        double_iter(R1, codomain, R1, left);
        double_iter(R2, codomain, R2, left);

        theta_chain_comput_rec(
            out, codomain, R1, R2, right, index, advance, P1, P2, stacklen + 1, total_length);
        R1[right * advance] = P1[stacklen];
        R2[right * advance] = P2[stacklen];
        theta_chain_comput_rec(out,
                               codomain,
                               R1 + right * advance,
                               R2 + right * advance,
                               left,
                               right + index,
                               advance,
                               P1,
                               P2,
                               stacklen,
                               total_length);
    }
}

void
theta_chain_comput_balanced(theta_chain_t *out,
                            int n,
                            theta_couple_curve_t *E12,
                            const theta_couple_point_t *T1,
                            const theta_couple_point_t *T2,
                            const theta_couple_point_t *T1m2)
{

    theta_couple_point_t P1, P2, P1m2;
    theta_point_t Q1, Q2, R1, R2;
    theta_isogeny_t steps[n - 1];
    theta_structure_t codomain;

    ibz_t a, b;
    ibz_init(&a);
    ibz_init(&b);

    // init of the isogeny chain
    out->domain = *E12;
    out->length = n;
    out->T1 = *T1;
    out->T2 = *T2;
    out->steps = malloc((n - 1) * sizeof(theta_isogeny_t));

    clock_t t = tic();
    // lift the basis
    theta_couple_jac_point_t xyT1, xyT2, xyK1, xyK2;
    ec_basis_t bas1, bas2;
    copy_point(&bas1.P, &T1->P1);
    copy_point(&bas1.Q, &T2->P1);
    copy_point(&bas1.PmQ, &T1m2->P1);

    copy_point(&bas2.P, &T1->P2);
    copy_point(&bas2.Q, &T2->P2);
    copy_point(&bas2.PmQ, &T1m2->P2);

    lift_basis(&xyT1.P1, &xyT2.P1, &bas1, &E12->E1);
    lift_basis(&xyT1.P2, &xyT2.P2, &bas2, &E12->E2);
    TOC_clock(t, "lifting the two basis");

    t = tic();

    // prepare the kernel of the first step
    double_couple_jac_point_iter(&xyK1, n - 1, E12, &xyT1);
    double_couple_jac_point_iter(&xyK2, n - 1, E12, &xyT2);

    TOC_clock(t, "4 xyz doubling");

    // compute the gluing isogeny
    t = tic();
    gluing_comput(&out->first_step, E12, &xyK1, &xyK2);
    TOC_clock(t, "gluing comput");

    // set-up the theta_structure for the first codomain
    t = tic();
    codomain.null_point = out->first_step.codomain;
    codomain.precomputation = 0;
    theta_precomputation(&codomain);
    TOC_clock(t, "precomputation for the first codomain");

    // push the kernel through the gluing isogeny
    // need to setup the input before
    ibz_pow(&a, &ibz_const_two, n);
    ibz_set(&b, 0);
    t = tic();
    gluing_eval_basis(&Q1,
                      &Q2,
                      // T1,T2,T1m2,&a,&b,
                      &xyT1,
                      &xyT2,
                      E12,
                      &out->first_step);
    TOC_clock(t, "gluing eval");
    t = tic();

    // now we launch the evaluation of the n-3 with a strategy other steps
    // setting_up the kernel
    double_iter(&R1, &codomain, &Q1, 2);
    double_iter(&R2, &codomain, &Q2, 2);
    // setting up other parameters
    long log, len = n - 3;
    for (log = 0; len > 1; len >>= 1)
        log++;
    theta_point_t stack1[10 * log + 1];
    theta_point_t stack2[10 * log + 1];
    stack1[0] = Q1;
    stack2[0] = Q2;
    theta_chain_comput_rec(out, &codomain, &R1, &R2, n - 3, 0, false, stack1, stack2, 1, n);
    Q1 = stack1[0];
    Q2 = stack2[0];

    TOC_clock(t, "middle steps");
    t = tic();

    // and now we do the remaining steps
    for (int i = n - 3; i < n - 1; i++) {
        // computing the kernel of the next step
        double_iter(&R1, &codomain, &Q1, n - i - 2);
        double_iter(&R2, &codomain, &Q2, n - i - 2);

        // computing the next step
        if (i == n - 3) {
            theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 0, 0);
        } else if (i == n - 2) {
            theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 1, 0);
        } else {
            theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 0, 1);
        }

        // updating the codomain
        codomain = out->steps[i].codomain;

        // pushing the kernel
        if (i < n - 2) {
            theta_isogeny_eval(&Q1, &out->steps[i], &Q1);
            theta_isogeny_eval(&Q2, &out->steps[i], &Q2);
        }
    }
    TOC_clock(t, "last two steps");

    t = tic();

    // final splitting step
    int is_split = splitting_comput(&out->last_step, &out->steps[n - 2].codomain);
    assert(is_split);

    // computing the curves of the codomain
    theta_product_structure_to_elliptic_product(&out->codomain, &out->last_step.B);
    TOC_clock(t, "splitting");

    ibz_finalize(&a);
    ibz_finalize(&b);
}

int
theta_chain_comput_strategy(theta_chain_t *out,
                            int n,
                            theta_couple_curve_t *E12,
                            const theta_couple_point_t *T1,
                            const theta_couple_point_t *T2,
                            const theta_couple_point_t *T1m2,
                            int *strategy,
                            int eight_above)
{

    theta_couple_point_t P1, P2, P1m2;
    theta_point_t R1, R2;
    theta_isogeny_t steps[n - 1];
    theta_structure_t codomain;

    // init of the isogeny chain
    out->domain = *E12;
    out->length = n;
    out->T1 = *T1;
    out->T2 = *T2;
    out->steps = malloc((n - 1) * sizeof(theta_isogeny_t));

    clock_t t = tic();
    // lift the basis
    theta_couple_jac_point_t xyT1, xyT2, xyK1, xyK2;
    ec_basis_t bas1, bas2;
    copy_point(&bas1.P, &T1->P1);
    copy_point(&bas1.Q, &T2->P1);
    copy_point(&bas1.PmQ, &T1m2->P1);

    copy_point(&bas2.P, &T1->P2);
    copy_point(&bas2.Q, &T2->P2);
    copy_point(&bas2.PmQ, &T1m2->P2);

    if (eight_above) {
        assert(test_point_order_twof(&bas2.P, &E12->E2, n + 2));
        assert(test_point_order_twof(&bas2.Q, &E12->E2, n + 2));
    } else {
        assert(test_point_order_twof(&bas2.P, &E12->E2, n));
        assert(test_point_order_twof(&bas2.Q, &E12->E2, n));
    }

    // t = tic();
    lift_basis(&xyT1.P1, &xyT2.P1, &bas1, &E12->E1);
    lift_basis(&xyT1.P2, &xyT2.P2, &bas2, &E12->E2);
    if (eight_above) {
        assert(test_jac_order_twof(&xyT1.P2, &E12->E2, n + 2));
    } else {
        assert(test_jac_order_twof(&xyT1.P2, &E12->E2, n));
    }

    // TOC_clock(t,"lifting the two basis");

    int adjusting = 2 * (1 - eight_above);
    assert((eight_above && adjusting == 0) || (!eight_above && adjusting == 2));

    // t = tic();

    // prepare the points for the first step
    // first we must compute the length of the list
    int len_count = 0;
    int index = 0;
    while (len_count != n - 1 - adjusting && index < n + 10) {
        len_count = len_count + strategy[index];
        index++;
    }
    int len_list = index + 1;

    theta_couple_jac_point_t points1[n];
    theta_couple_jac_point_t points2[n];
    theta_point_t Q1[n];
    theta_point_t Q2[n];
    int level[n];

    // the first point
    copy_jac_point(&points1[0].P1, &xyT1.P1);
    copy_jac_point(&points1[0].P2, &xyT1.P2);
    copy_jac_point(&points2[0].P1, &xyT2.P1);
    copy_jac_point(&points2[0].P2, &xyT2.P2);
    level[0] = 0;
    // and then the rest
    // t= tic();
    for (int i = 1; i < len_list; i++) {
        double_couple_jac_point_iter(&points1[i], strategy[i - 1], E12, &points1[i - 1]);
        double_couple_jac_point_iter(&points2[i], strategy[i - 1], E12, &points2[i - 1]);
        level[i] = strategy[i - 1];
    }

    // prepare the kernel of the first step
    copy_jac_point(&xyK1.P1, &points1[len_list - 1].P1);
    copy_jac_point(&xyK1.P2, &points1[len_list - 1].P2);
    copy_jac_point(&xyK2.P1, &points2[len_list - 1].P1);
    copy_jac_point(&xyK2.P2, &points2[len_list - 1].P2);
    // TOC_clock(t,"lifting + xyz doubling");

    assert(test_jac_order_twof(&xyK1.P2, &E12->E2, 3));
    assert(test_jac_order_twof(&xyK2.P2, &E12->E2, 3));

    // compute the gluing isogeny
    t = tic();
    gluing_comput(&out->first_step, E12, &xyK1, &xyK2);
    // TOC_clock(t,"gluing comput");

    // set-up the theta_structure for the first codomain
    codomain.null_point = out->first_step.codomain;
    codomain.precomputation = 0;
    theta_precomputation(&codomain);

    len_list--;

    // push the kernel through the gluing isogeny
    // need to setup the input before
    // t = tic();
    for (int i = 0; i < len_list; i++) {
        gluing_eval_basis(&Q1[i], &Q2[i], &points1[i], &points2[i], E12, &out->first_step);
    }

    // TOC_clock(t,"gluing eval");
    // t = tic();
    // and now we do the remaining steps

    for (int i = 0; i < n - 1 - adjusting; i++) {

        len_count = 0;
        for (int j = 0; j < len_list; j++) {
            len_count = len_count + level[j];
        }
        while (len_count != n - i - 2 - adjusting) {
            len_count = len_count + strategy[index];
            double_iter(&Q1[len_list], &codomain, &Q1[len_list - 1], strategy[index]);
            double_iter(&Q2[len_list], &codomain, &Q2[len_list - 1], strategy[index]);
            level[len_list] = strategy[index];
            index++;
            len_list++;
        }

        // computing the kernel of the next step
        // double_iter(&R1,&codomain,&Q1,n-i-2);
        // double_iter(&R2,&codomain,&Q2,n-i-2);
        // TODO : proper copying here ?
        R1 = Q1[len_list - 1];
        R2 = Q2[len_list - 1];

        // computing the next step
        if (i == n - 3) {
            if (eight_above) {
                theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 0, 0);
            }

        } else if (i == n - 2) {
            if (eight_above) {
                theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 1, 0);
            }

        } else {
            theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 0, 1);
        }
        // updating the codomain
        codomain = out->steps[i].codomain;

        len_list--;

        // pushing the kernel
        if (i < n - 2) {
            for (int j = 0; j < len_list; j++) {
                theta_isogeny_eval(&Q1[j], &out->steps[i], &Q1[j]);
                theta_isogeny_eval(&Q2[j], &out->steps[i], &Q2[j]);
            }
        }
    }

    if (!eight_above) {
        // the last two steps are done here
        R1 = Q1[0];
        R2 = Q2[0];
        theta_isogeny_eval(&R1, &out->steps[n - 4], &R1);
        theta_isogeny_eval(&R2, &out->steps[n - 4], &R2);
        theta_isogeny_comput4(&out->steps[n - 3], &codomain, &R1, &R2, 0, 0);
        codomain = out->steps[n - 3].codomain;
        theta_isogeny_eval(&R1, &out->steps[n - 3], &R1);
        theta_isogeny_eval(&R2, &out->steps[n - 3], &R2);
        theta_isogeny_comput2(&out->steps[n - 2], &codomain, &R1, &R2, 1, 0);
        codomain = out->steps[n - 2].codomain;
    }

    // TOC_clock(t,"middle steps");

    t = tic();

    // final splitting step
    int is_split = splitting_comput(&out->last_step, &out->steps[n - 2].codomain);
    if (!is_split) {
        return 0;
    }

    // computing the curves of the codomain
    theta_product_structure_to_elliptic_product(&out->codomain, &out->last_step.B);
    // TOC_clock(t,"splitting");
    return 1;
}

int
theta_chain_comput_strategy_superglue(theta_chain_t *out,
                            int n,
                            theta_couple_curve_t *E12,
                            const theta_couple_point_t *T1,
                            const theta_couple_point_t *T2,
                            const theta_couple_point_t *T1m2,
                            int *strategy,
                            int eight_above)
{

    theta_couple_point_t P1, P2, P1m2;
    theta_point_t R1, R2;
    theta_isogeny_t steps[n - 1];
    theta_structure_t codomain;

    // init of the isogeny chain
    out->domain = *E12;
    out->length = n;
    out->T1 = *T1;
    out->T2 = *T2;
    out->steps = malloc((n - 1) * sizeof(theta_isogeny_t));

    clock_t t = tic();
    // lift the basis
    theta_couple_jac_point_t xyT1, xyT2, xyK1, xyK2;
    ec_basis_t bas1, bas2;
    copy_point(&bas1.P, &T1->P1);
    copy_point(&bas1.Q, &T2->P1);
    copy_point(&bas1.PmQ, &T1m2->P1);

    copy_point(&bas2.P, &T1->P2);
    copy_point(&bas2.Q, &T2->P2);
    copy_point(&bas2.PmQ, &T1m2->P2);

    if (eight_above) {
        assert(test_point_order_twof(&bas2.P, &E12->E2, n + 2));
        assert(test_point_order_twof(&bas2.Q, &E12->E2, n + 2));
    } else {
        assert(test_point_order_twof(&bas2.P, &E12->E2, n));
        assert(test_point_order_twof(&bas2.Q, &E12->E2, n));
    }

    // t = tic();
    lift_basis(&xyT1.P1, &xyT2.P1, &bas1, &E12->E1);
    lift_basis(&xyT1.P2, &xyT2.P2, &bas2, &E12->E2);
    if (eight_above) {
        assert(test_jac_order_twof(&xyT1.P2, &E12->E2, n + 2));
    } else {
        assert(test_jac_order_twof(&xyT1.P2, &E12->E2, n));
    }

    // TOC_clock(t,"lifting the two basis");

    int adjusting = 2 * (1 - eight_above);
    assert((eight_above && adjusting == 0) || (!eight_above && adjusting == 2));

    // t = tic();

    // prepare the points for the first step
    // first we must compute the length of the list
    int len_count = 0;
    int index = 0;
    while (len_count != n - 1 - adjusting && index < n + 10) {
        len_count = len_count + strategy[index];
        index++;
    }
    int len_list = index + 1;

    theta_couple_jac_point_t points1[n];
    theta_couple_jac_point_t points2[n];
    theta_point_t Q1[n];
    theta_point_t Q2[n];
    int level[n];

    // the first point
    copy_jac_point(&points1[0].P1, &xyT1.P1);
    copy_jac_point(&points1[0].P2, &xyT1.P2);
    copy_jac_point(&points2[0].P1, &xyT2.P1);
    copy_jac_point(&points2[0].P2, &xyT2.P2);
    level[0] = 0;
    // and then the rest
    // t= tic();
    for (int i = 1; i < len_list; i++) {
        double_couple_jac_point_iter(&points1[i], strategy[i - 1], E12, &points1[i - 1]);
        double_couple_jac_point_iter(&points2[i], strategy[i - 1], E12, &points2[i - 1]);
        level[i] = strategy[i - 1];
    }

    // prepare the kernel of the first step
    copy_jac_point(&xyK1.P1, &points1[len_list - 1].P1);
    copy_jac_point(&xyK1.P2, &points1[len_list - 1].P2);
    copy_jac_point(&xyK2.P1, &points2[len_list - 1].P1);
    copy_jac_point(&xyK2.P2, &points2[len_list - 1].P2);
    // TOC_clock(t,"lifting + xyz doubling");

    assert(test_jac_order_twof(&xyK1.P2, &E12->E2, 3));
    assert(test_jac_order_twof(&xyK2.P2, &E12->E2, 3));

    // compute the gluing isogeny
    t = tic();
    superglue_gluing_compute(&out->first_step, E12, &xyK1, &xyK2, 1);
    // gluing_comput(&out->first_step, E12, &xyK1, &xyK2);
    // TOC_clock(t,"gluing comput");

    // set-up the theta_structure for the first codomain
    codomain.null_point = out->first_step.codomain;
    codomain.precomputation = 0;
    theta_precomputation(&codomain);

    len_list--;

    // push the kernel through the gluing isogeny
    // need to setup the input before
    // t = tic();
    fp2_t FP2_ONE;
    fp2_set_one(&FP2_ONE);  // TODO Why not precompute this constant ? Why ?
    for (int i = 0; i < len_list; i++) {
        superglue_gluing_eval_point_jac(&Q1[i], &points1[i], &out->first_step, 0, &FP2_ONE);
        superglue_gluing_eval_point_jac(&Q2[i], &points2[i], &out->first_step, 0, &FP2_ONE);
        // gluing_eval_basis(&Q1[i], &Q2[i], &points1[i], &points2[i], E12, &out->first_step);
    }

    // TOC_clock(t,"gluing eval");
    // t = tic();
    // and now we do the remaining steps

    for (int i = 0; i < n - 1 - adjusting; i++) {

        len_count = 0;
        for (int j = 0; j < len_list; j++) {
            len_count = len_count + level[j];
        }
        while (len_count != n - i - 2 - adjusting) {
            len_count = len_count + strategy[index];
            double_iter(&Q1[len_list], &codomain, &Q1[len_list - 1], strategy[index]);
            double_iter(&Q2[len_list], &codomain, &Q2[len_list - 1], strategy[index]);
            level[len_list] = strategy[index];
            index++;
            len_list++;
        }

        // computing the kernel of the next step
        // double_iter(&R1,&codomain,&Q1,n-i-2);
        // double_iter(&R2,&codomain,&Q2,n-i-2);
        // TODO : proper copying here ?
        R1 = Q1[len_list - 1];
        R2 = Q2[len_list - 1];

        // computing the next step
        if (i == n - 3) {
            if (eight_above) {
                theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 0, 0);
            }

        } else if (i == n - 2) {
            if (eight_above) {
                theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 1, 0);
            }

        } else {
            theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 0, 1);
        }
        // updating the codomain
        codomain = out->steps[i].codomain;

        len_list--;

        // pushing the kernel
        if (i < n - 2) {
            for (int j = 0; j < len_list; j++) {
                theta_isogeny_eval(&Q1[j], &out->steps[i], &Q1[j]);
                theta_isogeny_eval(&Q2[j], &out->steps[i], &Q2[j]);
            }
        }
    }

    if (!eight_above) {
        // the last two steps are done here
        R1 = Q1[0];
        R2 = Q2[0];
        theta_isogeny_eval(&R1, &out->steps[n - 4], &R1);
        theta_isogeny_eval(&R2, &out->steps[n - 4], &R2);
        theta_isogeny_comput4(&out->steps[n - 3], &codomain, &R1, &R2, 0, 0);
        codomain = out->steps[n - 3].codomain;
        theta_isogeny_eval(&R1, &out->steps[n - 3], &R1);
        theta_isogeny_eval(&R2, &out->steps[n - 3], &R2);
        theta_isogeny_comput2(&out->steps[n - 2], &codomain, &R1, &R2, 1, 0);
        codomain = out->steps[n - 2].codomain;
    }

    // TOC_clock(t,"middle steps");

    t = tic();

    // final splitting step
    int is_split = splitting_comput(&out->last_step, &out->steps[n - 2].codomain);
    // printf("good = %d", is_split);
    if (!is_split) {
        printf("THE CURVE IS NOT SPLIT! \n");
        return 0;
    }
    // out->last_step.split_type = is_split;

    // computing the curves of the codomain
    theta_product_structure_to_elliptic_product(&out->codomain, &out->last_step.B);
    // TOC_clock(t,"splitting");
    return 1;
}

void
theta_chain_comput_strategy_faster_no_eval(theta_chain_t *out,
                                           int n,
                                           theta_couple_curve_t *E12,
                                           const theta_couple_point_t *T1,
                                           const theta_couple_point_t *T2,
                                           const theta_couple_point_t *T1m2,
                                           int *strategy,
                                           int eight_above)
{

    theta_point_t R1, R2;
    theta_isogeny_t steps[n - 1];
    theta_structure_t codomain;

    // init of the isogeny chain
    out->domain = *E12;
    out->length = n;
    out->T1 = *T1;
    out->T2 = *T2;
    out->steps = malloc((n - 1) * sizeof(theta_isogeny_t));

    clock_t t = tic();
    // lift the basis
    theta_couple_jac_point_t xyT1, xyT2, xyK1, xyK2;
    ec_basis_t bas1, bas2;
    copy_point(&bas1.P, &T1->P1);
    copy_point(&bas1.Q, &T2->P1);
    copy_point(&bas1.PmQ, &T1m2->P1);

    copy_point(&bas2.P, &T1->P2);
    copy_point(&bas2.Q, &T2->P2);
    copy_point(&bas2.PmQ, &T1m2->P2);

    if (eight_above) {
        assert(test_point_order_twof(&bas2.P, &E12->E2, n + 2));
    } else {
        assert(test_point_order_twof(&bas2.P, &E12->E2, n));
    }

    // t = tic();
    lift_basis(&xyT1.P1, &xyT2.P1, &bas1, &E12->E1);
    lift_basis(&xyT1.P2, &xyT2.P2, &bas2, &E12->E2);
    if (eight_above) {
        assert(test_jac_order_twof(&xyT1.P2, &E12->E2, n + 2));
    } else {
        assert(test_jac_order_twof(&xyT1.P2, &E12->E2, n));
    }

    int adjusting = 2 * (1 - eight_above);
    assert((eight_above && adjusting == 0) || (!eight_above && adjusting == 2));

    // t = tic();

    // prepare the points for the first step
    // first we must compute the length of the list
    int len_count = 0;
    int index = 0;
    while (len_count != n - 1 - adjusting && index < n + 10) {
        len_count = len_count + strategy[index];
        index++;
    }
    int len_list = index + 1;

    theta_couple_jac_point_t points1[n];
    theta_couple_jac_point_t points2[n];
    theta_point_t Q1[n];
    theta_point_t Q2[n];
    int level[n];

    // the first point
    copy_jac_point(&points1[0].P1, &xyT1.P1);
    copy_jac_point(&points1[0].P2, &xyT1.P2);
    copy_jac_point(&points2[0].P1, &xyT2.P1);
    copy_jac_point(&points2[0].P2, &xyT2.P2);
    level[0] = 0;
    // and then the rest
    // t= tic();
    for (int i = 1; i < len_list; i++) {
        double_couple_jac_point_iter(&points1[i], strategy[i - 1], E12, &points1[i - 1]);
        double_couple_jac_point_iter(&points2[i], strategy[i - 1], E12, &points2[i - 1]);
        level[i] = strategy[i - 1];
    }

    // prepare the kernel of the first step
    copy_jac_point(&xyK1.P1, &points1[len_list - 1].P1);
    copy_jac_point(&xyK1.P2, &points1[len_list - 1].P2);
    copy_jac_point(&xyK2.P1, &points2[len_list - 1].P1);
    copy_jac_point(&xyK2.P2, &points2[len_list - 1].P2);

    assert(test_jac_order_twof(&xyK1.P2, &E12->E2, 3));
    assert(test_jac_order_twof(&xyK1.P1, &E12->E1, 3));

    // compute the gluing isogeny
    t = tic();
    gluing_comput(&out->first_step, E12, &xyK1, &xyK2);
    // TOC_clock(t,"gluing comput");

    // set-up the theta_structure for the first codomain
    codomain.null_point = out->first_step.codomain;
    codomain.precomputation = 0;
    theta_precomputation(&codomain);

    len_list--;

    // push the kernel through the gluing isogeny
    // need to setup the input before
    // t = tic();
    for (int i = 0; i < len_list; i++) {
        gluing_eval_basis(&Q1[i], &Q2[i], &points1[i], &points2[i], E12, &out->first_step);
    }
    // TOC_clock(t,"gluing eval");
    // t = tic();
    // and now we do the remaining steps

    for (int i = 0; i < n - 1 - adjusting; i++) {

        len_count = 0;
        for (int j = 0; j < len_list; j++) {
            len_count = len_count + level[j];
        }
        while (len_count != n - i - 2 - adjusting) {
            len_count = len_count + strategy[index];
            double_iter(&Q1[len_list], &codomain, &Q1[len_list - 1], strategy[index]);
            double_iter(&Q2[len_list], &codomain, &Q2[len_list - 1], strategy[index]);
            level[len_list] = strategy[index];
            index++;
            len_list++;
        }

        // computing the kernel of the next step
        // double_iter(&R1,&codomain,&Q1,n-i-2);
        // double_iter(&R2,&codomain,&Q2,n-i-2);
        // TODO : proper copying here ?
        R1 = Q1[len_list - 1];
        R2 = Q2[len_list - 1];

        // computing the next step
        if (i == n - 3) {
            if (eight_above) {
                theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 0, 0);
            }

        } else if (i == n - 2) {
            if (eight_above) {
                theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 1, 0);
            }

        } else {
            theta_isogeny_comput(&out->steps[i], &codomain, &R1, &R2, 0, 1);
        }
        // updating the codomain
        codomain = out->steps[i].codomain;

        len_list--;

        // pushing the kernel
        if (i < n - 2) {
            for (int j = 0; j < len_list; j++) {
                theta_isogeny_eval(&Q1[j], &out->steps[i], &Q1[j]);
                theta_isogeny_eval(&Q2[j], &out->steps[i], &Q2[j]);
            }
        }
    }

    if (!eight_above) {
        // the last two steps are done here
        R1 = Q1[0];
        R2 = Q2[0];
        theta_isogeny_eval(&R1, &out->steps[n - 4], &R1);
        theta_isogeny_eval(&R2, &out->steps[n - 4], &R2);
        theta_isogeny_comput4(&out->steps[n - 3], &codomain, &R1, &R2, 0, 0);
        codomain = out->steps[n - 3].codomain;
        theta_isogeny_eval(&R1, &out->steps[n - 3], &R1);
        theta_isogeny_eval(&R2, &out->steps[n - 3], &R2);
        theta_isogeny_comput2(&out->steps[n - 2], &codomain, &R1, &R2, 1, 0);
        codomain = out->steps[n - 2].codomain;
    }

    // TOC_clock(t,"middle steps");

    t = tic();

    // final splitting step
    int is_split = splitting_comput(&out->last_step, &out->steps[n - 2].codomain);
    assert(is_split);

    // computing the curves of the codomain
    theta_product_structure_to_elliptic_product(&out->codomain, &out->last_step.B);
    // TOC_clock(t,"splitting");
}

void
theta_point_to_montgomery_point(theta_couple_point_t *P12,
                                const theta_point_t *P,
                                const theta_structure_t *A)
{

    fp2_t temp;

    // P1.X = A.null_point.y * P.x + A.null_point.x * P.y
    // P1.Z = - A.null_point.y * P.x + A.null_point.x * P.y
    fp2_mul(&P12->P1.x, &A->null_point.y, &P->x);
    fp2_mul(&temp, &A->null_point.x, &P->y);
    fp2_sub(&P12->P1.z, &temp, &P12->P1.x);
    fp2_add(&P12->P1.x, &P12->P1.x, &temp);

    // P2.X = A.null_point.t * P.y + A.null_point.y * P.t
    // P2.Z = - A.null_point.t * P.y + A.null_point.y * P.t
    fp2_mul(&P12->P2.x, &A->null_point.t, &P->y);
    fp2_mul(&temp, &A->null_point.y, &P->t);
    fp2_sub(&P12->P2.z, &temp, &P12->P2.x);
    fp2_add(&P12->P2.x, &P12->P2.x, &temp);
}

/**
 * @brief Evaluate a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the image point
 * @param phi : the (2,2) isogeny chain of domain E12
 * @param P12 a couple point on E12,
 * @param Help a couple point on E12
 *
 * phi : E1xE2 -> E3xE4 of kernel
 * P12 in E1xE2
 * out = phi(P12) in E3xE4
 * Help is equal to phi.first_step.K1_4 + P12
 *
 */
void
theta_chain_eval(theta_couple_point_t *out,
                 theta_chain_t *phi,
                 theta_couple_point_t *P12,
                 const theta_couple_point_t *Help)
{

    theta_point_t temp;

    theta_couple_jac_point_t jac_help;

    // first, we apply the gluing
    gluing_eval_point(&temp, P12, Help, &phi->first_step);

    // then, we apply the successive isogenies
    for (int i = 0; i < phi->length - 1; i++) {
        theta_isogeny_eval(&temp, &phi->steps[i], &temp);
    }

    // we send the result to the theta product structure of the codomain
    apply_isomorphism(&temp, &phi->last_step, &temp);

    // finaly the send the result to the elliptic product
    theta_point_to_montgomery_point(out, &temp, &phi->last_step.B);
}

/**
 * @brief Evaluate a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the image point
 * @param phi : the (2,2) isogeny chain of domain E12
 * @param P12 a couple jac point on E12,
 *
 * phi : E1xE2 -> E3xE4 of kernel
 * P12 in E1xE2
 * out = phi(P12) in E3xE4
 *
 */
void
theta_chain_eval_no_help(theta_couple_point_t *out,
                         theta_chain_t *phi,
                         theta_couple_jac_point_t *P12,
                         const theta_couple_curve_t *E12)
{

    theta_point_t temp;

    // first, we apply the gluing
    gluing_eval_point_no_help(&temp, P12, E12, &phi->first_step);

    // then, we apply the successive isogenies
    for (int i = 0; i < phi->length - 1; i++) {
        theta_isogeny_eval(&temp, &phi->steps[i], &temp);
    }

    // we send the result to the theta product structure of the codomain
    apply_isomorphism(&temp, &phi->last_step, &temp);

    // finaly we send the result to the elliptic product
    theta_point_to_montgomery_point(out, &temp, &phi->last_step.B);
}

/**
 * @brief Evaluate a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the image point
 * @param phi : the (2,2) isogeny chain of domain E12
 * @param P12 a couple point on E12 where one of the two point is zero,
 *
 * phi : E1xE2 -> E3xE4 of kernel
 * P12 in E1xE2 with P12 = (P1,0) or (0,P2)
 * out = phi(P12) in E3xE4
 *
 */
void
theta_chain_eval_special_case(theta_couple_point_t *out,
                              theta_chain_t *phi,
                              theta_couple_point_t *P12,
                              const theta_couple_curve_t *E12)
{
    theta_point_t temp;

#ifndef NDEBUG
    assert(ec_is_zero(&P12->P1) || ec_is_zero(&P12->P2));
    assert(phi->first_step.zero_idx == 3);
#endif

    // first, we apply the gluing
    gluing_eval_point_special_case(&temp, P12, &phi->first_step);

    // then, we apply the successive isogenies
    for (int i = 0; i < phi->length - 1; i++) {
        theta_isogeny_eval(&temp, &phi->steps[i], &temp);
    }

    // we send the result to the theta product structure of the codomain
    apply_isomorphism(&temp, &phi->last_step, &temp);

    // finaly we send the result to the elliptic product
    theta_point_to_montgomery_point(out, &temp, &phi->last_step.B);
}

void theta_chain_eval_special_case_superglue(
    theta_couple_point_t *out,
    theta_chain_t *phi,
    theta_couple_proj_point_t *P12,
    const fp2_t* twist_P)
{
    theta_point_t temp;

#ifndef NDEBUG
    // assert(ec_is_zero(&P12->P1) || ec_is_zero(&P12->P2));
    assert((fp2_is_zero(&P12->P1.x) && fp2_is_zero(&P12->P1.z)) || (fp2_is_zero(&P12->P2.x) && fp2_is_zero(&P12->P2.z)));
    assert(phi->first_step.zero_idx == 3);
#endif

    // first, we apply the gluing
    superglue_gluing_eval_point_proj(&temp, P12, &phi->first_step, 1, twist_P);
    // gluing_eval_point_special_case(&temp, P12, &phi->first_step);

    // then, we apply the successive isogenies
    for (int i = 0; i < phi->length - 1; i++) {
        theta_isogeny_eval(&temp, &phi->steps[i], &temp);
    }

    // we send the result to the theta product structure of the codomain
    apply_isomorphism(&temp, &phi->last_step, &temp);

    // finaly we send the result to the elliptic product
    theta_point_to_montgomery_point(out, &temp, &phi->last_step.B);
}

void theta_chain_eval_special_case_superglue_vsq(
    theta_couple_point_t *out,
    theta_chain_t *phi,
    theta_couple_point_t *P12)
{
    theta_point_t temp;

#ifndef NDEBUG
    assert(ec_is_zero(&P12->P1) || ec_is_zero(&P12->P2));
    // assert((fp2_is_zero(&P12->P1.x) && fp2_is_zero(&P12->P1.z)) || (fp2_is_zero(&P12->P2.x) && fp2_is_zero(&P12->P2.z)));
    assert(phi->first_step.zero_idx == 3);
#endif

    // first, we apply the gluing
    superglue_gluing_eval_point_ec(&temp, P12, &phi->first_step);
    // gluing_eval_point_special_case(&temp, P12, &phi->first_step);

    // then, we apply the successive isogenies
    for (int i = 0; i < phi->length - 1; i++) {
        theta_isogeny_eval(&temp, &phi->steps[i], &temp);
    }

    // we send the result to the theta product structure of the codomain
    apply_isomorphism(&temp, &phi->last_step, &temp);

    // finaly we send the result to the elliptic product
    theta_point_to_montgomery_point(out, &temp, &phi->last_step.B);
}
