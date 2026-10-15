#include "curve_extras.h"
#include "ec.h"
#include "fp2.h"
#include "gf_constants.h"

void
xTPL(ec_point_t *Q, const ec_point_t *P, const ec_point_t *A3)
{
    /* ----------------------------------------------------------------------------- *
     * Differential point tripling given the montgomery coefficient A3 = (A+2C:A-2C)
     * ----------------------------------------------------------------------------- */

    fp2_t t0, t1, t2, t3, t4;
    fp2_sub(&t0, &P->x, &P->z);
    fp2_sqr(&t2, &t0);
    fp2_add(&t1, &P->x, &P->z);
    fp2_sqr(&t3, &t1);
    fp2_add(&t4, &t1, &t0);
    fp2_sub(&t0, &t1, &t0);
    fp2_sqr(&t1, &t4);
    fp2_sub(&t1, &t1, &t3);
    fp2_sub(&t1, &t1, &t2);
    fp2_mul(&Q->x, &t3, &A3->x);
    fp2_mul(&t3, &Q->x, &t3);
    fp2_mul(&Q->z, &t2, &A3->z);
    fp2_mul(&t2, &t2, &Q->z);
    fp2_sub(&t3, &t2, &t3);
    fp2_sub(&t2, &Q->x, &Q->z);
    fp2_mul(&t1, &t2, &t1);
    fp2_add(&t2, &t3, &t1);
    fp2_sqr(&t2, &t2);
    fp2_mul(&Q->x, &t2, &t4);
    fp2_sub(&t1, &t3, &t1);
    fp2_sqr(&t1, &t1);
    fp2_mul(&Q->z, &t1, &t0);
}

int
ec_is_on_curve(const ec_curve_t *curve, const ec_point_t *P)
{

    fp2_t t0, t1, t2;

    // Check if xz*(C^2x^2+zACx+z^2C^2) is a square
    fp2_mul(&t0, &curve->C, &P->x);
    fp2_mul(&t1, &t0, &P->z);
    fp2_mul(&t1, &t1, &curve->A);
    fp2_mul(&t2, &curve->C, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sqr(&t2, &t2);
    fp2_add(&t0, &t0, &t1);
    fp2_add(&t0, &t0, &t2);
    fp2_mul(&t0, &t0, &P->x);
    fp2_mul(&t0, &t0, &P->z);
    return fp2_is_square(&t0);
}

