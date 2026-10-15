/**
 * @file theta_structure.c
 *
 * @note Modified by Sai-Yu Wang on 2026/09/14.
 * Modified the function `theta_precomputation` to compute the null point
 * and store it in the `null_point` member of the `theta_structure_t`.
 *
 * Add function:
 * - theta_precomputation_dual
 */
#include "theta_structure.h"
#include <assert.h>

void
theta_extract_dual_null_point(theta_point_t *dual_null_point, theta_structure_t *A)
{
    // Reminder, inv_dual_null_point is
    // A^{-1} = CzBxDw = (BCD)(xzw)
    // C^{-1} = AzBxDw = (ABD)(xzw)
    // B^{-1} = AxCzDw = (ACD)(xzw)
    // D^{-1} = AxCzBw = (ABC)(xzw)
    // and dbl_data := (Az, Ax, Bx, Bw), so
    // C = AxCzBw = C(ABxzw) = D^{-1}
    // D = AzBxDw = D(ABxzw) = C^{-1}
    // A = AzAxBw = A(ABxzw) = Az * Ax * Bw
    // B = AzBxBw = B(ABxzw) = Az * Bx * Bw

    fp2_t tmp;
    fp2_copy(&dual_null_point->t, &A->inv_dual_null_point.z);
    fp2_copy(&dual_null_point->z, &A->inv_dual_null_point.t);
    fp2_mul(&tmp, &A->dbl_data.x, &A->dbl_data.t);
    fp2_mul(&dual_null_point->x, &tmp, &A->dbl_data.y);
    fp2_mul(&dual_null_point->y, &tmp, &A->dbl_data.z);
}

void
theta_precomputation(theta_structure_t *A)
{
    if (A->precomputation) {
        return;
    }

    theta_point_t dual_null_point;
    theta_extract_dual_null_point(&dual_null_point, A);
    theta_hadamard(&A->null_point, &dual_null_point);
    to_squared_theta(&dual_null_point, &dual_null_point);
    theta_invert_point(&A->inv_sqr_null_point, &dual_null_point);

    A->precomputation = true;
}

void
theta_precomputation_dual(theta_structure_t *A)
{
    if (A->precomputation_dual) {
        return;
    }

    if (!A->precomputation) {
        theta_extract_dual_null_point(&A->null_point, A);
        theta_hadamard(&A->null_point, &A->null_point);
    }

    if (!fp2_is_zero(&A->null_point.t)) {
        theta_invert_point(&A->inv_null_point, &A->null_point);
    } else {
        fp2_mul(&A->inv_null_point.x, &A->null_point.y, &A->null_point.z);
        fp2_mul(&A->inv_null_point.y, &A->null_point.z, &A->null_point.x);
        fp2_mul(&A->inv_null_point.z, &A->null_point.x, &A->null_point.y);
        fp2_set_zero(&A->inv_null_point.t);
    }

    A->precomputation_dual = true;
}

void
theta_DBL(theta_point_t *out, theta_structure_t *A, const theta_point_t *in)
{
    if (!A->precomputation) {
        theta_precomputation(A);
    }
    to_squared_theta(out, in);
    pointwise_square(out, out);

    pointwise_product(out, out, &A->inv_sqr_null_point);
    theta_hadamard(out, out);
    pointwise_product(out, out, &A->inv_dual_null_point);
}

void
theta_DBL_iter(theta_point_t *out, theta_structure_t *A, const theta_point_t *in, uint16_t exp)
{
    if (exp == 0) {
        *out = *in;
    } else {
        theta_DBL(out, A, in);
        for (uint16_t i = 1; i < exp; i++) {
            theta_DBL(out, A, out);
        }
    }
}

uint32_t
theta_is_product_theta_point(const theta_point_t *P)
{
    fp2_t t1, t2;
    fp2_mul(&t1, &P->x, &P->t);
    fp2_mul(&t2, &P->y, &P->z);
    return fp2_is_equal(&t1, &t2);
}
