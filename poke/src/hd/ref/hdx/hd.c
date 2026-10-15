#include <hd.h>
#include <assert.h>

/**
 * @brief Compute the double of the theta couple point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_point
 * @param E12 an elliptic product
 * @param in the theta couple point in the elliptic product
 * in = (P1,P2)
 * out = [2] (P1,P2)
 *
 */
void
double_couple_point(theta_couple_point_t *out,
                    const theta_couple_curve_t *A,
                    const theta_couple_point_t *in)
{
    ec_dbl(&out->P1, &A->E1, &in->P1);
    ec_dbl(&out->P2, &A->E2, &in->P2);
}

/**
 * @brief Compute the iterated double of the theta couple point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_point
 * @param n : the number of iteration
 * @param E12 an elliptic product
 * @param in the theta couple point in the elliptic product
 * in = (P1,P2)
 * out = [2^n] (P1,P2)
 *
 */
void
double_couple_point_iter(theta_couple_point_t *out,
                         int n,
                         const theta_couple_curve_t *A,
                         const theta_couple_point_t *in)
{
    if (n == 0) {
        *out = *in;
    } else {
        double_couple_point(out, A, in);
        for (int i = 0; i < n - 1; i++) {
            double_couple_point(out, A, out);
        }
    }
}

/**
 * @brief Compute the double of the theta couple jac point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_jac_point
 * @param E12 an elliptic product
 * @param in the theta couple jac point in the elliptic product
 * in = (P1,P2)
 * out = [2] (P1,P2)
 *
 */
void
double_couple_jac_point(theta_couple_jac_point_t *out,
                        const theta_couple_curve_t *A,
                        const theta_couple_jac_point_t *in)
{
    DBL(&out->P1, &in->P1, &A->E1);
    DBL(&out->P2, &in->P2, &A->E2);
}

/**
 * @brief Compute the iterated double of the theta couple jac point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_jac_point
 * @param n : the number of iteration
 * @param E12 an elliptic product
 * @param in the theta couple jac point in the elliptic product
 * in = (P1,P2)
 * out = [2^n] (P1,P2)
 *
 */
void
double_couple_jac_point_iter(theta_couple_jac_point_t *out,
                             int n,
                             const theta_couple_curve_t *A,
                             const theta_couple_jac_point_t *in)
{
    if (n == 0) {
        *out = *in;
    } else {
        double_couple_jac_point(out, A, in);
        for (int i = 0; i < n - 1; i++) {
            double_couple_jac_point(out, A, out);
        }
    }
}

void
couple_jac_to_xz(theta_couple_point_t *P, const theta_couple_jac_point_t *xyP)
{
    jac_to_xz(&P->P1, &xyP->P1);
    jac_to_xz(&P->P2, &xyP->P2);
}

/**
 * @brief Compute the addition of theta couple points in the elliptic product E12
 *
 * @param out Output: the theta_couple_point
 * @param E12 an elliptic product
 * @param T1 a theta couple point in the elliptic product
 * @param T2 a theta couple point in the elliptic product
 * T1 = (P1,P2) T2 = (Q1,Q2)
 * out = (P1+Q1,P2+Q2)
 *
 */
void
add_couple_point(theta_couple_point_t *out,
                 const theta_couple_curve_t *A,
                 const theta_couple_point_t *T1,
                 const theta_couple_point_t *T2)
{

    // perform the lifting
    // TODETERMINE : do we care about signs or not ?
    assert(0);
}

void jac_to_xz_add_components(
    add_components_t *add_comp,
    const jac_point_t *P,
    const jac_point_t *Q,
    const ec_curve_t *AC,
    const fp2_t *twist_P,
    const fp2_t *twist_Q)
{
    // Take P and Q in E distinct, two jac_point_t, return three components u,v and w in Fp2 such
    // that the xz coordinates of P+Q are (u-v:w) and of P-Q are (u+v:w)

    assert(fp2_is_zero(&AC->A) || fp2_is_one(&AC->C));
    fp2_t t0, t1, t2, t3, t4, t5, t6;

    // assert(fp2_is_equal(&twist_P, &twist_Q));

    fp2_sqr(&t0, &P->z);             // t0 = z1^2
    fp2_sqr(&t1, &Q->z);             // t1 = z2^2
    fp2_mul(&t2, &P->x, &t1);        // t2 = x1z2^2
    fp2_mul(&t3, &t0, &Q->x);        // t3 = z1^2x2
    fp2_mul(&t4, &P->y, &Q->z);      // t4 = y1z2
    fp2_mul(&t4, &t4, &t1);          // t4 = y1z2^3
    fp2_mul(&t5, &P->z, &Q->y);      // t5 = z1y2
    fp2_mul(&t5, &t5, &t0);          // t5 = z1^3y2
    fp2_mul(&t0, &t0, &t1);          // t0 = (z1z2)^2
    fp2_mul(&t6, &t4, &t5);          // t6 = (z1z_2)^3y1y2
    fp2_add(&add_comp->v, &t6, &t6); // v  = 2(z1z_2)^3y1y2
    fp2_sqr(&t4, &t4);               // t4 = y1^2z2^6
    fp2_mul(&t4, &t4, twist_P);
    fp2_sqr(&t5, &t5);               // t5 = z1^6y_2^2
    fp2_mul(&t5, &t5, twist_Q);
    fp2_add(&t4, &t4, &t5);          // t4 = z1^6y_2^2 + y1^2z2^6
    fp2_add(&t5, &t2, &t3);          // t5 = x1z2^2 +z_1^2x2
    fp2_add(&t6, &t3, &t3);          // t6 = 2z_1^2x2
    fp2_sub(&t6, &t5, &t6);          // t6 = lambda = x1z2^2 - z_1^2x2
    fp2_sqr(&t6, &t6);               // t6 = lambda^2 = (x1z2^2 - z_1^2x2)^2
    fp2_mul(&t1, &AC->A, &t0);       // t1 = A*(z1z2)^2
    fp2_add(&t1, &t5, &t1);          // t1 = gamma =A*(z1z2)^2 + x1z2^2 +z_1^2x2
    fp2_mul(&t1, &t1, &t6);          // t1 = gamma*lambda^2
    fp2_sub(&add_comp->u, &t4, &t1); // u  = z1^6y_2^2 + y1^2z2^6 - gamma*lambda^2
    fp2_mul(&add_comp->w, &t6, &t0); // w  = (z1z2)^2(lambda)^2
}

void proj_to_xz_add_components(
    add_components_t *add_comp,
    const proj_point_t *P,
    const proj_point_t *Q,
    const ec_curve_t *AC,
    const fp2_t *twist_P,
    const fp2_t *twist_Q)
{
    assert(fp2_is_zero(&AC->A) || fp2_is_one(&AC->C));
    fp2_t t1, t2, t3, t4, t5, t6, t7;

    fp2_mul(&t1, &P->z, &Q->z);         // t1 = z1 * z2
    fp2_mul(&t2, &P->x, &Q->z);         // t2 = x1 * z2
    fp2_mul(&t3, &P->z, &Q->x);         // t3 = z1 * x2
    fp2_mul(&t4, &P->y, &Q->z);         // t4 = y1 * z2
    fp2_mul(&t5, &P->z, &Q->y);         // t5 = z1 * y2
    fp2_mul(&t6, &t4, &t5);             // t6 = y1y2 * z1z2
    fp2_mul(&t6, &t6, &t1);             // t6 = y1y2 * (z1z2)^2
    fp2_add(&add_comp->v, &t6, &t6);    // v = 2 * y1y2 * (z1z2)^2
    fp2_mul(&t7, &AC->A, &t1);          // t7 = A * z1z2
    fp2_add(&t2, &t2, &t3);             // t2 = x1z2 + z1x2
    fp2_add(&t7, &t7, &t2);             // t7 = Az1z2 + x1z2 + z1x2
    fp2_add(&t3, &t3, &t3);             // t3 = 2 * z1x2
    fp2_sub(&t2, &t2, &t3);             // t2 = gamma = x1z2 - z1x2
    fp2_sqr(&t2, &t2);                  // t2 = gamma^2
    fp2_mul(&t7, &t7, &t2);             // t7 = (Az1z2 + x1z2 + z1x2) * gamma^2
    fp2_sqr(&t4, &t4);                  // t4 = y1^2 * z2^2
    fp2_mul(&t4, &t4, twist_P);
    fp2_sqr(&t5, &t5);                  // t5 = z1^2 * y2^2
    fp2_mul(&t5, &t5, twist_Q);
    fp2_add(&t4, &t4, &t5);             // t4 = y1^2 * z2^2 + z1^2 * y2^2
    fp2_mul(&t4, &t4, &t1);             // t4 = (y1^2 * z2^2 + z1^2 * y2^2) * z1z2
    fp2_sub(&add_comp->u, &t4, &t7);    // u = (y1^2 * z2^2 + z1^2 * y2^2) * z1z2 - (Az1z2 + x1z2 + z1x2) * gamma^2
    fp2_mul(&add_comp->w, &t2, &t1);    // w = gamma^2 * z1z2
}

void ec_to_xz_add_components(
    add_components_t *add_comp,
    const ec_point_t *P,
    const ec_point_t *Q,
    const ec_curve_t *AC)
{
    assert(fp2_is_zero(&AC->A) || fp2_is_one(&AC->C));
    fp2_t t1, t2;
    fp2_t alpha, beta, gamma, delta;

    fp2_mul(&t1, &P->x, &Q->x);             // t1 = x1 * x2
    fp2_mul(&t2, &P->z, &Q->z);             // t2 = z1 * z2
    fp2_add(&alpha, &t1, &t2);              // alpha = x1x2 + z1z2
    fp2_mul(&t1, &P->x, &Q->z);             // t1 = x1 * z2
    fp2_mul(&t2, &P->z, &Q->x);             // t2 = x2 * z1
    fp2_add(&beta, &t1, &t2);               // beta = x1z2 + x2z1
    fp2_sub(&gamma, &t1, &t2);              // gamma = x1z2 - x2z1
    fp2_sqr(&add_comp->w, &gamma);          // w = gamma^2
    fp2_mul(&t1, &t1, &t2);                 // t1 = xxzz
    fp2_mul(&t2, &t1, &AC->A);              // t2 = A*xxzz
    fp2_mul(&delta, &alpha, &beta);
    fp2_add(&delta, &delta, &t2);           // delta = alpha * beta + A*xxzz
    fp2_add(&add_comp->u, &delta, &t2);     // u = delta + A*xxzz
    fp2_mul(&t2, &delta, &AC->A);           // t2 = A * delta
    fp2_sqr(&add_comp->v, &alpha);
    fp2_add(&add_comp->v, &add_comp->v, &add_comp->w);
    fp2_add(&add_comp->v, &add_comp->v, &t2);
    fp2_mul(&add_comp->v, &add_comp->v, &t1);
    fp2_add(&add_comp->v, &add_comp->v, &add_comp->v);
    fp2_add(&add_comp->v, &add_comp->v, &add_comp->v);      // v^2 = 4 * xxzz * (alpha^2 + gamma^2 + A * delta)
}