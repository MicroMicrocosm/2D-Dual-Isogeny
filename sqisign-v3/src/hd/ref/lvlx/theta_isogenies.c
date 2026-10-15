/**
 * @file theta_isogenies.c
 *
 * @note Modified by Sai-Yu Wang on 2026/09/13.
 * Add parameter `Phi` to the function `theta_chain_compute_impl` to store the isogeny chain.
 *
 * Add functions:
 * - splitting_matrix_inv
 * - gluing_matrix_inv
 * - theta_isogeny_eval_dual
 * - theta_chain_compute_steps
 * - theta_chain_eval_special_forward
 * - theta_chain_eval_special_dual
 * - theta_chain_eval_special_forward_then_dual
 * - theta_chain_eval_special_dual_then_forward
 *
 * @note Modified by Sai-Yu Wang on 2026/09/18.
 * Rename function `theta_chain_compute_and_eval_verify` to `theta_chain_compute_and_eval_verify_E1`
 * and add function `theta_chain_compute_and_eval_verify_E2`
 * to compute and evaluate the theta chain for verification of E1 and E2 respectively.
 */

#include <hd.h>
#include "theta_structure.h"
#include "gluing.h"
#include "splitting.h"

/**
 * @brief Compute the projective inverse of a 2x2 change of coordinate matrix for splitting.
 * @param M_inv The inverse of the change of coordinate matrix M
 * @param M The change of coordinate matrix
 */
static inline void
splitting_matrix_inv(ec_change_coord_matrix_t *M_inv, const ec_change_coord_matrix_t *M)
{
    fp2_copy(&M_inv->a, &M->d);
    fp2_neg(&M_inv->b, &M->b);
    fp2_neg(&M_inv->c, &M->c);
    fp2_copy(&M_inv->d, &M->a);
}

/**
 * @brief Compute the projective inverse of a 4x4 change of coordinate matrix for gluing.
 * @param M_inv The inverse of the change of coordinate matrix M
 * @param M The change of coordinate matrix
 */
static void
gluing_matrix_inv(gluing_change_coord_matrix_t *M_inv, const gluing_change_coord_matrix_t *M)
{
    fp2_t x0, y0, z0, t0;
    fp2_t det_01_01, det_01_02, det_01_03, det_01_12, det_01_13, det_01_23;
    fp2_t det_23_01, det_23_02, det_23_03, det_23_12, det_23_13, det_23_23;

    fp2_mul(&x0, &M->m[0][0], &M->m[1][1]);
    fp2_mul(&y0, &M->m[1][0], &M->m[0][1]);
    fp2_sub(&det_01_01, &x0, &y0);
    fp2_mul(&x0, &M->m[0][0], &M->m[1][2]);
    fp2_mul(&y0, &M->m[1][0], &M->m[0][2]);
    fp2_sub(&det_01_02, &x0, &y0);
    fp2_mul(&x0, &M->m[0][0], &M->m[1][3]);
    fp2_mul(&y0, &M->m[1][0], &M->m[0][3]);
    fp2_sub(&det_01_03, &x0, &y0);
    fp2_mul(&x0, &M->m[0][1], &M->m[1][2]);
    fp2_mul(&y0, &M->m[1][1], &M->m[0][2]);
    fp2_sub(&det_01_12, &x0, &y0);
    fp2_mul(&x0, &M->m[0][1], &M->m[1][3]);
    fp2_mul(&y0, &M->m[1][1], &M->m[0][3]);
    fp2_sub(&det_01_13, &x0, &y0);
    fp2_mul(&x0, &M->m[0][2], &M->m[1][3]);
    fp2_mul(&y0, &M->m[1][2], &M->m[0][3]);
    fp2_sub(&det_01_23, &x0, &y0);

    fp2_mul(&x0, &M->m[2][0], &M->m[3][1]);
    fp2_mul(&y0, &M->m[3][0], &M->m[2][1]);
    fp2_sub(&det_23_01, &x0, &y0);
    fp2_mul(&x0, &M->m[2][0], &M->m[3][2]);
    fp2_mul(&y0, &M->m[3][0], &M->m[2][2]);
    fp2_sub(&det_23_02, &x0, &y0);
    fp2_mul(&x0, &M->m[2][0], &M->m[3][3]);
    fp2_mul(&y0, &M->m[3][0], &M->m[2][3]);
    fp2_sub(&det_23_03, &x0, &y0);
    fp2_mul(&x0, &M->m[2][1], &M->m[3][2]);
    fp2_mul(&y0, &M->m[3][1], &M->m[2][2]);
    fp2_sub(&det_23_12, &x0, &y0);
    fp2_mul(&x0, &M->m[2][1], &M->m[3][3]);
    fp2_mul(&y0, &M->m[3][1], &M->m[2][3]);
    fp2_sub(&det_23_13, &x0, &y0);
    fp2_mul(&x0, &M->m[2][2], &M->m[3][3]);
    fp2_mul(&y0, &M->m[3][2], &M->m[2][3]);
    fp2_sub(&det_23_23, &x0, &y0);

    fp2_mul(&x0, &M->m[1][1], &det_23_23);
    fp2_mul(&y0, &M->m[1][2], &det_23_13);
    fp2_mul(&z0, &M->m[1][3], &det_23_12);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[0][0], &t0, &y0);
    fp2_mul(&x0, &M->m[1][0], &det_23_23);
    fp2_mul(&y0, &M->m[1][2], &det_23_03);
    fp2_mul(&z0, &M->m[1][3], &det_23_02);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[1][0], &y0, &t0);
    fp2_mul(&x0, &M->m[1][0], &det_23_13);
    fp2_mul(&y0, &M->m[1][1], &det_23_03);
    fp2_mul(&z0, &M->m[1][3], &det_23_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[2][0], &t0, &y0);
    fp2_mul(&x0, &M->m[1][0], &det_23_12);
    fp2_mul(&y0, &M->m[1][1], &det_23_02);
    fp2_mul(&z0, &M->m[1][2], &det_23_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[3][0], &y0, &t0);

    fp2_mul(&x0, &M->m[0][1], &det_23_23);
    fp2_mul(&y0, &M->m[0][2], &det_23_13);
    fp2_mul(&z0, &M->m[0][3], &det_23_12);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[0][1], &y0, &t0);
    fp2_mul(&x0, &M->m[0][0], &det_23_23);
    fp2_mul(&y0, &M->m[0][2], &det_23_03);
    fp2_mul(&z0, &M->m[0][3], &det_23_02);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[1][1], &t0, &y0);
    fp2_mul(&x0, &M->m[0][0], &det_23_13);
    fp2_mul(&y0, &M->m[0][1], &det_23_03);
    fp2_mul(&z0, &M->m[0][3], &det_23_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[2][1], &y0, &t0);
    fp2_mul(&x0, &M->m[0][0], &det_23_12);
    fp2_mul(&y0, &M->m[0][1], &det_23_02);
    fp2_mul(&z0, &M->m[0][2], &det_23_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[3][1], &t0, &y0);

    fp2_mul(&x0, &M->m[3][1], &det_01_23);
    fp2_mul(&y0, &M->m[3][2], &det_01_13);
    fp2_mul(&z0, &M->m[3][3], &det_01_12);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[0][2], &t0, &y0);
    fp2_mul(&x0, &M->m[3][0], &det_01_23);
    fp2_mul(&y0, &M->m[3][2], &det_01_03);
    fp2_mul(&z0, &M->m[3][3], &det_01_02);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[1][2], &y0, &t0);
    fp2_mul(&x0, &M->m[3][0], &det_01_13);
    fp2_mul(&y0, &M->m[3][1], &det_01_03);
    fp2_mul(&z0, &M->m[3][3], &det_01_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[2][2], &t0, &y0);
    fp2_mul(&x0, &M->m[3][0], &det_01_12);
    fp2_mul(&y0, &M->m[3][1], &det_01_02);
    fp2_mul(&z0, &M->m[3][2], &det_01_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[3][2], &y0, &t0);

    fp2_mul(&x0, &M->m[2][1], &det_01_23);
    fp2_mul(&y0, &M->m[2][2], &det_01_13);
    fp2_mul(&z0, &M->m[2][3], &det_01_12);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[0][3], &y0, &t0);
    fp2_mul(&x0, &M->m[2][0], &det_01_23);
    fp2_mul(&y0, &M->m[2][2], &det_01_03);
    fp2_mul(&z0, &M->m[2][3], &det_01_02);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[1][3], &t0, &y0);
    fp2_mul(&x0, &M->m[2][0], &det_01_13);
    fp2_mul(&y0, &M->m[2][1], &det_01_03);
    fp2_mul(&z0, &M->m[2][3], &det_01_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[2][3], &y0, &t0);
    fp2_mul(&x0, &M->m[2][0], &det_01_12);
    fp2_mul(&y0, &M->m[2][1], &det_01_02);
    fp2_mul(&z0, &M->m[2][2], &det_01_01);
    fp2_add(&t0, &x0, &z0);
    fp2_sub(&M_inv->m[3][3], &t0, &y0);
}

/**
 * @brief Compute a (2,2) isogeny in dimension 2 in the theta_model
 *
 * @param out Output: the theta_isogeny
 * @param T1_8 a point in A[8]
 * @param T2_8 a point in A[8]
 * @param dual_domain_bool a boolean used for the last two steps of the chain
 *
 * out : A -> B of kernel [4](T1_8,T2_8)
 * theta_hadamard_bool_1 controls if the domain is in standard or dual coordinates
 * theta_hadamard_bool_2 controls if the codomain is in standard or dual coordinates
 * verify: add extra sanity check to ensure our 8-torsion points are coherent with the isogeny
 *
 */
static int
theta_isogeny_compute(theta_isogeny_t *out,
                      const theta_point_t *T1_8,
                      const theta_point_t *T2_8,
                      const bool dual_domain_bool,
                      const bool verify)
{
    out->codomain.precomputation = false;
    out->codomain.precomputation_dual = false;

    theta_point_t TT1, TT2;

    if (dual_domain_bool) {
        theta_hadamard(&TT1, T1_8);
        to_squared_theta(&TT1, &TT1);
        theta_hadamard(&TT2, T2_8);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT1, T1_8);
        to_squared_theta(&TT2, T2_8);
    }

    // Test that our projective factor ABCDxzw is non zero, where
    // TT1=(Ax, Bx, Cy, Dy), TT2=(Az, Bw, Cz, Dw)
    // But ABCDxzw=0 can only happen if we had an unexpected splitting in
    // the isogeny chain. In either case reject
    // (this is not strictly necessary, we could just return (0:0:0:0))
    if (verify && (fp2_is_zero(&TT2.x) | fp2_is_zero(&TT2.y) | fp2_is_zero(&TT2.z) | fp2_is_zero(&TT2.t) |
                   fp2_is_zero(&TT1.x) | fp2_is_zero(&TT1.y))) {
        return 0;
    }
    fp2_t t1, t2;

    // compute inverse dual theta null point
    fp2_mul(&t1, &TT1.y, &TT2.t); // BxDw
    fp2_mul(&t2, &TT1.x, &TT2.z); // AxCz

    fp2_mul(&out->codomain.inv_dual_null_point.x, &t1, &TT2.z); // A^{-1} = CzBxDw = (BCD)(xzw)
    fp2_mul(&out->codomain.inv_dual_null_point.z, &t1, &TT2.x); // C^{-1} = AzBxDw = (ABD)(xzw)
    fp2_mul(&out->codomain.inv_dual_null_point.y, &t2, &TT2.t); // B^{-1} = AxCzDw = (ACD)(xzw)
    fp2_mul(&out->codomain.inv_dual_null_point.t, &t2, &TT2.y); // D^{-1} = AxCzBw = (ABC)(xzw)

    // Keep in memory intermediary values to compute dual theta null point faster if necessary.
    // C = AxCzBw = C(ABxzw) = inv_dual_null_point.t
    // D = AzBxDw = D(ABxzw) = inv_dual_null_point.z
    // A = AzAxBw = A(ABxzw) = TT2.x * TT1.x * TT2.y
    // B = AzBxBw = B(ABxzw) = TT2.x * TT1.y * TT2.y
    fp2_copy(&out->codomain.dbl_data.x, &TT2.x);
    fp2_copy(&out->codomain.dbl_data.y, &TT1.x);
    fp2_copy(&out->codomain.dbl_data.z, &TT1.y);
    fp2_copy(&out->codomain.dbl_data.t, &TT2.y);

    // If T1_8 and T2_8 are our 8-torsion points, this ensures that the
    // 4-torsion points 2T1_8 and 2T2_8 are isotropic.
    // All 3 others equalities are a natural consequence of the interpolation formula.

    if (verify) {
        fp2_mul(&t1, &TT1.z, &out->codomain.inv_dual_null_point.z);
        fp2_mul(&t2, &TT1.t, &out->codomain.inv_dual_null_point.t);
        if (!fp2_is_equal(&t1, &t2))
            return 0;
    }

    return 1;
}

static void
theta_isogeny_eval(theta_point_t *out, const theta_isogeny_t *phi, const theta_point_t *P, const bool dual_domain_bool)
{
    if (dual_domain_bool) {
        theta_hadamard(out, P);
        to_squared_theta(out, out);
    } else {
        to_squared_theta(out, P);
    }
    fp2_mul(&out->x, &out->x, &phi->codomain.inv_dual_null_point.x);
    fp2_mul(&out->y, &out->y, &phi->codomain.inv_dual_null_point.y);
    fp2_mul(&out->z, &out->z, &phi->codomain.inv_dual_null_point.z);
    fp2_mul(&out->t, &out->t, &phi->codomain.inv_dual_null_point.t);
}

static int
theta_isogeny_eval_dual(theta_point_t *out,
                        theta_isogeny_t *phi,
                        const theta_point_t *P,
                        const bool dual_domain_bool,
                        const bool t_is_zero)
{
    if (dual_domain_bool) {
        to_squared_theta(out, P);
    } else {
        theta_hadamard(out, P);
        to_squared_theta(out, out);
    }

    if (t_is_zero && !fp2_is_zero(&out->t)) {
        debug_print("unexpected non-zero t coordinate in theta_isogeny_eval_dual");
        return 0;
    }

    if (!phi->codomain.precomputation_dual) {
        theta_precomputation_dual(&phi->codomain);
    }

    if (t_is_zero && !fp2_is_zero(&phi->codomain.inv_null_point.t)) {
        debug_print("unexpected non-zero t coordinate in phi->codomain.inv_null_point");
        return 0;
    }

    pointwise_product(out, out, &phi->codomain.inv_null_point);

    return 1;
}

static int
theta_chain_compute_impl(uint16_t n,
                         theta_couple_curve_t *E12,
                         const theta_kernel_couple_points_t *ker,
                         theta_couple_curve_t *E34,
                         theta_chain_t *Phi,
                         theta_couple_point_t *P12,
                         uint8_t numP,
                         int8_t mode)
{
    ec_curve_normalize_A24(&E12->E1);
    ec_curve_normalize_A24(&E12->E2);

#ifndef NDEBUG
    if (!test_point_order_twof(&ker->T1.P2, &E12->E2, n + HD_EXTRA_TORSION))
        debug_print("T1.P2 does not have correct order");

    if (!test_point_order_twof(&ker->T1.P1, &E12->E1, n + HD_EXTRA_TORSION))
        debug_print("T1.P1 does not have correct order");
#endif

    // points to evaluate throughout the chain
    theta_point_t pts[numP ? numP : 1];

    // init chain
    if (Phi) {
        Phi->length = n;
        Phi->steps = malloc(n * sizeof(theta_isogeny_t));
    }
    uint8_t space = 1;
    for (uint16_t i = 1; i < n; i *= 2)
        ++space;

    uint16_t todo[space];
    todo[0] = n;

    // kernel points for the remaining isogeny steps
    theta_point_t thetaQ1[space], thetaQ2[space];

    //-------------- GLUING STEP --------------------
    // start new gluing
    theta_gluing_t first_step;

    int8_t current = 0;
    current = gluing_start_chain(&first_step, thetaQ1, thetaQ2, pts, todo, space, numP, E12, P12, ker, mode);
    if (Phi) {
        Phi->first_step = first_step;
    }

    //-------------- GENERIC STEPS --------------------
    theta_isogeny_t step;
    uint32_t ret;
    bool is_not_step2;

    // set-up the generic theta_structure for the chain
    // here, we are not in dual but in regular mode (done to avoid zero in DBL formula)
    theta_invert_point(&step.codomain.inv_dual_null_point, &first_step.codomain);
    to_squared_theta(&step.codomain.inv_sqr_null_point, &first_step.codomain);
    theta_invert_point(&step.codomain.inv_sqr_null_point, &step.codomain.inv_sqr_null_point);
    step.codomain.precomputation = true;
    theta_copy(&step.codomain.null_point, &first_step.codomain);
    theta_copy(&step.codomain.inv_null_point, &step.codomain.inv_dual_null_point);
    step.codomain.precomputation_dual = true;
    if (Phi) {
        Phi->steps[0] = step;
    }

    // The generic steps
    for (uint16_t i = 1; current >= 0 && todo[current]; ++i) {
        assert(current < space);
        is_not_step2 = (i != 1);
        while (todo[current] != 1) {
            assert(todo[current] >= 2);
            ++current;
            assert(current < space);
            const uint16_t num_dbls = todo[current - 1] / 2;
            assert(num_dbls && num_dbls < todo[current - 1]);
            theta_DBL_iter(&thetaQ1[current], &step.codomain, &thetaQ1[current - 1], num_dbls);
            theta_DBL_iter(&thetaQ2[current], &step.codomain, &thetaQ2[current - 1], num_dbls);
            todo[current] = todo[current - 1] - num_dbls;
        }

        // computing the next step
        if (i == n - 1) // ultimate step
        {
            ret = theta_isogeny_compute(&step, &thetaQ1[current], &thetaQ2[current], is_not_step2, false);
        } else {
            ret = theta_isogeny_compute(&step,
                                        &thetaQ1[current],
                                        &thetaQ2[current],
                                        is_not_step2,
                                        (mode == THETA_CHAIN_MODE_VERIFY_E1) || (mode == THETA_CHAIN_MODE_VERIFY_E2));
        }
        if (Phi) {
            Phi->steps[i] = step;
        }
        if (!ret) {
            debug_print("failure of the codomain computation during the HD chain");
            return 0;
        }
        // pushing the points
        for (int j = 0; j < numP; ++j) {
            theta_isogeny_eval(&pts[j], &step, &pts[j], is_not_step2);
        }

        // pushing the kernel
        assert(todo[current] == 1);
        for (int j = 0; j < current; ++j) {
            theta_isogeny_eval(&thetaQ1[j], &step, &thetaQ1[j], is_not_step2);
            theta_isogeny_eval(&thetaQ2[j], &step, &thetaQ2[j], is_not_step2);
            assert(todo[j]);
            --todo[j];
        }
        --current;
    }

    assert(current == -1);

    //-------------- SPLITTING STEP --------------------
    theta_point_t dual_null_pt;
    theta_extract_dual_null_point(&dual_null_pt, &step.codomain);
    if (Phi) {
        ret = splitting_to_elliptic_product(E34, P12, &Phi->M_spl_E1, &Phi->M_spl_E2, &dual_null_pt, pts, numP, mode);
    } else {
        ret = splitting_to_elliptic_product(E34, P12, NULL, NULL, &dual_null_pt, pts, numP, mode);
    }
    if (!ret)
        return 0;

    return 1;
}

int
theta_chain_compute_and_eval(uint16_t n,
                             /*const*/ theta_couple_curve_t *E12,
                             const theta_kernel_couple_points_t *ker,
                             theta_couple_curve_t *E34,
                             theta_couple_point_t *P12,
                             uint8_t numP)
{
    return theta_chain_compute_impl(n, E12, ker, E34, NULL, P12, numP, THETA_CHAIN_MODE_BOTH);
}

int
theta_chain_compute_and_eval_E1(uint16_t n,
                                /*const*/ theta_couple_curve_t *E12,
                                const theta_kernel_couple_points_t *ker,
                                theta_couple_curve_t *E34,
                                theta_couple_point_t *P12,
                                uint8_t numP)
{
    return theta_chain_compute_impl(n, E12, ker, E34, NULL, P12, numP, THETA_CHAIN_MODE_E1);
}

int
theta_chain_compute_and_eval_E2(uint16_t n,
                                /*const*/ theta_couple_curve_t *E12,
                                const theta_kernel_couple_points_t *ker,
                                theta_couple_curve_t *E34,
                                theta_couple_point_t *P12,
                                uint8_t numP)
{
    return theta_chain_compute_impl(n, E12, ker, E34, NULL, P12, numP, THETA_CHAIN_MODE_E2);
}

// Like theta_chain_compute_and_eval, adding extra verification checks;
// used in the signature verification
int
theta_chain_compute_and_eval_verify_E1(uint16_t n,
                                       /*const*/ theta_couple_curve_t *E12,
                                       const theta_kernel_couple_points_t *ker,
                                       theta_couple_curve_t *E34,
                                       theta_couple_point_t *P12,
                                       uint8_t numP)
{
    return theta_chain_compute_impl(n, E12, ker, E34, NULL, P12, numP, THETA_CHAIN_MODE_VERIFY_E1);
}

// Like theta_chain_compute_and_eval, adding extra verification checks;
// used in the signature verification
int
theta_chain_compute_and_eval_verify_E2(uint16_t n,
                                       /*const*/ theta_couple_curve_t *E12,
                                       const theta_kernel_couple_points_t *ker,
                                       theta_couple_curve_t *E34,
                                       theta_couple_point_t *P12,
                                       uint8_t numP)
{
    return theta_chain_compute_impl(n, E12, ker, E34, NULL, P12, numP, THETA_CHAIN_MODE_VERIFY_E2);
}

int
theta_chain_compute_steps(uint16_t n,
                          /*const*/ theta_couple_curve_t *E12,
                          const theta_kernel_couple_points_t *ker,
                          theta_couple_curve_t *E34,
                          theta_chain_t *Phi)
{
    return theta_chain_compute_impl(n, E12, ker, E34, Phi, NULL, 0, 0);
}

/**
 * @brief Evaluate at a list of points of the form (P1,0) or (0,P2).
 *
 * @param Phi a isogeny chain
 * @param P12 Input/Output: pointer to points to be pushed through the isogeny (in-place)
 * @param numP: length of the list of points given in P12 (can be zero)
 * @return 1 on success, 0 on failure
 *
 * @note `P12` are on the DOMAIN of the isogeny chain `Phi`.
 */
static int
theta_chain_eval_special_forward(theta_chain_t *Phi, theta_couple_point_t *P12, uint8_t numP)
{
    theta_point_t Ptheta;

    for (int j = 0; j < numP; ++j) {
        // evaluating the image points through the first step of the chain
        assert(ec_is_zero(&P12[j].P1) || ec_is_zero(&P12[j].P2));
        if (!gluing_eval_point_special_case(&Ptheta, &P12[j], &Phi->first_step, 0)) {
            debug_print("failure of the gluing evaluation during the HD chain");
            return 0;
        }

        // evaluating the image points through the generic steps of the chain
        for (int i = 1; i < Phi->length; ++i) {
            theta_isogeny_eval(&Ptheta, &Phi->steps[i], &Ptheta, (i != 1));
        }

        // splitting to couple points on the elliptic product
        if (!splitting_eval_point(&P12[j], &Phi->M_spl_E1, &Phi->M_spl_E2, &Ptheta)) {
            debug_print("failure of the splitting evaluation during the HD chain");
            return 0;
        }
    }

    return 1;
}

/**
 * @brief Evaluate at a list of points of the form (P1,0) or (0,P2).
 *
 * @param Phi a isogeny chain
 * @param P12 Input/Output: pointer to points to be pushed through the isogeny (in-place)
 * @param numP: length of the list of points given in P12 (can be zero)
 * @return 1 on success, 0 on failure
 *
 * @note `P12` are on the CODOMAIN of the isogeny chain `Phi`.
 */
static int
theta_chain_eval_special_dual(theta_chain_t *Phi, theta_couple_point_t *P12, uint8_t numP)
{
    ec_change_coord_matrix_t M_spl_E1_inv, M_spl_E2_inv;
    gluing_change_coord_matrix_t M_glu_inv;
    splitting_matrix_inv(&M_spl_E1_inv, &Phi->M_spl_E1);
    splitting_matrix_inv(&M_spl_E2_inv, &Phi->M_spl_E2);
    gluing_matrix_inv(&M_glu_inv, &Phi->first_step.M);

    theta_point_t Ptheta;

    for (int j = 0; j < numP; ++j) {
        // evaluating the image points through the dual of the splitting step
        assert(ec_is_zero(&P12[j].P1) || ec_is_zero(&P12[j].P2));
        splitting_eval_point_dual(&Ptheta, &M_spl_E1_inv, &M_spl_E2_inv, &P12[j]);

        // evaluating the image points through the dual of the generic steps of the chain
        if (!theta_isogeny_eval_dual(&Ptheta, &Phi->steps[Phi->length - 2], &Ptheta, true, true)) {
            debug_print("failure of the first dual evaluation during the HD chain");
            return 0;
        }
        for (int i = Phi->length - 3; i >= 0; --i) {
            theta_isogeny_eval_dual(&Ptheta, &Phi->steps[i], &Ptheta, false, false);
        }

        // evaluating the image points through the dual of the first step of the chain
        if (!gluing_eval_point_dual(&P12[j], &Ptheta, &Phi->first_step, &M_glu_inv, false)) {
            debug_print("failure of the dual gluing evaluation during the HD chain");
            return 0;
        }
    }

    return 1;
}

int
theta_chain_eval_special_forward_then_dual(theta_chain_t *Phi, theta_couple_point_t *P12, uint8_t numP)
{
    if (!theta_chain_eval_special_forward(Phi, P12, numP)) {
        debug_print("failure during the forward evaluation of the HD chain");
        return 0;
    }

    for (int j = 0; j < numP; ++j) {
        ec_point_init(&P12[j].P2);
    }

    if (!theta_chain_eval_special_dual(Phi, P12, numP)) {
        debug_print("failure during the dual evaluation of the HD chain");
        return 0;
    }

    return 1;
}

int
theta_chain_eval_special_dual_then_forward(theta_chain_t *Phi, theta_couple_point_t *P12, uint8_t numP)
{
    if (!theta_chain_eval_special_dual(Phi, P12, numP)) {
        debug_print("failure during the dual evaluation of the HD chain");
        return 0;
    }

    for (int j = 0; j < numP; ++j) {
        ec_point_init(&P12[j].P2);
    }

    if (!theta_chain_eval_special_forward(Phi, P12, numP)) {
        debug_print("failure during the forward evaluation of the HD chain");
        return 0;
    }

    return 1;
}
