/** @file
 *
 * @authors Damien Robert, Max Duparc
 *
 * @brief the splitting isogeny header
 *
 * @note Modified by Sai-Yu Wang on 2026/09/13.
 * Add parameters `M_spl_E1` and `M_spl_E2` to the function `splitting_eval_point`
 * to store the change of coordinate matrices for E1 and E2.
 *
 * Add function declaration:
 * - splitting_tensor_product_to_couple_point
 *
 * Add functions:
 * - splitting_eval_point
 * - splitting_eval_point_dual
 */

#ifndef SPLITTING_ISOGENY_H
#define SPLITTING_ISOGENY_H

#include <sqisign_namespace.h>

#include <fp2.h>
#include <ec.h>
#include <hd.h>
#include <assert.h>

#include "theta_structure.h"

/**
 * @ingroup hd_module
 * @defgroup hd_split Splitting isogenies
 * @{
 */

/**
 * @brief Convert a tensor product of theta_point to a couple_point.
 * @param res The resulting couple point
 * @param P The input tensor product of theta points
 * @return 1 on success, 0 on failure
 */
bool splitting_tensor_product_to_couple_point(theta_couple_point_t *res, const theta_point_t *P);

/**
 * @brief Convert the codomain to a product of Montgomery elliptic curves
 * @param E12 The two curves
 * @param P12 The image of the points on both curves
 * @param M_spl_E1 The change of coordinate matrix for E1
 * @param M_spl_E2 The change of coordinate matrix for E2
 * @param dual_null_point The codomain null point in dual form.
 * @param pts The dim 2 theta points to convert
 * @param numP The number of points
 * @param index the mode we use the splitting in:
 *                 0: compute both E1 and E2, together with the associated points.
 *                 1: compute only E1, together with the associated points.
 *                 2: compute only E2, together with the associated points.
 *                -1: Same as 1, but also checking validity.
 *
 */
bool splitting_to_elliptic_product(theta_couple_curve_t *E12,
                                   theta_couple_point_t *P12,
                                   ec_change_coord_matrix_t *M_spl_E1,
                                   ec_change_coord_matrix_t *M_spl_E2,
                                   theta_point_t *dual_null_point,
                                   theta_point_t *pts,
                                   size_t numP,
                                   int index);

/**
 * @brief Evaluate a list of points through the splitting isogeny.
 * @param P12 The image of the points on both curves
 * @param M_spl_E1 The change of coordinate matrix for E1
 * @param M_spl_E2 The change of coordinate matrix for E2
 * @param Ptheta The dim2 theta points to evaluate
 */
bool splitting_eval_point(theta_couple_point_t *P12,
                          ec_change_coord_matrix_t *M_spl_E1,
                          ec_change_coord_matrix_t *M_spl_E2,
                          theta_point_t *Ptheta);

/**
 * @brief Evaluate a couple_point through the dual of the splitting isogeny.
 * @param Ptheta The image of the couple_point
 * @param M_spl_E1_inv The inverse of the change of coordinate matrix for E1
 * @param M_spl_E2_inv The inverse of the change of coordinate matrix for E2
 * @param P12 The couple_point to evaluate
 */
void splitting_eval_point_dual(theta_point_t *Ptheta,
                               ec_change_coord_matrix_t *M_spl_E1_inv,
                               ec_change_coord_matrix_t *M_spl_E2_inv,
                               theta_couple_point_t *P12);

/**
 * @}
 */

#endif
