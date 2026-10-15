#include "pokeb.h"
#include <hd.h>
#include <endomorphism_action.h>
#include <torsion_constants.h>
#include <klpt.h>
#include <quaternion.h>
#include <intbig.h>
#include <ec.h>
#include <rng.h>
#include <id2iso.h>
#include <biextension.h>
#include <fips202.h>
#include "b_isog.h"

int eval_dimtwo_isog_generic_dual(theta_chain_t *phi, ec_basis_t *evalPQ, const ec_basis_t *PQ)
{
    theta_couple_point_t input_points;
    theta_couple_point_t output_points;

    input_points.P1 = PQ->P;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case_superglue_vsq(&output_points, phi, &input_points);
    ec_set_zero(&output_points.P2);
    dual_theta_chain_eval_special_case(&output_points, phi, &output_points);
    copy_point(&evalPQ->P, &output_points.P1);

    input_points.P1 = PQ->Q;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case_superglue_vsq(&output_points, phi, &input_points);
    ec_set_zero(&output_points.P2);
    dual_theta_chain_eval_special_case(&output_points, phi, &output_points);
    copy_point(&evalPQ->Q, &output_points.P1);
    
    input_points.P1 = PQ->PmQ;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case_superglue_vsq(&output_points, phi, &input_points);
    ec_set_zero(&output_points.P2);
    dual_theta_chain_eval_special_case(&output_points, phi, &output_points);
    copy_point(&evalPQ->PmQ, &output_points.P1);
    
    return 1;
}

int eval_dimtwo_isog_single_generic_dual(theta_chain_t *phi, ec_point_t *evalP, const ec_point_t *P)
{
    theta_couple_point_t input_points;
    theta_couple_point_t output_points;

    input_points.P1 = *P;
    ec_set_zero(&input_points.P2);
    theta_chain_eval_special_case_superglue_vsq(&output_points, phi, &input_points);
    ec_set_zero(&output_points.P2);
    dual_theta_chain_eval_special_case(&output_points, phi, &output_points);
    copy_point(evalP, &output_points.P1);

    return 1;
}

int ibz_random_unit(ibz_t *q, const ibz_t *modulus, shake256ctx *state) {
    ibz_t gcd;
    ibz_init(&gcd);
    while (1) {
        if (state != NULL) {
            ibz_rand_interval_with_state(q, &ibz_const_zero, modulus, state);
        } else {
            ibz_rand_interval(q, &ibz_const_zero, modulus);
        }
        ibz_gcd(&gcd, q, modulus);
        if (ibz_is_one(&gcd)) {
            break;
        }
    }
    ibz_finalize(&gcd);
    return 1;
}

void compute_unreduced_matrix(ibz_mat_2x2_t *mat_unreduced, quat_alg_elem_t *theta) {
    ibz_t tmp, content;
    ibz_vec_4_t coeffs;
    ibz_init(&tmp); ibz_init(&content); ibz_vec_4_init(&coeffs);

    quat_alg_make_primitive(&coeffs, &content, theta, &MAXORD_O0, &QUATALG_PINFTY);

    ibz_set(&(*mat_unreduced)[0][0], 0); ibz_set(&(*mat_unreduced)[0][1], 0);
    ibz_set(&(*mat_unreduced)[1][0], 0); ibz_set(&(*mat_unreduced)[1][1], 0);

    for (unsigned i = 0; i < 2; ++i) {
        ibz_add(&(*mat_unreduced)[i][i], &(*mat_unreduced)[i][i], &coeffs[0]);
        for (unsigned j = 0; j < 2; ++j) {
            ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
            ibz_add(&(*mat_unreduced)[i][j], &(*mat_unreduced)[i][j], &tmp);
            
            ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
            ibz_add(&(*mat_unreduced)[i][j], &(*mat_unreduced)[i][j], &tmp);
            
            ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
            ibz_add(&(*mat_unreduced)[i][j], &(*mat_unreduced)[i][j], &tmp);
            
            ibz_mul(&(*mat_unreduced)[i][j], &(*mat_unreduced)[i][j], &content);
        }
    }
    ibz_finalize(&tmp); ibz_finalize(&content); ibz_vec_4_finalize(&coeffs);
}

void reduce_matrix_mod(ibz_mat_2x2_t *mat_out, const ibz_mat_2x2_t *mat_unreduced, const ibz_t *order) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_mod(&(*mat_out)[i][j], &(*mat_unreduced)[i][j], order);
        }
    }
}

void evaluate_linear_combination_generator(ec_point_t *out_gen, const ec_basis_t *bas, const ec_curve_t *E, const ibz_mat_2x2_t *mat, const ibz_t *order) {
    ibz_t X, Y, gcd_val, tmp_gcd;
    ibz_t k_ibz, tmp_mul;
    
    ibz_init(&X); ibz_init(&Y); ibz_init(&gcd_val); ibz_init(&tmp_gcd);
    ibz_init(&k_ibz); ibz_init(&tmp_mul);

    digit_t scalars[2][NWORDS_ORDER] = { 0 };
    long k = 0;

    while (1) {
        ibz_set(&k_ibz, k);
        ibz_mul(&tmp_mul, &(*mat)[0][1], &k_ibz);
        ibz_add(&X, &(*mat)[0][0], &tmp_mul);
        ibz_mod(&X, &X, order);

        ibz_mul(&tmp_mul, &(*mat)[1][1], &k_ibz);
        ibz_add(&Y, &(*mat)[1][0], &tmp_mul);
        ibz_mod(&Y, &Y, order);

        ibz_gcd(&gcd_val, &X, order);
        ibz_gcd(&tmp_gcd, &Y, &gcd_val);

        if (ibz_cmp(&tmp_gcd, &ibz_const_one) == 0) {
            break;
        }

        k++;
    }

    ibz_to_digits(scalars[0], &X);
    ibz_to_digits(scalars[1], &Y);

    ec_biscalar_mul(out_gen, E, scalars[0], scalars[1], bas);

    ibz_finalize(&X); ibz_finalize(&Y); ibz_finalize(&gcd_val); ibz_finalize(&tmp_gcd);
    ibz_finalize(&k_ibz); ibz_finalize(&tmp_mul);
}

int keygen(pokeb_sk_t *sk, pokeb_pk_t *pk) {
    ibz_t q, deg_isog, rhs, gcd, C, q_new;
    ibz_t A, q_bound;
    ibz_t alpha, beta, gamma1, gamma2, delta;
    
    ibz_init(&q); ibz_init(&deg_isog); ibz_init(&rhs); ibz_init(&gcd); ibz_init(&C); ibz_init(&q_new);
    ibz_init(&A); ibz_init(&q_bound);
    ibz_init(&alpha); ibz_init(&beta); ibz_init(&gamma1); ibz_init(&gamma2); ibz_init(&delta);
    
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_div_2exp(&q_bound, &TORSION_PLUS_2POWER, 4);
    ibz_mul(&C, &TORSION_C1, &TORSION_C2); 

    // Find candidate secret prime q
    for (int i = 0; i < 1000; i++) {
        ibz_rand_interval(&q, &ibz_const_zero, &q_bound);
        ibz_sub(&deg_isog, &A, &q); 
        
        if (ibz_divides(&q, &ibz_const_two) == 0) {
            ibz_gcd(&gcd, &q, &TORSION_C1);
            if (ibz_is_one(&gcd)) {
                ibz_gcd(&gcd, &q, &TORSION_C2);
                if (ibz_is_one(&gcd)) {
                    ibz_gcd(&gcd, &q, &TORSION_D);
                    if (ibz_is_one(&gcd)) {
                        ibz_gcd(&gcd, &deg_isog, &TORSION_C1);
                        if (ibz_is_one(&gcd)) {
                            ibz_gcd(&gcd, &deg_isog, &TORSION_C2);
                            if (ibz_is_one(&gcd)) {
                                ibz_gcd(&gcd, &deg_isog, &TORSION_D);
                                if (ibz_is_one(&gcd)) {
                                    break; 
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Compute quaternion norm
    ibz_mul(&rhs, &deg_isog, &q);
    ibz_mul(&rhs, &rhs, &C);

    // Sample and store private key
    ibz_random_unit(&alpha, &TORSION_PLUS_2POWER, NULL);
    ibz_random_unit(&beta,  &TORSION_PLUS_2POWER, NULL);
    ibz_random_unit(&gamma1, &TORSION_C1, NULL);
    ibz_random_unit(&gamma2, &TORSION_C2, NULL);
    ibz_random_unit(&delta, &TORSION_D, NULL);

    memset(&sk->deg, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->alpha, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->beta, 0, NWORDS_ORDER * RADIX / 8);
    memset(&sk->delta, 0, NWORDS_ORDER * RADIX / 8);

    ibz_to_digits(sk->deg, &q);
    ibz_to_digits(sk->alpha, &alpha);
    ibz_to_digits(sk->beta, &beta);
    ibz_to_digits(sk->delta, &delta);

    ec_curve_t curve = CURVE_E0;

    // Compute quaternion endomorphism matrix
    quat_alg_elem_t tau;
    quat_alg_elem_init(&tau);
    int rep_ok = represent_integer(&tau, &rhs, &QUATALG_PINFTY);
    assert(rep_ok);

    ibz_mat_2x2_t mat_unreduced, mat_c1, mat_c2, mat_2power;
    ibz_mat_2x2_init(&mat_unreduced); ibz_mat_2x2_init(&mat_c1); ibz_mat_2x2_init(&mat_c2);
    ibz_mat_2x2_init(&mat_2power);

    compute_unreduced_matrix(&mat_unreduced, &tau);
    reduce_matrix_mod(&mat_c1, &mat_unreduced, &TORSION_C1);
    reduce_matrix_mod(&mat_c2, &mat_unreduced, &TORSION_C2);
    reduce_matrix_mod(&mat_2power, &mat_unreduced, &TORSION_PLUS_2POWER);

    ec_basis_t eval_c1 = BASIS_C1;
    ec_point_t ker_C1;
    evaluate_linear_combination_generator(&ker_C1, &eval_c1, &curve, &mat_c1, &TORSION_C1);

    ec_basis_t eval_c2 = BASIS_C2;
    ec_point_t ker_C2;
    evaluate_linear_combination_generator(&ker_C2, &eval_c2, &curve, &mat_c2, &TORSION_C2);

    ec_basis_t E0_two = BASIS_EVEN;
    matrix_application_even_basis(&E0_two, &curve, &mat_2power, POWER_OF_2);

    ec_point_t E0_point_D = BASIS_D.P;

    // Evaluate 1D mixed isogeny
    ec_curve_t E1;
    ec_point_t mixed_pts[3];
    mixed_pts[0] = E0_two.P;   
    mixed_pts[1] = E0_two.Q;   
    mixed_pts[2] = E0_two.PmQ;

    evaluate_mixed_isogeny_pipeline(&E1, &curve, ker_C1, ker_C2, mixed_pts, 3);

    // Compute 2D Theta chain
    theta_couple_curve_t E01;
    theta_couple_point_t T1, T2, T1m2;
    theta_chain_t hd_isog;

    E01.E1 = curve; E01.E2 = E1;

    ibz_mul(&q_new, &q, &C);
    ibz_mod(&q_new, &q_new, &TORSION_PLUS_2POWER);

    digit_t q_new_dig[NWORDS_ORDER] = {0};
    ibz_to_digits(q_new_dig, &q_new);

    copy_point(&T1.P1, &BASIS_EVEN.P);
    copy_point(&T2.P1, &BASIS_EVEN.Q);
    copy_point(&T1m2.P1, &BASIS_EVEN.PmQ);

    ec_curve_normalize_A24(&curve);

    ec_mul_normalized_A24(&T1.P1, &curve, q_new_dig, &T1.P1, POWER_OF_2);
    ec_mul_normalized_A24(&T2.P1, &curve, q_new_dig, &T2.P1, POWER_OF_2);
    ec_mul_normalized_A24(&T1m2.P1, &curve, q_new_dig, &T1m2.P1, POWER_OF_2);

    T1.P2   = mixed_pts[0];
    T2.P2   = mixed_pts[1];
    T1m2.P2 = mixed_pts[2];

    theta_chain_comput_strategy_superglue(&hd_isog, POWER_OF_2 - 2, &E01, &T1, &T2, &T1m2, strategies[2], 1);
    pk->EA = E01.E2; 

    // Evaluate 2D Theta image points
    ec_basis_t evalPQC1;
    ec_basis_t evalPQC2;
    ec_point_t imX;

    eval_dimtwo_isog_generic_dual(&hd_isog, &evalPQC1, &BASIS_C1);
    eval_dimtwo_isog_generic_dual(&hd_isog, &evalPQC2, &BASIS_C2);
    eval_dimtwo_isog_single_generic_dual(&hd_isog, &imX, &E0_point_D);

    // Mask final public key
    digit_t alpha_dig[NWORDS_ORDER] = {0};
    digit_t beta_dig[NWORDS_ORDER] = {0};
    digit_t gamma1_dig[NWORDS_ORDER] = {0};
    digit_t gamma2_dig[NWORDS_ORDER] = {0};
    digit_t delta_dig[NWORDS_ORDER] = {0};
    
    ibz_to_digits(alpha_dig, &alpha);
    ibz_to_digits(beta_dig,  &beta);
    ibz_to_digits(gamma1_dig, &gamma1);
    ibz_to_digits(gamma2_dig, &gamma2);
    ibz_to_digits(delta_dig, &delta);

    ec_curve_normalize_A24(&pk->EA);

    // Mask 2^a torsion basis
    ec_mul_normalized_A24(&pk->PQ2.P, &pk->EA, alpha_dig, &mixed_pts[0], POWER_OF_2);
    ec_mul_normalized_A24(&pk->PQ2.Q, &pk->EA, beta_dig,  &mixed_pts[1], POWER_OF_2);
    
    ec_basis_t tmp_2;
    copy_point(&tmp_2.P, &mixed_pts[0]);
    copy_point(&tmp_2.Q, &mixed_pts[1]);
    copy_point(&tmp_2.PmQ, &mixed_pts[2]);
    xADD(&tmp_2.PmQ, &tmp_2.P, &tmp_2.Q, &tmp_2.PmQ);
    ec_biscalar_mul_bounded(&pk->PQ2.PmQ, &pk->EA, alpha_dig, beta_dig, &tmp_2, POWER_OF_2);

    // Mask C1 torsion basis
    ec_mul_normalized_A24(&pk->PQC1.P,   &pk->EA, gamma1_dig, &evalPQC1.P, TORSION_C1_BITLENGTH);
    ec_mul_normalized_A24(&pk->PQC1.Q,   &pk->EA, gamma1_dig, &evalPQC1.Q, TORSION_C1_BITLENGTH);
    ec_mul_normalized_A24(&pk->PQC1.PmQ, &pk->EA, gamma1_dig, &evalPQC1.PmQ, TORSION_C1_BITLENGTH);

    // Mask C2 torsion basis
    ec_mul_normalized_A24(&pk->PQC2.P,   &pk->EA, gamma2_dig, &evalPQC2.P, TORSION_C2_BITLENGTH);
    ec_mul_normalized_A24(&pk->PQC2.Q,   &pk->EA, gamma2_dig, &evalPQC2.Q, TORSION_C2_BITLENGTH);
    ec_mul_normalized_A24(&pk->PQC2.PmQ, &pk->EA, gamma2_dig, &evalPQC2.PmQ, TORSION_C2_BITLENGTH);

    // Mask point D
    ec_mul_normalized_A24(&pk->imX, &pk->EA, delta_dig, &imX, TORSION_D_BITLENGTH);

    // Clean up
    quat_alg_elem_finalize(&tau); 
    ibz_mat_2x2_finalize(&mat_unreduced); ibz_mat_2x2_finalize(&mat_c1); ibz_mat_2x2_finalize(&mat_c2);
    ibz_finalize(&q); ibz_finalize(&deg_isog); ibz_finalize(&rhs); ibz_finalize(&gcd); ibz_finalize(&C); ibz_finalize(&q_new);
    ibz_finalize(&A); ibz_finalize(&q_bound);
    ibz_finalize(&alpha); ibz_finalize(&beta); ibz_finalize(&gamma1); ibz_finalize(&gamma2); ibz_finalize(&delta);

    return 1;
}

int encrypt(pokeb_ct_t *ct, const pokeb_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len) {
    ibz_t beta1, beta2, omega, omega_inv, mask_D, A;
    ec_point_t pointT, eval_points_0[4], eval_points_A[4];
    ec_point_t EAB_imX; 
    shake256ctx state;

    digit_t beta1_scalar[NWORDS_ORDER] = {0};
    digit_t beta2_scalar[NWORDS_ORDER] = {0};
    digit_t omega_scalar[NWORDS_ORDER] = {0};
    digit_t omega_inv_scalar[NWORDS_ORDER] = {0};
    digit_t mask_D_scalar[NWORDS_ORDER] = {0};
    digit_t one_scalar[NWORDS_ORDER] = {0};
    one_scalar[0] = 1;

    ibz_init(&beta1); ibz_init(&beta2); ibz_init(&omega); ibz_init(&omega_inv);
    ibz_init(&mask_D); ibz_init(&A);

    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);

    // Sample random values
    if (seed != NULL) {
        shake256_absorb(&state, seed, seed_len);
        ibz_random_unit(&beta1, &TORSION_C1, &state); 
        ibz_random_unit(&beta2, &TORSION_C2, &state); 
        ibz_random_unit(&omega, &A, &state);          
        ibz_random_unit(&mask_D, &TORSION_D, &state); 
        shake256_ctx_release(&state);
    } else {
        ibz_random_unit(&beta1, &TORSION_C1, NULL);
        ibz_random_unit(&beta2, &TORSION_C2, NULL);
        ibz_random_unit(&omega, &A, NULL);
        ibz_random_unit(&mask_D, &TORSION_D, NULL);
    }

    ibz_invmod(&omega_inv, &omega, &A);

    ibz_to_digits(omega_scalar, &omega);
    ibz_to_digits(omega_inv_scalar, &omega_inv);
    ibz_to_digits(beta1_scalar, &beta1);
    ibz_to_digits(beta2_scalar, &beta2);
    ibz_to_digits(mask_D_scalar, &mask_D);

    // Path 1: E0 -> EB
    ec_point_t ker_C1_0, ker_C2_0;
    
    ec_biscalar_mul_bounded(&ker_C1_0, &CURVE_E0, one_scalar, beta1_scalar, &BASIS_C1, TORSION_C1_BYTES * 8);
    ec_biscalar_mul_bounded(&ker_C2_0, &CURVE_E0, one_scalar, beta2_scalar, &BASIS_C2, TORSION_C2_BYTES * 8);

    copy_point(&eval_points_0[0], &BASIS_EVEN.P);
    copy_point(&eval_points_0[1], &BASIS_EVEN.Q);
    copy_point(&eval_points_0[2], &BASIS_EVEN.PmQ);
    copy_point(&eval_points_0[3], &BASIS_D.P);

    evaluate_mixed_isogeny_pipeline(&ct->EB, &CURVE_E0, ker_C1_0, ker_C2_0, eval_points_0, 4);

    ec_curve_normalize_A24(&ct->EB);
    xMULv2_normalized(&ct->PQ2_B.P, &eval_points_0[0], omega_scalar, POWER_OF_2, &ct->EB.A24);
    xMULv2_normalized(&ct->PQ2_B.Q, &eval_points_0[1], omega_inv_scalar, POWER_OF_2, &ct->EB.A24);
    xADD(&pointT, &eval_points_0[0], &eval_points_0[1], &eval_points_0[2]); 
    xDBLMUL_bounded_normalized(&ct->PQ2_B.PmQ, &eval_points_0[0], omega_scalar, &eval_points_0[1], omega_inv_scalar, &pointT, &ct->EB, POWER_OF_2);

    xMULv2_normalized(&ct->imX_B, &eval_points_0[3], mask_D_scalar, TORSION_D_BITLENGTH, &ct->EB.A24);

    // Path 2: EA -> EAB
    ec_point_t ker_C1_A, ker_C2_A;
    
    ec_biscalar_mul_bounded(&ker_C1_A, &pk->EA, one_scalar, beta1_scalar, &pk->PQC1, TORSION_C1_BYTES * 8);
    ec_biscalar_mul_bounded(&ker_C2_A, &pk->EA, one_scalar, beta2_scalar, &pk->PQC2, TORSION_C2_BYTES * 8);

    copy_point(&eval_points_A[0], &pk->PQ2.P);
    copy_point(&eval_points_A[1], &pk->PQ2.Q);
    copy_point(&eval_points_A[2], &pk->PQ2.PmQ);
    copy_point(&eval_points_A[3], &pk->imX); 

    evaluate_mixed_isogeny_pipeline(&ct->EAB, &pk->EA, ker_C1_A, ker_C2_A, eval_points_A, 4);

    ec_curve_normalize_A24(&ct->EAB);
    xMULv2_normalized(&ct->PQ2_AB.P, &eval_points_A[0], omega_scalar, POWER_OF_2, &ct->EAB.A24);
    xMULv2_normalized(&ct->PQ2_AB.Q, &eval_points_A[1], omega_inv_scalar, POWER_OF_2, &ct->EAB.A24);
    xADD(&pointT, &eval_points_A[0], &eval_points_A[1], &eval_points_A[2]); 
    xDBLMUL_bounded_normalized(&ct->PQ2_AB.PmQ, &eval_points_A[0], omega_scalar, &eval_points_A[1], omega_inv_scalar, &pointT, &ct->EAB, POWER_OF_2);

    xMULv2_normalized(&EAB_imX, &eval_points_A[3], mask_D_scalar, TORSION_D_BITLENGTH, &ct->EAB.A24);

    // Hash and XOR to encapsulate shared secret
    unsigned char hash_input[2 * NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};

    ec_normalize_point(&EAB_imX);

    memcpy(hash_input, &EAB_imX.x.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + NWORDS_FIELD * RADIX / 8, &EAB_imX.x.im[0], NWORDS_FIELD * RADIX / 8);

    SHAKE256(hash_output, sizeof(hash_output), hash_input, sizeof(hash_input));

    memset(ct->ct, 0, sizeof(ct->ct));
    for (size_t i = 0; i < 32; i++) {
        if (i >= m_len) {
            ct->ct[i] = hash_output[i];
        } else {
            ct->ct[i] = m[i] ^ hash_output[i];
        }
    }

    ibz_finalize(&beta1); ibz_finalize(&beta2); ibz_finalize(&omega); ibz_finalize(&omega_inv);
    ibz_finalize(&mask_D); ibz_finalize(&A);

    return 1;
}

int decrypt(unsigned char *m, size_t *m_len, const pokeb_ct_t *ct, const pokeb_sk_t *sk) {
    unsigned char hash_input[2 * NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};
    digit_t T1_scalar[NWORDS_ORDER] = {0}, T2_scalar[NWORDS_ORDER] = {0};
    theta_chain_t hd_isog;
    theta_couple_curve_t EBAB;
    theta_couple_point_t T1, T2, T1m2;
    ec_basis_t PQ2_AB;
    ec_point_t imX_eval;
    ibz_t alpha_inv, beta_inv, deg, A;

    ibz_init(&alpha_inv);
    ibz_init(&beta_inv);
    ibz_init(&deg);
    ibz_init(&A);

    ibz_copy_digits(&deg, sk->deg, NWORDS_ORDER);
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_copy_digits(&alpha_inv, sk->alpha, NWORDS_ORDER);
    ibz_copy_digits(&beta_inv, sk->beta, NWORDS_ORDER);
    
    ibz_invmod(&alpha_inv, &alpha_inv, &TORSION_PLUS_2POWER);
    ibz_invmod(&beta_inv, &beta_inv, &TORSION_PLUS_2POWER);

    copy_point(&PQ2_AB.P, &ct->PQ2_AB.P);
    copy_point(&PQ2_AB.Q, &ct->PQ2_AB.Q);
    copy_point(&PQ2_AB.PmQ, &ct->PQ2_AB.PmQ);

    EBAB.E1 = ct->EB;
    EBAB.E2 = ct->EAB;

    ibz_t q_new, C_val;
    ibz_init(&q_new); 
    ibz_init(&C_val);

    // 1. C = TORSION_C1 * TORSION_C2
    ibz_mul(&C_val, &TORSION_C1, &TORSION_C2);
    
    // 2. q_new = (deg * C) mod 2^a
    ibz_mul(&q_new, &deg, &C_val);
    ibz_mod(&q_new, &q_new, &TORSION_PLUS_2POWER);

    // 3. Scalar multiplication with q_new
    digit_t q_new_dig[NWORDS_ORDER] = {0};
    ibz_to_digits(q_new_dig, &q_new);

    ec_curve_normalize_A24(&EBAB.E1);
    xMULv2_normalized(&T1.P1, &ct->PQ2_B.P, q_new_dig, POWER_OF_2, &EBAB.E1.A24);
    xMULv2_normalized(&T2.P1, &ct->PQ2_B.Q, q_new_dig, POWER_OF_2, &EBAB.E1.A24);
    xMULv2_normalized(&T1m2.P1, &ct->PQ2_B.PmQ, q_new_dig, POWER_OF_2, &EBAB.E1.A24);

    // Build 2D Theta isogeny core
    // ec_mul_ibz(&T1.P1, &EBAB.E1, &deg, &ct->PQ2_B.P);
    // ec_mul_ibz(&T2.P1, &EBAB.E1, &deg, &ct->PQ2_B.Q);
    // ec_mul_ibz(&T1m2.P1, &EBAB.E1, &deg, &ct->PQ2_B.PmQ);
    
    ibz_to_digits(T1_scalar, &alpha_inv);
    ibz_to_digits(T2_scalar, &beta_inv);
    
    ec_curve_normalize_A24(&EBAB.E2);
    xMULv2_normalized(&T1.P2, &ct->PQ2_AB.P, T1_scalar, POWER_OF_2, &EBAB.E2.A24);
    xMULv2_normalized(&T2.P2, &ct->PQ2_AB.Q, T2_scalar, POWER_OF_2, &EBAB.E2.A24);
    xADD(&PQ2_AB.PmQ, &PQ2_AB.P, &PQ2_AB.Q, &PQ2_AB.PmQ);
    xDBLMUL_bounded_normalized(&T1m2.P2, &ct->PQ2_AB.P, T1_scalar, &ct->PQ2_AB.Q, T2_scalar, &PQ2_AB.PmQ, &EBAB.E2, POWER_OF_2);

    if(!theta_chain_comput_strategy_superglue(&hd_isog, POWER_OF_2 - 2, &EBAB, &T1, &T2, &T1m2, strategies[2], 1)) {
        ibz_finalize(&alpha_inv); ibz_finalize(&beta_inv); ibz_finalize(&deg); ibz_finalize(&A); ibz_finalize(&q_new); ibz_finalize(&C_val);
        return 0;
    }
    copy_point(&imX_eval, &ct->imX_B);
    if (!eval_dimtwo_isog_single_generic_dual(&hd_isog, &imX_eval, &imX_eval)) {
        printf("Failed to evaluate the 2-dim isogeny\n");
        ibz_finalize(&alpha_inv); ibz_finalize(&beta_inv); ibz_finalize(&deg); ibz_finalize(&A); ibz_finalize(&q_new); ibz_finalize(&C_val);
        return 0;
    }

    xMULv2_normalized(&imX_eval, &imX_eval, sk->delta, TORSION_D_BITLENGTH, &EBAB.E2.A24);
    ec_normalize_point(&imX_eval);

    memcpy(hash_input, &imX_eval.x.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + NWORDS_FIELD * RADIX / 8, &imX_eval.x.im[0], NWORDS_FIELD * RADIX / 8);

    SHAKE256(hash_output, sizeof(hash_output), hash_input, sizeof(hash_input));

    memset(m, 0, sizeof(hash_output));
    for (size_t i = 0; i < sizeof(hash_output); i++) {
        m[i] = ct->ct[i] ^ hash_output[i];
    }
    *m_len = sizeof(hash_output);

    ibz_finalize(&alpha_inv);
    ibz_finalize(&beta_inv);
    ibz_finalize(&deg);
    ibz_finalize(&A);
    ibz_finalize(&q_new);
    ibz_finalize(&C_val);
    return 1;
}

#define POKEB_CT_ENCODED_BYTES (NWORDS_FIELD * RADIX * 18 / 8)

const unsigned char G_hash_str[9] = "encrypt_";
const size_t G_hash_str_len = 8;

int check_ct(const pokeb_ct_t *ct, const pokeb_pk_t *pk, const unsigned char *seed, const size_t seed_len) {
    ibz_t beta1, beta2, omega, omega_inv, mask_D, A;
    ec_curve_t EB, EAB;
    ec_basis_t PQ2_B, PQ2_AB;
    ec_point_t pointT, eval_points_0[4], eval_points_A[3], imX_B;
    bool validation_check;
    shake256ctx state;

    digit_t beta1_scalar[NWORDS_ORDER] = {0};
    digit_t beta2_scalar[NWORDS_ORDER] = {0};
    digit_t omega_scalar[NWORDS_ORDER] = {0};
    digit_t omega_inv_scalar[NWORDS_ORDER] = {0};
    digit_t mask_D_scalar[NWORDS_ORDER] = {0};
    digit_t one_scalar[NWORDS_ORDER] = {0};
    one_scalar[0] = 1;

    ibz_init(&beta1); ibz_init(&beta2); ibz_init(&omega); ibz_init(&omega_inv);
    ibz_init(&mask_D); ibz_init(&A);

    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);

    if (seed != NULL) {
        shake256_absorb(&state, seed, seed_len);
        ibz_random_unit(&beta1, &TORSION_C1, &state);
        ibz_random_unit(&beta2, &TORSION_C2, &state);
        ibz_random_unit(&omega, &A, &state);
        ibz_random_unit(&mask_D, &TORSION_D, &state);
        shake256_ctx_release(&state);
    } else {
        ibz_finalize(&beta1); ibz_finalize(&beta2); ibz_finalize(&omega); 
        ibz_finalize(&omega_inv); ibz_finalize(&mask_D); ibz_finalize(&A);
        return 0; 
    }

    ibz_invmod(&omega_inv, &omega, &A);
    ibz_to_digits(omega_scalar, &omega);
    ibz_to_digits(omega_inv_scalar, &omega_inv);
    ibz_to_digits(beta1_scalar, &beta1);
    ibz_to_digits(beta2_scalar, &beta2);
    ibz_to_digits(mask_D_scalar, &mask_D);

    // Path 1: E0 -> EB
    ec_point_t ker_C1_0, ker_C2_0;
    ec_biscalar_mul_bounded(&ker_C1_0, &CURVE_E0, one_scalar, beta1_scalar, &BASIS_C1, TORSION_C1_BYTES * 8);
    ec_biscalar_mul_bounded(&ker_C2_0, &CURVE_E0, one_scalar, beta2_scalar, &BASIS_C2, TORSION_C2_BYTES * 8);
    
    copy_point(&eval_points_0[0], &BASIS_EVEN.P);
    copy_point(&eval_points_0[1], &BASIS_EVEN.Q);
    copy_point(&eval_points_0[2], &BASIS_EVEN.PmQ);
    copy_point(&eval_points_0[3], &BASIS_D.P); 

    evaluate_mixed_isogeny_pipeline(&EB, &CURVE_E0, ker_C1_0, ker_C2_0, eval_points_0, 4);

    ec_curve_normalize_A24(&EB);
    xMULv2_normalized(&PQ2_B.P, &eval_points_0[0], omega_scalar, POWER_OF_2, &EB.A24);
    xMULv2_normalized(&PQ2_B.Q, &eval_points_0[1], omega_inv_scalar, POWER_OF_2, &EB.A24);
    xADD(&pointT, &eval_points_0[0], &eval_points_0[1], &eval_points_0[2]); 
    xDBLMUL_bounded_normalized(&PQ2_B.PmQ, &eval_points_0[0], omega_scalar, &eval_points_0[1], omega_inv_scalar, &pointT, &EB, POWER_OF_2);

    xMULv2_normalized(&imX_B, &eval_points_0[3], mask_D_scalar, TORSION_D_BITLENGTH, &EB.A24);

    validation_check = ec_is_equal(&ct->PQ2_B.P, &PQ2_B.P) & 
                       ec_is_equal(&ct->PQ2_B.Q, &PQ2_B.Q) & 
                       ec_is_equal(&ct->PQ2_B.PmQ, &PQ2_B.PmQ) &
                       ec_is_equal(&ct->imX_B, &imX_B);

    fp2_t j1, j2;
    ec_j_inv(&j1, &EB);
    ec_j_inv(&j2, &ct->EB);
    validation_check &= (fp2_is_equal(&j1, &j2) != 0);

    // Path 2: EA -> EAB
    ec_point_t ker_C1_A, ker_C2_A;
    ec_biscalar_mul_bounded(&ker_C1_A, &pk->EA, one_scalar, beta1_scalar, &pk->PQC1, TORSION_C1_BYTES * 8);
    ec_biscalar_mul_bounded(&ker_C2_A, &pk->EA, one_scalar, beta2_scalar, &pk->PQC2, TORSION_C2_BYTES * 8);

    copy_point(&eval_points_A[0], &pk->PQ2.P);
    copy_point(&eval_points_A[1], &pk->PQ2.Q);
    copy_point(&eval_points_A[2], &pk->PQ2.PmQ);

    evaluate_mixed_isogeny_pipeline(&EAB, &pk->EA, ker_C1_A, ker_C2_A, eval_points_A, 3);

    ec_curve_normalize_A24(&EAB);
    xMULv2_normalized(&PQ2_AB.P, &eval_points_A[0], omega_scalar, POWER_OF_2, &EAB.A24);
    xMULv2_normalized(&PQ2_AB.Q, &eval_points_A[1], omega_inv_scalar, POWER_OF_2, &EAB.A24);
    xADD(&pointT, &eval_points_A[0], &eval_points_A[1], &eval_points_A[2]); 
    xDBLMUL_bounded_normalized(&PQ2_AB.PmQ, &eval_points_A[0], omega_scalar, &eval_points_A[1], omega_inv_scalar, &pointT, &EAB, POWER_OF_2);

    validation_check &= ec_is_equal(&ct->PQ2_AB.P, &PQ2_AB.P) & 
                        ec_is_equal(&ct->PQ2_AB.Q, &PQ2_AB.Q) & 
                        ec_is_equal(&ct->PQ2_AB.PmQ, &PQ2_AB.PmQ);

    fp2_t j3, j4;
    ec_j_inv(&j3, &EAB);
    ec_j_inv(&j4, &ct->EAB);
    validation_check &= (fp2_is_equal(&j3, &j4) != 0);

    ibz_finalize(&A); ibz_finalize(&beta1); ibz_finalize(&beta2);
    ibz_finalize(&mask_D); ibz_finalize(&omega); ibz_finalize(&omega_inv);

    return validation_check;
}

int ct_encode(unsigned char *encoded_ct, pokeb_ct_t *ct) {
    ec_curve_normalize_A24(&ct->EB);
    ec_curve_normalize_A24(&ct->EAB);

    fp2_encode(encoded_ct, &ct->EB.A);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 2 / 8, &ct->PQ2_B.P.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 4 / 8, &ct->PQ2_B.Q.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 6 / 8, &ct->PQ2_B.PmQ.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 8 / 8, &ct->imX_B.x);
    
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 10 / 8, &ct->EAB.A);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 12 / 8, &ct->PQ2_AB.P.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 14 / 8, &ct->PQ2_AB.Q.x);
    fp2_encode(encoded_ct + NWORDS_FIELD * RADIX * 16 / 8, &ct->PQ2_AB.PmQ.x);

    return 1;
}

int encaps(unsigned char *key, pokeb_ct_t *ct, const pokeb_pk_t *pk) {
    unsigned char m[32];
    unsigned char tt[32 + G_hash_str_len];
    unsigned char gm[32];
    unsigned char encoded_ct[32 + POKEB_CT_ENCODED_BYTES];

    randombytes(m, 32);
    memcpy(tt, G_hash_str, G_hash_str_len);
    memcpy(tt + G_hash_str_len, m, 32);

    SHAKE256(gm, 32, tt, 32 + G_hash_str_len);
    
    encrypt(ct, pk, m, 32, gm, 32);
    
    memcpy(encoded_ct, m, 32);
    ct_encode(encoded_ct + 32, ct);
    
    SHAKE256(key, 32, encoded_ct, 32 + POKEB_CT_ENCODED_BYTES);
    
    return 1;
}

int decaps(unsigned char *key, pokeb_ct_t *ct, const pokeb_pk_t *pk, const pokeb_sk_t *sk, unsigned char *dummy_m) {
    unsigned char m[32];
    unsigned char tt[32 + G_hash_str_len];
    unsigned char gm[32];
    unsigned char ct_bytes[32 + POKEB_CT_ENCODED_BYTES];
    size_t m_len;

    ct_encode(ct_bytes + 32, ct);
    
    if(!decrypt(m, &m_len, ct, sk)) {
        memcpy(ct_bytes, dummy_m, 32);
    } else {
        memcpy(tt, G_hash_str, G_hash_str_len);
        memcpy(tt + G_hash_str_len, m, 32);
        
        SHAKE256(gm, 32, tt, 32 + G_hash_str_len);
        
        if (check_ct(ct, pk, gm, 32) != 1) {
            memcpy(ct_bytes, dummy_m, 32); 
        } else {
            memcpy(ct_bytes, m, 32);       
        }
    }
    
    SHAKE256(key, 32, ct_bytes, 32 + POKEB_CT_ENCODED_BYTES);

    return 1;
}