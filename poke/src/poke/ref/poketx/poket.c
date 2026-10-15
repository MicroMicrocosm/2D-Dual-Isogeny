#include "poket.h"
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

void evaluate_generator(ec_point_t *out_gen, const ec_basis_t *bas, const ec_curve_t *E, const ibz_mat_2x2_t *mat, const ibz_t *order) {
    ibz_t X, Y, gcd_val, tmp_gcd;
    ibz_t k_ibz, tmp_mul;
    
    ibz_init(&X); ibz_init(&Y); ibz_init(&gcd_val); ibz_init(&tmp_gcd);
    ibz_init(&k_ibz); ibz_init(&tmp_mul);

    digit_t scalars[2][NWORDS_ORDER] = { 0 };

    ibz_copy(&X, &(*mat)[0][0]);
    ibz_copy(&Y, &(*mat)[1][0]);
    ibz_gcd(&gcd_val, &X, order);
    ibz_gcd(&tmp_gcd, &Y, &gcd_val);

    if (ibz_cmp(&tmp_gcd, &ibz_const_one) != 0) {
        ibz_copy(&X, &(*mat)[0][1]);
        ibz_copy(&Y, &(*mat)[1][1]);
    }

    ibz_to_digits(scalars[0], &X);
    ibz_to_digits(scalars[1], &Y);

    ec_biscalar_mul_bounded(out_gen, E, scalars[0], scalars[1], bas, ibz_bitsize(order));

    ibz_finalize(&X); ibz_finalize(&Y); ibz_finalize(&gcd_val); ibz_finalize(&tmp_gcd);
    ibz_finalize(&k_ibz); ibz_finalize(&tmp_mul);
}


int keygen(poket_sk_t *sk, poket_pk_t *pk) {
    ibz_t q, deg_isog, rhs, gcd, q_new;
    ibz_t A, q_bound;
    ibz_t alpha, beta, gamma, delta;
    
    ibz_init(&q); ibz_init(&deg_isog); ibz_init(&rhs); ibz_init(&gcd); ibz_init(&q_new);
    ibz_init(&A); ibz_init(&q_bound);
    ibz_init(&alpha); ibz_init(&beta); ibz_init(&gamma); ibz_init(&delta);
    
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    ibz_div_2exp(&q_bound, &TORSION_PLUS_2POWER, 4);

    // ============================================================
    // Step 1: Find a candidate secret prime q satisfying coprimality conditions
    // ============================================================
    for (int i = 0; i < 1000; i++) {
        ibz_rand_interval(&q, &ibz_const_zero, &q_bound);
        ibz_sub(&deg_isog, &A, &q); 
        
        if (ibz_divides(&q, &ibz_const_two) == 0) {
            ibz_gcd(&gcd, &q, &ibz_const_three);
            if (ibz_is_one(&gcd)) {
                ibz_gcd(&gcd, &q, &TORSION_D);
                if (ibz_is_one(&gcd)) {
                    ibz_gcd(&gcd, &deg_isog, &ibz_const_three);
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

    // Compute the quaternion norm: rhs = q * (2^a - q) * 3^b
    ibz_mul(&rhs, &deg_isog, &q);
    ibz_mul(&rhs, &rhs, &TORSION_MINUS_3POWER);

    // ============================================================
    // Step 2: Sample and archive the private key
    // ============================================================
    ibz_random_unit(&alpha, &TORSION_PLUS_2POWER, NULL);
    ibz_random_unit(&beta,  &TORSION_PLUS_2POWER, NULL);
    ibz_random_unit(&gamma, &TORSION_MINUS_3POWER, NULL);
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

    // ============================================================
    // Step 3: Compute the quaternion endomorphism large matrix (using tau action)
    // ============================================================
    quat_alg_elem_t tau;
    quat_alg_elem_init(&tau);
    int rep_ok = represent_integer(&tau, &rhs, &QUATALG_PINFTY);
    assert(rep_ok);

    ibz_mat_2x2_t mat_unreduced, mat_3power, mat_2power;
    ibz_mat_2x2_init(&mat_unreduced); ibz_mat_2x2_init(&mat_3power); ibz_mat_2x2_init(&mat_2power);

    compute_unreduced_matrix(&mat_unreduced, &tau);
    reduce_matrix_mod(&mat_3power, &mat_unreduced, &TORSION_MINUS_3POWER);
    reduce_matrix_mod(&mat_2power, &mat_unreduced, &TORSION_PLUS_2POWER);

    ec_basis_t eval_three = BASIS_THREE;
    ec_point_t ker_3power;
    evaluate_generator(&ker_3power, &eval_three, &curve, &mat_3power, &TORSION_MINUS_3POWER);

    ec_basis_t E0_two = BASIS_EVEN;
    matrix_application_even_basis(&E0_two, &curve, &mat_2power, POWER_OF_2);

    ec_point_t E0_point_D = BASIS_D.P;

    // ============================================================
    // Step 4: Evaluate 1D isogeny pipeline
    // ============================================================
    ec_isog_odd_t isog;
    isog.curve = curve;
    isog.degree[0] = POWER_OF_3;
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isog.degree[i] = 0;
    }
    isog.ker_minus = ker_3power;
    ec_set_zero(&isog.ker_plus);

    ec_curve_t E1;
    ec_point_t mixed_pts[3];
    mixed_pts[0] = E0_two.P;
    mixed_pts[1] = E0_two.Q;
    mixed_pts[2] = E0_two.PmQ;

    ec_eval_three_twist(&E1, &isog, mixed_pts, 3);

    // ============================================================
    // Step 5: Setup and compute the 2D Theta chain
    // ============================================================
    theta_couple_curve_t E01;
    theta_couple_point_t T1, T2, T1m2;
    theta_chain_t hd_isog;

    E01.E1 = curve; E01.E2 = E1;

    ibz_mul(&q_new, &q, &TORSION_MINUS_3POWER);
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

    // ============================================================
    // Step 6: Evaluate 2D Theta image points (3^b and single point D)
    // ============================================================
    ec_basis_t evalPQ3;
    ec_point_t imX;

    eval_dimtwo_isog_generic_dual(&hd_isog, &evalPQ3, &BASIS_THREE);
    eval_dimtwo_isog_single_generic_dual(&hd_isog, &imX, &E0_point_D);

    // ============================================================
    // Step 7: Apply private key mask to the final public key (strictly using pure ec_mul)
    // ============================================================
    digit_t alpha_dig[NWORDS_ORDER] = {0};
    digit_t beta_dig[NWORDS_ORDER] = {0};
    digit_t gamma_dig[NWORDS_ORDER] = {0};
    digit_t delta_dig[NWORDS_ORDER] = {0};
    
    ibz_to_digits(alpha_dig, &alpha);
    ibz_to_digits(beta_dig,  &beta);
    ibz_to_digits(gamma_dig, &gamma);
    ibz_to_digits(delta_dig, &delta);

    ec_curve_normalize_A24(&pk->EA);

    // --- Mask the 2^a torsion basis (PQ2) ---
    ec_mul_normalized_A24(&pk->PQ2.P, &pk->EA, alpha_dig, &mixed_pts[0], POWER_OF_2);
    ec_mul_normalized_A24(&pk->PQ2.Q, &pk->EA, beta_dig,  &mixed_pts[1], POWER_OF_2);
    
    ec_basis_t tmp_2;
    copy_point(&tmp_2.P, &mixed_pts[0]);
    copy_point(&tmp_2.Q, &mixed_pts[1]);
    copy_point(&tmp_2.PmQ, &mixed_pts[2]);
    xADD(&tmp_2.PmQ, &tmp_2.P, &tmp_2.Q, &tmp_2.PmQ);
    ec_biscalar_mul_bounded(&pk->PQ2.PmQ, &pk->EA, alpha_dig, beta_dig, &tmp_2, POWER_OF_2);

    // --- Mask the 3^b torsion basis (PQ3) (using separate scalar multiplications) ---
    int three_bits = TORSION_3POWER_BYTES * 8;
    ec_mul_normalized_A24(&pk->PQ3.P,   &pk->EA, gamma_dig, &evalPQ3.P, three_bits);
    ec_mul_normalized_A24(&pk->PQ3.Q,   &pk->EA, gamma_dig, &evalPQ3.Q, three_bits);
    ec_mul_normalized_A24(&pk->PQ3.PmQ, &pk->EA, gamma_dig, &evalPQ3.PmQ, three_bits);

    // --- Mask the single point D twist image ---
    ec_mul_normalized_A24(&pk->imX, &pk->EA, delta_dig, &imX, TORSION_D_BITLENGTH);

    // ============================================================
    // Step 8: Clean up resources
    // ============================================================
    quat_alg_elem_finalize(&tau); 
    ibz_mat_2x2_finalize(&mat_unreduced); ibz_mat_2x2_finalize(&mat_3power); ibz_mat_2x2_finalize(&mat_2power);
    ibz_finalize(&q); ibz_finalize(&deg_isog); ibz_finalize(&rhs); ibz_finalize(&gcd); ibz_finalize(&q_new);
    ibz_finalize(&A); ibz_finalize(&q_bound);
    ibz_finalize(&alpha); ibz_finalize(&beta); ibz_finalize(&gamma); ibz_finalize(&delta);

    return 1;
}

int encrypt(poket_ct_t *ct, const poket_pk_t *pk, const unsigned char *m, const size_t m_len, const unsigned char *seed, const size_t seed_len) {
    ibz_t beta, omega, omega_inv, mask_D, A;
    ec_isog_odd_t isogB, isogB_prime;
    ec_point_t pointT, eval_points_0[4], eval_points_A[4];
    ec_point_t EAB_imX; 
    shake256ctx state;

    digit_t beta_scalar[NWORDS_ORDER] = {0};
    digit_t omega_scalar[NWORDS_ORDER] = {0};
    digit_t omega_inv_scalar[NWORDS_ORDER] = {0};
    digit_t mask_D_scalar[NWORDS_ORDER] = {0};
    digit_t one_scalar[NWORDS_ORDER] = {0};
    one_scalar[0] = 1;

    ibz_init(&beta); ibz_init(&omega); ibz_init(&omega_inv);
    ibz_init(&mask_D); ibz_init(&A);

    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);

    // ==========================================
    // Random number sampling (seed mechanism)
    // ==========================================
    if (seed != NULL) {
        shake256_absorb(&state, seed, seed_len);
        ibz_random_unit(&beta, &TORSION_MINUS_3POWER, &state); // 3^b side isogeny kernel scalar
        ibz_random_unit(&omega, &A, &state);                  // 2^a torsion image point mask
        ibz_random_unit(&mask_D, &TORSION_D, &state);         // D torsion uniform mask
        shake256_ctx_release(&state);
    } else {
        ibz_random_unit(&beta, &TORSION_MINUS_3POWER, NULL);
        ibz_random_unit(&omega, &A, NULL);
        ibz_random_unit(&mask_D, &TORSION_D, NULL);
    }

    // ibz_invmod(&omega_inv, &omega, &A);
    ibz_invmod(&omega_inv, &omega, &TORSION_PLUS_2POWER);

    ibz_to_digits(omega_scalar, &omega);
    ibz_to_digits(omega_inv_scalar, &omega_inv);
    ibz_to_digits(beta_scalar, &beta);
    ibz_to_digits(mask_D_scalar, &mask_D);

    // ==========================================
    // Path 1: E0 -> EB (\psi)
    // ==========================================
    isogB.curve = CURVE_E0;
    isogB.degree[0] = POWER_OF_3;
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isogB.degree[i] = 0;
    }
    
    // Construct isogeny kernel on E0
    ec_set_zero(&isogB.ker_plus);
    ec_biscalar_mul_bounded(&isogB.ker_minus, &CURVE_E0, one_scalar, beta_scalar, &BASIS_THREE, TORSION_3POWER_BYTES * 8);

    // Transport the 2^a basis and the D torsion point P
    copy_point(&eval_points_0[0], &BASIS_EVEN.P);
    copy_point(&eval_points_0[1], &BASIS_EVEN.Q);
    copy_point(&eval_points_0[2], &BASIS_EVEN.PmQ);
    copy_point(&eval_points_0[3], &BASIS_D.P);

    // Evaluate through the 3-isogeny twist pipeline
    ec_eval_three_twist(&ct->EB, &isogB, eval_points_0, 4);

    // Apply scalar masking to the image points on E_B
    ec_curve_normalize_A24(&ct->EB);
    xMULv2_normalized(&ct->PQ2_B.P, &eval_points_0[0], omega_scalar, POWER_OF_2, &ct->EB.A24);
    xMULv2_normalized(&ct->PQ2_B.Q, &eval_points_0[1], omega_inv_scalar, POWER_OF_2, &ct->EB.A24);
    xADD(&pointT, &eval_points_0[0], &eval_points_0[1], &eval_points_0[2]); 
    xDBLMUL_bounded_normalized(&ct->PQ2_B.PmQ, &eval_points_0[0], omega_scalar, &eval_points_0[1], omega_inv_scalar, &pointT, &ct->EB, POWER_OF_2);

    // Generate ciphertext element imX_B
    xMULv2_normalized(&ct->imX_B, &eval_points_0[3], mask_D_scalar, TORSION_D_BITLENGTH, &ct->EB.A24);

    // ==========================================
    // Path 2: EA -> EAB (\psi')
    // ==========================================
    isogB_prime.curve = pk->EA;
    isogB_prime.degree[0] = POWER_OF_3;
    for(int i = 1; i < P_LEN + M_LEN; i++) {
        isogB_prime.degree[i] = 0;
    }
    
    // Construct isogeny kernel on EA (using the basis images published by Alice)
    ec_set_zero(&isogB_prime.ker_plus);
    ec_biscalar_mul_bounded(&isogB_prime.ker_minus, &pk->EA, one_scalar, beta_scalar, &pk->PQ3, TORSION_3POWER_BYTES * 8);

    // Transport the 2^a basis on EA and the D torsion image imX
    copy_point(&eval_points_A[0], &pk->PQ2.P);
    copy_point(&eval_points_A[1], &pk->PQ2.Q);
    copy_point(&eval_points_A[2], &pk->PQ2.PmQ);
    copy_point(&eval_points_A[3], &pk->imX); 

    // Evaluate through the 3-isogeny twist pipeline again
    ec_eval_three_twist(&ct->EAB, &isogB_prime, eval_points_A, 4);

    // Apply scalar masking to the image points on E_AB
    ec_curve_normalize_A24(&ct->EAB);
    xMULv2_normalized(&ct->PQ2_AB.P, &eval_points_A[0], omega_scalar, POWER_OF_2, &ct->EAB.A24);
    xMULv2_normalized(&ct->PQ2_AB.Q, &eval_points_A[1], omega_inv_scalar, POWER_OF_2, &ct->EAB.A24);
    xADD(&pointT, &eval_points_A[0], &eval_points_A[1], &eval_points_A[2]); 
    xDBLMUL_bounded_normalized(&ct->PQ2_AB.PmQ, &eval_points_A[0], omega_scalar, &eval_points_A[1], omega_inv_scalar, &pointT, &ct->EAB, POWER_OF_2);

    // Generate the final shared secret point using mask_D
    xMULv2_normalized(&EAB_imX, &eval_points_A[3], mask_D_scalar, TORSION_D_BITLENGTH, &ct->EAB.A24);

    // ==========================================
    // Extract shared secret via KDF and XOR
    // ==========================================
    unsigned char hash_input[2 * NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};

    ec_normalize_point(&EAB_imX);

    // Copy the affine x-coordinate real and imaginary parts into the Hash
    memcpy(hash_input, &EAB_imX.x.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + NWORDS_FIELD * RADIX / 8, &EAB_imX.x.im[0], NWORDS_FIELD * RADIX / 8);

    SHAKE256(hash_output, sizeof(hash_output), hash_input, sizeof(hash_input));

    // Final XOR to encapsulate the message
    memset(ct->ct, 0, sizeof(ct->ct));
    for (size_t i = 0; i < 32; i++) {
        if (i >= m_len) {
            ct->ct[i] = hash_output[i];
        } else {
            ct->ct[i] = m[i] ^ hash_output[i];
        }
    }

    // Memory cleanup
    ibz_finalize(&beta); ibz_finalize(&omega); ibz_finalize(&omega_inv);
    ibz_finalize(&mask_D); ibz_finalize(&A);

    return 1;
}

int decrypt(unsigned char *m, size_t *m_len, const poket_ct_t *ct, const poket_sk_t *sk) {
    unsigned char hash_input[2 * NWORDS_FIELD * RADIX / 8] = {0};
    unsigned char hash_output[32] = {0};
    digit_t T1_scalar[NWORDS_ORDER] = {0}, T2_scalar[NWORDS_ORDER] = {0};
    
    theta_chain_t hd_isog;
    theta_couple_curve_t EBAB;
    theta_couple_point_t T1, T2, T1m2;
    ec_basis_t PQ2_AB;
    ec_point_t imX_eval;
    ibz_t alpha_inv, beta_inv, deg, A;
    ibz_t q_new;

    ibz_init(&alpha_inv);
    ibz_init(&beta_inv);
    ibz_init(&deg);
    ibz_init(&A);
    ibz_init(&q_new);

    // ==========================================
    // 1. Extract private key and compute scalars
    // ==========================================
    ibz_copy_digits(&deg, sk->deg, NWORDS_ORDER);
    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);
    
    ibz_copy_digits(&alpha_inv, sk->alpha, NWORDS_ORDER);
    ibz_copy_digits(&beta_inv, sk->beta, NWORDS_ORDER);
    
    // Calculate the inverses for eliminating the alpha and beta masks on E_AB
    ibz_invmod(&alpha_inv, &alpha_inv, &TORSION_PLUS_2POWER);
    ibz_invmod(&beta_inv, &beta_inv, &TORSION_PLUS_2POWER);

    copy_point(&PQ2_AB.P, &ct->PQ2_AB.P);
    copy_point(&PQ2_AB.Q, &ct->PQ2_AB.Q);
    copy_point(&PQ2_AB.PmQ, &ct->PQ2_AB.PmQ);

    EBAB.E1 = ct->EB;
    EBAB.E2 = ct->EAB;

    // ==========================================
    // 2. Compute q_new = deg * 3^b mod 2^a
    // ==========================================
    ibz_mul(&q_new, &deg, &TORSION_MINUS_3POWER);
    ibz_mod(&q_new, &q_new, &TORSION_PLUS_2POWER);

    digit_t q_new_dig[NWORDS_ORDER] = {0};
    ibz_to_digits(q_new_dig, &q_new);

    // ==========================================
    // 3. Construct the kernels for the 2D isogeny chain
    // ==========================================
    // Left side (E_B): eliminate the cofactor generated by the isogeny branch
    ec_curve_normalize_A24(&EBAB.E1);
    xMULv2_normalized(&T1.P1, &ct->PQ2_B.P, q_new_dig, POWER_OF_2, &EBAB.E1.A24);
    xMULv2_normalized(&T2.P1, &ct->PQ2_B.Q, q_new_dig, POWER_OF_2, &EBAB.E1.A24);
    xMULv2_normalized(&T1m2.P1, &ct->PQ2_B.PmQ, q_new_dig, POWER_OF_2, &EBAB.E1.A24);
    
    // Right side (E_AB): eliminate the \alpha and \beta masks
    ibz_to_digits(T1_scalar, &alpha_inv);
    ibz_to_digits(T2_scalar, &beta_inv);
    
    ec_curve_normalize_A24(&EBAB.E2);
    xMULv2_normalized(&T1.P2, &ct->PQ2_AB.P, T1_scalar, POWER_OF_2, &EBAB.E2.A24);
    xMULv2_normalized(&T2.P2, &ct->PQ2_AB.Q, T2_scalar, POWER_OF_2, &EBAB.E2.A24);
    xADD(&PQ2_AB.PmQ, &PQ2_AB.P, &PQ2_AB.Q, &PQ2_AB.PmQ);
    xDBLMUL_bounded_normalized(&T1m2.P2, &ct->PQ2_AB.P, T1_scalar, &ct->PQ2_AB.Q, T2_scalar, &PQ2_AB.PmQ, &EBAB.E2, POWER_OF_2);

    // Compute the 2D isogeny strategy using superglue optimization
    if(!theta_chain_comput_strategy_superglue(&hd_isog, POWER_OF_2 - 2, &EBAB, &T1, &T2, &T1m2, strategies[2], 1)) {
        ibz_finalize(&alpha_inv); ibz_finalize(&beta_inv); ibz_finalize(&deg); ibz_finalize(&A); ibz_finalize(&q_new);
        return 0; // Isogeny chain computation failed
    }

    // ==========================================
    // 4. Evaluate the D torsion image point to E_AB
    // ==========================================
    copy_point(&imX_eval, &ct->imX_B);
    if (!eval_dimtwo_isog_single_generic_dual(&hd_isog, &imX_eval, &imX_eval)) {
        ibz_finalize(&alpha_inv); ibz_finalize(&beta_inv); ibz_finalize(&deg); ibz_finalize(&A); ibz_finalize(&q_new);
        return 0; // Point evaluation failed
    }

    // Generate the final shared secret point using private key delta
    xMULv2_normalized(&imX_eval, &imX_eval, sk->delta, TORSION_D_BITLENGTH, &EBAB.E2.A24);
    ec_normalize_point(&imX_eval);

    // ==========================================
    // 5. KDF and message recovery
    // ==========================================
    memcpy(hash_input, &imX_eval.x.re[0], NWORDS_FIELD * RADIX / 8);
    memcpy(hash_input + NWORDS_FIELD * RADIX / 8, &imX_eval.x.im[0], NWORDS_FIELD * RADIX / 8);

    SHAKE256(hash_output, sizeof(hash_output), hash_input, sizeof(hash_input));

    memset(m, 0, sizeof(hash_output));
    for (size_t i = 0; i < sizeof(hash_output); i++) {
        m[i] = ct->ct[i] ^ hash_output[i];
    }
    *m_len = sizeof(hash_output);

    // Resource cleanup
    ibz_finalize(&alpha_inv); ibz_finalize(&beta_inv);
    ibz_finalize(&deg); ibz_finalize(&A); ibz_finalize(&q_new);
    
    return 1;
}

// CT encodes 9 fp2_t elements: EB.A, PQ2_B(3), imX_B(1), EAB.A, PQ2_AB(3)
#define POKET_CT_ENCODED_BYTES (NWORDS_FIELD * RADIX * 18 / 8)

const unsigned char G_hash_str[9] = "encrypt_";
const size_t G_hash_str_len = 8;

int check_ct(const poket_ct_t *ct, const poket_pk_t *pk, const unsigned char *seed, const size_t seed_len) {
    ec_isog_odd_t isogB, isogB_prime;
    ibz_t beta, omega, omega_inv, mask_D, A;
    ec_curve_t EB, EAB;
    ec_basis_t PQ2_B, PQ2_AB;
    ec_point_t pointT, eval_points_0[4], eval_points_A[3];
    ec_point_t imX_B;
    bool validation_check;
    shake256ctx state;

    digit_t beta_scalar[NWORDS_ORDER] = {0};
    digit_t omega_scalar[NWORDS_ORDER] = {0};
    digit_t omega_inv_scalar[NWORDS_ORDER] = {0};
    digit_t mask_D_scalar[NWORDS_ORDER] = {0};
    digit_t one_scalar[NWORDS_ORDER] = {0};
    one_scalar[0] = 1;

    ibz_init(&beta); ibz_init(&omega); ibz_init(&omega_inv);
    ibz_init(&mask_D); ibz_init(&A);

    ibz_div_2exp(&A, &TORSION_PLUS_2POWER, 2);

    if (seed != NULL) {
        shake256_absorb(&state, seed, seed_len);
        ibz_random_unit(&beta, &TORSION_MINUS_3POWER, &state);
        ibz_random_unit(&omega, &A, &state);
        ibz_random_unit(&mask_D, &TORSION_D, &state);
        shake256_ctx_release(&state);
    } else {
        ibz_finalize(&beta); ibz_finalize(&omega); ibz_finalize(&omega_inv);
        ibz_finalize(&mask_D); ibz_finalize(&A);
        return 0; 
    }

    ibz_invmod(&omega_inv, &omega, &TORSION_PLUS_2POWER);
    
    ibz_to_digits(omega_scalar, &omega);
    ibz_to_digits(omega_inv_scalar, &omega_inv);
    ibz_to_digits(beta_scalar, &beta);
    ibz_to_digits(mask_D_scalar, &mask_D);

    // Path 1: Verify E0 -> EB
    isogB.curve = CURVE_E0;
    isogB.degree[0] = POWER_OF_3;
    for(int i = 1; i < P_LEN + M_LEN; i++) isogB.degree[i] = 0;
    
    ec_set_zero(&isogB.ker_plus);
    ec_biscalar_mul_bounded(&isogB.ker_minus, &CURVE_E0, one_scalar, beta_scalar, &BASIS_THREE, TORSION_3POWER_BYTES * 8);
    
    copy_point(&eval_points_0[0], &BASIS_EVEN.P);
    copy_point(&eval_points_0[1], &BASIS_EVEN.Q);
    copy_point(&eval_points_0[2], &BASIS_EVEN.PmQ);
    copy_point(&eval_points_0[3], &BASIS_D.P); 

    ec_eval_three_twist(&EB, &isogB, eval_points_0, 4);

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

    // Path 2: Verify EA -> EAB
    isogB_prime.curve = pk->EA;
    isogB_prime.degree[0] = POWER_OF_3;
    for(int i = 1; i < P_LEN + M_LEN; i++) isogB_prime.degree[i] = 0;
    
    ec_set_zero(&isogB_prime.ker_plus);
    ec_biscalar_mul_bounded(&isogB_prime.ker_minus, &pk->EA, one_scalar, beta_scalar, &pk->PQ3, TORSION_3POWER_BYTES * 8);

    // Optimized: Only 3 points are evaluated as imX_AB is not included in the ciphertext
    copy_point(&eval_points_A[0], &pk->PQ2.P);
    copy_point(&eval_points_A[1], &pk->PQ2.Q);
    copy_point(&eval_points_A[2], &pk->PQ2.PmQ);

    ec_eval_three_twist(&EAB, &isogB_prime, eval_points_A, 3);

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

    ibz_finalize(&beta); ibz_finalize(&omega); ibz_finalize(&omega_inv);
    ibz_finalize(&mask_D); ibz_finalize(&A);

    return validation_check;
}

int ct_encode(unsigned char *encoded_ct, poket_ct_t *ct) {
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

int encaps(unsigned char *key, poket_ct_t *ct, const poket_pk_t *pk) {
    unsigned char m[32];
    unsigned char tt[32 + G_hash_str_len];
    unsigned char gm[32];
    unsigned char encoded_ct[32 + POKET_CT_ENCODED_BYTES];

    randombytes(m, 32);
    memcpy(tt, G_hash_str, G_hash_str_len);
    memcpy(tt + G_hash_str_len, m, 32);  

    SHAKE256(gm, 32, tt, 32 + G_hash_str_len);
    
    encrypt(ct, pk, m, 32, gm, 32);
    
    memcpy(encoded_ct, m, 32);
    ct_encode(encoded_ct + 32, ct);
    
    SHAKE256(key, 32, encoded_ct, 32 + POKET_CT_ENCODED_BYTES);
    
    return 1;
}

int decaps(unsigned char *key, poket_ct_t *ct, const poket_pk_t *pk, const poket_sk_t *sk, unsigned char *dummy_m) {
    unsigned char m[32];
    unsigned char tt[32 + G_hash_str_len];
    unsigned char gm[32];
    unsigned char ct_bytes[32 + POKET_CT_ENCODED_BYTES];
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
    
    SHAKE256(key, 32, ct_bytes, 32 + POKET_CT_ENCODED_BYTES);

    return 1;
}