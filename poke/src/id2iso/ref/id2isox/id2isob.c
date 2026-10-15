#include <quaternion.h>
#include <ec.h>
#include <endomorphism_action.h>
#include <id2iso.h>
#include <inttypes.h>
#include <locale.h>
#include <bench.h>
#include <curve_extras.h>
#include <biextension.h>

static __inline__ uint64_t
rdtsc(void)
{
    return (uint64_t)cpucycles();
}

void
id2iso_long_two_isog_init(id2iso_long_two_isog_t *isog, const size_t length)
{
    isog->length = length;
    isog->chain = malloc(length * sizeof(*isog->chain));
}

void
id2iso_long_two_isog_finalize(id2iso_long_two_isog_t *isog)
{
    free(isog->chain);
}

void
id2iso_compressed_long_two_isog_init(id2iso_compressed_long_two_isog_t *zip, const size_t length)
{
    zip->length = length;
    zip->zip_chain = malloc(length * sizeof(*zip->zip_chain));
    for (size_t i = 0; i < length; ++i)
        ibz_init(&zip->zip_chain[i]);
}

void
id2iso_compressed_long_two_isog_finalize(id2iso_compressed_long_two_isog_t *zip)
{
    for (size_t i = 0; i < zip->length; ++i)
        ibz_finalize(&zip->zip_chain[i]);
    free(zip->zip_chain);
}

// untested
void
matrix_of_endomorphism_even(ibz_mat_2x2_t *mat, const quat_alg_elem_t *alpha)
{
    ibz_t tmp;
    ibz_init(&tmp);

    // construct the matrix of alpha on the 2^f-torsion
    ibz_vec_4_t coeffs;
    ibz_vec_4_init(&coeffs);
    from_1ijk_to_O0basis(&coeffs, alpha);

    ibz_mat_2x2_set(mat, 0, 0, 0, 0);

    for (unsigned i = 0; i < 2; ++i) {
        ibz_add(&((*mat)[i][i]), &((*mat)[i][i]), &coeffs[0]);
        ibz_mod(&((*mat)[i][i]), &((*mat)[i][i]), &TORSION_PLUS_2POWER);

        for (unsigned j = 0; j < 2; ++j) {
            ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
            ibz_add(&((*mat)[i][j]), &((*mat)[i][j]), &tmp);
            ibz_mod(&((*mat)[i][j]), &((*mat)[i][j]), &TORSION_PLUS_2POWER);

            ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
            ibz_add(&((*mat)[i][j]), &((*mat)[i][j]), &tmp);
            ibz_mod(&((*mat)[i][j]), &((*mat)[i][j]), &TORSION_PLUS_2POWER);

            ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
            ibz_add(&((*mat)[i][j]), &((*mat)[i][j]), &tmp);
            ibz_mod(&((*mat)[i][j]), &((*mat)[i][j]), &TORSION_PLUS_2POWER);
        }
    }

#ifndef NDEBUG
    // compate determinant and norm
    ibz_t det, norm;
    ibq_t norm_q;
    ibz_init(&det);
    ibz_init(&norm);
    ibq_init(&norm_q);

    ibz_mul(&det, &((*mat)[0][0]), &((*mat)[1][1]));
    ibz_mul(&tmp, &((*mat)[0][1]), &((*mat)[1][0]));
    ibz_sub(&det, &det, &tmp);
    ibz_mod(&det, &det, &TORSION_PLUS_2POWER);

    quat_alg_norm(&norm_q, alpha, &QUATALG_PINFTY);
    ibq_to_ibz(&norm, &norm_q);
    ibz_mod(&norm, &norm, &TORSION_PLUS_2POWER);

    assert(ibz_cmp(&det, &norm) == 0);
    ibz_finalize(&det);
    ibz_finalize(&norm);
    ibq_finalize(&norm_q);
#endif

    ibz_vec_4_finalize(&coeffs);

    ibz_finalize(&tmp);
}

void
id2iso_ideal_to_kernel_dlogs_even(ibz_vec_2_t *vec, const quat_left_ideal_t *lideal)
{
    ibz_t tmp;
    ibz_init(&tmp);

    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);

    // construct the matrix of the dual of alpha on the 2^f-torsion
    {
        quat_alg_elem_t alpha;
        quat_alg_elem_init(&alpha);

        int lideal_generator_ok;
        lideal_generator_ok = quat_lideal_generator(&alpha, lideal, &QUATALG_PINFTY, 0);
        assert(lideal_generator_ok);
        quat_alg_conj(&alpha, &alpha);

        ibz_vec_4_t coeffs;
        ibz_vec_4_init(&coeffs);
        from_1ijk_to_O0basis(&coeffs, &alpha);

        for (unsigned i = 0; i < 2; ++i) {
            ibz_add(&mat[i][i], &mat[i][i], &coeffs[0]);
            for (unsigned j = 0; j < 2; ++j) {
                ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
                ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
                ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
            }
        }

        ibz_vec_4_finalize(&coeffs);
        quat_alg_elem_finalize(&alpha);
    }

    // find the kernel of alpha modulo the norm of the ideal
    {
        ibz_t const *const norm = &lideal->norm;

        ibz_mod(&(*vec)[0], &mat[0][0], norm);
        ibz_mod(&(*vec)[1], &mat[1][0], norm);
        ibz_gcd(&tmp, &(*vec)[0], &(*vec)[1]);
        if (!(ibz_get(&tmp) & 1)) {
            ibz_mod(&(*vec)[0], &mat[0][1], norm);
            ibz_mod(&(*vec)[1], &mat[1][1], norm);
        }
#ifndef NDEBUG
        ibz_gcd(&tmp, &(*vec)[0], norm);
        ibz_gcd(&tmp, &(*vec)[1], &tmp);
        assert(!ibz_cmp(&tmp, &ibz_const_one));
#endif
    }

    ibz_mat_2x2_finalize(&mat);
    ibz_finalize(&tmp);
}

void
id2iso_ideal_to_isogeny_even(ec_isog_even_t *isog, const quat_left_ideal_t *lideal_input)
{
    // compute length
    isog->length = 0;
    ibz_t norm;
    ibz_init(&norm);
    ibz_copy(&norm, &lideal_input->norm);
    while (!ibz_is_one(&norm)) {
        assert(!ibz_is_zero(&norm) && !(ibz_get(&norm) & 1));
        ibz_div_2exp(&norm, &norm, 1);
        ++isog->length;
    }
    ibz_finalize(&norm);
    assert(isog->length <= TORSION_PLUS_EVEN_POWER);

    digit_t scalars[2][NWORDS_FIELD];
    {
        ibz_vec_2_t vec;
        ibz_vec_2_init(&vec);

        id2iso_ideal_to_kernel_dlogs_even(&vec, lideal_input);

        // multiply out unnecessary cofactor from 2^f-torsion basis
        for (size_t i = isog->length; i < TORSION_PLUS_EVEN_POWER; ++i) {
            ibz_add(&vec[0], &vec[0], &vec[0]);
            ibz_add(&vec[1], &vec[1], &vec[1]);
        }

        ibz_to_digit_array(scalars[0], &vec[0]);
        ibz_to_digit_array(scalars[1], &vec[1]);

        ibz_vec_2_finalize(&vec);
    }

    isog->curve = CURVE_E0;
    ec_biscalar_mul(&isog->kernel, &isog->curve, scalars[0], scalars[1], &BASIS_EVEN);
}

void
ec_mul_ibz(ec_point_t *res, const ec_curve_t *curve, const ibz_t *scalarP, const ec_point_t *P)
{

    digit_t scalars[NWORDS_FIELD];
    if (ibz_cmp(scalarP, &ibz_const_one) == 0) {
        copy_point(res, P);
    } else {
        ibz_to_digit_array(scalars, scalarP);
        ec_mul(res, curve, scalars, P);
    }
}

void
ec_biscalar_mul_ibz(ec_point_t *res,
                    const ec_curve_t *curve,
                    const ibz_t *scalarP,
                    const ibz_t *scalarQ,
                    const ec_basis_t *PQ,
                    int f)
{

    digit_t scalars[2][NWORDS_FIELD];

    if (ibz_cmp(scalarP, &ibz_const_zero) == 0) {
        ec_mul_ibz(res, curve, scalarQ, &PQ->Q);
    } else if (ibz_cmp(scalarQ, &ibz_const_zero) == 0) {
        ec_mul_ibz(res, curve, scalarP, &PQ->P);
    } else {
        ibz_to_digit_array(scalars[0], scalarP);
        ibz_to_digit_array(scalars[1], scalarQ);
        ec_biscalar_mul_bounded(res, curve, scalars[0], scalars[1], PQ, f);
    }
}

// helper function to apply a matrix to a basis of E[2^f]
// works in place
void
matrix_application_even_basis(ec_basis_t *bas, const ec_curve_t *E, ibz_mat_2x2_t *mat, int f)
{
    digit_t scalars[2][NWORDS_FIELD] = { 0 };

    ibz_t tmp, pow_two;
    ibz_init(&tmp);
    ibz_init(&pow_two);
    ibz_pow(&pow_two, &ibz_const_two, f);

    ec_basis_t tmp_bas;
    copy_point(&tmp_bas.P, &bas->P);
    copy_point(&tmp_bas.Q, &bas->Q);
    copy_point(&tmp_bas.PmQ, &bas->PmQ);

    // reduction mod 2f
    ibz_mod(&(*mat)[0][0], &(*mat)[0][0], &pow_two);
    ibz_mod(&(*mat)[0][1], &(*mat)[0][1], &pow_two);
    ibz_mod(&(*mat)[1][0], &(*mat)[1][0], &pow_two);
    ibz_mod(&(*mat)[1][1], &(*mat)[1][1], &pow_two);

    // first basis element
    ibz_to_digit_array(scalars[0], &(*mat)[0][0]);
    // ibz_set(&mat[0][1],0);
    ibz_to_digit_array(scalars[1], &(*mat)[1][0]);
    ec_biscalar_mul_bounded(&bas->P, E, scalars[0], scalars[1], &tmp_bas, f);
    ibz_to_digit_array(scalars[0], &(*mat)[0][1]);
    ibz_to_digit_array(scalars[1], &(*mat)[1][1]);
    ec_biscalar_mul_bounded(&bas->Q, E, scalars[0], scalars[1], &tmp_bas, f);

    ibz_sub(&tmp, &(*mat)[0][0], &(*mat)[0][1]);
    ibz_mod(&tmp, &tmp, &pow_two);
    ibz_to_digit_array(scalars[0], &tmp);
    ibz_sub(&tmp, &(*mat)[1][0], &(*mat)[1][1]);
    ibz_mod(&tmp, &tmp, &pow_two);
    ibz_to_digit_array(scalars[1], &tmp);
    ec_biscalar_mul_bounded(&bas->PmQ, E, scalars[0], scalars[1], &tmp_bas, f);

    ibz_finalize(&tmp);
    ibz_finalize(&pow_two);
}

// helper function to apply some endomorphism of E0 on the precomputed basis of E[2^f]
// works in place
void
endomorphism_application_even_basis(ec_basis_t *bas,
                                    const ec_curve_t *E,
                                    quat_alg_elem_t *theta,
                                    int f)
{
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_vec_4_t coeffs;
    ibz_vec_4_init(&coeffs);
    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);

    ibz_t content;
    ibz_init(&content);

    // // decomposing theta on the basis
    quat_alg_make_primitive(&coeffs, &content, theta, &MAXORD_O0, &QUATALG_PINFTY);
    assert(ibz_get(&content) % 2 == 1);

    ibz_set(&mat[0][0], 0);
    ibz_set(&mat[0][1], 0);
    ibz_set(&mat[1][0], 0);
    ibz_set(&mat[1][1], 0);

    // computing the matrix
    for (unsigned i = 0; i < 2; ++i) {
        ibz_add(&mat[i][i], &mat[i][i], &coeffs[0]);
        for (unsigned j = 0; j < 2; ++j) {
            ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&mat[i][j], &mat[i][j], &content);
            // ibz_mod(&mat[i][j],&mat[i][j],&twopow);
        }
    }

    // and now we apply it
    matrix_application_even_basis(bas, E, &mat, f);

    ibz_vec_4_finalize(&coeffs);
    ibz_mat_2x2_finalize(&mat);
    ibz_finalize(&content);
}

// helper function to apply some endomorphism of E0 on the precomputed basis of E0[2^f]
// works in place
void
endomorphism_application_even_jac_basis(jac_point_t *P,
                                        jac_point_t *Q,
                                        quat_alg_elem_t *theta,
                                        int f)
{
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_vec_4_t coeffs;
    ibz_vec_4_init(&coeffs);
    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);
    ibz_t twopow;
    ibz_init(&twopow);
    ibz_pow(&twopow, &ibz_const_two, f);

    ibz_t content;
    ibz_init(&content);

    digit_t scalars[2][NWORDS_ORDER] = { 0 };

    // // decomposing theta on the basis
    quat_alg_make_primitive(&coeffs, &content, theta, &MAXORD_O0, &QUATALG_PINFTY);
    assert(ibz_get(&content) % 2 == 1);

    ibz_set(&mat[0][0], 0);
    ibz_set(&mat[0][1], 0);
    ibz_set(&mat[1][0], 0);
    ibz_set(&mat[1][1], 0);

    // computing the matrix
    for (unsigned i = 0; i < 2; ++i) {
        ibz_add(&mat[i][i], &mat[i][i], &coeffs[0]);
        for (unsigned j = 0; j < 2; ++j) {
            ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mod(&mat[i][j], &mat[i][j], &twopow);
        }
    }

    // and now we apply it
    jac_point_t tmp1, tmp2, tmp3;
    copy_jac_point(&tmp1, P);
    copy_jac_point(&tmp2, Q);

    // first basis element
    ibz_to_digit_array(scalars[0], &mat[0][0]);
    ibz_to_digit_array(scalars[1], &mat[1][0]);
    DBLMUL_generic(P, &tmp1, scalars[0], &tmp2, scalars[1], &CURVE_E0, NWORDS_ORDER);

    ibz_to_digit_array(scalars[0], &mat[0][1]);
    ibz_to_digit_array(scalars[1], &mat[1][1]);
    DBLMUL_generic(Q, &tmp1, scalars[0], &tmp2, scalars[1], &CURVE_E0, NWORDS_ORDER);

    ibz_finalize(&tmp);
    ibz_vec_4_finalize(&coeffs);
    ibz_mat_2x2_finalize(&mat);
    ibz_finalize(&twopow);
    ibz_finalize(&content);
}

// helper function to apply some endomorphism of E0 on a single jacobian point P of an arbitrary order
// works in place
void endomorphism_application_single_point(ec_point_t *R, ec_curve_t *E, const quat_alg_elem_t *theta, const ibz_t *order) {
    jac_point_t Rt, Rtt, R_jac;
    fp2_t fp2_i;
    digit_t tmp_digit[NWORDS_ORDER] = {0};
    ibz_t a, b, c, d, t;
    ec_point_t tmp;

    fp2_set_zero(&fp2_i);
    fp_set_one(&fp2_i.im);
    lift_point(&R_jac, R, E);
    // jac_init(&Rt);
    ibz_init(&a); ibz_init(&b); ibz_init(&c); ibz_init(&d); ibz_init(&t);
    ibz_invmod(&t, &theta->denom, order);
    ibz_mul(&a, &theta->coord[0], &t); ibz_mod(&a, &a, order);
    ibz_mul(&b, &theta->coord[1], &t); ibz_mod(&b, &b, order);
    ibz_mul(&c, &theta->coord[2], &t); ibz_mod(&c, &c, order);
    ibz_mul(&d, &theta->coord[3], &t); ibz_mod(&d, &d, order);

    // quat_alg_elem_print(theta);

    ibz_to_digits(tmp_digit, &a);
    MUL_generic(&Rt, &R_jac, tmp_digit, (*order)->_mp_size, E);

    // (x, y) -> (-x, yi)
    memset(tmp_digit, 0, NWORDS_ORDER * RADIX / 8);
    copy_jac_point(&Rtt, &R_jac);
    ibz_to_digits(tmp_digit, &b);
    fp2_neg(&Rtt.x, &Rtt.x);
    fp2_mul(&Rtt.y, &Rtt.y, &fp2_i);
    MUL_generic(&Rtt, &Rtt, tmp_digit, (*order)->_mp_size, E);
    ADD(&Rt, &Rt, &Rtt, E);

    // (x, y) -> (x^p, y^p)
    memset(tmp_digit, 0, NWORDS_ORDER * RADIX / 8);
    copy_jac_point(&Rtt, &R_jac);
    ibz_to_digits(tmp_digit, &c);
    fp2_conj(&Rtt.x, &Rtt.x);
    fp2_conj(&Rtt.y, &Rtt.y);
    MUL_generic(&Rtt, &Rtt, tmp_digit, (*order)->_mp_size, E);
    ADD(&Rt, &Rt, &Rtt, E);

    // (x, y) -> (-x^p, y^p * i)
    memset(tmp_digit, 0, NWORDS_ORDER * RADIX / 8);
    copy_jac_point(&Rtt, &R_jac);
    ibz_to_digits(tmp_digit, &d);
    fp2_conj(&Rtt.x, &Rtt.x);
    fp2_neg(&Rtt.x, &Rtt.x);
    fp2_conj(&Rtt.y, &Rtt.y);
    fp2_mul(&Rtt.y, &Rtt.y, &fp2_i);
    MUL_generic(&Rtt, &Rtt, tmp_digit, (*order)->_mp_size, E);
    ADD(&Rt, &Rt, &Rtt, E);

    jac_to_xz(R, &Rt);

    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&t);
}
