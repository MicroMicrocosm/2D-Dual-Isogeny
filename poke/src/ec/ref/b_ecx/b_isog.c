#include "ec.h"
#include "b_isog.h"
#include <assert.h>


void ec_eval_strategy_core(
    ec_curve_t *image,
    const ec_curve_t *start_curve,
    ec_point_t kernel,
    const int *strategy,
    int start_p_idx,   
    int end_p_idx,     
    ec_point_t *points,
    size_t num_points)
{
    
    int leaf_array[MAX_STRATEGY_LEAVES];
    int total_leaves = 0;
    for (int p = end_p_idx; p >= start_p_idx; p--) {
        for (int k = 0; k < TORSION_ODD_POWERS[p]; k++) {
            leaf_array[total_leaves++] = p;
        }
    }

    ec_point_t R, A24, B24;
    ec_point_t stack[MAX_STRATEGY_LEAVES];
    unsigned int stack_index[MAX_STRATEGY_LEAVES];
    int index = 0, nstack = 0, ii = 0;
    int m;

    copy_point(&R, &kernel);
    AC_to_A24(&A24, start_curve);

    for (int row = 1; row < total_leaves; row++) {
        while (index < total_leaves - row) {
            copy_point(&stack[nstack], &R);
            stack_index[nstack++] = index;
            m = strategy[ii++];
            
            for (int j = 0; j < m; j++) {
                int p_idx = leaf_array[index + j];
                DAC_prime(&R, &R, &A24, p_idx); 
            }
            index += m;
        }
        
        int current_p_idx = leaf_array[index];
        kps(current_p_idx, R, A24);
        xisog(&B24, current_p_idx, A24);

        for(int i = 0; i < nstack; i++) {
            xeval(&stack[i], current_p_idx, stack[i], A24);
        }
        
        for(int i = 0; i < num_points; i++) {
            xeval(&points[i], current_p_idx, points[i], A24);
        }

        copy_point(&A24, &B24);
        copy_point(&R, &stack[nstack - 1]);
        index = stack_index[nstack - 1];
        nstack--;
    }
    
    int current_p_idx = leaf_array[index];
    kps(current_p_idx, R, A24);
    xisog(&B24, current_p_idx, A24);

    for(int i = 0; i < num_points; i++) {
        xeval(&points[i], current_p_idx, points[i], A24);
    }

    A24_to_AC(image, &B24);
    image->is_A24_computed_and_normalized = 0;
}

void ec_eval_three_core_isolated(
    ec_curve_t *image, 
    const ec_curve_t *start_curve, 
    ec_point_t kernel, 
    int degree, 
    ec_point_t *points, 
    size_t num_points)
{
    ec_point_t R, A24, A3, K;
    ec_point_t stack[POWER_OF_3];
    unsigned int stack_index[POWER_OF_3];
    int index = 0, nstack = 0, ii = 0;
    int m;

    copy_point(&R, &kernel);

    AC_to_A24(&A24, start_curve);
    fp2_sub(&A3.z, &A24.x, &A24.z);
    fp2_copy(&A3.x, &A24.x);

    for (int row = 1; row < degree; row ++) {
        while (index < degree - row) {
            copy_point(&stack[nstack], &R);
            stack_index[nstack++] = index;
            m = STRATEGY3[ii++];
            for (int j = 0; j < m; j++) {
                xTPL(&R, &R, &A3);
            }
            index += m;
        }
        kps_3(R, &K);
        xisog_3(&A24, K);
        fp2_sub(&A3.z, &A24.x, &A24.z);
        fp2_copy(&A3.x, &A24.x);

        for(int i = 0; i < nstack; i++) xeval_3(&stack[i], stack[i], K);
        for(int i = 0; i < num_points; i++) xeval_3(&points[i], points[i], K);

        copy_point(&R, &stack[nstack - 1]);
        index = stack_index[nstack - 1];
        nstack--;
    }
    kps_3(R, &K);
    xisog_3(&A24, K);

    for(int i = 0; i < num_points; i++) xeval_3(&points[i], points[i], K);

    A24_to_AC(image, &A24);
    image->is_A24_computed_and_normalized = 0;
}

static void split_kernel_with_3(
    ec_point_t *K_3, 
    ec_point_t *K_other, 
    const ec_point_t *ker, 
    const ec_curve_t *curve, 
    int start_idx, 
    int end_idx) 
{
    ec_point_t A24, A3;
    copy_point(K_3, ker);
    copy_point(K_other, ker);
    
    AC_to_A24(&A24, curve);
    fp2_sub(&A3.z, &A24.x, &A24.z);
    fp2_copy(&A3.x, &A24.x);

    for (int p = end_idx; p >= start_idx; p--) {
        if (TORSION_ODD_POWERS[p] > 0) {
            for (size_t k = 0; k < TORSION_ODD_POWERS[p]; k++) {
                DAC_prime(K_3, K_3, &A24, p);
            }
        }
    }
    
    for (size_t k = 0; k < TORSION_ODD_POWERS[IDX_OF_3]; k++) {
        xTPL(K_other, K_other, &A3);
    }
}

void evaluate_mixed_isogeny_pipeline(
    ec_curve_t *out_curve,
    const ec_curve_t *start_curve,
    ec_point_t ker_C1,  
    ec_point_t ker_C2,  
    ec_point_t *points, 
    size_t num_points) 
{
    ec_curve_t current_curve = *start_curve;

#if C1_FIRST == 1

    #if THREE_on_MAIN == 1
        ec_point_t K_3, K_C1_other;
        split_kernel_with_3(&K_3, &K_C1_other, &ker_C1, &current_curve, C1_START_IDX, C1_END_IDX);

        size_t push_len_3 = num_points + 2; 
        ec_point_t push_pts_3[push_len_3];
        push_pts_3[0] = K_C1_other; 
        push_pts_3[1] = ker_C2;     
        for (size_t i = 0; i < num_points; i++) push_pts_3[2 + i] = points[i];

        ec_curve_t curve_after_3;
        ec_eval_three_core_isolated(&curve_after_3, &current_curve, K_3, TORSION_ODD_POWERS[IDX_OF_3], push_pts_3, push_len_3);
        
        K_C1_other = push_pts_3[0];
        ker_C2     = push_pts_3[1];

        size_t push_len_c1 = num_points + 1;
        ec_point_t push_pts_c1[push_len_c1];
        push_pts_c1[0] = ker_C2;
        for (size_t i = 0; i < num_points; i++) push_pts_c1[1 + i] = push_pts_3[2 + i];

        ec_curve_t curve_after_c1;
        ec_eval_strategy_core(&curve_after_c1, &curve_after_3, K_C1_other, STRATEGY_C1, C1_START_IDX, C1_END_IDX, push_pts_c1, push_len_c1);

        ker_C2 = push_pts_c1[0];
        ec_eval_strategy_core(out_curve, &curve_after_c1, ker_C2, STRATEGY_C2, C2_START_IDX, C2_END_IDX, push_pts_c1 + 1, num_points);
        for (size_t i = 0; i < num_points; i++) points[i] = push_pts_c1[1 + i];

    #else
        size_t push_len_c1 = num_points + 1;
        ec_point_t push_pts_c1[push_len_c1];
        push_pts_c1[0] = ker_C2;
        for (size_t i = 0; i < num_points; i++) push_pts_c1[1 + i] = points[i];

        ec_curve_t curve_after_c1;
        ec_eval_strategy_core(&curve_after_c1, &current_curve, ker_C1, STRATEGY_C1, C1_START_IDX, C1_END_IDX, push_pts_c1, push_len_c1);

        ker_C2 = push_pts_c1[0];
        ec_point_t K_3, K_C2_other;
        split_kernel_with_3(&K_3, &K_C2_other, &ker_C2, &curve_after_c1, C2_START_IDX, C2_END_IDX);

        size_t push_len_3 = num_points + 1;
        ec_point_t push_pts_3[push_len_3];
        push_pts_3[0] = K_C2_other;
        for (size_t i = 0; i < num_points; i++) push_pts_3[1 + i] = push_pts_c1[1 + i];

        ec_curve_t curve_after_3;
        ec_eval_three_core_isolated(&curve_after_3, &curve_after_c1, K_3, TORSION_ODD_POWERS[IDX_OF_3], push_pts_3, push_len_3);
        
        K_C2_other = push_pts_3[0];
        ec_eval_strategy_core(out_curve, &curve_after_3, K_C2_other, STRATEGY_C2, C2_START_IDX, C2_END_IDX, push_pts_3 + 1, num_points);
        for (size_t i = 0; i < num_points; i++) points[i] = push_pts_3[1 + i];
    #endif

#else 

    #if THREE_on_MAIN == 1
        size_t push_len_c2 = num_points + 1;
        ec_point_t push_pts_c2[push_len_c2];
        push_pts_c2[0] = ker_C1; 
        for (size_t i = 0; i < num_points; i++) push_pts_c2[1 + i] = points[i];

        ec_curve_t curve_after_c2;
        ec_eval_strategy_core(&curve_after_c2, &current_curve, ker_C2, STRATEGY_C2, C2_START_IDX, C2_END_IDX, push_pts_c2, push_len_c2);
        
        ker_C1 = push_pts_c2[0]; 

        ec_point_t K_3, K_C1_other;
        split_kernel_with_3(&K_3, &K_C1_other, &ker_C1, &curve_after_c2, C1_START_IDX, C1_END_IDX);

        size_t push_len_3 = num_points + 1;
        ec_point_t push_pts_3[push_len_3];
        push_pts_3[0] = K_C1_other; 
        for (size_t i = 0; i < num_points; i++) push_pts_3[1 + i] = push_pts_c2[1 + i];

        ec_curve_t curve_after_3;
        ec_eval_three_core_isolated(&curve_after_3, &curve_after_c2, K_3, TORSION_ODD_POWERS[IDX_OF_3], push_pts_3, push_len_3);
        
        K_C1_other = push_pts_3[0];
        ec_eval_strategy_core(out_curve, &curve_after_3, K_C1_other, STRATEGY_C1, C1_START_IDX, C1_END_IDX, push_pts_3 + 1, num_points);
        for (size_t i = 0; i < num_points; i++) points[i] = push_pts_3[1 + i];

    #else
        ec_point_t K_3, K_C2_other;
        split_kernel_with_3(&K_3, &K_C2_other, &ker_C2, &current_curve, C2_START_IDX, C2_END_IDX);

        size_t push_len_3 = num_points + 2; 
        ec_point_t push_pts_3[push_len_3];
        push_pts_3[0] = K_C2_other; 
        push_pts_3[1] = ker_C1;     
        for (size_t i = 0; i < num_points; i++) push_pts_3[2 + i] = points[i];

        ec_curve_t curve_after_3;
        ec_eval_three_core_isolated(&curve_after_3, &current_curve, K_3, TORSION_ODD_POWERS[IDX_OF_3], push_pts_3, push_len_3);
        
        K_C2_other = push_pts_3[0];
        ker_C1     = push_pts_3[1];

        size_t push_len_c2 = num_points + 1;
        ec_point_t push_pts_c2[push_len_c2];
        push_pts_c2[0] = ker_C1;
        for (size_t i = 0; i < num_points; i++) push_pts_c2[1 + i] = push_pts_3[2 + i];

        ec_curve_t curve_after_c2;
        ec_eval_strategy_core(&curve_after_c2, &curve_after_3, K_C2_other, STRATEGY_C2, C2_START_IDX, C2_END_IDX, push_pts_c2, push_len_c2);

        ker_C1 = push_pts_c2[0];
        ec_eval_strategy_core(out_curve, &curve_after_c2, ker_C1, STRATEGY_C1, C1_START_IDX, C1_END_IDX, push_pts_c2 + 1, num_points);
        for (size_t i = 0; i < num_points; i++) points[i] = push_pts_c2[1 + i];
    #endif

#endif
}