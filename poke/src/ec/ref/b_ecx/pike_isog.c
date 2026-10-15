#include "ec.h"
#include "isog.h"
#include <assert.h>


void
ec_eval_three_twist(ec_curve_t *image,
              const ec_isog_odd_t *phi,
              ec_point_t *points,
              unsigned short length)
{

    ec_point_t R, A24, A3, K;
    ec_point_t stack[POWER_OF_3];
    unsigned int stack_index[POWER_OF_3];
    int index = 0, nstack = 0, ii = 0;
    digit_t m;

    copy_point(&R, &phi->ker_minus);

    AC_to_A24(&A24, &phi->curve);
    fp2_sub(&A3.z, &A24.x, &A24.z);
    fp2_copy(&A3.x, &A24.x);

    // Traverse Tree
    for (int row = 1; row < phi->degree[0]; row ++) {
        while (index < phi->degree[0] - row) {
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

        for(int i = 0; i < nstack; i++) {
            xeval_3(&stack[i], stack[i], K);
        }

        for(int i = 0; i < length; i++) {
            xeval_3(&points[i], points[i], K);
        }

        copy_point(&R, &stack[nstack - 1]);
        index = stack_index[nstack - 1];
        nstack--;
    }
    kps_3(R, &K);
    xisog_3(&A24, K);

    for(int i = 0; i < length; i++) {
        xeval_3(&points[i], points[i], K);
    }

    A24_to_AC(image, &A24);
    image->is_A24_computed_and_normalized = 0;
}
