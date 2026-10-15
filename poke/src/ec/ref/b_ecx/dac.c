#include "curve_extras.h"
#include <ec_params.h>
#include <assert.h>
#include <stdio.h>
#include "differential_addition_chain.h"
#include "ec.h"

void DAC_prime(ec_point_t *Q, const ec_point_t *P, const ec_point_t *A24, uint8_t const i)
{
    ec_point_t R[3], T;

    copy_point(&R[0], P);
    xDBL_A24(&R[1], &R[0], A24); 
    xADD(&R[2], &R[1], &R[0], &R[0]); 
    
    uint8_t j;
    uint32_t tmp = ADDITION_CHAIN[i];
    for(j = 0; j < ADDITION_CHAIN_LENGTH[i]; j++)
    {
        xADD(&T, &R[2], &R[(tmp & 0x1) ^ 0x1], &R[tmp & 0x1]);
        
        copy_point(&R[0], &R[(tmp & 0x1) ^ 0x1]);
        copy_point(&R[1], &R[2]);
        copy_point(&R[2], &T);
        
        tmp >>= 1;
    }
    copy_point(Q, &R[2]); 
}