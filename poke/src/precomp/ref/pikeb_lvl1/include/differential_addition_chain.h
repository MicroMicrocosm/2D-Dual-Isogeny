#ifndef DIFFERENTIAL_ADDITION_CHAIN_H
#define DIFFERENTIAL_ADDITION_CHAIN_H

#include <stdint.h>
#include "ec_params.h"

// Differential addition chains for odd primes (aligned with P_LEN and M_LEN)
static const uint32_t ADDITION_CHAIN[P_LEN + M_LEN] = 
{ 0, 2, 10, 32, 160, 1216, 3744, 2144, 0, 4, 4, 16, 8, 34, 44, 48, 96, 192, 280, 264, 2400, 208, 2058, 4288, 5192, 2912 };

static const uint8_t ADDITION_CHAIN_LENGTH[P_LEN + M_LEN] = 
{ 0, 2, 4, 6, 10, 11, 12, 13, 1, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 10, 12, 12, 12, 13, 13, 14 };

#endif
