#ifndef DIFFERENTIAL_ADDITION_CHAIN_H
#define DIFFERENTIAL_ADDITION_CHAIN_H

#include <stdint.h>
#include "ec_params.h"

// Differential addition chains for odd primes (aligned with P_LEN and M_LEN)
static const uint32_t ADDITION_CHAIN[P_LEN + M_LEN] = 
{ 2, 192, 592, 3712, 2608, 4944, 5696, 0, 0, 4, 0, 0, 354, 328, 1290, 1032, 1700, 6272, 3080, 4610 };

static const uint8_t ADDITION_CHAIN_LENGTH[P_LEN + M_LEN] = 
{ 2, 8, 10, 12, 13, 14, 14, 0, 1, 3, 3, 9, 10, 10, 11, 11, 12, 13, 13, 14 };

#endif
