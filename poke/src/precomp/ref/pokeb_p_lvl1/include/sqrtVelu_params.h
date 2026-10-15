#ifndef SQRTVELU_PARAMS_H
#define SQRTVELU_PARAMS_H

#include "ec_params.h"

// Execution order determination based on total xEVAL costs (including factor 3)
// Total xEVAL cost for C1: 20409 Fp-mul
// Total xEVAL cost for C2: 28771 Fp-mul
#define C1_FIRST 1 // 1 if C1 should be executed first (lower xEVAL cost), else 0

// Sqrt-Velu parameters (aligned with P_LEN and M_LEN, includes factor 3)
static const int sizeI[P_LEN + M_LEN] = 
	{ 2, 5, 8, 13, 15, 20, 26, 1, 1, 2, 2, 7, 7, 10, 10, 10, 11, 15, 22, 26 };
static const int sizeJ[P_LEN + M_LEN] = 
	{ 0, 4, 8, 9, 16, 18, 15, 0, 1, 1, 1, 8, 8, 7, 10, 11, 12, 16, 12, 18 };
static const int sizeK[P_LEN + M_LEN] = 
	{ 3, 8, 3, 5, 8, 3, 3, 1, 0, 1, 2, 4, 7, 0, 4, 13, 14, 5, 3, 8 };

#define sI_max 26
#define sJ_max 18
#define sK_max 191

#define ceil_log_sI_max 5
#define ceil_log_sJ_max 5

#define MAX_STRATEGY_LEAVES 100

// After stripping 3, the start and end indices (inclusive) of C1 and C2 in the TORSION_ODD_PRIMES array.
#define IDX_OF_3 P_LEN

#define C1_START_IDX 0
#define C1_END_IDX (P_LEN - 1)
#define C2_START_IDX (P_LEN + 1)
#define C2_END_IDX (P_LEN + M_LEN - 1)

// Velu strategies for large odd primes (sorted from largest to smallest, excluding 3)
#define STRATEGY_C1_LEN 6
static const int STRATEGY_C1[STRATEGY_C1_LEN] = 
	{ 4, 1, 1, 3, 2, 1 };

#define STRATEGY_C2_LEN 23
static const int STRATEGY_C2[STRATEGY_C2_LEN] = 
	{ 8, 4, 4, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 7, 6, 5, 4, 3, 2, 1 };

#endif
