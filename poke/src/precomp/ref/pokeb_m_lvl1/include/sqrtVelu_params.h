#ifndef SQRTVELU_PARAMS_H
#define SQRTVELU_PARAMS_H

#include "ec_params.h"

// Execution order determination based on total xEVAL costs (including factor 3)
// Total xEVAL cost for C1: 22212 Fp-mul
// Total xEVAL cost for C2: 29978 Fp-mul
#define C1_FIRST 1 // 1 if C1 should be executed first (lower xEVAL cost), else 0

// Sqrt-Velu parameters (aligned with P_LEN and M_LEN, includes factor 3)
static const int sizeI[P_LEN + M_LEN] = 
	{ 1, 2, 2, 4, 11, 10, 12, 20, 1, 2, 3, 3, 3, 4, 4, 6, 6, 7, 10, 9, 14, 14, 14, 15, 14, 20 };
static const int sizeJ[P_LEN + M_LEN] = 
	{ 0, 0, 2, 2, 7, 9, 9, 13, 1, 1, 1, 2, 2, 2, 3, 4, 4, 6, 6, 8, 10, 12, 12, 16, 18, 17 };
static const int sizeK[P_LEN + M_LEN] = 
	{ 1, 3, 0, 7, 1, 3, 3, 4, 0, 1, 3, 2, 3, 5, 6, 3, 5, 2, 5, 9, 4, 0, 18, 11, 2, 6 };

#define sI_max 20
#define sJ_max 18
#define sK_max 191

#define ceil_log_sI_max 5
#define ceil_log_sJ_max 5

#define MAX_STRATEGY_LEAVES 100

// After stripping 3, the start and end indices (inclusive) of C1 and C2 in the TORSION_ODD_PRIMES array.
#define IDX_OF_3 0

#define C1_START_IDX 1
#define C1_END_IDX (P_LEN - 1)
#define C2_START_IDX P_LEN
#define C2_END_IDX (P_LEN + M_LEN - 1)

// Velu strategies for large odd primes (sorted from largest to smallest, excluding 3)
#define STRATEGY_C1_LEN 15
static const int STRATEGY_C1[STRATEGY_C1_LEN] = 
	{ 7, 4, 2, 1, 1, 1, 2, 1, 1, 6, 5, 4, 3, 2, 1 };

#define STRATEGY_C2_LEN 20
static const int STRATEGY_C2[STRATEGY_C2_LEN] = 
	{ 8, 5, 2, 2, 2, 1, 1, 1, 1, 3, 1, 2, 1, 7, 6, 5, 4, 3, 2, 1 };

#endif
