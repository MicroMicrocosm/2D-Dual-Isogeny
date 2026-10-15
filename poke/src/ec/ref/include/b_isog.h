#ifndef B_ISOG_H
#define B_ISOG_H

#include <ec.h>
#include "curve_extras.h"
#include "sqrtVelu_params.h"
#include "differential_addition_chain.h"
#include "poly.h"

extern int sI, sJ, sK;	

extern fp2_t I[sI_max][2],		
			EJ_0[sJ_max][3], EJ_1[sJ_max][3];	

extern ec_point_t J[sJ_max], K[sK_max];		
extern fp2_t XZJ4[sJ_max],		
    rtree_A[(1 << (ceil_log_sI_max+2)) - 1],		
    A0;			

extern poly ptree_hI[(1 << (ceil_log_sI_max+2)) - 1],		
     rtree_hI[(1 << (ceil_log_sI_max+2)) - 1],		
     ptree_EJ[(1 << (ceil_log_sJ_max+2)) - 1];		
     
extern fp2_t R0[2*sJ_max + 1];		

extern int deg_ptree_hI[(1 << (ceil_log_sI_max+2)) - 1],	
    deg_ptree_EJ[(1 << (ceil_log_sJ_max+2)) - 1];	

extern fp2_t leaves[sI_max];		



typedef struct
{
    ec_point_t K;
} ec_kps2_t;
typedef struct
{
    ec_point_t K[3];
} ec_kps4_t;

extern ec_point_t K[]; 

void kps_t(uint64_t const i, ec_point_t const P, ec_point_t const A); 
void kps_3(const ec_point_t P, ec_point_t *K);

void xisog_3(ec_point_t* B, const ec_point_t K);
void xeval_3(ec_point_t* R, const ec_point_t Q, const ec_point_t K);

void xisog_2(ec_kps2_t *kps, ec_point_t *B, const ec_point_t P); 
void xisog_2_singular(ec_kps2_t *kps, ec_point_t *B24, ec_point_t A24);

void xisog_4(ec_kps4_t *kps, ec_point_t *B, const ec_point_t P); 
void xisog_4_singular(ec_kps4_t *kps, ec_point_t *B24, const ec_point_t P, ec_point_t A24);

void xeval_2(ec_point_t *R, ec_point_t *const Q, const int lenQ, const ec_kps2_t *kps);
void xeval_2_singular(ec_point_t *R, const ec_point_t *Q, const int lenQ, const ec_kps2_t *kps);

void xeval_4(ec_point_t *R, const ec_point_t *Q, const int lenQ, const ec_kps4_t *kps);
void xeval_4_singular(ec_point_t *R,
                      const ec_point_t *Q,
                      const int lenQ,
                      const ec_point_t P,
                      const ec_kps4_t *kps);


static void ec_eval_even_strategy(ec_curve_t *image,
                                  ec_point_t *points,
                                  unsigned short points_len,
                                  ec_point_t *A24,
                                  const ec_point_t *kernel,
                                  const int isog_len);

void xisog_t(ec_point_t *B, uint64_t const i, ec_point_t const A); 
void xeval_t(ec_point_t *Q, uint64_t const i, ec_point_t const P); 

void kps_s(uint64_t const i, ec_point_t const P, ec_point_t const A);	
void xisog_s(ec_point_t* B, uint64_t const i, ec_point_t const A);	
void xeval_s(ec_point_t* Q, uint64_t const i, ec_point_t const P, ec_point_t const A);	


void kps_clear(int i);	

static inline void kps(uint64_t const i, ec_point_t const P, ec_point_t const A)	
{
	if (TORSION_ODD_PRIMES[i] <= gap)
		kps_t(i, P, A);
	else
		kps_s(i, P, A);
}

static inline void xisog(ec_point_t* B, uint64_t const i, ec_point_t const A)
{
	
	if (TORSION_ODD_PRIMES[i] <= gap)
		xisog_t(B, i, A);
	else
		xisog_s(B, i, A);
}

static inline void xeval(ec_point_t* Q, uint64_t const i, ec_point_t const P, ec_point_t const A)
{
	
	if (TORSION_ODD_PRIMES[i] <= gap)
		xeval_t(Q, i, P);
	else
		xeval_s(Q, i, P, A);
}

void ec_eval_strategy_core(
    ec_curve_t *image,
    const ec_curve_t *start_curve,
    ec_point_t kernel,
    const int *strategy,
    int start_p_idx,   
    int end_p_idx,     
    ec_point_t *points,
    size_t num_points);
void ec_eval_three_core_isolated(
    ec_curve_t *image, 
    const ec_curve_t *start_curve, 
    ec_point_t kernel, 
    int degree, 
    ec_point_t *points, 
    size_t num_points);
void evaluate_mixed_isogeny_pipeline(
    ec_curve_t *out_curve,
    const ec_curve_t *start_curve,
    ec_point_t ker_C1,  
    ec_point_t ker_C2,  
    ec_point_t *points, 
    size_t num_points) ;

#endif