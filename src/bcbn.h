/*
 * bcbn.h
 *
 */

#ifndef BCBN_H_
#define BCBN_H_

#include "dt_cbn.h"

/**
 * computes the likelihood of the data
 * @param eps
 * @param M
 * @param theta
 * @param D
 * @param N_u
 * @param Prob
 * @param cond_err_prob
 * @return likelihood
 */
long double compute_likelihood( double eps, model* M, double* theta, data* D, int N_u, double* Prob, double** cond_err_prob );


/**
 * unnorm P[M,theta,eps|D]
 * @param eps
 * @param M
 * @param theta
 * @param D
 * @param N_u
 * @return unnormalized posterior
 */
long double compute_unnomr_posterior( double eps, model*M, double* theta, data* D, int N_u );
void relocate_whole_theta( double* theta, double* theta_p, int n );
void relocate_theta_i( double* theta, double* theta_p, int n, int i );
double compute_theta_transition_prob( double* theta, double* theta_p, int n );
void valid_new_cover_relation_matrix( model *M, int** V );
int propose_new_cover_relation( model *M, double* tp, data* D, int N_u );
double get_tp_for_remove_cover_move( model* M );
int propose_delete_cover_relation( model *M, double* tp );
double get_tp_for_new_cover_move( model* M, data* D, int N_u );
void propose_event_exchange_move( model* M, double* tp, int *ti, int *tj );
int propose_reincarnation_move( model* M, data* D, int N_u );
void relocate_epsilon( double eps, double* epsilon_p );
double get_tp_epsilon_relocation( double epsilon_p );
void valid_new_transitive_closer_relation_matrix( model *M, int** V );
int propose_new_transitive_closure_relation( model *M, double* tp, data* D, int N_u );
void valid_delete_transitive_closure_relation_matrix( model *M, int** V );
int propose_delete_transitive_closure_relation( model *M, double* tp, data* D, int N_u );
double get_tp_for_delete_transitive_closure_relation_move( model* M );
double get_tp_for_new_transitive_closure_relation_move( model* M );
// void initialize_random(uint32_t new_seed);
// double random_double();

/**
 *
 * @param eps IN epsilon
 * @param M
 * @param theta
 * @param D
 * @param N_u
 * @param burn_in
 * @param number_samples
 */
void start_MH( double eps, model* M, double* theta, data* D, int N_u, int burn_in, int number_samples, int record_ith );

/**
 * Expectation over theta (i.e. P[M|D)
 * @param eps
 * @param M
 * @param theta
 * @param D
 * @param N_u
 * @param burn_in
 * @param number_samples
 * @param record_ith
 * @return
 */
long double start_Exp_theta_MH( double eps, model* M, double* theta, data* D, int N_u, int burn_in, int number_samples, int record_ith );
void start_nested_MH(double eps, model* M, double* theta, data* D, int N_u, int burn_in, int number_samples, int record_ith );
void run_MH_sampler( model* M, double *theta_in, double epsilon_in, data* D, int N_u, int number_samples, int thinout, double **theta_matrix_out, int ***edges_cube_out, double *epsilon_out, double *log_posterior_out );

/**
 *
 * @param theta_in array[n] of starting thetas
 * @param nevents number of events
 * @param epsilon_in starting epsilon
 * @param edges_in array[n*n] of starting edges (poset)
 * @param number_smaples number of desired returned samples
 * @param thinout number of samples to discard between recording one
 * @param data array[n*ncases] data for inference
 * @param ncases number of cases i.e. data rows
 * @param theta_out array[n*number_samples] thetas sampled
 * @param epsilon_out array[number_samples] epsilons sampled
 * @param edges_out array[n*n*k] edges sampled (poset)
 */

void sample_full_cbn(double *theta_in, int *nevents, double *epsilon_in, int * edges_in, int *number_samples, int *thinout, int *patdata, int *number_cases, double *theta_out, double *epsilon_out, int* edges_out, double* log_posterior_out );

int main(int argc, char **argv);

#endif /* BCBN_H_ */
