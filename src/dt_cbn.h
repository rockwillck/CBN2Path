/*
 * dt_cbn.h
 *
 */

#ifndef DT_CBN_H_
#define DT_CBN_H_

#include "queue.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
// #include <omp.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>

#define DOUBLE_FORMAT "%.10g"

#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

int verbose;

gsl_rng *RNG;  // random number generator

int** GENOTYPE;  // GENOTYPE[i] is the integer i in binary (as int array)


typedef struct {
  int** P;  // cover relations of the event poset
  int n;  // number of events
  int* lin_ext;  // linear extension of the poset
  
  int* J_P;  // lattice of order ideals (genotypes)
  int m;  // lattice size
  
  double eps;  // error tolerance used for constructing P
  
  int* N_pa;  // number of parents for each genotype
  int** pa;  // list of parents
  int** pa_diff;  // index in which parent and child differ
  
  int* N_ch; // number of children in lattice
  int** ch_diff;  // index in which parent and child differ
} model;


typedef struct {
  int *g;  // genotype
  int is_compatible;  // compatible with model?
  int** Q;  // induced refinement of poset
  int* J_Q;  // corresponding sublattice
  int count;  // number of observations of this type
} data;

typedef struct {
  model* M; //model
  data* D; //data
  int* N_u; //No obs
} gsl_parameters;


int* bcbn_get_int_array(const int n);
unsigned int* bcbn_get_uint_array(const int n);
double* bcbn_get_double_array(const int n);
int** bcbn_get_int_matrix(const int m, const int n);
double** bcbn_get_double_matrix(const int m, const int n);
double*** bcbn_get_double_cube(const int m, const int n, const int l);
int*** get_int_cube(const int m, const int n, const int l);
void bcbn_print_int_array(int* x, int n);
void bcbn_print_int_matrix(int** X, int m, int n);
void bcbn_print_double_array(double* x, int n);
void bcbn_print_double_matrix(double** X, int m, int n);
void bcbn_write_poset(int k, char* filestem, int** P, int n, int b);
void bcbn_write_patterns(char* filestem, int** pat, int N, int n);
inline void genotype_of(int index, int* x, int n);
int bcbn_pow2(int k);
void bcbn_precompute_binary(const int n);
int bcbn_index_of(int* x, int n);
int bcbn_is_equal_int_matrix(int** A, int** B, int n);
void bcbn_boolean_matrix_sum(int** A, int** B, int** C, int n);
void bcbn_boolean_matrix_product(int** A, int** B, int** C, int n);
void free_double_matrix( double** A, int n);
void free_int_matrix( int** A, int n);
void free_int_cube( int*** A, int n, int m);
void bcbn_free_poset(model* M);
void bcbn_free_lattice(model* M);
void bcbn_free_lattice_children(model* M);
void bcbn_free_model(model* M);
void print_model(model* M);
void bcbn_print_data(data* D, int N_u, int n, int m);
void bcbn_free_data(data* D, int N_u, int n);
int** bcbn_read_patterns(char* filestem, int* N, int n);
void bcbn_read_poset(char* filestem, model* M);
void bcbn_print_genotype(int* x, int n);
int* bcbn_bfs_order_ideals(int** poset, const int len, int* count, int* lin_ext);
int bcbn_norm1(const int g_idx, const int n);
double bcbn_power(double m, int n);
unsigned bcbn_hamdist(unsigned x, unsigned y);
int bcbn_hamming_distance(int g_idx, int h_idx, int* diff_idx, int n);
void parents_dt(model* M);
void children_dt(model* M);
void bcbn_read_poset_dt(char* filestem, model* M);
data* bcbn_make_data_set(int** pat, int N, int n, int* N_u, int* pat_idx);
void bcbn_compatibility(data* D, int N_u, model* M);
/**
 * Pr(X|M,theta) for all X in J(P)
 * @param theta IN mutation probabilities; size n
 * @param M IN the model
 * @param Prob OUT the probabilities of all valid (according to the poset) genotypes; size m
 * @param theta_exit IN the exit probabilities for all elements in the genotype lattice; size: m
 */
void compute_all_cbn_prob( double* theta, model* M, double* Prob, double* theta_exit );
/*
 * IN double *theta : mutation probabilities; size: n
 * IN model *M : the model
 * OUT double * theta_exit: the exit probabilities for all elements in the genotype lattice; size: m
 */
void compute_theta_exit(double *theta, model* M, double* theta_exit);
// Compute the conditional probability Prob[Z|X,epsilon]
void compute_cond_error_prob(model* M, data* D, int N_u, double epsilon, double** condprob);
void clone_poset( model* M_src, model* M_dest );
void bcbn_transitive_closure(int** A, int**T, int n);
/**
 *
 * @param P transitive closure to cover relations
 * @param n
 * @param C changes
 * @return status
 */
int bcbn_reduce_to_cover_relations ( int** P, int n, int** C );
#endif /* DT_CBN_H_ */
