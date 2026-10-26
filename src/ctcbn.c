#include <stdio.h>

#include "ct-cbn.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdint.h>

// Function to convert char* to SEXP
SEXP char_to_sexp(const char* input) {
    // Create a character vector of length 1
    SEXP result = PROTECT(allocVector(STRSXP, 1));

    // Set the value of the character vector
    SET_STRING_ELT(result, 0, mkChar(input));

    // Unprotect the SEXP object to manage R's garbage collection
    UNPROTECT(1);

    // Return the SEXP object
    return result;
}

// Caller is responsible for freeing the returned buffer.
char* read_file(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return NULL;
    }

    // Determine the file size
    if (fseek(file, 0, SEEK_END) != 0) {
        fclose(file);
        return NULL;
    }
    long file_size = ftell(file);
    if (file_size < 0) {
        fclose(file);
        return NULL;
    }
    rewind(file);

    // Allocate memory for the file content (+1 for the null-terminator)
    char *buffer = (char *)malloc(((size_t)file_size + 1) * sizeof(char));
    if (buffer == NULL) {
        perror("Memory allocation error");
        fclose(file);
        return NULL;
    }

    // Read the file into the buffer
    size_t read_size = fread(buffer, sizeof(char), (size_t)file_size, file);
    if (read_size != (size_t)file_size) {
        perror("Error reading file");
        free(buffer);
        fclose(file);
        return NULL;
    }

    // Null-terminate the string
    buffer[file_size] = '\0';

    // Close the file
    fclose(file);

    return buffer;
}

SEXP ctcbn_(SEXP ofs, SEXP fs1, SEXP fs2, SEXP mb, SEXP bs, SEXP rs, SEXP sr, SEXP epsilon, SEXP nd, SEXP emr)
{
  const char *ofilestemRaw = CHAR(STRING_ELT(ofs, 0));
  char *ofilestem = strdup(ofilestemRaw);
  const char *filestem1Raw = CHAR(STRING_ELT(fs1, 0));
  char *filestem1 = strdup(filestem1Raw);
  const char *filestem2Raw = CHAR(STRING_ELT(fs2, 0));
  char *filestem2 = strdup(filestem2Raw);
  if (ofilestem == NULL || filestem1 == NULL || filestem2 == NULL)
  {
    free(ofilestem);
    free(filestem1);
    free(filestem2);
    Rf_error("Error: Out of memory!");
  }

  // defaults:
  double eps = REAL(epsilon)[0];                         // e
  int R = INTEGER(emr)[0];                                // # of EM runs
  double S = REAL(sr)[0];                           // sampling rate \lambda_s
  // unsigned int seed = (unsigned) time(NULL);  // r, random seed
  unsigned int seed = (unsigned) INTEGER(rs)[0]; // r, random seed
  int verbose = 0;
  int N_draw = INTEGER(nd)[0]; // # of samples to draw
  int error_flag = 0;
  int f_flag = 1;
  int e_flag = 0;
  int GPS = 0;
  int mode = LEARN_PARAM;
  int B = INTEGER(bs)[0];  // bootstrap samples
  int BM = INTEGER(mb)[0]; // bootstrap mode
  int c = 0;
  char* output;

  // // no epsilon provided from R
  if (eps > 1.0) {
    e_flag = 0;
    eps = 0.0;
  } else {
    e_flag = 1;
  }

  // if ((outputFileObj = fopen(outputFile, "w")) == NULL) {
  //   fprintf(stderr, "ERROR: Could not create file, %s\n", outputFile);;
  // }

  if (B > 0)
    if ((eps < 0.0) || (N_draw > 0) || ((BM == LEARN_POSET) && (!e_flag)))
      error_flag++;

  if (error_flag)
  {
    free(ofilestem);
    free(filestem1);
    free(filestem2);
    return(char_to_sexp("Error: Bad parameter values!"));
    // fprintf(stderr, "Error: Bad parameter values!\n");
    // exit(1);
  }

  if (!f_flag)
  {
    // fprintf(stderr, "Error: No input file specified!\n");
    // exit(1);
  }

  if (B > 0)
  {
    mode = BM;
  }
  else
  {
    if (e_flag) {
      mode = LEARN_BOTH;
    }
    else {
      mode = LEARN_PARAM;
    }
  }

  pcg_srand(seed);
  RNG = gsl_rng_alloc(gsl_rng_taus); // global variable
  gsl_rng_set(RNG, seed);            // seed rng

  int i, j, k;

  model M = {0};
  read_poset(filestem1, &M);

  // precompute binary expansions
  precompute_binary(M.n+1);

  M.lin_ext = get_int_array(M.n);             // a linear extension of the poset
  double *lambda = get_double_array(M.n+1); // Exp rates
  lambda[0] = S;


  if (N_draw == 0) // learn model
  {
    int N, N_u;
    int **pat = read_patterns(filestem2, &N, M.n);
    int *pat_idx = get_int_array(N);
    data *D = make_data_set(pat, N, M.n, &N_u, pat_idx);
    for (k = 0; k < N; k++)
      free(pat[k]);
    free(pat);
    free(pat_idx);

    if (eps >= 0.0) // fixed epsilon
    {
      int b;
      eps = MIN(eps, 1.0);

      // single run:
      output = select_poset(0, eps, &M, lambda, D, N_u, R, mode, 1);
      // if (e_flag)
      write_poset(0, ofilestem, M.P, M.n, -1);
      write_lambda(ofilestem, lambda, M.n);

      // bootstrap runs:
      double *p_orig = get_double_array(N_u); // frequencies of original data
      for (k = 0; k < N_u; k++)
        p_orig[k] = (double) D[k].count;

      int **bootstrap_count = get_int_matrix(M.n + 1, M.n + 1);       // all relations
      int **bootstrap_cover_count = get_int_matrix(M.n + 1, M.n + 1); // cover relations only
      for (i = 0; i <= M.n; i++)
        for (j = 0; j <= M.n; j++)
          bootstrap_count[i][j] = bootstrap_cover_count[i][j] = 0;

      int **T = get_int_matrix(M.n + 1, M.n + 1);
      for (b = 1; b <= B; b++)
      {
        resample(D, p_orig, N_u);
        if (output != NULL && strcmp(output, "N") != 0)
          free(output);
        output = select_poset(b, eps, &M, lambda, D, N_u, R, mode, 1); // e_flag

        int_matrix_sum(bootstrap_cover_count, M.P, bootstrap_cover_count, M.n + 1);
        transitive_closure(M.P, T, M.n + 1);
        int_matrix_sum(bootstrap_count, T, bootstrap_count, M.n + 1);

        // if (e_flag)
        write_poset(0, ofilestem, M.P, M.n, b);
      }

      // if ((B > 0) && (BM == 0))
      // {
      //   fprintf(outputFileObj, "\nbootstrap counts: matrix entry (i,j) counts edge i-->j\n");
      //   fprintf(outputFileObj, "\nall poset relations =\n");
      //   print_int_matrix(bootstrap_count, M.n + 1, M.n + 1);
      //   fprintf(outputFileObj, "\ncover relations =\n");
      //   print_int_matrix(bootstrap_cover_count, M.n + 1, M.n + 1);
      // }

      for (i = 0; i <= M.n; i++)
      {
        free(T[i]);
        free(bootstrap_count[i]);
        free(bootstrap_cover_count[i]);
      }
      free(T);
      free(bootstrap_count);
      free(bootstrap_cover_count);
      free(p_orig);
    }
    else // sample epsilon
    {
      double epsilon;
      int k;

      FILE *summaryFile;

      char summaryFileName[255];
      snprintf(summaryFileName, sizeof(summaryFileName), "%s.summary", ofilestem);

      summaryFile = fopen(summaryFileName, "w");

      output = "multi";

      for (k = 0; k < N; k++)
      {
        epsilon = (double)k / (2.0 * (double)N);
        char *poset_out = select_poset(k, epsilon, &M, lambda, D, N_u, R, LEARN_BOTH, 1);
        if (summaryFile != NULL) {
          fprintf(summaryFile, "%s\n", poset_out);
        }
        if (poset_out != NULL && strcmp(poset_out, "N") != 0) {
          free(poset_out);
        }

        if (e_flag) {
          write_poset(k, ofilestem, M.P, M.n, -1);
        }
      }

      if (summaryFile != NULL)
        fclose(summaryFile);
    }
    free_data(D, N_u, M.n);
    if (M.J_P != NULL)
      free_lattice(&M);
    free_poset(&M);
  }
  else // sample from given model:
  {
    // construct model:
    M.J_P = bfs_order_ideals(M.P, M.n + 1, &(M.m), M.lin_ext);
    parents(&M);
    children(&M);

    int **pat_draw = get_int_matrix(N_draw, M.n + 1);
    double **t_draw = get_double_matrix(N_draw, M.n + 1);

    read_lambda(filestem2, lambda, M.n);

    draw_samples(M.P, lambda, M.lin_ext, M.n, pat_draw, t_draw, N_draw);
    write_patterns(ofilestem, pat_draw, N_draw, M.n);
    write_times(ofilestem, t_draw, N_draw, M.n);

    for (k = 0; k < N_draw; k++)
    {
      free(pat_draw[k]);
      free(t_draw[k]);
    }
    free(pat_draw);
    free(t_draw);

    int n_genotypes = pow2(M.n + 1);
    for (i = 0; i < n_genotypes; i++)
      free(GENOTYPE[i]);
    free(GENOTYPE);
    GENOTYPE = NULL;

    free_model(&M);
    free(lambda);
    free(ofilestem);
    free(filestem1);
    free(filestem2);
    gsl_rng_free(RNG);
    RNG = NULL;
    return char_to_sexp("haiii >_<");
  }

  int n_genotypes = pow2(M.n + 1);
  for (i = 0; i < n_genotypes; i++)
    free(GENOTYPE[i]);
  free(GENOTYPE);
  GENOTYPE = NULL;

  free(lambda);
  free(ofilestem);
  free(filestem1);
  free(filestem2);
  gsl_rng_free(RNG);
  RNG = NULL;

  SEXP res = char_to_sexp(output);
  if (output != NULL && strcmp(output, "N") != 0 && strcmp(output, "multi") != 0)
    free(output);

  return res;
}

SEXP hcbn_(SEXP ofs, SEXP fs1, SEXP fs2, SEXP fs3, SEXP s, SEXP temp, SEXP n, SEXP repsilon)
{
  const char *ofilestemRaw = CHAR(STRING_ELT(ofs, 0));
  char *ofilestem = strdup(ofilestemRaw);
  const char *filestem1Raw = CHAR(STRING_ELT(fs1, 0));
  char *filestem1 = strdup(filestem1Raw);
  const char *filestem2Raw = CHAR(STRING_ELT(fs2, 0));
  char *filestem2 = strdup(filestem2Raw);
  const char *filestem3Raw = CHAR(STRING_ELT(fs3, 0));
  char *filestem3 = strdup(filestem3Raw);
  if (ofilestem == NULL || filestem1 == NULL || filestem2 == NULL || filestem3 == NULL)
  {
    free(ofilestem);
    free(filestem1);
    free(filestem2);
    free(filestem3);
    Rf_error("Error: Out of memory!");
  }

  // defaults:
  double eps = REAL(repsilon)[0];  // e
  int R = 1;  // # of EM runs
  double S = 1.0;  // sampling rate \lambda_s
  // char *filestem = CHAR(STRING_ELT(fs, 0));
  unsigned int seed = (unsigned) time(NULL);  // r, random seed
  verbose = 0;
  int N_draw = 0;  // # of samples to draw
  int error_flag = 0;
  int f_flag = 0;
  int e_flag = 0;
  int gps_flag = 0;
  int l_flag = 0;
  int s_flag = INTEGER(s)[0];
  int t_flag = 1;
  // n = 0, then equivalent to no argument supplied
  int n_flag = 0;
  int m_flag = 0;
  int w_flag = 1;
  double T = REAL(temp)[0];
  int N_iter = 0;
  int only_falsepos = 0;
  int p_flag = 0;
  int mode = LEARN_PARAM;
  double epsilon = 0.0;
  int c = 0;
  char* rOutput = "N";

  // // no epsilon provided from R
  if (eps > 1.0) {
    e_flag = 0;
    eps = 0.0;
  } else {
    e_flag = 1;
  }

  int n_temp = INTEGER(n)[0];
  if (n_temp > 0) {
    n_flag = 1;
    N_iter = n_temp;
  }

  pcg_srand(seed);
  RNG = gsl_rng_alloc(gsl_rng_taus);  // global variable
  gsl_rng_set(RNG, seed);  // seed rng

  int i, k;

  model M = {0};
  read_poset(filestem1, &M);

  // precompute binary expansions
  precompute_binary(M.n+1);

  M.lin_ext = get_int_array(M.n);  // a linear extension of the poset
  double* lambda = get_double_array(M.n+1);  // Exp rates
  lambda[0] = S;

  double total_loglik = 0.0;

  {

    if (N_draw == 0)  // learn model
    {
      int N, N_u;
      int** pat = read_patterns(filestem2, &N, M.n);
      int* pat_idx = get_int_array(N);
      data* D = make_data_set(pat, N, M.n, &N_u, pat_idx);
      for (k=0; k<N; k++)
        free(pat[k]);
      free(pat);
      free(pat_idx);

      if (eps >= 0.0)  // fixed epsilon
      {
        eps = MIN(eps, 1.0);

        /*  single run: */
        char *sel_out = select_poset(0 /* poset value in output (k in select_poset) */, eps, &M, lambda, D, N_u, R, mode,1);
        if (sel_out != NULL && strcmp(sel_out, "N") != 0)
          free(sel_out);

        if(eps > 0){
          read_lambda(filestem3, lambda, M.n);
        }

        /* Compute variables */
        double* Prob = get_double_array(M.m);
        double** condprob = get_double_matrix(M.m, N_u);
        double* lambda_exit = get_double_array(M.m);
        compute_lambda_exit(lambda, &M, lambda_exit);
        int* lattice_index = get_int_array(pow2(M.n+1));
        for (i=0; i < M.m; i++)
          lattice_index[M.J_P[i]]=i;
        compute_all_prob(lambda, &M, lambda_exit, Prob, lattice_index);

        /* Compute prob of all observations for ct-cbn;
         needed to be done here, befor params are re-estimated. */
        if (p_flag){
          double* ProbY_ctcbn = compute_ProbY_ctcbn(&M, Prob);
          write_double_array(ofilestem, ".prY_ctcbn", ProbY_ctcbn, pow2(M.n));
          free(ProbY_ctcbn);
        }


        if (eps == 0){ // ie, epsilon not specified
          epsilon=0.5; // initial value

          /* Estimate epsilon | lambda */
          total_loglik = EM_epsilon(&M, D, N_u, Prob, condprob, &epsilon);

          /* Estimate epsilon and lambda */
          int N_compatible = 0;
          int N = 0;
          double alpha;
          for (k=0; k<N_u; k++)
          {
            N += D[k].count;
            N_compatible += D[k].is_compatible * D[k].count;
          }
          alpha = (double) N_compatible / (double) N;

          epsilon = 0.00001;
          total_loglik = EM_EM(&M, D, N_u, lambda, &epsilon);

          char* buffer = (char*)malloc(512);
          if (buffer == NULL) {
            rOutput = "N";
          } else {
            char* lambda_str = print_double_array(lambda, M.n+1);
            snprintf(buffer, 512, "0.0 %f %g %g %s",  epsilon, alpha, total_loglik, lambda_str ? lambda_str : "");
            free(lambda_str);
            rOutput = buffer;
          }
        }
        else{ //ie, given epsilon
          epsilon = eps;
          total_loglik = compute_total_loglik(&M, D, N_u, Prob, condprob, &epsilon);

          char* buffer = (char*)malloc(512);
          if (buffer == NULL) {
            rOutput = "N";
          } else {
            char* lambda_str = print_double_array(lambda, M.n+1);
            snprintf(buffer, 512, "0 %f %g %g %s",  epsilon, pow(1-epsilon,M.n), total_loglik, lambda_str ? lambda_str : "");
            free(lambda_str);
            rOutput = buffer;
          }
        }

        /*Write output?*/
        if (w_flag){
          write_poset(0, ofilestem, M.P, M.n, -1);
          write_lambda(ofilestem, lambda, M.n);
        }

        /* Output probability of each observation Y */
        if (p_flag){
          double* ProbY = compute_ProbY(&M, Prob, only_falsepos, epsilon);
          write_double_array(ofilestem, ".prY", ProbY, pow2(M.n));
          free(ProbY);
        }

        /* Local search */
        if (s_flag)
        {
          if (!n_flag)
            N_iter =  2 * M.n * M.n;
          // printf("\n---Local search---\n");
          // local_search(&M, D, N_u, lambda, &epsilon, total_loglik, T, N_iter, filestem, l_flag, w_flag);
        }

        for (i=0; i < M.m; i++)
          free(condprob[i]);
        free(condprob);
        free(Prob);
        free(lambda_exit);
        free(lattice_index);
      }
      free_data(D, N_u, M.n);
      /* Print final poset */
      // for (i = 0; i< pow2(M.n); i++)
      // print_int_array(GENOTYPE[i + pow2(M.n)], M.n+1);
      if(m_flag) ML_path(&M, lambda);
      if (M.J_P != NULL)
        free_lattice(&M);
      free_poset(&M);
    }
  }

  int n_genotypes = pow2(M.n + 1);
  for (i = 0; i < n_genotypes; i++)
    free(GENOTYPE[i]);
  free(GENOTYPE);
  GENOTYPE = NULL;

  free(lambda);
  free(ofilestem);
  free(filestem1);
  free(filestem2);
  free(filestem3);
  gsl_rng_free(RNG);
  RNG = NULL;

  SEXP res = char_to_sexp(rOutput);
  if (rOutput != NULL && strcmp(rOutput, "N") != 0)
    free(rOutput);

  return res;
}
