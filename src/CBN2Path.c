#include <stdio.h>

#include "ct-cbn.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

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

char* read_file(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return NULL;
    }

    // Determine the file size
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    // Allocate memory for the file content (+1 for the null-terminator)
    char *buffer = (char *)malloc((file_size + 1) * sizeof(char));
    if (buffer == NULL) {
        perror("Memory allocation error");
        fclose(file);
        return NULL;
    }

    // Read the file into the buffer
    size_t read_size = fread(buffer, sizeof(char), file_size, file);
    if (read_size != file_size) {
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
  const char *ofilestem = CHAR(STRING_ELT(ofs, 0));
  const char *filestem1 = CHAR(STRING_ELT(fs1, 0));
  const char *filestem2 = CHAR(STRING_ELT(fs2, 0));

  // defaults:
  double eps = REAL(epsilon)[0];                         // e
  int R = INTEGER(emr)[0];                                // # of EM runs
  double S = REAL(sr)[0];                           // sampling rate \lambda_s
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

  // no epsilon provided from R
  if (eps >= 1.0) {
    e_flag = 0;
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
    if (e_flag)
      mode = LEARN_BOTH;
    else
      mode = LEARN_PARAM;
  }

  srand(seed);
  RNG = gsl_rng_alloc(gsl_rng_taus); // global variable
  gsl_rng_set(RNG, seed);            // seed rng

  int i, j, k;

  model M;
  read_poset(filestem1, &M);

  // precompute binary expansions
  precompute_binary(M.n + 1);

  M.lin_ext = get_int_array(M.n);             // a linear extension of the poset
  double *lambda = get_double_array(M.n + 1); // Exp rates
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

    if (eps >= 0.0) // fixed epsilon
    {
      int b;
      eps = MIN(eps, 1.0);

      // single run:
      output = select_poset(0, eps, &M, lambda, D, N_u, R, mode, 1);
      if (e_flag)
        write_poset(0, ofilestem, M.P, M.n, -1);
      write_lambda(ofilestem, lambda, M.n);

      // bootstrap runs:
      double *p_orig = get_double_array(N_u); // frequencies of original data
      for (k = 0; k < N_u; k++)
        p_orig[k] = (double)D[k].count;

      int **bootstrap_count = get_int_matrix(M.n + 1, M.n + 1);       // all relations
      int **bootstrap_cover_count = get_int_matrix(M.n + 1, M.n + 1); // cover relations only
      for (i = 0; i <= M.n; i++)
        for (j = 0; j <= M.n; j++)
          bootstrap_count[i][j] = bootstrap_cover_count[i][j] = 0;

      int **T = get_int_matrix(M.n + 1, M.n + 1);
      for (b = 1; b <= B; b++)
      {
        resample(D, p_orig, N_u);
        output = select_poset(b, eps, &M, lambda, D, N_u, R, mode, 1); // e_flag

        int_matrix_sum(bootstrap_cover_count, M.P, bootstrap_cover_count, M.n + 1);
        transitive_closure(M.P, T, M.n + 1);
        int_matrix_sum(bootstrap_count, T, bootstrap_count, M.n + 1);

        if (e_flag)
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
      }
      free(T);
      free(bootstrap_count);
      free(p_orig);
    }
    else // sample epsilon
    {
      double epsilon;
      int k;
      for (k = 0; k < N; k++)
      {
        epsilon = (double)k / (2.0 * (double)N);
        output = select_poset(k, epsilon, &M, lambda, D, N_u, R, LEARN_BOTH, 1);
        if (e_flag)
          write_poset(k, ofilestem, M.P, M.n, -1);
      }
    }
    free_data(D, N_u, M.n);
    free_poset(&M);
  }
  else // sample from given model:
  {
    // construct model:
    M.lin_ext = get_int_array(M.n); // a linear extension of the poset
    M.J_P = bfs_order_ideals(M.P, M.n + 1, &(M.m), M.lin_ext);
    parents(&M);
    children(&M);
    if (verbose)
      print_model(&M);

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

    free_model(&M);
  }

  free(lambda);

  return char_to_sexp(output);
}

SEXP hcbn_(SEXP ofs, SEXP fs1, SEXP fs2, SEXP s, SEXP temp, SEXP n)
{
  const char *ofilestem = CHAR(STRING_ELT(ofs, 0));
  const char *filestem1 = CHAR(STRING_ELT(fs1, 0));
  const char *filestem2 = CHAR(STRING_ELT(fs2, 0));
  
  // defaults:
  double eps = 0.0;  // e
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
  
  int n_temp = INTEGER(n)[0];
  if (n_temp > 0) {
    n_flag = 1;
    N_iter = n_temp;
  }
  
  srand(seed);
  RNG = gsl_rng_alloc(gsl_rng_taus);  // global variable
  gsl_rng_set(RNG, seed);  // seed rng
  
  int i, k;
  
  model M;
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

      if (eps >= 0.0)  // fixed epsilon
      {
        eps = MIN(eps, 1.0);

        /*  single run: */
        select_poset(0, eps, &M, lambda, D, N_u, R, mode,0);

        if(eps > 0){
          read_lambda(filestem2, lambda, M.n);
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
          print_double_array(lambda, M.n+1);
        }
        else{ //ie, given epsilon
          epsilon = eps;
          total_loglik = compute_total_loglik(&M, D, N_u, Prob, condprob, &epsilon);
          print_double_array(lambda, M.n+1);
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
        }


        /* GPS */
        if(gps_flag)
        {
          //printf("\n+++GPS+++\n");

          double* all_GPS = get_double_array(M.m);
          double* cond_GPS = get_double_array(N_u);
          double* loglik = get_double_array(N_u);

          GPS(&M, D, N_u, lambda, epsilon, all_GPS, cond_GPS);
          compute_loglik(&M, D, N_u, Prob, condprob, &epsilon, loglik);

          int* ml_pat = get_int_array(N_u);
          double** exp_pat = get_double_matrix(M.n,N_u);
          compute_hidden_patterns(&M, D, N_u, lambda, epsilon, condprob, ml_pat, exp_pat);

          // Write to file
          write_gps(ofilestem, cond_GPS, N, pat_idx);

          // Write MAP to file
          char suffix[15] = ".map";
          char *filename = (char *) calloc(strlen(ofilestem) + strlen(suffix) + 1, sizeof(char));
          strcat(filename, ofilestem);
          strcat(filename, suffix);

          FILE *output;
          if ( (output = fopen(filename, "w")) == NULL)
          {
            // fprintf(stderr, "Error:  Could not write to file %s\n", filename);
            // exit(1);
          }

          for (k=0; k<N;k++)
          {
            c = pat_idx[k]; //Unique observation index
            int ml_gen = M.J_P[ml_pat[c]]; //Lattice entry
            for(i=0;i < M.n + 1;i++)
              fprintf(output, "%i ", GENOTYPE[ml_gen][i]);
            //printf(" %i ", hamdist(index_of(D[c].g, M.n+1), ml_gen));
            //printf(" ");
            //for(i=0; i < M.n; i++)
            //	printf("%.5f ", exp_pat[i][c]);
            //printf(" %.5f ", cond_GPS[c]);
            //printf(" %.5f", all_GPS[ml_pat[c]]);
            //printf(" %.5f\n", loglik[c]);
            fprintf(output, "\n");

          }
          fclose(output);
        }

        /* Local search */
        if (s_flag)
        {
          if (!n_flag)
            N_iter =  2 * M.n * M.n;
          // printf("\n---Local search---\n");
          // local_search(&M, D, N_u, lambda, &epsilon, total_loglik, T, N_iter, filestem, l_flag, w_flag);
        }

      }
      free_data(D, N_u, M.n);
      /* Print final poset */
      // for (i = 0; i< pow2(M.n); i++)
      // print_int_array(GENOTYPE[i + pow2(M.n)], M.n+1);
      if(m_flag) ML_path(&M, lambda);
      free_poset(&M);
    }
  }
  free(lambda);
  
  return char_to_sexp(ofilestem);
}
