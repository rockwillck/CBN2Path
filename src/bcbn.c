/*
 * bcbn.c
 *
 */

#include "dt_cbn.h"
#include "bcbn.h"
#include "bcbn_rand.h"

int bcbn_verbose;
gsl_rng *bcbn_RNG;  // random number generator
int** bcbn_GENOTYPE;  // GENOTYPE[i] is the integer i in binary (as int array)

int* bcbn_get_int_array(const int n)
{
  int* x = calloc(n, sizeof(int));

  if (x == NULL)
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  return x;
}


unsigned int* bcbn_get_uint_array(const int n)
{
  unsigned int* x = calloc(n, sizeof(unsigned int));

  if (x == NULL)
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  return x;
}


double* bcbn_get_double_array(const int n)
{
  double* x = calloc(n, sizeof(double));

  if (x == NULL)
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  return x;
}


int** bcbn_get_int_matrix(const int m, const int n)
{
  int** x = malloc(m * sizeof(int *));

  if (x == NULL)
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  int i;
  for (i=0; i<m; i++)
  {
    x[i] = calloc(n, sizeof(int));
    if (x[i] == NULL)
    {
      // fprintf(stderr, "Error: Out of memory!\n");
      // exit(1);
    }
  }

  return x;
}


double** bcbn_get_double_matrix(const int m, const int n)
{
  double** x = malloc(m * sizeof(double *));

  if (x == NULL)
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  int i;
  for (i=0; i<m; i++)
  {
    x[i] = calloc(n, sizeof(double));
    if (x[i]  == NULL)
    {
      // fprintf(stderr, "Error: Out of memory!\n");
      // exit(1);
    }
  }

  return x;
}

double*** bcbn_get_double_cube(const int m, const int n, const int l)
{
  double*** x = malloc(m * sizeof(double **));

  if (x == NULL)
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  int i;
  for (i=0; i<m; i++)
  {
    x[i] = bcbn_get_double_matrix(n,l);
    if (x[i]  == NULL)
    {
      // fprintf(stderr, "Error: Out of memory!\n");
      // exit(1);
    }
  }

  return x;
}

int*** get_int_cube(const int m, const int n, const int l)
{
  int*** x = malloc(m * sizeof(int **));

  if (x == NULL)
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  int i;
  for (i=0; i<m; i++)
  {
    x[i] = bcbn_get_int_matrix(n,l);
    if (x[i]  == NULL)
    {
      // fprintf(stderr, "Error: Out of memory!\n");
      // exit(1);
    }
  }

  return x;
}

void bcbn_print_int_array(int* x, int n)
{
  // int j;
  //
  // if (n > 0)
  // {
  // 	for (j=0; j<n-1; j++)
  // 		printf("%d ", x[j]);
  // 	printf("%d", x[n-1]);
  // }
  // printf("\n");
}



void bcbn_print_int_matrix(int** X, int m, int n)
{
  int i;

  for (i=0; i<m; i++)
    bcbn_print_int_array(X[i], n);

}


void bcbn_print_double_array(double* x, int n)
{
  // int j;
  //
  // if (n > 0)
  // {
  // 	for (j=0; j<n; j++)
  // 	{
  // 		printf(DOUBLE_FORMAT, x[j]);
  // 		if (j < n-1)
  // 			printf("\t");
  // 	}
  // }
  // printf("\n");
}



void bcbn_print_double_matrix(double** X, int m, int n)
{
  int i;

  for (i=0; i<m; i++)
    bcbn_print_double_array(X[i], n);

}



void bcbn_write_poset(int k, char* filestem, int** P, int n, int b)
{
  int i, j;

  char filename[255];

  if (b >= 0)
    snprintf(filename, sizeof(filename), "%s/b%09d.poset", filestem, b);
  else
    snprintf(filename, sizeof(filename), "%s/%09d.poset", filestem, k);

  FILE *output;
  if ( (output = fopen(filename, "w")) == NULL )
  {
    // fprintf(stderr, "Error:  Could not write to file %s\n", filename);
    // fprintf(stderr, "        Make sure the directory '%s' exists.\n", filestem);
    // exit(1);
  }

  fprintf(output, "%d\n", n);
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      if (P[i][j])
        fprintf(output, "%d %d\n", i+1, j+1);

      fprintf(output, "0 0\n");
      fclose(output);
}



void bcbn_write_patterns(char* filestem, int** pat, int N, int n)
{
  int i, j;

  char suffix[15] = ".sim.pat";
  //char suffix[15] = ".pat";
  char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
  strcat(filename, filestem);
  strcat(filename, suffix);

  FILE *output;
  if ( (output = fopen(filename, "w")) == NULL )
  {
    // fprintf(stderr, "Error:  Could not write to file %s\n", filename);
    // exit(1);
  }

  fprintf(output, "%d %d\n", N, n+1);
  for (i=0; i<N; i++)
  {
    for (j=0; j<n; j++)
      fprintf(output, "%d ", pat[i][j]);
    fprintf(output, "%d\n", pat[i][n]);
  }

  fclose(output);
}




inline void genotype_of(int index, int* x, int n)
{
  int i;

  for (i=n-1; i>=0; i--)
  {
    x[i] = index % 2;
    index = index / 2;
  }

}



int bcbn_pow2(int k)
{
  int i, pow = 1;

  for (i=0; i<k; i++)
    pow *= 2;

  return(pow);
}



void bcbn_precompute_binary(const int n)
{
  int i, j;
  int m = bcbn_pow2(n);
  int *g = bcbn_get_int_array(n);

  bcbn_GENOTYPE = bcbn_get_int_matrix(m, n);

  for (i=0; i<m; i++)
  {
    genotype_of(i, g, n);
    for (j=0; j<n; j++)
      bcbn_GENOTYPE[i][j] = g[j];
  }
  free(g);
}



int bcbn_index_of(int* x, int n)
{
  int i, index = 0;

  for (i=0; i<n; i++)
  {
    if (x[i] == 1)
      index += bcbn_pow2(n-1-i);
  }

  return index;
}

int bcbn_is_equal_int_matrix(int** A, int** B, int n)
{
  int i, j;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      if (A[i][j] != B[i][j])
        return 0;

      return 1;
}



void bcbn_boolean_matrix_sum(int** A, int** B, int** C, int n)
{
  /*
   Boolean matrix sum  A + B = C
   */

  int i, j;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
    {
      C[i][j] = A[i][j] + B[i][j];
      C[i][j] = C[i][j] ? 1 : 0;
    }

}

void bcbn_boolean_matrix_product(int** A, int** B, int** C, int n)
{
  /*
   Boolean matrix product  A * B = C
   */

  int i, j, k;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
    {
      C[i][j] = 0;
      for (k=0; k<n; k++)
        C[i][j] += A[i][k] * B[k][j];
      C[i][j] = C[i][j] ? 1 : 0;
    }

}

void free_int_matrix( int** A, int n) {
  int j;

  for (j=0; j<n; j++)
    free(A[j]);
  free(A);
}

void free_double_matrix( double** A, int n) {
  int j;

  for (j=0; j<n; j++)
    free(A[j]);
  free(A);
}

void free_int_cube( int*** A, int n, int m) {
  int j;

  for (j=0; j<n; j++) {
    //printf("before free int mat\n");
    free_int_matrix(A[j], m );
  }
  //printf("before free A\n");
  //free(A);
}



void bcbn_free_poset(model* M)
{
  free_int_matrix( M->P, M->n );
  free(M->lin_ext);
}



void bcbn_free_lattice(model* M)
{
  int i;

  free(M->J_P);
  for (i=0; i<M->m; i++)
  {
    free(M->pa[i]);
    free(M->pa_diff[i]);
    free(M->ch_diff[i]);
  }
  free(M->pa);
  free(M->pa_diff);
  free(M->N_pa);
  free(M->ch_diff);
}

void bcbn_free_lattice_children(model* M)
{
  int i;

  free(M->J_P);
  free_int_matrix( M->ch_diff, M->m );
  free(M->N_ch);
}



void bcbn_free_model(model* M)
{
  bcbn_free_poset(M);
  bcbn_free_lattice(M);
}



void bcbn_print_model(model* M)
{
  // int i;
  // int m = M->m;
  // int n = M->n;
  // int* g = bcbn_get_int_array(n);
  //
  // printf("\n\nMODEL\n\n");
  // printf("\nP =\n");
  // bcbn_print_int_matrix(M->P, n, n);
  // printf("\n");
  //
  // printf("lattice size, m = %d\nsorted lattice = \n", m);
  // for (i=0; i<m; i++)
  // {
  // 	printf("%d\t", i);
  // 	genotype_of(M->J_P[i], g, n);
  // 	bcbn_print_int_array(g, n);
  // }
  // //bcbn_print_int_array(M->J_P, m);
  //
  // printf("\nlinear extension of the poset = ");
  // bcbn_print_int_array(M->lin_ext, n);
  //
  // for (i=0; i<m; i++)
  // {
  // 	printf("\nparents of %5d = ", i);
  // 	//genotype_of(M->J_P[i], g, n+1);
  // 	//bcbn_print_int_array(g, n+1);
  // 	bcbn_print_int_array(M->pa[i], M->N_pa[i]);
  // 	printf("differing events = ");
  // 	bcbn_print_int_array(M->pa_diff[i], M->N_pa[i]);
  // }
  // printf("\n\n");
  //
  // for (i=0; i<m; i++)
  // {
  // 	printf("differing events to children of %5d = ", i);
  // 	//genotype_of(M->J_P[i], g, n+1);
  // 	//bcbn_print_int_array(g, n+1);
  // 	bcbn_print_int_array(M->ch_diff[i], M->N_ch[i]);
  // }
  //
  // printf("\n");
  //
  // free(g);
}



void bcbn_print_data(data* D, int N_u, int n, int m)
{
  // int k;
  //
  // printf("\n\nDATA\n\n");
  // for (k=0; k<N_u; k++)
  // {
  // 	printf("# %d\n", k);
  // 	printf("g = ");  bcbn_print_int_array(D[k].g, n);
  // 	//printf("Q =\n");  bcbn_print_int_matrix(D[k].Q, n, n);
  // 	//printf("J_Q = ");  bcbn_print_int_array(D[k].J_Q, m);
  // 	printf("count = %d\n", D[k].count);
  // 	printf("is compatible = %d\n", D[k].is_compatible);
  // 	printf("--------------------\n");
  // }

}



void bcbn_free_data(data* D, int N_u, int n)
{
  int k, j;

  for (k=0; k<N_u; k++)
  {
    free(D[k].g);
    for (j=0; j<=n; j++)
      free(D[k].Q[j]);
    free(D[k].Q);
    free(D[k].J_Q);
  }
  free(D);
}



int** bcbn_read_patterns(char* filestem, int* N, int n)
{
  int j, k, p;

  char suffix[15] = ".pat";
  char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
  strcat(filename, filestem);
  strcat(filename, suffix);

  FILE *input;
  if ( (input = fopen(filename, "r")) == NULL)
  {
    // fprintf(stderr, "Error:  Could not read %s\n", filename);
    // exit(1);
  }

  /* Read dimensions */
  if(fscanf(input, "%d %d", N, &p));
  // if (bcbn_verbose) printf("\nreading data from file %s :  %d samples, %d events ...\n\n", filename, *N, p-1);
  if (*N < 1)
  {
    // fprintf(stderr, "Error:  Less than one data point!\n");
    // exit(1);
  }
  if (n != p)
  {
    // fprintf(stderr, "Error:  Number of events in poset and data do not match!\n");
    // exit(1);
  }

  int** pat = bcbn_get_int_matrix(*N, p);

  /* Read patterns */
  int x;
  for (k=0; k<*N; k++)
  {
    for (j=0; j<n; j++)
    {
      if (fscanf(input,"%d ", &x) == 1)
      {
        if ((x != 0) && (x != 1))
          x = -1;
        pat[k][j] = x;
      }
      else
      {
        // fprintf(stderr, "Error reading data from %s!\n", filename);
        // exit(1);
      }

    }
  }

  return pat;
}

void bcbn_read_poset(char* filestem, model* M)
{
  int left, right;

  char suffix[15] = ".poset";
  char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
  strcat(filename, filestem);
  strcat(filename, suffix);

  FILE *input;
  if ( (input = fopen(filename, "r")) == NULL)
  {
    // fprintf(stderr, "Error:  Could not read %s\n", filename);
    // exit(1);
  }

  /* Read number of relations */
  int n;
  if(fscanf(input, "%d", &n));
  // if (bcbn_verbose)  printf("n = %d events\n\n", n);
  if ((n < 1) || (n > 25))
  {
    // fprintf(stderr, "Error:  Number of events is %d.  Supported range is {1, ..., 14}.\n", n);
    // exit(1);
  }

  M->n = n;
  M->P = bcbn_get_int_matrix(n+1, n+1);

  /* Read partial orderings from file */
  if(fscanf(input,"%d %d", &left, &right));
  while (left != 0)
  {
    // if (bcbn_verbose)  printf("%d --> %d\n", left, right);  // i.e., left < right
    if ((left > n) || (right > n) || (left < 0) || (right < 1))
    {
      // fprintf(stderr, "Error:  Undefined event in %s!\n", filename);
      // exit(1);
    }
    M->P[left][right] = 1;
    if(fscanf(input,"%d %d", &left, &right));
  }

  fclose(input);
}

void bcbn_print_genotype(int* x, int n)
{
  int i;

  // for (i=0; i<n; i++)
  // 	printf("%d", x[i]);

}



int* bcbn_bfs_order_ideals(int** poset, const int len, int* count, int* lin_ext)
{
  // implements BFS in the genotype lattice

  int i, j, k, is_compatible;
  int lin_ext_size = 0;

  queue q;
  int g_idx, new_idx;  // current, new genotype index
  int* g = bcbn_get_int_array(len);  // current genotype

  int* lattice = malloc(sizeof(int));

  int* added = bcbn_get_int_array(bcbn_pow2(len));  // records added genotypes
  /* TODO: make this dynamic, e.g., using a list */

  init_queue(&q);
  new_idx = 0;  // wild type 0...0
  enqueue(&q, new_idx);
  added[new_idx] = 1;

  *count = 0;
  while (empty(&q) == FALSE) {
    g_idx = dequeue(&q);
    genotype_of(g_idx, g, len);

    // visit...
    (*count)++;
    if ((lattice = (int*) realloc(lattice, (*count) * sizeof(int))) == NULL)
    {
      // fprintf(stderr, "Error: Out of memory!\n");
      // exit(1);
    }
    lattice[(*count)-1] = g_idx;

    // linear extension
    for (i=1; i<len; i++)  // exclude 0
    {
      if (g[i] == 1)
      {
        int is_in = 0;
        for (j=0; j<lin_ext_size; j++)
          if (lin_ext[j] == i)
          {
            is_in = 1;
            break;
          }
          if (! is_in)  // add to linear extension:
          {
            lin_ext[lin_ext_size] = i;
            lin_ext_size++;
          }
      }
    }


    // generate children:
    for (i=0; i<len; i++)
    {
      if (g[i] == 0)
      {
        g[i] = 1;  // try this new event
        new_idx = bcbn_index_of(g, len);

        // check bcbn_compatibility:
        is_compatible = 1;
        for (k=0; k<len; k++)
        {
          if ((poset[k][i] == 1) && (g[k] == 0))
          {
            is_compatible = 0;
            break;
          }
        }

        if ((is_compatible) && (! added[new_idx]))
          /* add if compatible and really new */
        {
          enqueue(&q, new_idx);
          added[new_idx] = 1;
        }

        g[i] = 0;  // undo event
      }
    }
  }

  free(g);
  free(added);

  return(lattice);
}



int bcbn_norm1(const int g_idx, const int n)
{
  int i;
  int *g = bcbn_get_int_array(n);

  genotype_of(g_idx, g, n);

  int bcbn_norm1 = 0;
  for(i=0; i<n; i++)
    bcbn_norm1 += g[i];

  free(g);

  return bcbn_norm1;
}


double bcbn_power(double m, int n)
{
  int i;
  double x=1;
  for(i=0;i<n;i++)
    x*=m;
  return x;
}

unsigned bcbn_hamdist(unsigned x, unsigned y)
{
  unsigned dist = 0, val = x ^ y;

  while(val)
  {
    ++dist;
    val &= val - 1;
  }

  return dist;
}

int bcbn_hamming_distance(int g_idx, int h_idx, int* diff_idx, int n)
{
  int i;
  int *g = bcbn_get_int_array(n);
  int *h = bcbn_get_int_array(n);

  genotype_of(g_idx, g, n);
  genotype_of(h_idx, h, n);

  int dist = 0;
  for(i=0; i<n; i++)
    if (g[i] != h[i])
    {
      dist++;
      *diff_idx = i;
    }

    free(g);
    free(h);

    return dist;
}

void parents_dt(model* M)
{
  int i, j, k, c;
  int m = M->m;
  int n = M->n;

  M->N_pa = bcbn_get_int_array(m);
  M->pa = malloc(m * sizeof(int*));
  M->pa_diff = malloc(m * sizeof(int*));
  if ((M->pa == NULL) || (M->pa_diff == NULL))
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  // count parents:
  for (i=0; i<m; i++)
  {
    M->N_pa[i] = 0;  // number of parents
    k = i-1;
    while ((k >=0) && (bcbn_norm1(M->J_P[k], n) >= (bcbn_norm1(M->J_P[i], n) - 1)))
    {
      //if (bcbn_hamming_distance(M->J_P[i], M->J_P[k], &j, n+1) == 1)
      if (bcbn_hamdist(M->J_P[i], M->J_P[k]) == 1)
      {
        M->N_pa[i]++;  // found a parent!
      }
      k--;
    }
  }

  // list parents:
  for (i=0; i<m; i++)
  {
    M->pa[i] = bcbn_get_int_array(M->N_pa[i]);
    M->pa_diff[i] = bcbn_get_int_array(M->N_pa[i]);
    // generate parents in sublattice
    k = i-1;
    c = 0;
    while ((k >=0) && (bcbn_norm1(M->J_P[k], n) >= (bcbn_norm1(M->J_P[i], n) - 1)))
    {
      if(bcbn_hamming_distance(M->J_P[i], M->J_P[k], &j, n) == 1)
      {
        M->pa[i][c] = k;  // save index of parent in J_P[]
        M->pa_diff[i][c] = j;  // save index by which parent and child differ
        c++;
      }
      k--;
    }
  }

}



void children_dt(model* M)
{
  int i, j, k, c;
  int m = M->m;
  int n = M->n;

  M->N_ch = bcbn_get_int_array(m);
  M->ch_diff = malloc(m * sizeof(int*));
  if (M->ch_diff == NULL)
  {
    // fprintf(stderr, "Error: Out of memory!\n");
    // exit(1);
  }

  // count children:
  for (i=0; i<m; i++)
  {
    M->N_ch[i] = 0;
    /* generate children in lattice */
    k = i + 1;
    while ( (k < m) && (bcbn_norm1(M->J_P[k], n) <= (bcbn_norm1(M->J_P[i], n) + 1)) )
    {
      //if (bcbn_hamming_distance(M->J_P[i], M->J_P[k], &j, n+1) == 1)
      if (bcbn_hamdist(M->J_P[i], M->J_P[k]) == 1)
      {
        M->N_ch[i]++;
      }
      k++;
    }
  }

  // list children:
  for (i=0; i<m; i++)
  {
    M->ch_diff[i] = bcbn_get_int_array(M->N_ch[i]);
    /* generate children in lattice */
    k = i + 1;
    c = 0;
    while ( (k < m) && (bcbn_norm1(M->J_P[k], n) <= (bcbn_norm1(M->J_P[i], n) + 1)) )
    {
      if (bcbn_hamming_distance(M->J_P[i], M->J_P[k], &j, n) == 1)
      {
        M->ch_diff[i][c] = j;
        c++;
      }
      k++;
    }
  }

}

void bcbn_read_poset_dt(char* filestem, model* M)
{
  int left, right;

  char suffix[15] = ".poset";
  char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
  strcat(filename, filestem);
  strcat(filename, suffix);

  FILE *input;
  if ( (input = fopen(filename, "r")) == NULL)
  {
    // fprintf(stderr, "Error:  Could not read %s\n", filename);
    // exit(1);
  }

  /* Read number of relations */
  int n;
  if(fscanf(input, "%d", &n));
  // if (bcbn_verbose)  printf("n = %d events\n\n", n);
  if ((n < 1) || (n > 25))
  {
    // fprintf(stderr, "Error:  Number of events is %d.  Supported range is {1, ..., 14}.\n", n);
    // exit(1);
  }

  M->n = n;
  M->P = bcbn_get_int_matrix(n, n);

  /* Read partial orderings from file */
  if(fscanf(input,"%d %d", &left, &right));
  while (left != 0)
  {
    // if (bcbn_verbose)  printf("%d --> %d\n", left, right);  // i.e., left < right
    if ((left > n) || (right > n) || (left < 0) || (right < 1))
    {
      // fprintf(stderr, "Error:  Undefined event in %s!\n", filename);
      // exit(1);
    }
    M->P[left-1][right-1] = 1;
    if(fscanf(input,"%d %d", &left, &right));
  }

  fclose(input);
}

data* bcbn_make_data_set(int** pat, int N, int n, int* N_u, int* pat_idx)
{
  int j, k, l, skip;

  int* idx = bcbn_get_int_array(N);
  int* count = bcbn_get_int_array(N);

  for (k=0; k<N; k++)
    idx[k] = bcbn_index_of(pat[k], n);

  // count unique patterns:
  int c = 0;
  for (k=0; k<N; k++)
  {
    skip = 0;
    for (l=0; l<k; l++)
    {
      if (idx[k] == idx[l]) // observed before?
      {
        count[l]++;
        pat_idx[k] = pat_idx[l];
        skip = 1;
        break;
      }
    }
    if (! skip) // Unique!
    {
      count[k]++;
      pat_idx[k] = c;
      c++;
    }
  }

  *N_u = 0;
  for (k=0; k<N; k++)
    *N_u += (count[k] > 0);

  // if (bcbn_verbose)  printf("N_u = %d unique patterns\n", *N_u);

  data* D = calloc(*N_u, sizeof(data));

  c = 0;
  for (k=0; k<N; k++)
  {
    //    pat_idx[k] = c;
    if (count[k] > 0)
    {
      //	  pat_idx[k]=c;
      D[c].count = count[k];
      D[c].g = bcbn_get_int_array(n);
      for (j=0; j<n; j++)
        D[c].g[j] = pat[k][j];
      c++;
    }
  }

  return D;
}

void bcbn_compatibility(data* D, int N_u, model* M)
{
  int i, j, k;

  for (k=0; k<N_u; k++)
  {
    D[k].is_compatible = 1;
    i = 0;
    while ((i < M->n) && D[k].is_compatible)
    {
      j = 0;
      while ((j < M->n) && D[k].is_compatible)
      {
        if (M->P[i][j] && (! D[k].g[i]) && D[k].g[j])
          /* It is sufficient to check the cover relations! */
          //if (lt_poset(i, j, M->P, M->n) && (! D[k].g[i]) && D[k].g[j])
          D[k].is_compatible = 0;
        j++;
      }
      i++;
    }
  }

}


/**
 * Pr(X|M,theta) for all X in J(P)
 * @param theta IN mutation probabilities; size n
 * @param M IN the model
 * @param Prob OUT the probabilities of all valid (according to the poset) genotypes; size m
 * @param theta_exit IN the exit probabilities for all elements in the genotype lattice; size: m
 */
void compute_all_cbn_prob( double* theta, model* M, double* Prob, double* theta_exit ) {

  int m = M->m;
  int n = M->n;

  int i,j,c,k;

  Prob[0] = 1.0;

  for (i=0; i < m; i++)
  {
    //printf("%d ",M->J_P[i]);
    //bcbn_print_int_array(bcbn_GENOTYPE[M->J_P[i]],n);
    double t = 1.0;
    for (c=0; c<n; c++)
    {
      if( bcbn_GENOTYPE[M->J_P[i]][c] ) {
        t *= theta[c];
      }
    }
    Prob[i] = t * theta_exit[i];

  }

  double p = 0;
  for(i=0;i<m;i++) {
    p += Prob[i];
  }
  //printf(DOUBLE_FORMAT ,p);
}

/*
 * IN double *theta : mutation probabilities; size: n
 * IN model *M : the model
 * OUT double * theta_exit: the exit probabilities for all elements in the genotype lattice; size: m
 */
void compute_theta_exit(double *theta, model* M, double* theta_exit)
{
  int i, j, c;

  for (i=0; i<M->m; i++)
  {
    theta_exit[i] = 1.0;

    for (c=0; c<M->N_ch[i]; c++)
    {
      j = M->ch_diff[i][c];  // index of differing event
      theta_exit[i] *= (1-theta[j]);

    }
  }
}

// Compute the conditional probability Prob[Z|X,epsilon]
void compute_cond_error_prob(model* M, data* D, int N_u, double epsilon, double** condprob)
{
  int n = M->n;
  int m = M->m;

  int i,k,d,index;


  for (k=0; k<N_u; k++)   // All patients
  {
    //prob_tmp = 0;
    index = bcbn_index_of(D[k].g, n);

    for (i=1; i<m; i++) // All patterns
    {
      condprob[i][k] = 0.0;
      d = bcbn_hamdist(M->J_P[i], index);
      condprob[i][k] = bcbn_power(epsilon, d) * bcbn_power(1 - epsilon, n - d);
      //prob_tmp += condprob[i][k];
    }
    // Now divide by Prob[Y_k]
    //for (i=1; i<m; i++) // All patterns
    //condprob[i][k] /= prob_tmp;
  }
}


void clone_poset( model* M_src, model* M_dest ) {

  int i,j;
  int n = M_src->n;
  //int m = M_src->m;
  M_dest->n = n;
  M_dest->P = bcbn_get_int_matrix(n, n);
  M_dest->lin_ext = bcbn_get_int_array( n );
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      M_dest->P[i][j] = M_src->P[i][j];
    }
    M_dest->lin_ext[i] = M_src->lin_ext[i];
  }
  /*
   M_dest->m = m;
   M_dest->J_P = bcbn_get_int_array( m );
   for( i=0;i<m;i++ ) {
   M_dest->J_P[i] = M_src->J_P[i];
   }
   M_dest->eps = M_src->eps;
   */

}

void bcbn_transitive_closure(int** A, int**T, int n)
{
  /*
   T is the transitive closure of the relation A
   */

  int** R = bcbn_get_int_matrix(n, n);
  int** S = bcbn_get_int_matrix(n, n);

  int i, j, k = 0;

  // initialize T = A
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      T[i][j] = A[i][j];

  while ((k == 0) || (! bcbn_is_equal_int_matrix(S, T, n)))
  {
    // S = T
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        S[i][j] = T[i][j];

    // T = A*S + A
    bcbn_boolean_matrix_product(A, S, R, n);
    bcbn_boolean_matrix_sum(R, A, T, n);

    k++;
  }

  free_int_matrix(R, n);
  free_int_matrix(S, n);

}

/**
 *
 * @param P transitive closure to cover relations
 * @param n
 * @param C changes
 * @return status
 */
int bcbn_reduce_to_cover_relations ( int** P, int n, int** C )
{
  int i,j,k;
  queue q;
  int stat = 0;

  // Sample all nodes
  for (i=0; i<n; i++)
  {
    int* visit = bcbn_get_int_array(n);
    init_queue(&q);

    // Fill queue with children of i
    for (j=0; j<n; j++)
      if (P[i][j])
      {
        enqueue(&q,j);
      }

      // Walk through grandchildren
      while (empty(&q) == FALSE)
      {
        j = dequeue(&q);
        for (k=0; k<n; k++)
        {
          if (P[j][k] && !(visit[k]))
          {
            visit[k] = 1;
            enqueue(&q,k);

            // Remove non-cover relations
            if (P[i][k])
            {
              P[i][k] = 0;
              C[i][k] = 1;
              stat = 1; // Report changes
            }

            // Check if cyclic
            if (P[k][i]) {
              free(visit);
              return 2; // Fatal
            }
          }
        }
      }
      free(visit);
  }
  return stat;
}

#define LCG_MODULUS 4294967296U // 2^32
#define LCG_MULTIPLIER 1103515245U
#define LCG_INCREMENT 12345U
// static uint32_t seed = 1; // Seed value
// void srand(uint32_t new_seed) {
//   seed = new_seed;
// }

// int rand() {
//   // Update the seed using LCG formula
//   seed = (LCG_MULTIPLIER * seed + LCG_INCREMENT) % LCG_MODULUS;
//   // Convert the seed to a double between 0 and 1
//   return seed;
// }

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
long double compute_likelihood( double eps, model* M, double* theta, data* D, int N_u, double* Prob, double** cond_err_prob ) {

  int i,j;
  int m = M->m;

  long double likelihood, likelihood_d;

  //likelihood = 1;
  long double loglik = 0;
  //product over observed genotypes (i.e. data)
  for( i=0;i<N_u;i++ ) {
    //sum over valid genotypes
    likelihood_d = 0;
    for( j=0;j<m;j++ ) {
      long double tempr = (long double) (Prob[j] * cond_err_prob[j][i]);
      likelihood_d += tempr;
    }
    //likelihood *= bcbn_power( likelihood_d, D[i].count );
    long double temprr = D[i].count;
    loglik += temprr * log10l( likelihood_d );
  }
  //printf("\n likelihood: %20Lg\n", likelihood);
  //printf("\n loglik: %20Lg\n", logl( likelihood ) );
  //return likelihood;
  return(loglik);
}


/**
 * unnorm P[M,theta,eps|D]
 * @param eps
 * @param M
 * @param theta
 * @param D
 * @param N_u
 * @return unnormalized posterior
 */
long double compute_unnomr_posterior( double eps, model*M, double* theta, data* D, int N_u ) {

  int m = M->m;
  int n = M->n;

  double* Prob = bcbn_get_double_array(m);
  double* theta_exit = bcbn_get_double_array(m);
  compute_theta_exit( theta, M, theta_exit);
  compute_all_cbn_prob( theta, M, Prob, theta_exit );
  free(theta_exit);
  double** cond_err_prob = bcbn_get_double_matrix(m, N_u);
  compute_cond_error_prob( M, D, N_u, eps, cond_err_prob );
  long double lik;
  lik = compute_likelihood(eps, M, theta, D, N_u, Prob, cond_err_prob );
  free(Prob);
  int i;
  free_double_matrix( cond_err_prob, m );

  double t=1;
  double alpha, beta;
  for( i=0;i<n;i++ ) {
    //TODO: treat alpha and beta like parameters not like f***-given constants
    alpha = 0.00001;
    beta = 0.00001;
    alpha=beta=1;
    t *= gsl_ran_beta_pdf( theta[i], alpha, beta);
  }
  //TODO: treat alpha and beta like parameters not like f***-given constants
  double eps_pdf = gsl_ran_beta_pdf( eps, 5, 30 );
  //long double posterior;
  //posterior = lik*t*eps_pdf;

  long double temp1 = (long double) log10(t);
  long double temp2 = (long double) log10(eps_pdf);
  long double lposterior = (lik+temp1+temp2);
  //return posterior;
  return lposterior;
}

void relocate_whole_theta( double* theta, double* theta_p, int n ) {
  int i;
  double alpha,beta,var,x;
  for( i=0;i<n;i++ ) {
    x = theta[i];
    //TODO: set alpha and beta!!!
    beta = 1;
    var = x * (1 - x) / ( beta / ( 1 - x) + 1 );
    alpha = x * ( x * ( 1 - x ) / var - 1 );
    alpha = beta = 1;
    theta_p[i] = gsl_ran_beta (bcbn_RNG, alpha, beta);
  }
}

void relocate_theta_i( double* theta, double* theta_p, int n, int i ) {
  double alpha,beta,var,x;
  x = theta[i];
  //TODO: set alpha and beta!!!
  //beta = 1;
  //var = x * (1 - x) / ( beta / ( 1 - x) + 1 );
  //alpha = x * ( x * ( 1 - x ) / var - 1 );
  alpha = beta = 1;
  theta_p[i] = gsl_ran_beta (bcbn_RNG, alpha, beta);
}

double compute_theta_transition_prob( double* theta, double* theta_p, int n ) {
  double t = 1;
  int i;
  double alpha,beta,var,x;
  for( i=0;i<n;i++ ) {
    x = theta[i];
    //TODO: set alpha and beta!!!
    //beta = 1;
    //var = x * (1 - x) / ( beta / ( 1 - x) + 1 );
    //alpha = x * ( x * ( 1 - x ) / var - 1 );
    alpha = beta = 1;
    t *= gsl_ran_beta_pdf(theta_p[i], alpha, beta);
  }
  return t;
}


void valid_new_cover_relation_matrix( model *M, int** V ) {
  int n = M->n;
  int i,j,k,l;
  int** T = bcbn_get_int_matrix(n, n);
  int** C = bcbn_get_int_matrix(n, n);

  for( i=0;i<n;i++ ) {
    for( k=0;k<n;k++ ) {
      C[i][k] = 0;
      T[i][k] = 0;
    }
  }
  int stat;
  bcbn_transitive_closure( M->P, T, n );

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( !T[i][j] && i!=j ) {
        T[i][j] = 1;
        stat = bcbn_reduce_to_cover_relations( T, n, C );
        //check changes
        if( stat == 1 && !C[i][j] ) {
          for( k=0;k<n;k++ ) {
            for( l=0;l<n;l++ ) {
              if( C[k][l] == 1 && M->P[k][l] == 1 ) {
                T[i][j] = 0;
              }
            }
          }
          V[i][j] = T[i][j];
        }
        //cycle!!
        else if( stat == 2) {
          T[i][j] = 0;
        }
        else if( stat == 0 ){
          V[i][j] = 1;
        }
      }

      for( l=0;l<n;l++ ) {
        for( k=0;k<n;k++ ) {
          T[l][k] = 0;
          C[l][k] = 0;
        }
      }
      bcbn_transitive_closure( M->P, T, n );
    }
  }
  free_int_matrix( T, n );
  free_int_matrix( C, n );
}

int propose_new_cover_relation( model *M, double* tp, data* D, int N_u ) {
  int n = M->n;
  int i,j,k,N_compatible,N_all_comp;
  int** V = bcbn_get_int_matrix(n, n);
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ )
      V[i][j] = 0;
  }
  valid_new_cover_relation_matrix( M, V );
  int valid_c = 0;
  int c = 0;
  N_all_comp = 0;

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( V[i][j] ) {
        valid_c++;
        //M->P[i][j] = 1;
        //bcbn_compatibility( D, N_u, M );
        //N_compatible = 0;
        //for (k=0; k<N_u; k++) {
        //N_compatible += D[k].is_compatible * D[k].count;
        //N_compatible = 1;
        //}
        N_compatible = 1;
        V[i][j] = N_compatible;
        N_all_comp += N_compatible;
        //M->P[i][j] = 0;
      }
    }
  }
  //TODO: what if N_compatible is 0 for a specific instert???
  if (valid_c) {
    //(*tp) = (double)1/valid_c;
    int new_index = bcbn_pcg_rand()%N_all_comp;
    for( i=0;i<n;i++ ) {
      for( j=0;j<n;j++ ) {
        if( V[i][j] ) {
          c += V[i][j];
          if( new_index < c ) {
            M->P[i][j] = 1;
            (*tp) = (double)V[i][j]/N_all_comp;
            new_index = N_all_comp + 1;
          }
        }
      }
    }
  }
  else {
    free_int_matrix( V, n );
    return 0;
  }
  free_int_matrix( V, n );
  return 1;
}

double get_tp_for_remove_cover_move( model* M ) {
  int i,j;
  int delete_c = 0;
  int n = M->n;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if ( M->P[i][j] ) {
        delete_c++;
      }
    }
  }
  return (double)1/delete_c;
}

int propose_delete_cover_relation( model *M, double* tp ) {
  int valid_c = 0;
  int c = 0;
  int i,j;
  int n = M->n;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( M->P[i][j] ) {
        valid_c++;
      }
    }
  }
  if (valid_c) {
    (*tp) = (double)1/valid_c;
    int new_index = bcbn_pcg_rand()%valid_c;
    for( i=0;i<n;i++ ) {
      for( j=0;j<n;j++ ) {
        if( M->P[i][j] ) {
          if( new_index == c ) {
            M->P[i][j] = 0;
          }
          c++;
        }
      }
    }
  }
  else {
    return 0;
  }
  return 1;
}

double get_tp_for_new_cover_move( model* M, data* D, int N_u ) {
  int i,j,N_all_comp,N_compatible,k;
  N_all_comp = 0;
  int c = 0;
  int n = M->n;
  double tp = 0;
  int** V = bcbn_get_int_matrix(n, n);
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      V[i][j] = 0;
    }
  }
  valid_new_cover_relation_matrix( M, V );
  int valid_c = 0;
  N_all_comp = 0;

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( V[i][j] ) {
        valid_c++;
        //M->P[i][j] = 1;
        //bcbn_compatibility( D, N_u, M );
        //N_compatible = 0;
        //for (k=0; k<N_u; k++) {
        //N_compatible += D[k].is_compatible * D[k].count;
        //N_compatible = 1;
        //}
        //N_compatible = 1;
        //V[i][j] = N_compatible;
        //N_all_comp += N_compatible;
        //M->P[i][j] = 0;
      }
    }
  }

  if (valid_c) {
    //(*tp) = (double)1/valid_c;
    int new_index = bcbn_pcg_rand()%valid_c;
    for( i=0;i<n;i++ ) {
      for( j=0;j<n;j++ ) {
        if( V[i][j] ) {
          c += V[i][j];
          if( new_index < c ) {
            tp = (double)V[i][j]/valid_c;
            new_index = valid_c + 1;
          }
        }
      }
    }
  }
  free_int_matrix( V, n );
  return tp;
}

void propose_event_exchange_move( model* M, double* tp, int *ti, int *tj ) {
  int n = M->n;
  int i,j,k;
  int* h_in = bcbn_get_int_array( n );
  int* h_out = bcbn_get_int_array( n );
  int h=0;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      h += M->P[i][j];
    }
  }

  i = j = bcbn_pcg_rand()%n;
  while( i==j ) {
    j = bcbn_pcg_rand()%n;
  }

  int ij,ji;
  ij = M->P[i][j];
  ji = M->P[j][i];

  for( k=0;k<n;k++ ) {
    h_out[k] = M->P[i][k];
    h_in[k] = M->P[k][i];
  }
  for( k=0;k<n;k++ ) {
    if ( i != k ) {
      M->P[i][k] = M->P[j][k];
      M->P[k][i] = M->P[k][j];
    }
  }
  for( k=0;k<n;k++ ) {
    if ( j != k ) {
      M->P[j][k] = h_out[k];
      M->P[k][j] = h_in[k];
    }
  }
  M->P[i][j] = ji;
  M->P[j][i] = ij;

  (*ti) = i;
  (*tj) = j;


  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      h -= M->P[i][j];
    }
  }

  if ( h != 0 ) {
    // printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
  }

  (*tp) = 2/(double)(n*(n-1));

}


//death move followed by birth move
int propose_reincarnation_move( model* M, data* D, int N_u ) {
  double tp;
  if (! propose_delete_cover_relation( M, &tp) ) {
    return 0;
  }
  else {
    propose_new_cover_relation( M, &tp, D, N_u );
  }
  return 1;
}

void relocate_epsilon( double eps, double* epsilon_p ) {
  double alpha,beta,var,x;

  //BetaDist(2,20) is a bit broader than (5,30) used for the prior
  (*epsilon_p) = gsl_ran_beta (bcbn_RNG, 2, 20);
}

double get_tp_epsilon_relocation( double epsilon_p ) {

  return gsl_ran_beta_pdf( epsilon_p, 2, 20);
}

void valid_new_transitive_closer_relation_matrix( model *M, int** V ) {
  int n = M->n;
  int i,j,k,l;
  int** T = bcbn_get_int_matrix(n, n);
  int** T2 = bcbn_get_int_matrix(n, n);
  int** C = bcbn_get_int_matrix(n, n);
  int new_n_edges;
  int orig_n_edges = 0;
  int changes;

  for( i=0;i<n;i++ ) {
    for( k=0;k<n;k++ ) {
      C[i][k] = 0;
      T[i][k] = 0;
      T2[i][k] = 0;
    }
  }
  int stat;
  bcbn_transitive_closure( M->P, T, n );
  for( l=0;l<n;l++ ) {
    for( k=0;k<n;k++ ) {
      orig_n_edges += T[l][k];
    }
  }

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( !T[i][j] && i!=j ) {
        T[i][j] = 1;
        bcbn_transitive_closure( T, T2, n );
        new_n_edges = 0;
        for( l=0;l<n;l++ ) {
          for( k=0;k<n;k++ ) {
            new_n_edges += T2[l][k];
          }

        }
        if ((new_n_edges - orig_n_edges) == 1) {
          stat = bcbn_reduce_to_cover_relations( T2, n, C );
          //check changes
          if( stat == 1 && !C[i][j] ) {
            changes = 0;
            for( l=0;l<n;l++ ) {
              for( k=0;k<n;k++ ) {
                changes += abs(T2[l][k]-M->P[l][k]);
              }
            }
            if ( changes > 2 ) {
              V[i][j] = 1;
            }
          }
        }
        T[i][j] = 0;
        for( l=0;l<n;l++ ) {
          for( k=0;k<n;k++ ) {
            T2[l][k] = 0;
            C[l][k] = 0;
          }
        }
      }
    }
  }
  free_int_matrix( T, n );
  free_int_matrix( T2, n );
  free_int_matrix( C, n );
}

int propose_new_bcbn_transitive_closure_relation( model *M, double* tp, data* D, int N_u ) {
  int n = M->n;
  int i,j,k,N_compatible,N_all_comp;
  int** V = bcbn_get_int_matrix(n, n);
  int** T = bcbn_get_int_matrix(n, n);
  int** C = bcbn_get_int_matrix(n, n);
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      V[i][j] = 0;
      T[i][j] = 0;
    }
  }
  valid_new_transitive_closer_relation_matrix( M, V );
  int valid_c = 0;
  int c = 0;
  N_all_comp = 0;
  bcbn_transitive_closure( M->P, T, n );

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( V[i][j] ) {
        valid_c++;
        //M->P[i][j] = 1;
        //bcbn_compatibility( D, N_u, M );
        //N_compatible = 0;
        //for (k=0; k<N_u; k++) {
        //N_compatible += D[k].is_compatible * D[k].count;
        //N_compatible = 1;
        //}
        N_compatible = 1;
        V[i][j] = N_compatible;
        N_all_comp += N_compatible;
        //M->P[i][j] = 0;
      }
    }
  }
  //TODO: what if N_compatible is 0 for a specific instert???
  if (valid_c) {
    //(*tp) = (double)1/valid_c;
    int new_index = bcbn_pcg_rand()%N_all_comp;
    for( i=0;i<n;i++ ) {
      for( j=0;j<n;j++ ) {
        if( V[i][j] ) {
          c += V[i][j];
          if( new_index < c ) {
            T[i][j] = 1;
            bcbn_reduce_to_cover_relations( T, n, C );
            (*tp) = (double)V[i][j]/N_all_comp;
            new_index = N_all_comp + 1;
          }
        }
      }
    }
    for( i=0;i<n;i++ ) {
      for( j=0;j<n;j++ ) {
        M->P[i][j] = T[i][j];
      }
    }
  }
  else {
    free_int_matrix( V, n );
    free_int_matrix( C, n );
    free_int_matrix( T, n );
    return 0;
  }
  free_int_matrix( V, n );
  free_int_matrix( C, n );
  free_int_matrix( T, n );
  return 1;
}
void valid_delete_bcbn_transitive_closure_relation_matrix( model *M, int** V ) {
  int n = M->n;
  int i,j,k,l;
  int** T = bcbn_get_int_matrix(n, n);
  int** T2 = bcbn_get_int_matrix(n, n);
  int** C = bcbn_get_int_matrix(n, n);
  int changes = 0;

  for( i=0;i<n;i++ ) {
    for( k=0;k<n;k++ ) {
      C[i][k] = 0;
      T[i][k] = 0;
      T2[i][k] = 0;
    }
  }
  int stat;
  bcbn_transitive_closure( M->P, T, n );
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( T[i][j] && i!=j ) {
        T[i][j] = 0;
        bcbn_transitive_closure( T, T2, n );
        if ( !T2[i][j] ) {
          stat = bcbn_reduce_to_cover_relations( T2, n, C );
          //check changes
          if( stat == 1 && !C[i][j] ) {
            changes = 0;
            for( l=0;l<n;l++ ) {
              for( k=0;k<n;k++ ) {
                changes += abs(T2[l][k]-M->P[l][k]);
              }
            }
            if ( changes > 2 ) {
              V[i][j] = 1;
            }
          }
        }
        T[i][j] = 1;
        for( l=0;l<n;l++ ) {
          for( k=0;k<n;k++ ) {
            T2[l][k] = 0;
            C[l][k] = 0;
          }
        }
      }
    }
  }
  free_int_matrix( T, n );
  free_int_matrix( T2, n );
  free_int_matrix( C, n );
}

int propose_delete_bcbn_transitive_closure_relation( model *M, double* tp, data* D, int N_u ) {

  int n = M->n;
  int i,j,k,N_compatible,N_all_comp;
  int** V = bcbn_get_int_matrix(n, n);
  int** T = bcbn_get_int_matrix(n, n);
  int** C = bcbn_get_int_matrix(n, n);
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ )
      V[i][j] = 0;
  }
  valid_delete_bcbn_transitive_closure_relation_matrix( M, V );
  int valid_c = 0;
  int c = 0;
  N_all_comp = 0;

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( V[i][j] ) {
        valid_c++;
        //M->P[i][j] = 1;
        //bcbn_compatibility( D, N_u, M );
        //N_compatible = 0;
        //for (k=0; k<N_u; k++) {
        //N_compatible += D[k].is_compatible * D[k].count;
        //N_compatible = 1;
        //}
        //N_compatible = 1;
        //V[i][j] = N_compatible;
        //N_all_comp += N_compatible;
        //M->P[i][j] = 0;
      }
    }
  }
  //TODO: what if N_compatible is 0 for a specific instert???
  if (valid_c) {
    //(*tp) = (double)1/valid_c;
    int new_index = bcbn_pcg_rand()%valid_c;
    for( i=0;i<n;i++ ) {
      for( j=0;j<n;j++ ) {
        if( V[i][j] ) {
          c += V[i][j];
          if( new_index < c ) {
            T[i][j] = 1;
            bcbn_reduce_to_cover_relations( T, n, C );
            (*tp) = (double)V[i][j]/valid_c;
            new_index = valid_c + 1;
          }
        }
      }
    }
    for( i=0;i<n;i++ ) {
      for( j=0;j<n;j++ ) {
        M->P[i][j] = T[i][j];
      }
    }
  }
  else {
    free_int_matrix( V, n );
    free_int_matrix( C, n );
    free_int_matrix( T, n );
    return 0;
  }
  free_int_matrix( V, n );
  free_int_matrix( C, n );
  free_int_matrix( T, n );
  return 1;
}

double get_tp_for_delete_bcbn_transitive_closure_relation_move( model* M ) {
  int n = M->n;
  int i,j,k,N_compatible,N_all_comp;
  int** V = bcbn_get_int_matrix(n, n);
  int** T = bcbn_get_int_matrix(n, n);
  int** C = bcbn_get_int_matrix(n, n);
  double tp = 0;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ )
      V[i][j] = 0;
  }
  valid_delete_bcbn_transitive_closure_relation_matrix( M, V );
  int valid_c = 0;
  int c = 0;
  N_all_comp = 0;

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( V[i][j] ) {
        valid_c++;
        //M->P[i][j] = 1;
        //bcbn_compatibility( D, N_u, M );
        //N_compatible = 0;
        //for (k=0; k<N_u; k++) {
        //N_compatible += D[k].is_compatible * D[k].count;
        //N_compatible = 1;
        //}
        //N_compatible = 1;
        //V[i][j] = N_compatible;
        //N_all_comp += N_compatible;
        //M->P[i][j] = 0;
      }
    }
  }
  //TODO: what if N_compatible is 0 for a specific instert???
  if (valid_c) {
    tp = (double)1/valid_c;

  }
  else {
    free_int_matrix( V, n );
    free_int_matrix( C, n );
    free_int_matrix( T, n );
    return 0;
  }
  free_int_matrix( V, n );
  free_int_matrix( C, n );
  free_int_matrix( T, n );
  return tp;
}

double get_tp_for_new_bcbn_transitive_closure_relation_move( model* M ) {
  int n = M->n;
  int i,j,k,N_compatible,N_all_comp;
  int** V = bcbn_get_int_matrix(n, n);
  int** T = bcbn_get_int_matrix(n, n);
  int** C = bcbn_get_int_matrix(n, n);
  double tp = 0;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ )
      V[i][j] = 0;
  }
  valid_new_transitive_closer_relation_matrix( M, V );
  int valid_c = 0;
  int c = 0;
  N_all_comp = 0;

  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( V[i][j] ) {
        valid_c++;
        //M->P[i][j] = 1;
        //bcbn_compatibility( D, N_u, M );
        //N_compatible = 0;
        //for (k=0; k<N_u; k++) {
        //N_compatible += D[k].is_compatible * D[k].count;
        //N_compatible = 1;
        //}
        //N_compatible = 1;
        //V[i][j] = N_compatible;
        //N_all_comp += N_compatible;
        //M->P[i][j] = 0;
      }
    }
  }
  //TODO: what if N_compatible is 0 for a specific instert???
  if (valid_c) {
    tp = (double)1/valid_c;

  }
  else {
    free_int_matrix( V, n );
    free_int_matrix( C, n );
    free_int_matrix( T, n );
    return 0;
  }
  free_int_matrix( V, n );
  free_int_matrix( C, n );
  free_int_matrix( T, n );
  return tp;
}

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
void start_MH( double eps, model* M, double* theta, data* D, int N_u, int burn_in, int number_samples, int record_ith ) {

  /*
   * suffixes of variables in this method
   * *_p ... proposal
   * *_o ... old value
   * *_r ... reverse
   */

  double cm = 0.5;
  double fraction_exchange = 0.03;
  double fraction_reincarnation = 0.05;
  double fraction_birth = 0.2;
  double fraction_death = 0.1;
  double fraction_tr_birth = 0.01;
  double fraction_tr_death = 0.01;
  long double alpha, MH_ratio;
  double epsilon_p;
  double u,ms;
  int n = M->n;
  double* theta_p = bcbn_get_double_array(n);
  memcpy( theta_p, theta, sizeof(double)*n);
  model* M_p = malloc( sizeof(model) );
  int i,j,k = 0;
  long double posterior_o, posterior_p;
  double msp, msp_r, tp, tp_r;
  int accepted = 0;

  int si;
  int n_edges = 0;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( M->P[i][j] ) {
        n_edges++;
      }
    }
  }

  posterior_o = compute_unnomr_posterior(eps, M, theta, D, N_u );
  // printf("posterior_o %.20Lg\n", posterior_o);

  clone_poset( M, M_p );
  M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
  //parents_dt(&M);
  children_dt(M_p);
  epsilon_p = eps;

  FILE *samplelog;
  if ( (samplelog = fopen("data/output/SYNrun10.3sample", "w")) == NULL ) {
    // exit(1);
  }

  for( k=0;k<burn_in+number_samples;k++) {

    ms = gsl_ran_flat (bcbn_RNG, 0, 1);

    if( ms < cm ) {
      //relocate_whole_theta( theta, theta_p, n );

      si = bcbn_pcg_rand()%n;;
      //si = k%n;
      relocate_theta_i( theta, theta_p, n, si );
      tp = compute_theta_transition_prob( theta, theta_p, n );
      tp_r = compute_theta_transition_prob( theta_p, theta, n );
      msp = cm;
      msp_r = cm;
    }
    else {
      if ( ms < cm+fraction_exchange ) {

        int ti,tj;
        propose_event_exchange_move( M_p, &tp, &ti, &tj );
        msp = msp_r = fraction_exchange;
        relocate_theta_i( theta, theta_p, n, ti );
        tp = compute_theta_transition_prob( theta, theta_p, n );
        tp_r = compute_theta_transition_prob( theta_p, theta, n );
        relocate_theta_i( theta, theta_p, n, tj );
        tp = tp*compute_theta_transition_prob( theta, theta_p, n );
        tp_r = tp_r*compute_theta_transition_prob( theta_p, theta, n );
        //tp = tp_r = 1;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation ) {
        if ( !propose_reincarnation_move( M_p, D, N_u ) ) {
          k--;
          continue;
        }
        msp = msp_r = fraction_reincarnation;
        tp = tp_r = 1;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth ) {
        //TODO: ugly: handle the case where no more cover relations are possible nicer!
        if( !propose_new_cover_relation(M_p, &tp, D, N_u ) ) {
          k--;
          continue;
        }
        tp_r = get_tp_for_remove_cover_move(M_p);
        msp = fraction_birth;
        msp_r = fraction_death;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death ) {
        if( !propose_delete_cover_relation(M_p, &tp ) ) {
          k--;
          continue;
        }
        tp_r = get_tp_for_new_cover_move(M_p, D, N_u );
        msp_r = fraction_birth;
        msp = fraction_death;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death+fraction_tr_birth ) {
        if( !propose_new_bcbn_transitive_closure_relation(M_p, &tp, D, N_u ) ) {
          k--;
          continue;
        }
        tp_r = get_tp_for_delete_bcbn_transitive_closure_relation_move(M_p);
        msp = fraction_tr_birth;
        msp_r = fraction_tr_death;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death+fraction_tr_birth+fraction_tr_death ) {
        if( !propose_delete_bcbn_transitive_closure_relation(M_p, &tp, D, N_u ) ) {
          k--;
          continue;
        }
        tp_r = get_tp_for_new_bcbn_transitive_closure_relation_move(M_p);
        msp = fraction_tr_death;
        msp_r = fraction_tr_birth;
      }
      else {
        relocate_epsilon( eps, &epsilon_p );
        tp = tp_r = get_tp_epsilon_relocation( epsilon_p );
        msp = msp_r = 1 - (cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death);
      }
      //just for testing perturb theta after structure move
      for ( i=0;i<2;i++ ) {
        si = bcbn_pcg_rand()%n;
        //si = k%n;
        relocate_theta_i( theta, theta_p, n, si );
        tp = tp*compute_theta_transition_prob( theta, theta_p, n );
        tp_r = tp_r*compute_theta_transition_prob( theta_p, theta, n );
      }
      bcbn_free_lattice_children( M_p );
      M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
      //parents_dt(&M);
      children_dt(M_p);
    }

    posterior_p = compute_unnomr_posterior(epsilon_p, M_p, theta_p, D, N_u );
    //printf("posterior_p %.20Lg\n", posterior_p);


    MH_ratio = posterior_p * msp_r * tp_r / (posterior_o * msp * tp );

    alpha = MIN(1,MH_ratio);

    u = gsl_ran_flat (bcbn_RNG, 0, 1);

    if ( alpha > u ) {
      posterior_o = posterior_p;
      if( ms < cm ) {
        memcpy( theta, theta_p, sizeof(double)*n);
      }
      else {
        int h=0;
        for( i=0;i<n;i++ ) {
          for( j=0;j<n;j++ ) {
            if( M->P[i][j] != M_p->P[i][j] ) {
              if( M_p->P[i][j] ) {
                // printf("%d+%d  ", i+1, j+1);
              }
              else {
                // printf("%d-%d  ", i+1, j+1);
              }
            }
            h += M->P[i][j];
          }
        }
        if ( h != n_edges ) {
          // printf("AAAALLLLEEEERRRRTTTTT");
        }
        bcbn_free_poset( M );
        bcbn_free_lattice_children( M );
        clone_poset( M_p, M );
        M->J_P = bcbn_bfs_order_ideals(M->P, M->n, &(M->m), M->lin_ext);
        //parents_dt(M);
        children_dt(M);

        eps = epsilon_p;
        if ( ms < cm+fraction_exchange ) {
          // printf("successful poset event EXCHANGE move!!!! posterior %.10Lg number edges %d epsilon: %.3f\n", posterior_o, n_edges, eps);
        }
        else if ( ms < cm+fraction_exchange+fraction_reincarnation ) {
          // printf("successful poset REINCARNATION move!!!! posterior %.10Lg number edges %d epsilon: %.3f\n", posterior_o, n_edges, eps);
        }
        else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth ) {
          n_edges++;
          // printf("successful poset BIRTH move!!!! posterior %.10Lg number edges %d epsilon: %.3f\n", posterior_o, n_edges, eps);
        }
        else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death ) {
          n_edges--;
          // printf("successful poset DEATH move!!!! posterior %.10Lg number edges %d epsilon: %.3f\n", posterior_o, n_edges, eps);
        }
        else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death+fraction_tr_birth ) {
          n_edges = 0;
          for( i=0;i<n;i++ ) {
            for( j=0;j<n;j++ ) {
              if( M->P[i][j] ) {
                n_edges++;
              }
            }
          }
          // printf("successful transitive closure BIRTH move!!!! posterior %.10Lg number edges %d epsilon: %.3f\n", posterior_o, n_edges, eps);
        }
        else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death+fraction_tr_birth+fraction_tr_death ) {
          n_edges = 0;
          for( i=0;i<n;i++ ) {
            for( j=0;j<n;j++ ) {
              if( M->P[i][j] ) {
                n_edges++;
              }
            }
          }
          // printf("successful transitive closure DEATH move!!!! posterior %.10Lg number edges %d epsilon: %.3f\n", posterior_o, n_edges, eps);
        }
      }


      if( k>burn_in ) {
        accepted ++;
        //if ( !(accepted%record_ith) ) {
        //if ( !(k%1000) ) {
        //printf("sample number %d: posterior %.10Lg number edges %d\n", k,posterior_o, n_edges);
        //if( ms < cm ) {
        //bcbn_print_model( M );
        //printf("theta: ");
        //bcbn_print_double_array( theta, n );
        //}
        //else
        //{
        //printf("poset:\n");
        //bcbn_print_int_matrix( M->P, n, n);
        //}
        //bcbn_write_poset(0, "data/output/SYNrun1", M->P, n, k);
        //int** T = bcbn_get_int_matrix(n, n);
        //bcbn_transitive_closure( M->P, T, n );
        //bcbn_bcbn_write_poset(0, "data/output_tc/SYNrun1", T, n, k);
        //free_int_matrix( T, n );

        //}
      }
    }
    else {
      memcpy( theta_p, theta, sizeof(double)*n);

      bcbn_free_poset( M_p );
      bcbn_free_lattice_children( M_p );
      clone_poset( M, M_p );
      M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
      //parents_dt(M_p);
      children_dt(M_p);
      epsilon_p = eps;
    }

    //log samples
    if( k>burn_in && !(k%record_ith) ) {
      for (j=0; j<n; j++) {
        fprintf(samplelog, DOUBLE_FORMAT, theta[j]);
        fprintf(samplelog, "\t");
      }
      fprintf(samplelog, "posterior %.10Lg log-posterior %.5Lg number edges %d epsilon: %.3f accepted: %d, samplenr: %d u: %.5f ms: %.5f\n", posterior_o, log10l(posterior_o), n_edges, eps, accepted, k, u, ms );
      //bcbn_write_poset(0, "data/output/SYNrun8.1", M->P, n, k);
      //int** T = bcbn_get_int_matrix(n, n);
      //bcbn_transitive_closure( M->P, T, n );
      //bcbn_write_poset(0, "data/output_tc/SYNrun8.1", T, n, k);
      //free_int_matrix( T, n );
      if ( !(k%1000) ) {
        // printf("sample number %d: posterior %.10Lg number edges %d number order ideals %d\n", k,posterior_o, n_edges, M_p->m);
      }
    }

  }

  fclose(samplelog);

  free(M_p);
  // printf("poset:\n");
  bcbn_print_int_matrix( M->P, n, n);
  // printf("\n accepted: %d, ratio %f\n", accepted, (float)accepted/number_samples);
}

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
long double start_Exp_theta_MH( double eps, model* M, double* theta, data* D, int N_u, int burn_in, int number_samples, int record_ith ) {

  int i,j,k = 0;
  long double posterior_o, posterior_p, posterior_sum;
  posterior_sum = 0;
  double msp, msp_r, tp, tp_r, u;
  int accepted = 0;
  int recorded = 0;
  int si;
  int n = M->n;
  double* theta_p = bcbn_get_double_array(n);
  memcpy( theta_p, theta, sizeof(double)*n);
  long double alpha, MH_ratio;

  posterior_o = compute_unnomr_posterior(eps, M, theta, D, N_u );

  for( k=0;k<burn_in+number_samples;k++) {

    si = bcbn_pcg_rand()%n;;
    //si = k%n;
    relocate_theta_i( theta, theta_p, n, si );
    tp = compute_theta_transition_prob( theta, theta_p, n );
    tp_r = compute_theta_transition_prob( theta_p, theta, n );
    msp = 1;
    msp_r = 1;
    posterior_p = compute_unnomr_posterior(eps, M, theta_p, D, N_u );
    //printf("posterior_p %.20Lg\n", posterior_p);


    MH_ratio = posterior_p * msp_r * tp_r / (posterior_o * msp * tp );

    alpha = MIN(1,MH_ratio);

    u = gsl_ran_flat (bcbn_RNG, 0, 1);

    if ( alpha > u ) {
      memcpy( theta, theta_p, sizeof(double)*n);

      posterior_o = posterior_p;

    }
    else {
      memcpy( theta_p, theta, sizeof(double)*n);

    }
    if( k>burn_in ) {
      if ( !(k%record_ith) ) {
        recorded++;
        posterior_sum += posterior_o;
      }
    }

  }
  return posterior_sum/recorded;
}



void start_nested_MH(double eps, model* M, double* theta, data* D, int N_u, int burn_in, int number_samples, int record_ith ) {

  double fraction_exchange = 0.03;
  double fraction_reincarnation = 0.1;
  double fraction_birth = 0.5;
  double fraction_death = 0.22;
  long double alpha, MH_ratio;
  double epsilon_p;
  double u,ms;
  int n = M->n;
  double* theta_p = bcbn_get_double_array(n);
  memcpy( theta_p, theta, sizeof(double)*n);
  model* M_p = malloc( sizeof(model) );
  int i,j,k = 0;
  long double posterior_o, posterior_p;
  double msp, msp_r, tp, tp_r;
  int accepted = 0;

  int si;
  int n_edges = 0;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( M->P[i][j] ) {
        n_edges++;
      }
    }
  }

  FILE *samplelog;
  if ( (samplelog = fopen("data/output/SYNrun1sampleNest2", "w")) == NULL ) {
    // exit(1);
  }
  epsilon_p = eps;
  posterior_o = start_Exp_theta_MH(eps, M, theta, D, N_u, 20000, 1000, 20 );;
  // printf("posterior_o %.20Lg\n", posterior_o);

  clone_poset( M, M_p );
  M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
  //parents_dt(&M);
  children_dt(M_p);

  for( k=0;k<burn_in+number_samples;k++) {
    ms = gsl_ran_flat (bcbn_RNG, 0, 1);

    if ( ms < fraction_exchange ) {

      int ti,tj;
      propose_event_exchange_move( M_p, &tp, &ti, &tj );
      msp = msp_r = fraction_exchange;
      relocate_theta_i( theta, theta_p, n, ti );
      tp = compute_theta_transition_prob( theta, theta_p, n );
      tp_r = compute_theta_transition_prob( theta_p, theta, n );
      relocate_theta_i( theta, theta_p, n, tj );
      tp = tp*compute_theta_transition_prob( theta, theta_p, n );
      tp_r = tp_r*compute_theta_transition_prob( theta_p, theta, n );
      //tp = tp_r = 1;
    }
    else if ( ms < fraction_exchange+fraction_reincarnation ) {
      if ( !propose_reincarnation_move( M_p, D, N_u ) ) {
        continue;
      }
      msp = msp_r = fraction_reincarnation;
      tp = tp_r = 1;
    }
    else if ( ms < fraction_exchange+fraction_reincarnation+fraction_birth ) {
      //TODO: ugly: handle the case where no more cover relations are possible nicer!
      if( !propose_new_cover_relation(M_p, &tp, D, N_u ) ) {
        continue;
      }
      tp_r = get_tp_for_remove_cover_move(M_p);
      msp = fraction_birth;
      msp_r = 1-(fraction_exchange+fraction_reincarnation+fraction_birth);
    }
    else if ( ms < fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death ){
      if( !propose_delete_cover_relation(M_p, &tp ) ) {
        continue;
      }
      tp_r = get_tp_for_new_cover_move(M_p, D, N_u );
      msp_r = fraction_birth;
      msp = fraction_death;
    }
    else {
      relocate_epsilon( eps, &epsilon_p );
      tp = tp_r = get_tp_epsilon_relocation( epsilon_p );
      msp = msp_r = 1 - (fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death);
    }

    //just for testing perturb theta after structure move
    for ( i=0;i<3;i++ ) {
      si = bcbn_pcg_rand()%n;
      //si = k%n;
      relocate_theta_i( theta, theta_p, n, si );
      tp = tp*compute_theta_transition_prob( theta, theta_p, n );
      tp_r = tp_r*compute_theta_transition_prob( theta_p, theta, n );
    }
    bcbn_free_lattice_children( M_p );
    M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
    //parents_dt(&M);
    children_dt(M_p);

    posterior_p = start_Exp_theta_MH(epsilon_p, M_p, theta_p, D, N_u, 1000, 1000, 20 );
    //printf("posterior_p %.20Lg\n", posterior_p);


    MH_ratio = posterior_p * msp_r * tp_r / (posterior_o * msp * tp );

    alpha = MIN(1,MH_ratio);

    u = gsl_ran_flat (bcbn_RNG, 0, 1);

    if ( alpha > u ) {
      posterior_o = posterior_p;

      for( i=0;i<n;i++ ) {
        for( j=0;j<n;j++ ) {
          if( M->P[i][j] != M_p->P[i][j] ) {
            // printf("%d %d  ", i+1, j+1);
          }
        }
      }
      bcbn_free_poset( M );
      bcbn_free_lattice_children( M );
      clone_poset( M_p, M );
      M->J_P = bcbn_bfs_order_ideals(M->P, M->n, &(M->m), M->lin_ext);
      //parents_dt(M);
      children_dt(M);

      eps = epsilon_p;
      if ( ms < fraction_exchange ) {
        // printf("successful poset event EXCHANGE move!!!! posterior %.10Lg number edges %d\n", posterior_o, n_edges);
      }
      else if ( ms < fraction_exchange+fraction_reincarnation ) {
        // printf("successful poset REINCARNATION move!!!! posterior %.10Lg number edges %d\n", posterior_o, n_edges);
      }
      else if ( ms < fraction_exchange+fraction_reincarnation+fraction_birth ) {
        n_edges++;
        // printf("successful poset BIRTH move!!!! posterior %.10Lg number edges %d\n", posterior_o, n_edges);
      }
      else if ( ms < fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death ) {
        n_edges--;
        // printf("successful poset DEATH move!!!! posterior %.10Lg number edges %d\n", posterior_o, n_edges);
      }


      if( k>burn_in ) {
        accepted ++;
        if ( !(accepted%record_ith) ) {
          // printf("sample number %d: posterior %.10Lg number edges %d\n", k,posterior_o, n_edges);

          //printf("poset:\n");
          //bcbn_print_int_matrix( M->P, n, n);
          //bcbn_write_poset(0, "data/output/run16", M->P, n, k);
          //int** T = bcbn_get_int_matrix(n, n);
          //bcbn_transitive_closure( M->P, T, n );
          //bcbn_write_poset(0, "data/output_tc/run16", T, n, k);
          //free_int_matrix( T, n );

        }
      }
    }
    else {

      bcbn_free_poset( M_p );
      bcbn_free_lattice_children( M_p );
      clone_poset( M, M_p );
      M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
      //parents_dt(M_p);
      children_dt(M_p);
      epsilon_p = eps;

    }

    //log samples
    if( k>burn_in && !(k%record_ith) ) {
      fprintf(samplelog, "posterior %.10Lg log-posterior %.5Lg number edges %d epsilon: %.3f\n", posterior_o, log10l(posterior_o), n_edges, eps );
      bcbn_write_poset(0, "data/output/SYNrun1Nest2", M->P, n, k);
      int** T = bcbn_get_int_matrix(n, n);
      bcbn_transitive_closure( M->P, T, n );
      bcbn_write_poset(0, "data/output_tc/SYNrun1Nest2", T, n, k);
      free_int_matrix( T, n );
    }
  }
  free(M_p);
  // printf("poset:\n");
  bcbn_print_int_matrix( M->P, n, n);
  // printf("\n accepted: %d, ratio %f\n", accepted, (float)accepted/number_samples);

}

void run_MH_sampler( model* M, double *theta_in, double epsilon_in, data* D, int N_u, int number_samples, int thinout,
                     double **theta_matrix_out, int ***edges_cube_out, double *epsilon_out, double *log_posterior_out ) {


  /*
   * suffixes of variables in this method
   * *_p ... proposal
   * *_o ... old value
   * *_r ... reverse
   */
  double cm = 0.5;
  double fraction_exchange = 0.03;
  double fraction_reincarnation = 0.05;
  double fraction_birth = 0.2;
  double fraction_death = 0.1;
  double fraction_tr_birth = 0.01;
  double fraction_tr_death = 0.01;
  long double alpha, MH_ratio;
  double epsilon_p;
  double u,ms;
  int n = M->n;
  double* theta_p = bcbn_get_double_array(n);
  memcpy( theta_p, theta_in, sizeof(double)*n);
  model* M_p = malloc( sizeof(model) );
  int i,j,k = 0;
  long double posterior_o, posterior_p;
  double msp, msp_r, tp, tp_r;
  int accepted = 0;

  int si;
  int n_edges = 0;
  for( i=0;i<n;i++ ) {
    for( j=0;j<n;j++ ) {
      if( M->P[i][j] ) {
        n_edges++;
      }
    }
  }
  posterior_o = compute_unnomr_posterior(epsilon_in, M, theta_in, D, N_u );

  clone_poset( M, M_p );
  M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
  //parents_dt(&M);
  children_dt(M_p);
  epsilon_p = epsilon_in;

  for( k=0;k<number_samples*thinout;k++) {
    ms = gsl_ran_flat (bcbn_RNG, 0, 1);

    if( ms < cm ) {
      //relocate_whole_theta( theta, theta_p, n );
      //printf("propose relocate theta\n");
      si = bcbn_pcg_rand()%n;;
      //si = k%n;
      relocate_theta_i( theta_in, theta_p, n, si );
      tp = compute_theta_transition_prob( theta_in, theta_p, n );
      tp_r = compute_theta_transition_prob( theta_p, theta_in, n );
      msp = cm;
      msp_r = cm;
    }
    else {
      if ( ms < cm+fraction_exchange ) {
        //printf("propose event exchange\n");
        int ti,tj;
        propose_event_exchange_move( M_p, &tp, &ti, &tj );
        msp = msp_r = fraction_exchange;
        relocate_theta_i( theta_in, theta_p, n, ti );
        relocate_theta_i( theta_in, theta_p, n, tj );
        tp = compute_theta_transition_prob( theta_in, theta_p, n );
        tp_r = compute_theta_transition_prob( theta_p, theta_in, n );
        //tp = tp_r = 1;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation ) {
        //printf("propose reincarnation\n");
        if ( !propose_reincarnation_move( M_p, D, N_u ) ) {
          k--;
          continue;
        }
        msp = msp_r = fraction_reincarnation;
        tp = tp_r = 1;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth ) {
        //TODO: ugly: handle the case where no more cover relations are possible nicer!
        //printf("propose new cover\n");
        if( !propose_new_cover_relation(M_p, &tp, D, N_u ) ) {
          k--;
          continue;
        }
        tp_r = get_tp_for_remove_cover_move(M_p);
        msp = fraction_birth;
        msp_r = fraction_death;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death ) {
        //printf("propose delete cover\n");
        if( !propose_delete_cover_relation(M_p, &tp ) ) {
          k--;
          continue;
        }
        tp_r = get_tp_for_new_cover_move(M_p, D, N_u );
        msp_r = fraction_birth;
        msp = fraction_death;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death+fraction_tr_birth ) {
        //printf("propose new tc\n");
        if( !propose_new_bcbn_transitive_closure_relation(M_p, &tp, D, N_u ) ) {
          k--;
          continue;
        }
        tp_r = get_tp_for_delete_bcbn_transitive_closure_relation_move(M_p);
        msp = fraction_tr_birth;
        msp_r = fraction_tr_death;
      }
      else if ( ms < cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death+fraction_tr_birth+fraction_tr_death ) {
        //printf("propose delete tc\n");
        if( !propose_delete_bcbn_transitive_closure_relation(M_p, &tp, D, N_u ) ) {
          k--;
          continue;
        }
        tp_r = get_tp_for_new_bcbn_transitive_closure_relation_move(M_p);
        msp = fraction_tr_death;
        msp_r = fraction_tr_birth;
      }
      else {
        //printf("propose relocate epsilon\n");
        relocate_epsilon( epsilon_in, &epsilon_p );
        tp = tp_r = get_tp_epsilon_relocation( epsilon_p );
        msp = msp_r = 1 - (cm+fraction_exchange+fraction_reincarnation+fraction_birth+fraction_death);
      }
      if ( ms > cm+fraction_exchange ) {
        //just for testing perturb theta after structure move
        for ( i=0;i<2;i++ ) {
          si = bcbn_pcg_rand()%n;
          //si = k%n;
          relocate_theta_i( theta_in, theta_p, n, si );
          tp = tp*compute_theta_transition_prob( theta_in, theta_p, n );
          tp_r = tp_r*compute_theta_transition_prob( theta_p, theta_in, n );
        }
      }
      bcbn_free_lattice_children( M_p );
      M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
      //parents_dt(&M);
      children_dt(M_p);
    }

    posterior_p = compute_unnomr_posterior(epsilon_p, M_p, theta_p, D, N_u );


    long double tempor = (posterior_p-posterior_o);

    long double temporr = (long double) (log10(msp_r)+log10(tp_r)-log10(msp)-log10(tp));
    long double temp_ratio = (tempor+temporr);
    //MH_ratio = posterior_p * msp_r * tp_r / (posterior_o * msp * tp );
    long double MH_ratio = (long double) pow(10,temp_ratio);

    //MH_ratio = posterior_p * msp_r * tp_r / (posterior_o * msp * tp );

    alpha = MIN(1,MH_ratio);

    u = gsl_ran_flat (bcbn_RNG, 0, 1);

    if ( alpha > u ) {
      posterior_o = posterior_p;
      if( ms < cm ) {
        memcpy( theta_in, theta_p, sizeof(double)*n);
      }
      else {

        bcbn_free_poset( M );
        //printf("before free lattice children M\n");
        bcbn_free_lattice_children( M );
        clone_poset( M_p, M );
        M->J_P = bcbn_bfs_order_ideals(M->P, M->n, &(M->m), M->lin_ext);
        //parents_dt(M);
        children_dt(M);

        epsilon_in = epsilon_p;
      }
      accepted ++;
    }
    else {
      memcpy( theta_p, theta_in, sizeof(double)*n);

      bcbn_free_poset( M_p );
      //printf("before free lattice children M_p\n");
      bcbn_free_lattice_children( M_p );
      clone_poset( M, M_p );
      M_p->J_P = bcbn_bfs_order_ideals(M_p->P, M_p->n, &(M_p->m), M_p->lin_ext);
      //parents_dt(M_p);
      children_dt(M_p);
      epsilon_p = epsilon_in;
    }

    //log samples
    if( !(k%thinout) ) {
      memcpy( theta_matrix_out[k/thinout], theta_in, sizeof(double)*n);
      for(i=0;i<n;i++ ) {
        memcpy( edges_cube_out[k/thinout][i], M->P[i], sizeof(int)*n );
      }
      epsilon_out[k/thinout] = epsilon_in;
      //log_posterior_out[k/thinout] = (double)log10l(posterior_o);
      log_posterior_out[k/thinout] = (double) posterior_o;
      //printf("posterior %.10Lg log-posterior %.5Lg epsilon: %.3f accepted: %d, samplenr: %d u: %.5f ms: %.5f\n", posterior_o, log10l(posterior_o), epsilon_in, accepted, k, u, ms );
      /*n_edges = 0;
       for( i=0;i<n;i++ ) {
       for( j=0;j<n;j++ ) {
       if( M->P[i][j] ) {
       n_edges++;
       }
       }
       }
       fprintf(samplelog, "posterior %.10Lg log-posterior %.5Lg number edges %d epsilon: %.3f accepted: %d, samplenr: %d u: %.5f ms: %.5f\n", posterior_o, log10l(posterior_o), n_edges, epsilon_in, accepted, k, u, ms );
       */
    }

  }
  bcbn_free_poset( M_p );
  bcbn_free_lattice_children( M_p );
  free(M_p);
}


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

void sample_full_cbn(double *theta_in, int *nevents, double *epsilon_in, int * edges_in, int *number_samples, int *thinout, int *patdata,
                     int *number_cases, double *theta_out, double *epsilon_out, int* edges_out, double* log_posterior_out ) {

  unsigned int seed = (unsigned) time(NULL);  // r, random seed
  bcbn_pcg_srand(seed);
  bcbn_RNG = gsl_rng_alloc (gsl_rng_taus);  // global variable
  gsl_rng_set(bcbn_RNG, seed);  // seed rng

  int i,j,k,n;
  n=(*nevents);
  int n_samples = (*number_samples);
  int n_cases = (*number_cases);


  model M;
  M.n = n;
  M.P = bcbn_get_int_matrix(n, n);

  for( i=0;i<n;i++){
    for( j=0;j<n;j++) {
      M.P[i][j] = edges_in[i*n+j];
    }
  }
  bcbn_precompute_binary(M.n);

  M.lin_ext = bcbn_get_int_array(M.n);  // a linear extension of the poset
  int** pat = bcbn_get_int_matrix( n_cases,n );
  for( i=0;i<n_cases;i++) {
    for( j=0;j<n;j++ ) {
      pat[i][j] = patdata[i*n+j];
    }
  }

  int* pat_idx = bcbn_get_int_array(n_cases);
  int N_u;
  data* D = bcbn_make_data_set(pat, n_cases, M.n, &N_u, pat_idx);

  for (k=0; k<n_cases; k++)
    free(pat[k]);
  free(pat);

  double ** theta_matrix_out = bcbn_get_double_matrix( n_samples ,n );
  int *** edges_cube_out = get_int_cube( n_samples, n, n );
  for( i=0;i<n_samples;i++ ) {
    for( j=0;j<n;j++ ) {
      for( k=0;k<n;k++ ) {
        edges_cube_out[i][j][k] = 0;
      }
    }
  }

  M.J_P = bcbn_bfs_order_ideals(M.P, M.n, &(M.m), M.lin_ext);
  //parents_dt(&M);
  children_dt(&M);
  run_MH_sampler( &M, theta_in, (*epsilon_in), D, N_u, n_samples, (*thinout), theta_matrix_out, edges_cube_out, epsilon_out, log_posterior_out );

  for( i=0;i<n_samples;i++) {
    for( j=0;j<n;j++ ) {
      theta_out[i*n+j] = theta_matrix_out[i][j];
    }
  }
  for( i=0;i<n_samples;i++ ) {
    for( j=0;j<n;j++ ) {
      for( k=0;k<n;k++ ) {
        //printf("%d ",i*n*n+j*n+k);
        edges_out[i*n*n+j*n+k] = edges_cube_out[i][j][k];
        //printf(".");
        //printf("%d ", edges_cube_out[i][j][k]);
      }
    }
  }
  free_double_matrix( theta_matrix_out, n_samples );
  free_int_cube( edges_cube_out, n_samples, n );
  //bcbn_free_data( D, N_u, n);
}
