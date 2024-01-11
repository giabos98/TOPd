/* -------------------------------------------------------------------- */
/*      Example program to show the use of the "PARDISO" routine        */
/*      for evaluating the Schur-complement                             */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Drosos Kourounis, Institute of Computational Science            */
/*      Universita della Svizzera Italiana.                             */
/*      Email: drosos.kourounis@usi.ch                                  */
/* -------------------------------------------------------------------- */

// C++ compatible

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

using namespace std;

#include "pardiso.h"

#if 0
/* PARDISO prototype. */
extern "C" void pardisoinit_d(void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso_d(void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix_d(int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec_d(int *, int *, double *, int *);
extern "C" void pardiso_printstats_d(int *, int *, double *, int *, int *, int *, double *, int *);
extern "C" void pardiso_get_schur_d(void*, int*, int*, int*, double*, int*, int*);
#endif

void dumpCSR(const char* filename, int n, int* ia, int* ja, double* a)
{
  fstream fout(filename, ios::out);
  fout << n << endl;
  fout << n << endl;
  fout << ia[n] << endl;
  
  for (int i = 0; i <= n; i++)
  {
    fout << ia[i] << endl;
  }

  for (int i = 0; i < ia[n]; i++)
  {
    fout << ja[i] << endl;
  }

  for (int i = 0; i < ia[n]; i++)
  {
    fout << a[i] << endl;
  }

  fout.close();
}


void printCSR(int n, int nnz, int* ia, int* ja, double* a)
{
  cout << "rows: " << setw(10) << n << endl;
  cout << "nnz : " << setw(10) << nnz << endl;

  if (nnz == n*n)
  {
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        cout << setw(10) << a[i*n + j];
      }
      cout << endl;
    }
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      for (int index = ia[i]; index < ia[i+1]; index++)
      {
        int j = ja[index];
        cout << setw(10) << "(" << i << ", " << j << ") " << a[index];
      }
      cout << endl;
    }
  }
}


void shiftIndices(int n, int nonzeros, int* ia, int* ja, int value)
{
  int i;
    for (i = 0; i < n+1; i++) 
    {
        ia[i] += value;
    }
    for (i = 0; i < nonzeros; i++) 
    {
        ja[i] += value;
    }
}




int main( void ) 
{
    /* Matrix data.
     */

    int    n =10;
    int    ia[11] = { 0, 6, 11, 15, 19, 22, 27, 31, 36, 45, 54};

    int    ja[54] = { 0,    2,       5, 6,    8, 9,
                         1, 2,    4,          8, 9,
                            2,             7, 8, 9, 
                               3,       6,    8, 9,
                         1,                   8, 9,
                            2,       5,    7, 8, 9,
                         1,             6,    8, 9,
                            2,          6, 7, 8, 9,
                      0, 1, 2, 3, 4, 5, 6, 7, 8,
                      0, 1, 2, 3, 4, 5, 6, 7,    9};
   
    double  a[54] = { 7.0,      1.0,           2.0, 7.0,     17.0, 8.5,
                          -4.0, 8.0,      2.0,                6.0, 3.0,
                                1.0,                     5.0, 6.0, 3.0,
                                     7.0,           9.0,     16.0, 8.0,
                          -4.0,                              -4.0,-2.0,
                                7.0,           3.0,      8.0,18.0, 9.0,
                           1.0,                    11.0,     12.0, 6.0,
                               -3.0,                2.0, 5.0, 4.0, 2.0,
                 17.0, 6.0, 6.0, 16.0, -4.0, 18.0, 12.0, 4.0, 0.0,
                  8.5, 3.0, 3.0,  8.0, -2.0,  9.0,  6.0, 2.0,      0.0};



    int      nnz = ia[n];
    int      mtype = 11;        /* Real unsymmetric matrix */

    dumpCSR("matrix_sample.csr", n, ia, ja, a);
    printCSR(n, nnz, ia, ja, a);

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[65];
    double   dparm[64];
    int      solver;
    int      maxfct, mnum, phase, error, msglvl;

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */

/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters and initialize the solvers      */
/*     internal adress pointers. This is only necessary for the FIRST   */
/*     call of the PARDISO solver.                                      */
/* ---------------------------------------------------------------------*/
      
    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit_d(pt,  &mtype, &solver, &iparm[1], dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         return 1;
    }
    else
        printf("[PARDISO]: License check was successful ... \n");
 

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else 
    {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[3]  = num_procs;
    iparm[11] = 1;
    iparm[13] = 0;
   
    
    maxfct = 1;         /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */


/* -------------------------------------------------------------------- */    
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */ 
    shiftIndices(n, nnz, ia, ja, 1);


/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
    pardiso_chkmatrix_d(&mtype, &n, a, ia, ja, &error);
    if (error != 0) 
    {
        printf("\nERROR in consistency of matrix: %d", error);
        exit(1);
    }

 
/* -------------------------------------------------------------------- */    
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */ 
    
    int nrows_S = 2;
    phase       = 12;
    iparm[38]   = nrows_S; 

    int nb = 0;

    pardiso_d(pt, &maxfct, &mnum, &mtype, &phase,
               &n, a, ia, ja, &idum, &nb,
               &iparm[1], &msglvl, &ddum, &ddum, &error,  dparm);
    
    if (error != 0) 
    {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    printf("\nReordering completed ...\n");
    printf("Number of nonzeros in factors  = %d\n", iparm[18]);
    printf("Number of factorization GFLOPS = %d\n", iparm[19]);
    printf("Number of nonzeros is   S      = %d\n", iparm[39]);
   
/* -------------------------------------------------------------------- */    
/* ..  allocate memory for the Schur-complement and copy it there.      */
/* -------------------------------------------------------------------- */ 
    int nonzeros_S = iparm[39];
 
    int* iS        = new int[nrows_S+1];
    int* jS        = new int[nonzeros_S];
    double* S      = new double[nonzeros_S];
  
    pardiso_get_schur_d(pt, &maxfct, &mnum, &mtype, S, iS, jS);

    printCSR(nrows_S, nonzeros_S, iS, jS, S);

/* -------------------------------------------------------------------- */    
/* ..  Convert matrix from 1-based Fortan notation to 0-based C         */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */ 
    shiftIndices(n, nnz, ia, ja, -1);
    
/* -------------------------------------------------------------------- */    
/* ..  Convert Schur complement from Fortran notation held internally   */
/*     to 0-based C notation                                            */
/* -------------------------------------------------------------------- */ 
    shiftIndices(nrows_S, nonzeros_S, iS, jS, -1);


/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */ 
    phase = -1;                 /* Release internal memory. */

    pardiso_d(pt, &maxfct, &mnum, &mtype, &phase,
            &n, &ddum, ia, ja, &idum, &idum,
            &iparm[1], &msglvl, &ddum, &ddum, &error,  dparm);


    delete[] iS;
    delete[] jS;
    delete[] S;

    printf ("EXIT: Completed\n");

    return 0;
} 



