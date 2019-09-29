#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

#include <stdlib.h>
//#ifndef WINDOWS
//#include <stdio.h>
//#else
#include "pfwstdio.h"
//#endif
#include <math.h>
#ifdef ANSIPROTO
#include <float.h>
#endif

#include "constant.h"
#include "param.h"
#include "sparse.h"

#ifdef ANSIPROTO
void InitializeWorkingRow(void);
void LoadUnloadWorkingRow1(SparseMatrixElement *Ptr2,BOOLEAN SetReset);
void Normalize(INDEX I);
int DoFactorization(void);
int factor(SparseMatrix *Mptr);
#else
void InitializeWorkingRow();
void LoadUnloadWorkingRow1();
void Normalize();
int DoFactorization();
int factor();
#endif

/* Global variable definitions */
  SparseMatrix *Matrix2;
  INDEX Nstop2;
  LONGINT Nmult2;
  ELEMENTVALUETYPE DiagVal2;
  ELEMENTVALUETYPE *Swr2;
  SparseMatrixElement **DiagPtr2;


/* ============================ InitializeWorkingRow ==================== */
    void InitializeWorkingRow()
    {
      INDEX I,J;
      SparseMatrixElement *Ptr1;
      /* BEGIN */
#ifdef WINDOWS
     Swr2 = new ELEMENTVALUETYPE[Matrix2->n1+1];
     DiagPtr2 = new SparseMatrixElement*[Matrix2->n1+1];
#else
     Swr2 = (ELEMENTVALUETYPE *) calloc((Matrix2->n1+1),sizeof(ELEMENTVALUETYPE));
     DiagPtr2 = (struct SparseMatrixElement **)
             calloc((Matrix2->n1+1),sizeof(struct SparseMatrixElement *));
#endif
     for (I=1; I<=Matrix2->n1; I++) {
      Ptr1 = Matrix2->RowHead[I];
      while (Ptr1 != NULL) {
        J = Ptr1->Col;
        if (I == J) DiagPtr2[I] = Ptr1;
        Ptr1 = Ptr1->RowNext;
      }
     }
     for (I=1; I<=Matrix2->n1; I++) {
      if (DiagPtr2[I] == NULL) {
        fCustomPrint(stderr,"\nError: Diagonal element missing in the Jacobian.\n");
        stopExecute(ERROREXIT);
      }
     }
    }

/* ========================== LoadUnloadWorkingRow ====================== */
#ifdef ANSIPROTO
    void LoadUnloadWorkingRow1(SparseMatrixElement *Ptr2,BOOLEAN SetReset)
#else
    void LoadUnloadWorkingRow1(Ptr2,SetReset)
    SparseMatrixElement *Ptr2;
    BOOLEAN SetReset;
#endif
    {
      /* BEGIN */
      while (Ptr2 != NULL) {
        if (SetReset) {
        Swr2[Ptr2->Col] = Ptr2->Value;
        } else {
        Ptr2->Value = Swr2[Ptr2->Col];
        Swr2[Ptr2->Col] = 0.0;
        }
        Ptr2 = Ptr2->RowNext;
      }
    }

/* ========================= Normalize ================================ */
#ifdef ANSIPROTO
    void Normalize(INDEX I)
#else
    void Normalize(I)
    INDEX I;
#endif
    {
      struct SparseMatrixElement *Ptr1;
      INDEX J;
      /* BEGIN */
      Ptr1 = Matrix2->RowHead[I];
      while (Ptr1 != NULL) {
        J = Ptr1->Col;
        if (J == I) {
        Ptr1->Value = DiagVal2;
        /* DiagPtr2[I] = Ptr1; */
        } else if (Ptr1->Col > I) {
        Ptr1->Value = DiagVal2 * Ptr1->Value;
        Nmult2++;
        }
        Ptr1 = Ptr1->RowNext;
      }
    }

/* ========================== DoFactorization ========================== */
int DoFactorization()
  {
    INDEX I,J,K;
    struct SparseMatrixElement *Ptr1, *Ptr2;
    /* BEGIN {DoFactorization} */
    InitializeWorkingRow();
    for (I=1; I<=Matrix2->n1; I++) {
      Ptr1 = Matrix2->RowHead[I];
      LoadUnloadWorkingRow1(Ptr1,TRUE);
      while (Ptr1 != NULL) {
        J = Ptr1->Col;
      if ((J < I) && (J <= Nstop2)) {
        Ptr2 = Matrix2->RowHead[J];
          while (Ptr2 != NULL) {
            K = Ptr2->Col;
            if (K > J) {
            Swr2[K] = Swr2[K] - Swr2[J]*Ptr2->Value;
              /*
            MultiplyValues(TempVal,Swr2[J],Ptr2^.Value);
            SubtractValues(Swr2[K],Swr2[K],TempVal);
              */
            }
            Ptr2 = Ptr2->RowNext;
          }
        Swr2[J] = Swr2[J] * DiagPtr2[J]->Value;
        Nmult2++;
          /*
        MultiplyValues(Swr2[J],Swr2[J],DiagPtr2[J]^.Value);
          */
      } else if (J == I) DiagPtr2[I] = Ptr1;
        Ptr1 = Ptr1->RowNext;
      }
      LoadUnloadWorkingRow1(Matrix2->RowHead[I],FALSE);
      DiagVal2 = DiagPtr2[I]->Value;
      /*
      EquateValues(DiagVal2,DiagPtr2[I]^.Value);
      */
      if (I <= Nstop2) {
      if (NearZero(DiagVal2)) return(WARNINGEXIT);
      DiagVal2 = 1.0 / DiagVal2;
      Normalize(I);
      }
    }
    return(0);
  } /* END {DoFactorization} */

/* ============================== factor ================================= */
#ifdef ANSIPROTO
int factor(SparseMatrix *Mptr)
#else
int factor(Mptr)
SparseMatrix *Mptr;
#endif
{
  /* BEGIN FactorMatrix2 */
  Matrix2 = Mptr;
  Nmult2 = 0;
  Nstop2 = Matrix2->n1; 
  if (DoFactorization()==WARNINGEXIT) return(WARNINGEXIT);
#ifdef WINDOWS
  delete[] Swr2;
  delete[] DiagPtr2;
#else
  free(Swr2);
  free(DiagPtr2);
#endif
  return(0);
}
