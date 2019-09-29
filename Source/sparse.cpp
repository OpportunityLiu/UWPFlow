#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

/* sparse.c 030590 */

#include <stdlib.h>
//#ifndef WINDOWS
//#include <stdio.h>
//#else
#include "pfwstdio.h"
//#endif
#include <math.h>
#include <string.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"

/* ------- Global Variables ------ */
BOOLEAN InputError;

/* ===================== ErrorHalt ================================ */
#ifdef ANSIPROTO
void ErrorHalt(char *Msg)
#else
void ErrorHalt(Msg)
char *Msg;
#endif
{
  fCustomPrint(stderr,"ERROR: %s\n",Msg);
  InputError=TRUE;
  /* stopExecute(ERROREXIT) */
}


/* ============================== TransposeMatrix ==================== */
#ifdef ANSIPROTO
  void TransposeMatrix(SparseMatrix *Matrix)
#else
  void TransposeMatrix(Matrix)
  SparseMatrix *Matrix;
#endif
  {
    INDEX I,RowI,ColJ,N0;
    SparseMatrixElement **TempHead; /* array[1..maxN] */
    SparseMatrixElement *Ptr1,*PtrR,*PtrC;

    /* BEGIN */
#ifdef WINDOWS
    TempHead = new SparseMatrixElement*[Matrix->n1+1];
#else
    TempHead = (SparseMatrixElement **)
               calloc(Matrix->n1+1,sizeof(SparseMatrixElement *));
#endif
    for (I=1; I<=Matrix->n1; I++) {
      TempHead[I] = Matrix->RowHead[I];
      Ptr1 = Matrix->RowHead[I];
      while (Ptr1 != NULL) {
        /* WITH Ptr1^ DO BEGIN */
        PtrR = Ptr1->RowNext;
        PtrC = Ptr1->ColNext;
        Ptr1->RowNext = PtrC;
        Ptr1->ColNext = PtrR;
        RowI = Ptr1->Row;
        ColJ = Ptr1->Col;
        Ptr1->Row = ColJ;
        Ptr1->Col = RowI;
        Ptr1 = PtrR;
      }
    }
    for (I=1; I<=Matrix->n2; I++) Matrix->RowHead[I] = Matrix->ColHead[I];
    for (I=1; I<=Matrix->n2; I++) Matrix->ColHead[I] = TempHead[I];
    N0 = Matrix->n1;
    Matrix->n1 = Matrix->n2;
    Matrix->n2 = N0;
#ifdef WINDOWS
    delete[] TempHead;
#else
    free(TempHead);
#endif
  }


/* ===================== AllocatePermutation ========================= */
#ifdef ANSIPROTO
struct IntegerVector * AllocatePermutation(INDEX Size)
#else
struct IntegerVector * AllocatePermutation(Size)
INDEX Size;
#endif
{
  struct IntegerVector * Vector;
  int i;

  /* BEGIN */
#ifdef WINDOWS
  Vector = new IntegerVector;
  Vector->p = new INDEX[Size+1];
#else
  Vector = (IntegerVector *) malloc(sizeof(IntegerVector));
  if (Vector == NULL) {ErrorHalt("Insufficient memory to allocate integer vector"); stopExecute(ERROREXIT);}
  Vector->p = (INDEX *) calloc((Size+1),sizeof(INDEX));
  if (Vector->p == NULL) {ErrorHalt("Insufficient memory to allocate permutation vector"); stopExecute(ERROREXIT);}
#endif
  Vector->N = Size;
  for(i=0;i<Size+1;i++) Vector->p[i]=0;
  return(Vector);
}


/* ===================== SortRowsColumns ======================= */
#ifdef ANSIPROTO
void SortRowsColumns(SparseMatrix *Matrix)
#else
void SortRowsColumns(Matrix)
SparseMatrix *Matrix;
#endif
{
  SparseMatrixElement *TempP1;
  SparseMatrixElement *TempP2;
  int I;
  /* BEGIN */
  for (I=1; I<=Matrix->n2; I++) Matrix->ColHead[I] = NULL;
  I = Matrix->n1;
  while (I>0) {
    TempP1 = Matrix->RowHead[I];
    while (TempP1 != NULL) {
      TempP2 = TempP1->RowNext;
      TempP1->ColNext = Matrix->ColHead[TempP1->Col];
      Matrix->ColHead[TempP1->Col] = TempP1;
      TempP1 = TempP2;
    } /* END WHILE */
    I--;
  } /* END WHILE */
  for (I=1; I<=Matrix->n1; I++) Matrix->RowHead[I] = NULL;
  I = Matrix->n2;
  while (I > 0) {
    TempP1 = Matrix->ColHead[I];
    while (TempP1 != NULL) {
      TempP2 = TempP1->ColNext;
      TempP1->RowNext = Matrix->RowHead[TempP1->Row];
      Matrix->RowHead[TempP1->Row] = TempP1;
      TempP1 = TempP2;
    } /* END WHILE */
    I--;
  } /* END WHILE */
} /* END RowSort */


/* ============================ NearZero ===========================*/
#ifdef ANSIPROTO
BOOLEAN NearZero(ELEMENTVALUETYPE Value)
#else
BOOLEAN NearZero(Value)
ELEMENTVALUETYPE Value;
#endif
{
  return(fabs(Value) < SINGULARITYZERO);
}

