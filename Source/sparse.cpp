/* sparse.c 030590 */

#include "sparse.h"
#include "constant.h"
#include "param.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------- Global Variables ------ */
bool InputError;

/* ===================== ErrorHalt ================================ */
void ErrorHalt(char *Msg) {
  fprintf(stderr, "ERROR: %s\n", Msg);
  InputError = true;
  /* exit(ERROREXIT) */
}

/* ============================== TransposeMatrix ==================== */
void TransposeMatrix(SparseMatrix *Matrix) {
  INDEX I, RowI, ColJ, N0;
  SparseMatrixElement **TempHead; /* array[1..maxN] */
  SparseMatrixElement *Ptr1, *PtrR, *PtrC;

  /* BEGIN */
  TempHead = new SparseMatrixElement *[Matrix->n1 + 1];
  for (I = 1; I <= Matrix->n1; I++) {
    TempHead[I] = Matrix->RowHead[I];
    Ptr1 = Matrix->RowHead[I];
    while (Ptr1 != nullptr) {
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
  for (I = 1; I <= Matrix->n2; I++)
    Matrix->RowHead[I] = Matrix->ColHead[I];
  for (I = 1; I <= Matrix->n2; I++)
    Matrix->ColHead[I] = TempHead[I];
  N0 = Matrix->n1;
  Matrix->n1 = Matrix->n2;
  Matrix->n2 = N0;
  delete[] TempHead;
}

/* ===================== AllocatePermutation ========================= */
struct IntegerVector *AllocatePermutation(INDEX Size) {
  struct IntegerVector *Vector;
  int i;

  /* BEGIN */
  Vector = new IntegerVector;
  Vector->p = new INDEX[Size + 1];
  Vector->N = Size;
  for (i = 0; i < Size + 1; i++)
    Vector->p[i] = 0;
  return (Vector);
}

/* ===================== SortRowsColumns ======================= */
void SortRowsColumns(SparseMatrix *Matrix) {
  SparseMatrixElement *TempP1;
  SparseMatrixElement *TempP2;
  int I;
  /* BEGIN */
  for (I = 1; I <= Matrix->n2; I++)
    Matrix->ColHead[I] = nullptr;
  I = Matrix->n1;
  while (I > 0) {
    TempP1 = Matrix->RowHead[I];
    while (TempP1 != nullptr) {
      TempP2 = TempP1->RowNext;
      TempP1->ColNext = Matrix->ColHead[TempP1->Col];
      Matrix->ColHead[TempP1->Col] = TempP1;
      TempP1 = TempP2;
    } /* END WHILE */
    I--;
  } /* END WHILE */
  for (I = 1; I <= Matrix->n1; I++)
    Matrix->RowHead[I] = nullptr;
  I = Matrix->n2;
  while (I > 0) {
    TempP1 = Matrix->ColHead[I];
    while (TempP1 != nullptr) {
      TempP2 = TempP1->ColNext;
      TempP1->RowNext = Matrix->RowHead[TempP1->Row];
      Matrix->RowHead[TempP1->Row] = TempP1;
      TempP1 = TempP2;
    } /* END WHILE */
    I--;
  } /* END WHILE */
} /* END RowSort */

/* ============================ NearZero ===========================*/
bool NearZero(ELEMENTVALUETYPE Value) {
  return (fabs(Value) < SINGULARITYZERO);
}
