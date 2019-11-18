#pragma once

/* sparse.h 030990 */

#include "constant.h"

typedef struct IntegerVector {
  INDEX N;
  INDEX *p;
} IntegerVector;

typedef struct SparseMatrixElement {
  INDEX Row;
  INDEX Col;
  ELEMENTVALUETYPE Value;
  struct SparseMatrixElement *RowNext;
  struct SparseMatrixElement *ColNext;
} SparseMatrixElement;

typedef struct /* SparseMatrix */ {
  INDEX n1;
  INDEX n2;
  SparseMatrixElement **RowHead;
  SparseMatrixElement **ColHead;
} SparseMatrix;


#ifdef ANSIPROTO
void ErrorHalt(char *Msg);
void	TransposeMatrix(SparseMatrix *Matrix);
struct IntegerVector *AllocatePermutation(INDEX Size);
void SortRowsColumns(SparseMatrix *Matrix);
bool NearZero(ELEMENTVALUETYPE Value);
#else
void ErrorHalt();
void	TransposeMatrix();
struct IntegerVector *AllocatePermutation();
void SortRowsColumns();
bool NearZero();
#endif
