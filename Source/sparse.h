#pragma once

/* sparse.h 030990 */

#include "constant.h"

struct IntegerVector {
  INDEX N;
  INDEX *p;
};

struct SparseMatrixElement {
  INDEX Row;
  INDEX Col;
  ELEMENTVALUETYPE Value;
  SparseMatrixElement *RowNext;
  SparseMatrixElement *ColNext;
};

struct SparseMatrix {
  INDEX n1;
  INDEX n2;
  SparseMatrixElement **RowHead;
  SparseMatrixElement **ColHead;
};

void ErrorHalt(char *Msg);
void TransposeMatrix(SparseMatrix *Matrix);
struct IntegerVector *AllocatePermutation(INDEX Size);
void SortRowsColumns(SparseMatrix *Matrix);
bool NearZero(ELEMENTVALUETYPE Value);
