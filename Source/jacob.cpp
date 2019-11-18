/* Initialize variables. */

#include "jacob.h"

/* ------- Global Variables ------ */
IntegerVector *RowPer, *ColPer;

/* -------------------- JacElement ---------------------- */
void JacElement(SparseMatrix *Mptr, INDEX I, INDEX J, VALUETYPE val)
/* Jacobian Element */
{
  SparseMatrixElement *Ptr;
  INDEX i, j;

  /*  if (val==0) return;*/
  i = I;
  j = J;
  if (RowPer->p[I] != 0)
    I = RowPer->p[I];
  if (ColPer->p[J] != 0)
    J = ColPer->p[J];
  for (Ptr = Mptr->RowHead[I]; Ptr != nullptr; Ptr = Ptr->RowNext)
    if (Ptr->Row == I && Ptr->Col == J) {
      if (flagReducedContinuation && (DxZero[i] || DxZero[j])) {
        if (i == j)
          Ptr->Value = 1.;
        else
          Ptr->Value = 0;
      } else
        Ptr->Value = Ptr->Value + val;
      break;
    }
  if (Ptr == nullptr &&
      (i == j || !(flagReducedContinuation && (DxZero[i] || DxZero[j])))) {
    Ptr = new SparseMatrixElement;
    Ptr->Row = I;
    Ptr->Col = J;
    Ptr->RowNext = Mptr->RowHead[I];
    Mptr->RowHead[I] = Ptr;
    Ptr->ColNext = Mptr->ColHead[J];
    Mptr->ColHead[J] = Ptr;
    if (flagReducedContinuation && (DxZero[i] || DxZero[j]))
      Ptr->Value = 1.;
    else
      Ptr->Value = val;
  }
}

/* ------------------- Jacobian --------------------- */
void Jacobian(void)
/* Allocate memory for the Jacobian. */
{
  ACbusData *ACptr;
  INDEX i, N, N1, N2;

  Dx = x0 = x0p = dx = dF = nullptr;
  DxZero = nullptr;
  Vlist = Vlistp = nullptr;
  ACvar = new INDEX[Nac + 1];
  N = 1;
  for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    ACvar[ACptr->N] = N;
    if (PQcont)
      N1 = ACptr->Ncont;
    else
      N1 = 0;
    if (Acont && strpbrk(ACptr->Type, "A"))
      N = N + 3 + N1;
    else
      N = N + 2 + N1;
    if (ACptr->Gen != nullptr) {
      ACptr->Gen->Nvar = N - 1;
      N = N + 11;
    }
  }
  NacVar = N - 1;
  N1 = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar + 7 * Nstatcom; /* FACTS  */
  RowPartition = new IntegerVector;
  ColPartition = new IntegerVector;
  if (flagH) {
    N1++;
    N2 = N1;
    RowPartition->N = 2;
    ColPartition->N = 1;
    RowPartition->p = new INDEX[3];
    ColPartition->p = new INDEX[3];
    RowPartition->p[1] = N1 - 1;
    RowPartition->p[2] = N1;
    ColPartition->p[1] = N1;
    RowPartition->p[0] = ColPartition->p[0] = ColPartition->p[2] = 0;
    Dx = new VALUETYPE[N1 + 1];
    x0 = new VALUETYPE[N1 + 1];
    x0p = new VALUETYPE[N1 + 1];
    for (i = 0; i < N1 + 1; i++)
      Dx[i] = x0[i] = x0p[i] = 0.;
  } else if (flagPoC) {
    N2 = N1;
    N1 = 2 * N1 + 1;
    RowPartition->N = ColPartition->N = 1;
    RowPartition->p = new INDEX[2];
    ColPartition->p = new INDEX[2];
    RowPartition->p[1] = ColPartition->p[1] = N2;
    RowPartition->p[0] = ColPartition->p[0] = 0;
    Dx = new VALUETYPE[N2 + 1];
    x0 = new VALUETYPE[N2 + 1];
    for (i = 0; i < N2 + 1; i++)
      Dx[i] = x0[i] = 0.;
  } else {
    N2 = N1;
    RowPartition->N = ColPartition->N = 1;
    RowPartition->p = new INDEX[2];
    ColPartition->p = new INDEX[2];
    RowPartition->p[1] = ColPartition->p[1] = N2;
    RowPartition->p[0] = ColPartition->p[0] = 0;
  }
  dx = new VALUETYPE[N1 + 2];
  dF = new VALUETYPE[N1 + 2];
  Jac = new SparseMatrix;
  for (i = 0; i < N1 + 2; i++)
    dx[i] = dF[i] = 0.;
  Jac->n1 = N2;
  Jac->n2 = N2;
  Jac->RowHead = new SparseMatrixElement *[N1 + 2];
  Jac->ColHead = new SparseMatrixElement *[N1 + 2];
  for (i = 0; i < N1 + 2; i++)
    Jac->RowHead[i] = Jac->ColHead[i] = nullptr;
  NewRow = AllocatePermutation(N1);
  OldRow = AllocatePermutation(N1);
  NewCol = AllocatePermutation(N1);
  OldCol = AllocatePermutation(N1);
  NewRow->N = NewCol->N = N2;
  OldRow->N = OldCol->N = N2;
}
