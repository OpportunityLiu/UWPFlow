/*  Perform Repeat Solution.
    Remarks: Matrix1 must be already factored (see FACTOR or FACTORNS).
             The input vector is not permuted, this is done internally
             in this routine using the row permutation vector created
             by the ordering routine (see FACTORNS).                  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"

#ifdef ANSIPROTO
void ForwardSubstitution(void);
void DiagonalScaling(void);
void BackSubstitution(void);
void CreateDiagonalPointer(void);
void repsolp(SparseMatrix *Mptr,VALUETYPE *Vptr,
             IntegerVector *PermR,IntegerVector *PermC);
#else
void ForwardSubstitution();
void DiagonalScaling();
void BackSubstitution();
void CreateDiagonalPointer();
void repsolp();
#endif


/* ==================== Global definitions ============================= */
SparseMatrix *Matrix1;
INDEX Nstop1;
LONGINT Nmult1;
SparseMatrixElement **DiagPtr;
VALUETYPE *FullVector;
int DetSign;



/* ======================== ForwardSubstition =========================== */
  void ForwardSubstitution()
  {
    INDEX I,J;
    SparseMatrixElement *Ptr1;
    /* BEGIN ForwardSubstitution */
    I = 1;
    while (I <= Matrix1->n1) {
      Ptr1 = Matrix1->RowHead[I];
      while (Ptr1 != nullptr) {
        J = Ptr1->Col;
        if ((J < I) && (J <= Nstop1)) {
          FullVector[I] = FullVector[I] - FullVector[J] * Ptr1->Value;
          Nmult1++;
        }
        Ptr1 = Ptr1->RowNext;
      }
      I++;
    }
  } /* END ForwardSubstitution */

/* =========================== DiagonalScaling =========================== */
  void DiagonalScaling()
  {
    INDEX I;
    /* BEGIN DiagonalScaling */
    DetSign=1;
    for (I=1; I<=Nstop1; I++) {
      FullVector[I] = FullVector[I] * DiagPtr[I]->Value;
      if (DiagPtr[I]->Value<0) DetSign=-DetSign;
      Nmult1++;
    }
  } /* END DiagonalScaling */

/* =============================== BackSubstitution ==================== */
  void BackSubstitution()
  {
    INDEX I,J;
    SparseMatrixElement *Ptr1;
    /* BEGIN BackSubstitution */
    I = Nstop1;
    while (I > 0) {
      Ptr1 = Matrix1->RowHead[I];
      while (Ptr1 != nullptr) {
        J = Ptr1->Col;
        if (J > I) {
          FullVector[I] = FullVector[I] - FullVector[J] * Ptr1->Value;
          Nmult1++;
        }
        Ptr1 = Ptr1->RowNext;
      }
      I--;
    }
  } /* END BackSubstitution */

/* ========================= CreateDiagonalPointer ======================= */
  void CreateDiagonalPointer()
  {
    INDEX i;
    SparseMatrixElement *Ptr1;

    /* BEGIN */
#ifdef WINDOWS
    DiagPtr = new SparseMatrixElement*[Matrix1->n1+1];
#else
    DiagPtr = (SparseMatrixElement **)
              calloc((Matrix1->n1+1),sizeof(SparseMatrixElement *));
#endif
    for(i=0;i<Matrix1->n1+1;i++) DiagPtr[i]=nullptr;
    for (i=1; i<=Matrix1->n1; i++) {
      Ptr1 = Matrix1->RowHead[i];
      while ((Ptr1 != nullptr) && (DiagPtr[i] == nullptr)) {
        if (Ptr1->Col == Ptr1->Row) DiagPtr[i] = Ptr1;
        Ptr1 = Ptr1->RowNext;
      }
     }
  }


/* =========================== repsolp ================================== */
#ifdef ANSIPROTO
void repsolp(SparseMatrix *Mptr,VALUETYPE *Vptr,
             IntegerVector *PermR,IntegerVector *PermC)
#else
void repsolp(Mptr,Vptr,PermR,PermC)
SparseMatrix *Mptr;
VALUETYPE *Vptr;
IntegerVector *PermR,*PermC;
#endif
{
  INDEX i;
  /* BEGIN RepeatSolution */
  Matrix1 = Mptr;
  CreateDiagonalPointer();
  Nstop1 = Matrix1->n1;
  Nmult1 = 0;
#ifdef WINDOWS
  FullVector= new VALUETYPE[Nstop1+1];
#else
  FullVector=(VALUETYPE *) malloc((Nstop1+1)*sizeof(VALUETYPE));
  if (FullVector==nullptr) {ErrorHalt("Insufficient memory for solution vector"); exit(ERROREXIT);}
#endif
  for (i=1;i<=Nstop1;i++) FullVector[i]=Vptr[PermR->p[i]];
  ForwardSubstitution();
  DiagonalScaling();
  BackSubstitution();
  for (i=1;i<=Nstop1;i++) Vptr[i]=FullVector[PermC->p[i]];
#ifdef WINDOWS
  delete[] DiagPtr;
  delete[] FullVector;
#else
  free(DiagPtr);
  free(FullVector);
#endif
/*  fprintf(stderr,"  Repeat Solution Multiplications (Tau+N): %ld\n",Nmult1);*/
}
