/* Factor unsymmetric square or rectangular matrix, add fills.
   Purpose: Factor Unsymmetric Matrix, add fills and
            Generate Permutation Vectors.
            Order according to min Degree rows / max Value OR min degree cols.
            Order within row partitions if partition vector is given. */

#include "constant.h"
#include "param.h"
#include "sparse.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DeltaValence 3
/* This is the only "trick": If row seems singular,
           defer factorization by increasing its valence.
           For square matrices, this happens only if singular.
           For rectangular factorization, this can happen. */

void InitializeLSWR(void);
void AllocateOrderArrays(void);
void LoadUnloadWorkingRow(INDEX i, bool SetReset);
SparseMatrixElement *PutMatrix(INDEX PutRow, INDEX PutCol, ELEMENTVALUETYPE PutVal, SparseMatrixElement *PutRNext, SparseMatrixElement *PutCNext);
void SimElim(INDEX Ipivot, INDEX Jpivot, ELEMENTVALUETYPE *DiagVal);
void LinkVal(INDEX Ival, VALENCE valVal);
void UNLinkVal(INDEX Iunl);
void FindVal(INDEX Ifnd);
void FormClass(INDEX Ibeg, INDEX Iend);
INDEX GetLowval(void);
VALENCE ColumnValence(INDEX Col);
INDEX GetHighest(INDEX Ipivot);
void UpdateVals(INDEX Jpivot);
int OrderMatrix(INDEX Ibeg, INDEX Iend);
void InvertNormalize(void);
void FinishVectors(void);
int factorns(SparseMatrix *Mptr, double Param, IntegerVector *PartRow, IntegerVector *PartCol, IntegerVector *P1Row, IntegerVector *P1Col, IntegerVector *P2Row,
             IntegerVector *P2Col);

/* Global array and variable definitions */
SparseMatrix *Matrix;
LONGINT Tau, Nfills, Nmult;
VALENCE maxVal = maxN;

bool *Lswr;            /* array 0..maxN */
ELEMENTVALUETYPE *Swr; /* array 0..maxN */
IntegerVector *RowNew, *RowOld, *ColNew, *ColOld;
INDEX *PartitionAt, *PartitionCol; /* array 0..maxPart */

double Alpha;
int LinkCount, Iorder, Nclasses;
INDEX Nstop, Ibeg, Iend, Ipart, Nparts;
INDEX CurrentColPart, NpartCol;

/* Following definitions global only within OrderMatrix routines */
VALENCE *Valence; /* array 0..maxN */
INDEX *FwdLink;   /* array 0..maxN */
INDEX *BckLink;   /* array 0..maxN */
INDEX *ValClass;  /* array 0..maxVal */
ELEMENTVALUETYPE PivotV, DiagVal;

/* ================= InitializeLSWR =================================== */
void InitializeLSWR() {
  int i;

  Lswr = new bool[Matrix->n1 + 1];
  Swr = new ELEMENTVALUETYPE[Matrix->n1 + 1];
  for (i = 1; i <= Matrix->n1; i++)
    Lswr[i] = false;
}

/* =================== AllocateOrderArrays ========================== */
void AllocateOrderArrays() {
  INDEX i, n0;
  /* BEGIN */
  if (Matrix->n1 > Matrix->n2)
    n0 = Matrix->n1;
  else
    n0 = Matrix->n2;
  maxVal = n0;                      /* This could be larger yet */
  Valence = new VALENCE[n0 + 1];    /* array 0..maxN */
  FwdLink = new INDEX[n0 + 1];      /* array 0..maxN */
  BckLink = new INDEX[n0 + 1];      /* array 0..maxN */
  ValClass = new INDEX[maxVal + 1]; /* array 0..maxVal */
  for (i = 0; i < n0 + 1; i++) {
    Valence[i] = 0;
    FwdLink[i] = BckLink[i] = 0;
  }
  for (i = 0; i < maxVal + 1; i++)
    ValClass[i] = 0;
}

/* ====================== LoadUnloadWorkingRow ======================== */
void LoadUnloadWorkingRow(INDEX i, bool SetReset) {
  SparseMatrixElement *Ptr2;
  /* BEGIN */
  Ptr2 = Matrix->RowHead[i];
  while (Ptr2 != nullptr) {
    if (SetReset)
      Tau++;
    Lswr[Ptr2->Col] = SetReset;
    Ptr2 = Ptr2->RowNext;
  }
}

/* ========================== PutMatrix =============================== */
SparseMatrixElement *PutMatrix(INDEX PutRow, INDEX PutCol, ELEMENTVALUETYPE PutVal, SparseMatrixElement *PutRNext,
                               SparseMatrixElement *PutCNext) { /* PutMatrix */
  SparseMatrixElement *PutPtr;
  /* BEGIN */
  PutPtr = new SparseMatrixElement;
  if (PutPtr != nullptr) {
    PutPtr->Row = PutRow;
    PutPtr->Col = PutCol;
    PutPtr->Value = PutVal;
    PutPtr->RowNext = PutRNext;
    PutPtr->ColNext = PutCNext;
  } else {
    ErrorHalt("Insufficient memory for fills");
    exit(ERROREXIT);
  }
  return (PutPtr);
} /* PutMatrix */

/* =========================== SimElim ============================== */
void SimElim(INDEX Ipivot, INDEX Jpivot, ELEMENTVALUETYPE *DiagVal) {
  INDEX Jsim, Ksim;
  SparseMatrixElement *OrdP1;
  SparseMatrixElement *OrdP2;
  /* BEGIN SimElim */
  OrdP1 = Matrix->ColHead[Jpivot];
  while (OrdP1 != nullptr) {
    Jsim = OrdP1->Row;
    if (RowNew->p[Jsim] == 0) {
      OrdP2 = Matrix->RowHead[Jsim];
      while (OrdP2 != nullptr) {
        Lswr[OrdP2->Col] = true;
        Swr[OrdP2->Col] = OrdP2->Value;
        OrdP2 = OrdP2->RowNext;
      }

      /* DivideValues(&Swr[Jpivot],Swr[Jpivot],DiagVal); */
      Swr[Jpivot] = Swr[Jpivot] / *DiagVal;

      OrdP2 = Matrix->RowHead[Ipivot];
      while (OrdP2 != nullptr) {
        Ksim = OrdP2->Col;
        if (ColNew->p[Ksim] == 0) {
          if (Lswr[Ksim] != true) {
            Nfills++;
            Matrix->RowHead[Jsim] = PutMatrix(Jsim, Ksim, 0.0, Matrix->RowHead[Jsim], Matrix->ColHead[Ksim]);
            Matrix->ColHead[Ksim] = Matrix->RowHead[Jsim];
            Lswr[Ksim] = true;
            Swr[Ksim] = 0.0;
          } /* END IF */
          /*
                  MultiplyValues(TempVal,Swr[Jpivot],OrdP2^.Value);
                  SubtractValues(Swr[Ksim],Swr[Ksim],TempVal);
                  */
          Swr[Ksim] = Swr[Ksim] - Swr[Jpivot] * OrdP2->Value;
          Nmult++;
        } /* END IF */
        OrdP2 = OrdP2->RowNext;
      } /* END WHILE */;
      OrdP2 = Matrix->RowHead[Jsim];
      while (OrdP2 != nullptr) {
        Lswr[OrdP2->Col] = false;
        OrdP2->Value = Swr[OrdP2->Col];
        Swr[OrdP2->Col] = 0.0;
        OrdP2 = OrdP2->RowNext;
      } /* END WHILE */
    }   /* END IF */
    OrdP1 = OrdP1->ColNext;
  } /* END WHILE */
} /* END SimElim */

/* ============== START OF ROUTINES LOCAL TO OrderMatrix ============= */

void LinkVal(INDEX Ival, VALENCE valVal) {
  if (ValClass[valVal] != 0)
    BckLink[ValClass[valVal]] = Ival;
  FwdLink[Ival] = ValClass[valVal];
  BckLink[Ival] = 0;
  ValClass[valVal] = Ival;
  LinkCount++;
} /* END LinkVal */

void UNLinkVal(INDEX Iunl) {
  if (BckLink[Iunl] != 0)
    FwdLink[BckLink[Iunl]] = FwdLink[Iunl];
  else
    ValClass[Valence[Iunl]] = FwdLink[Iunl];
  if (FwdLink[Iunl] != 0)
    BckLink[FwdLink[Iunl]] = BckLink[Iunl];
  LinkCount--;
} /* END UNLinkVal */

void FindVal(INDEX Ifnd) {
  SparseMatrixElement *OrdP1;
  /* BEGIN */
  Valence[Ifnd] = 0;
  OrdP1 = Matrix->RowHead[Ifnd];
  while (OrdP1 != nullptr) {
    if ((Valence[Ifnd] < maxVal) && (ColNew->p[OrdP1->Col] == 0)) {
      Valence[Ifnd]++;
    } /* END IF */
    OrdP1 = OrdP1->RowNext;
  } /* END WHILE */
} /* END FindVal */

void FormClass(INDEX Ibeg, INDEX Iend) {
  INDEX I;
  /* BEGIN */
  for (I = Ibeg; I <= Iend; I++) {
    RowOld->p[I] = 0;
    ColOld->p[I] = 0;
  }
  for (I = Ibeg; I <= Iend; I++)
    FindVal(I);
  for (I = 0; I <= maxVal; I++)
    ValClass[I] = 0;
  for (I = Ibeg; I <= Iend; I++)
    LinkVal(I, Valence[I]);
} /* END FormClass */

INDEX GetLowval() {
  INDEX Iget;
  VALENCE ValGet;
  /* BEGIN */
  ValGet = 0;
  do {
    Iget = ValClass[ValGet];
    ValGet++;
  } while (Iget == 0);
  UNLinkVal(Iget);
  return (Iget);
}

VALENCE ColumnValence(INDEX Col) {
  VALENCE TempCount;
  SparseMatrixElement *ColValPtr;
  /* BEGIN */
  TempCount = 0;
  ColValPtr = Matrix->ColHead[Col];
  while (ColValPtr != nullptr) {
    TempCount++;
    ColValPtr = ColValPtr->ColNext;
  }
  return (TempCount);
}

INDEX GetHighest(INDEX Ipivot) {
  INDEX LargestI, SmallestJ;
  VALENCE ColVal;
  SparseMatrixElement *GetPtr;
  ELEMENTVALUETYPE LargestV;
  float AlphaV;
  /* BEGIN */
  LargestV = 0;
  LargestI = 0;
  /* Get the largest pivot on row: */
  GetPtr = Matrix->RowHead[Ipivot];
  while ((GetPtr != nullptr) && (GetPtr->Col <= PartitionCol[CurrentColPart])) {
    if (ColNew->p[GetPtr->Col] == 0) {
      if (LargestV <= fabs(GetPtr->Value)) {
        LargestV = fabs(GetPtr->Value);
        PivotV = LargestV;
        DiagVal = GetPtr->Value;
        LargestI = GetPtr->Col;
      }
    }
    GetPtr = GetPtr->RowNext;
  }
  if (Alpha == 1.0)
    return (LargestI); /* The most robust case */
  SmallestJ = maxN + 1;
  AlphaV = LargestV * Alpha;
  GetPtr = Matrix->RowHead[Ipivot];
  while ((GetPtr != nullptr) && (GetPtr->Col <= PartitionCol[CurrentColPart])) {
    if (ColNew->p[GetPtr->Col] == 0) {
      ColVal = ColumnValence(GetPtr->Col);
      if ((ColVal <= SmallestJ) && (fabs(GetPtr->Value) >= AlphaV)) {
        SmallestJ = ColVal;
        PivotV = fabs(GetPtr->Value);
        DiagVal = GetPtr->Value;
        LargestI = GetPtr->Col;
      }
    }
    GetPtr = GetPtr->RowNext;
  }
  return (LargestI);
} /* END GetHighest */

void UpdateVals(INDEX Jpivot) {
  SparseMatrixElement *UpdP1;
  INDEX Jord;
  /* BEGIN */
  UpdP1 = Matrix->ColHead[Jpivot];
  while (UpdP1 != nullptr) {
    Jord = UpdP1->Row;
    if ((RowNew->p[Jord] == 0) && (Jord >= Ibeg) && (Jord <= Iend)) {
      UNLinkVal(Jord);
      FindVal(Jord);
      LinkVal(Jord, Valence[Jord]);
    } /* END IF */
    UpdP1 = UpdP1->ColNext;
  } /* END WHILE */
} /* END UpdateVals */

int OrderMatrix(INDEX Ibeg, INDEX Iend) {
  int LoopLimit;
  INDEX Istop, Ipivot, Jpivot;
  /* BEGIN OrderMatrix */
  FormClass(Ibeg, Iend);
  Istop = Iend;
  if (Istop > Matrix->n1)
    Istop = Matrix->n1;
  if (Istop > Matrix->n2)
    Istop = Matrix->n2;
  if (Istop > Nstop)
    Istop = Nstop;
  while (Iorder < Istop) {
    LoopLimit = Matrix->n1;
    Iorder++;

    if (Iorder > PartitionCol[CurrentColPart])
      CurrentColPart++;

    do {
      LoopLimit--;
      Ipivot = GetLowval();
      Jpivot = GetHighest(Ipivot);
      if (Jpivot == 0) {
        UNLinkVal(Ipivot);
        if (Matrix->n1 == Matrix->n2) {
          fprintf(stderr, "\nUnable to find a nonzero pivot for row %d\n", RowOld->p[Ipivot]);
          fprintf(stderr, "Matrix is probably numerically singular\n");
          return (WARNINGEXIT);
        }
        Valence[Ipivot] = Valence[Ipivot] + DeltaValence;
        LinkVal(Ipivot, Valence[Ipivot]);
      }
    } while (((Jpivot <= 0) || (PivotV <= SINGULARITYZERO)) && (LinkCount > 0) && (LoopLimit >= 0));

    if (Jpivot <= 0) {
      fprintf(stderr, "\nUnable to find an available pivot for row %d\n", RowOld->p[Ipivot]);
      fprintf(stderr, "Matrix is probably topologically singular\n");
      return (WARNINGEXIT);
      /* Iorder = Istop; */
    } else {
      RowNew->p[Ipivot] = Iorder;
      RowOld->p[Iorder] = Ipivot;
      if (NearZero(PivotV)) {
        fprintf(stderr, "\nEquation %d Depends on Other Equations\n", RowOld->p[Ipivot]);
        fprintf(stderr, "Matrix is probably numerically singular\n");
        return (WARNINGEXIT);
        /* Iorder = Istop; */
      } else {
        ColNew->p[Jpivot] = Iorder;
        ColOld->p[Iorder] = Jpivot;
        SimElim(Ipivot, Jpivot, &DiagVal);
        UpdateVals(Jpivot);
      }
    }
  } /* END WHILE */
  return (0);
} /* END OrderMatrix */

/* ============================ InvertNormalize ======================== */
void InvertNormalize() {
  INDEX I, Inew;
  SparseMatrixElement *Ptr1;
  ELEMENTVALUETYPE DiagValue;
  /* BEGIN {InvertNormalize} */
  for (Inew = 1; Inew <= Nstop; Inew++) {
    I = RowOld->p[Inew];
    Ptr1 = Matrix->RowHead[I];
    while (Ptr1 != nullptr) {
      if (RowNew->p[I] == ColNew->p[Ptr1->Col]) {
        /* InvertValue(Ptr1->Value); */
        /* EquateValues(DiagValue,Ptr1->Value); */
        Ptr1->Value = 1.0 / Ptr1->Value;
        DiagValue = Ptr1->Value;
      }
      Ptr1 = Ptr1->RowNext;
    }
    Ptr1 = Matrix->RowHead[I];
    while (Ptr1 != nullptr) {
      if (ColNew->p[Ptr1->Col] > RowNew->p[I]) {
        /* MultiplyValues(Ptr1^.Value,Ptr1^.Value,DiagValue); */
        Ptr1->Value = Ptr1->Value * DiagValue;
        Nmult++;
      }
      Ptr1 = Ptr1->RowNext;
    }
  }
} /* END {InvertNormalize} */

/* ========================== FinishVectors =========================== */
void FinishVectors() {
  INDEX I, K;
  /* BEGIN */
  K = 1;
  while ((RowOld->p[K] != 0) && (K <= Matrix->n1))
    K++;
  for (I = 1; I <= Matrix->n1; I++) {
    if (RowNew->p[I] == 0) {
      RowNew->p[I] = K;
      K++;
    }
  }
  for (I = 1; I <= Matrix->n1; I++) {
    K = RowNew->p[I];
    RowOld->p[K] = I;
  }
  K = 1;
  while ((ColOld->p[K] != 0) && (K <= Matrix->n2))
    K++;
  for (I = 1; I <= Matrix->n2; I++) {
    if (ColNew->p[I] == 0) {
      ColNew->p[I] = K;
      K++;
    }
  }
  for (I = 1; I <= Matrix->n2; I++)
    ColOld->p[ColNew->p[I]] = I;
}

/* ============================= factorns ================================ */
int factorns(SparseMatrix *Mptr, double Param, IntegerVector *PartRow, IntegerVector *PartCol, IntegerVector *P1Row, IntegerVector *P1Col, IntegerVector *P2Row,
             IntegerVector *P2Col) { /* FactorNS */
  double MinAlpha = 0.0;
  double MaxAlpha = 1.0;
  SparseMatrixElement *Ptr;
  INDEX i;
  FILE *OutFile;

  /* BEGIN */
  RowNew = P1Row;
  ColNew = P1Col;
  RowOld = P2Row;
  ColOld = P2Col;
  Alpha = Param;
  if (Alpha < MinAlpha)
    Alpha = MinAlpha;
  if (Alpha > MaxAlpha)
    Alpha = MaxAlpha;
  Matrix = Mptr;
  InitializeLSWR();
  AllocateOrderArrays();
  if (Matrix->n1 < Matrix->n2)
    Nstop = Matrix->n1;
  else
    Nstop = Matrix->n2;
  PartitionAt = PartRow->p;
  Nparts = PartRow->N;
  PartitionCol = PartCol->p;
  NpartCol = PartCol->N;
  CurrentColPart = 1;
  Iorder = 0;
  LinkCount = 0;
  Nfills = 0;
  Nmult = 0;
  for (Ipart = 1; Ipart <= Nparts; Ipart++) {
    Ibeg = PartitionAt[Ipart - 1] + 1;
    Iend = PartitionAt[Ipart];
    if (OrderMatrix(Ibeg, Iend))
      return (1);
  }
  /*  fprintf(stderr,"  Factorization Fills: %d",Nfills);*/
  FinishVectors();   /* For the case of rectangular matrix factorization */
  InvertNormalize(); /* To convert to true LDU factors, as expected by REPSOL */
                     /*  fprintf(stderr,"    Multiplications: %d\n",Nmult);*/
  for (i = 1; i <= Matrix->n1; i++) {
    Ptr = Matrix->RowHead[i];
    while (Ptr != nullptr) {
      Ptr->Row = RowNew->p[Ptr->Row];
      Ptr->Col = ColNew->p[Ptr->Col];
      Ptr = Ptr->RowNext;
    }
  }
  delete[] Lswr;
  delete[] Swr;
  delete[] Valence;
  delete[] FwdLink;
  delete[] BckLink;
  delete[] ValClass;
  /*
  if (ExistParameter('d')) {
    OutFile=(FILE *) OpenOutput("prow.dat");
    fprintf(OutFile,"%d\n",RowOld->N);
    for (i=1; i<=RowOld->N; i++) fprintf(OutFile,"%d ",RowOld->p[i]);
    fprintf(OutFile,"\n");
    fclose(OutFile);
    OutFile=(FILE *) OpenOutput("pcol.dat");
    fprintf(OutFile,"%d\n",ColOld->N);
    for (i=1; i<=ColOld->N; i++) fprintf(OutFile,"%d ",ColOld->p[i]);
    fprintf(OutFile,"\n");
    fclose(OutFile);
  }
*/
  return (0);
} /* FactorNS */
