#pragma once

#include "constant.h"
#include "param.h"
#include "pflow.h"
#include "sparse.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void JacElement(SparseMatrix *Mptr, INDEX I, INDEX J, VALUETYPE val);
AreaData *ACFunJac(SparseMatrix *Mptr, int *val, bool flagF, bool flagJ, bool flagP);
bool DCFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ);
void Jacobian(void);

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX Nac, Ndc, NacVar, Narea, Bl, Ngen, Nsvc, Ntcsc, NtcscVar, Nstatcom; /* FACTS */
extern INDEX *ACvar;
extern VALUETYPE *dx, *dF, Sn, lambda, *Dx, *x0, *x0p;
extern VALUETYPE K1, K2;
extern IntegerVector *NewRow, *OldRow, *NewCol, *OldCol, *RowPartition, *ColPartition;
extern bool Acont, PQcont, QRcont, Rcont, flagH, flagPoC, flagBS, flag2Vcontrol, flagPgMax, flagSmax, flagReducedContinuation, *DxZero;
extern bool *DxZero;
extern AClist *Vlist, *Vlistp;
