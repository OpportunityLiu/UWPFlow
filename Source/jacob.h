#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"

void JacElement(SparseMatrix *Mptr,INDEX I,INDEX J,VALUETYPE val);
AreaData *ACFunJac(SparseMatrix *Mptr,int *val,BOOLEAN flagF,BOOLEAN flagJ,BOOLEAN flagP);
BOOLEAN DCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);
void Jacobian(void);

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX Nac,Ndc,NacVar,Narea,Bl,Ngen,
             Nsvc,Ntcsc,NtcscVar,Nstatcom; /* FACTS */
extern INDEX *ACvar;
extern VALUETYPE *dx,*dF,Sn,lambda,*Dx,*x0,*x0p;
extern VALUETYPE K1,K2;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern BOOLEAN Acont,PQcont,QRcont,Rcont,flagH,flagPoC,flagBS,flag2Vcontrol,
	              flagPgMax,flagSmax,flagReducedContinuation,*DxZero;
extern BOOLEAN *DxZero;
extern AClist *Vlist,*Vlistp;

