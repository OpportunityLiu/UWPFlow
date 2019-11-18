#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"

void Output(INDEX Iter,char *File1,char *str);
void Print(FILE *File,int spaces,int width,int decimals,VALUETYPE val);
void IEEE(void);
void DeleteJac(SparseMatrix *Mptr,IntegerVector *P1Row,IntegerVector *P1Col,
	       IntegerVector *P2Row,IntegerVector *P2Col);
AreaData *ACFunJac(SparseMatrix *Mptr,int *val,bool flagF,bool flagJ,bool flagP);
bool DCFunJac(SparseMatrix *Mptr,bool flagF,bool flagJ);
int HFunJac(bool flagF,bool flagJ,AreaData *Aptr,VALUETYPE *vec);
void ACFunHes(bool flagF,bool flagJ);
bool DCFunHes(bool flagF,bool flagJ);
void SVCFunJac(SparseMatrix *Mptr,bool flagF,bool flagJ);          /* FACTS */
void TCSCFunJac(SparseMatrix *Mptr,bool flagF,bool flagJ);          /* FACTS */
void STATCOMFunJac(SparseMatrix *Mptr,bool flagF,bool flagJ);       /* FACTS */
void SVCFunHes(bool flagF,bool flagJ);        /* FACTS */
void TCSCFunHes(bool flagF,bool flagJ);       /* FACTS */
void STATCOMFunHes(bool flagF,bool flagJ);    /* FACTS */
void WriteJac(void);
void WriteSolution(INDEX Iter,char *File,char *str);
void WriteQmatrix(INDEX count,VALUETYPE *vec);

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX MaxIter,Nac,NacVar,NacEl,NregPQ,NregV,Ndc,Narea,Nslack,Nvolt,Bl,NZvolt,NXvolt,
             Nsvc,Ntcsc,NtcscVar,Nstatcom; /* FACTS */
extern VALUETYPE Tol,Sn,lambda,lambda_o,*dF,MaxdFi,*x0,*Dx,K3;
extern int GlobalArgc;
extern char **GlobalArgv;
extern bool Acont,PQcont,QRcont,Rcont,flagH,flagPoC,flagBS,flagPgMax,flagSmax;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern IntegerVector *RowPer,*ColPer;
extern FILE *OutputHomot;


