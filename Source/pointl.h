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
void WriteSolution(INDEX Iter,char *File1,char *str);
void DeleteJac(SparseMatrix *Mptr,IntegerVector *P1Row,IntegerVector *P1Col,
	       IntegerVector *P2Row,IntegerVector *P2Col);
int factorns(SparseMatrix *Mptr,double Param,IntegerVector *PartRow,IntegerVector *PartCol,
	     IntegerVector *P1Row,IntegerVector *P1Col,IntegerVector *P2Row,IntegerVector *P2Col);
void repsolp(SparseMatrix *Mptr,VALUETYPE *Vptr,
	     IntegerVector *PermR,IntegerVector *PermC);
int factor(SparseMatrix *Mptr);
bool DCsetup(void);
int Pflow(int iter,bool flagF,bool flagD,bool flagP);
AreaData *ACFunJac(SparseMatrix *Mptr,int *val,bool flagF,bool flagJ,bool flagP);
bool DCFunJac(SparseMatrix *Mptr,bool flagF,bool flagJ);
void SVCFunJac(SparseMatrix *Mptr,bool flagF,bool flagJ);         /* FACTS */
void TCSCFunJac(SparseMatrix *Mptr,bool flagF,bool flagJ);        /* FACTS */
void STATCOMFunJac(SparseMatrix *Mptr,bool flagF,bool flagJ);     /* FACTS */
void ACFunHes(bool flagF,bool flagJ);
bool DCFunHes(bool flagF,bool flagJ);
void SVCFunHes(bool flagF,bool flagJ);      /* FACTS */
void TCSCFunHes(bool flagF,bool flagJ);     /* FACTS */
void STATCOMFunHes(bool flagF,bool flagJ);  /* FACTS */
void UpdateEvector(VALUETYPE cons);
bool CheckRlimits(void);
bool CheckVlimits(void);
bool CheckQlimits(void);
bool ChangeDCmode(void);
bool ChangeSVCmode(void);       /* FACTS */
bool ChangeTCSCmode(void);      /* FACTS */
bool ChangeSTATCOMmode(void);   /* FACTS */
int Direction(SparseMatrix *Mptr,VALUETYPE *vec,bool flag);
VALUETYPE LoadX0(bool flagX,bool flagV,bool flag);
int Evector(int M,int iter,VALUETYPE tol,bool RightEvector,VALUETYPE *EigenValue);
int PoCPoint(void);

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX Nac,NacEl,NregPQ,NregV,Ndc,Nslack,Nvolt,NZvolt,NXvolt,Narea,NacVar,MaxIter,
             Nsvc,Ntcsc,NtcscVar,Nstatcom; /* FACTS */
extern INDEX *ACvar;
extern VALUETYPE *dx,*dF,Sn,lambda,param0,Dparam,*x0,*Dx;
extern VALUETYPE K1,K2,Tol,alpha;
extern bool Acont,PQcont,QRcont,Rcont,PQlim,Tlim,Qlim,Vlim,Zlim,Xlim,
               flagH,flagPoC,flagL,flagR,flagPgMax,flagSmax;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern IntegerVector *RowPer,*ColPer;
