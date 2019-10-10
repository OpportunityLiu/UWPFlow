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
void DeleteJac(SparseMatrix *Mptr,IntegerVector *P1Row,IntegerVector *P1Col,
	       IntegerVector *P2Row,IntegerVector *P2Col);
AreaData *ACFunJac(SparseMatrix *Mptr,int *val,BOOLEAN flagF,BOOLEAN flagJ,BOOLEAN flagP);
BOOLEAN DCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);
void SVCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);         /* FACTS */
void TCSCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);        /* FACTS */
void STATCOMFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);     /* FACTS */
int HFunJac(BOOLEAN FlagFunction,BOOLEAN FlagJacobian,AreaData *Aptr,VALUETYPE *vec);
VALUETYPE LoadX0(BOOLEAN FlagLoadX0,BOOLEAN FlagUpdateVar,BOOLEAN FlagMakeDxZero);
int factorns(SparseMatrix *Mptr,double Param,IntegerVector *PartRow,IntegerVector *PartCol,
	     IntegerVector *P1Row,IntegerVector *P1Col,IntegerVector *P2Row,IntegerVector *P2Col);
int factor(SparseMatrix *Mptr);
void repsolp(SparseMatrix *Mptr,VALUETYPE *Vptr,
	     IntegerVector *PermR,IntegerVector *PermC);
void WriteSolution(INDEX Iter,char *File1,char *str);
int Pflow(int iter,BOOLEAN flagF,BOOLEAN flagD,BOOLEAN flagP);
void MakeVlist(FILE *Out);
void VoltProf(BOOLEAN flag,FILE *Out);
int Direction(SparseMatrix *Mptr,VALUETYPE *vec,BOOLEAN flag);
BOOLEAN ChangeParam(void);
int Homot(void);
BOOLEAN CheckRlimits(void);
BOOLEAN CheckVlimits(void);
BOOLEAN CheckQlimits(void);
BOOLEAN ChangeDCmode(void);
BOOLEAN ChangeSVCmode(void);      /* FACTS */
BOOLEAN ChangeTCSCmode(void);     /* FACTS */
BOOLEAN ChangeSTATCOMmode(void);  /* FACTS */
BOOLEAN DCsetup(void);
void Print(FILE *File,int spaces,int width,int decimals,VALUETYPE val);
void VoltageSensFactor(VALUETYPE *vec,BOOLEAN first);
void TEFac(BOOLEAN flag);
void TEFdc(FILE *Out);
void MatlabV(FILE *Out);
void TEFMatlabFiles(void);
BOOLEAN InList(ACbusData *ACptr,AClist *Vptr);
void PrintDirection(char Option,VALUETYPE *vector,VALUETYPE Max);
void WriteQmatrix(INDEX count,VALUETYPE *vec);
void IndicesMatlab(INDEX count);

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX Nac,NacEl,Ndc,Narea,NacVar,Nvolt,NregV,Nslack,
             Nsvc,Ntcsc,NtcscVar,Nstatcom; /* FACTS */
extern INDEX *ACvar,Bl;
extern VALUETYPE *dF,lambda,alpha,Vac,Sn,*dx,Tol;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern BOOLEAN Acont,PQcont,QRcont,Rcont,PQlim,Tlim,Qlim,Vlim,Elim,Ilim,Xlim,Zlim,
	       flagH,flagL,flagR,flagReducedContinuation,flagPgMax,flagSmax;
extern ACbusData *BlPtr;
extern INDEX TFnum,TFbus;
extern int DetSign;
extern char TFname[13];


typedef struct ACranklist {
    struct ACbusData *AC;
    VALUETYPE val;
    struct ACranklist *Next;
  } ACranklist;

