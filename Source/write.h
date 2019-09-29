#include <stdlib.h>
//#ifndef WINDOWS
//#include <stdio.h>
//#else
#include "pfwstdio.h"
//#endif
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"


#ifdef ANSIPROTO
void Output(INDEX Iter,char *File1,char *str);
void Print(FILE *File,int spaces,int width,int decimals,VALUETYPE val);
void IEEE(void);
void DeleteJac(SparseMatrix *Mptr,IntegerVector *P1Row,IntegerVector *P1Col,
	       IntegerVector *P2Row,IntegerVector *P2Col);
AreaData *ACFunJac(SparseMatrix *Mptr,int *val,BOOLEAN flagF,BOOLEAN flagJ,BOOLEAN flagP);
BOOLEAN DCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);
int HFunJac(BOOLEAN flagF,BOOLEAN flagJ,AreaData *Aptr,VALUETYPE *vec);
void ACFunHes(BOOLEAN flagF,BOOLEAN flagJ);
BOOLEAN DCFunHes(BOOLEAN flagF,BOOLEAN flagJ);
void SVCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);          /* FACTS */
void TCSCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);          /* FACTS */
void STATCOMFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);       /* FACTS */
void SVCFunHes(BOOLEAN flagF,BOOLEAN flagJ);        /* FACTS */
void TCSCFunHes(BOOLEAN flagF,BOOLEAN flagJ);       /* FACTS */
void STATCOMFunHes(BOOLEAN flagF,BOOLEAN flagJ);    /* FACTS */
void WriteJac(void);
void WriteSolution(INDEX Iter,char *File,char *str);
void WriteQmatrix(INDEX count,VALUETYPE *vec);
#else
void Output();
void Print();
void IEEE();
void DeleteJac();
AreaData *ACFunJac();
BOOLEAN DCFunJac();
int HFunJac();
void ACFunHes();
BOOLEAN DCFunHes();
void SVCFunJac();      // FACTS 
void TCSCFunJac();     // FACTS 
void STATCOMFunJac();  // FACTS 
void SVCFunHes();      // FACTS 
void TCSCFunHes();     //. FACTS 
void STATCOMFunHes();  // FACTS 
void WriteJac();
void WriteSolution();
void WriteQmatrix();
#endif

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX MaxIter,Nac,NacVar,NacEl,NregPQ,NregV,Ndc,Narea,Nslack,Nvolt,Bl,NZvolt,NXvolt,
             Nsvc,Ntcsc,NtcscVar,Nstatcom; /* FACTS */
extern VALUETYPE Tol,Sn,lambda,lambda_o,*dF,MaxdFi,*x0,*Dx,K3;
extern int GlobalArgc;
extern char **GlobalArgv;
extern BOOLEAN Acont,PQcont,QRcont,Rcont,flagH,flagPoC,flagBS,flagPgMax,flagSmax;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern IntegerVector *RowPer,*ColPer;
extern FILE *OutputHomot;


