#include <stdlib.h>
//#ifndef WINDOWS
//#include <stdio.h>
//#else
#include "pfwstdio.h"
//#endif
#include <math.h>
#include <string.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"

#ifdef ANSIPROTO
void JacElement(SparseMatrix *Mptr,INDEX I,INDEX J,VALUETYPE val);
void WriteSolution(INDEX Iter,char *File1,char *str);
void DeleteJac(SparseMatrix *Mptr,IntegerVector *P1Row,IntegerVector *P1Col,
	       IntegerVector *P2Row,IntegerVector *P2Col);
int factorns(SparseMatrix *Mptr,double Param,IntegerVector *PartRow,IntegerVector *PartCol,
	     IntegerVector *P1Row,IntegerVector *P1Col,IntegerVector *P2Row,IntegerVector *P2Col);
void repsolp(SparseMatrix *Mptr,VALUETYPE *Vptr,
	     IntegerVector *PermR,IntegerVector *PermC);
int factor(SparseMatrix *Mptr);
BOOLEAN DCsetup(void);
int Pflow(int iter,BOOLEAN flagF,BOOLEAN flagD,BOOLEAN flagP);
AreaData *ACFunJac(SparseMatrix *Mptr,int *val,BOOLEAN flagF,BOOLEAN flagJ,BOOLEAN flagP);
BOOLEAN DCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);
void SVCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);         /* FACTS */
void TCSCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);        /* FACTS */
void STATCOMFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);     /* FACTS */
void ACFunHes(BOOLEAN flagF,BOOLEAN flagJ);
BOOLEAN DCFunHes(BOOLEAN flagF,BOOLEAN flagJ);
void SVCFunHes(BOOLEAN flagF,BOOLEAN flagJ);      /* FACTS */
void TCSCFunHes(BOOLEAN flagF,BOOLEAN flagJ);     /* FACTS */
void STATCOMFunHes(BOOLEAN flagF,BOOLEAN flagJ);  /* FACTS */
void UpdateEvector(VALUETYPE cons);
BOOLEAN CheckRlimits(void);
BOOLEAN CheckVlimits(void);
BOOLEAN CheckQlimits(void);
BOOLEAN ChangeDCmode(void);
BOOLEAN ChangeSVCmode(void);       /* FACTS */
BOOLEAN ChangeTCSCmode(void);      /* FACTS */
BOOLEAN ChangeSTATCOMmode(void);   /* FACTS */
int Direction(SparseMatrix *Mptr,VALUETYPE *vec,BOOLEAN flag);
VALUETYPE LoadX0(BOOLEAN flagX,BOOLEAN flagV,BOOLEAN flag);
int Evector(int M,int iter,VALUETYPE tol,BOOLEAN RightEvector,VALUETYPE *EigenValue);
int PoCPoint(void);
#else
void JacElement();
void WriteSolution();
void DeleteJac();
int factorns();
void repsolp();
int factor();
BOOLEAN DCsetup();
int Pflow();
AreaData *ACFunJac();
BOOLEAN DCFunJac();
void SVCFunJac();       //FACTS 
void TCSCFunJac();      //FACTS 
void STATCOMFunJac();   //FACTS 
void ACFunHes();
BOOLEAN DCFunHes();
void SVCFunHes();      //FACTS 
void TCSCFunHes();     //FACTS 
void STATCOMFunHes();  //FACTS 
void UpdateEvector();
BOOLEAN CheckRlimits();
BOOLEAN CheckVlimits();
BOOLEAN CheckQlimits();
BOOLEAN ChangeDCmode();
BOOLEAN ChangeSVCmode();     //FACTS 
BOOLEAN ChangeTCSCmode();    //FACTS 
BOOLEAN ChangeSTATCOMmode(); //FACTS 
int Direction();
VALUETYPE LoadX0();
int Evector();
int PoCPoint();
#endif

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX Nac,NacEl,NregPQ,NregV,Ndc,Nslack,Nvolt,NZvolt,NXvolt,Narea,NacVar,MaxIter,
             Nsvc,Ntcsc,NtcscVar,Nstatcom; /* FACTS */
extern INDEX *ACvar;
extern VALUETYPE *dx,*dF,Sn,lambda,param0,Dparam,*x0,*Dx;
extern VALUETYPE K1,K2,Tol,alpha;
extern BOOLEAN Acont,PQcont,QRcont,Rcont,PQlim,Tlim,Qlim,Vlim,Zlim,Xlim,
               flagH,flagPoC,flagL,flagR,flagPgMax,flagSmax;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern IntegerVector *RowPer,*ColPer;
