#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

/* Power Flow. */

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
void ErrorStop(char *Msg);
AreaData *ACFunJac(SparseMatrix *Mptr,int *val,BOOLEAN flagF,BOOLEAN flagJ,BOOLEAN flagFirst);
BOOLEAN DCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);
void UpdateSVCvar(VALUETYPE cons,INDEX j);                               /* FACTS */
void SVCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);          /* FACTS */
void UpdateTCSCvar(VALUETYPE cons,INDEX j);                              /* FACTS */
void TCSCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);         /* FACTS */
void UpdateSTATCOMvar(VALUETYPE cons,INDEX j);                           /* FACTS */
void STATCOMFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ);      /* FACTS */
int HFunJac(BOOLEAN flagF,BOOLEAN flagJ,AreaData *Aptr,VALUETYPE *vec);
void ACFunHes(BOOLEAN flagF,BOOLEAN flagJ);
BOOLEAN DCFunHes(BOOLEAN flagF,BOOLEAN flagJ);
void SVCFunHes(BOOLEAN flagF,BOOLEAN flagJ);       /* FACTS */
void TCSCFunHes(BOOLEAN flagF,BOOLEAN flagJ);      /* FACTS */
void STATCOMFunHes(BOOLEAN flagF,BOOLEAN flagJ);   /* FACTS */
int factorns(SparseMatrix *Mptr,double Param,IntegerVector *PartRow,IntegerVector *PartCol,
             IntegerVector *P1Row,IntegerVector *P1Col,IntegerVector *P2Row,IntegerVector *P2Col);
void repsolp(SparseMatrix *Mptr,VALUETYPE *Vptr,
             IntegerVector *PermR,IntegerVector *PermC);
VALUETYPE Norm(VALUETYPE *Vptr,INDEX N,INDEX *N1);
void WriteSolution(INDEX Iter,char *File1,char *str);
int factor(SparseMatrix *Mptr);
void UpdateACvar(VALUETYPE cons,INDEX j,BOOLEAN Limits,BOOLEAN Recover);
void UpdateDCvar(VALUETYPE cons,INDEX j,BOOLEAN Limits);
void UpdateEvector(VALUETYPE cons);
BOOLEAN CheckRlimits(void);
BOOLEAN CheckVlimits(void);
BOOLEAN CheckQlimits(void);
BOOLEAN CheckDClimits(void);
BOOLEAN ChangeSVCmode(void);     /* FACTS */
BOOLEAN ChangeTCSCmode(void);    /* FACTS */
BOOLEAN ChangeSTATCOMmode(void); /* FACTS */
BOOLEAN ChangeDCmode(void);
ACbusData *GetACbus(INDEX N);
void PrintMismatch(VALUETYPE val,INDEX j,INDEX N1);
void DeleteJac(SparseMatrix *Mptr,IntegerVector *P1Row,IntegerVector *P1Col,
               IntegerVector *P2Row,IntegerVector *P2Col);
void WriteJac(void);
int Pflow(int iter,BOOLEAN flagF,BOOLEAN flagD,BOOLEAN flagFirst);
void InitializeLoad(void);

#else
void ErrorStop();
AreaData *ACFunJac();
BOOLEAN DCFunJac();
void UpdateSVCvar();      // FACTS 
void SVCFunJac();         // FACTS 
void UpdateTCSCvar();     // FACTS 
void TCSCFunJac();        // FACTS 
void UpdateSTATCOMvar();  // FACTS 
void STATCOMFunJac();     // FACTS 
int HFunJac();
void ACFunHes();
BOOLEAN DCFunHes();
void SVCFunHes();       // FACTS 
void TCSCFunHes();      // FACTS 
void STATCOMFunHes();   // FACTS 
int factorns();
void repsolp();
VALUETYPE Norm();
void WriteSolution();
int factor();
void UpdateACvar();
void UpdateDCvar();
void UpdateEvector();
BOOLEAN CheckRlimits();
BOOLEAN CheckVlimits();
BOOLEAN CheckQlimits();
BOOLEAN CheckDClimits();
BOOLEAN ChangeSVCmode();     // FACTS 
BOOLEAN ChangeTCSCmode();    // FACTS 
BOOLEAN ChangeSTATCOMmode(); // FACTS 
BOOLEAN ChangeDCmode();
ACbusData *GetACbus();
void PrintMismatch();
void DeleteJac();
void WriteJac();
int Pflow();
void InitializeLoad();
#endif

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX MaxIter,Nac,NacEl,NregPQ,NregV,Ndc,Nslack,Nvolt,Narea,NacVar,Bl,
             Nsvc,Ntcsc,NtcscVar,Nstatcom; /* FACTS */
extern INDEX *ACvar;
extern VALUETYPE *dx,*dF,tol,Tol,Sn,lambda,*x0,*Dx;
extern VALUETYPE K1,K2,MaxdFi,alpha;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern IntegerVector *RowPer,*ColPer;
extern BOOLEAN Acont,PQcont,QRcont,Rcont,PQlim,Tlim,Qlim,Vlim,flagH,flagPoC,flagL,flagR;
extern INDEX *InvRowPerm,*InvColPerm;
extern BOOLEAN *MarkRowPerm,*MarkColPerm;
extern int SD0;

/* --------------- Norm ---------------------- */
#ifdef ANSIPROTO
VALUETYPE Norm(VALUETYPE *Vptr,INDEX N,INDEX *N1)
#else
VALUETYPE Norm(Vptr,N,N1)
VALUETYPE *Vptr;
INDEX N,*N1;
#endif
/* Find the norm (max. value) of a vector. */
{
  VALUETYPE val;
  INDEX i;

  val=-0.1;
  for (i=1;i<=N;i++) if (fabs(Vptr[i])>val) {val=fabs(Vptr[i]); *N1=i;}
  return(val);
}


/* ----------------- GetACbus ------------------------ */
#ifdef ANSIPROTO
ACbusData *GetACbus(INDEX N)
#else
ACbusData *GetACbus(N)
INDEX N;
#endif
{
  ACbusData *ACptr;

  for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
     if (ACptr->N==N) return(ACptr);
  }
  return(NULL);
}


/* ----------------- PrintMismatch ------------------------ */
#ifdef ANSIPROTO
void PrintMismatch(VALUETYPE val,INDEX j,INDEX N1)
#else
void PrintMismatch(val,j,N1)
VALUETYPE val;
INDEX j,N1;
#endif
{
  INDEX k,l,N,N2,N3; 
  INDEX m,n,o,N4,N5,N6; /* FACTS */
  ACbusData *ACptr;

  N2=N3=N4=N5=N6=0; l=NacVar;
  m=NacVar+11*Ndc/2;
  n=m+3*Nsvc;      /* FACTS */
  o=n+NtcscVar;    /* FACTS */
  N=o+7*Nstatcom;  /* FACTS */
  if (((!flagPoC || Jac->n1==N) && !flagH) || N1!=Jac->n1) {
     if (flagPoC && N1>N) N1=N1-N;
     else N=0;
     if (N1>NacVar && N1<=NacVar+11*Ndc/2) for (k=1; k<=Ndc/2; k++){
        N3=k;
        if (N1<=l+11) {N1=N1-l; break;}
        l=l+11;
     }
                              /* FACTS */
     else if (N1>NacVar+11*Ndc/2 && N1<=NacVar+11*Ndc/2+3*Nsvc) for (k=1;k<=Nsvc;k++){
        N4=k;
        if (N1<=m+3){N1=N1-m;break;}
        m=m+3;
     }
     else if (N1>NacVar+11*Ndc/2+3*Nsvc && N1<=NacVar+11*Ndc/2+3*Nsvc+NtcscVar) for (k=1;k<=Ntcsc;k++){
        N5=k;
        if (N1<=n+7){N1=N1-n;break;}
        n=n+7;
     }
     else if (N1>NacVar+11*Ndc/2+3*Nsvc+NtcscVar && N1<=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom) for (k=1;k<=Nstatcom;k++){
        N6=k;
        if (N1<=o+7){N1=N1-o;break;}
        o=o+7;
     }
                              /* END OF FACTS */
     else for(k=1;k<=Nac;k++){
        N2=k;
        if (k+1>Nac || N1<ACvar[k+1]) {N1=N1-ACvar[k]+1; break;}
     }
  } else N=0;
  if (j!=0) fCustomPrint(stderr,"%15s","");
  fCustomPrint(stderr,"Maximum mismatch: %8.4lg  ",val);
  if (N) fCustomPrint(stderr,"PoC-");
  fCustomPrint(stderr,"Equation: %d  ",N1);
  if (N2) {
    ACptr=(ACbusData *) GetACbus(N2);
    if (ACptr!=NULL) fCustomPrint(stderr,"AC bus: %d\n",ACptr->Num);
  }
  else if(N3) fCustomPrint(stderr,"DC link: %d\n",N3);
  else if(N4) fCustomPrint(stderr,"SVC: %d\n",N4);          /* FACTS */
  else if(N5) fCustomPrint(stderr,"TCSC: %d\n",N5);         /* FACTS */
  else if(N6) fCustomPrint(stderr,"STATCOM: %d\n",N6);      /* FACTS */
  else if (flagH) fCustomPrint(stderr,"Continuation Equation\n");
  else fCustomPrint(stderr,"PoC Eigenvector Equation\n");
}

/* ------------------------ DeleteJac --------------------------------- */
#ifdef ANSIPROTO
void DeleteJac(SparseMatrix *Mptr,IntegerVector *P1Row,IntegerVector *P1Col,
               IntegerVector *P2Row,IntegerVector *P2Col)
#else
void DeleteJac(Mptr,P1Row,P1Col,P2Row,P2Col)
SparseMatrix *Mptr;
IntegerVector *P1Row,*P1Col,*P2Row,*P2Col;
#endif
{
  INDEX k;
  SparseMatrixElement *Jptr,*Jptrp;

  for (k=1;k<=Mptr->n1;k++) {
    Jptr=Mptr->RowHead[k];
    while (Jptr!=NULL) {
      Jptrp=Jptr->RowNext;
#ifdef WINDOWS
      delete Jptr;
#else
      free(Jptr);
#endif
      Jptr=Jptrp;
    }
    Mptr->ColHead[k]=Mptr->RowHead[k]=NULL;
    P1Row->p[k]=P1Col->p[k]=0;
    P2Row->p[k]=P2Col->p[k]=0;
  }
}

/* --------------------------- Pflow --------------------------------- */
#ifdef ANSIPROTO
int Pflow(int iter,BOOLEAN flagDCLimits,BOOLEAN flagLimits,BOOLEAN flagFirst)
#else
int Pflow(iter,flagDCLimits,flagLimits,flagFirst)
int iter;
BOOLEAN flagDCLimits,flagLimits,flagFirst;
#endif
/* Power Flow Routine. */
{
  SparseMatrixElement *Jptr;
  AreaData *Aptr;
  int i,m=0,PgMax,PgMaxH;
  INDEX j,k,N,N1,N2,N3,MaxIterRun;
  VALUETYPE MaxdFi1,cons,val;
  BOOLEAN flag=FALSE,flagp=FALSE,flags=FALSE;


  N3=iter;
  if (flagH) MaxIterRun=MaxIter+iter-1;
  else MaxIterRun=MaxIter;
  RowPer=NewRow; ColPer=NewCol;
  NewtonRapson:	
  ACFunJac(Jac,&PgMax,TRUE,FALSE,flagFirst);
  if(DCFunJac(Jac,TRUE,FALSE)) return(-iter);
  SVCFunJac(Jac,TRUE,FALSE);                  /*  FACTS  */
  TCSCFunJac(Jac,TRUE,FALSE);                 /*  FACTS  */
  STATCOMFunJac(Jac,TRUE,FALSE);              /*  FACTS  */
  if (flagH) HFunJac(TRUE,FALSE,NULL,Dx);
  else if (iter!=1 && flagPoC) {
    ACFunHes(TRUE,FALSE);
    if(DCFunHes(TRUE,FALSE)) return(-iter);
    SVCFunHes(TRUE,FALSE);                   /* FACTS  */
    TCSCFunHes(TRUE,FALSE);                  /* FACTS  */
    STATCOMFunHes(TRUE,FALSE);               /* FACTS  */
  }
  MaxdFi=Norm(dF,Jac->n1,&N1);
  fCustomPrint(stderr,"\nIteration: %2d  ",N3);
  PrintMismatch(MaxdFi,0,N1);
  for (i=N3;i<=MaxIterRun;i++){
    N=Jac->n1;
    /* if (ExistParameter('d')) fCustomPrint(stderr,"i %d   MaxdFi %lf   Tol %lf\n",i,MaxdFi,Tol); */
    if ((MaxdFi>Tol) || flagDCLimits) {
      if ((i>1 || flagR || flagL) && !flagDCLimits && (flagLimits || (i==N3))) {
        if(ExistParameter('d')) fCustomPrint(stderr,"Factor Jacobian.\n");
        for(k=1; k<=N; k++) for(Jptr=Jac->RowHead[k];Jptr!=NULL;Jptr->Value=0,Jptr=Jptr->RowNext);
        Aptr=ACFunJac(Jac,&PgMax,TRUE,TRUE,flagFirst);
        if(DCFunJac(Jac,TRUE,TRUE)) return(-iter);
        SVCFunJac(Jac,TRUE,TRUE);                    /* FACTS */
        TCSCFunJac(Jac,TRUE,TRUE);                   /* FACTS */
        STATCOMFunJac(Jac,TRUE,TRUE);                /* FACTS */
        if (flagH) PgMaxH=HFunJac(TRUE,TRUE,Aptr,Dx);
        else if (iter!=1 && flagPoC) {
           ACFunHes(TRUE,TRUE);
           if(DCFunHes(TRUE,TRUE)) return(-iter);
           SVCFunHes(TRUE,FALSE);                   /* FACTS  */
           TCSCFunHes(TRUE,FALSE);                  /* FACTS  */
           STATCOMFunHes(TRUE,FALSE);               /* FACTS  */
        }
        if (PgMax<0 && (!flagH || (flagH && PgMaxH<0)) ) {
           if (Aptr!=NULL) {
              fCustomPrint(stderr,"\nError: Area %d %s does not have any spinning reserves.\n",Aptr->N,Aptr->Name);
              fCustomPrint(stderr,"       Increase the maximum P generation in this area, otherwise\n");
           } else {
              fCustomPrint(stderr,"\nError: The system does not have any spinning reserves.\n");
              fCustomPrint(stderr,"       Increase the maximum P generation in this system, otherwise\n");
           }
           fCustomPrint(stderr,"       the Jacobian matrix becomes singular.\n");
           if (flagFirst) InitializeLoad();
           WriteSolution(--i,TrueParamStr(2),"Pg Max. Problems:");
           stopExecute(1);
        }
        m=factor(Jac);
      }
      if ((i==1 && !flagR && !flagL) || m==WARNINGEXIT || flagDCLimits) {
        if(ExistParameter('d')) fCustomPrint(stderr,"Order and factor Jacobian.\n");
        flagDCLimits=FALSE;
        if (i>1 || flagR || flagL) DeleteJac(Jac,NewRow,NewCol,OldRow,OldCol);
        Aptr=ACFunJac(Jac,&PgMax,TRUE,TRUE,flagFirst);
        if(DCFunJac(Jac,TRUE,TRUE)) return(-iter);
        SVCFunJac(Jac,TRUE,TRUE);                    /* FACTS */
        TCSCFunJac(Jac,TRUE,TRUE);                   /* FACTS */
        STATCOMFunJac(Jac,TRUE,TRUE);                /* FACTS */
        if (flagH) PgMaxH=HFunJac(TRUE,TRUE,Aptr,Dx);
        else if (iter!=1 && flagPoC) {
           ACFunHes(TRUE,TRUE);
           if(DCFunHes(TRUE,TRUE)) return(-iter);
           SVCFunHes(TRUE,FALSE);                   /* FACTS  */
           TCSCFunHes(TRUE,FALSE);                  /* FACTS  */
           STATCOMFunHes(TRUE,FALSE);               /* FACTS  */
        }
        if (PgMax<0 && (!flagH || (flagH && PgMaxH<0)) ) {
           if (Aptr!=NULL) {
              fCustomPrint(stderr,"\nError: Area %d %s does not have any spinning reserves.\n",Aptr->N,Aptr->Name);
              fCustomPrint(stderr,"       Increase the maximum P generation in this area, otherwise\n");
           } else {
              fCustomPrint(stderr,"\nError: The system does not have any spinning reserves.\n");
              fCustomPrint(stderr,"       Increase the maximum P generation in this system, otherwise\n");
           }
           fCustomPrint(stderr,"       the Jacobian matrix becomes singular.\n");
           if (flagFirst) InitializeLoad();
           WriteSolution(--i,TrueParamStr(2),"Pg Max. Problems:");
           stopExecute(1);
        }
        SortRowsColumns(Jac);
        if(factorns(Jac,alpha,RowPartition,ColPartition,NewRow,NewCol,OldRow,OldCol)){
           fCustomPrint(stderr,"*** Singular Jacobian (possible voltage collapse, contol or limit problems).\n");
           fCustomPrint(stderr,"    Try changing the load levels, controls or limits, or use the -F option.\n");
           if (flagFirst) InitializeLoad();
           WriteSolution(--i,TrueParamStr(2),"Singular Jacobian:");
           stopExecute(1);
        }
        SortRowsColumns(Jac);
      }
      fCustomPrint(stderr,"Iteration: %2d  ",i);
      for(j=1;j<=N;j++) dx[j]=dF[j];
      repsolp(Jac,dx,OldRow,NewCol);
      if (m==WARNINGEXIT) SD0=0;
      N2=10; k=j=0;
      while(j<=N2){
        cons= j;
        if(j==0) cons= -1; else cons=1.0/pow(2.0,cons);
        UpdateACvar(cons,j,TRUE,!ExistParameter('G'));
        UpdateDCvar(cons,j,!flag);
        UpdateSVCvar(cons,j);                        /* FACTS */
        UpdateTCSCvar(cons,j);                       /* FACTS */
        UpdateSTATCOMvar(cons,j);                    /* FACTS */
        if (iter!=1 && flagPoC) UpdateEvector(cons);
        ACFunJac(Jac,&PgMax,TRUE,FALSE,flagFirst);
        flags=DCFunJac(Jac,TRUE,FALSE);
        SVCFunJac(Jac,TRUE,FALSE);                   /* FACTS */
        TCSCFunJac(Jac,TRUE,FALSE);                  /* FACTS */
        STATCOMFunJac(Jac,TRUE,FALSE);               /* FACTS */
        if(!flags) {
            if (flagH) HFunJac(TRUE,FALSE,NULL,Dx);
            else if (iter!=1 && flagPoC) {
              ACFunHes(TRUE,FALSE);
              flags=DCFunHes(TRUE,FALSE);
              SVCFunHes(TRUE,FALSE);                   /* FACTS  */
              TCSCFunHes(TRUE,FALSE);                  /* FACTS  */
              STATCOMFunHes(TRUE,FALSE);               /* FACTS  */
           }
        }
        if (!flags) {
           val=Norm(dF,N,&N1);
           PrintMismatch(val,k++,N1);
           if (j==10 && val>2*MaxdFi) N2=20;
           else if (MaxdFi>val) {
              MaxdFi1=MaxdFi;
              MaxdFi=val;
              break;
           }
        }
        j++;
      }
      /* if (ExistParameter('d')) fCustomPrint(stderr,"j %d   MaxdFi %lf   Tol %lf   MaxdFi1 %lf   tol %lf\n",j,MaxdFi,Tol,MaxdFi1,tol); */
      if ((MaxdFi>Tol) && (j>N2 || (fabs(MaxdFi1-MaxdFi)/MaxdFi)<=tol)) {
        if (flagL) return(-(++i));
        if (!flagR) {
          flagp=CheckRlimits();
          /*  Apply Q limits after convergence
          if (!ExistParameter('G') && !flagp) flagp=CheckQlimits();
          */
          if (!flagp) flagp=CheckVlimits();
          if (!flagp) flagp=CheckQlimits();
          if (!flagp) flagp=flag=CheckDClimits();
          if (!flagp) flagp=flagDCLimits=ChangeSVCmode();       /* FACTS */
          if (!flagp) flagp=flagDCLimits=ChangeTCSCmode();      /* FACTS */
          if (!flagp) flagp=flagDCLimits=ChangeSTATCOMmode();  /* FACTS */
        } else {
          flagp=CheckRlimits();
          if (!flagp) flagp=CheckVlimits(); else CheckVlimits();
          if (!flagp) flagp=CheckQlimits(); else CheckQlimits();
          if (!flagp) flagp=flagDCLimits=ChangeDCmode(); else flagDCLimits=ChangeDCmode();
          if (!flagp) flagp=flagDCLimits=ChangeSVCmode(); else flagDCLimits=ChangeSVCmode();          /* FACTS */
          if (!flagp) flagp=flagDCLimits=ChangeTCSCmode(); else flagDCLimits=ChangeTCSCmode();        /* FACTS */
          if (!flagp) flagp=flagDCLimits=ChangeSTATCOMmode(); else flagDCLimits=ChangeSTATCOMmode();  /* FACTS */
          flagL=TRUE;
       }
       if(!flagp) {
          if (flagR) return(-(++i));
          MaxdFi=val;
          fCustomPrint(stderr,"\n *** The case diverges (possible voltage collapse or AC/DC/FACTS control\n");
          fCustomPrint(stderr,"     problems).  Try changing the load levels or AC/DC/FACTS controls, or\n");
          fCustomPrint(stderr,"     use the -F option or decrease the tolerance between two\n");
          fCustomPrint(stderr,"     consecutive iterations with the -t option.\n");
          if (flagFirst) InitializeLoad();
          WriteSolution(i,TrueParamStr(2),"Divergence:");
          stopExecute(1);
        }
      }
    } else break;
  }
  /* if (ExistParameter('d')) fCustomPrint(stderr,"i %d   MaxdFi %lf   Tol %lf\n",i,MaxdFi,Tol); */
  if (i>MaxIterRun && MaxdFi>Tol) {
     if (flagR) return(-i);
     fCustomPrint(stderr,"\n *** The case has not been solved (possible voltage collapse, AC/DC/FACTS control\n");
     fCustomPrint(stderr,"     problems, or too few iterations).  Try running the case using the -F\n");
     fCustomPrint(stderr,"     option, or change load levels, AC/DC/FACTS controls, or increase the maximum\n");
     fCustomPrint(stderr,"     number of iterations with the -M option.\n");
     if (flagFirst) InitializeLoad();
     WriteSolution(--i,TrueParamStr(2),"Unsolved case:");
     stopExecute(1);
  }
  /*    Apply Q limits after convergence
  N3=i;
  if (ExistParameter('G') && !flagR) {
     flagp=CheckQlimits();
     if (flagp) goto NewtonRapson;
  }
  */
  if (flagFirst) InitializeLoad();
  return(i);
}
