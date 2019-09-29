#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

/* AC/DC Power Flow Program.

*/

#include <stdlib.h>
#include "pfwstdio.h"
#include <math.h>
#include <string.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"

#ifdef WINDOWS
#include "Win UWPflow.h"
#include "GraphDLG.h"
#endif

#ifdef ANSIPROTO
void ErrorStop(char *Msg);
void DCinit(void);
void SVCinit(void);     /* FACTS */
void TCSCinit(void);    /* FACTS */
void STATCOMinit(void); /* FACTS */
void Jacobian(void);
void ReadData(char *Name);
BOOLEAN ReadInit(void);
void WriteSolution(INDEX Iter,char *File1,char *str);
int  Pflow(int iter,BOOLEAN flagF,BOOLEAN flagD,BOOLEAN flagP);
int Homot(void);
int PoCPoint(void);
#ifdef WINDOWS
  int pfw_main(int argc, char **argv);
  void initFile();
#else
  int main(int argc,char **argv);
#endif
void InitializeLoad(void);
int Evector(int M,int iter,VALUETYPE tol,BOOLEAN RightEvector,VALUETYPE *EigenValue);
void PrintLeftEvector(INDEX N,FILE *Out);
void PrintDirection(char Option,VALUETYPE *vector,VALUETYPE Max);
void CleanUp(void);
void DeleteJac(SparseMatrix *Mptr,IntegerVector *P1Row,IntegerVector *P1Col,
               IntegerVector *P2Row,IntegerVector *P2Col);
#else
void ErrorStop();
void DCinit();
void SVCinit();      // FACTS 
void TCSCinit();     // FACTS 
void STATCOMinit();  // FACTS 
void Jacobian();
void ReadData();
BOOLEAN ReadInit();
void WriteSolution();
int  Pflow();
int Homot();
int PoCPoint();
int main();
void InitializeLoad();
int Evector();
void PrintLeftEvector();
void PrintDirection();
void CleanUp();
void DeleteJac();
#endif

/* ------- Global Variables ------ */
Data *dataPtr=NULL;
SparseMatrix *Jac=NULL;
INDEX MaxIter,Nac,NacEl,NregPQ,NregV,Ngen,Ndc,
      Nslack,Nvolt,Narea,NacVar,Bl,*ACvar,NZvolt,NXvolt,
      Nsvc,Ntcsc,NtcscVar,Nstatcom;  /* FACTS */
ACbusData *BlPtr;
VALUETYPE *dx,*dF,tol,Tol,Sn,lambda;
VALUETYPE K1,K2,K3=PI/180.0,MaxdFi,alpha;
IntegerVector *NewRow=NULL,*OldRow=NULL,*NewCol=NULL,*OldCol=NULL,
              *RowPartition=NULL,*ColPartition=NULL;
BOOLEAN Acont,PQcont,QRcont,Rcont,Xcont,
        PQlim,Tlim,Qlim,Vlim,Elim,Ilim,Zlim,Xlim,flag2Vcontrol,
        flagH,flagPoC,flagL,flagR,flagPgMax,flagSmax,flagReducedContinuation,
        flagVloads,flagKdirection,flagBS;
extern VALUETYPE *x0,*x0p,*Dx;
extern BOOLEAN *DxZero;
extern AClist *Vlist,*Vlistp;

#ifdef WINDOWS
extern CString dir;
extern GraphDLG* GraphDlg;
#endif

/* ------------------- DCinit -------------------- */
#ifdef ANSIPROTO
void DCinit(void)
#else
void DCinit()
#endif
/*    Intialize the DC system variables. */
{
  ACbusData *ACptrR,*ACptrI;
  DCbusData *DCptrR,*DCptrI;
  VALUETYPE KVr,KVi,Vr,Vi,ar,ai,cosar,cosgi,Xcr,Xci,Id,Rd,Ld,Vn,In;

  for(DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if (!strcmp(DCptrR->Type,"R")){
      KVr=DCptrR->Vn;  KVi=DCptrI->Vn;
      ACptrR=DCptrR->AC;
      ACptrI=DCptrI->AC;
      Vr=KVr*ACptrR->V;  Vi=KVi*ACptrI->V;
      ar=DCptrR->Tap*DCptrR->Ntrf;  ai=DCptrI->Tap*DCptrI->Ntrf;
      Xcr=DCptrR->Xc;    Xci=DCptrI->Xc;
      cosar=cos(DCptrR->Alfa*K3);  cosgi=cos(DCptrI->Gamma*K3);
      Id=DCptrR->Id/1000;
      Rd=DCptrR->Rd;   Ld=DCptrR->Ld;
      if (DCptrR->Vd<=0) DCptrR->Vd=fabs(K1*ar*Vr*cosar-K2*Xcr*Id);
      if (DCptrI->Vd<=0) DCptrI->Vd=fabs(K1*ai*Vi*cosgi-K2*Xci*Id);
      if (Id<=0) {
        Id=fabs(DCptrR->Vd-DCptrI->Vd)/Rd;
        DCptrR->Id=DCptrI->Id=Id;
      }
      if (DCptrR->P<=0) DCptrR->P=DCptrR->Vd*Id;
      if (DCptrI->P<=0) DCptrI->P=DCptrI->Vd*Id;
      if (DCptrR->Q<=0) DCptrR->Q=DCptrR->Vd*Id*tan(acos(cosar));
      if (DCptrI->Q<=0) DCptrI->Q=DCptrI->Vd*Id*tan(acos(cosgi));
  /* ------------- Per unit conversion ----------------------------- */
      Vn=KVr;  In=Sn/Vn;
      DCptrR->Vn=DCptrI->Vn=Vn;
      DCptrR->Rd=DCptrI->Rd=Rd*In/Vn;
      DCptrR->Ld=DCptrI->Ld=Ld*In/Vn;
      DCptrR->Id=DCptrI->Id=Id/In;
      DCptrR->Vd=DCptrR->Vd/Vn;
      DCptrI->Vd=DCptrI->Vd/Vn;
      DCptrR->Xc=DCptrR->Xc*In/Vn;
      DCptrI->Xc=DCptrI->Xc*In/Vn;
      DCptrR->Alfa=DCptrR->Alfa*K3;
      DCptrI->Alfa=DCptrI->Alfa*K3;
      DCptrR->Gamma=DCptrR->Gamma*K3;
      DCptrI->Gamma=DCptrI->Gamma*K3;
      DCptrR->AlfaMin=DCptrR->AlfaMin*K3;
      DCptrI->AlfaMin=DCptrI->AlfaMin*K3;
      DCptrR->AlfaMax=DCptrR->AlfaMax*K3;
      DCptrI->AlfaMax=DCptrI->AlfaMax*K3;
      DCptrR->GammaMin=DCptrR->GammaMin*K3;
      DCptrI->GammaMin=DCptrI->GammaMin*K3;
      DCptrI->Ntrf=DCptrI->Ntrf*KVi/KVr;
      DCptrR->P= -DCptrR->P/Sn;
      DCptrI->P=DCptrI->P/Sn;
      DCptrR->Q= -DCptrR->Q/Sn;
      DCptrI->Q= -DCptrI->Q/Sn;
      DCptrR->MVA=sqrt(DCptrR->P*DCptrR->P+DCptrR->Q*DCptrR->Q);
      DCptrI->MVA=sqrt(DCptrI->P*DCptrI->P+DCptrI->Q*DCptrI->Q);
    }
  }
}



/* --------------------------- InitializeLoad ------------------------------ */
#ifdef ANSIPROTO
void InitializeLoad(void)
#else
void InitializeLoad()
#endif
/* Initialize load parameters for OH load model. */
{
  ACbusData *ACptr;
  VALUETYPE Pn,Qn,Pz,Qz,val;

  for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
         if (ACptr->Pz==0 && ACptr->a==0) ACptr->Pn=ACptr->Pl;
         else {
                Pn=ACptr->Pn;
                Pz=1-Pn;
                val=ACptr->Pl/(Pn*pow(ACptr->V,ACptr->a)+Pz*ACptr->V*ACptr->V);
                ACptr->Pn=Pn*val;  ACptr->Pz=Pz*val;
         }
         if (ACptr->Qz==0 && ACptr->b==0) ACptr->Qn=ACptr->Ql;
         else {
                Qn=ACptr->Qn;
                Qz=1-Qn;
                val=ACptr->Ql/(Qn*pow(ACptr->V,ACptr->b)+Qz*ACptr->V*ACptr->V);
                ACptr->Qn=Qn*val;  ACptr->Qz=Qz*val;
         }
  }
}


/* --------------------------- Main Program  ------------------------------ */
#ifdef WINDOWS
  int pfw_main(int argc, char **argv)
#else
#ifdef ANSIPROTO
  int main(int argc,char **argv)
#else
  int main(argc,argv)
  int argc;
  char **argv;
#endif
#endif
/* Main Routine. */
{
  char *Name;
  BOOLEAN flagv=FALSE,flagD=FALSE,PrintEvalue=FALSE;
  VALUETYPE V,EigenValue;
  INDEX N1;
  FILE *Out;
  int i;

  SetArguments(argc,argv);

#ifdef WINDOWS
  if (setjmp(exit_main)) return(1);	 
  
#endif

  if (HelpRequested())
  {
	  //display error message
	  ErrorStop("");
	  //exit from main
	  return(1);
  }
  if (ExistParameter('l')) {
         Name=NameParameter('l');
#ifdef WINDOWS
		 char charDir[300];
		 strcpy_s(charDir, dir.GetBuffer(dir.GetLength()));
         if (!NullName(Name)) freopen(strcat(charDir, Name),"wt",stderr);
#else
		 if (!NullName(Name)) freopen(Name,"wt",stderr);
#endif
  }

  fCustomPrint(stderr,"UW Continuation Power Flow (c)1992,1996,1999, 2006 C. Canizares, F. Alvarado and S. Zhang.\n");
  if (ExistParameter('A') || ExistParameter('6')) Acont=FALSE;
  else Acont=TRUE;
  if (ExistParameter('P')) PQcont=FALSE;
  else PQcont=TRUE;
  if (ExistParameter('Q') && NullName(NameParameter('Q'))) QRcont=FALSE;
  else QRcont=TRUE;
  if (ExistParameter('Q') && !strcmp(NameParameter('Q'),"X")) Xcont=FALSE;
  else Xcont=TRUE;
  if (ExistParameter('R')) Rcont=FALSE;
  else Rcont=TRUE;
  if (ExistParameter('N')) Acont=PQcont=QRcont=Rcont=Xcont=FALSE;
  if (ExistParameter('a')) Tlim=FALSE;
  else Tlim=TRUE;
  if (ExistParameter('P')||ExistParameter('p')) PQlim=FALSE;
  else PQlim=TRUE;
  if (ExistParameter('q') && NullName(NameParameter('q'))) Qlim=FALSE;
  else Qlim=TRUE;
  if (ExistParameter('q') && !strcmp(NameParameter('q'),"z")) Zlim=FALSE;
  else Zlim=TRUE;
  if (ExistParameter('q') && !strcmp(NameParameter('q'),"x")) Xlim=FALSE;
  else Xlim=TRUE;
  if (ExistParameter('r')) Vlim=FALSE;
  else Vlim=TRUE;
  if (ExistParameter('3')) {
    if (ExistParameter('4')) Elim=FALSE;
    else Elim=TRUE;
    if (ExistParameter('5')) Ilim=FALSE;
    else Ilim=TRUE;
  } else Elim=Ilim=FALSE;
  flagPgMax=ExistParameter('X');
  flagSmax=ExistParameter('9');
  if (ExistParameter('n')) {Tlim=PQlim=Qlim=Vlim=Elim=Ilim=Xlim=Zlim=FALSE; flagPgMax=flagSmax=TRUE;}
  flagVloads=flagKdirection=FALSE;

  //Call ReadData with the name of the inpute file
  ReadData(TrueParamStr(1));
  Tol=1e-4;
  RealParameter('T',&Tol,1e-10,1e-3);
  tol=0.0001;
  RealParameter('t',&tol,Tol,0.1);
  alpha=0.01;
  RealParameter('F',&alpha,0.001,1.);
  MaxIter=(INDEX) IntegerParameter('M',MaxIter,0,1000);
  if (ExistParameter('d')) {
         fCustomPrint(stderr,"Tol %lf   tol %lf   alpha %lf\n",Tol,tol,alpha);
  }
  flagL=flagR=FALSE;
  if ((ExistParameter('H') || ExistParameter('c')) && (!NullName(NameParameter('K')) || flagKdirection)) flagH=TRUE;
  else flagH=FALSE;
  if (ExistParameter('C') && (!NullName(NameParameter('K')) || flagKdirection)) flagPoC=TRUE;
  else flagPoC=FALSE;
  if ((ExistParameter('D') && (!NullName(NameParameter('D')))) || flagVloads) flagD=TRUE;
  else flagD=FALSE;
  flagv=ReadInit();
  flagReducedContinuation=FALSE;
  flagBS = ExistParameter('x');
  K1=sqrt(2.0)*3/PI;
  K2=3/PI;
  lambda=0;
  DCinit();
  SVCinit();      /* FACTS */
  TCSCinit();     /* FACTS */
  STATCOMinit();  /* FACTS */
  Jacobian();
  if (flagH) {
         i=Homot();
         if(i<0) {
                i= -i;
                fCustomPrint(stderr,"\n *** The DC equations have a square root of a negative number.\n");
                fCustomPrint(stderr,"     Try changing the DC controls.\n");
                fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
                WriteSolution(--i,TrueParamStr(2),"DC problems:");
                stopExecute(1);
         }
         fCustomPrint(stderr,"**** Voltage Profile Case Solved **** \n\n");
         WriteSolution(--i,TrueParamStr(2),"U.E.P. Solution:");
         if (ExistParameter('y') && !NullName(NameParameter('y'))) {
                if(ExistParameter('d')) fCustomPrint(stderr,"Write left e-vector for base case (see file 'evect.dat').\n");
                Evector(40,0,0.00001,FALSE,&EigenValue);
                fCustomPrint(stderr,"Minimum |e_value| -> %-10.6lg\n\n",EigenValue);
                PrintEvalue=TRUE;
                Out=OpenOutput(NameParameter('y'));
                N1=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;   /* FACTS */
                PrintLeftEvector(N1,Out);
         }
         if (ExistParameter('Y') && !NullName(NameParameter('Y'))){
                if(ExistParameter('d')) fCustomPrint(stderr,"Write right e-vector for base case (see file 'evect.dat').\n");
                Evector(40,0,0.00001,TRUE,&EigenValue);
                if (!PrintEvalue) fCustomPrint(stderr,"Minimum |e_value| -> %-10.6lg\n\n",EigenValue);
                PrintDirection('Y',x0,1.0);
         }
  } else if (flagPoC) {
         i=PoCPoint();
         if(i<0) {
                i= -i;
                fCustomPrint(stderr,"\n *** The DC equations have a square root of a negative number.\n");
                fCustomPrint(stderr,"     Try changing the DC controls.\n");
                fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
                WriteSolution(--i,TrueParamStr(2),"DC problems:");
                stopExecute(1);
         }
         fCustomPrint(stderr,"**** Point of Collapse Case Solved ****   ");
         fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
         WriteSolution(--i,TrueParamStr(2),"PoC Solution:");
  } else if (flagv) {
         i=1;
         V=BlPtr->V;
         if (flagD) {
                i=Pflow(i,FALSE,TRUE,TRUE);
                if(i<0) {
                  i= -i;
                  fCustomPrint(stderr,"\n *** The DC equations have a square root of a negative number.\n");
                  fCustomPrint(stderr,"     Try changing the DC controls.\n");
                  WriteSolution(--i,TrueParamStr(2),"DC problems:");
                  stopExecute(1);
                }
                fCustomPrint(stderr,"**** Base Case Solved (to calculate OH load parameters) ****\n\n");
         } else InitializeLoad();
         RealParameter('L',&lambda,-1e6,1e6);
         strcpy_s(BlPtr->Type,"BL");
         if(BlPtr->Area!=NULL && BlPtr->Area->Slack==BlPtr) strcat_s(BlPtr->Type,"A");
         Bl=BlPtr->N;
         BlPtr->V=V;
         BlPtr->Cont=NULL;
         i=Pflow(i,flagD,TRUE,FALSE);
         if(i<0) {
                i= -i;
                fCustomPrint(stderr,"\n *** The DC equations have a square root of a negative number.\n");
                fCustomPrint(stderr,"     Try changing the DC controls.\n");
                fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
                WriteSolution(--i,TrueParamStr(2),"DC problems:");
                stopExecute(1);
         }
         fCustomPrint(stderr,"**** Voltage/Lambda Case Solved ****    ");
         fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
         WriteSolution(--i,TrueParamStr(2),"Voltage/Lambda Solution:");
  } else if (ExistParameter('L') && (ExistParameter('K') || flagKdirection)) {
         i=1;
         if (ExistParameter('b') || flagD ) {
                i=Pflow(i,FALSE,TRUE,TRUE);
                if(i<0) {
                  i= -i;
                  fCustomPrint(stderr,"\n *** The DC equations have a square root of a negative number.\n");
                  fCustomPrint(stderr,"     Try changing the DC controls.\n");
                  WriteSolution(--i,TrueParamStr(2),"DC problems:");
                  stopExecute(1);
                }
                fCustomPrint(stderr,"**** Base Case Solved (to initialize power flow) ****\n\n");
         } else InitializeLoad();
         RealParameter('L',&lambda,-1e6,1e6);
         if (lambda==0 || (NullName(NameParameter('K')) && !flagKdirection)) {
                fCustomPrint(stderr,"***Warning: The program has detected the -L option but either lambda is zero\n");
                fCustomPrint(stderr,"            or there is no gen./load variations defined.\n");
                if (!ExistParameter('b') && !flagD  ) {
                  fCustomPrint(stderr,"            The program will just solve the base case.\n");
                  i=Pflow(i,FALSE,TRUE,TRUE);
                  if(i<0) {
                         i= -i;
                         fCustomPrint(stderr,"\n *** The DC equations have a square root of a negative number.\n");
                         fCustomPrint(stderr,"     Try changing the DC controls.\n");
                         WriteSolution(--i,TrueParamStr(2),"DC problems:");
                         stopExecute(1);
                  }
                  fCustomPrint(stderr,"**** Base Case Solved ****\n\n");
                }
         } else {
                i=Pflow(i,flagD,TRUE,FALSE);
                if(i<0) {
                  i= -i;
                  fCustomPrint(stderr,"\n *** The DC equations have a square root of a negative number.\n");
                  fCustomPrint(stderr,"     Try changing the DC controls.\n");
                  fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
                  WriteSolution(--i,TrueParamStr(2),"DC problems:");
                  stopExecute(1);
                }
                fCustomPrint(stderr,"**** Lambda Case Solved ****    ");
                fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
         }
         WriteSolution(--i,TrueParamStr(2),"Lambda Solution:");
         if (ExistParameter('y') && !NullName(NameParameter('y'))) {
                N1=NacVar+11*Ndc/2+1+3*Nsvc+NtcscVar+7*Nstatcom;   /* FACTS */
#ifdef WINDOWS
                x0= new VALUETYPE[N1];
#else
                x0=(VALUETYPE *) calloc(N1,sizeof(VALUETYPE));
                if (x0==NULL) {ErrorHalt("Insufficient memory to allocate approx. left e-vector."); stopExecute(ERROREXIT);}
#endif
                if(ExistParameter('d')) fCustomPrint(stderr,"Write left e-vector for base case (see file 'evect.dat').\n");
                Evector(40,0,0.00001,FALSE,&EigenValue);
                fCustomPrint(stderr,"Minimum |e_value| -> %-10.6lg\n\n",EigenValue);
                PrintEvalue=TRUE;
                Out=OpenOutput(NameParameter('y'));
                PrintLeftEvector(N1,Out);
         }
         if (ExistParameter('Y') && !NullName(NameParameter('Y'))){
                N1=NacVar+11*Ndc/2+1+3*Nsvc+NtcscVar+7*Nstatcom;   /* FACTS */
                if (x0 == NULL) {
#ifdef WINDOWS
                  x0= new VALUETYPE[N1];
#else
                  x0=(VALUETYPE *) calloc(N1,sizeof(VALUETYPE));
                  if (x0==NULL) {ErrorHalt("Insufficient memory to allocate approx. right e-vector."); stopExecute(ERROREXIT);}
#endif
                }
                if(ExistParameter('d')) fCustomPrint(stderr,"Write right e-vector for base case (see file 'evect.dat').\n");
                Evector(40,0,0.00001,TRUE,&EigenValue);
                if (!PrintEvalue) fCustomPrint(stderr,"Minimum |e_value| -> %-10.6lg\n\n",EigenValue);
                PrintDirection('Y',x0,1.0);
         }
  } else {
         i=Pflow(1,FALSE,TRUE,TRUE);
         if(i<0) {
                i= -i;
                fCustomPrint(stderr,"\n *** The DC equations have a square root of a negative number.\n");
                fCustomPrint(stderr,"     Try changing the DC controls.\n");
                WriteSolution(--i,TrueParamStr(2),"DC problems:");
                stopExecute(1);
         }
         fCustomPrint(stderr,"**** Base Case Solved ****\n\n");
         WriteSolution(--i,TrueParamStr(2),"Base Solution:");
         if (ExistParameter('y') && !NullName(NameParameter('y'))) {
                N1=NacVar+11*Ndc/2+1+3*Nsvc+NtcscVar+7*Nstatcom;   /* FACTS */
#ifdef WINDOWS
                x0= new VALUETYPE[N1];
#else
                x0=(VALUETYPE *) calloc(N1,sizeof(VALUETYPE));
                if (x0==NULL) {ErrorHalt("Insufficient memory to allocate approx. left e-vector."); stopExecute(ERROREXIT);}
#endif
                if(ExistParameter('d')) fCustomPrint(stderr,"Write left e-vector for base case (see file 'evect.dat').\n");
                Evector(40,0,0.00001,FALSE,&EigenValue);
        fCustomPrint(stderr,"Minimum |e_value| -> %-10.6lg\n\n",EigenValue);
        PrintEvalue=TRUE;
                Out=OpenOutput(NameParameter('y'));
                PrintLeftEvector(N1,Out);
         }
         if (ExistParameter('Y') && !NullName(NameParameter('Y'))){
                N1=NacVar+11*Ndc/2+1+3*Nsvc+NtcscVar+7*Nstatcom;   /* FACTS */
                if (x0 == NULL) {
#ifdef WINDOWS
                  x0= new VALUETYPE[N1];
#else
                  x0=(VALUETYPE *) calloc(N1,sizeof(VALUETYPE));
                  if (x0==NULL) {ErrorHalt("Insufficient memory to allocate approx. right e-vector."); stopExecute(ERROREXIT);}
#endif
                }
                if(ExistParameter('d')) fCustomPrint(stderr,"Write right e-vector for base case (see file 'evect.dat').\n");
                Evector(40,0,0.00001,TRUE,&EigenValue);
                if (!PrintEvalue) fCustomPrint(stderr,"Minimum |e_value| -> %-10.6lg\n\n",EigenValue);
                PrintDirection('Y',x0,1.0);
         }
  }
#ifdef WINDOWS

  CleanUp();
#endif
  return(0);
}


/* -------------------- Clean up all data structures ----------------------*/
#ifdef WINDOWS
void CleanUp()
/* Clean up all Pflow data structures */
{

  INDEX k;
  SparseMatrixElement *Jptr,*Jptrp;
  AreaData *Aptr,*Aptrp;
  AClist *LACptr,*LACptrp;
  DClist *LDCptr,*LDCptrp;
  ElementList *LEptr,*LEptrp;
  ACbusData *ACptr,*ACptrp;
  SVClist *LSVCptr,*LSVCptrp;                 /* FACTS */
  TCSClist *LTCSCptr,*LTCSCptrp;              /* FACTS */
  STATCOMlist *LSTATCOMptr,*LSTATCOMptrp;     /* FACTS */
  DCbusData *DCptr,*DCptrp;
  ElementData *Eptr,*Eptrp;
  SVCbusData *SVCptr,*SVCptrp;                /* FACTS */
  TCSCbusData *TCSCptr,*TCSCptrp;             /* FACTS */
  STATCOMbusData *STATCOMptr,*STATCOMptrp;    /* FACTS */

  if (Jac!=NULL) {
    for (k=1;k<=Jac->n1;k++) {
      Jptr=Jac->RowHead[k];
      while (Jptr!=NULL) {
        Jptrp=Jptr->RowNext;
        delete Jptr;
        Jptr=Jptrp;
      }
    }
    delete[] Jac->ColHead;
    delete[] Jac->RowHead;
    delete Jac;
    delete[] NewRow->p;
    delete NewRow;
    delete[] NewCol->p;
    delete NewCol;
    delete[] OldRow->p;
    delete OldRow;
    delete[] OldCol->p;
    delete OldCol;
    delete[] RowPartition->p;
    delete RowPartition;
    delete[] ColPartition->p;
    delete ColPartition;
    Jac=NULL;
  }
  if (ACvar!=NULL)  { delete[] ACvar;  ACvar=NULL; }
  if (Dx!=NULL) { delete[] Dx;  Dx=NULL; }
  if (x0!=NULL) { delete[] x0;  x0=NULL; }
  if (x0p!=NULL) { delete[] x0p;  x0p=NULL; }
  if (dx!=NULL) { delete[] dx;  dx=NULL; }
  if (dF!=NULL) {  delete[] dF;  dF=NULL; }
  if (DxZero!=NULL) {  delete[] DxZero;  DxZero=NULL; }
  if (Vlist!=NULL) {
    if (Vlist!=Vlistp) {
      for(LACptr=Vlist;LACptr!=NULL;LACptrp=LACptr->Next,delete LACptr,LACptr=LACptrp);
      for(LACptr=Vlistp;LACptr!=NULL;LACptrp=LACptr->Next,delete LACptr,LACptr=LACptrp);
    } else
      for(LACptr=Vlist;LACptr!=NULL;LACptrp=LACptr->Next,delete LACptr,LACptr=LACptrp);
    Vlist=Vlistp=NULL;
  }
  if (dataPtr!=NULL) {
    for (Aptr=dataPtr->Area;Aptr!=NULL;){
      Aptrp=Aptr->Next;
      for(LACptr=Aptr->AC;LACptr!=NULL;LACptrp=LACptr->Next,delete LACptr,LACptr=LACptrp);
      for(LDCptr=Aptr->DC;LDCptr!=NULL;LDCptrp=LDCptr->Next,delete LDCptr,LDCptr=LDCptrp);
      for(LEptr=Aptr->Elem;LEptr!=NULL;LEptrp=LEptr->Next,delete LEptr,LEptr=LEptrp);
      delete Aptr;
      Aptr=Aptrp;
    }
    for(ACptr=dataPtr->ACbus;ACptr!=NULL;) {
      ACptrp=ACptr->Next;
      for(LEptr=ACptr->Reg;LEptr!=NULL;LEptrp=LEptr->Next,delete LEptr,LEptr=LEptrp);
      for(LEptr=ACptr->Elem;LEptr!=NULL;LEptrp=LEptr->Next,delete LEptr,LEptr=LEptrp);
      for(LDCptr=ACptr->DC;LDCptr!=NULL;LDCptrp=LDCptr->Next,delete LDCptr,LDCptr=LDCptrp);
      for(LSVCptr=ACptr->SVC;LSVCptr!=NULL;LSVCptrp=LSVCptr->Next,delete LSVCptr,LSVCptr=LSVCptrp);                                  /* FACTS */
      for(LTCSCptr=ACptr->TCSC;LTCSCptr!=NULL;LTCSCptrp=LTCSCptr->Next,delete LTCSCptr,LTCSCptr=LTCSCptrp);                          /* FACTS */
      for(LSTATCOMptr=ACptr->STATCOM;LSTATCOMptr!=NULL;LSTATCOMptrp=LSTATCOMptr->Next,delete LSTATCOMptr,LSTATCOMptr=LSTATCOMptrp);  /* FACTS */
      for(LACptr=ACptr->ContBus;LACptr!=NULL;LACptrp=LACptr->Next,delete LACptr,LACptr=LACptrp);
      if (ACptr->Gen!=NULL) delete ACptr->Gen;
      delete ACptr;
      ACptr=ACptrp;
    }
    for(DCptr=dataPtr->DCbus;DCptr!=NULL;) {
      DCptrp=DCptr->Next;
      delete DCptr;
      DCptr=DCptrp;
    }
    for(Eptr=dataPtr->Element;Eptr!=NULL;) {
      Eptrp=Eptr->Next;
      delete Eptr;
      Eptr=Eptrp;
    }

                           /* FACTS */
    for(SVCptr=dataPtr->SVCbus;SVCptr!=NULL;) {
      SVCptrp=SVCptr->Next;
      delete SVCptr;
      SVCptr=SVCptrp;
    }
    for(TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;) {
      TCSCptrp=TCSCptr->Next;
      delete TCSCptr;
      TCSCptr=TCSCptrp;
    }
    for(STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;) {
      STATCOMptrp=STATCOMptr->Next;
      delete STATCOMptr;
      STATCOMptr=STATCOMptrp;
    }
                        /* END FACTS */

    for(LACptr=dataPtr->KGbus;LACptr!=NULL;) {
      LACptrp=LACptr->Next;
      delete LACptr;
      LACptr=LACptrp;
    }
    delete dataPtr;
    dataPtr=NULL;

	
  }
  //close standard streams
  fclose(stdout);
  fclose(stderr);
  //fclose(stdin);

  
  //update graph
  //set range to around the max and min points on the graph
  GraphDlg->m_GraphCtrl.SetRange(GraphDlg->minX*0.9, GraphDlg->maxX*1.2, GraphDlg->minY*0.9, GraphDlg->maxY*1.1);
	  
  for (int i = 1; i<=GraphDlg->m_GraphCtrl.GetAnnoCount(); i++){
	//set the legend
	  GraphDlg->m_GraphCtrl.SetAnnotation(i-1);
	  GraphDlg->m_GraphCtrl.SetAnnoLabelX(GraphDlg->maxX*1.05);
	  GraphDlg->m_GraphCtrl.SetAnnoLabelY(GraphDlg->maxY * 1.1 *(1-i*0.1) + GraphDlg->minY* 0.9*0.1*i);
  }

  
  //set the x axis label
  GraphDlg->m_GraphCtrl.SetXLabel("L.F.");


}
#endif

