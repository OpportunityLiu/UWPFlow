#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

/* Point of Collapse: -Build the AC/DC Hessian.
                      -Find initial guess for the left e-vector,

        Main. */

#include "pointl.h"


#ifdef ANSIPROTO
void InitializeLoad(void);
void PrintLeftEvector(INDEX N,FILE *Out);
#else
void InitializeLoad();
void PrintLeftEvector();
#endif


extern BOOLEAN flagVloads;

/* -----------------UpdateEvector ---------------------------- */
#ifdef ANSIPROTO
void UpdateEvector(VALUETYPE cons)
#else
void UpdateEvector(cons)
VALUETYPE cons;
#endif
{
  INDEX i,N;

  N=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;   /* FACTS */
  for(i=1;i<=N;i++) x0[i]=x0[i]+cons*dx[i+N];
  Dparam=Dparam+cons*dx[Jac->n1];
  lambda=param0+Dparam;
}

/* ----------------- Evector ---------------------------- */
#ifdef ANSIPROTO
int Evector(int M,int iter,VALUETYPE tol,BOOLEAN RightEvector,VALUETYPE *EigenValue)
#else
int Evector(M,iter,tol,RightEvector,EigenValue)
int M,iter;
VALUETYPE tol;
BOOLEAN RightEvector;
VALUETYPE *EigenValue;
#endif
/* Find an inital guess for the left eigenvector. */
{
  AreaData *Aptr;
  VALUETYPE val,Nx0;
  INDEX i,N,count;
  int PgMax;
  FILE *Out;

  if(ExistParameter('d')) fCustomPrint(stderr,"Order and factor Jacobian.\n");
  DeleteJac(Jac,NewRow,NewCol,OldRow,OldCol);
  N=Jac->n2=Jac->n1=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;    /* FACTS */
  NewRow->N=NewCol->N=N;
  OldRow->N=OldCol->N=N;
  RowPartition->N=ColPartition->N=1;
  RowPartition->p[1]=ColPartition->p[1]=N;
  Aptr=ACFunJac(Jac,&PgMax,FALSE,TRUE,FALSE);
  if(DCFunJac(Jac,FALSE,TRUE)) return(-iter);
  SVCFunJac(Jac,FALSE,TRUE);      /*  FACTS  */
  TCSCFunJac(Jac,FALSE,TRUE);     /*  FACTS  */
  STATCOMFunJac(Jac,FALSE,TRUE);  /*  FACTS  */
  if (PgMax<0) {
    if (Aptr!=NULL) {
      fCustomPrint(stderr,"\nError: Area %d %s does not have any spinning reserves.\n",Aptr->N,Aptr->Name);
      fCustomPrint(stderr,"       Increase the maximum P generation in this area, otherwise\n");
    } else {
      fCustomPrint(stderr,"\nError: The system does not have any spinning reserves.\n");
      fCustomPrint(stderr,"       Increase the maximum P generation in this system, otherwise\n");
    }
    fCustomPrint(stderr,"       the Jacobian matrix becomes singular.\n");
    WriteSolution(--iter,TrueParamStr(2),"Pg Max. Problems:");
    stopExecute(1);
  }
  if (!RightEvector) TransposeMatrix(Jac);
  SortRowsColumns(Jac);
  if(factorns(Jac,alpha,RowPartition,ColPartition,NewRow,NewCol,OldRow,OldCol)){
    fCustomPrint(stderr,"*** Singular Jacobian (possible voltage collapse, contol or limit problems).\n");
    fCustomPrint(stderr,"    Try changing the load levels, controls or limits, or use the -F option.\n");
    WriteSolution(--iter,TrueParamStr(2),"Singular Jacobian:");
    stopExecute(1);
  }
  SortRowsColumns(Jac);
  for(i=1;i<=N;i++) x0[i]=dF[i]=1;
  Nx0=1;
  val=0;
  count=0;
  while(fabs((Nx0-val)/Nx0)>tol) {
    val=Nx0;
    repsolp(Jac,x0,OldRow,NewCol);
    for(Nx0=0,i=1;i<=N;i++) {
      if(x0[i]==dF[i]) x0[i]=0;
      dF[i]=x0[i];
      if(fabs(x0[i])>Nx0) Nx0=fabs(x0[i]);
    }
    for(i=1;i<=N;i++) x0[i]=x0[i]/Nx0;
    count++;
    if(ExistParameter('d')) fCustomPrint(stderr,"Evect iter.=%d  Max_x0i=%lf  E_value=%lf\n",count,Nx0,1./Nx0);
    if(count>M) break;
  }
  *EigenValue=1./Nx0;
  if (ExistParameter('d')) {
    Out=OpenOutput("evect.dat");
    fCustomPrint(Out,"%d 1\n",N);
    for(i=1;i<=N;i++) fCustomPrint(Out,"%d 1 %-lg\n",i,x0[i]);
    fCustomPrint(Out,"0 0 0.0\n");
    fclose(Out);
  }
  return(iter);
}


/* --------------------------- PrintLeftEvector --------------------------------- */
#ifdef ANSIPROTO
void PrintLeftEvector(INDEX N,FILE *Out)
#else
void PrintLeftEvector(N,Out)
INDEX N;
FILE *Out;
#endif
/* Print "zero" left e-vector */
{
  INDEX i,l,k,I,J;
  ACbusData *ACptr;
  DCbusData *DCptrR;
  SVCbusData *SVCptr;                  /* FACTS */
  TCSCbusData *TCSCptr;                /* FACTS */
  STATCOMbusData *STATCOMptr;          /* FACTS */
  ElementData *Eptr;
  ElementList *ELptr;
  char str[80];

  fCustomPrint(Out,"%d 1\n",N);
  for (i=0,ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
    sprintf_s(str,"dP%-d",ACptr->Num); i++;
    fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
    sprintf_s(str,"dQ%-d",ACptr->Num);  i++;
    fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
    if(Acont && strpbrk(ACptr->Type,"A")){
      sprintf_s(str,"dPA%-d",ACptr->Area->N); i++;
      fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
    }
    if (PQcont) for(ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next) {
      Eptr=ELptr->Eptr;
      if(strpbrk(Eptr->Type,"PQNM")) {
         if (Eptr->From==ACptr) {
           I=Eptr->From->Num;
           J=Eptr->To->Num;
         } else {
           J=Eptr->From->Num;
           I=Eptr->To->Num;
         }
         if(!strcmp(Eptr->Type,"RP") || strpbrk(Eptr->Type,"PM")){
           sprintf_s(str,"dP%-d_%-d",I,J); i++;
           fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
         } else {
           sprintf_s(str,"dQ%-d_%-d",I,J); i++;
           fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
         }
      }
    }
    if (ACptr->Gen!=NULL) {
      i=ACptr->Gen->Nvar;
      sprintf_s(str,"dPg%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dQg%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dEq%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dEd%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dVd%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dVq%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dId%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dIq%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dVr%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dVi%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
      sprintf_s(str,"dIa%-d",ACptr->Num); fCustomPrint(Out,"%4d %8s %-11.5g\n",++i,str,x0[i]);
    }
  }
  for(k=0,DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next) if(!strcmp(DCptrR->Type,"R")){
    for (k++,l=1;l<=11;l++){
      sprintf_s(str,"Fdc%-d_%-d",k,l); i++;
      fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
    }
  }
                                       /* FACTS */
  for(k=0,SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
    for (k++,l=1;l<=3;l++){
      sprintf_s(str,"Fsvc%-d_%-d",k,l); i++;
      fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
    }
  }
  for(k=0,TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;TCSCptr=TCSCptr->Next){
    for (k++,l=1;l<=7;l++){
      sprintf_s(str,"Ftcsc%-d_%-d",k,l); i++;
      fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
    }
  }
  for(k=0,STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
    for (k++,l=1;l<=7;l++){
      sprintf_s(str,"Fstat%-d_%-d",k,l); i++;
      fCustomPrint(Out,"%4d %8s %-11.5g\n",i,str,x0[i]);
    }
  }
                                    /* END OF FACTS */

  fCustomPrint(Out,"%4d %8s %-11.5g\n",0,"0",0.);
  fclose(Out);
}

/* ------------------- PoCPoint --------------------- */
#ifdef ANSIPROTO
int PoCPoint(void)
#else
int PoCPoint()
#endif
/* Find PoC. */
{
  int count,iter,M;
  INDEX i,N,J;
  FILE *Out;
  BOOLEAN flag,flagp,flags,flagD;
  VALUETYPE NDx,Nlim,Nlimp,EigenValue;

  RealParameter('L',&lambda,-1e6,1e6);
  param0=lambda;
  if ((ExistParameter('D') && (!NullName(NameParameter('D')))) || flagVloads) flagD=TRUE;
  else flagD=FALSE;

  /* Solve base load  */
  if (lambda!=0 && !flagD) {
    InitializeLoad();
    iter=Pflow(1,FALSE,TRUE,FALSE);
    fCustomPrint(stderr,"Loading factor -> %-10.6lg  ",lambda);
  } else {
    iter=Pflow(1,FALSE,TRUE,TRUE);
    fCustomPrint(stderr,"Loading factor -> %-10.6lg  ",0.);
  }
  if(iter<0) return(iter);
  fCustomPrint(stderr,"**** Base Case Solved ****\n\n");
  flagp=DCsetup();

  /* Initial loading of system to get it closer to bifurcation */
  if(Direction(Jac,Dx,flagp)<0) return(-iter);
  for(NDx=0,i=1;i<=Jac->n1;i++) if(fabs(Dx[i])>NDx) NDx=fabs(Dx[i]);
  flagR=TRUE;
  if (NDx) {
    if (Nvolt!=0) Dparam=sqrt((double) Nvolt)/NDx;
    else Dparam=1/NDx;
    for(i=1;i<=Jac->n1;i++) Dx[i]=Dparam*Dx[i];
    Nlim=LoadX0(TRUE,TRUE,FALSE);
    if (ExistParameter('d')) fCustomPrint(stderr,"Dparam=%lf   Nlim=%lf\n",Dparam,Nlim);
    Nlimp=0;
    if (Qlim) Nlimp+=Nvolt;
    if (Zlim) Nlimp+=NZvolt;
    if (Tlim) Nlimp+=NregV;
    if (PQlim) Nlimp+=NregPQ;
    while((Nlim-1)>Nlimp*0.1) {
      Dparam=0.8*Dparam;
      for(i=1;i<=Jac->n1;i++) Dx[i]=0.8*Dx[i];
      Nlim=LoadX0(FALSE,TRUE,FALSE);
      if (ExistParameter('d')) fCustomPrint(stderr,"Dparam=%lf   Nlim=%lf\n",Dparam,Nlim);
    }
    iter=Pflow(1,flagp,TRUE,FALSE);
    if(iter<0) iter= -iter;
    else { param0=lambda; Dparam=0;}
    fCustomPrint(stderr,"****  Initial Loading  ****   ");
    fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
  }
  J=4;
  Evector(J,iter,0.001,FALSE,&EigenValue);

  /* Direct method */
  for(flagL=flags=TRUE,M=5,count=0;;) {
    N=Jac->n2=Jac->n1=2*(NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom)+1;    /* FACTS */
    NewRow->N=NewCol->N=N;
    OldRow->N=OldCol->N=N;
    RowPartition->p[1]=ColPartition->p[1]=N;
    iter=Pflow(iter,TRUE,TRUE,FALSE);
    if(ExistParameter('d')) fCustomPrint(stderr,"Loading factor -> %-10.6lg\n",lambda);
    if (count==M) break;
    /* If no convergence, recalculate e-vector */
    if (iter<0) {
      iter= -iter;
      if (iter>MaxIter) {
        fCustomPrint(stderr,"\n *** The PoC case has not been solved (possible initial guess problems, or\n");
        fCustomPrint(stderr,"     AC/DC/FACTS problems, or too few iterations). Try running the case using\n");
        fCustomPrint(stderr,"     the -F option, or change load levels, AC/DC controls, or increase the\n");
        fCustomPrint(stderr,"     maximum number of iterations with the -M option.\n");
        fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
        WriteSolution(iter,TrueParamStr(2),"Unsolved PoC case:");
        stopExecute(1);
      }
      flag=ChangeDCmode();
      if (!flag) flag=ChangeSVCmode(); else ChangeSVCmode();            /* FACTS */
      if (!flag) flag=ChangeTCSCmode(); else ChangeTCSCmode();          /* FACTS */
      if (!flag) flag=ChangeSTATCOMmode(); else ChangeSTATCOMmode();    /* FACTS */
      if (!flag) {
        flag=CheckRlimits();
        if(!flag) flag=CheckVlimits(); else CheckVlimits();
        if(!flag) flag=CheckQlimits(); else CheckQlimits(); }
      if (flag) count=0;
      if (count>0 && flags) {J=0; flags=FALSE;}
      else {J=4; flags=TRUE;}
      if((iter=Evector(J,iter,0.001,FALSE,&EigenValue))<0) return(-iter);
    } else { count=0; break;}
    if (!flag) count++;
  }

  /* -------------------------- Print left e-vector ----------------------- */
  if (!NullName(NameParameter('C'))) {
    Out=OpenOutput(NameParameter('C'));
    N=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;     /* FACTS */
    PrintLeftEvector(N,Out);
  }

  if (count==M) {
    fCustomPrint(stderr,"\n *** The PoC case has not been solved (possible initial guess problems, or\n");
    fCustomPrint(stderr,"     AC/DC/FACTS problems, or too few iterations). Try running the case using\n");
    fCustomPrint(stderr,"     the -F option, or change load levels, AC/DC controls, or increase the\n");
    fCustomPrint(stderr,"     maximum number of iterations with the -M option.\n");
    fCustomPrint(stderr,"Loading factor -> %-10.6lg\n\n",lambda);
    WriteSolution(iter,TrueParamStr(2),"Unsolved PoC case:");
    stopExecute(1);
  }
  return(iter);
}
