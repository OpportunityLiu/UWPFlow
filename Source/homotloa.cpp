/* Homotopy Continuation Method: Load variables */

#include "homot.h"

/* ------- Global Variables ------ */
extern VALUETYPE *Dx,Dparam,param0,*x0,*x0p,Kh,Htol,SD0,AngTr,
                 VoltageSum,VoltageSum0,DxiMax,VSF,SF,ZeroDx,Tol;
extern INDEX NewNumEq,CountSteps,NumSteps;
extern AClist *Vlist,*Vlistp;
extern BOOLEAN flagReducedContinuation,flagReduceSystem,*DxZero,flagBS;


/* ------------------ LoadX0 ----------------------------- */
#ifdef ANSIPROTO
VALUETYPE LoadX0(BOOLEAN FlagLoadX0,BOOLEAN FlagUpdateVar,BOOLEAN FlagMakeDxZero)
#else
VALUETYPE LoadX0(FlagLoadX0,FlagUpdateVar,FlagMakeDxZero)
BOOLEAN FlagLoadX0,FlagUpdateVar,FlagMakeDxZero;
#endif
/* Load x0 vector and update AC/DC variables for power flow. */
{
  ACbusData *ACptr,*ACptrp,*Ptr,*BEptr;
  AClist *ALptr;
  DCbusData *DCptrR,*DCptrI,*DCptr;
  SVCbusData *SVCptr;                    /* FACTS */
  TCSCbusData *TCSCptr;                  /* FACTS */
  STATCOMbusData *STATCOMptr;            /* FACTS */
  ElementData *Eptr;
  ElementList *ELptr;
  VALUETYPE DPg,val,valp,vals,count,consp=0,Pg,Pmax,PgMax,Q,Qm;
  INDEX i,j;
  BOOLEAN Recover=TRUE;
  char Qmax[5],Qmin[5];

  if (ExistParameter('G')) Recover=FALSE;
  NewNumEq=Jac->n1 - 1 ;
  if (FlagLoadX0) VoltageSum0=Tol;
  if (FlagUpdateVar) VoltageSum=0;
  if (!Bl) {
    if (FlagLoadX0) param0=lambda;
    if (FlagUpdateVar) {
      lambda=param0+Dparam;
      if (lambda<=0) {
        val=fabs(lambda/Dparam);
        if(ExistParameter('d')) fprintf(stderr,"lambdaMin  %lf\n",val);
        lambda=0;}
      else val=0;
      if (val>consp) 
		  consp=val;
    }
  }
  for(ALptr=dataPtr->KGbus;ALptr!=NULL;ALptr=ALptr->Next) {
    ACptr=ALptr->AC;
    if (strpbrk(ACptr->Type,"S")) {
      i=ACvar[ACptr->N];
      if (FlagLoadX0)  {
        x0[i]=ACptr->Kg;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Kg=x0[i]+Dx[i];
      BEptr=ACptr;
    }
    else if(Acont){
      i=ACvar[ACptr->N]+2;
      if (FlagLoadX0) {
        x0[i]=ACptr->Kg;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Kg=x0[i]+Dx[i];
    }
  }
  for (Ptr=NULL,valp=0,count=0,ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
    if (FlagUpdateVar) {
      if (ACptr->Area!=NULL) BEptr=ACptr->Area->Slack;
      DPg=ACptr->DPg;
      if (flagBS) {
        if (ACptr!=BEptr) Pg=ACptr->Pg+lambda*DPg;
        else              Pg=ACptr->Pg+BEptr->Kg;
      }
      else  Pg=ACptr->Pg+BEptr->Kg*DPg;
      Pmax=ACptr->Smax*ACptr->Smax-ACptr->Qg*ACptr->Qg;
      if (Pmax>0) Pmax=sqrt(Pmax);
      else        Pmax=99999999.;
      if (!flagSmax && Pmax<ACptr->PgMax) PgMax=Pmax;
      else if (!flagPgMax)                PgMax=ACptr->PgMax;
      else                                PgMax=99999999.;
      if (Pg>PgMax) {
        Pg=PgMax;
        DPg=0;
      }
      ACptr->PG=Pg;
      ACptr->DPG=DPg;
      ACptr->Pmax=PgMax;
      if (ACptr->Qmax==ACptr->Max) strcpy_s(Qmax,"Qmax");
      else                         strcpy_s(Qmax,"Smax");
      if (ACptr->Qmin==ACptr->Min) strcpy_s(Qmin,"Qmin");
      else                         strcpy_s(Qmin,"Smin");
    }
    if (!strcmp(ACptr->Type,"B")||!strcmp(ACptr->Type,"BA")) {
      i=ACvar[ACptr->N];
      if (FlagLoadX0) {
        x0[i]=ACptr->Ang;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->V;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        if (InList(ACptr,Vlistp)) VoltageSum0+=ACptr->V;
      }
      if (FlagUpdateVar) {
        ACptr->V=x0[i]+Dx[i];
        if (ACptr->V<=0){
          count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=0.00001;}
        else val=0;
        if (val>consp) 
			consp=val;
        if (InList(ACptr,Vlistp)) VoltageSum+=ACptr->V;
      }
      if(!Bl && FlagLoadX0 && FlagUpdateVar &&  FlagMakeDxZero && fabs(Dx[i])>valp) {
        valp=fabs(Dx[i]*param0/x0[i]);
        Ptr=ACptr;
      }
    }
    else if(strpbrk(ACptr->Type,"L")){
      i=ACvar[ACptr->N];
      if (FlagLoadX0) {
        x0[i]=ACptr->Ang;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      i++; if (FlagLoadX0) x0[i]=lambda;
      if (FlagUpdateVar) {
        lambda=x0[i]+Dx[i];
        if (lambda<=0) {
          val=fabs(lambda/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"lambdaMin  %lf\n",val);
          lambda=0;}
        else val=0;
        if (val>consp) 
			consp=val;
      }
      if (FlagLoadX0) {
        param0=ACptr->V;
        if (InList(ACptr,Vlistp)) VoltageSum0+=ACptr->V;
      }
      if (FlagUpdateVar) {
        ACptr->V=param0+Dparam;
        if (ACptr->V<=0){
          count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=0.00001;}
        else val=0;
        if (val>consp) 
			consp=val;
        if (InList(ACptr,Vlistp)) VoltageSum+=ACptr->V;
        if(fabs(Dparam)>param0) {
          val=fabs(1-param0/fabs(Dparam)/8.0);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s lambdaMin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);}
        else val=0;
        if (val>consp) 
			consp=val;
      }
      if(FlagLoadX0 && FlagUpdateVar &&  FlagMakeDxZero && fabs(Dx[i])>valp) {
        valp=fabs(Dx[i]*param0/x0[i]);
        Ptr=ACptr;
      }
    }
    else if(strpbrk(ACptr->Type,"C")){
      i=ACvar[ACptr->N];
      if (FlagLoadX0) {
        x0[i]=ACptr->Ang;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      i++;
      if (!QRcont) {
        if (FlagLoadX0) {
          x0[i]=ACptr->V;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
          if (InList(ACptr,Vlistp)) VoltageSum0+=ACptr->V;
        }
        if (FlagUpdateVar) {
          ACptr->V=x0[i]+Dx[i];
          if (ACptr->V<=0){
            count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->V=0.00001;}
          else val=0;
          if (val>consp) 
			  consp=val;
          if (InList(ACptr,Vlistp)) VoltageSum+=ACptr->V;
        }
        if(FlagLoadX0 && FlagUpdateVar &&  FlagMakeDxZero && fabs(Dx[i])>valp) {
          valp=fabs(Dx[i]*param0/x0[i]);
          Ptr=ACptr;
        }
      } else if (ACptr->Kbg<1) {
        ACptrp=ACptr->ContBus->AC;
        if (FlagLoadX0) {
          x0[i]=ACptr->V;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->V=x0[i]+Dx[i];
          if (ACptr->V<=0){
            count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->V=0.00001;}
          else if (strpbrk(ACptrp->cont,"Q")) {
            if (Recover && ACptrp->Qg>=ACptrp->Max && ACptr->V>=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s Recover%s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmax,val);
              ACptr->V=ACptr->VCont;}
            else if (Recover && ACptrp->Qg<=ACptrp->Min && ACptr->V<=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s Recover%s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmin,val);
              ACptr->V=ACptr->VCont;}
            else val=0;
          }
          else if (strpbrk(ACptrp->cont,"E")) {
            if (Recover && ACptrp->Gen->Eq>=ACptrp->Gen->EqMax && ACptr->V>=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverEqMax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->V=ACptr->VCont;}
            else if (Recover && ACptrp->Gen->Eq<=ACptrp->Gen->EqMin && ACptr->V<=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverEqMin%lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->V=ACptr->VCont;}
            else val=0;
          }
          else if (strpbrk(ACptrp->cont,"I")) {
            if (Recover && ACptrp->Gen->Ia>=ACptrp->Gen->IaMax && ACptr->V>=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverIaMax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->V=ACptr->VCont;}
            else val=0;
          }
          else val=0;
          if (val>consp) 
			  consp=val;
        }
      } else {
        if (FlagLoadX0) {
          x0[i]=ACptr->Qr;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) ACptr->Qr=x0[i]+Dx[i];
      }
    }
    else if(strpbrk(ACptr->Type,"T")){
      i=ACvar[ACptr->N];
      if (FlagLoadX0) {
        x0[i]=ACptr->Ang;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      i++;
      if (Rcont) {
        for(ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next){
          Eptr=ELptr->Eptr;
          if(!strcmp(Eptr->Type,"R")) {
            if (FlagLoadX0) {
              x0[i]=Eptr->Tap;
              if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
            }
            if (FlagUpdateVar) {
              Eptr->Tap=x0[i]+Dx[i];
              if (Tlim && Eptr->Tap<=1/Eptr->Tmax+0.00001){
                count++; val=fabs((Eptr->Tap-1/Eptr->Tmax)/Dx[i]);
                if(ExistParameter('d')) fprintf(stderr,"%s %d %s %d %s Tmax %lf\n",Eptr->Type,Eptr->From->Num,Eptr->From->Name,
                                                Eptr->To->Num,Eptr->To->Name,val);
                Eptr->Tap=1/Eptr->Tmax;}
              else if (Tlim && Eptr->Tap>=1/Eptr->Tmin-0.00001){
                count++; val=fabs((Eptr->Tap-1/Eptr->Tmin)/Dx[i]);
                if(ExistParameter('d')) fprintf(stderr,"%s %d %s %d %s Tmin %lf\n",Eptr->Type,Eptr->From->Num,Eptr->From->Name,
                                                Eptr->To->Num,Eptr->To->Name,val);
                Eptr->Tap=1/Eptr->Tmin;}
              else val=0;
              if (val>consp) 
				  consp=val;
            }
          } else {
            if (FlagLoadX0) {
              x0[i]=ACptr->V;
              if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
            }
            if (FlagUpdateVar) {
              ACptr->V=x0[i]+Dx[i];
              if (ACptr->V<=0){
                count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
                if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
                ACptr->V=0.00001;}
              else val=0;
              if (val>consp) 
				  consp=val;
            }
          }
        }
      }
    }
    else if(strpbrk(ACptr->Type,"R")){
      i=ACvar[ACptr->N];
      if (FlagLoadX0) {
        x0[i]=ACptr->Ang;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->V;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) {
        ACptr->V=x0[i]+Dx[i];
        if (Rcont && ACptr->V>=ACptr->Vmax){
          count++; val=fabs((ACptr->V-ACptr->Vmax)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vmax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=ACptr->Vmax;}
        else if (Rcont && ACptr->V<=ACptr->Vmin){
          count++; val=fabs((ACptr->V-ACptr->Vmin)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vmin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=ACptr->Vmin;}
        else if (ACptr->V<=0){
          count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=0.00001;}
        else val=0;
        if (val>consp) 
			consp=val;
      }
    }
    else if(strpbrk(ACptr->Type,"V")){
      i=ACvar[ACptr->N];
      if (FlagLoadX0) {
        x0[i]=ACptr->Ang;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      i++;
      if (!strpbrk(ACptr->cont,"V")) {
        if (FlagLoadX0) {
          x0[i]=ACptr->V;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->V=x0[i]+Dx[i];
          if (ACptr->V<=0){
            count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->V=0.00001;}
          else if (Vlim && ACptr->V>=ACptr->Vmax){
            count++; val=fabs((ACptr->V-ACptr->Vmax)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vmax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->V=ACptr->Vmax;}
          else if (Vlim && ACptr->V<=ACptr->Vmin){
            count++; val=fabs((ACptr->V-ACptr->Vmin)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vmin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->V=ACptr->Vmin;}
          else val=0;
          if (val>consp) 
			  consp=val;
        }
      }
      else {
        if (FlagLoadX0) {
          x0[i]=ACptr->Qg;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->Qg=x0[i]+Dx[i];
          if (Recover && ACptr->V>=ACptr->Vmax && ACptr->Qg>=ACptr->VCont){
            count++; val=fabs((ACptr->Qg-ACptr->VCont)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverVmax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->Qg=ACptr->VCont;}
          else if (Recover && ACptr->V<=ACptr->Vmin && ACptr->Qg<=ACptr->VCont){
            count++; val=fabs((ACptr->Qg-ACptr->VCont)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverVmin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->Qg=ACptr->VCont;}
          else val=0;
          if (val>consp) 
			  consp=val;
        }
      }
    }
    else if(strpbrk(ACptr->Type,"X")){
      i=ACvar[ACptr->N];
      if (FlagLoadX0) {
        x0[i]=ACptr->Ang;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->V;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) {
        ACptr->V=x0[i]+Dx[i];
        if (ACptr->V<=0){
          count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=0.00001;}
        else if (Xlim && strpbrk(ACptr->cont,"X") && ACptr->V>=ACptr->Vmax){
          count++; val=fabs((ACptr->V-ACptr->Vmax)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vmax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=ACptr->Vmax;}
        else if (Xlim && strpbrk(ACptr->cont,"X") && ACptr->V<=ACptr->Vmin){
          count++; val=fabs((ACptr->V-ACptr->Vmin)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vmin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=ACptr->Vmin;}
        else if (Recover && strpbrk(ACptr->cont,"M") && ACptr->V<=ACptr->Vmax){
          count++; val=fabs((ACptr->V-ACptr->Vmax)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vmax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=ACptr->Vmax;}
        else if (Recover && strpbrk(ACptr->cont,"m") && ACptr->V>=ACptr->Vmin){
          count++; val=fabs((ACptr->V-ACptr->Vmin)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vmin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=ACptr->Vmin;}
        else val=0;
        if (val>consp) 
			consp=val;
      }
    }
    else if(QRcont && strpbrk(ACptr->Type,"G")){
      i=ACvar[ACptr->N];
      if (!strpbrk(ACptr->Type,"S")) {
        if (FlagLoadX0) {
          x0[i]=ACptr->Ang;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      }
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->V;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) {
        ACptr->V=x0[i]+Dx[i];
        if (ACptr->V<=0){
          count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
          ACptr->V=0.00001;}
        else val=0;
        if (val>consp) 
			consp=val;
      }
    }
    else if(strpbrk(ACptr->Type,"Q,S")||(!QRcont && strpbrk(ACptr->Type,"G"))){
      i=ACvar[ACptr->N];
      if (!strpbrk(ACptr->Type,"S")) {
        if (FlagLoadX0) {
          x0[i]=ACptr->Ang;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      }
      i++;
      if (strpbrk(ACptr->cont,"V")) {
        if (FlagLoadX0) {
          x0[i]=ACptr->Qg;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->Qg=x0[i]+Dx[i];
          if (Qlim && ACptr->Qg>=ACptr->Max){///
            count++; val=fabs((ACptr->Qg-ACptr->Max)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmax,val);
            ACptr->Qg=ACptr->Max;}
          else if (Qlim && ACptr->Qg<=ACptr->Min){
            count++; val=fabs((ACptr->Qg-ACptr->Min)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmin,val);
            ACptr->Qg=ACptr->Min;}
          else val=0;
          if (val>consp) //
			  consp=val;
        }
      } else {
        if (FlagLoadX0) {
          x0[i]=ACptr->V;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->V=x0[i]+Dx[i];
          if (ACptr->V<=0){
            count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->V=0.00001;}
          else if (strpbrk(ACptr->cont,"Q")) {
            if (Recover && ACptr->Qg>=ACptr->Max && ACptr->V>=ACptr->VCont){///
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s Recover%s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmax,val);
              ACptr->V=ACptr->VCont;}
            else if (Recover && ACptr->Qg<=ACptr->Min && ACptr->V<=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s Recover%s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmin,val);
              ACptr->V=ACptr->VCont;}
            else val=0;
          }
          else if (strpbrk(ACptr->cont,"E")) {
            if (Recover && ACptr->Gen->Eq>=ACptr->Gen->EqMax && ACptr->V>=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverEqMax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->V=ACptr->VCont;}
            else if (Recover && ACptr->Gen->Eq<=ACptr->Gen->EqMin && ACptr->V<=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverEqMin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->V=ACptr->VCont;}
            else val=0;
          }
          else if (strpbrk(ACptr->cont,"I")) {
            if (Recover && ACptr->Gen->Ia>=ACptr->Gen->IaMax && ACptr->V>=ACptr->VCont){
              count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverIaMax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->V=ACptr->VCont;}
            else val=0;
          }
          else val=0;//
          if (val>consp) 
			  consp=val;
        }
      }
    }
    else if(strpbrk(ACptr->Type,"Z")){
      i=ACvar[ACptr->N];
      if (FlagLoadX0) {
        x0[i]=ACptr->Ang;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Ang=x0[i]+Dx[i];
      i++;
      if (strpbrk(ACptr->cont,"V")) {
        if (FlagLoadX0) {
          x0[i]=ACptr->Qg;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->Qg=x0[i]+Dx[i];
          if (Zlim && ACptr->Qg>=ACptr->Qmax){
            count++; val=fabs((ACptr->Qg-ACptr->Qmax)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmax,val);
            ACptr->Qg=ACptr->Qmax;}
          else if (Zlim && ACptr->Qg<=ACptr->Qmin){
            count++; val=fabs((ACptr->Qg-ACptr->Qmin)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmin,val);
            ACptr->Qg=ACptr->Qmin;}
          else val=0;
          if (val>consp) 
			  consp=val;
        }
      } else {
        if (FlagLoadX0) {
          x0[i]=ACptr->V;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->V=x0[i]+Dx[i];
          ACptr->Qg=ACptr->V*ACptr->V*ACptr->Bz;
          if (ACptr->V<=0){
            count++; val=fabs((ACptr->V-0.00001)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Vzero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->V=0.00001;}
          else if (Recover && ACptr->Qg>=ACptr->Qmax && ACptr->V>=ACptr->VCont){
            count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Recover%s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmax,val);
            ACptr->V=ACptr->VCont;}
          else if (Recover && ACptr->Qg<=ACptr->Qmin && ACptr->V<=ACptr->VCont){
            count++; val=fabs((ACptr->V-ACptr->VCont)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s Recover%s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmin,val);
            ACptr->V=ACptr->VCont;}
          else val=0;
          if (val>consp) 
			  consp=val;
          ACptr->Qg=ACptr->V*ACptr->V*ACptr->Bz;
        }
      }
    }
    if (PQcont) for(ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next) {
      Eptr=ELptr->Eptr;
      if(strpbrk(Eptr->Type,"PQNM")) {
        i=ACvar[ACptr->N]+1+ACptr->Ncont-Eptr->Ncont;
        if (Acont && strpbrk(ACptr->Type,"A")) i++;
        if(!strcmp(Eptr->Type,"RP")) {
          if (FlagLoadX0) {
            x0[i]=Eptr->Ang;
            if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
          }
          if (FlagUpdateVar) {
            Eptr->Ang=x0[i]+Dx[i];
            if (Tlim && Eptr->Ang>=Eptr->Tmax){
              count++; val=fabs((Eptr->Ang-Eptr->Tmax)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %d %s Max %lf\n",Eptr->Type,Eptr->From->Num,Eptr->From->Name,
                                              Eptr->To->Num,Eptr->To->Name,val);
              Eptr->Ang=Eptr->Tmax;}
            else if (Tlim && Eptr->Ang<=Eptr->Tmin){
              count++; val=fabs((Eptr->Ang-Eptr->Tmin)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %d %s Min %lf\n",Eptr->Type,Eptr->From->Num,Eptr->From->Name,
                                              Eptr->To->Num,Eptr->To->Name,val);
              Eptr->Ang=Eptr->Tmin;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
        }
        else if(!strcmp(Eptr->Type,"RQ")){
          if (FlagLoadX0) {
            x0[i]=Eptr->Tap;
            if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
          }
          if (FlagUpdateVar) {
            Eptr->Tap=x0[i]+Dx[i];
            if (Tlim && Eptr->Tap<=1/Eptr->Tmax+0.00001){
              count++; val=fabs((Eptr->Tap-1/Eptr->Tmax)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %d %s Max %lf\n",Eptr->Type,Eptr->From->Num,Eptr->From->Name,
                                              Eptr->To->Num,Eptr->To->Name,val);
              Eptr->Tap=1/Eptr->Tmax;}
            else if (Tlim && Eptr->Tap>=1/Eptr->Tmin-0.00001){
              count++; val=fabs((Eptr->Tap-1/Eptr->Tmin)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %d %s Min %lf\n",Eptr->Type,Eptr->From->Num,Eptr->From->Name,
                                              Eptr->To->Num,Eptr->To->Name,val);
              Eptr->Tap=1/Eptr->Tmin;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
        }
        else if(strpbrk(Eptr->Type,"MN")) {
          if (FlagLoadX0) {
            x0[i]=Eptr->Cvar;
            if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
          }
          if (FlagUpdateVar) {
            Eptr->Cvar=x0[i]+Dx[i];
            if (PQlim && Eptr->Cvar>=Eptr->Max){
              count++; val=fabs((Eptr->Cvar-Eptr->Max)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %d %s Max %lf\n",Eptr->Type,Eptr->From->Num,Eptr->From->Name,
                                              Eptr->To->Num,Eptr->To->Name,val);
              Eptr->Cvar=Eptr->Max;}
            else if (PQlim && Eptr->Cvar<=Eptr->Min){
              count++; val=fabs((Eptr->Cvar-Eptr->Min)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %d %s Min %lf\n",Eptr->Type,Eptr->From->Num,Eptr->From->Name,
                                              Eptr->To->Num,Eptr->To->Name,val);
              Eptr->Cvar=Eptr->Min;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
        }
        else {
          if (FlagLoadX0) {
            x0[i]=Eptr->Cvar;
            if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
          }
          if (FlagUpdateVar) Eptr->Cvar=x0[i]+Dx[i];
        }
      }
    }
    if (ACptr->Gen!=NULL) {
      i=ACptr->Gen->Nvar+1;
      if (!strpbrk(ACptr->cont,"E")) {
        if (FlagLoadX0) {
          x0[i]=ACptr->Gen->Eq;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->Gen->Eq=x0[i]+Dx[i];
          if (Elim && ACptr->Gen->Eq>=ACptr->Gen->EqMax) {
            count++; val=fabs((ACptr->Gen->Eq-ACptr->Gen->EqMax)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s EqMax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->Gen->Eq=ACptr->Gen->EqMax;}
          else if (Elim && ACptr->Gen->Eq<=ACptr->Gen->EqMin) {
            count++; val=fabs((ACptr->Gen->Eq-ACptr->Gen->EqMin)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s EqMin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->Gen->Eq=ACptr->Gen->EqMin;}
          else if (ACptr->Gen->Eq<=0) {
            count++; val=fabs((ACptr->Gen->Eq-0.00001)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s EqZero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->Gen->Eq=0.00001;}
          else val=0;
          if (val>consp) 
			  consp=val;
        }
      }
      else {
        if (FlagLoadX0) {
          x0[i]=ACptr->Qg;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->Qg=x0[i]+Dx[i];
          if (strpbrk(ACptr->Type,"V")) {
            if (Recover && ACptr->Gen->Eq>=ACptr->Gen->EqMax && ACptr->Qg>=ACptr->VCont){
              count++; val=fabs((ACptr->Qg-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverEqMax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->Qg=ACptr->VCont;}
            else if (Recover && ACptr->Gen->Eq<=ACptr->Gen->EqMin && ACptr->Qg<=ACptr->VCont){
              count++; val=fabs((ACptr->Qg-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverEqMin %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->Qg=ACptr->VCont;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
          else {
            if (Qlim && ACptr->Qg>=ACptr->Max){
              count++; val=fabs((ACptr->Qg-ACptr->Max)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmax,val);
              ACptr->Qg=ACptr->Max;}
            else if (Qlim && ACptr->Qg<=ACptr->Min){
              count++; val=fabs((ACptr->Qg-ACptr->Min)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmin,val);
              ACptr->Qg=ACptr->Min;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
        }
      }
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->dg;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->dg=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->Vr;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->Vr=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->Vi;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->Vi=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->Ir;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->Ir=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->Ii;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->Ii=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->Vq;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->Vq=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->Vd;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->Vd=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->Iq;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->Iq=x0[i]+Dx[i];
      i++;
      if (FlagLoadX0) {
        x0[i]=ACptr->Gen->Id;
        if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
      }
      if (FlagUpdateVar) ACptr->Gen->Id=x0[i]+Dx[i];
      i++;
      if (!strpbrk(ACptr->cont,"I")) {
        if (FlagLoadX0) {
          x0[i]=ACptr->Gen->Ia;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->Gen->Ia=x0[i]+Dx[i];
          if (Ilim && ACptr->Gen->Ia>=ACptr->Gen->IaMax) {
            count++; val=fabs((ACptr->Gen->Ia-ACptr->Gen->IaMax)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s IaMax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->Gen->Ia=ACptr->Gen->IaMax;}
          else if (ACptr->Gen->Ia<=0) {
            count++; val=fabs((ACptr->Gen->Ia-0.00001)/Dx[i]);
            if(ExistParameter('d')) fprintf(stderr,"%s %d %s IaZero %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
            ACptr->Gen->Ia=0.00001;}
          else val=0;
          if (val>consp) 
			  consp=val;
        }
      }
      else {
        if (FlagLoadX0) {
          x0[i]=ACptr->Qg;
          if (FlagMakeDxZero && (DxZero!=NULL && (DxZero[i] || (flagReduceSystem && x0[i]!=0 && CountSteps>=NumSteps && fabs(Dx[i]/x0[i])<ZeroDx)))) { Dx[i]=0;  DxZero[i]=TRUE; NewNumEq--; }
        }
        if (FlagUpdateVar) {
          ACptr->Qg=x0[i]+Dx[i];
          if (strpbrk(ACptr->Type,"V")) {
            if (Recover && ACptr->Gen->Ia>=ACptr->Gen->IaMax && ACptr->Qg>=ACptr->VCont){
              count++; val=fabs((ACptr->Qg-ACptr->VCont)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s RecoverIaMax %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,val);
              ACptr->Qg=ACptr->VCont;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
          else {
            if (Qlim && ACptr->Qg>=ACptr->Max){
              count++; val=fabs((ACptr->Qg-ACptr->Max)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmax,val);
              ACptr->Qg=ACptr->Max;}
            else if (Qlim && ACptr->Qg<=ACptr->Min){
              count++; val=fabs((ACptr->Qg-ACptr->Min)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmin,val);
              ACptr->Qg=ACptr->Min;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
        }
      }
    }
  }

  if (FlagUpdateVar)  for (ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
    if(QRcont && strpbrk(ACptr->Type,"G") && strpbrk(ACptr->cont,"V")){
      i=ACvar[ACptr->Cont->N]+1;
      ACptr->Qg=(x0[i]+Dx[i])*ACptr->Kbg;
      if (Qlim && ACptr->Qg>=ACptr->Max){
        count++; val=fabs((ACptr->Qg-ACptr->Max)/(Dx[i]*ACptr->Kbg));
        if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmax,val);
        ACptr->Qg=ACptr->Max;}
      else if (Qlim && ACptr->Qg<=ACptr->Min){
        count++; val=fabs((ACptr->Qg-ACptr->Min)/(Dx[i]*ACptr->Kbg));
        if(ExistParameter('d')) fprintf(stderr,"%s %d %s %s %lf\n",ACptr->Type,ACptr->Num,ACptr->Name,Qmin,val);
        ACptr->Qg=ACptr->Min;}
      else val=0;
      if (val>consp) 
		  consp=val;
    }
  }

  i=NacVar;
  for(DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if(!strcmp(DCptrR->Type,"R")){
      for (j=1;j<=2;j++) {
        if (j==1) DCptr=DCptrR;
        else DCptr=DCptrI;
        if(strcmp(DCptr->Cont1,"VD")&&strcmp(DCptr->Cont2,"VD")) {
          i++; if (FlagLoadX0) x0[i]=DCptr->Vd;
          if (FlagUpdateVar) {
            DCptr->Vd=x0[i]+Dx[i];
            if (j==2 && ((DCptr->Tap>=DCptr->TapMax && DCptr->Vd>=DCptr->VdN)||
                (DCptr->Tap<=DCptr->TapMin && DCptr->Vd<=DCptr->VdN) ||
                (strpbrk(DCptr->Cont2,"IP")&&DCptr->Tap>DCptr->TapMin&&DCptr->Vd>=DCptr->VdN))){
              val=fabs((DCptr->Vd-DCptr->VdN)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %s VdN %lf\n",DCptr->Type,DCptr->Name,val);
              DCptr->Vd=DCptr->VdN;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
        }
        if(strcmp(DCptr->Cont1,"AT")&&strcmp(DCptr->Cont2,"AT")) {
          i++; if (FlagLoadX0) x0[i]=DCptr->Tap*DCptr->Ntrf;
          if (FlagUpdateVar) {
            DCptr->Tap=(x0[i]+Dx[i])/DCptr->Ntrf;
            if (DCptr->Tap>=DCptr->TapMax){
              val=fabs(DCptr->Ntrf*(DCptr->Tap-DCptr->TapMax)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %s Tmax %lf\n",DCptr->Type,DCptr->Name,val);
              DCptr->Tap=DCptr->TapMax;}
            else if (DCptr->Tap<=DCptr->TapMin){
              val=fabs(DCptr->Ntrf*(DCptr->Tap-DCptr->TapMin)/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %s Tmin %lf\n",DCptr->Type,DCptr->Name,val);
              DCptr->Tap=DCptr->TapMin;}
            else val=0;
            if (val>consp) 
				consp=val;
          }
        }
        if(strcmp(DCptr->Cont1,"AL")&&strcmp(DCptr->Cont2,"AL")) {
          i++; if (FlagLoadX0) x0[i]=cos(DCptr->Alfa);
          if (FlagUpdateVar) {
            DCptr->val[0]=x0[i]+Dx[i];
            if (DCptr->val[0]<=cos(DCptr->AlfaMax)){
              val=fabs((DCptr->val[0]-cos(DCptr->AlfaMax))/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %s Amax %lf\n",DCptr->Type,DCptr->Name,val);
              DCptr->Alfa=DCptr->AlfaMax;}
            if (DCptr->val[0]>=cos(DCptr->AlfaMin)){
              val=fabs((DCptr->val[0]-cos(DCptr->AlfaMin))/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %s Amin %lf\n",DCptr->Type,DCptr->Name,val);
              DCptr->Alfa=DCptr->AlfaMin;}
            else if (j==1 && ((DCptr->Tap>=DCptr->TapMax && DCptr->val[0]<=cos(DCptr->AlfaN))||
                     (DCptr->Tap<=DCptr->TapMin && DCptr->val[0]>=cos(DCptr->AlfaN)))){
              val=fabs((DCptr->val[0]-cos(DCptr->AlfaN))/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %s Anom %lf\n",DCptr->Type,DCptr->Name,val);
              DCptr->Alfa=DCptr->AlfaN;}
            else { val=0; DCptr->Alfa=acos(DCptr->val[0]);}
            if (val>consp) 
				consp=val;
          }
        }
        if(strcmp(DCptr->Cont1,"GA")&&strcmp(DCptr->Cont2,"GA")) {
          i++; if (FlagLoadX0) x0[i]=cos(DCptr->Gamma);
          if (FlagUpdateVar) {
            DCptr->val[0]=x0[i]+Dx[i];
            if (DCptr->val[0]>=cos(DCptr->GammaMin)){
              val=fabs((DCptr->val[0]-cos(DCptr->GammaMin))/Dx[i]);
              if(ExistParameter('d')) fprintf(stderr,"%s %s Gmin %lf\n",DCptr->Type,DCptr->Name,val);
              DCptr->Gamma=DCptr->GammaMin;}
            else { val=0; DCptr->Gamma=acos(DCptr->val[0]);}
            if (val>consp) 
				consp=val;
          }
        }
        i++; if (FlagLoadX0) x0[i]=DCptr->MVA;
        if (FlagUpdateVar) DCptr->MVA=x0[i]+Dx[i];
        if(strcmp(DCptr->Cont1,"PA")&&strcmp(DCptr->Cont2,"PA")) {
          i++; if (FlagLoadX0) x0[i]=DCptr->P;
          if (FlagUpdateVar) DCptr->P=x0[i]+Dx[i];
        }
        if(strcmp(DCptr->Cont1,"QA")&&strcmp(DCptr->Cont2,"QA")) {
          i++; if (FlagLoadX0) x0[i]=DCptr->Q;
          if (FlagUpdateVar) DCptr->Q=x0[i]+Dx[i];
        }
      }
      if(strcmp(DCptrR->Cont1,"ID")&&strcmp(DCptrR->Cont2,"ID")&&
         strcmp(DCptrI->Cont1,"ID")&&strcmp(DCptrI->Cont2,"ID")) {
        i++; if (FlagLoadX0) x0[i]=DCptrR->Id;
        if (FlagUpdateVar) DCptrR->Id=DCptrI->Id=x0[i]+Dx[i];
      }
    }
  }
                            /* FACTS */
  i=NacVar+11*Ndc/2;
  for(SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
   i++; if(FlagLoadX0) x0[i]=SVCptr->Qsvc;
   if(FlagUpdateVar) SVCptr->Qsvc=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=SVCptr->Bv;
   if(FlagUpdateVar) SVCptr->Bv=x0[i]+Dx[i];
   if(!strcmp(SVCptr->Cont,"AL")){
     i++; if(FlagLoadX0) x0[i]=SVCptr->alpha_svc;
     if(FlagUpdateVar){
        SVCptr->val=x0[i]+Dx[i];
        if(SVCptr->val>=SVCptr->AlphaMax){
          val=fabs((SVCptr->val-SVCptr->AlphaMax)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s AlphaMax %lf\n",SVCptr->Name,val);
          SVCptr->alpha_svc=SVCptr->AlphaMax;}
        else if(SVCptr->val<=SVCptr->AlphaMin){
          val=fabs((SVCptr->val-SVCptr->AlphaMin)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s AlphaMin %lf\n",SVCptr->Name,val);
          SVCptr->alpha_svc=SVCptr->AlphaMin;}
        else {val=0;SVCptr->alpha_svc=SVCptr->val;}
        if (val>consp) 
			consp=val;
     }
   }
   else if(!strcmp(SVCptr->Cont,"MN")){
     i++; if(FlagLoadX0) x0[i]=SVCptr->Vvar;
     if(FlagUpdateVar){
        SVCptr->val=x0[i]+Dx[i];
        if(SVCptr->val<=SVCptr->Vref){
          val=fabs((SVCptr->val-SVCptr->Vref)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s Vref_min %lf\n",SVCptr->Name,val);
          SVCptr->Vvar=SVCptr->Vref;}
        else {val=0;SVCptr->Vvar=SVCptr->val;}
        if (val>consp) 
			consp=val;
     }
   }
   else {
     i++; if(FlagLoadX0) x0[i]=SVCptr->Vvar;
     if(FlagUpdateVar){
        SVCptr->val=x0[i]+Dx[i];
        if(SVCptr->val>=SVCptr->Vref){
          val=fabs((SVCptr->val-SVCptr->Vref)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s Vref_max %lf\n",SVCptr->Name,val);
          SVCptr->Vvar=SVCptr->Vref;}
        else {val=0;SVCptr->Vvar=SVCptr->val;}
        if (val>consp) 
			consp=val;
     }
   }
  }

  i=NacVar+11*Ndc/2+3*Nsvc;
  for(TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;TCSCptr=TCSCptr->Next){
   i++; if(FlagLoadX0) x0[i]=TCSCptr->Ptcsc;
   if(FlagUpdateVar) TCSCptr->Ptcsc=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=TCSCptr->Qtcsck;
   if(FlagUpdateVar) TCSCptr->Qtcsck=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=TCSCptr->Qtcscm;
   if(FlagUpdateVar) TCSCptr->Qtcscm=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=TCSCptr->Be;
   if(FlagUpdateVar) TCSCptr->Be=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=TCSCptr->alpha_tcsc;
   if(FlagUpdateVar){
     TCSCptr->val=x0[i]+Dx[i];
     if(TCSCptr->val>=TCSCptr->AlphaMax){
       val=fabs((TCSCptr->val-TCSCptr->AlphaMax)/Dx[i]);
       if(ExistParameter('d')) fprintf(stderr,"%s AlphaMax %lf\n",TCSCptr->Name,val);
       TCSCptr->alpha_tcsc=TCSCptr->AlphaMax;}
     else if(TCSCptr->val<=TCSCptr->AlphaMin){
       val=fabs((TCSCptr->val-TCSCptr->AlphaMin)/Dx[i]);
       if(ExistParameter('d')) fprintf(stderr,"%s AlphaMin %lf\n",TCSCptr->Name,val);
       TCSCptr->alpha_tcsc=TCSCptr->AlphaMin;}
     else {val=0; TCSCptr->alpha_tcsc=TCSCptr->val;}
     if (val>consp) 
		 consp=val;
   }
   i++; if(FlagLoadX0) x0[i]=TCSCptr->Itcsc;
   if(FlagUpdateVar) TCSCptr->Itcsc=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=TCSCptr->delta_t;
   if(FlagUpdateVar) TCSCptr->delta_t=x0[i]+Dx[i];
  }

  i=NacVar+11*Ndc/2+3*Nsvc+NtcscVar;
  for(STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
   Q=STATCOMptr->Q;
   if(!strcmp(STATCOMptr->Cont,"PW") || !strcmp(STATCOMptr->Cont,"AL")){
     i++; if(FlagLoadX0) x0[i]=STATCOMptr->I;
     if(FlagUpdateVar){
        STATCOMptr->val=x0[i]+Dx[i];
        if(STATCOMptr->val<=0){
          val=fabs((STATCOMptr->val-0.0001)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s Izero %lf\n",STATCOMptr->Name,val);
          STATCOMptr->I=0.0001;}
        else if (STATCOMptr->val>=STATCOMptr->Imax && Q>0){
          val=fabs((STATCOMptr->val-STATCOMptr->Imax)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s Imax %lf\n",STATCOMptr->Name,val);
          STATCOMptr->I=STATCOMptr->Imax;}
        else if(STATCOMptr->val>=STATCOMptr->Imin && Q<0){
          val=fabs((STATCOMptr->val-STATCOMptr->Imin)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s Imin %lf\n",STATCOMptr->Name,val);
          STATCOMptr->I=STATCOMptr->Imin;}
        else {val=0;STATCOMptr->I=STATCOMptr->val;}
        if (val>consp) 
			consp=val;
     }
   }
   else if(!strcmp(STATCOMptr->Cont,"MX")){
     i++; if(FlagLoadX0) x0[i]=STATCOMptr->Vvar;
     if(FlagUpdateVar){
        STATCOMptr->val=x0[i]+Dx[i];
        if(STATCOMptr->val<=STATCOMptr->Vref){
          val=fabs((STATCOMptr->val-STATCOMptr->Vref)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s Vref_max %lf\n",STATCOMptr->Name,val);
          STATCOMptr->Vvar=STATCOMptr->Vref;}
        else {val=0;STATCOMptr->Vvar=STATCOMptr->val;}
        if (val>consp) 
			consp=val;
     }
   }
   else {
     i++; if(FlagLoadX0) x0[i]=STATCOMptr->Vvar;
     if(FlagUpdateVar){
        STATCOMptr->val=x0[i]+Dx[i];
        if(STATCOMptr->val>=STATCOMptr->Vref){
          val=fabs((STATCOMptr->val-STATCOMptr->Vref)/Dx[i]);
          if(ExistParameter('d')) fprintf(stderr,"%s Vref_min %lf\n",STATCOMptr->Name,val);
          STATCOMptr->Vvar=STATCOMptr->Vref;}
        else {val=0;STATCOMptr->Vvar=STATCOMptr->val;}
        if (val>consp) 
			consp=val;
     }
   }
   i++; if(FlagLoadX0) x0[i]=STATCOMptr->theta;
   if(FlagUpdateVar) STATCOMptr->theta=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=STATCOMptr->Vdc;
   if(FlagUpdateVar) STATCOMptr->Vdc=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=STATCOMptr->k;
   if(FlagUpdateVar) STATCOMptr->k=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=STATCOMptr->alpha;
   if(FlagUpdateVar) STATCOMptr->alpha=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=STATCOMptr->P;
   if(FlagUpdateVar) STATCOMptr->P=x0[i]+Dx[i];
   i++; if(FlagLoadX0) x0[i]=STATCOMptr->Q;
   if(FlagUpdateVar) STATCOMptr->Q=x0[i]+Dx[i];
  }
                        /* END FACTS */

  if (consp>1) consp=1;

  if(Ptr!=NULL) {
    BlPtr=Ptr; DxiMax=fabs(valp/Dparam);
    if(ExistParameter('d')) fprintf(stderr,"%s %d %s Dx=%lf Max.dx_i=%lf\n",Ptr->Type,Ptr->Num,Ptr->Name,Dx[ACvar[Ptr->N]+1],DxiMax);
  }
/*  if (FlagLoadX0 && ExistParameter('d') && DxZero!=NULL && NewNumEq<Jac->n1-1) {
      for(i=1;i<=Jac->n1-1;i++) {if (DxZero[i]) fprintf(stderr,"%d ",i);}
      fprintf(stderr,"\n",i);
    } */
  if (FlagMakeDxZero) return(consp);
  else return(count);
}

