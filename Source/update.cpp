/* Update AC/DC variables. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"

void UpdateACvar(VALUETYPE cons,INDEX j,bool Limits,bool Recover);
void UpdateDCvar(VALUETYPE cons,INDEX j,bool Limits);

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX MaxIter,Nac,NacEl,NregPQ,NregV,Ndc,Nslack,Nvolt,Narea,NacVar,Bl,NZbuses,NXbuses;
extern INDEX *ACvar;
extern VALUETYPE *dx,*dF,tol,Tol,Sn,lambda,*x0;
extern VALUETYPE K1,K2,MaxdFi,alpha;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern IntegerVector *RowPer,*ColPer;
extern bool Acont,PQcont,QRcont,Rcont,Xcont,
               PQlim,Tlim,Qlim,Vlim,Elim,Ilim,Zlim,Xlim,
               flagH,flagPoC,flagL,flagR,flagBS,flagPgMax,flagSmax;


/* -----------------UpdateACvar ---------------------------- */
void UpdateACvar(VALUETYPE cons,INDEX j,bool Limits,bool Recover)
{
  ACbusData *ACptr,*ACptrp,*BEptr;
  AClist *ALptr;
  ElementList *ELptr;
  ElementData *Eptr;
  VALUETYPE Pg,DPg,Qm,Pmax,PgMax;
  INDEX k;

  for(ALptr=dataPtr->KGbus;ALptr!=nullptr;ALptr=ALptr->Next) {
    ACptr=ALptr->AC;
    k=ACvar[ACptr->N];
    if (strpbrk(ACptr->Type,"S")) {
      ACptr->Kg=ACptr->Kg+cons*dx[k];
      BEptr=ACptr;
    }
    else if (Acont) ACptr->Kg=ACptr->Kg+cons*dx[k+2];
  }
  for (ACptr=dataPtr->ACbus; ACptr!=nullptr; ACptr=ACptr->Next){
    if (ACptr->Area!=nullptr) BEptr=ACptr->Area->Slack;
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
      if (ACptr->flagPgMax!=1) ACptr->flagPgMax++;
    } else ACptr->flagPgMax=0;
    ACptr->PG=Pg;
    ACptr->DPG=DPg;
    ACptr->Pmax=PgMax;
    Qm=ACptr->Smax*ACptr->Smax-Pg*Pg;
    if (!flagSmax && Qm>0) {
      Qm=sqrt(Qm);
      if (ACptr->Qmax<Qm) ACptr->Max=ACptr->Qmax;
      else ACptr->Max=Qm;
      if (ACptr->Qmin>-Qm) ACptr->Min=ACptr->Qmin;
      else ACptr->Min=-Qm;
    } else {
      ACptr->Max=ACptr->Qmax;
      ACptr->Min=ACptr->Qmin;
    }
    if (!strcmp(ACptr->Type,"B") || !strcmp(ACptr->Type,"BA")) {
      k=ACvar[ACptr->N];
      ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (j==0) ACptr->val=ACptr->V;
      ACptr->val=ACptr->val+cons*dx[k+1];
      if (ACptr->val<=0) ACptr->V=0.00001;
      else ACptr->V=ACptr->val;
    }
    else if (strpbrk(ACptr->Type,"V")) {
      k=ACvar[ACptr->N];
      ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (!strpbrk(ACptr->cont,"V")) {
        if (j==0) ACptr->val=ACptr->V;
        ACptr->val=ACptr->val+cons*dx[k+1];
        if (Vlim && ACptr->val>ACptr->Vmax && Limits)      ACptr->V=ACptr->Vmax;
        else if (Vlim && ACptr->val<ACptr->Vmin && Limits) ACptr->V=ACptr->Vmin;
        else if (ACptr->val<=0)                            ACptr->V=0.00001;
        else                                               ACptr->V=ACptr->val;
      } else {
        if (j==0) ACptr->val=ACptr->Qg;
        ACptr->val=ACptr->val+cons*dx[k+1];
        if (Recover && ACptr->V>=ACptr->Vmax && ACptr->val>ACptr->VCont && Limits)      ACptr->Qg=ACptr->VCont;
        else if (Recover && ACptr->V<=ACptr->Vmin && ACptr->val<ACptr->VCont && Limits) ACptr->Qg=ACptr->VCont;
        else                                                                            ACptr->Qg=ACptr->val;
      }
    }
    else if (strpbrk(ACptr->Type,"L")){
      k=ACvar[ACptr->N];
      ACptr->Ang=ACptr->Ang+cons*dx[k];
      lambda=lambda+cons*dx[k+1];
      if (flagH) {
        if (j==0) ACptr->val=ACptr->V;
        ACptr->val=ACptr->val+cons*dx[Jac->n1];
        if (ACptr->val<=0) ACptr->V=0.00001;
        else ACptr->V=ACptr->val;
      }
    }
    else if (QRcont && strpbrk(ACptr->Type,"G")){
      k=ACvar[ACptr->N];
      if (!strpbrk(ACptr->Type,"S")) ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (j==0) ACptr->val=ACptr->V;
      ACptr->val=ACptr->val+cons*dx[k+1];
      if (ACptr->val<=0) ACptr->V=0.00001;
      else ACptr->V=ACptr->val;
    }
    else if (strpbrk(ACptr->Type,"Q,S") || (!QRcont && strpbrk(ACptr->Type,"G"))){
      k=ACvar[ACptr->N];
      if (!strpbrk(ACptr->Type,"S")) ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (strpbrk(ACptr->cont,"V")) {
        if (j==0) ACptr->val=ACptr->Qg;
        ACptr->val=ACptr->val+cons*dx[k+1];
        if (Qlim && ACptr->val>ACptr->Max && Limits)      ACptr->Qg=ACptr->Max;
        else if (Qlim && ACptr->val<ACptr->Min && Limits) ACptr->Qg=ACptr->Min;
        else                                              ACptr->Qg=ACptr->val;
      }
      else {
        if (j==0) ACptr->val=ACptr->V;
        ACptr->val=ACptr->val+cons*dx[k+1];
        if (ACptr->val<=0) ACptr->V=0.00001;
        else if (strpbrk(ACptr->cont,"Q")) {
          if (Recover && ACptr->Qg>=ACptr->Max && ACptr->val>ACptr->VCont && Limits)      ACptr->V=ACptr->VCont;
          else if (Recover && ACptr->Qg<=ACptr->Min && ACptr->val<ACptr->VCont && Limits) ACptr->V=ACptr->VCont;
          else                                                                            ACptr->V=ACptr->val;
        }
        else if (strpbrk(ACptr->cont,"E")) {
          if (Recover && ACptr->Gen->Eq>=ACptr->Gen->EqMax && ACptr->val>ACptr->VCont && Limits)      ACptr->V=ACptr->VCont;
          else if (Recover && ACptr->Gen->Eq<=ACptr->Gen->EqMin && ACptr->val<ACptr->VCont && Limits) ACptr->V=ACptr->VCont;
          else                                                                                        ACptr->V=ACptr->val;
        }
        else if (strpbrk(ACptr->cont,"I")) {
          if (Recover && ACptr->Gen->Ia>=ACptr->Gen->IaMax && ACptr->val>ACptr->VCont && Limits) ACptr->V=ACptr->VCont;
          else                                                                                   ACptr->V=ACptr->val;
        }
        else ACptr->V=ACptr->val;
      }
    }
    else if (strpbrk(ACptr->Type,"X")) {
      k=ACvar[ACptr->N];
      ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (j==0) ACptr->val=ACptr->V;
      ACptr->val=ACptr->val+cons*dx[k+1];
      if (ACptr->val<=0)                                   ACptr->V=0.00001;
      else if (strpbrk(ACptr->cont,"X")) {
        if (Xlim && ACptr->val>ACptr->Vmax && Limits)      ACptr->V=ACptr->Vmax;
        else if (Xlim && ACptr->val<ACptr->Vmin && Limits) ACptr->V=ACptr->Vmin;
        else                                               ACptr->V=ACptr->val;
      }
      else if (strpbrk(ACptr->cont,"M")) {
        if (Recover && ACptr->val<ACptr->Vmax && Limits)   ACptr->V=ACptr->Vmax;
        else                                               ACptr->V=ACptr->val;
      }
      else {
        if (Recover && ACptr->val>ACptr->Vmin && Limits)   ACptr->V=ACptr->Vmin;
        else                                               ACptr->V=ACptr->val;
      }
    }
    else if (strpbrk(ACptr->Type,"Z")){
      k=ACvar[ACptr->N];
      ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (strpbrk(ACptr->cont,"V")) {
        if (j==0) ACptr->val=ACptr->Qg;
        ACptr->val=ACptr->val+cons*dx[k+1];
        if (Zlim && ACptr->val>ACptr->Qmax && Limits)      ACptr->Qg=ACptr->Qmax;
        else if (Zlim && ACptr->val<ACptr->Qmin && Limits) ACptr->Qg=ACptr->Qmin;
        else                                               ACptr->Qg=ACptr->val;
      }
      else {
        if (j==0) ACptr->val=ACptr->V;
        ACptr->val=ACptr->val+cons*dx[k+1];
        ACptr->Qg=ACptr->val*ACptr->val*ACptr->Bz;
        if (ACptr->val<=0)                                                               ACptr->V=0.00001;
        else if (Recover && ACptr->Qg>=ACptr->Qmax && ACptr->val>ACptr->VCont && Limits) ACptr->V=ACptr->VCont;
        else if (Recover && ACptr->Qg<=ACptr->Qmin && ACptr->val<ACptr->VCont && Limits) ACptr->V=ACptr->VCont;
        else                                                                             ACptr->V=ACptr->val;
        ACptr->Qg=ACptr->V*ACptr->V*ACptr->Bz;
      }
    }
    else if (strpbrk(ACptr->Type,"C")) {
      k=ACvar[ACptr->N];
      ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (!QRcont) {
        if (j==0) ACptr->val=ACptr->V;
        ACptr->val=ACptr->val+cons*dx[k+1];
        if (ACptr->val<=0) ACptr->V=0.00001;
        else ACptr->V=ACptr->val;
      }
      else if (ACptr->Kbg<1) {
        ACptrp=ACptr->ContBus->AC;
        if (j==0) ACptr->val=ACptr->V;
        ACptr->val=ACptr->val+cons*dx[k+1];
        if (ACptr->val<=0) ACptr->V=0.00001;
        else if (strpbrk(ACptrp->cont,"Q")) {
          if (Recover && ACptrp->Qg>=ACptrp->Max && ACptr->val>ACptr->VCont && Limits)      ACptr->V=ACptr->VCont;
          else if (Recover && ACptrp->Qg<=ACptrp->Min && ACptr->val<ACptr->VCont && Limits) ACptr->V=ACptr->VCont;
          else                                                                              ACptr->V=ACptr->val;
        }
        else if (strpbrk(ACptrp->cont,"E")) {
          if (Recover && ACptrp->Gen->Eq>=ACptrp->Gen->EqMax && ACptr->val>ACptr->VCont && Limits)      ACptr->V=ACptr->VCont;
          else if (Recover && ACptrp->Gen->Eq<=ACptrp->Gen->EqMin && ACptr->val<ACptr->VCont && Limits) ACptr->V=ACptr->VCont;
          else                                                                                          ACptr->V=ACptr->val;
        }
        else if (strpbrk(ACptrp->cont,"I")) {
          if (Recover && ACptrp->Gen->Ia>=ACptrp->Gen->IaMax && ACptr->val>ACptr->VCont && Limits) ACptr->V=ACptr->VCont;
          else                                                                                     ACptr->V=ACptr->val;
        }
        else ACptr->V=ACptr->val;
      }
      else
        ACptr->Qr=ACptr->Qr+cons*dx[k+1];
    }
    else if (strpbrk(ACptr->Type,"R")) {
      k=ACvar[ACptr->N];
      ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (j==0) ACptr->val=ACptr->V;
      ACptr->val=ACptr->val+cons*dx[k+1];
      if (Rcont && ACptr->val>ACptr->Vmax && Limits)      ACptr->V=ACptr->Vmax;
      else if (Rcont && ACptr->val<ACptr->Vmin && Limits) ACptr->V=ACptr->Vmin;
      else if (ACptr->val<=0)                             ACptr->V=0.00001;
      else                                                ACptr->V=ACptr->val;
    }
    else if (strpbrk(ACptr->Type,"T")){
      k=ACvar[ACptr->N];
      ACptr->Ang=ACptr->Ang+cons*dx[k];
      if (Rcont) {
        for(ELptr=ACptr->Reg;ELptr!=nullptr;ELptr=ELptr->Next){
          Eptr=ELptr->Eptr;
          if (!strcmp(Eptr->Type,"R")) {
            if (j==0) Eptr->val=Eptr->Tap;
            Eptr->val=Eptr->val+cons*dx[k+1];
            if (Tlim && Eptr->val<1/Eptr->Tmax+0.00001 && Limits) Eptr->Tap=1/Eptr->Tmax;
            else if (Tlim && Eptr->val>1/Eptr->Tmin-0.00001 && Limits) Eptr->Tap=1/Eptr->Tmin;
            else Eptr->Tap=Eptr->val;
          }
        }
      } else {
        if (j==0) ACptr->val=ACptr->V;
        ACptr->val=ACptr->val+cons*dx[k+1];
        if (ACptr->val<=0) ACptr->V=0.00001;
        else ACptr->V=ACptr->val;
      }
    }
    if (PQcont) for(ELptr=ACptr->Reg;ELptr!=nullptr;ELptr=ELptr->Next) {
      Eptr=ELptr->Eptr;
      if(strpbrk(Eptr->Type,"PQMN")) {
        k=ACvar[ACptr->N]+1;
        if(Acont && strpbrk(ACptr->Type,"A")) k=k+1;
        k=k+Eptr->Cont->Ncont-Eptr->Ncont;
        if (!strcmp(Eptr->Type,"RP")) {
          if (j==0) Eptr->val=Eptr->Ang;
          Eptr->val=Eptr->val+cons*dx[k];
          if (Tlim && Eptr->val>Eptr->Tmax && Limits) Eptr->Ang=Eptr->Tmax;
          else if (Tlim && Eptr->val<Eptr->Tmin && Limits) Eptr->Ang=Eptr->Tmin;
          else Eptr->Ang=Eptr->val;
        }
        else if (!strcmp(Eptr->Type,"RQ")) {
          if (j==0) Eptr->val=Eptr->Tap;
          Eptr->val=Eptr->val+cons*dx[k];
          if (Tlim && Eptr->val<1/Eptr->Tmax+0.00001 && Limits) Eptr->Tap=1/Eptr->Tmax;
          else if (Tlim && Eptr->val>1/Eptr->Tmin-0.00001 && Limits) Eptr->Tap=1/Eptr->Tmin;
          else Eptr->Tap=Eptr->val;
        }
        else if (strpbrk(Eptr->Type,"MN")) {
          if (j==0) Eptr->val=Eptr->Cvar;
          Eptr->val=Eptr->val+cons*dx[k];
          if (PQlim && Eptr->val>Eptr->Max && Limits) Eptr->Cvar=Eptr->Max;
          else if (PQlim && Eptr->val<Eptr->Min && Limits) Eptr->Cvar=Eptr->Min;
          else Eptr->Cvar=Eptr->val;
        }
        else Eptr->Cvar=Eptr->Cvar+cons*dx[k];
      }
    }
    if (ACptr->Gen!=nullptr) {
      k=ACptr->Gen->Nvar;
      if (!strpbrk(ACptr->cont,"E")) {
        if (j==0) ACptr->vals=ACptr->Gen->Eq;
        ACptr->vals=ACptr->vals+cons*dx[k+1];
        if (Elim && ACptr->vals>ACptr->Gen->EqMax && Limits)      ACptr->Gen->Eq=ACptr->Gen->EqMax;
        else if (Elim && ACptr->vals<ACptr->Gen->EqMin && Limits) ACptr->Gen->Eq=ACptr->Gen->EqMin;
        else if (ACptr->vals<=0)                                  ACptr->Gen->Eq=0.00001;
        else                                                      ACptr->Gen->Eq=ACptr->vals;
      }
      else {
        if (j==0) ACptr->vals=ACptr->Qg;
        ACptr->vals=ACptr->vals+cons*dx[k+1];
        if (strpbrk(ACptr->Type,"V")) {
          if (Recover && ACptr->Gen->Eq>=ACptr->Gen->EqMax && ACptr->vals>ACptr->VCont && Limits)      ACptr->Qg=ACptr->VCont;
          else if (Recover && ACptr->Gen->Eq<=ACptr->Gen->EqMin && ACptr->vals<ACptr->VCont && Limits) ACptr->Qg=ACptr->VCont;
          else                                                                                         ACptr->Qg=ACptr->vals;
        }
        else {
          if (Qlim && ACptr->vals>ACptr->Max && Limits)      ACptr->Qg=ACptr->Max;
          else if (Qlim && ACptr->vals<ACptr->Min && Limits) ACptr->Qg=ACptr->Min;
          else                                               ACptr->Qg=ACptr->vals;
        }
      }
      ACptr->Gen->dg=ACptr->Gen->dg+cons*dx[k+2];
      ACptr->Gen->Vr=ACptr->Gen->Vr+cons*dx[k+3];
      ACptr->Gen->Vi=ACptr->Gen->Vi+cons*dx[k+4];
      ACptr->Gen->Ir=ACptr->Gen->Ir+cons*dx[k+5];
      ACptr->Gen->Ii=ACptr->Gen->Ii+cons*dx[k+6];
      ACptr->Gen->Vq=ACptr->Gen->Vq+cons*dx[k+7];
      ACptr->Gen->Vd=ACptr->Gen->Vd+cons*dx[k+8];
      ACptr->Gen->Iq=ACptr->Gen->Iq+cons*dx[k+9];
      ACptr->Gen->Id=ACptr->Gen->Id+cons*dx[k+10];
      if (!strpbrk(ACptr->cont,"I")) {
        if (j==0) ACptr->valt=ACptr->Gen->Ia;
        ACptr->valt=ACptr->valt+cons*dx[k+11];
        if (Ilim && ACptr->valt>ACptr->Gen->IaMax && Limits) ACptr->Gen->Ia=ACptr->Gen->IaMax;
        else if (ACptr->valt<=0)                             ACptr->Gen->Ia=0.00001;
        else                                                 ACptr->Gen->Ia=ACptr->valt;
      }
      else {
        if (j==0) ACptr->valt=ACptr->Qg;
        ACptr->valt=ACptr->valt+cons*dx[k+11];
        if (strpbrk(ACptr->Type,"V")) {
          if (Recover && ACptr->Gen->Ia>=ACptr->Gen->IaMax && ACptr->valt>ACptr->VCont && Limits) ACptr->Qg=ACptr->VCont;
          else                                                                                    ACptr->Qg=ACptr->valt;
        }
        else {
          if (Qlim && ACptr->valt>ACptr->Max && Limits)      ACptr->Qg=ACptr->Max;
          else if (Qlim && ACptr->valt<ACptr->Min && Limits) ACptr->Qg=ACptr->Min;
          else                                               ACptr->Qg=ACptr->valt;
        }
      }
    }
  }

  for (ACptr=dataPtr->ACbus; ACptr!=nullptr; ACptr=ACptr->Next) {
    if (QRcont && strpbrk(ACptr->Type,"G")  && strpbrk(ACptr->cont,"V")){
      k=ACvar[ACptr->Cont->N];
      if (ACptr->Cont->Qr>=0) ACptr->Kbg=ACptr->Kbg1;
      else                    ACptr->Kbg=ACptr->Kbg2;
      ACptr->valp=ACptr->Kbg*ACptr->Cont->Qr;
      if (Qlim && ACptr->valp>ACptr->Max && Limits)      ACptr->Qg=ACptr->Max;
      else if (Qlim && ACptr->valp<ACptr->Min && Limits) ACptr->Qg=ACptr->Min;
      else                                               ACptr->Qg=ACptr->valp;
    }
  }

  if (flagH && !Bl) lambda=lambda+cons*dx[Jac->n1];
}


/* -----------------UpdateDCvar ---------------------------- */
void UpdateDCvar(VALUETYPE cons,INDEX j,bool Limits)
{
  DCbusData *DCptrR,*DCptrI,*DCptr;
  INDEX i,k,m;

  k=NacVar;
  for (DCptrR=dataPtr->DCbus;DCptrR!=nullptr;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if (!strcmp(DCptrR->Type,"R")){
      for (m=1;m<=2;m++) {
        i= -1;
        if(m==1) DCptr=DCptrR;
        else DCptr=DCptrI;
        if (strcmp(DCptr->Cont1,"VD") && strcmp(DCptr->Cont2,"VD")){
          i++;
          if (j==0) DCptr->val[i]=DCptr->Vd;
          DCptr->val[i]=DCptr->val[i]+cons*dx[++k];
          if (DCptr->VdN && ((DCptr->Tap>=DCptr->TapMax && DCptr->val[i]>DCptr->VdN)||
              (DCptr->Tap<=DCptr->TapMin && DCptr->val[i]<DCptr->VdN) ||
              (strpbrk(DCptr->Cont1,"IP")&&DCptr->Tap>DCptr->TapMin&&DCptr->val[i]>DCptr->VdN)))
             DCptr->Vd=DCptr->VdN;
          else DCptr->Vd=DCptr->val[i];}
        if (strcmp(DCptr->Cont1,"AT") && strcmp(DCptr->Cont2,"AT")){
          i++;
          if (j==0) DCptr->val[i]=DCptr->Tap;
          DCptr->val[i]=DCptr->val[i]+cons*dx[++k]/DCptr->Ntrf;
          if (DCptr->val[i]>DCptr->TapMax && Limits) DCptr->Tap=DCptr->TapMax;
          else if (DCptr->val[i]<DCptr->TapMin && Limits) DCptr->Tap=DCptr->TapMin;
          else DCptr->Tap=DCptr->val[i];}
        if (strcmp(DCptr->Cont1,"AL") && strcmp(DCptr->Cont2,"AL")){
          i++;
          if (j==0) DCptr->val[i]=cos(DCptr->Alfa);
          DCptr->val[i]=DCptr->val[i]+cons*dx[++k];
          if (DCptr->val[i]<cos(DCptr->AlfaMax)) DCptr->Alfa=DCptr->AlfaMax;
          else if (DCptr->val[i]>cos(DCptr->AlfaMin)) DCptr->Alfa=DCptr->AlfaMin;
          else if (DCptr->AlfaN && ((DCptr->Tap>=DCptr->TapMax && DCptr->val[i]<cos(DCptr->AlfaN))||
                   (DCptr->Tap<=DCptr->TapMin && DCptr->val[i]>cos(DCptr->AlfaN)))) DCptr->Alfa=DCptr->AlfaN;
          else DCptr->Alfa=acos(DCptr->val[i]);}
        if (strcmp(DCptr->Cont1,"GA") && strcmp(DCptr->Cont2,"GA")){
          i++;
          if (j==0) DCptr->val[i]=cos(DCptr->Gamma);
          DCptr->val[i]=DCptr->val[i]+cons*dx[++k];
          if (fabs(DCptr->val[i])>cos(DCptr->GammaMin)) DCptr->Gamma=DCptr->GammaMin;
          else DCptr->Gamma=acos(DCptr->val[i]);}
        DCptr->MVA=DCptr->MVA+cons*dx[++k];
        if (strcmp(DCptr->Cont1,"PA") && strcmp(DCptr->Cont2,"PA"))
          DCptr->P=DCptr->P+cons*dx[++k];
        if (strcmp(DCptr->Cont1,"QA") && strcmp(DCptr->Cont2,"QA"))
          DCptr->Q=DCptr->Q+cons*dx[++k];
      }
      if (strcmp(DCptrR->Cont1,"ID") && strcmp(DCptrR->Cont2,"ID") &&
          strcmp(DCptrI->Cont1,"ID") && strcmp(DCptrI->Cont2,"ID"))
        DCptrI->Id=DCptrR->Id=DCptrR->Id+cons*dx[++k];
    }
  }
}

