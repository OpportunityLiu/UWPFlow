/* Homotopy Continuation Method: Jacobian */

#include "homot.h"

/* ------- Global Variables ------ */
extern VALUETYPE *Dx,Dparam,param0,*x0,*x0p,Kh,Htol,SD0,AngTr,
                 VoltageSum,VoltageSum0,DxiMax,VSF,SF,ZeroDx,Tol;
extern INDEX NewNumEq,CountSteps,NumSteps;
extern AClist *Vlist,*Vlistp;
extern bool flagReducedContinuation,flagReduceSystem,*DxZero;


/* ------------------ HFunJac ----------------------------- */
int HFunJac(bool FlagFunction,bool FlagJacobian,AreaData *Aptr,VALUETYPE *vec)
/* Add a row to the Jacobian. */
{
  ACbusData *ACptr,*ACptrp;
  DCbusData *DCptrR,*DCptrI,*DCptr;
  SVCbusData *SVCptr;                  /* FACTS */
  TCSCbusData *TCSCptr;                /* FACTS */
  STATCOMbusData *STATCOMptr;          /* FACTS */
  ElementData *Eptr;
  ElementList *ELptr;
  INDEX i,j,l;
  VALUETYPE valp;
  int val;

  val=0;
  l=Jac->n1;
  if (FlagFunction) {
    if (Bl) dF[l]=0;
    else dF[l]=Dparam*(lambda-param0-Dparam);
  }
  if (FlagJacobian) {
    if (Dparam) JacElement(Jac,l,l,Dparam);
    else JacElement(Jac,l,l,1.0);
  }
  for(ACptr=dataPtr->ACbus;ACptr!=nullptr;ACptr=ACptr->Next){
    if (ACptr->Cont!=nullptr) {
      i=ACvar[ACptr->N];
      if (strpbrk(ACptr->Type,"S")) {
        if (FlagFunction) {
          if (flagReducedContinuation && DxZero[i]) dF[i]=0;
          dF[l] +=vec[i]*(ACptr->Kg-x0[i]-vec[i]);
        }
        if (FlagJacobian) {
          if (Aptr==nullptr && fabs(vec[i])<1e-6) val=-1;
          JacElement(Jac,l,i,vec[i]);
        }
      } else {
        if (FlagFunction) {
          if (flagReducedContinuation && DxZero[i]) dF[i]=0;
          dF[l] +=vec[i]*(ACptr->Ang-x0[i]-vec[i]);
        }
        if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      }
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i+1]) dF[i+1]=0;
        dF[l] +=vec[i+1]*(ACptr->V-x0[i+1]-vec[i+1]);
      }
      if (FlagJacobian) JacElement(Jac,l,i+1,vec[i+1]);
    }
    else if(strpbrk(ACptr->Type,"L")) {
      i=ACvar[ACptr->N];
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Ang-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      if (FlagFunction) dF[l] +=vec[i+1]*(lambda-x0[i+1]-vec[i+1]);
      if (FlagJacobian) JacElement(Jac,l,i+1,vec[i+1]);
      if (FlagFunction) dF[l] +=Dparam*(ACptr->V-param0-Dparam);
    }
    else if(QRcont && strpbrk(ACptr->Type,"C")){
      i=ACvar[ACptr->N];
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Ang-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i+1]) dF[i+1]=0;
        dF[l] +=vec[i+1]*(ACptr->Qr-x0[i+1]-vec[i+1]);
      }
      if (FlagJacobian) JacElement(Jac,l,i+1,vec[i+1]);
    }
    else if(Rcont && strpbrk(ACptr->Type,"T")){
      i=ACvar[ACptr->N];
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Ang-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      for(ELptr=ACptr->Reg;ELptr!=nullptr;ELptr=ELptr->Next){
        Eptr=ELptr->Eptr;
        if(!strcmp(Eptr->Type,"R")) break;
      }
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i+1]) dF[i+1]=0;
        dF[l] +=vec[i+1]*(Eptr->Tap-x0[i+1]-vec[i+1]);
      }
      if (FlagJacobian) JacElement(Jac,l,i+1,vec[i+1]);
    }
    else if(strpbrk(ACptr->Type,"Q,S,V,Z")|| (!QRcont && strpbrk(ACptr->Type,"G"))){
      i=ACvar[ACptr->N];
      if (strpbrk(ACptr->Type,"S")) {
        if (FlagFunction) {
          if (flagReducedContinuation && DxZero[i]) dF[i]=0;
          dF[l] +=vec[i]*(ACptr->Kg-x0[i]-vec[i]);
        }
        if (FlagJacobian) {
          if (Aptr==nullptr && fabs(vec[i])<1e-6) val=-1;
          JacElement(Jac,l,i,vec[i]);
        }
      } else {
        if (FlagFunction) {
          if (flagReducedContinuation && DxZero[i]) dF[i]=0;
          dF[l] +=vec[i]*(ACptr->Ang-x0[i]-vec[i]);
        }
        if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      }
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i+1]) dF[i+1]=0;
        dF[l] +=vec[i+1]*(ACptr->Qg-x0[i+1]-vec[i+1]);
      }
      if (FlagJacobian) JacElement(Jac,l,i+1,vec[i+1]);
    }
    if(Acont && strpbrk(ACptr->Type,"A")){
      i=ACvar[ACptr->N]+2;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Kg-x0[i]-vec[i]);
      }
      if (FlagJacobian) {
        if (Aptr!=nullptr && ACptr->Area==Aptr && fabs(vec[i])<1e-6) val=-1;
        JacElement(Jac,l,i,vec[i]);
      }
    }
    if (PQcont) for(ELptr=ACptr->Reg;ELptr!=nullptr;ELptr=ELptr->Next) {
      Eptr=ELptr->Eptr;
      if(strpbrk(Eptr->Type,"PQNM")) {
        i=ACvar[ACptr->N]+1+ACptr->Ncont-Eptr->Ncont;
        if (Acont && strpbrk(ACptr->Type,"A")) i++;
        if(!strcmp(Eptr->Type,"RP")){
          if (FlagFunction) {
            if (flagReducedContinuation && DxZero[i]) dF[i]=0;
            dF[l] +=vec[i]*(Eptr->Ang-x0[i]-vec[i]);
          }
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
        else if(!strcmp(Eptr->Type,"RQ")){
          if (FlagFunction) {
            if (flagReducedContinuation && DxZero[i]) dF[i]=0;
            dF[l] +=vec[i]*(Eptr->Tap-x0[i]-vec[i]);
          }
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
        else {
          if (FlagFunction) {
            if (flagReducedContinuation && DxZero[i]) dF[i]=0;
            dF[l] +=vec[i]*(Eptr->Cvar-x0[i]-vec[i]);
          }
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
      }
    }
    if (ACptr->Gen!=nullptr) {
      i=ACptr->Gen->Nvar+1;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        if (!strpbrk(ACptr->cont,"E")) dF[l] +=vec[i]*(ACptr->Gen->Eq-x0[i]-vec[i]);
        else                                dF[l] +=vec[i]*(ACptr->Qg-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->dg-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->Vr-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->Vi-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->Ir-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->Ii-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->Vq-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->Vd-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->Iq-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        dF[l] +=vec[i]*(ACptr->Gen->Id-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      i++;
      if (FlagFunction) {
        if (flagReducedContinuation && DxZero[i]) dF[i]=0;
        if (!strpbrk(ACptr->cont,"I")) dF[l] +=vec[i]*(ACptr->Gen->Ia-x0[i]-vec[i]);
        else                           dF[l] +=vec[i]*(ACptr->Qg-x0[i]-vec[i]);
      }
      if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
    }
  }
  i=NacVar;
  for(DCptrR=dataPtr->DCbus;DCptrR!=nullptr;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if(!strcmp(DCptrR->Type,"R")){
      for (j=1;j<=2;j++) {
        if (j==1) DCptr=DCptrR;
        else DCptr=DCptrI;
        if(strcmp(DCptr->Cont1,"VD")&&strcmp(DCptr->Cont2,"VD")) {
          i++; if (FlagFunction) dF[l] +=vec[i]*(DCptr->Vd-x0[i]-vec[i]);
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
        if(strcmp(DCptr->Cont1,"AT")&&strcmp(DCptr->Cont2,"AT")) {
          i++; if (FlagFunction) dF[l] +=vec[i]*(DCptr->Tap*DCptr->Ntrf-x0[i]-vec[i]);
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
        if(strcmp(DCptr->Cont1,"AL")&&strcmp(DCptr->Cont2,"AL")) {
          i++; if (FlagFunction) dF[l] +=vec[i]*(cos(DCptr->Alfa)-x0[i]-vec[i]);
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
        if(strcmp(DCptr->Cont1,"GA")&&strcmp(DCptr->Cont2,"GA")) {
          i++; if (FlagFunction) dF[l] +=vec[i]*(cos(DCptr->Gamma)-x0[i]-vec[i]);
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
        i++; if (FlagFunction) dF[l] +=vec[i]*(DCptr->MVA-x0[i]-vec[i]);
        if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        if(strcmp(DCptr->Cont1,"PA")&&strcmp(DCptr->Cont2,"PA")) {
          i++; if (FlagFunction) dF[l] +=vec[i]*(DCptr->P-x0[i]-vec[i]);
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
        if(strcmp(DCptr->Cont1,"QA")&&strcmp(DCptr->Cont2,"QA")) {
          i++; if (FlagFunction) dF[l] +=vec[i]*(DCptr->Q-x0[i]-vec[i]);
          if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
        }
      }
      if(strcmp(DCptrR->Cont1,"ID")&&strcmp(DCptrR->Cont2,"ID")&&
         strcmp(DCptrI->Cont1,"ID")&&strcmp(DCptrI->Cont2,"ID")) {
        i++; if (FlagFunction) dF[l] +=vec[i]*(DCptrR->Id-x0[i]-vec[i]);
        if (FlagJacobian) JacElement(Jac,l,i,vec[i]);
      }
    }
  }
                                     /* FACTS */
  i=NacVar+11*Ndc/2;
  for(SVCptr=dataPtr->SVCbus;SVCptr!=nullptr;SVCptr=SVCptr->Next){
     i++; if(FlagFunction) dF[l] +=vec[i]*(SVCptr->Qsvc-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(SVCptr->Bv-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     if(!strcmp(SVCptr->Cont,"AL")){
       i++; if(FlagFunction) dF[l] +=vec[i]*(SVCptr->alpha_svc-x0[i]-vec[i]);
       if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     } else {
       i++; if(FlagFunction) dF[l] +=vec[i]*(SVCptr->Vvar-x0[i]-vec[i]);
       if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     }
  }
  i=NacVar+11*Ndc/2+3*Nsvc;
  for(TCSCptr=dataPtr->TCSCbus;TCSCptr!=nullptr;TCSCptr=TCSCptr->Next){
     i++; if(FlagFunction) dF[l] +=vec[i]*(TCSCptr->Ptcsc-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(TCSCptr->Qtcsck-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(TCSCptr->Qtcscm-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(TCSCptr->Be-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(TCSCptr->alpha_tcsc-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(TCSCptr->Itcsc-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(TCSCptr->delta_t-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
  }
  i=NacVar+11*Ndc/2+3*Nsvc+NtcscVar;
  for(STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=nullptr;STATCOMptr=STATCOMptr->Next){
     if (!strcmp(STATCOMptr->Cont,"PW") || !strcmp(STATCOMptr->Cont,"AL")) {
       i++; if(FlagFunction) dF[l] +=vec[i]*(STATCOMptr->I-x0[i]-vec[i]);
       if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     } else {
       i++; if(FlagFunction) dF[l] +=vec[i]*(STATCOMptr->Vvar-x0[i]-vec[i]);
       if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     }
     i++; if(FlagFunction) dF[l] +=vec[i]*(STATCOMptr->theta-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(STATCOMptr->Vdc-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(STATCOMptr->k-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(STATCOMptr->alpha-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(STATCOMptr->P-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
     i++; if(FlagFunction) dF[l] +=vec[i]*(STATCOMptr->Q-x0[i]-vec[i]);
     if(FlagJacobian)JacElement(Jac,l,i,vec[i]);
  }
                                  /* END FACTS */
  return(val);
}


