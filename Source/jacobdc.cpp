/* DC Function and Jacobian. */

#include "jacob.h"



/* ------------------ DCFunJac ----------------------------- */
#ifdef ANSIPROTO
BOOLEAN DCFunJac(SparseMatrix *Mptr,BOOLEAN flagF,BOOLEAN flagJ)
#else
BOOLEAN DCFunJac(Mptr,flagF,flagJ)
BOOLEAN flagF,flagJ;
SparseMatrix *Mptr;
#endif
/* Construct the DC part of the Jacobian. */
{
  INDEX i,j,k,l,kp,lp,m,n;
  DCbusData *DCptrR,*DCptrI;
  ACbusData *BEptr;
  VALUETYPE Pa1,Pa2,Sa1,Sa2;
  VALUETYPE Vdr,Vr,ar,cosar,cosgr,Xcr,Sr,Pr,Qr,Dr,Id;
  VALUETYPE Vdi,Vi,ai,cosai,cosgi,Xci,Si,Pi,Qi,Di,Rd;
  BOOLEAN dVr=FALSE,dVi=FALSE;

  i=NacVar;
  for(DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if (!strcmp(DCptrR->Type,"R")){
      Id=DCptrR->Id;                Rd=DCptrR->Rd;
      Vdr=DCptrR->Vd;               Vdi=DCptrI->Vd;
      k=ACvar[DCptrR->AC->N];       l=ACvar[DCptrI->AC->N];
      kp=k+1;                       lp=l+1;
      if (DCptrR->AC->Cont!=NULL) dVr=TRUE;
      else dVr=FALSE;
      if (DCptrI->AC->Cont!=NULL) dVi=TRUE;
      else dVi=FALSE;
      Vr=DCptrR->AC->V;                 Vi=DCptrI->AC->V;
      ar=DCptrR->Tap*DCptrR->Ntrf;  ai=DCptrI->Tap*DCptrI->Ntrf;
      Xcr=DCptrR->Xc;               Xci=DCptrI->Xc;
      cosar=cos(DCptrR->Alfa);      cosai=cos(DCptrI->Alfa);
      cosgr=cos(DCptrR->Gamma);     cosgi=cos(DCptrI->Gamma);
      Pr=DCptrR->P;                 Pi=DCptrI->P;
      Qr=DCptrR->Q;                 Qi=DCptrI->Q;
      Sr=DCptrR->MVA;               Si=DCptrI->MVA;
      Dr=Sr*Sr-Pr*Pr;               Di=Si*Si-Pi*Pi;
      if (Dr<=0 || Di<=0) return(TRUE);
      if (Acont && DCptrR->Meter!=NULL) {
        if(DCptrR->Meter==DCptrR) {
          Pa1= -Pr;  Pa2=Pr;
          Sa1= -1;   Sa2=1;
        }
        else {
          Pa1=Pi;  Pa2= -Pi;
          Sa1=1;   Sa2= -1;
        }
        BEptr=DCptrR->Area->Slack;
        if (!strpbrk(BEptr->Type,"S"))  m=ACvar[BEptr->N]+2;
        else m=0;
        BEptr=DCptrI->Area->Slack;
        if (!strpbrk(BEptr->Type,"S"))  n=ACvar[BEptr->N]+2;
        else n=0;
      }
      else m=n=0;
      if (flagF) {
        dF[i+1]=Vdr-K1*ar*Vr*cosar+K2*Xcr*Id;
        dF[i+2]=Sr-K1*ar*Vr*Id;
        dF[i+3]=Pr+Vdr*Id;
        dF[i+4]=Qr+sqrt(Sr*Sr-Pr*Pr);
        dF[i+5]=Vdi-K1*ai*Vi*cosgi+K2*Xci*Id;
        dF[i+6]=Si-K1*ai*Vi*Id;
        dF[i+7]=Pi-Vdi*Id;
        dF[i+8]=Qi+sqrt(Si*Si-Pi*Pi);
        dF[i+9] =sqrt(2.0)*Xcr*Id-ar*Vr*(cosar+cosgr);
        dF[i+10]=sqrt(2.0)*Xci*Id-ai*Vi*(cosai+cosgi);
        dF[i+11]=Vdr-Vdi-Rd*Id;
        dF[k]=dF[k]+Pr;
        dF[k+1]=dF[k+1]+Qr;
        dF[l]=dF[l]+Pi;
        dF[l+1]=dF[l+1]+Qi;
        if (m!=0) dF[m]=dF[m]-Pa1;
        if (n!=0) dF[n]=dF[n]-Pa2;
      }
      if (flagJ) {
        j=i;
        if (dVr){
          JacElement(Mptr,i+1,kp,-K1*ar*cosar);
          JacElement(Mptr,i+2,kp,-K1*ar*Id);
          JacElement(Mptr,i+9,kp,-ar*(cosar+cosgr));
        } else if (flagH && strpbrk(DCptrR->AC->Type,"L")) {
          JacElement(Mptr,i+1,Mptr->n1,-K1*ar*cosar);
          JacElement(Mptr,i+2,Mptr->n1,-K1*ar*Id);
          JacElement(Mptr,i+9,Mptr->n1,-ar*(cosar+cosgr));
        } else {
          JacElement(Mptr,i+1,kp,0.);
          JacElement(Mptr,i+2,kp,0.);
          JacElement(Mptr,i+9,kp,0.);
        }
        if (strcmp(DCptrR->Cont1,"VD") && strcmp(DCptrR->Cont2,"VD")){
          JacElement(Mptr,i+1,++j,1.0);
          JacElement(Mptr,i+3,j,Id);
          JacElement(Mptr,i+11,j,1.0);
        }
        if (strcmp(DCptrR->Cont1,"AT") && strcmp(DCptrR->Cont2,"AT")){
          JacElement(Mptr,i+1,++j,-K1*Vr*cosar);
          JacElement(Mptr,i+2,j,-K1*Vr*Id);
          JacElement(Mptr,i+9,j,-Vr*(cosar+cosgr));
        }
        if (strcmp(DCptrR->Cont1,"AL") && strcmp(DCptrR->Cont2,"AL")){
          JacElement(Mptr,i+1,++j,-K1*ar*Vr);
          JacElement(Mptr,i+9,j,-ar*Vr);
        }
        if (strcmp(DCptrR->Cont1,"GA") && strcmp(DCptrR->Cont2,"GA")) JacElement(Mptr,i+9,++j,-ar*Vr);
        JacElement(Mptr,i+2,++j,1.0);
        JacElement(Mptr,i+4,j,Sr/sqrt(Sr*Sr-Pr*Pr));
        if (strcmp(DCptrR->Cont1,"PA") && strcmp(DCptrR->Cont2,"PA")){
          JacElement(Mptr,i+3,++j,1.0);
          JacElement(Mptr,i+4,j,-Pr/sqrt(Sr*Sr-Pr*Pr));
          JacElement(Mptr,k,j,1.0);
          if (m!=0 && DCptrR->Meter==DCptrR) JacElement(Mptr,m,j,-Sa1);
          if (n!=0 && DCptrR->Meter==DCptrR) JacElement(Mptr,n,j,-Sa2);
        }
        if (strcmp(DCptrR->Cont1,"QA") && strcmp(DCptrR->Cont2,"QA")){
          JacElement(Mptr,i+4,++j,1.0);
          JacElement(Mptr,k+1,j,1.0);
        }
        if (dVi){
          JacElement(Mptr,i+5,lp,-K1*ai*cosgi);
          JacElement(Mptr,i+6,lp,-K1*ai*Id);
          JacElement(Mptr,i+10,lp,-ai*(cosai+cosgi));
        } else if (flagH && strpbrk(DCptrI->AC->Type,"L")) {
          JacElement(Mptr,i+5,Mptr->n1,-K1*ai*cosgi);
          JacElement(Mptr,i+6,Mptr->n1,-K1*ai*Id);
          JacElement(Mptr,i+10,Mptr->n1,-ai*(cosai+cosgi));
        } else {
          JacElement(Mptr,i+5,lp,0.);
          JacElement(Mptr,i+6,lp,0.);
          JacElement(Mptr,i+10,lp,0.);
        }
        if (strcmp(DCptrI->Cont1,"VD") && strcmp(DCptrI->Cont2,"VD")){
          JacElement(Mptr,i+5,++j,1.0);
          JacElement(Mptr,i+7,j,-Id);
          JacElement(Mptr,i+11,j,-1.0);
        }
        if (strcmp(DCptrI->Cont1,"AT") && strcmp(DCptrI->Cont2,"AT")){
          JacElement(Mptr,i+5,++j,-K1*Vi*cosgi);
          JacElement(Mptr,i+6,j,-K1*Vi*Id);
          JacElement(Mptr,i+10,j,-Vi*(cosai+cosgi));
        }
        if (strcmp(DCptrI->Cont1,"AL") && strcmp(DCptrI->Cont2,"AL"))
          JacElement(Mptr,i+10,++j,-ai*Vi);
        if (strcmp(DCptrI->Cont1,"GA") && strcmp(DCptrI->Cont2,"GA")){
          JacElement(Mptr,i+5,++j,-K1*ai*Vi);
          JacElement(Mptr,i+10,j,-ai*Vi);
        }
        JacElement(Mptr,i+6,++j,1.0);
        JacElement(Mptr,i+8,j,Si/sqrt(Si*Si-Pi*Pi));
        if (strcmp(DCptrI->Cont1,"PA") && strcmp(DCptrI->Cont2,"PA")){
          JacElement(Mptr,i+7,++j,1.0);
          JacElement(Mptr,i+8,j,-Pi/sqrt(Si*Si-Pi*Pi));
          JacElement(Mptr,l,j,1.0);
          if (m!=0 && DCptrI->Meter==DCptrI) JacElement(Mptr,m,j,-Sa1);
          if (n!=0 && DCptrI->Meter==DCptrI) JacElement(Mptr,n,j,-Sa2);
        }
        if (strcmp(DCptrI->Cont1,"QA") && strcmp(DCptrI->Cont2,"QA")){
          JacElement(Mptr,i+8,++j,1.0);
          JacElement(Mptr,l+1,j,1.0);
        }
        if (strcmp(DCptrR->Cont1,"ID") && strcmp(DCptrR->Cont2,"ID") &&
            strcmp(DCptrI->Cont1,"ID") && strcmp(DCptrI->Cont2,"ID")){
          JacElement(Mptr,i+1,++j,K2*Xcr);
          JacElement(Mptr,i+2,j,-K1*ar*Vr);
          JacElement(Mptr,i+3,j,Vdr);
          JacElement(Mptr,i+5,j,K2*Xci);
          JacElement(Mptr,i+6,j,-K1*ai*Vi);
          JacElement(Mptr,i+7,j,-Vdi);
          JacElement(Mptr,i+9,j,sqrt(2.0)*Xcr);
          JacElement(Mptr,i+10,j,sqrt(2.0)*Xci);
          JacElement(Mptr,i+11,j,-Rd);
        }
      }
      i=i+11;
    }
  }
  return(FALSE);
}
