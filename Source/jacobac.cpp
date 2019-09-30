/* AC Function and Jacobian. */

#include "jacob.h"

/* ------------------ ACFunJac ----------------------------- */
#ifdef ANSIPROTO
AreaData *ACFunJac(SparseMatrix *Mptr,int *val,BOOLEAN flagF,BOOLEAN flagJ,BOOLEAN flagP)
#else
AreaData *ACFunJac(Mptr,val,flagF,flagJ,flagP)
SparseMatrix *Mptr;
int *val;
BOOLEAN flagF,flagJ,flagP;
#endif
/* Construct the AC part of the Jacobian. */
{
  ACbusData *ACptr,*To,*From,*BEptr;
  AClist *ALptr;
  AreaData *Aptr;
  ElementList *ELptr;
  ElementData *Eptr;
  VALUETYPE Pi,Qi,Vi,Vj,di,dj,SPij,dPijd,dPiid,dQijd,dQiid,val1;
  VALUETYPE gij,bij,gsij,bsij,SQij,dPijv,dPiiv,dQijv,dQiiv,val2;
  VALUETYPE Ra,Xd,Xq,Qg,Eq,dg,Vr,Vim,Ir,Iim,Vq,Vd,Iq,Id,Ia;
  VALUETYPE Pl,Ql,Pg,DPg,SPg=0;
  INDEX i,j,k,l;

  if (val!=NULL) *val=0;
  for(ALptr=dataPtr->KGbus;ALptr!=NULL;ALptr=ALptr->Next) {
    BEptr=ALptr->AC;
    if(strpbrk(ALptr->AC->Type,"S"))  break;
  }
  if (Acont) for(Aptr=dataPtr->Area;Aptr!=NULL;Aptr->SPg=0,Aptr=Aptr->Next);
  if (Bl) l=ACvar[Bl]+1;
  else if (flagH)  l=Mptr->n1;
  else l=0;
  for (ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
    if (ACptr->Area!=NULL) BEptr=ACptr->Area->Slack;
    i=ACvar[ACptr->N];

  /* ------------------- Power Mismatches ------------------- */
    Vi=ACptr->V;
    di=ACptr->Ang;
    if (flagP) {Pl=ACptr->Pl; Ql=ACptr->Ql;}
    else {
      Pl=(ACptr->Pn+lambda*ACptr->Pnl)*pow(Vi,ACptr->a)+
         (ACptr->Pz+lambda*ACptr->Pzl)*Vi*Vi;
      Ql=(ACptr->Qn+lambda*ACptr->Qnl)*pow(Vi,ACptr->b)+
         (ACptr->Qz+lambda*ACptr->Qzl)*Vi*Vi;
    }
    Pg=ACptr->PG;
    DPg=ACptr->DPG;
    if (Acont && Narea>1) ACptr->Area->Slack->Area->SPg+=DPg;
    else SPg+=DPg;
    Pi=Pg-Pl-Vi*Vi*ACptr->G;
    Qg=ACptr->Qg;
    Qi=Qg-Ql+Vi*Vi*ACptr->B;
    dPiid=dQiid=0;
    dPiiv=dQiiv=0;
    SPij=SQij=0;
    for(ELptr=ACptr->Elem; ELptr!=NULL; ELptr=ELptr->Next) {
      Eptr=ELptr->Eptr;
      if (Eptr->From==ACptr) To=Eptr->To;
      else To=Eptr->From;
      j=ACvar[To->N];
      Vj=To->V;
      dj=To->Ang;
      if(Eptr->From==ACptr) {
        gij=(Eptr->G*cos(Eptr->Ang)-Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
        bij=(Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
        gsij=(Eptr->G1+Eptr->G)*Eptr->Tap*Eptr->Tap-gij;
        bsij=(Eptr->B1+Eptr->B)*Eptr->Tap*Eptr->Tap-bij;
      } else {
        gij=(Eptr->G*cos(Eptr->Ang)+Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
        bij=(-Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
        gsij=Eptr->G+Eptr->G2-gij;
        bsij=Eptr->B+Eptr->B2-bij;
      }
      dPijd=Vi*Vj*(gij*sin(di-dj)-bij*cos(di-dj));
      dPiid=dPiid-dPijd;
      dPijv=Vi*(gij*cos(di-dj)+bij*sin(di-dj));
      dPiiv=dPiiv+Vj/Vi*dPijv-2*Vi*(gij+gsij);
      dQijd= -Vj*dPijv;
      dQiid=dQiid-dQijd;
      dQijv=dPijd/Vj;
      dQiiv=dQiiv+Vj/Vi*dQijv+2*Vi*(bij+bsij);
      SPij=SPij+Vi*Vi*(gij+gsij)+dQijd;
      SQij=SQij-Vi*Vi*(bij+bsij)-dPijd;
      if(flagJ){
        /* df/ddelta_j */
        if (!strpbrk(To->Type,"S")) {
          JacElement(Mptr,i,j,dPijd);
          JacElement(Mptr,i+1,j,dQijd);
        }
        /* df/dV_j */
        if (To->Cont!=NULL) {
          JacElement(Mptr,i,j+1,dPijv);
          JacElement(Mptr,i+1,j+1,dQijv);
        } else if (flagH && strpbrk(To->Type,"L")) {
          JacElement(Mptr,i,Mptr->n1,dPijv);
          JacElement(Mptr,i+1,Mptr->n1,dQijv);
        } else {
          JacElement(Mptr,i,j+1,0.);
          JacElement(Mptr,i+1,j+1,0.);
        }
        /* df/dControl_ij */
        if (PQcont && strpbrk(Eptr->Type,"PQMN")){
          if (Acont && strpbrk(Eptr->Cont->Type,"A")) k=1; else k=0;
          j=ACvar[Eptr->Cont->N]+1+k+Eptr->Cont->Ncont-Eptr->Ncont;
          if(!strcmp(Eptr->Type,"RQ")) {
            val1= -dQijd/Eptr->Tap;
            val2=dPijd/Eptr->Tap;
            if (Eptr->From==ACptr){
              val1=val1-2*Vi*Vi*(Eptr->G1+Eptr->G)*Eptr->Tap;
              val2=val2+2*Vi*Vi*(Eptr->B1+Eptr->B)*Eptr->Tap;
            }
          }
          else if(!strcmp(Eptr->Type,"RP")) {
            val1=dPijd;
            val2=dQijd;
            if (Eptr->From!=ACptr) {
              val1= -val1;
              val2= -val2;
            }
          } else val1=val2=0;
          JacElement(Mptr,i,j,val1);
          JacElement(Mptr,i+1,j,val2);
        }
        else if(Rcont && Eptr->Cont!=NULL) {
          j=ACvar[Eptr->Cont->N]+1;
          if(!strcmp(Eptr->Type,"R")) {
            val1= -dQijd/Eptr->Tap;
            val2=dPijd/Eptr->Tap;
            if (Eptr->From==ACptr){
              val1=val1-2*Vi*Vi*(Eptr->G1+Eptr->G)*Eptr->Tap;
              val2=val2+2*Vi*Vi*(Eptr->B1+Eptr->B)*Eptr->Tap;
            }
          } else val1=val2=0;
          JacElement(Mptr,i,j,val1);
          JacElement(Mptr,i+1,j,val2);
        }
      }
    }
    if (flagF) {
      dF[i]=Pi-SPij;
      dF[i+1]=Qi-SQij;
    }
    if (flagJ) {
      dPiiv=dPiiv-2*Vi*ACptr->G;
      dQiiv=dQiiv+2*Vi*ACptr->B;
      if (!flagP) {
        dPiiv=dPiiv-2*(ACptr->Pz+ACptr->Pzl*lambda)*Vi-
              ACptr->a*(ACptr->Pn+ACptr->Pnl*lambda)*pow(Vi,ACptr->a-1.0);
        dQiiv=dQiiv-2*(ACptr->Qz+ACptr->Qzl*lambda)*Vi-
              ACptr->b*(ACptr->Qn+ACptr->Qnl*lambda)*pow(Vi,ACptr->b-1.0);
      }
      /* df/dKg */
      j=ACvar[BEptr->N];
      if(DPg) {
        if (strpbrk(BEptr->Type,"S")) { JacElement(Mptr,i,j,DPg);}
        else if(Acont) { JacElement(Mptr,i,j+2,DPg);}
      }
      /* df/dlambda */
      if (l) {
        if (ACptr->Pzl || ACptr->Pnl) JacElement(Mptr,i,l,-ACptr->Pzl*Vi*Vi-ACptr->Pnl*pow(Vi,ACptr->a));
        if (ACptr->Qzl || ACptr->Qnl) JacElement(Mptr,i+1,l,-ACptr->Qzl*Vi*Vi-ACptr->Qnl*pow(Vi,ACptr->b));
      }
      /* df/ddelta_i */
      if (!strpbrk(ACptr->Type,"S")) {
        JacElement(Mptr,i,i,dPiid);
        JacElement(Mptr,i+1,i,dQiid);
      }
      /* df/dV_i */
      if (ACptr->Cont!=NULL) {
        JacElement(Mptr,i,i+1,dPiiv);
        JacElement(Mptr,i+1,i+1,dQiiv);
      } else if (flagH && strpbrk(ACptr->Type,"L")) {
        JacElement(Mptr,i,Mptr->n1,dPiiv);
        JacElement(Mptr,i+1,Mptr->n1,dQiiv);
      } else {
        JacElement(Mptr,i,i+1,0.);
        JacElement(Mptr,i+1,i+1,0.);
      }
      /*  df/dQg   */
      if(strpbrk(ACptr->Type,"V,Q,S,G,Z")) {
        if(QRcont && strpbrk(ACptr->Type,"G")) {
          j=ACvar[ACptr->Cont->N];
          val1=ACptr->Kbg;
        } else {
          j=i;
          val1=1.;
        }
        if (strpbrk(ACptr->cont,"V")) JacElement(Mptr,i+1,j+1,val1);
        else                          JacElement(Mptr,i+1,j+1,0.);
        if (ACptr->Gen!=NULL) {  /*Gen. Model BEGIN */
          j=ACptr->Gen->Nvar;
          if (strpbrk(ACptr->cont,"I")) JacElement(Mptr,i+1,j+11,1.);
          else                          JacElement(Mptr,i+1,j+11,0.);
          if (strpbrk(ACptr->cont,"E")) JacElement(Mptr,i+1,j+1,1.);
          else                          JacElement(Mptr,i+1,j+1,0.);
        } /* Gen. Model  END */
      }
    }

  /* ---------------------- Area Control --------------------- */
    if (Acont && strpbrk(ACptr->Type,"A")){
      i=i+2;
      SPij=0;
      for (ELptr=BEptr->Area->Elem; ELptr!=NULL; ELptr=ELptr->Next){
        Eptr=ELptr->Eptr;
        From=Eptr->From;
        To=Eptr->To;
        if (Eptr->From!=Eptr->Meter) {
          From=Eptr->To;
          To=Eptr->From;
        }
        val2=1.0;
        if (Eptr->Meter->Area!=ACptr->Area) val2= -1.0;
        Vi=From->V;
        di=From->Ang;
        Vj=To->V;
        dj=To->Ang;
        if(Eptr->From==Eptr->Meter) {
          gij=(Eptr->G*cos(Eptr->Ang)-Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
          bij=(Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
          gsij=(Eptr->G1+Eptr->G)*Eptr->Tap*Eptr->Tap-gij;
        } else {
          gij=(Eptr->G*cos(Eptr->Ang)+Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
          bij=(-Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
          gsij=Eptr->G+Eptr->G2-gij;
        }
        dPiid= val2*(-Vi*Vj*(gij*sin(di-dj)-bij*cos(di-dj)));
        dPiiv= val2*(-2*Vi*(gij+gsij)+Vj*(gij*cos(di-dj)+bij*sin(di-dj)));
        dPijd= -dPiid;
        dPijv=val2*(Vi*(gij*cos(di-dj)+bij*sin(di-dj)));
        SPij=SPij+val2*Vi*Vi*(gij+gsij)-Vj*dPijv;
        if(flagJ){
          j=ACvar[From->N];
          if (!strpbrk(From->Type,"S")) {  JacElement(Mptr,i,j,dPiid);}
          if (From->Cont!=NULL) JacElement(Mptr,i,j+1,dPiiv);
          else if (flagH && strpbrk(From->Type,"L")) JacElement(Mptr,i,Mptr->n1,dPiiv);
          else JacElement(Mptr,i,j+1,0.);
          j=ACvar[To->N];
          if (!strpbrk(To->Type,"S")) { JacElement(Mptr,i,j,dPijd);}
          if (To->Cont!=NULL) JacElement(Mptr,i,j+1,dPijv);
          else if (flagH && strpbrk(To->Type,"L")) JacElement(Mptr,i,Mptr->n1,dPijv);
          else JacElement(Mptr,i,j+1,0.);
          if (PQcont && strpbrk(Eptr->Type,"PQMN")){
            if (Acont && strpbrk(Eptr->Cont->Type,"A")) k=1; else k=0;
            j=ACvar[Eptr->Cont->N]+1+k+Eptr->Cont->Ncont-Eptr->Ncont;
            if (!strcmp(Eptr->Type,"RQ")){
              val1=Vj*dPijv/Eptr->Tap;
              if (Eptr->From==Eptr->Meter) val1=val1-2*Vi*Vi*(Eptr->G1+Eptr->G)*Eptr->Tap;
            }
            else if (!strcmp(Eptr->Type,"RP")){
              val1=dPiid;
              if (Eptr->From==Eptr->Meter) val1= -val1;
            } else val1=0;
            JacElement(Mptr,i,j,val1);
          }
          else if(Rcont && Eptr->Cont!=NULL) {
            j=ACvar[Eptr->Cont->N]+1;
            if (!strcmp(Eptr->Type,"R")){
              val1=Vj*dPijv/Eptr->Tap;
              if (Eptr->From==Eptr->Meter) val1=val1-2*Vi*Vi*(Eptr->G1+Eptr->G)*Eptr->Tap;
            } else val1=0;
            JacElement(Mptr,i,j,val1);
          }
        }
      }
      if (flagF) dF[i]=ACptr->Area->P-SPij;
    }

  /* -------------- Regulating Transf. ----------------------- */
    if(PQcont) for (ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next) {
      Eptr=ELptr->Eptr;
      if(strpbrk(Eptr->Type,"PQMN")){
        if (Acont && strpbrk(Eptr->Cont->Type,"A")) k=1; else k=0;
        i=ACvar[ACptr->N]+1+k+ACptr->Ncont-Eptr->Ncont;
        if(Eptr->From==ACptr) {
          From=Eptr->From;
          To=Eptr->To;
          gij=(Eptr->G*cos(Eptr->Ang)-Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
          bij=(Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
          gsij=(Eptr->G1+Eptr->G)*Eptr->Tap*Eptr->Tap-gij;
          bsij=(Eptr->B1+Eptr->B)*Eptr->Tap*Eptr->Tap-bij;
        }
        else {
          From=Eptr->To;
          To=Eptr->From;
          gij=(Eptr->G*cos(Eptr->Ang)+Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
          bij=(-Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
          gsij=Eptr->G+Eptr->G2-gij;
          bsij=Eptr->B+Eptr->B2-bij;
        }
        Vi=From->V;
        di=From->Ang;
        Vj=To->V;
        dj=To->Ang;
        dPiid= -Vi*Vj*(gij*sin(di-dj)-bij*cos(di-dj));
        dPiiv= -2*Vi*(gij+gsij)+Vj*(gij*cos(di-dj)+bij*sin(di-dj));
        dPijd= -dPiid;
        dPijv=Vi*(gij*cos(di-dj)+bij*sin(di-dj));
        SPij=Vi*Vi*(gij+gsij)-Vj*dPijv;
        dQiid=Vj*dPijv;
        dQijd= -dQiid;
        dQiiv= -dPiid/Vi+2*Vi*(bij+bsij);
        dQijv= -dPiid/Vj;
        SQij= -Vi*Vi*(bij+bsij)+dPiid;
        if (strpbrk(Eptr->Type,"PM")){
          if(flagJ){
            if(!strcmp(Eptr->Type,"RP")) {
              val1= -dPiid;
              if (Eptr->From!=ACptr) val1= -val1;
            } else val1=1;
            JacElement(Mptr,i,i,val1);
            j=ACvar[From->N];
            if (!strpbrk(From->Type,"S")) { JacElement(Mptr,i,j,dPiid);}
            if (From->Cont!=NULL) JacElement(Mptr,i,j+1,dPiiv);
            else if (flagH && strpbrk(From->Type,"L")) JacElement(Mptr,i,Mptr->n1,dPiiv);
            else JacElement(Mptr,i,j+1,0.);
            j=ACvar[To->N];
            if (!strpbrk(To->Type,"S")) { JacElement(Mptr,i,j,dPijd);}
            if (To->Cont!=NULL) JacElement(Mptr,i,j+1,dPijv);
            else if (flagH && strpbrk(To->Type,"L")) JacElement(Mptr,i,Mptr->n1,dPijv);
            else JacElement(Mptr,i,j+1,0.);
          }
          if (flagF) dF[i]=Eptr->Cvar-SPij;
        }
        else {
          if(flagJ){
            if(!strcmp(Eptr->Type,"RQ")) {
              val1= -dPiid/Eptr->Tap;
              if (Eptr->From==ACptr) val1=val1+2*Vi*Vi*(Eptr->B1+Eptr->B)*Eptr->Tap;
            } else val1=1;
            JacElement(Mptr,i,i,val1);
            j=ACvar[From->N];
            if (!strpbrk(From->Type,"S")) { JacElement(Mptr,i,j,dQiid);}
            if (From->Cont!=NULL) JacElement(Mptr,i,j+1,dQiiv);
            else if (flagH && strpbrk(From->Type,"L")) JacElement(Mptr,i,Mptr->n1,dQiiv);
            else JacElement(Mptr,i,j+1,0.);
            j=ACvar[To->N];
            if (!strpbrk(To->Type,"S"))  JacElement(Mptr,i,j,dQijd);
            if (To->Cont!=NULL) JacElement(Mptr,i,j+1,dQijv);
            else if (flagH && strpbrk(To->Type,"L")) JacElement(Mptr,i,Mptr->n1,dQijv);
            else JacElement(Mptr,i,j+1,0.);
          }
          if (flagF) dF[i]=Eptr->Cvar-SQij;
        }
      }
    }

  /* -------------- Generator Model ----------------------- */
    if (ACptr->Gen!=NULL) {
      i=ACptr->Gen->Nvar;
      Ra=ACptr->Gen->Ra;
      Xd=ACptr->Gen->Xd;
      Xq=ACptr->Gen->Xq;
      Eq=ACptr->Gen->Eq;
      dg=ACptr->Gen->dg;
      Vr=ACptr->Gen->Vr;
      Vim=ACptr->Gen->Vi;
      Ir=ACptr->Gen->Ir;
      Iim=ACptr->Gen->Ii;
      Vq=ACptr->Gen->Vq;
      Vd=ACptr->Gen->Vd;
      Iq=ACptr->Gen->Iq;
      Id=ACptr->Gen->Id;
      Ia=ACptr->Gen->Ia;
      if (flagF) {
        dF[i+1]=Pg-Vr*Ir-Vim*Iim;
        dF[i+2]=Qg-Vim*Ir+Vr*Iim;
        dF[i+3]=Eq-Vq-Ra*Iq+Xd*Id;
        dF[i+4]=Vd+Ra*Id+Xq*Iq;
        dF[i+5]=Vr-cos(dg)*Vq+sin(dg)*Vd;
        dF[i+6]=Vim-sin(dg)*Vq-cos(dg)*Vd;
        dF[i+7]=Ir-cos(dg)*Iq+sin(dg)*Id;
        dF[i+8]=Iim-sin(dg)*Iq-cos(dg)*Id;
        dF[i+9]=Vr-Vi*cos(di);
        dF[i+10]=Vim-Vi*sin(di);
        dF[i+11]=Ia*Ia-Ir*Ir-Iim*Iim;
      }
      if (flagJ) {
        /* df1/dKg */
        j=ACvar[BEptr->N];
        if(DPg) {
          if (strpbrk(BEptr->Type,"S")) JacElement(Mptr,i+1,j,DPg);
          else if(Acont) JacElement(Mptr,i+1,j+2,DPg);
        }
        /* df1/dVr, df1/dVi, df1/dIr, df1/dIi */
        JacElement(Mptr,i+1,i+3,-Ir);
        JacElement(Mptr,i+1,i+4,-Iim);
        JacElement(Mptr,i+1,i+5,-Vr);
        JacElement(Mptr,i+1,i+6,-Vim);
        /* df2/dQg */
        if (QRcont && strpbrk(ACptr->Type,"G")) {
          j=ACvar[ACptr->Cont->N];
          if (strpbrk(ACptr->cont,"V")) JacElement(Mptr,i+2,j+1,ACptr->Kbg);
          else                          JacElement(Mptr,i+2,j+1,0.);
        }
        else {
          j=ACvar[ACptr->N];
          if (strpbrk(ACptr->cont,"V")) JacElement(Mptr,i+2,j+1,1.);
          else                          JacElement(Mptr,i+2,j+1,0.);
        }
        if (strpbrk(ACptr->cont,"I")) JacElement(Mptr,i+2,i+11,1.);
        else                          JacElement(Mptr,i+2,i+11,0.);
        if (strpbrk(ACptr->cont,"E")) JacElement(Mptr,i+2,i+1,1.);
        else                          JacElement(Mptr,i+2,i+1,0.);
        /* df2/dVr, df2/dVi, df2/dIr, df2/dIi */
        JacElement(Mptr,i+2,i+3,Iim);
        JacElement(Mptr,i+2,i+4,-Ir);
        JacElement(Mptr,i+2,i+5,-Vim);
        JacElement(Mptr,i+2,i+6,Vr);
        /* df3/dx  */
        if (strpbrk(ACptr->cont,"E")) JacElement(Mptr,i+3,i+1,0.);
        else                          JacElement(Mptr,i+3,i+1,1.);
        JacElement(Mptr,i+3,i+7,-1.);
        JacElement(Mptr,i+3,i+9,-Ra);
        JacElement(Mptr,i+3,i+10,Xd);
        /* df4/dx  */
        JacElement(Mptr,i+4,i+8,1.);
        JacElement(Mptr,i+4,i+9,Xq);
        JacElement(Mptr,i+4,i+10,Ra);
        /* df5/dx  */
        JacElement(Mptr,i+5,i+2,sin(dg)*Vq+cos(dg)*Vd);
        JacElement(Mptr,i+5,i+3,1.);
        JacElement(Mptr,i+5,i+7,-cos(dg));
        JacElement(Mptr,i+5,i+8,sin(dg));
        /* df6/dx  */
        JacElement(Mptr,i+6,i+2,-cos(dg)*Vq+sin(dg)*Vd);
        JacElement(Mptr,i+6,i+4,1.);
        JacElement(Mptr,i+6,i+7,-sin(dg));
        JacElement(Mptr,i+6,i+8,-cos(dg));
        /* df7/dx  */
        JacElement(Mptr,i+7,i+2,sin(dg)*Iq+cos(dg)*Id);
        JacElement(Mptr,i+7,i+5,1.);
        JacElement(Mptr,i+7,i+9,-cos(dg));
        JacElement(Mptr,i+7,i+10,sin(dg));
        /* df8/dx  */
        JacElement(Mptr,i+8,i+2,-cos(dg)*Iq+sin(dg)*Id);
        JacElement(Mptr,i+8,i+6,1.);
        JacElement(Mptr,i+8,i+9,-sin(dg));
        JacElement(Mptr,i+8,i+10,-cos(dg));
        /* df9/dV  */
        j=ACvar[ACptr->N];
        if (ACptr->Cont!=NULL) JacElement(Mptr,i+9,j+1,-cos(di));
        else                   JacElement(Mptr,i+9,j+1,0.);
        /* df9/ddelta  */
        if (!strpbrk(ACptr->Type,"S")) {
          j=ACvar[ACptr->N];
          JacElement(Mptr,i+9,j,Vi*sin(di));
        }
        /* df9/dVr  */
        JacElement(Mptr,i+9,i+3,1.);
        /* df10/dV  */
        j=ACvar[ACptr->N];
        if (ACptr->Cont!=NULL) JacElement(Mptr,i+10,j+1,-sin(di));
        else                   JacElement(Mptr,i+10,j+1,0.);
        /* df10/ddelta  */
        if (!strpbrk(ACptr->Type,"S")) {
          j=ACvar[ACptr->N];
          JacElement(Mptr,i+10,j,-Vi*cos(di));
        }
        /* df10/dVi  */
        JacElement(Mptr,i+10,i+4,1.);
        /* df11/dx  */
        JacElement(Mptr,i+11,i+5,-2.*Ir);
        JacElement(Mptr,i+11,i+6,-2.*Iim);
        if (strpbrk(ACptr->cont,"I")) JacElement(Mptr,i+11,i+11,0.);
        else                          JacElement(Mptr,i+11,i+11,2.*Ia);
      }
    }
  }

  /* -------------- Detect Area/System Generation Errors ----------------------- */
  if (!flagPgMax && val!=NULL) {
    if (Acont && Narea>1) for(Aptr=dataPtr->Area;Aptr!=NULL;Aptr=Aptr->Next) {
      if (!Aptr->SPg) {
        if (!*val) *val=-1;
        else *val=-2;
        return(Aptr);
      }
    } else  if (!SPg) {
      if (!*val) *val=-1;
      else *val=-2;
    }
  }
  return(NULL);
}


