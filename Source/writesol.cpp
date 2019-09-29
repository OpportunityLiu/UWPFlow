#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

/* Write solution in files: -Printed output.
                            -Bus voltages and angles.
                            -IEEE common format for SSSP or other
                             P.F. programs.   */

#include "write.h"
#include <string.h> /* FACTS */

/* --------------- Output ----------------- */
#ifdef ANSIPROTO
void Output(INDEX Iter,char *File1,char *str)
#else
void Output(Iter,File1,str)
INDEX Iter;
char *File1,*str;
#endif
{
  ACbusData *ACptr,*From,*To;
  AClist *ACLptr;
  DClist *DCLptr;
  DCbusData *DCptrR,*DCptrI,*DCptr,*DCptrp;
  ElementData *Eptr;
  ElementList *ELptr;
  AreaData *Aptr;
  VALUETYPE P,Q,Pl,Ql,Pij,Qij,Pji,Qji,Vi,di,Vj,dj,Iij,Iji,I,ratio;
  VALUETYPE G,B,Gp,Bp,Gi,Bi,Gj,Bj,Pg,KVi,KVj,KV,Qm,delta,theta,vals;
  FILE *OutFile;
  INDEX i,j;
  SVCbusData *SVCptr;                             /* FACTS */
  TCSCbusData *TCSCptr;                           /* FACTS */
  STATCOMbusData *STATCOMptr;                     /* FACTS */
  VALUETYPE Xc,Max,Vn;                            /* FACTS */

  OutFile=OpenOutput(File1);
  fCustomPrint(OutFile,"\n%s\n",str);
  i=0;
  while(i<=2 && dataPtr->Title[0][0]!='\0'){
    fCustomPrint(OutFile,"%s",dataPtr->Title[i]);
    i++;
  }
  if (!flagH) lambda_o=0;
  fCustomPrint(OutFile,"      Loading factor -> %-10.6lg\n",lambda+lambda_o);
  fCustomPrint(OutFile,"            AC buses -> %d\n",Nac);
  fCustomPrint(OutFile,"            PV buses -> %d\n",Nvolt);
  fCustomPrint(OutFile,"            X buses  -> %d\n",NXvolt);
  fCustomPrint(OutFile,"            Z buses  -> %d\n",NZvolt);
  fCustomPrint(OutFile,"            AC elem. -> %d\n",NacEl);
  fCustomPrint(OutFile,"         V Reg. Trf. -> %d\n",NregV);
  fCustomPrint(OutFile,"        PQ Reg. Trf. -> %d\n",NregPQ);
  fCustomPrint(OutFile,"            DC buses -> %d\n",Ndc);
  fCustomPrint(OutFile,"            DC lines -> %d\n",Ndc/2);
  fCustomPrint(OutFile,"                SVCs -> %d\n",Nsvc);     /* FACTS */
  fCustomPrint(OutFile,"               TCSCs -> %d\n",Ntcsc);    /* FACTS */
  fCustomPrint(OutFile,"            STATCOMs -> %d\n",Nstatcom); /* FACTS */
  fCustomPrint(OutFile,"           No. Areas -> %d\n",Narea);
  fCustomPrint(OutFile,"           Iterations -> %d (Maximum = %d)\n",Iter,MaxIter);
  fCustomPrint(OutFile,"   Max. p.u. mismatch -> %-8.4lg (Tolerance = %-8.4lg)\n",MaxdFi,Tol);
  fCustomPrint(OutFile,"    Reference Bus(es) -> ");
  i=0;
  for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next)
    if(strpbrk(ACptr->Type,"S")){
      if (i>0)   fCustomPrint(OutFile,"                         ");
      fCustomPrint(OutFile,"%d %s (Angle=%6.2lf deg.)\n",ACptr->Num,ACptr->Name,ACptr->Ang/K3);
      i++;
    }
  if (i==0) fCustomPrint(OutFile,"\n");
  fCustomPrint(OutFile,"\n");

  if (Narea<2) ACptr=dataPtr->ACbus;
  else {
    Aptr=dataPtr->Area;
    ACLptr=Aptr->AC;
    ACptr=ACLptr->AC;
  }
  fCustomPrint(OutFile,"\n                                      ***** AC RESULTS *****\n");
  fCustomPrint(OutFile,"                   L=lower limit    H=higher limit      O=over limit    U=under limit\n");
  fCustomPrint(OutFile,"--|----|------------|------|-------|--------|--------|--------|----|------------|-|--------|--------|--------|-------|-|-----------------\n");
  fCustomPrint(OutFile," A    i Bus           V(pu)   V(kV)   Pg(MW)    Pload   Pshunt|   j Bus          C      Pij  Plosses |Iij|(A) kVi/kVj T  Controlled Bus  \n");
  fCustomPrint(OutFile," n      Name         d(deg)  d(rad) Qg(MVAR)    Qload   Qshunt|     Name         r      Qij  Qlosses           a(deg)      k Name        \n");
  while(ACptr!=NULL){
    fCustomPrint(OutFile,"--|----|------------|------|-------|--------|--------|--------|----|------------|-|--------|--------|--------|-------|-|----|------------\n");
    if(ACptr->Area!=NULL) {
      fCustomPrint(OutFile,"%2d ",ACptr->Area->N);
    } else fCustomPrint(OutFile,"%2d ",0);
    fCustomPrint(OutFile,"%4d ",ACptr->Num);
    fCustomPrint(OutFile,"%12s ",ACptr->Name);
    fCustomPrint(OutFile,"%6.4lf ",ACptr->V);
    KVi=ACptr->KV;
    if (KVi > 0) fCustomPrint(OutFile,"%7.2lf",KVi*ACptr->V);
    else fCustomPrint(OutFile,"%7s","");
    if(ACptr->Vmax==ACptr->Vmin) {
      if(ACptr->Vlmax==ACptr->Vlmin) fCustomPrint(OutFile," ");
      else if(ACptr->V==ACptr->Vlmin) fCustomPrint(OutFile,"L");
      else if(ACptr->V==ACptr->Vlmax) fCustomPrint(OutFile,"H");
      else if(ACptr->V<ACptr->Vlmin) fCustomPrint(OutFile,"U");
      else if(ACptr->V>ACptr->Vlmax) fCustomPrint(OutFile,"O");
      else fCustomPrint(OutFile," ");
    }
    else if(ACptr->V==ACptr->Vmin) fCustomPrint(OutFile,"L");
    else if(ACptr->V==ACptr->Vmax) fCustomPrint(OutFile,"H");
    else if(ACptr->V<ACptr->Vmin) fCustomPrint(OutFile,"U");
    else if(ACptr->V>ACptr->Vmax) fCustomPrint(OutFile,"O");
    else fCustomPrint(OutFile," ");
    Pg=ACptr->PG;
    fCustomPrint(OutFile,"%8.2lf",Pg*Sn);
    if(Pg==ACptr->Pmax) fCustomPrint(OutFile,"H");
    else if(Pg>ACptr->Pmax) fCustomPrint(OutFile,"O");
    else fCustomPrint(OutFile," ");
    Pl=(ACptr->Pn+lambda*ACptr->Pnl)*pow(ACptr->V,ACptr->a)+
       (ACptr->Pz+lambda*ACptr->Pzl)*ACptr->V*ACptr->V;
    Ql=(ACptr->Qn+lambda*ACptr->Qnl)*pow(ACptr->V,ACptr->b)+
       (ACptr->Qz+lambda*ACptr->Qzl)*ACptr->V*ACptr->V;
    fCustomPrint(OutFile,"%8.2lf ",Pl*Sn);
    fCustomPrint(OutFile,"%8.2lf|",ACptr->G*ACptr->V*ACptr->V*Sn);
    i=0;

  /* -------------------------------------- DC element results ------------------------------------------- */
    for (DCLptr=ACptr->DC;DCLptr!=NULL;DCLptr=DCLptr->Next) {
      if (i!=0) fCustomPrint(OutFile,"%62s|","");
      DCptr=DCLptr->DC;
      fCustomPrint(OutFile,"  DC %-12s %1s ",DCptr->Name,DCptr->Type);
      fCustomPrint(OutFile,"%8.2lf ",-DCptr->P*Sn);
      fCustomPrint(OutFile," Alpha=%6.2lf",DCptr->Alfa/K3);
      if(DCptr->Alfa<=DCptr->AlfaMin) fCustomPrint(OutFile,"L");
      else if(DCptr->Alfa>=DCptr->AlfaMax) fCustomPrint(OutFile,"H");
      else fCustomPrint(OutFile," ");
      fCustomPrint(OutFile," Tap=%6.4lf",DCptr->Tap);
      if(DCptr->Tap==DCptr->TapMin) fCustomPrint(OutFile,"L");
      else if(DCptr->Tap==DCptr->TapMax) fCustomPrint(OutFile,"H");
      else if(DCptr->Tap<DCptr->TapMin) fCustomPrint(OutFile,"U");
      else if(DCptr->Tap>DCptr->TapMax) fCustomPrint(OutFile,"O");
      else fCustomPrint(OutFile," ");
      fCustomPrint(OutFile,"\n");
      if (i==0) {
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
        i++; fCustomPrint(OutFile,"%21s","");
        delta=ACptr->Ang;
        if (delta>=0) vals=1.00;
        else          vals=-1.00;
        if (fabs(delta)>2*PI) delta=delta-vals*floor(fabs(delta)/(2*PI))*2*PI;
        if (fabs(delta)>PI) delta=delta-vals*2*PI;
        ACptr->Ang=delta;
        fCustomPrint(OutFile,"%6.2lf ",ACptr->Ang/K3);
        fCustomPrint(OutFile,"%7.4lf ",ACptr->Ang);
        fCustomPrint(OutFile,"%8.2lf",ACptr->Qg*Sn);
        if(ACptr->Max==ACptr->Min) fCustomPrint(OutFile," ");
        else if(ACptr->Qg==ACptr->Min) fCustomPrint(OutFile,"L");
        else if(ACptr->Qg==ACptr->Max) fCustomPrint(OutFile,"H");
        else if(ACptr->Qg<ACptr->Min) fCustomPrint(OutFile,"U");
        else if(ACptr->Qg>ACptr->Max) fCustomPrint(OutFile,"O");
        else fCustomPrint(OutFile," ");
        fCustomPrint(OutFile,"%8.2lf ",Ql*Sn);
        fCustomPrint(OutFile,"%8.2lf|",ACptr->B*ACptr->V*ACptr->V*Sn);
      } else fCustomPrint(OutFile,"%62s|","");
      fCustomPrint(OutFile,"%19s ","");
      fCustomPrint(OutFile,"%8.2lf ",-DCptr->Q*Sn);
      fCustomPrint(OutFile," Gamma=%6.2lf",DCptr->Gamma/K3);
      if(DCptr->Gamma<=DCptr->GammaMin) fCustomPrint(OutFile,"L");
      else fCustomPrint(OutFile," ");
      fCustomPrint(OutFile,"\n");
    }

  /* --------------------------------------- AC element results -------------------------------------------- */
    for (ELptr=ACptr->Elem;ELptr!=NULL;ELptr=ELptr->Next) { Eptr=ELptr->Eptr;
      Vi=Eptr->From->V;  di=Eptr->From->Ang;  KVi=Eptr->From->KV;
      Vj=Eptr->To->V;    dj=Eptr->To->Ang;    KVj=Eptr->To->KV;
      if (KVi>0 && KVj>0) ratio=KVi/KVj/Eptr->Tap;
      else                ratio=1/Eptr->Tap;
      G=(Eptr->G*cos(Eptr->Ang)-Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
      B=(Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
      Gi=(Eptr->G1+Eptr->G)*pow(Eptr->Tap,2.0)-G;
      Bi=(Eptr->B1+Eptr->B)*pow(Eptr->Tap,2.0)-B;
      Gp=(Eptr->G*cos(Eptr->Ang)+Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
      Bp=(-Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
      Gj=Eptr->G+Eptr->G2-Gp;
      Bj=Eptr->B+Eptr->B2-Bp;
      Pij=Vi*Vi*(Gi+G)-Vi*Vj*(G*cos(di-dj)+B*sin(di-dj));
      Qij= -Vi*Vi*(Bi+B)-Vi*Vj*(G*sin(di-dj)-B*cos(di-dj));
      Iij=sqrt(Pij*Pij+Qij*Qij)/Vi;
      Pji=Vj*Vj*(Gj+Gp)-Vi*Vj*(Gp*cos(dj-di)+Bp*sin(dj-di));
      Qji= -Vj*Vj*(Bj+Bp)-Vi*Vj*(Gp*sin(dj-di)-Bp*cos(dj-di));
      Iji=sqrt(Pji*Pji+Qji*Qji)/Vj;
      if(Eptr->From==ACptr)  To=Eptr->To;
      else {
        To=Eptr->From;
        P=Pij;   Pij=Pji;  Pji=P;
        Q=Qij;   Qij=Qji;  Qji=Q;
        I=Iij;   Iij=Iji;  Iji=I;
        KV=KVi;  KVi=KVj;  KVj=KV;
        ratio=1/ratio;
      }
      if (i!=0) fCustomPrint(OutFile,"%62s|","");
      fCustomPrint(OutFile,"%4d ",To->Num);
      fCustomPrint(OutFile,"%12s ",To->Name);
      fCustomPrint(OutFile,"%1s ",Eptr->Ckt);
      fCustomPrint(OutFile,"%8.2lf",Pij*Sn);
      if(ACptr!=Eptr->Cont || strcmp(Eptr->Ctype,"P") || Eptr->Max==Eptr->Min) fCustomPrint(OutFile," ");
      else if(Pij==Eptr->Min) fCustomPrint(OutFile,"L");
      else if(Pij==Eptr->Max) fCustomPrint(OutFile,"H");
      else if(Pij<Eptr->Min) fCustomPrint(OutFile,"U");
      else if(Pij>Eptr->Max) fCustomPrint(OutFile,"O");
      else fCustomPrint(OutFile," ");
      fCustomPrint(OutFile,"%8.2lf ",(Pij+Pji)*Sn);
      if (KVi > 0) fCustomPrint(OutFile,"%8.2lf",Iij*1000*Sn/(sqrt(3.0)*KVi));
      else fCustomPrint(OutFile,"%8s","");
      if(Eptr->Imax<=0) fCustomPrint(OutFile," ");
      else if(Iij==Eptr->Imax) fCustomPrint(OutFile,"H");
      else if(Iij>Eptr->Imax) fCustomPrint(OutFile,"O");
      else fCustomPrint(OutFile," ");
      if (strpbrk(Eptr->Type,"TR")) {
        fCustomPrint(OutFile,"%7.4lf",ratio);
        if(!strcmp(Eptr->Ctype,"P")||Eptr->Tmax==Eptr->Tmin) fCustomPrint(OutFile," ");
        else if(Eptr->Tap==1/Eptr->Tmin) fCustomPrint(OutFile,"L");
        else if(Eptr->Tap==1/Eptr->Tmax) fCustomPrint(OutFile,"H");
        else if(Eptr->Tap>1/Eptr->Tmin) fCustomPrint(OutFile,"U");
        else if(Eptr->Tap<1/Eptr->Tmax) fCustomPrint(OutFile,"O");
        else fCustomPrint(OutFile," ");
      }
      else  fCustomPrint(OutFile,"%8s","");
      fCustomPrint(OutFile,"%1s ",Eptr->Ctype);
      if (Eptr->Cont!=NULL) {
        fCustomPrint(OutFile,"%4d ",Eptr->Cont->Num);
        fCustomPrint(OutFile,"%12s ",Eptr->Cont->Name);
      }
      fCustomPrint(OutFile,"\n");
      if (i==0) {
        i++; fCustomPrint(OutFile,"%21s","");
        delta=ACptr->Ang;
        if (delta>=0) vals=1.00;
        else          vals=-1.00;
        if (fabs(delta)>2*PI) delta=delta-vals*floor(fabs(delta)/(2*PI))*2*PI;
        if (fabs(delta)>PI) delta=delta-vals*2*PI;
        ACptr->Ang=delta;
        fCustomPrint(OutFile,"%6.2lf ",ACptr->Ang/K3);
        fCustomPrint(OutFile,"%7.4lf ",ACptr->Ang);
        fCustomPrint(OutFile,"%8.2lf",ACptr->Qg*Sn);
        if(ACptr->Qmax==ACptr->Qmin) fCustomPrint(OutFile," ");
        else if(ACptr->Qg==ACptr->Qmin) fCustomPrint(OutFile,"L");
        else if(ACptr->Qg==ACptr->Qmax) fCustomPrint(OutFile,"H");
        else if(ACptr->Qg<ACptr->Qmin) fCustomPrint(OutFile,"U");
        else if(ACptr->Qg>ACptr->Qmax) fCustomPrint(OutFile,"O");
        else fCustomPrint(OutFile," ");
        fCustomPrint(OutFile,"%8.2lf ",Ql*Sn);
        fCustomPrint(OutFile,"%8.2lf|",ACptr->B*ACptr->V*ACptr->V*Sn);
      } else fCustomPrint(OutFile,"%62s|","");
      fCustomPrint(OutFile,"%19s ","");
      fCustomPrint(OutFile,"%8.2lf",Qij*Sn);
      if(ACptr!=Eptr->Cont || strcmp(Eptr->Ctype,"Q") || Eptr->Max==Eptr->Min)
        fCustomPrint(OutFile," ");
      else if(Qij==Eptr->Min) fCustomPrint(OutFile,"L");
      else if(Qij==Eptr->Max) fCustomPrint(OutFile,"H");
      else if(Qij<Eptr->Min) fCustomPrint(OutFile,"U");
      else if(Qij>Eptr->Max) fCustomPrint(OutFile,"O");
      else fCustomPrint(OutFile," ");
      fCustomPrint(OutFile,"%8.2lf ",(Qij+Qji)*Sn);
      fCustomPrint(OutFile,"%9s","");
      if (strpbrk(Eptr->Type,"TR")) {
        fCustomPrint(OutFile,"%7.3lf",Eptr->Ang/K3);
        if(strpbrk(Eptr->Ctype,"QV")||Eptr->Tmax==Eptr->Tmin) fCustomPrint(OutFile," ");
        else if(Eptr->Ang==Eptr->Tmin) fCustomPrint(OutFile,"L");
        else if(Eptr->Ang==Eptr->Tmax) fCustomPrint(OutFile,"H");
        else if(Eptr->Ang<Eptr->Tmin) fCustomPrint(OutFile,"U");
        else if(Eptr->Ang>Eptr->Tmax) fCustomPrint(OutFile,"O");
        else fCustomPrint(OutFile," ");
      }
      fCustomPrint(OutFile,"\n");
    }
    if (i==0) {
      i++; fCustomPrint(OutFile,"%21s","");
      delta=ACptr->Ang;
      if (delta>=0) vals=1.00;
      else          vals=-1.00;
      if (fabs(delta)>2*PI) delta=delta-vals*floor(fabs(delta)/(2*PI))*2*PI;
      if (fabs(delta)>PI) delta=delta-vals*2*PI;
      ACptr->Ang=delta;
      fCustomPrint(OutFile,"%6.2lf ",ACptr->Ang/K3);
      fCustomPrint(OutFile,"%7.4lf ",ACptr->Ang);
      fCustomPrint(OutFile,"%8.2lf",ACptr->Qg*Sn);
      if(ACptr->Qmax==ACptr->Qmin) fCustomPrint(OutFile," ");
      else if(ACptr->Qg==ACptr->Qmin) fCustomPrint(OutFile,"L");
      else if(ACptr->Qg==ACptr->Qmax) fCustomPrint(OutFile,"H");
      else if(ACptr->Qg<ACptr->Qmin) fCustomPrint(OutFile,"U");
      else if(ACptr->Qg>ACptr->Qmax) fCustomPrint(OutFile,"O");
      else fCustomPrint(OutFile," ");
      fCustomPrint(OutFile,"%8.2lf ",Ql*Sn);
      fCustomPrint(OutFile,"%8.2lf|\n",ACptr->B*ACptr->V*ACptr->V*Sn);
    }

  /* -------------------------------------- Generator results -------------------------------------------- */
    if (ACptr->Gen!=NULL){
      fCustomPrint(OutFile," Gen.-> Ia=%7.4lf[pu]",ACptr->Gen->Ia);
      if (ACptr->Gen->Ia>ACptr->Gen->IaMax) fCustomPrint(OutFile,"O");
      else if (ACptr->Gen->Ia==ACptr->Gen->IaMax) fCustomPrint(OutFile,"H");
      else fCustomPrint(OutFile," ");
      if (ACptr->Gen->Ii==0 && ACptr->Gen->Ir==0) fCustomPrint(OutFile," al=%7.2lf[deg]",0.0);
      else fCustomPrint(OutFile," al=%7.2lf[deg]",atan2(ACptr->Gen->Ii,ACptr->Gen->Ir)/K3);
      fCustomPrint(OutFile,"%23s|\n","");
      fCustomPrint(OutFile,"        Eq=%7.4lf[pu]",ACptr->Gen->Eq);
      if (ACptr->Gen->Eq>ACptr->Gen->EqMax) fCustomPrint(OutFile,"O");
      else if (ACptr->Gen->Eq==ACptr->Gen->EqMax) fCustomPrint(OutFile,"H");
      else if (ACptr->Gen->Eq<ACptr->Gen->EqMin) fCustomPrint(OutFile,"U");
      else if (ACptr->Gen->Eq==ACptr->Gen->EqMin) fCustomPrint(OutFile,"L");
      else fCustomPrint(OutFile," ");
      fCustomPrint(OutFile," dg=%7.2lf[deg]",ACptr->Gen->dg/K3);
      fCustomPrint(OutFile,"%23s|\n","");
    }


    if (Narea>1) {
      ACLptr=ACLptr->Next;
      if(ACLptr!=NULL) ACptr=ACLptr->AC;
      else {
        Aptr=Aptr->Next;
        if (Aptr!=NULL) { ACLptr=Aptr->AC; ACptr=ACLptr->AC;}
        else ACptr=NULL;
      }
    } else ACptr=ACptr->Next;
  }

  if (dataPtr->Area!=NULL) {
    fCustomPrint(OutFile,"\n\n                                ***** AC AREA RESULTS *****\n");
    fCustomPrint(OutFile,"---------------------------------|-----------------|----|------------|----|------------|-|--------|----\n");
    fCustomPrint(OutFile," Area                                Slack Bus     |   i Bus Name        j Bus          C  Pmeter   To \n");
    fCustomPrint(OutFile," n Name                              k Name        |     (Meter)           Name         k   (MW)   Area\n");
    for(Aptr=dataPtr->Area;Aptr!=NULL;Aptr=Aptr->Next){
    fCustomPrint(OutFile,"--|------------------------------|----|------------|----|------------|----|------------|-|--------|----\n");
      Aptr->SPg=0;
      fCustomPrint(OutFile,"%2d ",Aptr->N);
      fCustomPrint(OutFile,"%30s ",Aptr->Name);
      fCustomPrint(OutFile,"%4d ",Aptr->BSptr->Num);
      fCustomPrint(OutFile,"%12s|",Aptr->BSptr->Name);
      i=0;
      for (DCLptr=Aptr->DC;DCLptr!=NULL;DCLptr=DCLptr->Next) {
        DCptr=DCLptr->DC->Meter;
        DCptrp=DCptr->To;
        if (i!=0) fCustomPrint(OutFile,"%51s|",""); else i++;
        fCustomPrint(OutFile,"DC %1s %-12s ",DCptr->Type,DCptr->Name);
        fCustomPrint(OutFile,"   %1s %-12s   ",DCptrp->Type,DCptrp->Name);
        P= -DCptr->P;
        if(DCptr->Area!=Aptr) P= -P;
        Aptr->SPg=Aptr->SPg+P;
        fCustomPrint(OutFile,"%8.2lf ",P*Sn);
        if(DCptr->Area!=Aptr) fCustomPrint(OutFile,"%2d\n",DCptr->Area->N);
        else fCustomPrint(OutFile,"%2d\n",DCptrp->Area->N);
      }
      for(ELptr=Aptr->Elem;ELptr!=NULL;ELptr=ELptr->Next) { 
        Eptr=ELptr->Eptr;
        Vi=Eptr->From->V;  di=Eptr->From->Ang;
        Vj=Eptr->To->V;    dj=Eptr->To->Ang;
        G=(Eptr->G*cos(Eptr->Ang)-Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
        B=(Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
        Gi=(Eptr->G1+Eptr->G)*pow(Eptr->Tap,2.0)-G;
        Bi=(Eptr->B1+Eptr->B)*pow(Eptr->Tap,2.0)-B;
        Gp=(Eptr->G*cos(Eptr->Ang)+Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
        Bp=(-Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
        Gj=Eptr->G+Eptr->G2-Gp;
        Bj=Eptr->B+Eptr->B2-Bp;
        if (Eptr->From==Eptr->Meter) {
          From=Eptr->From;
          To=Eptr->To;
          P=Vi*Vi*(Gi+G)-Vi*Vj*(G*cos(di-dj)+B*sin(di-dj));
        } else {
          From=Eptr->To;
          To=Eptr->From;
          P=Vj*Vj*(Gj+Gp)-Vi*Vj*(Gp*cos(dj-di)+Bp*sin(dj-di));
        }
        if (i!=0) fCustomPrint(OutFile,"%51s|",""); else i++;
        fCustomPrint(OutFile,"%4d ",From->Num);
        fCustomPrint(OutFile,"%12s ",From->Name);
        fCustomPrint(OutFile,"%4d ",To->Num);
        fCustomPrint(OutFile,"%12s ",To->Name);
        fCustomPrint(OutFile,"%1s ",Eptr->Ckt);
        if(Eptr->Meter->Area!=Aptr) P= -P;
        Aptr->SPg=Aptr->SPg+P;
        fCustomPrint(OutFile,"%8.2lf ",P*Sn);
        if(From->Area!=Aptr) fCustomPrint(OutFile,"%2d\n",From->Area->N);
        else fCustomPrint(OutFile,"%2d\n",To->Area->N);
      }
      if(i==0) fCustomPrint(OutFile,"\n");
      else fCustomPrint(OutFile,"%51s| MW Export = %-10.2lf  MW Sched. = %-10.2lf\n","",Aptr->SPg*Sn,Aptr->P*Sn);
    }
  }

  if (dataPtr->DCbus!=NULL) {
    fCustomPrint(OutFile,"\n\n                                  ***** DC RESULTS *****\n");
    fCustomPrint(OutFile,"                   L=lower limit    H=higher limit      O=over limit    U=under limit\n");
    fCustomPrint(OutFile,"---|--|----|------------|--------|-|--|--|--------|--------|--------|--------|--------|--------|---------\n");
    fCustomPrint(OutFile,"  i  A    j AC Bus       DC Bus   C S1 S2    Vd       Id       P        Q      Alpha    Gamma      Tap   \n");
    fCustomPrint(OutFile,"     n      Name         Name     v         (KV)      (A)     (MW)    (MVAR)   (deg)    (deg)      (pu)  \n");
    i=0;
    for(DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
      DCptrI=DCptrR->To;
      if(!strcmp(DCptrR->Type,"R")){
    fCustomPrint(OutFile,"---|--|----|------------|--------|-|--|--|--------|--------|--------|--------|--------|--------|---------\n");
        fCustomPrint(OutFile,"%3d ",++i);
        for (j=1;j<=2;j++) {
          if(j==1) {
            DCptr=DCptrR;
            DCptrp=DCptrI;
          }
          else {
            DCptr=DCptrI;
            DCptrp=DCptrR;
            fCustomPrint(OutFile,"%4s","");
          }
          if(DCptr->Area!=NULL) fCustomPrint(OutFile,"%2d ",DCptr->Area->N);
          else fCustomPrint(OutFile,"%2d ",0);
          fCustomPrint(OutFile,"%4d %12s ",DCptr->AC->Num,DCptr->AC->Name);
          fCustomPrint(OutFile,"%8s ",DCptr->Name);
          fCustomPrint(OutFile,"%1s ",DCptr->Type);
          fCustomPrint(OutFile,"%2s ",DCptr->Cont1);
          fCustomPrint(OutFile,"%2s ",DCptr->Cont2);
          fCustomPrint(OutFile,"%8.2lf ",DCptr->Vd*DCptr->Vn);
          fCustomPrint(OutFile,"%8.2lf ",DCptr->Id*Sn/DCptr->Vn*1000);
          fCustomPrint(OutFile,"%8.2lf ",DCptr->P*Sn);
          fCustomPrint(OutFile,"%8.2lf ",DCptr->Q*Sn);
          fCustomPrint(OutFile,"%8.2lf",DCptr->Alfa/K3);
          if(DCptr->Alfa<=DCptr->AlfaMin) fCustomPrint(OutFile,"L");
          else if(DCptr->Alfa>=DCptr->AlfaMax) fCustomPrint(OutFile,"H");
          else fCustomPrint(OutFile," ");
          fCustomPrint(OutFile,"%8.2lf",DCptr->Gamma/K3);
          if(DCptr->Gamma<=DCptr->GammaMin) fCustomPrint(OutFile,"L");
          else fCustomPrint(OutFile," ");
          fCustomPrint(OutFile,"%8.4lf",DCptr->Tap);
          if(DCptr->Tap==DCptr->TapMin) fCustomPrint(OutFile,"L");
          else if(DCptr->Tap==DCptr->TapMax) fCustomPrint(OutFile,"H");
          else if(DCptr->Tap<DCptr->TapMin) fCustomPrint(OutFile,"U");
          else if(DCptr->Tap>DCptr->TapMax) fCustomPrint(OutFile,"O");
          else fCustomPrint(OutFile," ");
          fCustomPrint(OutFile,"\n");
        }
      }
    }
  }

                                        /* FACTS */
  if(dataPtr->SVCbus!=NULL){
    fCustomPrint(OutFile,"\n\n                        ***** SVC RESULTS *****\n\n");
    fCustomPrint(OutFile,"                         L=lower limit    H=higher limit\n\n");
    fCustomPrint(OutFile,"---|--|----|------------|------|------|----|------------|------|------|--------|--------|--------\n");
    fCustomPrint(OutFile,"  i  A    j Cont. Bus      Vj     dj      k SVC Bus        Vk    dk     Qsvc    Alpha     Bv   \n");
    fCustomPrint(OutFile,"     n      Na    me      (pu)   (deg)      Name          (pu)  (deg)  (MVar)   (deg)    (pu)  \n");
    i=0;
    for(SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
      fCustomPrint(OutFile,"---|--|----|------------|------|------|----|------------|------|------|--------|--------|--------\n");
      fCustomPrint(OutFile,"%3d ",++i);
      if(SVCptr->From->Area!=NULL) fCustomPrint(OutFile,"%2d ",SVCptr->From->Area->N);
      else fCustomPrint(OutFile,"%2d ",0);
      fCustomPrint(OutFile,"%4d %12s ",SVCptr->Ctrl->Num,SVCptr->Ctrl->Name);
      fCustomPrint(OutFile,"%6.4lf %6.2lf ",SVCptr->Ctrl->V,SVCptr->Ctrl->Ang/K3);
      fCustomPrint(OutFile,"%4d %12s ",SVCptr->From->Num,SVCptr->Name);
      fCustomPrint(OutFile,"%6.4lf %6.2lf ",SVCptr->From->V,SVCptr->From->Ang/K3);
      fCustomPrint(OutFile,"%8.2lf ",SVCptr->Qsvc*Sn);
      fCustomPrint(OutFile,"%8.2lf",SVCptr->alpha_svc/K3);
      if(SVCptr->alpha_svc<=SVCptr->AlphaMin)      fCustomPrint(OutFile,"L");
      else if(SVCptr->alpha_svc>=SVCptr->AlphaMax) fCustomPrint(OutFile,"H");
      else                                         fCustomPrint(OutFile," ");
      fCustomPrint(OutFile,"%8.4lf ",SVCptr->Bv);
      fCustomPrint(OutFile,"\n\n");
    }
  }

  if(dataPtr->TCSCbus!=NULL){
    fCustomPrint(OutFile,"\n\n                                             ***** TCSC RESULTS *****\n\n");
    fCustomPrint(OutFile,"                                              L=lower limit    H=higher limit\n\n");
    fCustomPrint(OutFile,"---|--|----|------------|------|------|----|------------|------|------|--------|--------|--------|--------|--------|--------|--------|----|--------\n");
    fCustomPrint(OutFile,"  i  A    k From Bus       Vk     dk      m To Bus         Vm     dm    Ptcsc    Qtcsck   Qtcscm     Be     Alpha    Itcsc    delta_t Ctrl   Ctrl  \n");
    fCustomPrint(OutFile,"     n      Name          (pu)   (deg)      Name          (pu)   (deg)  (MW)     (MVar)   (MVar)    (pu)    (deg)    (kA)      (deg)  mode    val  \n");
    i=0;
    for(TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;TCSCptr=TCSCptr->Next){
      fCustomPrint(OutFile,"---|--|----|------------|------|------|----|------------|------|------|--------|--------|--------|--------|--------|--------|--------|----|--------\n");
      Vn=TCSCptr->From->KV;
      Xc=TCSCptr->Xc;
      Max=TCSCptr->Max;
      fCustomPrint(OutFile,"%3d ",++i);
      if(TCSCptr->From->Area!=NULL) fCustomPrint(OutFile,"%2d ",TCSCptr->From->Area->N);
      else fCustomPrint(OutFile,"%2d ",0);
      fCustomPrint(OutFile,"%4d %12s ",TCSCptr->From->Num,TCSCptr->From->Name);
      fCustomPrint(OutFile,"%6.4lf %6.2lf ",TCSCptr->From->V,TCSCptr->From->Ang/K3);
      fCustomPrint(OutFile,"%4d %12s ",TCSCptr->To->Num,TCSCptr->To->Name);
      fCustomPrint(OutFile,"%6.4lf %6.2lf ",TCSCptr->To->V,TCSCptr->To->Ang/K3);
      fCustomPrint(OutFile,"%8.2lf ",TCSCptr->Ptcsc*Sn);
      fCustomPrint(OutFile,"%8.2lf ",TCSCptr->Qtcsck*Sn);
      fCustomPrint(OutFile,"%8.2lf ",TCSCptr->Qtcscm*Sn);
      fCustomPrint(OutFile,"%8.2lf ",TCSCptr->Be);
      fCustomPrint(OutFile,"%8.2lf",TCSCptr->alpha_tcsc/K3);
      if(TCSCptr->alpha_tcsc<=TCSCptr->AlphaMin)      fCustomPrint(OutFile,"L");
      else if(TCSCptr->alpha_tcsc>=TCSCptr->AlphaMax) fCustomPrint(OutFile,"H");
      else                                            fCustomPrint(OutFile," ");
      fCustomPrint(OutFile,"%8.3lf ",TCSCptr->Itcsc*Sn/(sqrt(3.0)*Vn));
      fCustomPrint(OutFile,"%8.4lf ",TCSCptr->delta_t/K3);
      fCustomPrint(OutFile," %c  ",TCSCptr->Cont[0]);
      if (!strcmp(TCSCptr->Cont,"X"))      fCustomPrint(OutFile," %8.3lf",100.0/(TCSCptr->Bset*Max*Xc));
      else if (!strcmp(TCSCptr->Cont,"P")) fCustomPrint(OutFile," %8.3lf",TCSCptr->Control*Sn);
      else if (!strcmp(TCSCptr->Cont,"I")) fCustomPrint(OutFile," %8.3lf",TCSCptr->Control*Sn/(sqrt(3.0)*Vn));
      else if (!strcmp(TCSCptr->Cont,"D")) fCustomPrint(OutFile," %8.3lf",TCSCptr->Control/K3);
      fCustomPrint(OutFile,"\n\n");
    }
  }

  if(dataPtr->STATCOMbus!=NULL){
    fCustomPrint(OutFile,"\n\n                                ***** STATCOM RESULTS *****\n\n");
    fCustomPrint(OutFile,"                                  L=lower limit    H=higher limit\n\n");
    fCustomPrint(OutFile,"---|--|----|------------|------|------|----|------------|------|--------|--------|--------|--------|--------|--------|--------------|--------\n");
    fCustomPrint(OutFile,"  i  A    j Cont. Bus      Vj     dj      l STATCOM Bus   Vinv     K       Vdc     Alpha      P        Q        I      I (pu w.r.t.    Theta \n");
    fCustomPrint(OutFile,"     n      Name          (pu)   (deg)      Name          (pu)    (pu)     (pu)    (deg)     (MW)    (MVar)    (kA)    STATCOM base)   (deg) \n");
    i=0;
    for(STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
      fCustomPrint(OutFile,"---|--|----|------------|------|------|----|------------|------|--------|--------|--------|--------|--------|--------|--------------|--------\n");
      fCustomPrint(OutFile,"%3d ",++i);
      if(STATCOMptr->From->Area!=NULL) fCustomPrint(OutFile,"%2d ",STATCOMptr->From->Area->N);
      else fCustomPrint(OutFile,"%2d ",0);
      fCustomPrint(OutFile,"%4d %12s ",STATCOMptr->Ctrl->Num,STATCOMptr->Ctrl->Name);
      fCustomPrint(OutFile,"%6.4lf %6.2lf ",STATCOMptr->Ctrl->V,STATCOMptr->Ctrl->Ang/K3);
      fCustomPrint(OutFile,"%4d %12s ",STATCOMptr->From->Num,STATCOMptr->Name);
      fCustomPrint(OutFile,"%6.4lf ",STATCOMptr->k*STATCOMptr->Vdc);
      fCustomPrint(OutFile,"%8.6lf ",STATCOMptr->k);
      KV=STATCOMptr->From->KV;
      fCustomPrint(OutFile,"%8.6lf ",STATCOMptr->Vdc);
      fCustomPrint(OutFile,"%8.3lf ",STATCOMptr->alpha/K3);
      fCustomPrint(OutFile,"%8.4lf ",STATCOMptr->P*Sn);
      fCustomPrint(OutFile,"%8.2lf ",STATCOMptr->Q*Sn);
      Q=STATCOMptr->Q;
      fCustomPrint(OutFile,"%8.3lf ",STATCOMptr->I*Sn/(sqrt(3.0)*KV));
      fCustomPrint(OutFile,"      %8.3lf",STATCOMptr->I*Sn/STATCOMptr->MVA);
      if(STATCOMptr->I>=STATCOMptr->Imin && Q<0)      fCustomPrint(OutFile,"L");
      else if(STATCOMptr->I>=STATCOMptr->Imax && Q>0) fCustomPrint(OutFile,"H");
      else                                            fCustomPrint(OutFile," ");
      theta=STATCOMptr->theta;
      if (theta>=0) vals=1.00;
      else          vals=-1.00;
      if (fabs(theta)>2*PI) theta=theta-vals*floor(fabs(theta)/(2*PI))*2*PI;
      if (fabs(theta)>PI) theta=theta-vals*2*PI;
      fCustomPrint(OutFile,"%8.3lf ",theta/K3);
      fCustomPrint(OutFile,"\n\n");
    }
  }

  if (OutFile!=NULL) fclose(OutFile);
#ifdef WINDOWS
  delete[] File1;
#else
  free(File1);
#endif

}


/* --------------- WriteSolution ----------------- */
#ifdef ANSIPROTO
void WriteSolution(INDEX Iter,char *File,char *str)
#else
void WriteSolution(Iter,File,str)
INDEX Iter;
char *File,*str;
#endif
{
  //fclose(OutputHomot);
  if(!ExistParameter('s')) Output(Iter,File,str);
  if(ExistParameter('w')||ExistParameter('W')) IEEE();
  if(ExistParameter('J')||ExistParameter('j')) WriteJac();
}
