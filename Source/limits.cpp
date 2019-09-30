/* Check AC/DC limits. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"

#ifdef ANSIPROTO
BOOLEAN CheckRlimits(void);
BOOLEAN CheckVlimits(void);
BOOLEAN CheckQlimits(void);
BOOLEAN CheckDClimits(void);
void WriteSolution(INDEX Iter,char *File1,char *str);
#else
BOOLEAN CheckRlimits();
BOOLEAN CheckVlimits();
BOOLEAN CheckQlimits();
BOOLEAN CheckDClimits();
void WriteSolution();
#endif

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX MaxIter,Nac,NacEl,NregPQ,NregV,Ndc,Nslack,Nvolt,Narea,NacVar,Bl,NZvolt,NXvolt;
extern INDEX *ACvar;
extern VALUETYPE *dx,*dF,tol,Tol,Sn,lambda,*x0;
extern VALUETYPE K1,K2,MaxdFi,alpha;
extern IntegerVector *NewRow,*OldRow,*NewCol,*OldCol,*RowPartition,*ColPartition;
extern IntegerVector *RowPer,*ColPer;
extern BOOLEAN Acont,PQcont,QRcont,Rcont,
               PQlim,Tlim,Qlim,Ilim,Elim,Vlim,Zlim,Xlim,
               flagH,flagPoC,flagL,flagR,flagBS,flagPgMax,flagSmax;


/* -------------------- CheckRlimits ----------------------- */
#ifdef ANSIPROTO
BOOLEAN CheckRlimits(void)
#else
BOOLEAN CheckRlimits()
#endif
{
  BOOLEAN flag=FALSE;
  char str[4];
  ACbusData *ACptr;
  ElementList *ELptr,*Rptr;
  ElementData *Eptr;
  INDEX i,j;

  for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next){
    for (ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next){
      Eptr=ELptr->Eptr;
      if(Tlim&&!strcmp(Eptr->Type,"R")&&(Eptr->Tap<=1/Eptr->Tmax+0.00001||
                                         Eptr->Tap>=1/Eptr->Tmin-0.00001)){
        flag=TRUE; i=j=0; NregV--;
        if(Eptr->Tap<=1/Eptr->Tmax+0.00001) {Eptr->Tap=1/Eptr->Tmax; strcpy_s(str,"max");}
        else {Eptr->Tap=1/Eptr->Tmin; strcpy_s(str,"min");}
        for(Rptr=ACptr->Reg;Rptr!=NULL;Rptr=Rptr->Next){
          if(!strcmp(Rptr->Eptr->Type,"R")) i++;
          else if(!strcmp(Rptr->Eptr->Type,"RV")) j++;
        }
        if (i==1 && j==0) {
          strcpy_s(ACptr->Type,"B");
          if(ACptr->Area!=NULL && ACptr->Area->Slack==ACptr) strcat_s(ACptr->Type,"A");
          ACptr->Cont=ACptr;
          if (flagH) x0[ACvar[ACptr->N]+1]=ACptr->V;
        }
        else if(i==1 && j!=0) {
          strcpy_s(ACptr->Type,"BR");
          if(ACptr->Area!=NULL && ACptr->Area->Slack==ACptr) strcat_s(ACptr->Type,"A");
          ACptr->Cont=ACptr;
          if (flagH) x0[ACvar[ACptr->N]+1]=ACptr->V;
        }
        strcpy_s(Eptr->Type,"TR");
        fprintf(stderr,"***Warning: %1s Reg. transf. from %d %s to %d %s\n",
                Eptr->Ctype,Eptr->From->Num,Eptr->From->Name,Eptr->To->Num,Eptr->To->Name);
        fprintf(stderr,"            has been changed to a fixed tap transf. at t%3s.\n",str);
      }
      else if(Rcont&&!strcmp(Eptr->Type,"RV")&&(ACptr->V>=Eptr->Max||ACptr->V<=Eptr->Min)){
        flag=TRUE;
        strcpy_s(Eptr->Type,"R");
        if(ACptr->V>=Eptr->Max) {ACptr->V=Eptr->Max; strcpy_s(str,"max");}
        else {ACptr->V=Eptr->Min;  strcpy_s(str,"min");}
        if(!strpbrk(ACptr->Type,"T")){
          strcpy_s(ACptr->Type,"BT");
          if(ACptr->Area!=NULL && ACptr->Area->Slack==ACptr) strcat_s(ACptr->Type,"A");
          ACptr->Cont=NULL;
          if (flagH) x0[ACvar[ACptr->N]+1]=Eptr->Tap;
        }
        fprintf(stderr,"***Warning: %1s Reg. transf. from %d %s to %d %s\n",
				Eptr->Ctype,Eptr->From->Num,Eptr->From->Name,Eptr->To->Num,Eptr->To->Name);
        fprintf(stderr,"            has been changed to an LTC at V%3s for bus %d %s.\n",
                str,ACptr->Num,ACptr->Name);
      }
      else if(Tlim&&!strcmp(Eptr->Type,"RP")&&(Eptr->Ang>=Eptr->Tmax || Eptr->Ang<=Eptr->Tmin)){
        flag=TRUE; NregPQ--;
        strcpy_s(Eptr->Type,"TP");
        if(Eptr->Ang>=Eptr->Tmax) {Eptr->Ang=Eptr->Tmax; strcpy_s(str,"max");}
        else {Eptr->Ang=Eptr->Tmin; strcpy_s(str,"min");}
        i=ACvar[ACptr->N]+1+ACptr->Ncont-Eptr->Ncont;
        if (Acont && strpbrk(ACptr->Type,"A")) i++;
        if (flagH) x0[i]=Eptr->Cvar;
        fprintf(stderr,"***Warning: %1s Reg. transf. from %d %s to %d %s\n",
                Eptr->Ctype,Eptr->From->Num,Eptr->From->Name,Eptr->To->Num,Eptr->To->Name);
        fprintf(stderr,"            has been changed to a fixed angle transf. at a%3s.\n",str);
      }
      else if(Tlim&&!strcmp(Eptr->Type,"RQ")&&(Eptr->Tap<=1/Eptr->Tmax+0.00001||
                                                 Eptr->Tap>=1/Eptr->Tmin-0.00001)){
        flag=TRUE; NregPQ--;
        strcpy_s(Eptr->Type,"TQ");
        if(Eptr->Tap<=1/Eptr->Tmax+0.00001) {Eptr->Tap=1/Eptr->Tmax; strcpy_s(str,"max");}
        else {Eptr->Tap=1/Eptr->Tmin; strcpy_s(str,"min");}
        i=ACvar[ACptr->N]+1+ACptr->Ncont-Eptr->Ncont;
        if (Acont && strpbrk(ACptr->Type,"A")) i++;
        if (flagH) x0[i]=Eptr->Cvar;
        fprintf(stderr,"***Warning: %1s Reg. transf. from %d %s to %d %s\n",
                Eptr->Ctype,Eptr->From->Num,Eptr->From->Name,Eptr->To->Num,Eptr->To->Name);
        fprintf(stderr,"            has been changed to a fixed tap at t%3s.\n",str);
      }
      else if(PQcont&&PQlim&&strpbrk(Eptr->Type,"M")&&(Eptr->Cvar>=Eptr->Max||Eptr->Cvar<=Eptr->Min)){
        flag=TRUE;
        strcpy_s(Eptr->Type,"RP");
        if(Eptr->Cvar>=Eptr->Max) {Eptr->Cvar=Eptr->Max; strcpy_s(str,"max");}
        else {Eptr->Cvar=Eptr->Min; strcpy_s(str,"min");}
        i=ACvar[ACptr->N]+1+ACptr->Ncont-Eptr->Ncont;
        if (Acont && strpbrk(ACptr->Type,"A")) i++;
        if (flagH) x0[i]=Eptr->Ang;
        fprintf(stderr,"***Warning: %1s Reg. transf. from %d %s to %d %s\n",
                Eptr->Ctype,Eptr->From->Num,Eptr->From->Name,Eptr->To->Num,Eptr->To->Name);
        fprintf(stderr,"            has been changed to control %1s at P%3s.\n",Eptr->Ctype,str);
      }
      else if(PQcont&&PQlim&&strpbrk(Eptr->Type,"N")&&(Eptr->Cvar>=Eptr->Max||Eptr->Cvar<=Eptr->Min)){
        flag=TRUE;
        strcpy_s(Eptr->Type,"RQ");
        if(Eptr->Cvar>=Eptr->Max) {Eptr->Cvar=Eptr->Max; strcpy_s(str,"max");}
        else {Eptr->Cvar=Eptr->Min; strcpy_s(str,"min");}
        i=ACvar[ACptr->N]+1+ACptr->Ncont-Eptr->Ncont;
        if (Acont && strpbrk(ACptr->Type,"A")) i++;
        if (flagH) x0[i]=Eptr->Tap;
        fprintf(stderr,"***Warning: %1s Reg. transf. from %d %s to %d %s\n",
                Eptr->Ctype,Eptr->From->Num,Eptr->From->Name,Eptr->To->Num,Eptr->To->Name);
        fprintf(stderr,"            has been changed to control %1s at Q%3s.\n",Eptr->Ctype,str);
      }
    }
  }
  return(flag);
}

/* -------------------- CheckVlimits ----------------------- */
#ifdef ANSIPROTO
BOOLEAN CheckVlimits(void)
#else
BOOLEAN CheckVlimits()
#endif
{
  ACbusData *ACptr,*BEptr;
  char str[4];
  BOOLEAN flag=FALSE,Recover=TRUE;
  INDEX i;

  if (Narea<2){
    for(BEptr=dataPtr->ACbus;BEptr!=NULL;BEptr=BEptr->Next)
      if(strpbrk(BEptr->Type,"S")) break;
  }
  if(ExistParameter('G')) Recover=FALSE;
  if (Vlim || Elim || Ilim || Xlim) for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next){
    if (strpbrk(ACptr->Type,"V")) {
      if (Vlim && strpbrk(ACptr->cont,"Q") && (ACptr->V>=ACptr->Vmax || ACptr->V<=ACptr->Vmin)) {
        flag=TRUE; Nvolt++;
        if(ACptr->V>=ACptr->Vmax) {ACptr->V=ACptr->Vmax; strcpy_s(str,"max");}
        else                      {ACptr->V=ACptr->Vmin; strcpy_s(str,"min");}
        strcpy_s(ACptr->cont,"V");
        ACptr->VCont=ACptr->Qg;
        ACptr->Cont=NULL;
        if (flagH) x0[ACvar[ACptr->N]+1]=ACptr->Qg;
        fprintf(stderr,"***Warning: Generator %d %s will start controlling the voltage\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            due to V_%3s problems.\n",str);
      }
      else if (Vlim && Recover && strpbrk(ACptr->cont,"V") &&
               ((ACptr->V==ACptr->Vmax && ACptr->Qg>=ACptr->VCont)||
                (ACptr->V==ACptr->Vmin && ACptr->Qg<=ACptr->VCont)) ){
        flag=TRUE; Nvolt--;
        strcpy_s(ACptr->cont,"Q");
        ACptr->Qg=ACptr->VCont;
        ACptr->Cont=ACptr;
        if (flagH) x0[ACvar[ACptr->N]+1]=ACptr->V;
        fprintf(stderr,"***Warning: Generator %d %s voltage is back\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            within limits.\n",str);
      }
      else if (ACptr->Gen!=NULL) {
        if (Ilim && strpbrk(ACptr->cont,"Q") && ACptr->Gen->Ia>=ACptr->Gen->IaMax) {
          flag=TRUE;  Nvolt++;
          ACptr->Gen->Ia=ACptr->Gen->IaMax;
          strcpy_s(ACptr->cont,"I");
          ACptr->VCont=ACptr->Qg;
          ACptr->Cont=ACptr;
          if (flagH) x0[ACptr->Gen->Nvar+11]=ACptr->Qg;
          fprintf(stderr,"***Warning: Generator %d %s has lost Q control\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            due to Ia_max limit problems.\n");
        }
        else if (Ilim && Recover && strpbrk(ACptr->cont,"I") &&
                 ACptr->Gen->Ia==ACptr->Gen->IaMax && ACptr->Qg>=ACptr->VCont) {
          flag=TRUE; Nvolt--;
          strcpy_s(ACptr->cont,"Q");
          ACptr->Cont=ACptr;
          if (flagH) x0[ACptr->Gen->Nvar+11]=ACptr->Gen->Ia;
          fprintf(stderr,"***Warning: Generator %d %s has recovered\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            Q control as Ia is again within limits.\n");
        }
        else if (Elim && strpbrk(ACptr->cont,"Q") &&
                 (ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)) {
          flag=TRUE; Nvolt++;
          if (ACptr->Gen->Eq>=ACptr->Gen->EqMax) {
            ACptr->Gen->Eq=ACptr->Gen->EqMax;
            strcpy_s(str,"max");
          } else {
            ACptr->Gen->Eq=ACptr->Gen->EqMin;
            strcpy_s(str,"min");
          }
          strcpy_s(ACptr->cont,"E");
          ACptr->VCont=ACptr->Qg;
          ACptr->Cont=ACptr;
          if (flagH) x0[ACptr->Gen->Nvar+1]=ACptr->Qg;
          fprintf(stderr,"***Warning: Generator %d %s has lost Q control\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            due to Eq_%3s limit problems.\n",str);
        }
        else if (Elim && Recover && strpbrk(ACptr->cont,"E") &&
                 ((ACptr->Gen->Eq==ACptr->Gen->EqMax && ACptr->Qg>=ACptr->VCont) ||
                  (ACptr->Gen->Eq==ACptr->Gen->EqMin && ACptr->Qg<=ACptr->VCont)) ){
          flag=TRUE; Nvolt--;
          strcpy_s(ACptr->cont,"Q");
          ACptr->Qg=ACptr->VCont;
          ACptr->Cont=ACptr;
          if (flagH) x0[ACptr->Gen->Nvar+1]=ACptr->Gen->Eq;
          fprintf(stderr,"***Warning: Generator %d %s has recovered\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            Q control as Eq is again within limits.\n");
        }
        else if (strpbrk(ACptr->cont,"V") &&
                 ((Ilim && ACptr->Gen->Ia>=ACptr->Gen->IaMax) ||
                  (Elim &&(ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at a V limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing Ia and/or Eq limits.\n");
          WriteSolution(0,TrueParamStr(2),"Ia/Eq Limit Problems:");
          exit(1);
        }
        else if (strpbrk(ACptr->cont,"I") &&
                 ((Vlim && (ACptr->V>=ACptr->Vmax || ACptr->V<=ACptr->Vmin)) ||
                  (Elim && (ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at an Ia limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing V and/or Eq limits.\n");
          WriteSolution(0,TrueParamStr(2),"V/Eq Limit Problems:");
          exit(1);
        }
        else if (strpbrk(ACptr->cont,"E") &&
                 ((Ilim && ACptr->Gen->Ia>=ACptr->Gen->IaMax) ||
                  (Vlim && (ACptr->V>=ACptr->Vmax || ACptr->V<=ACptr->Vmin)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at an Eq limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing V and/or Ia limits.\n");
          WriteSolution(0,TrueParamStr(2),"V/Ia Limit Problems:");
          exit(1);
        }
      }
    }
    else if (strpbrk(ACptr->Type,"X")) {
      if (Xlim && strpbrk(ACptr->cont,"X") && ACptr->V<=ACptr->Vmin) {
        flag=TRUE;
        if (ACptr->Cont->Bx[0]>0) {
          ACptr->Cont->step++;
          i=ACptr->Cont->step;
          if (i<=ACptr->Cont->steps) {
            ACptr->Cont->B=ACptr->Cont->B+ACptr->Cont->Bx[i];
            fprintf(stderr,"***Warning: Reactance-controlled bus %d %s requires a %d st/nd/rd/th\n",ACptr->Num,ACptr->Name,i);
            fprintf(stderr,"            increment of %4.0lf MVar due to V_min problems.\n",Sn*ACptr->Cont->Bx[i]);
          }
        } else {
          i=ACptr->Cont->step;
          if (i>=1) {
            ACptr->Cont->B=ACptr->Cont->B-ACptr->Cont->Bx[i];
            fprintf(stderr,"***Warning: Reactance-controlled bus %d %s requires a %d st/nd/rd/th\n",ACptr->Num,ACptr->Name,i);
            fprintf(stderr,"            reduction of %4.0lf MVar due to V_min problems.\n",-Sn*ACptr->Cont->Bx[i]);
            ACptr->Cont->step--;
          }
        }
        if (ACptr->Cont->step>ACptr->Cont->steps || ACptr->Cont->step<1) {
          NXvolt--;
          strcpy_s(ACptr->cont,"m");
          fprintf(stderr,"***Warning: Reactance-controlled bus %d %s has run out of \n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            reactance support due to V_min problems.\n");
        }
      }
      else if (Xlim && strpbrk(ACptr->cont,"X") && ACptr->V>=ACptr->Vmax) {
        flag=TRUE;
        if (ACptr->Cont->Bx[0]>0) {
          i=ACptr->Cont->step;
          if (i>=1) {
            ACptr->Cont->B=ACptr->Cont->B-ACptr->Cont->Bx[i];
            fprintf(stderr,"***Warning: Reactance-controlled bus %d %s requires a %d st/nd/rd/th\n",ACptr->Num,ACptr->Name,i);
            fprintf(stderr,"            reduction of %4.0lf MVar due to V_max problems.\n",Sn*ACptr->Cont->Bx[i]);
            ACptr->Cont->step--;
          }
        } else {
          ACptr->Cont->step++;
          i=ACptr->Cont->step;
          if (i<=ACptr->Cont->steps) {
            ACptr->Cont->B=ACptr->Cont->B+ACptr->Cont->Bx[i];
            fprintf(stderr,"***Warning: Reactance-controlled bus %d %s requires a %d st/nd/rd/th\n",ACptr->Num,ACptr->Name,i);
            fprintf(stderr,"            increment of %4.0lf MVar due to V_max problems.\n",-Sn*ACptr->Cont->Bx[i]);
          }
        }
        if (ACptr->Cont->step>ACptr->Cont->steps || ACptr->Cont->step<1) {
          NXvolt--;
          strcpy_s(ACptr->cont,"M");
          fprintf(stderr,"***Warning: Reactance-controlled bus %d %s has run out of \n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            reactance support due to V_max problems.\n");
        }
      }
      else if (Xlim && Recover && strpbrk(ACptr->cont,"M") && ACptr->V<=ACptr->Vmax) {
        flag=TRUE; NXvolt++;
        strcpy_s(ACptr->cont,"X");
        fprintf(stderr,"***Warning: Reactance-controlled bus %d %s has recovered \n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            reactance support as V < V_max.\n");
      }
      else if (Xlim && Recover && strpbrk(ACptr->cont,"m") && ACptr->V>=ACptr->Vmin) {
        flag=TRUE; NXvolt++;
        strcpy_s(ACptr->cont,"X");
        fprintf(stderr,"***Warning: Reactance-controlled bus %d %s has recovered \n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            reactance support as V > V_min.\n");
      }
    }
  }
  return(flag);
}

/* -------------------- CheckQlimits ----------------------- */
#ifdef ANSIPROTO
BOOLEAN CheckQlimits(void)
#else
BOOLEAN CheckQlimits()
#endif
{
  ACbusData *ACptr;
  char str[5],Qmax[5],Qmin[5];
  BOOLEAN flag=FALSE,Recover=TRUE,RemoteVlost=FALSE;

  if(ExistParameter('G')) Recover=FALSE;
  if (Qlim || Elim || Ilim || Zlim) for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
    if (ACptr->flagPgMax==1) {
      fprintf(stderr,"***Warning: Generating bus %d %s has reached the maximum P/S limit;\n",ACptr->Num,ACptr->Name);
      fprintf(stderr,"            hence, Pg will be fixed at its max. MW limit -> %lf.\n",ACptr->Pmax*Sn);
      ACptr->flagPgMax++;
    }
    if (ACptr->Qmax==ACptr->Max) strcpy_s(Qmax,"Qmax");
    else                         strcpy_s(Qmax,"Smax");
    if (ACptr->Qmin==ACptr->Min) strcpy_s(Qmin,"Qmin");
    else                         strcpy_s(Qmin,"Smin");
    if (QRcont && strpbrk(ACptr->Type,"G")) {
      if (Qlim && strpbrk(ACptr->cont,"V") && (ACptr->Qg>=ACptr->Max||ACptr->Qg<=ACptr->Min)) {
        flag=TRUE; Nvolt--;
        if(ACptr->Qg>=ACptr->Max) {ACptr->Qg=ACptr->Max; strcpy_s(str,Qmax);}
        else                      {ACptr->Qg=ACptr->Min; strcpy_s(str,Qmin);}
        strcpy_s(ACptr->cont,"Q");
        fprintf(stderr,"***Warning: Generator %d %s has lost remote voltage control\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            due to %4s problems.\n",str);
        ACptr->Cont->Kbg--;
        if (ACptr->Cont->Kbg<1) {
          RemoteVlost=TRUE;
          ACptr->Cont->Kbg=-1;
          ACptr->Cont->Cont=ACptr->Cont;
          ACptr->Cont->VCont=ACptr->Cont->V;
          ACptr->Cont->Qr=ACptr->Qg/ACptr->Kbg;
          if (flagH) x0[ACvar[ACptr->Cont->N]+1]=ACptr->Cont->V;
          fprintf(stderr,"***Warning: Remote voltage controlled bus %d %s\n",ACptr->Cont->Num,ACptr->Cont->Name);
          fprintf(stderr,"            has lost all generator support.\n");
        }
      }
      else if (Qlim && Recover && strpbrk(ACptr->cont,"Q") &&  ACptr->Cont->Kbg==0 &&
               ((ACptr->Qg>=ACptr->Max && ACptr->Cont->V>=ACptr->Cont->VCont)||
                (ACptr->Qg<=ACptr->Min && ACptr->Cont->V<=ACptr->Cont->VCont)) ) {
        flag=TRUE; Nvolt++;
        strcpy_s(ACptr->cont,"V");
        ACptr->Cont->Kbg++;
        fprintf(stderr,"***Warning: Generator %d %s has recovered\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            remote voltage control as Qg/Sg is again within limits.\n");
        if (ACptr->Cont->Kbg==1) {
          ACptr->Cont->Cont=NULL;
          ACptr->Cont->V=ACptr->Cont->VCont;
          ACptr->Cont->Qr=ACptr->Qg/ACptr->Kbg;
          if (flagH) x0[ACvar[ACptr->Cont->N]+1]=ACptr->Cont->Qr;
          fprintf(stderr,"***Warning: Remote voltage controlled bus %d %s\n",ACptr->Cont->Num,ACptr->Cont->Name);
          fprintf(stderr,"            has recovered generator support.\n");
        }
      }
      else if (ACptr->Gen!=NULL) {
        if (Ilim && strpbrk(ACptr->cont,"V") && ACptr->Gen->Ia>=ACptr->Gen->IaMax) {
          flag=TRUE; Nvolt--;
          ACptr->Gen->Ia=ACptr->Gen->IaMax;
          strcpy_s(ACptr->cont,"I");
          if (flagH) x0[ACptr->Gen->Nvar+11]=ACptr->Qg;
          fprintf(stderr,"***Warning: Generator %d %s has lost remote\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            V control due to Ia_max limit problems.\n");
          ACptr->Cont->Kbg--;
          if (ACptr->Cont->Kbg<1) {
            ACptr->Cont->Cont=ACptr->Cont;
            ACptr->Cont->VCont=ACptr->Cont->V;
            if (flagH) x0[ACvar[ACptr->Cont->N]+1]=ACptr->Cont->V;
            fprintf(stderr,"***Warning: Remote voltage controlled bus %d %s\n",ACptr->Cont->Num,ACptr->Cont->Name);
            fprintf(stderr,"            has lost all generator support.\n");
          }
        }
        else if (Ilim && Recover && strpbrk(ACptr->cont,"I") &&   ACptr->Cont->Kbg<1 &&
                 ACptr->Gen->Ia==ACptr->Gen->IaMax && ACptr->Cont->V>=ACptr->Cont->VCont) {
          flag=TRUE;  Nvolt++;
          strcpy_s(ACptr->cont,"V");
          if (flagH) x0[ACptr->Gen->Nvar+11]=ACptr->Gen->Ia;
          fprintf(stderr,"***Warning: Generator %d %s has recovered\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            remote V control as Ia is again within limits.\n");
          ACptr->Cont->Kbg++;
          if (ACptr->Cont->Kbg==1) {
            ACptr->Cont->Cont=NULL;
            ACptr->Cont->V=ACptr->Cont->VCont;
            if (flagH) x0[ACvar[ACptr->Cont->N]+1]=ACptr->Qg/ACptr->Kbg;
            fprintf(stderr,"***Warning: Remote voltage controlled bus %d %s\n",ACptr->Cont->Num,ACptr->Cont->Name);
            fprintf(stderr,"            has recovered generator support.\n");
          }
        }
        else if (Elim && strpbrk(ACptr->cont,"V") &&
                 (ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)) {
          flag=TRUE;  Nvolt--;
          if (ACptr->Gen->Eq>=ACptr->Gen->EqMax) {
            ACptr->Gen->Eq=ACptr->Gen->EqMax;
            strcpy_s(str,"max");
          } else {
            ACptr->Gen->Eq=ACptr->Gen->EqMin;
            strcpy_s(str,"min");
          }
          strcpy_s(ACptr->cont,"E");
          if (flagH) x0[ACptr->Gen->Nvar+1]=ACptr->Qg;
          fprintf(stderr,"***Warning: Generator %d %s has lost\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            remote V control due to Eq_%3s limit problems.\n",str);
          ACptr->Cont->Kbg--;
          if (ACptr->Cont->Kbg<1) {
            ACptr->Cont->Cont=ACptr->Cont;
            ACptr->Cont->VCont=ACptr->Cont->V;
            if (flagH) x0[ACvar[ACptr->Cont->N]+1]=ACptr->Cont->V;
            fprintf(stderr,"***Warning: Remote voltage controlled bus %d %s\n",ACptr->Cont->Num,ACptr->Cont->Name);
            fprintf(stderr,"            has lost all generator support.\n");
          }
        }
        else if (Elim && Recover && strpbrk(ACptr->cont,"E") &&   ACptr->Cont->Kbg<1  &&
                 ((ACptr->Gen->Eq==ACptr->Gen->EqMax && ACptr->Cont->V>=ACptr->Cont->VCont) ||
                  (ACptr->Gen->Eq==ACptr->Gen->EqMin && ACptr->Cont->V<=ACptr->Cont->VCont)) ){
          flag=TRUE;  Nvolt++;
          strcpy_s(ACptr->cont,"V");
          if (flagH) x0[ACptr->Gen->Nvar+1]=ACptr->Gen->Eq;
          fprintf(stderr,"***Warning: Generator %d %s has recovered\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            remote V control as Eq is again within limits.\n");
          ACptr->Cont->Kbg++;
          if (ACptr->Cont->Kbg==1) {
            ACptr->Cont->Cont=NULL;
            ACptr->Cont->V=ACptr->Cont->VCont;
            if (flagH) x0[ACvar[ACptr->Cont->N]+1]=ACptr->Qg/ACptr->Kbg;
            fprintf(stderr,"***Warning: Remote voltage controlled bus %d %s\n",ACptr->Cont->Num,ACptr->Cont->Name);
            fprintf(stderr,"            has recovered generator support.\n");
          }
        }
        else if (strpbrk(ACptr->cont,"Q") &&
                 ((Ilim && ACptr->Gen->Ia>=ACptr->Gen->IaMax) ||
                  (Elim && (ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at a Q limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing Ia and/or Eq limits.\n");
          WriteSolution(0,TrueParamStr(2),"Ia/Eq Limit Problems:");
          exit(1);
        }
        else if (strpbrk(ACptr->cont,"I") &&
                 ((Qlim && (ACptr->Qg>=ACptr->Max || ACptr->Qg<=ACptr->Min)) ||
                  (Elim && (ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at an Ia limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing Q and/or Eq limits.\n");
          WriteSolution(0,TrueParamStr(2),"Q/Eq Limit Problems:");
          exit(1);
        }
        else if (strpbrk(ACptr->cont,"E") &&
                 ((Ilim && ACptr->Gen->Ia>=ACptr->Gen->IaMax) ||
                  (Qlim && (ACptr->Qg>=ACptr->Max || ACptr->Qg<=ACptr->Min)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at an Eq limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing Q and/or Ia limits.\n");
          WriteSolution(0,TrueParamStr(2),"Q/Ia Limit Problems:");
          exit(1);
        }
      }
    }
    else if (strpbrk(ACptr->Type,"Q,S") || (!QRcont && strpbrk(ACptr->Type,"G")) ) {
      if (Qlim && strpbrk(ACptr->cont,"V") && (ACptr->Qg>=ACptr->Max || ACptr->Qg<=ACptr->Min)) {
        flag=TRUE; Nvolt--;
        if(ACptr->Qg>=ACptr->Max) {
          ACptr->Qg=ACptr->Max;
          strcpy_s(str,Qmax);
        }
        else {
          ACptr->Qg=ACptr->Min;
          strcpy_s(str,Qmin);
        }
        strcpy_s(ACptr->cont,"Q");
        ACptr->VCont=ACptr->V;
        if (strpbrk(ACptr->Type,"G")) ACptr->Cont=ACptr->ContBus->AC;
        else                          ACptr->Cont=ACptr;
        if (flagH) x0[ACvar[ACptr->N]+1]=ACptr->V;
        fprintf(stderr,"***Warning: Generator %d %s has lost voltage control\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            due to %4s problems.\n",str);
      }
      else if (Qlim && Recover && strpbrk(ACptr->cont,"Q") &&
               ((ACptr->Qg>=ACptr->Max && ACptr->V>=ACptr->VCont) ||
                (ACptr->Qg<=ACptr->Min && ACptr->V<=ACptr->VCont)) ){
        flag=TRUE; Nvolt++;
        ACptr->V=ACptr->VCont;
        strcpy_s(ACptr->cont,"V");
        ACptr->Cont=NULL;
        if (flagH) x0[ACvar[ACptr->N]+1]=ACptr->Qg;
        fprintf(stderr,"***Warning: Generator %d %s has recovered\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            voltage control as Qg/Sg is again within limits.\n");
      }
      else if (ACptr->Gen!=NULL) {
        if (Ilim && strpbrk(ACptr->cont,"V") && ACptr->Gen->Ia>=ACptr->Gen->IaMax) {
          flag=TRUE;  Nvolt--;
          ACptr->Gen->Ia=ACptr->Gen->IaMax;
          strcpy_s(ACptr->cont,"I");
          ACptr->VCont=ACptr->V;
          if (strpbrk(ACptr->Type,"G")) ACptr->Cont=ACptr->ContBus->AC;
          else                          ACptr->Cont=ACptr;
          if (flagH) {
            x0[ACvar[ACptr->N]+1]=ACptr->V;
            x0[ACptr->Gen->Nvar+11]=ACptr->Qg;
          }
          fprintf(stderr,"***Warning: Generator %d %s has lost V control\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            due to Ia_max limit problems.\n");
        }
        else if (Ilim && Recover && strpbrk(ACptr->cont,"I") &&
                 ACptr->Gen->Ia==ACptr->Gen->IaMax && ACptr->V>=ACptr->VCont) {
          flag=TRUE;  Nvolt++;
          strcpy_s(ACptr->cont,"V");
          ACptr->V=ACptr->VCont;
          ACptr->Cont=NULL;
          if (flagH) {
            x0[ACvar[ACptr->N]+1]=ACptr->Qg;
            x0[ACptr->Gen->Nvar+11]=ACptr->Gen->Ia;
          }
          fprintf(stderr,"***Warning: Generator %d %s has recovered\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            V control as Ia is again within limits.\n");
        }
        else if (Elim && strpbrk(ACptr->cont,"V") &&
                 (ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)) {
          flag=TRUE; Nvolt--;
          if (ACptr->Gen->Eq>=ACptr->Gen->EqMax) {
            ACptr->Gen->Eq=ACptr->Gen->EqMax;
            strcpy_s(str,"max");
          } else {
            ACptr->Gen->Eq=ACptr->Gen->EqMin;
            strcpy_s(str,"min");
          }
          strcpy_s(ACptr->cont,"E");
          ACptr->VCont=ACptr->V;
          if (strpbrk(ACptr->Type,"G")) ACptr->Cont=ACptr->ContBus->AC;
          else                          ACptr->Cont=ACptr;
          if (flagH) {
            x0[ACvar[ACptr->N]+1]=ACptr->V;
            x0[ACptr->Gen->Nvar+1]=ACptr->Qg;
          }
          fprintf(stderr,"***Warning: Generator %d %s has lost V control\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            due to Eq_%3s limit problems.\n",str);
        }
        else if (Elim && Recover && strpbrk(ACptr->cont,"E") &&
                 ((ACptr->Gen->Eq==ACptr->Gen->EqMax && ACptr->V>=ACptr->VCont) ||
                  (ACptr->Gen->Eq==ACptr->Gen->EqMin && ACptr->V<=ACptr->VCont)) ){
          flag=TRUE;  Nvolt++;
          strcpy_s(ACptr->cont,"V");
          ACptr->V=ACptr->VCont;
          ACptr->Cont=NULL;
          if (flagH) {
            x0[ACvar[ACptr->N]+1]=ACptr->Qg;
            x0[ACptr->Gen->Nvar+1]=ACptr->Gen->Eq;
          }
          fprintf(stderr,"***Warning: Generator %d %s has recovered\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"            V control as Eq is again within limits.\n");
        }
        else if (strpbrk(ACptr->cont,"Q") &&
                 ((Ilim && ACptr->Gen->Ia>=ACptr->Gen->IaMax) ||
                  (Elim && (ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at a Q limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing Ia and/or Eq limits.\n");
          WriteSolution(0,TrueParamStr(2),"Ia/Eq Limit Problems:");
          exit(1);
        }
        else if (strpbrk(ACptr->cont,"I") &&
                 ((Qlim && (ACptr->Qg>=ACptr->Max || ACptr->Qg<=ACptr->Min)) ||
                  (Elim && (ACptr->Gen->Eq>=ACptr->Gen->EqMax || ACptr->Gen->Eq<=ACptr->Gen->EqMin)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at an Ia limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing Q and/or Eq limits.\n");
          WriteSolution(0,TrueParamStr(2),"Q/Eq Limit Problems:");
          exit(1);
        }
        else if (strpbrk(ACptr->cont,"E") &&
                 ((Ilim && ACptr->Gen->Ia>=ACptr->Gen->IaMax) ||
                  (Qlim && (ACptr->Qg>=ACptr->Max || ACptr->Qg<=ACptr->Min)))) {
          fprintf(stderr,"***Error: Generator %d %s it's already at an Eq limit.\n",ACptr->Num,ACptr->Name);
          fprintf(stderr,"          Try changing Q and/or Ia limits.\n");
          WriteSolution(0,TrueParamStr(2),"Q/Ia Limit Problems:");
          exit(1);
        }
      }
    }
    else if (strpbrk(ACptr->Type,"Z")) {
      if (Zlim && strpbrk(ACptr->cont,"V") && (ACptr->Qg>=ACptr->Qmax || ACptr->Qg<=ACptr->Qmin)) {
        flag=TRUE; NZvolt--;
        if(ACptr->Qg>=ACptr->Qmax) {
          ACptr->Qg=ACptr->Qmax;
          strcpy_s(str,Qmax);
        }
        else {
          ACptr->Qg=ACptr->Qmin;
          strcpy_s(str,Qmin);
        }
        ACptr->Bz=ACptr->Qg/(ACptr->V*ACptr->V);
        strcpy_s(ACptr->cont,"Q");
        ACptr->VCont=ACptr->V;
        ACptr->Cont=ACptr;
        if (flagH) x0[ACvar[ACptr->N]+1]=ACptr->V;
        fprintf(stderr,"***Warning: Reactance-controlled bus %d %s has lost\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            voltage control due to %4s problems.\n",str);
      }
      else if (Zlim && strpbrk(ACptr->cont,"Q") &&
               ((ACptr->Qg>=ACptr->Qmax && ACptr->V>=ACptr->VCont) ||
                (ACptr->Qg<=ACptr->Qmin && ACptr->V<=ACptr->VCont)) ){
        flag=TRUE; NZvolt++;
        ACptr->V=ACptr->VCont;
        strcpy_s(ACptr->cont,"V");
        ACptr->Cont=NULL;
        if (flagH) x0[ACvar[ACptr->N]+1]=ACptr->Qg;
        fprintf(stderr,"***Warning: Reactance-controlled bus %d %s has recovered\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"            voltage control as Q is again within limits.\n");
      }
    }
  }
  
  if ((Qlim || Elim || Ilim || Zlim) && RemoteVlost && QRcont)
    for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
      if (strpbrk(ACptr->Type,"C") && ACptr->Kg<0) ACptr->Kbg=0;
    }

  return(flag);
}

/* -------------------- CheckDClimits ----------------------- */
#ifdef ANSIPROTO
BOOLEAN CheckDClimits(void)
#else
BOOLEAN CheckDClimits()
#endif
{
  DCbusData *DCptrR,*DCptrI;
  BOOLEAN flag=FALSE;

  for (DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next) if(!strcmp(DCptrR->Type,"R")) {
    DCptrI=DCptrR->To;
    if ((strcmp(DCptrR->Cont1,"AT") && strcmp(DCptrR->Cont2,"AT") &&
         (DCptrR->Tap>=DCptrR->TapMax || DCptrR->Tap<=DCptrR->TapMin)) ||
        (strcmp(DCptrI->Cont1,"AT") && strcmp(DCptrI->Cont2,"AT") &&
         (DCptrI->Tap>=DCptrI->TapMax || DCptrI->Tap<=DCptrI->TapMin))){
      flag=TRUE;
      fprintf(stderr,"***Warning: The program will release tap limits in the HVDC links.\n");
      break;
    }
  }
  return(flag);
}

