/* Read AC data in WSCC/EPRI/BPA/ANAREDE format
   (EPRI's P.F. V 5.0 input data formats). */

#include "readdata.h"

/* ------- Global Variables ------ */
VALUETYPE G,B,G1,B1,G2,B2,Tap,Ang;
extern FILE *InputDataFile;


/* ---------------- Multiply 2 Complex Numbers ----------------------------- */
#ifdef ANSIPROTO
void Multiply(VALUETYPE *a,VALUETYPE *b,VALUETYPE c,VALUETYPE d)
#else
void Multiply(a,b,c,d)
VALUETYPE *a,*b,c,d;
#endif
{
  VALUETYPE r,i;

  r=(*a)*c - (*b)*d;
  i=(*b)*c + (*a)*d;
  *a=r; *b=i;
}

/* ---------------- Divide 2 Complex Numbers ----------------------------- */
#ifdef ANSIPROTO
void Divide(VALUETYPE *a,VALUETYPE *b,VALUETYPE c,VALUETYPE d)
#else
void Divide(a,b,c,d)
VALUETYPE *a,*b,c,d;
#endif
{
  VALUETYPE D;

  D=c*c+d*d;
  if (!D) { c=d=0; D=1;}
  Multiply(a,b,c/D,-d/D);
}


/* ---------------- AddSection --------------------------- */
#ifdef ANSIPROTO
BOOLEAN AddSection(ACbusData *From,ACbusData *To,char *Line,char *Ckt,INDEX Sec)
#else
BOOLEAN AddSection(From,To,Line,Ckt,Sec)
ACbusData *From,*To;
char *Line,*Ckt;
INDEX Sec;
#endif
{
  ElementList *ELptr;
  ElementData *Eptr;
  VALUETYPE Ge,Be,Ge1,Be1,Ge2,Be2,Gs,Bs;

  for(ELptr=From->Elem;ELptr!=NULL;ELptr=ELptr->Next){
    Eptr=ELptr->Eptr;
    if (((From==Eptr->From&&To==Eptr->To)||(From==Eptr->To&&To==Eptr->From))
        &&!strcmp(Ckt,Eptr->Ckt)&&strcmp(Eptr->Zone,"")){
      if (!strncmp(Line,"T",1)&&!strpbrk(Eptr->Type,"RT")) strcpy(Eptr->Type,"T");
      if(Sec<Eptr->Sec) {
        Ge=Eptr->G;   Eptr->G=G/(Eptr->Tap*Eptr->Tap);   G=Ge;
        Be=Eptr->B;   Eptr->B=B/(Eptr->Tap*Eptr->Tap);   B=Be;
        Ge1=Eptr->G1; Eptr->G1=G1/(Eptr->Tap*Eptr->Tap); G1=Ge1;
        Be1=Eptr->B1; Eptr->B1=B1/(Eptr->Tap*Eptr->Tap); B1=Be1;
        Ge2=Eptr->G2; Eptr->G2=G2/(Eptr->Tap*Eptr->Tap); G2=Ge2;
        Be2=Eptr->B2; Eptr->B2=B2/(Eptr->Tap*Eptr->Tap); B2=Be2;
      } else {
        Eptr->G=Eptr->G/(Tap*Tap);
        Eptr->B=Eptr->B/(Tap*Tap);
        Eptr->G1=Eptr->G1/(Tap*Tap);
        Eptr->B1=Eptr->B1/(Tap*Tap);
        Eptr->G2=Eptr->G2/(Tap*Tap);
        Eptr->B2=Eptr->B2/(Tap*Tap);
      }
      if ((fabs(Eptr->G)<1e-7 && fabs(Eptr->B)<1e-7)||(fabs(G)<1e-7 && fabs(B)<1e-7)) {
         fprintf(stderr,"***Warning: A section of element %d %s - %d %s\n",
                  Eptr->From->Num,Eptr->From->Name,Eptr->To->Num,Eptr->To->Name);
         fprintf(stderr,"            has a zero impedance.  Check the section data for this element.\n");
      }
      if (!strncmp(Line,"T",1)&&!strpbrk(Eptr->Type,"RT")) strcpy(Eptr->Type,"T");
      Eptr->Tap*=Tap;
      Eptr->Ang+=Ang;
      Eptr->Sec=Sec;
      Gs=Eptr->G+G+Eptr->G2+G1;
      Bs=Eptr->B+B+Eptr->B2+B1;
      Ge=Ge1=Eptr->G;  Be=Be1=Eptr->B;
      Ge2=G;  Be2=B;
      Multiply(&Ge,&Be,G,B);
      Divide(&Ge,&Be,Gs,Bs);
      Multiply(&Ge1,&Be1,Eptr->G2+G1,Eptr->B2+B1);
      Divide(&Ge1,&Be1,Gs,Bs);
      Ge1+=Eptr->G1;  Be1+=Eptr->B1;
      Multiply(&Ge2,&Be2,Eptr->G2+G1,Eptr->B2+B1);
      Divide(&Ge2,&Be2,Gs,Bs);
      Ge2+=G2;       Be2+=B2;
      Eptr->G=Ge;    Eptr->B=Be;
      Eptr->G1=Ge1;  Eptr->B1=Be1;
      Eptr->G2=Ge2;  Eptr->B2=Be2;
      return(TRUE);
    }
  }
  return(FALSE);
}


/* ---------------- ReadWSCC ----------------------------- */
#ifdef ANSIPROTO
void ReadWSCC()
#else
void ReadWSCC()
#endif
/* Read Bus and Element data in WSCC format. */
{
  ACbusData *ACptr,*ACptrp,*ACptrs;
  ElementData *Eptr;
  AreaData *Aptr;
  char Line[BUFLEN],Name[31],str[6];
  VALUETYPE KV,KVp,KVs,R,X,Bx,Taps,Tap1,Tap2,Imax;
  int i,j,k,s;
  INDEX Sec;
  BOOLEAN flag=FALSE,flagPrint=TRUE;


  for(i=0;i<=2;strcpy(dataPtr->Title[i],"\n"),i++);
  Sn=100.0;
  RealParameter('$',&Sn,1.0,100000000.0);
  ACptr=NULL;
  for(;;){ /* Reading Loop */
    if (fgets(Line,BUFLEN,InputDataFile)==NULL){
      ErrorHalt("Missing END or Anarede's FIM card.");
      break;
    }
    LineNum++;
    if (!strncmp(Line,"C",1) || !strncmp(Line,"BAS",3)) continue;

  /* ------------------------ Title ----------------------------- */
    else if (!strncmp(Line,"HDG",3)||!strncmp(Line,"TITU",4)) {
      i=0;
      for(;;){ /* Title Reading Loop */
        if (fgets(Line,BUFLEN,InputDataFile)==NULL){
          ErrorHalt("Missing BAS or Anarede's D title cards.");
          break;
        }
        LineNum++;
        if (strncmp(Line,"BAS",3)&&strncmp(Line,"D",1)) {
          strcpy(dataPtr->Title[i],Line);
          i++;
          if (i>2) break;
        }
        else break;
      }
    }

  /* --------------- AC bus data -------------------------------- */
    else if (!strncmp(Line,"B  ",3) || !strncmp(Line,"BQ ",3) ||
        !strncmp(Line,"BV ",3) || !strncmp(Line,"BE ",3) ||
        !strncmp(Line,"BC ",3) || !strncmp(Line,"BG ",3) ||
        !strncmp(Line,"BT ",3) || !strncmp(Line,"BX ",3) ||
        !strncmp(Line,"BA ",3) || !strncmp(Line,"BF ",3) ||
        !strncmp(Line,"BJ ",3) || !strncmp(Line,"BK ",3) ||
        !strncmp(Line,"BL ",3) || !strncmp(Line,"BS ",3)){
      GetStr(Line,7,12,12,Name);
      KV=GetValue(Line,15,4,0);
      ACptr=ACbusInList(0,Name,KV,Nac,1);
      if (ACptr->V>0) {
         fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
         ErrorHalt("The AC bus was previously defined (check B cards).");
      }
      if (ACptr->N==0) { Nac++; ACptr->Num=ACptr->N=Nac;}
      GetStr(Line,19,2,2,ACptr->Zone);
      GetStr(Line,4,3,3,ACptr->Owner);
      ACptr->V=GetValue(Line,58,4,3);
      if (ACptr->V<=0) ACptr->V=1;
      ACptr->VCont=ACptr->V;
      ACptr->Pl=GetValue(Line,21,5,0)/Sn;
      ACptr->Ql=GetValue(Line,26,5,0)/Sn;
      ACptr->G=GetValue(Line,31,4,0)/Sn;
      ACptr->B=GetValue(Line,35,4,0)/Sn;
      ACptr->Pg=GetValue(Line,43,5,0)/Sn;
      ACptr->PgMax=GetValue(Line,39,4,0)/Sn;
      if (ACptr->PgMax==0) ACptr->PgMax=99999999.;
      if (ACptr->PgMax<ACptr->Pg) {
         ACptr->PgMax=99999999.;
         fprintf(stderr,"***Warning: Bus %d %s has its maximum generating power PgMax < Pg.\n",ACptr->N,ACptr->Name);
         fprintf(stderr,"            PgMax will be given value of 99999999.\n");
      }
      ACptr->Pmax=ACptr->PgMax;
      ACptr->Qg=GetValue(Line,48,5,0)/Sn;
      if ((!strncmp(Line,"B  ",3)||!strncmp(Line,"BJ ",3)||!strncmp(Line,"BF ",3))
          && strcmp(ACptr->Type,"BC")) {
         ACptr->Cont=ACptr;
         if (!strncmp(Line,"BJ ",3)){
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            fprintf(stderr,"***Warning: BJ bus has been transformed to a B bus.\n");
         }
         if (!strncmp(Line,"BF ",3)){
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            fprintf(stderr,"***Warning: BF bus has been transformed to a B bus.\n");
         }
      }
      else if (!strncmp(Line,"BC ",3)) strcpy(ACptr->Type,"BC");
      else if (!strncmp(Line,"BT ",3)) strcpy(ACptr->Type,"BT");
      else if (!strncmp(Line,"BE ",3) || !strncmp(Line,"BK ",3)){
         Nvolt++;
         ACptr->Qg=0;
         strcpy(ACptr->Type,"BQ");
         strcpy(ACptr->cont,"V");
         ACptr->Qmax=99999999.;
         ACptr->Qmin= -99999999.;
         if (!strncmp(Line,"BK ",3)){
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            fprintf(stderr,"***Warning: BK bus has been transformed to a BE bus.\n");
         }
      }
      if (!strncmp(Line,"BQ ",3) || !strncmp(Line,"BK ",3)){
         Nvolt++;
         ACptr->Qg=0;
         strcpy(ACptr->Type,"BQ");
         strcpy(ACptr->cont,"V");
         ACptr->Qmax=GetValue(Line,48,5,0)/Sn;
         ACptr->Qmin=GetValue(Line,53,5,0)/Sn;
         if (ACptr->Qmax<=ACptr->Qmin) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("AC bus Q limits are wrong: Qmin >= Qmax.");
         }
         if (!strncmp(Line,"BL ",3)){
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            fprintf(stderr,"***Warning: BL bus has been transformed to a BQ bus.\n");
         }
      }
      else if (!strncmp(Line,"BG ",3)){
         Nvolt++;
         ACptr->Qg=0;
         strcpy(ACptr->Type,"BG");
         strcpy(ACptr->cont,"V");
         ACptr->Qmax=GetValue(Line,48,5,0)/Sn;
         ACptr->Qmin=GetValue(Line,53,5,0)/Sn;
         if (ACptr->Qmax<=ACptr->Qmin) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("AC bus Q limits are wrong: Qmin >= Qmax.");
         }
         GetStr(Line,66,12,12,Name);
         KV=GetValue(Line,74,4,0);
         if (!strcmp(Name,"            ")) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("The controlled bus has not been defined.");
         }
         ACptrp=ACbusInList(0,Name,KV,Nac,1);
         if (ACptrp->N==0) { Nac++; ACptrp->Num=ACptrp->N=Nac;}
         if (!strcmp(ACptr->Name,"            ")) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("The controlled bus has not been defined.");
         }
         if(!strcmp(ACptrp->Type,"B")) strcpy(ACptrp->Type,"BC");
         else if(strcmp(ACptrp->Type,"BC")) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("The voltage controlled bus is not a PQ bus.");
         }
         ACptr->Cont=ACptrp;
         ACptrp->Cont=NULL;
         ACptrp->Kbg++;
         ACptr->Kbg=GetValue(Line,78,3,0)/100;
         if (flag2Vcontrol) {
           ACptrp->Kbg1=ACptrp->Kbg1+ACptr->Qmax;
           ACptrp->Kbg2=ACptrp->Kbg2+ACptr->Qmin;
         }
         if (ACptr->Kbg<=0) ACptr->Kbg=1.;
      }
      else if (!strncmp(Line,"BV ",3)||!strncmp(Line,"BA ",3)) {
         strcpy(ACptr->Type,"BV");
         strcpy(ACptr->cont,"Q");
         ACptr->VCont=ACptr->Qg;
         ACptr->Vmax=GetValue(Line,58,4,3);
         if (ACptr->Vmax<=0) ACptr->Vmax=10.;
         ACptr->Vmin=GetValue(Line,62,4,3);
         if (ACptr->Vmin<=0) ACptr->Vmin=0.001;
         if (ACptr->Vmax<=ACptr->Vmin) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("AC bus V limits are wrong: Vmin >= Vmax.");
         }
         ACptr->Qmax=99999999.;
         ACptr->Qmin= -99999999.;
         ACptr->Cont=ACptr;
         if (!strncmp(Line,"BA ",3)){
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            fprintf(stderr,"***Warning: BA bus has been transformed to a BV bus.\n");
         }
      }
      else if (!strncmp(Line,"BX ",3)) {
         ACptr->Vmax=GetValue(Line,58,4,3);
         if (ACptr->Vmax<=0) ACptr->Vmax=10.;
         ACptr->Vmin=GetValue(Line,62,4,3);
         if (ACptr->Vmin<=0) ACptr->Vmin=0.001;
         if (ACptr->Vmax<=ACptr->Vmin) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("AC bus V limits are wrong: Vmin >= Vmax.");
         }
         ACptr->Cont=ACptr;
      }
      else if (!strncmp(Line,"BS ",3)) {
         strcat_s(ACptr->Type,"S");
         strcpy(ACptr->cont,"V");
         Nslack++;
         ACptr->Qg=0;
         AngSlack=0;
         ACptr->Qmax=GetValue(Line,48,5,0)/Sn;
         ACptr->Qmin=GetValue(Line,53,5,0)/Sn;
         if (ACptr->Qmax<=ACptr->Qmin) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("AC bus Q limits are wrong: Qmin >= Qmax.");
         }
      }
    }

  /* --------------- BPA '+' Cards -------------------------------- */
  /* These cards are assumed to come right after a bus card.
     The program associates the card to the last bus defined in
     the data set.  */
    else if (!strncmp(Line,"+",1) && ACptr!=NULL) {
      ACptr->Pl+=GetValue(Line,21,5,0)/Sn;
      ACptr->Ql+=GetValue(Line,26,5,0)/Sn;
      ACptr->G+=GetValue(Line,31,4,0)/Sn;
      ACptr->B+=GetValue(Line,35,4,0)/Sn;
      ACptr->Pg+=GetValue(Line,43,5,0)/Sn;
      ACptr->Ql-=GetValue(Line,39,4,0)/Sn;
    }

  /* --------------- X data cards -------------------------------- */
    else if (!strncmp(Line,"X  ",3)){
      flagPrint=TRUE;
      GetStr(Line,7,12,12,Name);
      KV=GetValue(Line,15,4,0);
      ACptr=ACbusInList(0,Name,KV,Nac,1);
      if (ACptr->N==0) { Nac++; ACptr->Num=ACptr->N=Nac;}
      GetStr(Line,21,12,12,Name);
      KVp=GetValue(Line,29,4,0);
      if (KVp>0) {
        ACptrp=ACbusInList(0,Name,KVp,Nac,1);
        if (ACptrp->N==0) { Nac++; ACptrp->Num=ACptr->N=Nac;}
      }
      else ACptrp=ACptr;
      if (Xcont) ACptrp=ACptr;
      if (!strcmp(ACptrp->Type,"B")) {
        NXvolt++;
        strcpy(ACptrp->Type,"BX");
        if (ACptrp!=ACptr) strcpy(ACptr->Type,"B");
        ACptr->V=1;
        ACptrp->Cont=ACptr;
        strcpy(ACptrp->cont,"X");
        for(i=1,j=33; i<=8; i++,j=j+6) {
          k=GetInt(Line,j,1);
          Bx=GetValue(Line,j+1,5,0);
          if (k>0 && k<10 && Bx!=0) for(s=1;s<=k;s++) {
             ACptr->steps++;
             ACptr->Bx[ACptr->steps]=fabs(Bx)/Sn;
             if (ACptrp->Bx[0]==0) {
               if (Bx>0) ACptr->Bx[0]=1;
               else      ACptr->Bx[0]=-1;
             } else {
               if (ACptr->Bx[0]*Bx<0 && flagPrint) {
                 fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
                 fprintf(stderr,"***Warning: All MVAr compensation steps are assumed to be either positive or negative,\n");
                 fprintf(stderr,"            as defined by the first nonzero step on this card. \n");
                 flagPrint=FALSE;
               }
             }
          }
        }
      } else {
        fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
        fprintf(stderr,"***Warning: This bus has been previously defined as a type %s bus, and hence\n",ACptrp->Type);
        fprintf(stderr,"            it cannot be defined as a BX remote bus. The BX data will be ignored.\n");
      }

    }

  /* --------------- AC element data -------------------------------- */
    else if (!strncmp(Line,"L  ",3)||!strncmp(Line,"T  ",3)||
        !strncmp(Line,"E  ",3)){
      GetStr(Line,7,12,12,Name);
      KV=GetValue(Line,15,4,0);
      if (KV==0){
         fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
         ErrorHalt("Base voltage at bus 1 is zero.");
         KV=1;
      }
      ACptr=ACbusInList(0,Name,KV,Nac,1);
      if (ACptr->N==0) { Nac++; ACptr->Num=ACptr->N=Nac;}
      GetStr(Line,20,12,12,Name);
      KVp=GetValue(Line,28,4,0);
      if (KVp==0){
         fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
         ErrorHalt("Base voltage at bus 2 is zero.");
         KVp=1;
      }
      ACptrp=ACbusInList(0,Name,KVp,Nac,1);
      if (ACptrp->N==0) { Nac++; ACptrp->Num=ACptrp->N=Nac;}
      if (ACptr==ACptrp){
         fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
         ErrorHalt("Both AC element buses are the same.");
      }
      R=GetValue(Line,39,6,5);
      X=GetValue(Line,45,6,5);
      if (fabs(R)<0.0001 && fabs(X)<0.0001) {
         fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
         ErrorHalt("AC element is a short circuit. Try eliminating it.");
         G=B=0;
      } else {
         G= R/(R*R+X*X);
         B= -X/(R*R+X*X);
      }
      G1=GetValue(Line,51,6,5);
      B1=GetValue(Line,57,6,5);
      Tap=Taps=1;
      Ang=0;
      if (!strncmp(Line,"E  ",3)) {
         G2=GetValue(Line,63,6,5);
         B2=GetValue(Line,69,6,5);
         Imax=GetValue(Line,34,4,0)/(Sn*1000/(sqrt(3.0)*KV));
      }
      else if (!strncmp(Line, "L  ",3)) {
        G2=G1;
        B2=B1;
        Imax=GetValue(Line,34,4,0)/(Sn*1000/(sqrt(3.0)*KV));
      } else {
         G2=G1; B2=B1;
         G1=B1=0;
         Tap2=GetValue(Line,68,5,2);
         if (Tap2>0) {
            Tap1=GetValue(Line,63,5,2)/KV;
            if (Tap1<=0) Tap1=1;
            Taps=Tap2/KVp;
            Tap=Taps/Tap1;
         } else Ang=GetValue(Line,63,5,2)*K3;
         Imax=GetValue(Line,34,4,0)/Sn;
      }
      GetStr(Line,32,1,1,str);
      Sec=GetInt(Line,33,1);
      if (Sec) flag=AddSection(ACptr,ACptrp,Line,str,Sec);
      else flag=FALSE;
      if (!flag) {
          if (!strncmp(Line,"T  ",3)) {
             Eptr=ElemInList(ACptr,ACptrp,NacEl,1,"R",str);
             if (Eptr==NULL) Eptr=ElemInList(ACptr,ACptrp,NacEl,0,"",str);
          }
          else Eptr=ElemInList(ACptr,ACptrp,NacEl,0,"",str);
          Eptr->Sec=Sec;
          Eptr->G=G;
          Eptr->B=B;
          Eptr->G1=G1;
          Eptr->B1=B1;
          Eptr->G2=G2;
          Eptr->B2=B2;
          Eptr->Tap=Tap;
          Eptr->Taps=Taps;
          Eptr->Ang=Ang;
          Eptr->Imax=Imax;
          GetStr(Line,4,3,3,Eptr->Owner);
          strcpy(Eptr->Zone,"  ");
          NacEl++;
          GetStr(Line,19,1,1,str);
          if (!strcmp(str,"1")) Eptr->Meter=ACptr;
          else if (!strcmp(str,"2")) Eptr->Meter=ACptrp;
          if (!strncmp(Line,"L  ",3)) strcpy(Eptr->Type,"L");
          else if (!strncmp(Line,"E  ",3)) strcpy(Eptr->Type,"E");
          else if (!strncmp(Line,"T  ",3)&&strncmp(Eptr->Type,"R",1)) strcpy(Eptr->Type,"T");
      }
    }

  /* ------------ Regulating transformer data ------------------- */
    else if (!strncmp(Line,"R  ",3)||!strncmp(Line,"RP ",3)||
        !strncmp(Line,"RM ",3)||!strncmp(Line,"RQ ",3)||
        !strncmp(Line,"RN ",3)){
      GetStr(Line,7,12,12,Name);
      KV=GetValue(Line,15,4,0);
      ACptr=ACbusInList(0,Name,KV,Nac,1);
      if (ACptr->N==0) { Nac++; ACptr->Num=ACptr->N=Nac;}
      GetStr(Line,20,12,12,Name);
      KVp=GetValue(Line,28,4,0);
      ACptrp=ACbusInList(0,Name,KVp,Nac,1);
      if (ACptrp->N==0) { Nac++; ACptrp->Num=ACptrp->N=Nac;}
      if (ACptr==ACptrp){
         fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
         ErrorHalt("Both AC element buses are the same.");
      }
      GetStr(Line,32,1,1,str);
      Eptr=ElemInList(ACptr,ACptrp,NacEl,1,"T",str);
      if (Eptr==NULL) Eptr=ElemInList(ACptr,ACptrp,NacEl,0,"",str);
      NacEl++; NregPQ++;
      GetStr(Line,34,12,12,Name);
      KVs=GetValue(Line,42,4,0);
      ACptrs=ACbusInList(0,Name,KVs,Nac,1);
      if (ACptrs->N==0) { Nac++; ACptrs->Num=ACptrs->N=Nac;}
      if (!strncmp(Line,"R  ",3)){
         NregV++;
         strcpy(Eptr->Type,"R");
         Eptr->Tmax=GetValue(Line,46,5,2)/KV;
         Eptr->Tmin=GetValue(Line,51,5,2)/KV;
         if(Eptr->Tmax<=0) Eptr->Tmax=1.1;
         if(Eptr->Tmin<=0) Eptr->Tmin=0.9;
         if(Eptr->Tmax<=Eptr->Tmin) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("LTC limits are wrong: Tmin > Tmax.");
         }
         if (!strcmp(ACptrs->Type,"B")) strcpy(ACptrs->Type,"BT");
         ACptrs->Reg=AddElemToList(ACptrs->Reg,Eptr);
         Eptr->Cont=ACptrs;
      }
      else if (!strncmp(Line,"RP ",3) || !strncmp(Line,"RM ",3)){
         Eptr->Tmax=GetValue(Line,46,5,2)*K3;
         Eptr->Tmin=GetValue(Line,51,5,2)*K3;
         if(Eptr->Tmax<=Eptr->Tmin) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("Phase shifter limits are wrong: Amin >= Amax.");
         }
         if (!strncmp(Line,"RP ",3)) {
            strcpy(Eptr->Type,"RP");
            Eptr->Cvar=GetValue(Line,58,5,0)/Sn;
         } else {
            strcpy(Eptr->Type,"RM");
            Eptr->Max=GetValue(Line,58,5,0)/Sn;
            Eptr->Min=GetValue(Line,63,5,0)/Sn;
            if(Eptr->Max<=Eptr->Min) {
               fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
               ErrorHalt("Phase shifter limits are wrong: Pmin >= Pmax.");
            }
         }
         if (ACptrs!=ACptr && ACptrs!=ACptrp) ACptrs=ACptr;
         Eptr->Cont=ACptrs;
         Eptr->Ncont=ACptrs->Ncont;
         ACptrs->Reg=AddElemToList(ACptrs->Reg,Eptr);
         ACptrs->Ncont++;
      }
      else if (!strncmp(Line,"RQ ",3) || !strncmp(Line,"RN ",3)){
         Eptr->Tmax=GetValue(Line,46,5,2)/KV;
         Eptr->Tmin=GetValue(Line,51,5,2)/KV;
         if(Eptr->Tmax<0) Eptr->Tmax=1.1;
         if(Eptr->Tmin<0) Eptr->Tmin=0.9;
         if(Eptr->Tmax<=Eptr->Tmin) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("LTC limits are wrong: Tmin >= Tmax.");
         }
         if (!strncmp(Line,"RQ ",3)){
            strcpy(Eptr->Type,"RQ");
            Eptr->Cvar=GetValue(Line,58,5,0)/Sn;
         } else {
            strcpy(Eptr->Type,"RN");
            Eptr->Max=GetValue(Line,58,5,0)/Sn;
            Eptr->Min=GetValue(Line,63,5,0)/Sn;
            if(Eptr->Max<=Eptr->Min) {
               fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
               ErrorHalt("LTC limits are wrong: Qmin >= Qmax.");
            }
         }
         if (ACptrs!=ACptr && ACptrs!=ACptrp) ACptrs=ACptr;
         Eptr->Cont=ACptrs;
         Eptr->Ncont=ACptrs->Ncont;
         ACptrs->Reg=AddElemToList(ACptrs->Reg,Eptr);
         ACptrs->Ncont++;
      }
    }

  /* -------------------- DC data -------------------------------- */
    else if (!strncmp(Line,"BD ",3) || !strncmp(Line,"BZ ",3)||
        !strncmp(Line,"LD ",3)) ReadEPRIdc(Line);


  /* ------------------------ Area data -------------------------- */
    else if (!strncmp(Line,"A  ",3)) {
      GetStr(Line,4,10,30,Name);
      Aptr=(AreaData *) AreaInList(0,Name,Narea);
      if (Aptr->N==0) { Narea++; Aptr->N=Narea;}
      else {
         fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
         ErrorHalt("The Area was previously defined (check A cards).");
      }
      Aptr->P=GetValue(Line,27,8,0)/Sn;
      GetStr(Line,14,12,12,Name);
      KV=GetValue(Line,22,4,0);
      ACptr=ACbusInList(0,Name,KV,Nac,1);
      if (ACptr->N==0) { Nac++; ACptr->Num=ACptr->N=Nac;}
      ACptr->Area=Aptr;
      Aptr->Slack=Aptr->BSptr=ACptr;
      for (i=36, j=1; i<=63; i=i+3, j++) GetStr(Line,i,2,2,Aptr->Zone[j]);
    }

  /* ---------------------- Solution data ------------------------ */
    else if (!strncmp(Line,"SOL",3)) {
      GetStr(Line,31,12,12,Name);
      KV=GetValue(Line,39,4,0);
      AngSlack=GetValue(Line,46,10,4);
      AngSlack=AngSlack*K3;
      ACptr=ACbusInList(0,Name,KV,Nac,1);
      if (ACptr->N==0) { Nac++; ACptr->Num=ACptr->N=Nac;}
      if(!strpbrk(ACptr->Type,"S")){
         strcat_s(ACptr->Type,"S");
         Nslack++;
      }
      MaxIter=GetInt(Line,24,5);
    }
                         /* FACTS */

 /* ---------------------- SVC data ------------------------ */
    else if (!strncmp(Line,"FS ",3)) ReadSVC(Line);

 /* ---------------------- TCSC data ------------------------ */
    else if (!strncmp(Line,"FC ",3)) ReadTCSC(Line);

 /* ---------------------- STATCOM data ------------------------ */
    else if (!strncmp(Line,"FT ",3)) ReadSTATCOM(Line);

                      /* END FACTS */

    else if (!strncmp(Line,"END",3)||!strncmp(Line,"FIM",3)) break;
    else if (strncmp(Line,"ZZ",2)&&strncmp(Line,"D",1)&&strncmp(Line,"9999",4)) {
      fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
      fprintf(stderr,"***Warning: The program will ignore this line.\n");
    }
  }
  fclose(InputDataFile);
  NacEl-=NregPQ;
  NregPQ-=NregV;
  if (MaxIter==0) MaxIter=50;
}
