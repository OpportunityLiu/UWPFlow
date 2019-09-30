/* Write output in IEEE format. */

#include "write.h"
#include <time.h>
#include <string.h>  /* FACTS */


/* --------------- Print ----------------- */
#ifdef ANSIPROTO
void Print(FILE *File,int spaces,int width,int decimals,VALUETYPE val)
#else
void Print(File,spaces,width,decimals,val)
FILE *File;
int spaces,width,decimals;
VALUETYPE val;
#endif
{
  int i,j;
  double cons;
  BOOLEAN flagNegative=FALSE;
  char str[50],*ptr;

  if (spaces>0) fprintf(File,"%*s",spaces,"");
  width=abs(width); decimals=abs(decimals);
  if(val<0) {
    width--;
    val=fabs(val);
    flagNegative=TRUE;
  }
  if(width<=decimals) decimals=width-1;
  if(decimals) width--;
  if(!width) {width=1; decimals=0;}
  if(decimals) cons=5*pow(10.,(double)-width-1);
  else cons=5*pow(10.,(double)-width);
  if (val<cons) {
    if(flagNegative && width>1) width++;
    if(decimals) width++;
    if (width==1) strcpy_s(str,"0");
    else {
      strcpy_s(str,"0.");
      for(i=2;i<=width-1;i++) strcat_s(str,"0");
    }
    fprintf(File,"%*s",width,str);
    return;
  }
  if(width==1 && flagNegative) {
    fprintf(File,"*");
    return;
  }
  if(val<(pow(10.,(double)width)-0.5)){
    if(flagNegative) fprintf(File,"-");
    j=0;
    for(i=width-1;i>=0;i--) {
      if(val>(pow(10.,(double)i)-0.5*pow(10.,(double)-j))) break;
      j++;
    }
    if(decimals) width++;
    else if(i<width-1) j--;
    sprintf_s(str,"%-*.*lf",width,j,val);
    if(!strncmp(str,"0.",2)) ptr=strpbrk(str,".");
    else if(strlen(str)==width && strpbrk(str," ")) {
      str[width-1]='.';
      str[width]='\0';
      ptr=str;
    }
    else ptr=str;
    fprintf(File,"%*s",width,ptr);
    return;
  }
  else {
    if(flagNegative) width++;
    if(decimals) width++;
    for(i=0;i<=width-2;i++) {
      if (flagNegative && i==0) str[i]='-';
      else str[i]='9';
    }
    if(decimals) str[width-1]='.';
    else         str[width-1]='9';
    str[width]='\0';
    fprintf(File,"%s",str);
    return;
  }
}


/* --------------- IEEE ----------------- */
#ifdef ANSIPROTO
void IEEE(void)
#else
void IEEE()
#endif
{
  ACbusData *ACptr;
  DCbusData *DCptrR,*DCptrI,*DCptr;
  SVCbusData *SVCptr;           /* FACTS */
  TCSCbusData *TCSCptr;         /* FACTS */
  STATCOMbusData *STATCOMptr;   /* FACTS */
  ElementData *Eptr;
  ElementList *ELptr;
  AreaData *Aptr;
  char str[32],Num[5],Nump[5],Area[3],Zone[4];
  int Type,Nties,i,j;
  VALUETYPE R,X,Vn,In,Pl,Ql,Imax,delta,vals;
  VALUETYPE Ssvc,Stcsc,Sstatcom,G,B,Max,Xc; /* FACTS */
  BOOLEAN card=FALSE,flag=FALSE;
  FILE *OutFile;
  time_t t;
  struct tm *localt;
  int month,year,yearp;
  Type = 0;


  if (ExistParameter('W')) {
    if (NullName(NameParameter('W'))) return;
    OutFile=OpenOutput(NameParameter('W'));
  } else {
    if (NullName(NameParameter('w'))) return;
    OutFile=OpenOutput(NameParameter('w'));
    card=TRUE;
  }
  if(card) fprintf(OutFile,"CARD\n");
  else fprintf(OutFile,"TAPE\n");
  t = time(NULL);
  localt=localtime(&t);
  if (localt->tm_mday<10) fprintf(OutFile," 0%1d/",localt->tm_mday);
  else                    fprintf(OutFile," %2d/",localt->tm_mday);
  month=localt->tm_mon+1;
  if (month<10) fprintf(OutFile,"0%1d/",month);
  else          fprintf(OutFile,"%2d/",month);
  year=localt->tm_year;
  yearp=year+1900;
  if (year>=100) year=year-100;
  if (year<10) fprintf(OutFile,"0%1d",year);
  else         fprintf(OutFile,"%2d",year);
  fprintf(OutFile," Generated with PFLOW");
  Print(OutFile,1,6,0,Sn);
  fprintf(OutFile," %4d",yearp);
  if (month>3 && month<10) fprintf(OutFile," S ");
  else                     fprintf(OutFile," W ");
  strncpy_s(str,dataPtr->Title[0],29);
  for(i=0;i<=28;i++){
    if (str[i]=='\n') str[i]='\0';
  }
  str[28]='\0';
  fprintf(OutFile,"%s\n",str);

 /* --------------------- AC bus results -----------------------------*/
  fprintf(OutFile,"BUS DATA FOLLOWS                  %5d ITEMS\n",Nac);
  for (ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
    if (ACptr->Num<=9999) sprintf_s(Num,"%4d",ACptr->Num);
    else strcpy_s(Num,"****");
    if (Narea<2 || ACptr->Area==NULL) strcpy_s(Area,"0");
    else if (ACptr->Area->N<=99) sprintf_s(Area,"%2d",ACptr->Area->N);
    else                         strcpy_s(Area,"99");
    if ((isdigit(ACptr->Zone[0])||ACptr->Zone[0]==' ') && (isdigit(ACptr->Zone[1])||ACptr->Zone[1]==' '))
         sprintf_s(Zone,"%3d",atoi(ACptr->Zone));
    else sprintf_s(Zone,"%3d",toascii(ACptr->Zone[0])+toascii(ACptr->Zone[1]));
    if(!strcmp(ACptr->Type,"B")||!strcmp(ACptr->Type,"BA")||strpbrk(ACptr->Type,"L,T,C,R")) Type=0;
    else if(strpbrk(ACptr->Type,"V,X")) Type=1;
    else {
      if(strpbrk(ACptr->Type,"G,Q,Z")) Type=2;
      if(strpbrk(ACptr->Type,"S")) Type=3;
    }
    fprintf(OutFile,"%4s %12s %2s%3s %2d",Num,ACptr->Name,Area,Zone,Type);
    Print(OutFile,1,6,4,ACptr->V);
    delta=ACptr->Ang;
    if (delta>=0) vals=1.00;
    else          vals=-1.00;
    if (fabs(delta)>2*PI) delta=delta-vals*floor(fabs(delta)/(2*PI))*2*PI;
    if (fabs(delta)>PI) delta=delta-vals*2*PI;
    ACptr->Ang=delta;
    Print(OutFile,1,6,2,ACptr->Ang/K3);
    Pl=(ACptr->Pn+lambda*ACptr->Pnl)*pow(ACptr->V,ACptr->a)+
       (ACptr->Pz+lambda*ACptr->Pzl)*ACptr->V*ACptr->V;
    Print(OutFile,1,8,2,Pl*Sn);
    Ql=(ACptr->Qn+lambda*ACptr->Qnl)*pow(ACptr->V,ACptr->b)+
       (ACptr->Qz+lambda*ACptr->Qzl)*ACptr->V*ACptr->V;
    Print(OutFile,1,8,2,Ql*Sn);
    Print(OutFile,1,8,2,ACptr->PG*Sn);
    Print(OutFile,1,7,2,ACptr->Qg*Sn);
    if(card) fprintf(OutFile,"\n");
    Print(OutFile,1,7,2,ACptr->KV);
    if(strpbrk(ACptr->Type,"G") && ACptr->Cont!=NULL) Print(OutFile,1,6,4,ACptr->Cont->VCont);
    else if(strpbrk(ACptr->Type,"Q,S")) Print(OutFile,1,6,4,ACptr->VCont);
    else Print(OutFile,1,6,4,ACptr->V);
    if(strpbrk(ACptr->Type,"V,X")) {
      Print(OutFile,1,7,2,ACptr->Vmax);
      Print(OutFile,1,7,2,ACptr->Vmin);
    }
    else {
      Print(OutFile,1,7,2,ACptr->Qmax*Sn);
      Print(OutFile,1,7,2,ACptr->Qmin*Sn);
    }
    for(ELptr=ACptr->Elem;ELptr!=NULL;ELptr=ELptr->Next){
      Eptr=ELptr->Eptr;
      if(Eptr->From==ACptr){
        ACptr->G=ACptr->G+Eptr->G1*Eptr->Tap*Eptr->Tap;
        if(Eptr->B1!=Eptr->B2) ACptr->B=ACptr->B+Eptr->B1*Eptr->Tap*Eptr->Tap;
      }
      else {
        ACptr->G=ACptr->G+Eptr->G2;
        if(Eptr->B1!=Eptr->B2) ACptr->B=ACptr->B+Eptr->B2;
      }
    }
    Print(OutFile,1,7,4,ACptr->G);
    Print(OutFile,1,7,4,ACptr->B);
    if(strpbrk(ACptr->Type,"G") && ACptr->Cont!=NULL) {
      if (ACptr->Cont->Num<=9999) sprintf_s(Num,"%4d",ACptr->Cont->Num);
      else strcpy_s(Num,"****");
    }
         else strcpy_s(Num,"   0");
    fprintf(OutFile," %4s\n",Num);
  }
  fprintf(OutFile,"-999\n");

 /* --------------------- AC element results ---------------------------*/
  Nties=0;
  fprintf(OutFile,"BRANCH DATA FOLLOWS               %5d ITEMS\n",NacEl);
  for(Eptr=dataPtr->Element;Eptr!=NULL;Eptr=Eptr->Next){
    if(Narea>1 && Eptr->From->Area!=Eptr->To->Area) Nties++;
    if (Eptr->From->Num<=9999) sprintf_s(Num,"%4d",Eptr->From->Num);
    else strcpy_s(Num,"****");
    if (Eptr->To->Num<=9999) sprintf_s(Nump,"%4d",Eptr->To->Num);
    else strcpy_s(Nump,"****");
    if (Eptr->Area==NULL) {
      if (Eptr->From->Area==Eptr->To->Area) Eptr->Area=Eptr->From->Area;
      else Eptr->Area=Eptr->Meter->Area;
    }
    if (Narea<2 || Eptr->Area==NULL) strcpy_s(Area,"0");
    else if (Eptr->Area->N<=99) sprintf_s(Area,"%2d",Eptr->Area->N);
    else strcpy_s(Area,"**");
    if ((isdigit(Eptr->Zone[0])||Eptr->Zone[0]==' ') && (isdigit(Eptr->Zone[1])||Eptr->Zone[1]==' '))
         sprintf_s(Zone,"%3d",atoi(Eptr->Zone));
    else sprintf_s(Zone,"%3d",toascii(Eptr->Zone[0])+toascii(Eptr->Zone[1]));
    if(strpbrk(Eptr->Type,"LE")) Type=0;
    else if(strpbrk(Eptr->Type,"T")) Type=1;
    else if(!strcmp(Eptr->Type,"R")||strpbrk(Eptr->Type,"V")) Type=2;
    else if(strpbrk(Eptr->Type,"QN")) Type=3;
    else Type=4;
    if(!strcmp(Eptr->Ckt," ")) strcpy_s(Eptr->Ckt,"1");
    fprintf(OutFile,"%4s %4s %2s%3s %1s %1d",Num,Nump,Area,Zone,Eptr->Ckt,Type);
    R=Eptr->G/(Eptr->G*Eptr->G+Eptr->B*Eptr->B);
    X= -Eptr->B/(Eptr->G*Eptr->G+Eptr->B*Eptr->B);
    Print(OutFile,1,9,6,R);
    Print(OutFile,1,9,6,X);
    if (Eptr->B1==Eptr->B2) Print(OutFile,1,9,5,Eptr->B1+Eptr->B2);
    else Print(OutFile,1,9,5,0.);
    Imax=ceil(Eptr->Imax*Sn);
    if (Imax>99999) fprintf(OutFile,"99999");
    else fprintf(OutFile," %5.0lf",Imax);
    fprintf(OutFile,"%12s","");
    if(strpbrk(Eptr->Type,"R")) {
      if (Eptr->Cont->Num<=9999) sprintf_s(Num,"%4d",Eptr->Cont->Num);
      else strcpy_s(Num,"****");
      fprintf(OutFile," %4s",Num);
      if (strpbrk(Eptr->Type,"PQNM")) fprintf(OutFile," 0 ");
      else {
        if(Eptr->Cont==Eptr->From||Eptr->Cont==Eptr->To) fprintf(OutFile," 0 ");
        else fprintf(OutFile," 1 ");
      }
    } else fprintf(OutFile,"    0 0 ");
    if(strpbrk(Eptr->Type,"RT")) {
      if(card) fprintf(OutFile,"\n");
      if(Eptr->Tap) Print(OutFile,1,6,4,1/Eptr->Tap);
      else Print(OutFile,1,6,4,0.);
      Print(OutFile,1,7,2,Eptr->Ang/K3);
      if(!strcmp(Eptr->Type,"RP") || !strcmp(Eptr->Type,"RM")) {
        Print(OutFile,1,6,2,Eptr->Tmin/K3);
        Print(OutFile,1,6,2,Eptr->Tmax/K3);
      }
      else if(strpbrk(Eptr->Type,"R")) {
        Print(OutFile,1,6,2,Eptr->Tmin);
        Print(OutFile,1,6,2,Eptr->Tmax);
      }
      else {
        Print(OutFile,1,6,2,0.);
        Print(OutFile,1,6,2,0.);
      }
      Print(OutFile,1,6,5,0.);
      if(!strcmp(Eptr->Type,"RP") || !strcmp(Eptr->Type,"RQ")) {
        Print(OutFile,1,6,5,Eptr->Cvar*Sn);
        Print(OutFile,1,6,5,Eptr->Cvar*Sn);
      }
      else if(strpbrk(Eptr->Type,"NM")) {
        Print(OutFile,1,6,5,Eptr->Min*Sn);
        Print(OutFile,1,6,5,Eptr->Max*Sn);
      }
      else if(!strcmp(Eptr->Type,"R")) {
        Print(OutFile,1,6,5,Eptr->Cont->V);
        Print(OutFile,1,6,5,Eptr->Cont->V);
      }
      else if(strpbrk(Eptr->Type,"V")) {
        Print(OutFile,1,6,5,Eptr->Cont->Vmin);
        Print(OutFile,1,6,5,Eptr->Cont->Vmax);
      }
      else {
        Print(OutFile,1,6,5,0.);
        Print(OutFile,1,6,5,0.);
      }
    } else {
	  Print(OutFile,1,6,4,0.);
      Print(OutFile,1,7,2,0.);
      Print(OutFile,1,6,2,0.);
      Print(OutFile,1,6,2,0.);
      Print(OutFile,1,6,5,0.);
      Print(OutFile,1,6,5,0.);
      Print(OutFile,1,6,5,0.);
	}
    fprintf(OutFile,"\n");
  }
  fprintf(OutFile,"-999\n");

 /* ------------------------- Area results ---------------------------*/
  fprintf(OutFile,"LOSS ZONES FOLLOW                 %5d ITEMS\n",0);
  fprintf(OutFile,"-99\n");
  if(Narea>1) fprintf(OutFile,"INTERCHANGE DATA FOLLOWS          %5d ITEMS\n",Narea);
  else fprintf(OutFile,"INTERCHANGE DATA FOLLOWS          %5d ITEMS\n",0);
  if(Narea>1) for(Aptr=dataPtr->Area;Aptr!=NULL;Aptr=Aptr->Next){
    if (Aptr->N<=99) sprintf_s(Area,"%2d",Aptr->N);
    else strcpy_s(Area,"**");
    if (Aptr->BSptr->Num<=9999) sprintf_s(Num,"%4d",Aptr->BSptr->Num);
    else strcpy_s(Num,"****");
    fprintf(OutFile,"%2s %4s %12s",Area,Num,Aptr->BSptr->Name);
    Print(OutFile,1,7,1,Aptr->P*Sn);
    fprintf(OutFile,"%15s  %30s\n","",Aptr->Name);
  }
  fprintf(OutFile,"-9\n");
  fprintf(OutFile,"TIE LINES FOLLOW                  %5d ITEMS\n",Nties);
  if(Nties) for(Eptr=dataPtr->Element;Eptr!=NULL;Eptr=Eptr->Next){
    if(Eptr->From->Area!=Eptr->To->Area) {
      ACptr=Eptr->Meter;
      if (ACptr->Num<=9999) sprintf_s(Num,"%4d",ACptr->Num);
      else strcpy_s(Num,"****");
      if (ACptr->Area->N<=99) sprintf_s(Area,"%2d",ACptr->Area->N);
      else strcpy_s(Area,"**");
      fprintf(OutFile,"%4s  %2s",Num,Area);
      if (Eptr->Meter==Eptr->From) ACptr=Eptr->To;
      else ACptr=Eptr->From;
      if (ACptr->Num<=9999) sprintf_s(Num,"%4d",ACptr->Num);
      else strcpy_s(Num,"****");
      if (ACptr->Area->N<=99) sprintf_s(Area,"%2d",ACptr->Area->N);
      else strcpy_s(Area,"**");
      fprintf(OutFile,"  %4s  %2s  %1s\n",Num,Area,Eptr->Ckt);
    }
  }
  fprintf(OutFile,"-999\n");

 /* --------------------- DC bus results -----------------------------*/
  for(DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if(!strcmp(DCptrR->Type,"R")){
      Vn=DCptrR->Vn; In=Sn/Vn;
      DCptrR->Xc=DCptrR->Xc*Vn/In/DCptrR->Nbr;
      DCptrI->Xc=DCptrI->Xc*Vn/In/DCptrI->Nbr;
      DCptrI->Ntrf=DCptrI->Ntrf*Vn/DCptrI->AC->KV;
      DCptrR->Ntrf=DCptrR->Ntrf/DCptrR->Nbr;
      DCptrI->Ntrf=DCptrI->Ntrf/DCptrI->Nbr;
      for(i=1;i<=2;i++){
        if(i==1) DCptr=DCptrR;
        else DCptr=DCptrI;
        fprintf(OutFile,"BD %8sGROUND  %12s",DCptr->Name,DCptr->AC->Name);
        fprintf(OutFile,"    %2s%2.0lf",DCptr->Zone,DCptr->Nbr);
        Print(OutFile,0,5,3,DCptr->Xc);
        Print(OutFile,0,5,4,DCptr->Ntrf);
        Print(OutFile,0,5,3,0.);
        if(DCptr->Tap>=DCptr->TapMin) Print(OutFile,0,5,3,DCptr->TapMin);
        else {
          Print(OutFile,0,5,3,(floor(DCptr->Tap*1000.)*.001));
          flag=TRUE;
        }
        if(DCptr->Tap<=DCptr->TapMax) Print(OutFile,0,5,3,DCptr->TapMax);
        else {
          Print(OutFile,0,5,3,(ceil(DCptr->Tap*1000.)*.001));
                         flag=TRUE;
        }
        if(flag){ 
          flag=FALSE;
          if (!strcmp(DCptr->Cont1,"AL")||!strcmp(DCptr->Cont1,"GA")) strcpy_s(DCptr->Cont1,"AT");
          else if (!strcmp(DCptr->Cont2,"AL")||!strcmp(DCptr->Cont2,"GA")) strcpy_s(DCptr->Cont2,"AT");
          else strcpy_s(DCptr->Cont1,"AT");
        }
        Print(OutFile,0,4,1,DCptr->AlfaMin/K3);
        Print(OutFile,0,4,1,DCptr->AlfaMax/K3);
        Print(OutFile,0,4,1,DCptr->GammaMin/K3);
        fprintf(OutFile,"\nBZ %8sGROUND  %12s",DCptr->Name,DCptr->AC->Name);
        fprintf(OutFile," %2s%2s",DCptr->Cont1,DCptr->Cont2);
        for(j=1;j<=2;j++){
          if(j==1) strcpy_s(str,DCptr->Cont1);
          else strcpy_s(str,DCptr->Cont2);
          if(!strcmp(str,"VD")) Print(OutFile,0,6,2,DCptr->Vd*Vn);
          if(!strcmp(str,"ID")) Print(OutFile,0,6,2,DCptr->Id*In*1000);
          if(!strcmp(str,"PA")) Print(OutFile,0,6,2,fabs(DCptr->P)*Sn);
          if(!strcmp(str,"QA")) Print(OutFile,0,6,2,fabs(DCptr->Q)*Sn);
          if(!strcmp(str,"AT")) Print(OutFile,0,6,2,DCptr->Tap);
          if(!strcmp(str,"AL")) Print(OutFile,0,6,2,DCptr->Alfa/K3);
          if(!strcmp(str,"GA")) Print(OutFile,0,6,2,DCptr->Gamma/K3);
        }
        fprintf(OutFile,"%1s\n",DCptr->Type);
      }
      fprintf(OutFile,"LD %8s%8s",DCptrR->Name,DCptrI->Name);
      fprintf(OutFile,"    %2s",DCptrR->Lzone);
      Print(OutFile,3,6,2,DCptrR->Rd*Vn/In);
      Print(OutFile,0,6,2,DCptrR->Ld*Vn/In*1000);
      fprintf(OutFile,"\n");
    }
  }

                        /* FACTS */
/* --------------------- SVC results ----------------------------- */
  for(SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
     fprintf(OutFile,"%2s",SVCptr->Type);
     if(SVCptr->From->Area!=NULL) fprintf(OutFile," %2d ",SVCptr->From->Area->N);
     else fprintf(OutFile," %2d ",0);
     fprintf(OutFile,"%12s ",SVCptr->From->Name);
     fprintf(OutFile,"%12s",SVCptr->Ctrl->Name);
     Vn=SVCptr->Ctrl->KV;
     Ssvc=SVCptr->SVC_base;
     Print(OutFile,1,7,6,SVCptr->Xc*Ssvc/Sn);
     Print(OutFile,1,7,6,SVCptr->Xl*Ssvc/Sn);
     Print(OutFile,0,3,0,SVCptr->AlphaMin/K3);
     Print(OutFile,0,3,0,SVCptr->AlphaMax/K3);
     Print(OutFile,0,2,0,SVCptr->slope*100.0*Ssvc/Sn);
     Print(OutFile,0,4,1,SVCptr->SVC_base);
     Print(OutFile,1,3,2,SVCptr->Vref);
     fprintf(OutFile,"%12s","            ");
     Print(OutFile,1,4,1,SVCptr->alpha_svc/K3);
     fprintf(OutFile,"\n");
  }
/* --------------------- TCSC results ----------------------------- */
  for(TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;TCSCptr=TCSCptr->Next){
     fprintf(OutFile,"%2s",TCSCptr->Type);
     if(TCSCptr->From->Area!=NULL) fprintf(OutFile," %2d ",TCSCptr->From->Area->N);
     else fprintf(OutFile," %2d ",0);
     fprintf(OutFile,"%12s ",TCSCptr->From->Name);
     fprintf(OutFile,"%12s",TCSCptr->To->Name);
     Vn=TCSCptr->From->KV;
     Stcsc=TCSCptr->TCSC_base;
     Xc=TCSCptr->Xc;
     Max=TCSCptr->Max;
     Print(OutFile,1,7,6,TCSCptr->Xc*Stcsc/Sn);
     Print(OutFile,1,7,6,TCSCptr->Xl*Stcsc/Sn);
     Print(OutFile,0,3,0,TCSCptr->AlphaMin/K3);
     Print(OutFile,0,3,0,TCSCptr->AlphaMax/K3);
     if (!strcmp(TCSCptr->Cont,"X")) Print(OutFile,1,7,3,100.0/(TCSCptr->Bset*Max*Xc));
     else if (!strcmp(TCSCptr->Cont,"P")) Print(OutFile,1,7,3,TCSCptr->Control*Sn);
     else if (!strcmp(TCSCptr->Cont,"I")) Print(OutFile,1,7,3,TCSCptr->Control*Sn/Vn);
     else if (!strcmp(TCSCptr->Cont,"D")) Print(OutFile,1,7,3,TCSCptr->Control/K3);
     fprintf(OutFile,"%c",TCSCptr->Cont[0]);
     Print(OutFile,0,4,1,Stcsc);
     Print(OutFile,0,5,1,TCSCptr->alpha_tcsc/K3);
     Print(OutFile,0,3,1,TCSCptr->Max);
     fprintf(OutFile,"\n");
  }
/* --------------------- STATCOM results ----------------------------- */
  for(STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
     fprintf(OutFile,"%2s",STATCOMptr->Type);
     if(STATCOMptr->From->Area!=NULL) fprintf(OutFile," %2d ",STATCOMptr->From->Area->N);
     else fprintf(OutFile," %2d ",0);
     fprintf(OutFile,"%12s ",STATCOMptr->From->Name);
     fprintf(OutFile,"%12s",STATCOMptr->Ctrl->Name);
     Vn=STATCOMptr->From->KV;
     Sstatcom=STATCOMptr->MVA;
     G=STATCOMptr->G;
     B=STATCOMptr->B;
     X=-B/(G*G+B*B);
     Print(OutFile,1,5,4,STATCOMptr->R*Sstatcom/Sn);
     Print(OutFile,1,5,4,X*Sstatcom/Sn);
     Print(OutFile,1,5,4,STATCOMptr->Gc*Sn/Sstatcom);
     Print(OutFile,1,5,3,STATCOMptr->Imin*Sn/Sstatcom);
     Print(OutFile,1,5,3,STATCOMptr->Imax*Sn/Sstatcom);
     Print(OutFile,1,2,0,STATCOMptr->slope*100.0*Sstatcom/Sn);
     Print(OutFile,1,4,0,STATCOMptr->MVA);
     Print(OutFile,1,3,2,STATCOMptr->Vref);
     if (!strcmp(STATCOMptr->Cont1,"PW")) {
       fprintf(OutFile,"%1s","W");
       Print(OutFile,1,5,4,STATCOMptr->Vdc);
     } else {
       fprintf(OutFile,"%1s","P");
       Print(OutFile,1,5,4,STATCOMptr->k);
     }
     fprintf(OutFile,"\n");
  }
                     /* END OF FACTS */

  fprintf(OutFile,"END OF DATA\n");
  fclose(OutFile);
}
