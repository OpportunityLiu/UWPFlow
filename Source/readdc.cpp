/* Read DC data in WSCC/EPRI/BPA format   */

#include "readdata.h"

/* --------- Global Input File --------- */
extern FILE *InputDataFile;

/* ---------------- ReadEPRIdc ----------------------------- */
#ifdef ANSIPROTO
void ReadEPRIdc(const char *Line)
#else
void ReadEPRIdc(Line)
char *Line;
#endif
/* Read DC Bus and Element data in ETMSP format. */
{
  ACbusData *ACptr;
  DClist *ptr;
  DCbusData *DCptr,*DCptrp;
  char Name[13],Mode[2];
  VALUETYPE KV,KVp,Set1,Set2,R,L,P,V,alpha,gamma;
  BOOLEAN flag=FALSE,flagEPRI=TRUE;
  int i;

  /* --------------- DC bus data -------------------------------- */
    if (!strncmp(Line,"BD ",3)) {
      GetStr(Line,51,8,8,Name);
      for (i=0;i<8;i++) if (isalpha(Name[i])) {
        flagEPRI=FALSE;
        break;
      }

    /* ------------ EPRI Format Multiterminal Format ------------- */
      if (flagEPRI) {
        flag=FALSE;
        for(i=4;i<=12;i+=8) {
          GetStr(Line,i,8,8,Name);
          if(strncmp(Name,"GROUND",6)){
            DCptr=(DCbusData *) DCbusInList(Name,Ndc);
            if (DCptr->Nbr!=0) {
              fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
              ErrorHalt("The DC bus was previously defined (check BD cards).");
            }
            if (DCptr->N==0) { Ndc++; DCptr->N=Ndc;}
          } else flag=TRUE;
        }
        if (!flag){
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("One of the DC buses must be GROUND (2 terminal HVDC links).");
        }
        GetStr(Line,20,12,12,Name);
        KV=GetValue(Line,28,4,0);
        ACptr=ACbusInList(0,Name,KV,Nac,1);
        if (ACptr->N==0) { Nac++; ACptr->Num=ACptr->N=Nac;}
        DCptr->Vn=KV;
        ptr=ACptr->DC;
#ifdef WINDOWS
        ACptr->DC= new DClist;
#else
        ACptr->DC=(DClist *) malloc(sizeof(DClist));
        if(ACptr->DC==NULL) {
          fclose(InputDataFile);
          ErrorHalt("Insufficient memory to allocate DC elemet data");
          exit(ERROREXIT);
        }
#endif
        ACptr->DC->DC=DCptr;
        ACptr->DC->Next=ptr;
        DCptr->AC=ACptr;
        GetStr(Line,36,2,2,DCptr->Zone);
        DCptr->Nbr=GetValue(Line,38,2,0);
        DCptr->Xc=GetValue(Line,40,5,3);
        DCptr->Ntrf=GetValue(Line,45,5,4);
        if (DCptr->Nbr<1 || DCptr->Ntrf<0.0001) {
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("Wrong number of AC/DC bridges or transf. ratio too small.");
        }
        DCptr->Xc=DCptr->Xc*DCptr->Nbr;
        DCptr->Ntrf=DCptr->Ntrf*DCptr->Nbr;
        DCptr->TapMin=GetValue(Line,55,5,3);
        if (DCptr->TapMin<=0) DCptr->TapMin=0.9;
        DCptr->TapMax=GetValue(Line,60,5,3);
        if (DCptr->TapMax<=DCptr->TapMin) DCptr->TapMax=1.1;
        DCptr->AlfaMin=GetValue(Line,65,4,1);
        if (DCptr->AlfaMin<0) DCptr->AlfaMin=0.;
        DCptr->AlfaMax=GetValue(Line,69,4,1);
        if (DCptr->AlfaMax<=DCptr->AlfaMin) DCptr->AlfaMax=90.;
        DCptr->GammaMin=GetValue(Line,73,4,1);
        if (DCptr->GammaMin<0) DCptr->GammaMin=0.;
      }

    /* ---------------- WSC/BPA 2 Terminal Format ---------------- */
      else {
        GetStr(Line,7,12,12,Name);
        KVp=GetValue(Line,15,4,0);
        ACptr=ACbusInList(0,Name,KVp,Nac,1);
        if (ACptr->N==0) {
          Nac++;
          ACptr->Num=ACptr->N=Nac;
          ACptr->V=1.;
          GetStr(Line,19,2,2,ACptr->Zone);
          GetStr(Line,4,3,3,ACptr->Owner);
        }
        DCptr=(DCbusData *) DCbusInList(Name,Ndc);
        if (DCptr->Nbr!=0) {
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("The DC bus was previously defined (check BD cards).");
        }
        if (DCptr->N==0) { Ndc++; DCptr->N=Ndc;}
        GetStr(Line,19,2,2,DCptr->Zone);
        DCptr->Nbr=GetValue(Line,24,2,0);
        GetStr(Line,51,12,12,Name);
        KV=GetValue(Line,59,4,0);
        DCptr->Ntrf=(KVp/KV)*DCptr->Nbr;
        if (DCptr->Nbr<1 || DCptr->Ntrf<0.0001) {
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("Wrong number of AC/DC bridges or transf. ratio too small.");
        }
        ACptr=ACbusInList(0,Name,KV,Nac,1);
        if (ACptr->N==0) { Nac++; ACptr->Num=ACptr->N=Nac;}
        DCptr->Vn=KV;
        ptr=ACptr->DC;
#ifdef WINDOWS
        ACptr->DC= new DClist;
#else
        ACptr->DC=(DClist *) malloc(sizeof(DClist));
        if(ACptr->DC==NULL) {
          fclose(InputDataFile);
          ErrorHalt("Insufficient memory to allocate DC elemet data");
          exit(ERROREXIT);
        }
#endif
        ACptr->DC->DC=DCptr;
        ACptr->DC->Next=ptr;
        DCptr->AC=ACptr;
        DCptr->Ld+=GetValue(Line,26,5,1)/1000.;
        DCptr->AlfaMin=GetValue(Line,31,5,1);
        if (DCptr->AlfaMin<0) DCptr->AlfaMin=0.;
        DCptr->GammaMin=GetValue(Line,36,5,1);
        if (DCptr->GammaMin<0) DCptr->GammaMin=0.;
      }
    }

    /* ------------ EPRI Format Multiterminal Format ------------- */
    else if (!strncmp(Line,"BZ ",3)) {
      flag=FALSE;
      for(i=4;i<=12;i+=8) {
        GetStr(Line,i,8,8,Name);
        if(strncmp(Name,"GROUND",6)){
          DCptr=(DCbusData *) DCbusInList(Name,Ndc);
          if (strcmp(DCptr->Type,"")) {
            fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
            ErrorHalt("The DC bus was previously defined (check BD cards).");
          }
          if (DCptr->N==0) { Ndc++; DCptr->N=Ndc;}
        } else flag=TRUE;
      }
      if (!flag){
        fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
        ErrorHalt("One of the DC buses must be GROUND (2 terminal HVDC links).");
      }
      GetStr(Line,33,2,2,DCptr->Cont1);
      GetStr(Line,35,2,2,DCptr->Cont2);
      GetStr(Line,49,1,1,DCptr->Type);
      if ((strcmp(DCptr->Cont1,"ID") && strcmp(DCptr->Cont1,"VD") &&
           strcmp(DCptr->Cont1,"PA") && strcmp(DCptr->Cont1,"QA") &&
           strcmp(DCptr->Cont1,"AL") && strcmp(DCptr->Cont1,"GA") &&
           strcmp(DCptr->Cont1,"AT"))||
          (strcmp(DCptr->Cont2,"ID") && strcmp(DCptr->Cont2,"VD") &&
           strcmp(DCptr->Cont2,"PA") && strcmp(DCptr->Cont2,"QA") &&
           strcmp(DCptr->Cont2,"AL") && strcmp(DCptr->Cont2,"GA") &&
           strcmp(DCptr->Cont2,"AT"))||
          (strcmp(DCptr->Type,"R") && strcmp(DCptr->Type,"I"))){
        fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
        ErrorHalt("DC link has invalid type and/or control modes.");
      }
      if ((!strcmp(DCptr->Cont1,"ID") && !strcmp(DCptr->Cont2,"VD")) ||
          (!strcmp(DCptr->Cont1,"VD") && !strcmp(DCptr->Cont2,"ID"))){
        fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
        ErrorHalt("DC link has an invalid control mode: ID VD.");
      }
      if ((!strcmp(DCptr->Cont1,"PA") && !strcmp(DCptr->Cont2,"VD")) ||
          (!strcmp(DCptr->Cont1,"VD") && !strcmp(DCptr->Cont2,"PA"))){
        fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
        ErrorHalt("DC link has an invalid control mode: VD PA.");
      }
      if ((!strcmp(DCptr->Cont1,"ID") && !strcmp(DCptr->Cont2,"PA")) ||
          (!strcmp(DCptr->Cont1,"PA") && !strcmp(DCptr->Cont2,"ID"))){
        fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
        ErrorHalt("DC link has an invalid control mode: ID PA.");
      }
      Set1=fabs(GetValue(Line,37,6,2));
      Set2=fabs(GetValue(Line,43,6,2));
      if (Set1==0 || Set2==0){
        fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
        ErrorHalt("Zero DC control set point.");
      }
      if (!strcmp(DCptr->Cont1,"ID")) DCptr->Id=Set1;
      if (!strcmp(DCptr->Cont2,"ID")) DCptr->Id=Set2;
      if (!strcmp(DCptr->Cont1,"VD")) DCptr->Vd=Set1;
      if (!strcmp(DCptr->Cont2,"VD")) DCptr->Vd=Set2;
      if (!strcmp(DCptr->Cont1,"AL")) DCptr->Alfa=Set1;
      if (!strcmp(DCptr->Cont2,"AL")) DCptr->Alfa=Set2;
      if (!strcmp(DCptr->Cont1,"GA")) DCptr->Gamma=Set1;
      if (!strcmp(DCptr->Cont2,"GA")) DCptr->Gamma=Set2;
      if (!strcmp(DCptr->Cont1,"PA")) DCptr->P=Set1;
      if (!strcmp(DCptr->Cont2,"PA")) DCptr->P=Set2;
      if (!strcmp(DCptr->Cont1,"QA")) DCptr->Q=Set1;
      if (!strcmp(DCptr->Cont2,"QA")) DCptr->Q=Set2;
      if (!strcmp(DCptr->Cont1,"AT")) DCptr->Tap=Set1;
      if (!strcmp(DCptr->Cont2,"AT")) DCptr->Tap=Set2;
      DCptr->MVA=GetValue(Line,50,5,1);
    }

  /* --------------- DC element data -------------------------------- */
    else if (!strncmp(Line,"LD ",3)) {
      P=GetValue(Line,57,5,1);
      if (P!=0) flagEPRI=FALSE;

    /* ------------ EPRI Format Multiterminal Format ------------- */
      if (flagEPRI) {
        NdcEl++;
        GetStr(Line,4,8,8,Name);
        DCptr=(DCbusData *) DCbusInList(Name,Ndc);
        if (DCptr->N==0) { Ndc++; DCptr->N=Ndc;}
        GetStr(Line,12,8,8,Name);
        DCptrp=(DCbusData *) DCbusInList(Name,Ndc);
        if (DCptrp->N==0) { Ndc++; DCptrp->N=Ndc;}
        if (DCptr->N==DCptrp->N){
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("Both DC line buses are the same.");
        }
        R=GetValue(Line,29,6,2);
        L=GetValue(Line,35,6,2)/1000.;
        if (R<0.01) {
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("DC line is a short circuit or has negative R.");
        }
        if (ExistParameter('O') && L<=0) {
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("Program needs a DC L>0 to calculate ac/dc TEF.");
        }
        DCptr->Rd=DCptrp->Rd=R;
        DCptr->Ld=DCptrp->Ld=L;
        DCptr->To=DCptrp;
        DCptrp->To=DCptr;
        GetStr(Line,24,2,2,DCptr->Lzone);
        strcpy_s(DCptrp->Lzone,DCptr->Lzone);
        if (strcmp(DCptr->Zone,DCptr->Lzone)) DCptr->Meter=DCptrp->Meter=DCptr;
        else DCptr->Meter=DCptrp->Meter=DCptrp;
      }

    /* ---------------- WSC/BPA 2 Terminal Format ---------------- */
      else {
        NdcEl++;
        GetStr(Line,7,12,12,Name);
        DCptr=(DCbusData *) DCbusInList(Name,Ndc);
        if (DCptr->N==0) { Ndc++; DCptr->N=Ndc;}
        GetStr(Line,20,12,12,Name);
        DCptrp=(DCbusData *) DCbusInList(Name,Ndc);
        if (DCptrp->N==0) { Ndc++; DCptrp->N=Ndc;}
        if (DCptr->N==DCptrp->N){
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("Both DC line buses are the same.");
        }
        R=GetValue(Line,38,6,2);
        L=GetValue(Line,44,6,2)/1000.;
        if (R<0.01) {
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("DC line is a short circuit or has negative R.");
        }
        if (ExistParameter('O') && L<=0) {
          fprintf(stderr,"Input Line-> %d\n%s",LineNum,Line);
          ErrorHalt("Program needs a DC L>0 to calculate ac/dc TEF.");
        }
        DCptr->Rd=DCptrp->Rd=R;
        DCptrp->Ld+=L;
        DCptr->Ld=DCptrp->Ld;
        V=GetValue(Line,62,5,1);
        alpha=GetValue(Line,67,4,1);
        gamma=GetValue(Line,71,4,1);
        GetStr(Line,56,1,1,Mode);
        if (Mode[0]=='R' && P>0) {
          strcpy_s(DCptr->Type,"R");
          strcpy_s(DCptrp->Type,"I");
          strcpy_s(DCptr->Cont1,"PA");
          strcpy_s(DCptr->Cont2,"AL");
          strcpy_s(DCptrp->Cont1,"VD");
          strcpy_s(DCptrp->Cont2,"GA");
          DCptr->P=P;
          DCptr->Alfa=alpha;
          DCptrp->Vd=V-(P/V)*R;
          DCptrp->Gamma=gamma;
          DCptr->Meter=DCptrp->Meter=DCptr;
          DCptr->AlfaMax=90.;
          DCptrp->AlfaMax=180.;
        }
        else if (Mode[0]=='I' && P>0) {
          strcpy_s(DCptr->Type,"R");
          strcpy_s(DCptrp->Type,"I");
          strcpy_s(DCptr->Cont1,"VD");
          strcpy_s(DCptr->Cont2,"AL");
          strcpy_s(DCptrp->Cont1,"PA");
          strcpy_s(DCptrp->Cont2,"GA");
          DCptr->Vd=V;
          DCptr->Alfa=alpha;
          DCptrp->P=P;
          DCptrp->Gamma=gamma;
          DCptr->Meter=DCptrp->Meter=DCptrp;
          DCptr->AlfaMax=90.;
          DCptrp->AlfaMax=180.;
        }
        else if (Mode[0]=='R' && P<0) {
          strcpy_s(DCptr->Type,"I");
          strcpy_s(DCptrp->Type,"R");
          strcpy_s(DCptrp->Cont1,"PA");
          strcpy_s(DCptrp->Cont2,"AL");
          strcpy_s(DCptr->Cont1,"VD");
          strcpy_s(DCptr->Cont2,"GA");
          DCptrp->P=-P;
          DCptrp->Alfa=alpha;
          DCptr->Vd=V+(P/V)*R;
          DCptr->Gamma=gamma;
          DCptr->Meter=DCptrp->Meter=DCptr;
          DCptrp->AlfaMax=90.;
          DCptr->AlfaMax=180.;
        }
        else if (Mode[0]=='I' && P<0) {
          strcpy_s(DCptr->Type,"I");
          strcpy_s(DCptrp->Type,"R");
          strcpy_s(DCptrp->Cont1,"VD");
          strcpy_s(DCptrp->Cont2,"AL");
          strcpy_s(DCptr->Cont1,"PA");
          strcpy_s(DCptr->Cont2,"GA");
          DCptrp->Vd=V;
          DCptrp->Alfa=alpha;
          DCptr->P=-P;
          DCptr->Gamma=gamma;
          DCptr->Meter=DCptrp->Meter=DCptrp;
          DCptrp->AlfaMax=90.;
          DCptr->AlfaMax=180.;
        }
        DCptr->To=DCptrp;
        DCptrp->To=DCptr;
        strcpy_s(DCptr->Lzone,"  ");
        strcpy_s(DCptrp->Lzone,DCptr->Lzone);
      }
    }

}

