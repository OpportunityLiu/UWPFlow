/*  Read input data.
    Transform the input data into the data structures needed by
    the power flow, i.e., AC, DC and FACTS bus data and interconnecting
    elements (lines, transformers with fixed taps, and PI equivalents).

    Main. */

#include "readdata.h"


/* ------- Global Variables ------ */
VALUETYPE AngSlack;
INDEX NdcEl,LineNum;
FILE *InputDataFile;

/* ------------------- ErrorDetect -------------------------- */
void ErrorDetect(void)
/* Detect inconsistencies in input data. */
{
  ACbusData *ACptr,*ACptrp,*ACptrs;
  AClist *ptrp,*ptrs;
  DClist *ptr;
  DCbusData *DCptr,*DCptrp;
  ElementData *Eptr,*Eptrp;
  ElementList *ELptr,*ELptrp;
  AreaData *Aptr,*Aptrp;
  SVCbusData *SVCptr;          /* FACTS */
  STATCOMbusData *STATCOMptr;  /* FACTS */
  VALUETYPE KVs;
  BOOLEAN flag,ieee;
  int i,j;

  if (Nac==0||(NacEl+NdcEl)==0) ErrorHalt("No AC buses and/or elements in input data.");
  if (Ndc==0) dataPtr->DCbus=NULL;
  if (NacEl==0) dataPtr->Element=NULL;
  if (Nsvc==0) dataPtr->SVCbus=NULL;         /* FACTS */
  if (Ntcsc==0) dataPtr->TCSCbus=NULL;       /* FACTS */
  if (Nstatcom==0) dataPtr->STATCOMbus=NULL; /* FACTS */
  if (Narea<2) {dataPtr->Area=NULL; Narea=0;}
  if (Nslack==0) ErrorHalt("No angle reference Bus for AC system.");

 /* ------------------- AC buses ------------------------ */
  ACptr=dataPtr->ACbus;
  ieee=ExistParameter('I');
  while(ACptr!=NULL){
    if (ACptr->V==0){
      fprintf(stderr,"ERROR: AC/DC bus %d %s has not been defined.\n",ACptr->Num,ACptr->Name);
      fprintf(stderr,"       Check the bus input cards.\n");
      InputError=TRUE;
    }
                           /* FACTS */
    if (ACptr->Elem==NULL && ACptr->DC==NULL && ACptr->TCSC==NULL){
      fprintf(stderr,"ERROR: AC bus %d %s is isolated.\n",ACptr->Num,ACptr->Name);
      fprintf(stderr,"       Check the AC/DC/FACTS input cards.\n");
      InputError=TRUE;
    }
                       /* END OF FACTS */
    if (strpbrk(ACptr->Type,"G")){
#ifdef WINDOWS
      ptrp= new AClist;
#else
      ptrp= (AClist *) malloc(sizeof(AClist));
      if (ptrp==NULL) {ErrorHalt("Insufficient memory to allocate list of AC controlled buses."); exit(ERROREXIT);}
#endif
      ACptr->ContBus=ptrp;
      ACptr->ContBus->AC=ACptr->Cont;
      ACptr->ContBus->Next=ACptr->ContBus->Prev=NULL;
      if (ACptr->Cont!=NULL) {
        if (!strpbrk(ACptr->Cont->Type,"C")){
         fprintf(stderr,"ERROR: The voltage controlled bus %d %s of PV bus\n",ACptr->Cont->Num,ACptr->Cont->Name);
         fprintf(stderr,"       %d %s, is not a PQ bus.  Check the AC bus input data.\n",ACptr->Num,ACptr->Name);
         InputError=TRUE;
        } else {
#ifdef WINDOWS
         ptrp= new AClist;
#else
         ptrp= (AClist *) malloc(sizeof(AClist));
         if (ptrp==NULL) {ErrorHalt("Insufficient memory to allocate list of AC controlling buses."); exit(ERROREXIT);}
#endif
         ptrp->AC=ACptr;
         ptrp->Prev=NULL;
         if(ACptr->Cont->ContBus==NULL) {
              ACptr->Cont->ContBus=ptrp;
              ACptr->Cont->ContBus->Next=ACptr->Cont->ContBus->Prev=NULL;
         } else {
              ptrs=ACptr->Cont->ContBus;
              ACptr->Cont->ContBus=ptrp;
              ptrs->Prev=ptrp;
              ptrp->Next=ptrs;
         }
        }
        if (!QRcont) { ACptr->Cont->Cont=ACptr->Cont; ACptr->Cont=NULL; }
        else {
          if (flag2Vcontrol) {
            if (ACptr->Cont->Kbg1>0) ACptr->Kbg1=ACptr->Qmax/ACptr->Cont->Kbg1;
            else                     ACptr->Kbg1=1;
            if (ACptr->Cont->Kbg2<0) ACptr->Kbg2=ACptr->Qmin/ACptr->Cont->Kbg2;
            else                     ACptr->Kbg2=1;
            ACptr->Kbg=ACptr->Kbg1;
          }
		  else {
			ACptr->Kbg1=ACptr->Kbg;
			ACptr->Kbg2=ACptr->Kbg;
		  }
          ACptr->Cont->Qr=ACptr->Cont->Qr+ACptr->Qg;
        }
      } else {
         fprintf(stderr,"ERROR: The remote controlled bus of PV bus %d %s\n",ACptr->Num,ACptr->Name);
         fprintf(stderr,"       has not been defined.  Check the AC bus input data.\n");
         InputError=TRUE;
      }
    }
    if (strpbrk(ACptr->Type,"C") && ACptr->Kbg<1){
      fprintf(stderr,"ERROR: The voltage controlled bus %d %s does not have\n",ACptr->Num,ACptr->Name);
      fprintf(stderr,"       any generator controlling the voltage. Check the AC bus input data.\n");
      InputError=TRUE;
    }
    if (Rcont && strpbrk(ACptr->Type,"T")) ACptr->Cont=NULL;
    if (Narea>1){
      if(ACptr->Area==NULL) for (Aptr=dataPtr->Area; Aptr!=NULL; Aptr=Aptr->Next){
         for (i=1;i<=10;i++){
            if(!strcmp(ACptr->Zone,Aptr->Zone[i])) {
              ACptr->Area=Aptr;
              if (strpbrk(ACptr->Type,"S")) Aptr->i++;
              if (Aptr->i>1) {
                fprintf(stderr,"ERROR: Area %d %s has 2 slack buses.\n",Aptr->N,Aptr->Name);
                fprintf(stderr,"       Check AC area and bus input data.\n");
                InputError=TRUE;
              }
              break;
            }
         }
         if (ACptr->Area!=NULL) break;
      }
      if (ACptr->Area==NULL) {
         fprintf(stderr,"ERROR: AC/DC bus %d %s is not in any area.\n",ACptr->Num,ACptr->Name);
         fprintf(stderr,"       Check area and bus input data.\n");
         InputError=TRUE;
      }
      else {
        if (!strcmp(ACptr->Area->Name,"")) {
           fprintf(stderr,"ERROR: Area %d has not been defined.\n",ACptr->Area->N);
           fprintf(stderr,"       Check area and bus input data.\n");
           InputError=TRUE;
        }
        if(ACptr==ACptr->Area->Slack && !strpbrk(ACptr->Type,"S")) {
          strcat(ACptr->Type,"A");
          ACptr->Area->i++;
          if (ACptr->Area->i>1) {
            fprintf(stderr,"ERROR: Area %d %s has 2 slack buses.\n",ACptr->Area->N,ACptr->Area->Name);
            fprintf(stderr,"       Check AC area and bus input data.\n");
            InputError=TRUE;
          }
        }
        Aptr=ACptr->Area;
        ptrp=Aptr->AC;
#ifdef WINDOWS
        Aptr->AC= new AClist;
#else
        Aptr->AC=(AClist *) malloc(sizeof(AClist));
        if (Aptr->AC==NULL) {ErrorHalt("Insufficient memory to allocate area data."); exit(ERROREXIT);}
#endif
        Aptr->AC->AC=ACptr;
        Aptr->AC->Next=ptrp;
        Aptr->AC->Prev=NULL;
        if (ptrp!=NULL) ptrp->Prev=Aptr->AC;
        if (ExistParameter('6') && Aptr->BSptr==NULL) Aptr->Slack=Aptr->BSptr=ACptr;
      }
    } else ACptr->Area=NULL;
    if(ACptr->Ang>=1000) ACptr->Ang=AngSlack;
    if(ACptr->DPg==0 && strpbrk(ACptr->Type,"AS")) {
       if (ExistParameter('6')) {if (strpbrk(ACptr->Type,"S")) ACptr->DPg=1;}
       else ACptr->DPg=1;
    }
    if (ACptr->PgMax<=0) ACptr->PgMax=99999999.;
    if (ACptr->Smax<=0) ACptr->Smax=99999999.;
    ACptr->DPG=ACptr->DPg;
    if (strpbrk(ACptr->Type,"X")) {
      if (ACptr->Cont==NULL) {
        fprintf(stderr,"ERROR: The reactance controlled bus %d %s does not have\n",ACptr->Num,ACptr->Name);
        fprintf(stderr,"       a controlling bus. Check the AC input data.\n");
        InputError=TRUE;
      }
      else {
        ACptrp=ACptr->Cont;
        ACptr->Vmax=ACptrp->Vmax;
        ACptr->Vmin=ACptrp->Vmin;
        if (ACptr->Vmin>=ACptr->Vmax) {
          fprintf(stderr,"ERROR: The reactance controlling bus %d %s has inconsistent\n",ACptrp->Num,ACptrp->Name);
          fprintf(stderr,"       voltage limits. Check the AC input data.\n");
          InputError=TRUE;
        }
        else if (ACptrp->steps==0) {
          fprintf(stderr,"ERROR: The reactance controlling bus %d %s has zero MVAr steps.\n",ACptrp->Num,ACptrp->Name);
          fprintf(stderr,"       Check the AC input data.\n");
          InputError=TRUE;
        }
      }
    }

 /* -------------------------- AC elements --------------------------- */
    i=0;
    for (ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next){
      Eptr=ELptr->Eptr;
      if((!strcmp(Eptr->Type,"R") && !strpbrk(ACptr->Type,"T")) ||
         (!strcmp(Eptr->Type,"RV") && !strpbrk(ACptr->Type,"R"))){
         fprintf(stderr,"ERROR: LTC volt. controlled bus %d %s is not PQ.\n",ACptr->Num,ACptr->Name);
         fprintf(stderr,"       Check the AC/DC bus input cards.\n");
         InputError=TRUE;
      }
      if(!strcmp(Eptr->Type,"R") || !strcmp(Eptr->Type,"RV")) i++;
    }
    if((ACptr->Reg==NULL || (ACptr->Reg!=NULL && i==0)) && strpbrk(ACptr->Type,"T")) {
      fprintf(stderr,"ERROR: LTC volt. controlled bus %d %s does not have LTCs.\n",ACptr->Num,ACptr->Name);
      fprintf(stderr,"       Check the AC/DC bus and element input cards.\n");
      InputError=TRUE;
    }
    if (strpbrk(ACptr->Type,"S,A")) {
      ptrp=dataPtr->KGbus;
#ifdef WINDOWS
      ptrs= new AClist;
#else
      ptrs=(AClist *) malloc(sizeof(AClist));
      if (ptrs==NULL) {ErrorHalt("Insufficient memory to allocate area data."); exit(ERROREXIT);}
#endif
      ptrs->AC=ACptr;
      ptrs->Next=ptrp;
      ptrs->Prev=NULL;
      dataPtr->KGbus=ptrs;
      if (ptrp!=NULL) ptrp->Prev=ptrs;
    }
    ACptr->PG=ACptr->Pg;
    ACptr->PL=ACptr->Pl;
    ACptr->QL=ACptr->Ql;
    ACptrp=ACptr->Next;
    if(ieee) {
      ACptr->Next=ACptr->Prev;
      ACptr->Prev=ACptrp;
      dataPtr->ACbus=ACptr;
    }
    ACptr=ACptrp;
  }
  Eptr=dataPtr->Element;
  while(Eptr!=NULL){
    if(strpbrk(Eptr->Type,"R") && !strcmp(Eptr->Zone,"")) {
      fprintf(stderr,"ERROR: Reg. transf. from %d %s to %d %s\n",
               Eptr->From->Num,Eptr->From->Name,Eptr->To->Num,Eptr->To->Name);
      fprintf(stderr,"       has not been completely defined.  Check T cards on WSCC format.\n");
      InputError=TRUE;
    }
    if (!strcmp(Eptr->Type,"R") || !strcmp(Eptr->Type,"RV")) strcpy(Eptr->Ctype,"V");
    else if (strpbrk(Eptr->Type,"PM")) strcpy(Eptr->Ctype,"P");
    else if (strpbrk(Eptr->Type,"QN")) strcpy(Eptr->Ctype,"Q");
    if (strcmp(Eptr->Owner,"")){
      if (!strcmp(Eptr->Owner,Eptr->From->Owner)) strcpy(Eptr->Zone,Eptr->From->Zone);
      else strcpy(Eptr->Zone,Eptr->To->Zone);
    }
    if (Narea>1) {
      Aptr=Eptr->From->Area;
      Aptrp=Eptr->To->Area;
      if(Aptr!=Aptrp  && Aptr!=NULL && Aptrp!=NULL){
         Aptr->Elem=(ElementList *) AddElemToList(Aptr->Elem,Eptr);
         Aptrp->Elem=(ElementList *) AddElemToList(Aptrp->Elem,Eptr);
      }
      if (Eptr->Meter==NULL) {
         if (strcmp(Eptr->Owner,"")){
            if (strcmp(Eptr->Owner,Eptr->From->Owner)) Eptr->Meter=Eptr->From;
            else Eptr->Meter=Eptr->To;
         }
         else if (Eptr->Area==Aptr) Eptr->Meter=Eptr->To;
         else Eptr->Meter=Eptr->From;
      }
    }
    Eptrp=Eptr->Next;
    Eptr->Next=Eptr->Prev;
    Eptr->Prev=Eptrp;
    dataPtr->Element=Eptr;
    Eptr=Eptrp;
  }

 /* ------------------- DC buses ------------------------ */
  for(DCptr=dataPtr->DCbus;DCptr!=NULL;DCptr=DCptr->Next){
    ACptr=DCptr->AC;
    if (ACptr==NULL || !strcmp(DCptr->Type,"")){
      fprintf(stderr,"ERROR: DC bus %8s has not been fully defined in the input data.\n",DCptr->Name);
      fprintf(stderr,"       Check the BD, LD, and/or BZ cards.\n");
      InputError=TRUE;
    }

    /* --------- Fix data read in WSCC/BPA Format --------- */
    if (DCptr->Xc==0) {
      for (ELptrp=ELptr=ACptr->Elem;ELptr!=NULL;ELptrp=ELptr,ELptr=ELptr->Next) {
        Eptr=ELptr->Eptr;
        if (!strcmp(Eptr->From->Name,DCptr->Name)) { ACptrp=Eptr->From; break;}
        if (!strcmp(Eptr->To->Name,DCptr->Name))  { ACptrp=Eptr->To; break;}
      }
      DCptr->Name[8]='\0';
      if (ELptr==NULL) {
        fprintf(stderr,"ERROR: The tranformer for DC bus %8s has not been defined in the\n",DCptr->Name);
        fprintf(stderr,"       input data. Check the T data cards.\n");
        InputError=TRUE;
      } else {
        DCptr->Xc=-Eptr->B/(Eptr->G*Eptr->G+Eptr->B*Eptr->B);
        DCptr->Xc=(DCptr->Xc*ACptr->KV*ACptr->KV/Sn)*DCptr->Nbr;
        DCptr->Tap=Eptr->Tap;
        if (strpbrk(Eptr->Type,"P,Q,M,N")) {
           fprintf(stderr,"ERROR: DC bus %8s regulating transfomer may only control voltage.\n",DCptr->Name);
           fprintf(stderr,"       Check the related R transfomer input cards.\n");
           InputError=TRUE;
        }
        if (!strcmp(Eptr->Type,"R")) {
          if (Eptr->Cont!=ACptrp) {
            fprintf(stderr,"ERROR: DC bus %8s regulating transfomer must control voltage\n",DCptr->Name);
            fprintf(stderr,"       of DC bus.  Check the related R transfomer input cards.\n");
            InputError=TRUE;
          }
          NregV--;
          DCptr->TapMax=Eptr->Tmax;
          DCptr->TapMin=Eptr->Tmin;
        } else {
          DCptr->TapMax=1.1;
          DCptr->TapMin=0.9;
        }
      }
      /* Remove extra element from AC data base */
      if (ELptrp==ELptr) {
        ACptr->Elem=ELptr->Next;
        if (ACptr->Elem==NULL) {
           fprintf(stderr,"ERROR: AC bus %d %s is isolated from the rest of the system.\n",ACptr->Num,ACptr->Name);
           fprintf(stderr,"       Check the related AC element input cards.\n");
           InputError=TRUE;
        }
      }
      else ELptrp->Next=ELptr->Next;
#ifndef WINDOWS
      free(ELptr);
#else
      delete ELptr;
#endif
      Aptr=ACptr->Area;
      if (Aptr!=NULL) {
        for (ELptrp=ELptr=Aptr->Elem;ELptr!=NULL;ELptrp=ELptr,ELptr=ELptr->Next) if (Eptr==ELptr->Eptr) break;
        if (ELptr!=NULL) {
          if (ELptrp==ELptr) Aptr->Elem=ELptr->Next;
          else ELptrp->Next=ELptr->Next;
#ifndef WINDOWS
          free(ELptr);
#else
          delete ELptr;
#endif
        }
      }
      NacEl--;
      if (Eptr->Prev==NULL) {
        Eptrp=Eptr->Next;
        dataPtr->Element=Eptrp;
        Eptrp->Prev=NULL;
      } else {
        Eptrp=Eptr->Prev;
        Eptrp->Next=Eptr->Next;
        if (Eptr->Next!=NULL) Eptr->Next->Prev=Eptrp;
      }
#ifndef WINDOWS
      free(Eptr);
#else
      delete Eptr;
#endif
      /* Remove DC bus from AC data base */
      Aptr=ACptr->Area;
      if (Aptr!=NULL) {
        for (ptrp=Aptr->AC;ptrp!=NULL;ptrp=ptrp->Next) if (ACptrp==ptrp->AC) break;
        if (ptrp!=NULL) {
          if (ptrp->Prev==NULL) {
            ptrs=ptrp->Next;
            Aptr->AC=ptrs;
            ptrs->Prev=NULL;
          } else {
            ptrs=ptrp->Prev;
            ptrs->Next=ptrp->Next;
            if (ptrp->Next!=NULL) ptrp->Next->Prev=ptrs;
          }
#ifndef WINDOWS
          free(ptrp);
#else
          delete ptrp;
#endif
        }
      }
      Nac--;
      for(ACptrs=ACptrp->Next;ACptrs!=NULL;ACptrs->N--,ACptrs=ACptrs->Next);
      if (ACptrp->Prev==NULL) {
        ACptrs=ACptrp->Next;
        dataPtr->ACbus=ACptrs;
        ACptrs->Prev=NULL;
      } else {
        ACptrs=ACptrp->Prev;
        ACptrs->Next=ACptrp->Next;
        if (ACptrp->Next!=NULL) ACptrp->Next->Prev=ACptrs;
      }
#ifndef WINDOWS
      free(ACptrp);
#else
      delete ACptrp;
#endif
    }


    if (DCptr->MVA>0) {
      KVs=ACptr->KV*DCptr->Ntrf/DCptr->Nbr;
      DCptr->Xc=DCptr->Xc/DCptr->Nbr;
      DCptr->Xc=(DCptr->Xc*KVs*KVs/DCptr->MVA)*DCptr->Nbr;
    }
    if (DCptr->Xc<=0.1){
      fprintf(stderr,"ERROR: DC bus %8s has a negative or zero commutation reactance ->\n",DCptr->Name);
      fprintf(stderr,"       Xc=%6.3lf (Ohms)\n",DCptr->Xc);
      fprintf(stderr,"       Check the BD and/or BZ cards.\n");
      InputError=TRUE;
    }
    if (DCptr->To==NULL){
      fprintf(stderr,"ERROR: DC bus %8s is isolated.\n",DCptr->Name);
      fprintf(stderr,"       Check the DC line input cards.\n");
      InputError=TRUE;
    }
    DCptr->Area=ACptr->Area;
    if(!strcmp(DCptr->Cont1,"AL")||!strcmp(DCptr->Cont2,"AL")) {
      if (DCptr->Alfa>DCptr->AlfaMax || DCptr->Alfa<DCptr->AlfaMin) {
         fprintf(stderr,"ERROR: DC bus %8s has ALPHA outside its limits.\n",DCptr->Name);
         fprintf(stderr,"       Check the BD (limits) and/or BZ card.\n");
         InputError=TRUE;
      }
      DCptr->Gamma=180-DCptr->Alfa;
    }
    else if(!strcmp(DCptr->Cont1,"GA")||!strcmp(DCptr->Cont2,"GA")) {
      if (DCptr->Gamma<DCptr->GammaMin) {
         fprintf(stderr,"ERROR: DC bus %8s has GAMMA outside its limits.\n",DCptr->Name);
         fprintf(stderr,"       Check the BD (limits) and/or BZ card.\n");
         InputError=TRUE;
      }
      DCptr->Alfa=180-DCptr->Gamma;
    }
    else if(!strcmp(DCptr->Type,"R")) {
      DCptr->Alfa=DCptr->AlfaMin;
      DCptr->Gamma=180-DCptr->Alfa;
    }
    else {
      DCptr->Gamma=DCptr->GammaMin;
      DCptr->Alfa=180-DCptr->Gamma;
    }
    if (DCptr->Alfa>DCptr->AlfaMax || DCptr->Alfa<DCptr->AlfaMin || DCptr->Gamma<DCptr->GammaMin)
         fprintf(stderr,"***Warning: DC bus %8s could have wrong ALPHA or GAMMA limits.\n",DCptr->Name);
    if (DCptr->Tap<=0) DCptr->Tap=1;
    if (DCptr->Tap<DCptr->TapMin) DCptr->Tap=DCptr->TapMin;
    if (DCptr->Tap>DCptr->TapMax) DCptr->Tap=DCptr->TapMax;
  }

 /* ------------------- DC elements ------------------------ */
  if (Ndc!=(2*NdcEl)) {
    fprintf(stderr,"ERROR: There are inconsistencies between the DC bus and DC line input data.\n");
    fprintf(stderr,"       Check DC input data and remember that the program just allows for \n");
    fprintf(stderr,"       two-terminal HVDC links.\n");
    InputError=TRUE;
  }
  for(DCptr=dataPtr->DCbus;DCptr!=NULL;DCptr=DCptr->Next){
    DCptrp=DCptr->To;
    if (DCptr->N!=0 && DCptrp->N!=0){
      if ((strcmp(DCptr->Type,"R") || strcmp(DCptrp->Type,"I")) &&
         (strcmp(DCptr->Type,"I") || strcmp(DCptrp->Type,"R"))){
         fprintf(stderr,"ERROR: Both converters for the DC link between %8s and %8s\n",
                  DCptr->Name,DCptrp->Name);
         fprintf(stderr,"       are either rectifiers or inverters.\n");
         InputError=TRUE;
      }
      if ((!strcmp(DCptr->Cont1,"ID") || !strcmp(DCptr->Cont2,"ID")) &&
           (!strcmp(DCptrp->Cont1,"ID") || !strcmp(DCptrp->Cont2,"ID"))){
         fprintf(stderr,"ERROR: Both converters for the DC link between %8s and %8s\n",
                  DCptr->Name,DCptrp->Name);
         fprintf(stderr,"       are controlling the current.\n");
         InputError=TRUE;
      }
      if (DCptr->Area==DCptrp->Area||!strcmp(DCptr->Zone,DCptrp->Zone))   DCptr->Meter=DCptrp->Meter=NULL;
      if (DCptr->Area!=DCptrp->Area) for(i=1;i<=2;i++) {
         if (i==1) Aptr=DCptr->Area;
         else Aptr=DCptrp->Area;
         if (Aptr!=NULL) {
           ptr=Aptr->DC;
#ifdef WINDOWS
           Aptr->DC= new DClist;
#else
           Aptr->DC=(DClist *) malloc(sizeof(DClist));
           if (Aptr->DC==NULL) {ErrorHalt("Insufficient memory to allocate area data."); exit(ERROREXIT);}
#endif
           if (i==1) Aptr->DC->DC=DCptr;
           else Aptr->DC->DC=DCptrp;
           Aptr->DC->Next=ptr;
         }
      }
      if (DCptr->Id>0) DCptrp->Id=DCptr->Id;
      else if (DCptrp->Id>0) DCptr->Id=DCptrp->Id;
      DCptr->N=0;
      DCptrp->N=0;
    }
  }

 /* -------------------------- Areas --------------------------- */
  for (Aptr=dataPtr->Area; Aptr!=NULL; Aptr=Aptr->Next) {
    flag=TRUE;
    i=0; j=0;
    for(ELptr=Aptr->Elem;ELptr!=NULL;ELptr=ELptr->Next){
      Eptr=ELptr->Eptr; flag=FALSE; i++;
      if(!strcmp(Eptr->Type,"RP")) j++;
      if (Eptr->From->Area!=Aptr) Aptrp=Eptr->From->Area;
      else Aptrp=Eptr->To->Area;
      if (!Acont && strpbrk(Aptr->Slack->Type,"S") && !strpbrk(Aptrp->Slack->Type,"S")) ExpandSlack(Aptr->Slack,Aptrp);
    }
    if (!flag && i==j) {
      fprintf(stderr,"ERROR: All tie lines for area %d %s are P reg. transf.\n",Aptr->N,Aptr->Name);
      fprintf(stderr,"       Change at least one reg. transf. to a standard one.\n");
      InputError=TRUE;
    }
    if (flag && !strpbrk(Aptr->Slack->Type,"S")) {
      strcat(Aptr->Slack->Type,"S");
      Nslack++;
    }
    ptrp=Aptr->AC;
    while(ptrp!=NULL){
      ptrs=ptrp->Next;
      ptrp->Next=ptrp->Prev;
      ptrp->Prev=ptrs;
      if (ptrs==NULL) Aptr->AC=ptrp;
      ptrp=ptrs;
    }
  }

 /* -------------------------- FACTS --------------------------- */
  for(SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
    ACptr=SVCptr->Ctrl;
    if (ACptr->Cont==NULL) {
      fprintf(stderr,"ERROR: The SVC controlled bus %d %s is already controlled.\n",ACptr->N,ACptr->Name);
      fprintf(stderr,"       Check the AC bus cards.\n");
      InputError=TRUE;
    }
    ACptr=SVCptr->From;
    if (ACptr->Cont==NULL) {
      fprintf(stderr,"ERROR: The SVC bus %d %s is a voltage controlled bus.\n",ACptr->N,ACptr->Name);
      fprintf(stderr,"       Check the AC bus cards.\n");
      InputError=TRUE;
    }
    for (ELptr=ACptr->Elem,i=0;ELptr!=NULL;ELptr=ELptr->Next,i++);
    if (i>1) {
      fprintf(stderr,"ERROR: The SVC bus %d %s has more than one AC element connected to it.\n",ACptr->N,ACptr->Name);
      fprintf(stderr,"       Check the AC element cards.\n");
      InputError=TRUE;
    }
  }
  for(STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
    ACptr=STATCOMptr->Ctrl;
    if (ACptr->Cont==NULL) {
      fprintf(stderr,"ERROR: The STATCOM controlled bus %d %s is already controlled.\n",ACptr->N,ACptr->Name);
      fprintf(stderr,"       Check the AC bus cards.\n");
      InputError=TRUE;
    }
    ACptr=STATCOMptr->From;
    if (ACptr->Cont==NULL) {
      fprintf(stderr,"ERROR: The STATCOM bus %d %s is a voltage controlled bus.\n",ACptr->N,ACptr->Name);
      fprintf(stderr,"       Check the AC bus cards.\n");
      InputError=TRUE;
    }
  }
}

/* --------- ExpandSlack --------- */
void ExpandSlack(ACbusData *BSptr,AreaData *Aptr)
{
  AreaData *Aptrp;
  ElementList *ELptr;
  ElementData *Eptr;

  Aptr->Slack=BSptr;
  for(ELptr=Aptr->Elem;ELptr!=NULL;ELptr=ELptr->Next) {
    Eptr=ELptr->Eptr;
    if (Eptr->From->Area!=Aptr) Aptrp=Eptr->From->Area;
    else Aptrp=Eptr->To->Area;
    if (!strpbrk(Aptrp->Slack->Type,"S")) ExpandSlack(Aptr->Slack,Aptrp);
  }
}

/* --------- WriteSummary --------- */
void WriteSummary(void)
{
  ACbusData *ACptr;
  int i;

  fprintf(stderr,"Summary of input data for case:\n");
  i=0;
  while(i<=2 && dataPtr->Title[0][0]!='\0'){
    fprintf(stderr,"%s",dataPtr->Title[i]);
    i++;
  }
  fprintf(stderr,"            AC buses -> %d\n",Nac);
  fprintf(stderr,"            PV buses -> %d\n",Nvolt);
  fprintf(stderr,"            X buses  -> %d\n",NXvolt);
  fprintf(stderr,"            Z buses  -> %d\n",NZvolt);
  fprintf(stderr,"            AC elem. -> %d\n",NacEl);
  fprintf(stderr,"         V Reg. Trf. -> %d\n",NregV);
  fprintf(stderr,"        PQ Reg. Trf. -> %d\n",NregPQ);
  fprintf(stderr,"            DC buses -> %d\n",Ndc);
  fprintf(stderr,"            DC lines -> %d\n",NdcEl);
  fprintf(stderr,"                SVCs -> %d\n",Nsvc);     /* FACTS */
  fprintf(stderr,"               TCSCs -> %d\n",Ntcsc);    /* FACTS */
  fprintf(stderr,"            STATCOMs -> %d\n",Nstatcom); /* FACTS */
  fprintf(stderr,"           No. Areas -> %d\n",Narea);
  fprintf(stderr,"   Reference Bus(es) -> ");
  i=0;
  for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next)
    if(strpbrk(ACptr->Type,"S")){
      if (i>0)   fprintf(stderr,"                        ");
      fprintf(stderr,"%d %s (Angle=%6.2lf deg.)\n",ACptr->Num,ACptr->Name,ACptr->Ang*180/PI);
      i++;
    }
  if (i==0) fprintf(stderr,"\n");
  return;
}



/* ------------------------- ReadData ------------------------------- */
void ReadData(char *Name)
 /* Main routine. */
{
	InputDataFile=(FILE *) OpenInput(Name);  
#ifdef WINDOWS
	dataPtr= new Data;
#else
  dataPtr=(Data *) malloc(sizeof(Data));
  if (dataPtr==NULL) {
    fclose(InputDataFile);
    ErrorHalt("Insufficient memory to read input data.");
    exit(ERROREXIT);
  }
#endif
  dataPtr->Title[0][0]='\0';
  dataPtr->ACbus=NULL;
  dataPtr->DCbus=NULL;
  dataPtr->Element=NULL;
  dataPtr->Area=NULL;
  dataPtr->KGbus=NULL;
  dataPtr->SVCbus=NULL;     /* FACTS */
  dataPtr->TCSCbus=NULL;    /* FACTS */
  dataPtr->STATCOMbus=NULL; /* FACTS */
  Nac=0; Ndc=0; NacEl=0; NdcEl=0;
  Nsvc=Nstatcom=0; Ntcsc=NtcscVar=0;  /* FACTS */
  LineNum=0;  Nvolt=0; Nslack=0;
  Narea=0; NregPQ=NregV=0;
  NZvolt=NXvolt=0;
  InputError=FALSE;
  flag2Vcontrol=ExistParameter('#');
  if(ExistParameter('I')) ReadIEEE();
  else if(ExistParameter('6')) ReadITALY();
  else ReadWSCC();
  ErrorDetect();
  WriteSummary();
  if (InputError==TRUE) {
      fprintf(stderr,"*** The data has errors! Please review the input file. ***\n");
      exit(ERROREXIT);
  }
  else fprintf(stderr,"*** The data has been read successfully ***\n");
#ifdef WINDOWS
  delete[] Name;
#else
  free(Name);
#endif

}
