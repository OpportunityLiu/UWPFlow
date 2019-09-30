#include "readdata.h"

/* --------- Global Input File --------- */
extern FILE *InputDataFile;

/* --------------- SVCbusInList  ------------------- */
#ifdef ANSIPROTO
SVCbusData *SVCbusInList(char *BusName,INDEX N,ACbusData *ptrac,ACbusData *ptrac1)
#else
SVCbusData *SVCbusInList(BusName,N,ptrac,ptrac1)
char *BusName;
INDEX N;
ACbusData *ptrac,*ptrac1;
#endif
{
  SVCbusData *ptr,*ptrp,*ptrn;

  if(N==0) {
    ptrn=dataPtr->SVCbus;
    ptrp=NULL;
#ifdef WINDOWS
    dataPtr->SVCbus= new SVCbusData;
#else
    dataPtr->SVCbus=(SVCbusData *) malloc(sizeof(SVCbusData));
    if (dataPtr->SVCbus==NULL) {
      fclose(InputDataFile);
      ErrorHalt("Insufficient memory to allocate SVC bus data.");
      exit(ERROREXIT);
    }
#endif
    ptr=dataPtr->SVCbus;
    if(ptrn!=NULL)ptrn->Prev=ptr;
  } else {
    ptr=dataPtr->SVCbus;
    while (ptr!=NULL){
      if(!strcmp(ptr->Name,BusName)) return(ptr);
      ptrp=ptr;
      ptr=ptr->Next;
    }
  }
  strcpy(ptr->Name,BusName);
  ptr->N=0;
  strcpy(ptr->Type,"FS");
  strcpy(ptr->Cont,"AL");
  ptr->Vsvc=0;
  ptr->Xth_l=0;
  ptr->Xc=0;
  ptr->Xl=0;
  ptr->AlphaMin=0;
  ptr->AlphaMax=0;
  ptr->slope=0;
  ptr->SVC_base=0;
  ptr->Qsvc=0;
  ptr->Bv=0;
  ptr->alpha_svc=0;
  ptr->Vref=0;
  ptr->val=0;
  ptr->Vvar=0;
  ptr->From=ptrac;
  ptr->Ctrl=ptrac1;
  ptr->Prev=ptrp;
  ptr->Next=ptrn;
  return(ptr);
}


/* --------------- TCSCbusInList  ------------------- */
#ifdef ANSIPROTO
TCSCbusData *TCSCbusInList(char *BusName,INDEX N,ACbusData *ptrac,ACbusData *ptrac1)
#else
TCSCbusData *TCSCbusInList(BusName,N,ptrac,ptrac1)
char *BusName;
INDEX N;
ACbusData *ptrac,*ptrac1;
#endif
{
  TCSCbusData *ptr,*ptrp,*ptrn;

  if(N==0) {
    ptrn=dataPtr->TCSCbus;
    ptrp=NULL;
#ifdef WINDOWS
    dataPtr->TCSCbus=new TCSCbusData;
#else
    dataPtr->TCSCbus=(TCSCbusData *) malloc(sizeof(TCSCbusData));
    if (dataPtr->TCSCbus==NULL) {
      fclose(InputDataFile);
      ErrorHalt("Insufficient memory to allocate TCSC bus data.");
      exit(ERROREXIT);
    }
#endif
    ptr=dataPtr->TCSCbus;
    if(ptrn!=NULL)ptrn->Prev=ptr;
  } else {
    ptr=dataPtr->TCSCbus;
    while (ptr!=NULL){
      if(!strcmp(ptr->Name,BusName)) return(ptr);
      ptrp=ptr;
      ptr=ptr->Next;
    }
  }
  strcpy(ptr->Name,BusName);
  ptr->N=0;
  strcpy(ptr->Type,"FC");
  strcpy(ptr->Cont,"X");
  ptr->Xc=0;
  ptr->Xl=0;
  ptr->AlphaMin=0;
  ptr->AlphaMax=0;
  ptr->Control=0;
  ptr->Bset=0;
  ptr->TCSC_base=0;
  ptr->Ptcsc=0;
  ptr->Qtcsck=0;
  ptr->Qtcscm=0;
  ptr->Be=0;
  ptr->alpha_tcsc=0;
  ptr->Itcsc=0;
  ptr->delta_t=0;
  ptr->val=0;
  ptr->Max=0;
  ptr->From=ptrac;
  ptr->To=ptrac1;
  ptr->Prev=ptrp;
  ptr->Next=ptrn;
  return(ptr);
}


/* --------------- STATCOMbusInList  ------------------- */
#ifdef ANSIPROTO
STATCOMbusData *STATCOMbusInList(char *BusName,INDEX N,ACbusData *ptrac,ACbusData *ptrac1)
#else
STATCOMbusData *STATCOMbusInList(BusName,N,ptrac,ptrac1)
char *BusName;
INDEX N;
ACbusData *ptrac,*ptrac1;
#endif
{
  STATCOMbusData *ptr,*ptrp,*ptrn;

  if(N==0) {
    ptrn=dataPtr->STATCOMbus;
    ptrp=NULL;
#ifdef WINDOWS
    dataPtr->STATCOMbus= new STATCOMbusData;
#else
    dataPtr->STATCOMbus=(STATCOMbusData *) malloc(sizeof(STATCOMbusData));
    if (dataPtr->STATCOMbus==NULL) {
      fclose(InputDataFile);
      ErrorHalt("Insufficient memory to allocate STATCOM bus data.");
      exit(ERROREXIT);
    }
#endif
    ptr=dataPtr->STATCOMbus;
    if(ptrn!=NULL)ptrn->Prev=ptr;
  } else {
    ptr=dataPtr->STATCOMbus;
    while (ptr!=NULL){
      if(!strcmp(ptr->Name,BusName)) return(ptr);
      ptrp=ptr;
      ptr=ptr->Next;
    }
  }
  strcpy(ptr->Name,BusName);
  ptr->N=0;
  strcpy(ptr->Type,"FT");
  strcpy(ptr->Cont,"PW");
  strcpy(ptr->Cont1,"PW");
  ptr->I=0;
  ptr->theta=0;
  ptr->k=0;
  ptr->Vdc=0;
  ptr->alpha=0;
  ptr->R=0;
  ptr->G=0;
  ptr->B=0;
  ptr->Gc=0;
  ptr->Imin=0;
  ptr->Imax=0.;
  ptr->slope=0;
  ptr->P=0;
  ptr->Q=0;
  ptr->MVA=0;
  ptr->Vref=0;
  ptr->Contref=0;
  ptr->val=0;
  ptr->Vvar=0;
  ptr->From=ptrac;
  ptr->Ctrl=ptrac1;
  ptr->Prev=ptrp;
  ptr->Next=ptrn;
  return(ptr);
}

