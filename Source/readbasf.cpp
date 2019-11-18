#include "readdata.h"

/* --------- Global Input File --------- */
extern FILE *InputDataFile;

/* --------------- SVCbusInList  ------------------- */
SVCbusData *SVCbusInList(char *BusName, INDEX N, ACbusData *ptrac,
                         ACbusData *ptrac1) {
  SVCbusData *ptr, *ptrp, *ptrn;

  if (N == 0) {
    ptrn = dataPtr->SVCbus;
    ptrp = nullptr;
    dataPtr->SVCbus = new SVCbusData;
    ptr = dataPtr->SVCbus;
    if (ptrn != nullptr)
      ptrn->Prev = ptr;
  } else {
    ptr = dataPtr->SVCbus;
    while (ptr != nullptr) {
      if (!strcmp(ptr->Name, BusName))
        return (ptr);
      ptrp = ptr;
      ptr = ptr->Next;
    }
  }
  strcpy(ptr->Name, BusName);
  ptr->N = 0;
  strcpy(ptr->Type, "FS");
  strcpy(ptr->Cont, "AL");
  ptr->Vsvc = 0;
  ptr->Xth_l = 0;
  ptr->Xc = 0;
  ptr->Xl = 0;
  ptr->AlphaMin = 0;
  ptr->AlphaMax = 0;
  ptr->slope = 0;
  ptr->SVC_base = 0;
  ptr->Qsvc = 0;
  ptr->Bv = 0;
  ptr->alpha_svc = 0;
  ptr->Vref = 0;
  ptr->val = 0;
  ptr->Vvar = 0;
  ptr->From = ptrac;
  ptr->Ctrl = ptrac1;
  ptr->Prev = ptrp;
  ptr->Next = ptrn;
  return (ptr);
}

/* --------------- TCSCbusInList  ------------------- */
TCSCbusData *TCSCbusInList(char *BusName, INDEX N, ACbusData *ptrac,
                           ACbusData *ptrac1) {
  TCSCbusData *ptr, *ptrp, *ptrn;

  if (N == 0) {
    ptrn = dataPtr->TCSCbus;
    ptrp = nullptr;
    dataPtr->TCSCbus = new TCSCbusData;
    ptr = dataPtr->TCSCbus;
    if (ptrn != nullptr)
      ptrn->Prev = ptr;
  } else {
    ptr = dataPtr->TCSCbus;
    while (ptr != nullptr) {
      if (!strcmp(ptr->Name, BusName))
        return (ptr);
      ptrp = ptr;
      ptr = ptr->Next;
    }
  }
  strcpy(ptr->Name, BusName);
  ptr->N = 0;
  strcpy(ptr->Type, "FC");
  strcpy(ptr->Cont, "X");
  ptr->Xc = 0;
  ptr->Xl = 0;
  ptr->AlphaMin = 0;
  ptr->AlphaMax = 0;
  ptr->Control = 0;
  ptr->Bset = 0;
  ptr->TCSC_base = 0;
  ptr->Ptcsc = 0;
  ptr->Qtcsck = 0;
  ptr->Qtcscm = 0;
  ptr->Be = 0;
  ptr->alpha_tcsc = 0;
  ptr->Itcsc = 0;
  ptr->delta_t = 0;
  ptr->val = 0;
  ptr->Max = 0;
  ptr->From = ptrac;
  ptr->To = ptrac1;
  ptr->Prev = ptrp;
  ptr->Next = ptrn;
  return (ptr);
}

/* --------------- STATCOMbusInList  ------------------- */
STATCOMbusData *STATCOMbusInList(char *BusName, INDEX N, ACbusData *ptrac,
                                 ACbusData *ptrac1) {
  STATCOMbusData *ptr, *ptrp, *ptrn;

  if (N == 0) {
    ptrn = dataPtr->STATCOMbus;
    ptrp = nullptr;
    dataPtr->STATCOMbus = new STATCOMbusData;
    ptr = dataPtr->STATCOMbus;
    if (ptrn != nullptr)
      ptrn->Prev = ptr;
  } else {
    ptr = dataPtr->STATCOMbus;
    while (ptr != nullptr) {
      if (!strcmp(ptr->Name, BusName))
        return (ptr);
      ptrp = ptr;
      ptr = ptr->Next;
    }
  }
  strcpy(ptr->Name, BusName);
  ptr->N = 0;
  strcpy(ptr->Type, "FT");
  strcpy(ptr->Cont, "PW");
  strcpy(ptr->Cont1, "PW");
  ptr->I = 0;
  ptr->theta = 0;
  ptr->k = 0;
  ptr->Vdc = 0;
  ptr->alpha = 0;
  ptr->R = 0;
  ptr->G = 0;
  ptr->B = 0;
  ptr->Gc = 0;
  ptr->Imin = 0;
  ptr->Imax = 0.;
  ptr->slope = 0;
  ptr->P = 0;
  ptr->Q = 0;
  ptr->MVA = 0;
  ptr->Vref = 0;
  ptr->Contref = 0;
  ptr->val = 0;
  ptr->Vvar = 0;
  ptr->From = ptrac;
  ptr->Ctrl = ptrac1;
  ptr->Prev = ptrp;
  ptr->Next = ptrn;
  return (ptr);
}
