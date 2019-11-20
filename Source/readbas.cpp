/*  Read data: Basic routines. */

#include "readdata.h"

/* --------- Global Input File --------- */
extern FILE *InputDataFile;

/* --------- GetStr ------------------- */
char *GetStr(const char *ptr, int Pos, int Leng, int Tot, char *str) {
  int i, count;

  for (i = 1; i < Pos; i++) {
    ptr++;
    if (*ptr == '\0' || *ptr == '\n') {
      for (i = 1; i <= Tot; str[i - 1] = ' ', i++)
        ;
      str[Tot] = '\0';
      return (str);
    }
  }
  count = 1;
  for (i = 1; i <= Leng && *ptr != '\0' && *ptr != '\n'; str[i - 1] = *ptr, ptr++, i++, count++)
    ;
  for (i = count; i <= Tot; str[i - 1] = ' ', i++)
    ;
  str[Tot] = '\0';
  return (str);
}

/* ---------------------- GetValue ------------------- */
VALUETYPE GetValue(const char *ptr, int Pos, int Leng, int Dec) {
  int i, count;
  bool flag = false, flagp = false;
  char str[11];
  VALUETYPE val;

  for (i = 1; i < Pos; i++) {
    ptr++;
    if (*ptr == '\n') {
      return (0);
    }
  }
  count = 1;
  for (i = 1; i <= Leng && *ptr != '\n'; ptr++, i++, count++) {
    if (isdigit(*ptr) || *ptr == '.' || *ptr == '-' || *ptr == '+') {
      str[i - 1] = *ptr;
      flagp = true;
    } else if (*ptr == ' ') {
      if (flagp)
        str[i - 1] = '0';
      else
        str[i - 1] = ' ';
    } else
      return (0);
    if (*ptr == '.')
      flag = true;
  }
  for (i = count; i <= Leng; str[i - 1] = '0', i++)
    ;
  str[Leng] = '\0';
  val = atof(str);
  if (Dec != 0 && flag == false)
    for (i = 1; i <= Dec; i++, val = val / 10)
      ;
  return (val);
}

/* ---------------------- GetInt ------------------- */
INDEX GetInt(const char *ptr, int Pos, int Leng) {
  int i, count;
  bool flagp = false;
  char str[11];
  INDEX val;

  for (i = 1; i < Pos; i++) {
    ptr++;
    if (*ptr == '\n') {
      return (0);
    }
  }
  count = 1;
  for (i = 1; i <= Leng && *ptr != '\n'; ptr++, i++, count++) {
    if (isdigit(*ptr) || *ptr == '-' || *ptr == '+') {
      str[i - 1] = *ptr;
      flagp = true;
    } else if (*ptr == ' ') {
      if (flagp)
        str[i - 1] = '0';
      else
        str[i - 1] = ' ';
    } else
      return (0);
  }
  for (i = count; i <= Leng; str[i - 1] = '0', i++)
    ;
  str[Leng] = '\0';
  val = (INDEX)atoi(str);
  return (val);
}

/* --------------- ACbusInList  ------------------- */
ACbusData *ACbusInList(INDEX BusN, char *BusName, VALUETYPE V, INDEX N1, INDEX N2) {
  ACbusData *ptr, *ptrp, *ptrn;
  int i;

  if (N1 == 0 || N2 == 0) {
    if (N1 == 0)
      ptrn = nullptr;
    else
      ptrn = dataPtr->ACbus;
    dataPtr->ACbus = new ACbusData;
    ptr = dataPtr->ACbus;
    if (ptrn != nullptr)
      ptrn->Prev = ptr;
    ptrp = nullptr;
  } else {
    ptr = dataPtr->ACbus;
    while (ptr != nullptr) {
      if (BusN == ptr->Num || (BusN == 0 && !strcmp(ptr->Name, BusName)))
        return (ptr);
      ptrp = ptr;
      ptr = ptr->Next;
    }
    ptr = new ACbusData;
    ptrp->Next = ptr;
    ptrn = nullptr;
  }
  ptr->Num = 0;
  strcpy(ptr->Name, BusName);
  ptr->N = 0;
  ptr->KV = V;
  strcpy(ptr->Type, "B");
  strcpy(ptr->Zone, "");
  strcpy(ptr->Owner, "");
  strcpy(ptr->cont, "");
  ptr->Area = nullptr;
  ptr->Ncont = 0;
  ptr->Reg = nullptr;
  ptr->Elem = nullptr;
  ptr->Gen = nullptr;
  ptr->V = 0;
  ptr->VCont = 0;
  ptr->Ang = 1000.;
  ptr->Pg = 0;
  ptr->Qg = 0;
  ptr->Pl = 0;
  ptr->Ql = 0;
  ptr->G = 0;
  ptr->B = 0;
  ptr->Bz = 0;
  ptr->step = 0;
  ptr->steps = 0;
  for (i = 0; i <= 72; ptr->Bx[i] = 0, i++)
    ;
  ptr->PgMax = 99999999.;
  ptr->Smax = 99999999.;
  ptr->flagPgMax = 0;
  ptr->DPg = 0;
  ptr->Kg = 0;
  ptr->Pz = 0;
  ptr->Qz = 0;
  ptr->Pn = 0;
  ptr->Qn = 0;
  ptr->Pzl = 0;
  ptr->Qzl = 0;
  ptr->Pnl = 0;
  ptr->Qnl = 0;
  ptr->a = 0;
  ptr->b = 0;
  ptr->Qmax = 0;
  ptr->Qmin = 0;
  ptr->Vmax = 0;
  ptr->Vmin = 0;
  ptr->Vlmax = 0;
  ptr->Vlmin = 0;
  ptr->CheckVlimits = true;
  ptr->Qr = 0;
  ptr->Kbg = 0;
  ptr->Kbg1 = 0;
  ptr->Kbg2 = 0;
  ptr->val = 0;
  ptr->valp = 0;
  ptr->vals = 0;
  ptr->valt = 0;
  ptr->DC = nullptr;
  ptr->SVC = nullptr;     /* FACTS */
  ptr->TCSC = nullptr;    /* FACTS */
  ptr->STATCOM = nullptr; /* FACTS */
  ptr->Cont = nullptr;
  ptr->ContBus = nullptr;
  ptr->Next = ptrn;
  ptr->Prev = ptrp;
  return (ptr);
}

/* --------------- AreaInList  ------------------- */
AreaData *AreaInList(INDEX i, char *Name, INDEX N) {
  AreaData *ptr, *prevptr;
  int j;

  if (N == 0) {
    dataPtr->Area = new AreaData;
    ptr = dataPtr->Area;
  } else {
    ptr = dataPtr->Area;
    while (ptr != nullptr) {
      if (i == ptr->N || (i == 0 && !strcmp(ptr->Name, Name)))
        return (ptr);
      prevptr = ptr;
      ptr = ptr->Next;
    }
    ptr = new AreaData;
    prevptr->Next = ptr;
  }
  ptr->N = 0;
  ptr->i = 0;
  strcpy(ptr->Name, Name);
  for (j = 1; j <= 11; j++)
    strcpy(ptr->Zone[j], "");
  ptr->P = 0;
  ptr->SPg = 0;
  ptr->Slack = nullptr;
  ptr->BSptr = nullptr;
  ptr->Elem = nullptr;
  ptr->AC = nullptr;
  ptr->DC = nullptr;
  ptr->Next = nullptr;
  return (ptr);
}

/* --------------- AddElemToList  ------------------- */
ElementList *AddElemToList(ElementList *ELptr, ElementData *Eptr) {
  ElementList *ptr, *prevptr;

  prevptr = ELptr;
  ELptr = new ElementList;
  ptr = ELptr;
  ptr->Eptr = Eptr;
  ptr->Next = prevptr;
  return (ptr);
}

/* --------------- ElemInList  ------------------- */
ElementData *ElemInList(ACbusData *From, ACbusData *To, INDEX N1, INDEX N2, char *Type, char *Ckt) {
  ElementData *ptr, *ptrn;
  ElementList *ELptr;

  if (N1 == 0 || N2 == 0) {
    if (N1 != 0)
      ptrn = dataPtr->Element;
    else
      ptrn = nullptr;
    dataPtr->Element = new ElementData;
    ptr = dataPtr->Element;
    if (ptrn != nullptr)
      ptrn->Prev = ptr;
    ptr->From = From;
    ptr->To = To;
    From->Elem = (ElementList *)AddElemToList(From->Elem, ptr);
    To->Elem = (ElementList *)AddElemToList(To->Elem, ptr);
    strcpy(ptr->Ckt, Ckt);
    strcpy(ptr->Type, "");
    strcpy(ptr->Zone, "");
    strcpy(ptr->Owner, "");
    ptr->Area = nullptr;
    ptr->Meter = nullptr;
    ptr->Sec = 0;
    ptr->G = 0;
    ptr->B = 0;
    ptr->G1 = 0;
    ptr->B1 = 0;
    ptr->G2 = 0;
    ptr->B2 = 0;
    ptr->Tap = 1;
    ptr->Taps = 1;
    ptr->Ang = 0;
    ptr->Cont = nullptr;
    ptr->Ncont = 0;
    ptr->Cvar = 0;
    strcpy(ptr->Ctype, "");
    ptr->Tmin = 0;
    ptr->Tmax = 0;
    ptr->Min = 0;
    ptr->Max = 0;
    ptr->Imax = 0;
    ptr->CheckIlimits = true;
    ptr->Next = ptrn;
    ptr->Prev = nullptr;
    return (ptr);
  } else {
    for (ELptr = From->Elem; ELptr != nullptr; ELptr = ELptr->Next) {
      ptr = ELptr->Eptr;
      if (strpbrk(ptr->Type, Type) && (!strcmp(ptr->Ckt, Ckt) || !strcmp(Ckt, " ") || !strcmp(ptr->Ckt, " ")) &&
          ((ptr->From == From && ptr->To == To) || (ptr->To == From && ptr->From == To))) {
        if (strcmp(Type, "R") || !strcmp(ptr->Zone, "")) {
          if (!strcmp(ptr->Ckt, " "))
            strcpy(ptr->Ckt, Ckt);
          return (ptr);
        }
      }
    }
    return (nullptr);
  }
}

/* --------------- DCbusInList  ------------------- */
DCbusData *DCbusInList(char *BusName, INDEX N) {
  DCbusData *ptr, *prevptr;

  if (N == 0) {
    dataPtr->DCbus = new DCbusData;
    ptr = dataPtr->DCbus;
  } else {
    ptr = dataPtr->DCbus;
    while (ptr != nullptr) {
      if (!strcmp(ptr->Name, BusName))
        return (ptr);
      prevptr = ptr;
      ptr = ptr->Next;
    }
    ptr = new DCbusData;
    prevptr->Next = ptr;
  }
  strcpy(ptr->Name, BusName);
  ptr->N = 0;
  strcpy(ptr->Type, "");
  strcpy(ptr->Cont1, "");
  strcpy(ptr->Cont2, "");
  strcpy(ptr->Zone, "");
  ptr->Meter = nullptr;
  ptr->Area = nullptr;
  ptr->Xc = 0;
  ptr->Nbr = 0;
  ptr->Ntrf = 0;
  ptr->MVA = 0;
  ptr->Vd = 0;
  ptr->VdN = 0;
  ptr->Id = 0;
  ptr->P = 0;
  ptr->Q = 0;
  ptr->Alfa = 0;
  ptr->AlfaN = 0;
  ptr->Gamma = 0;
  ptr->AlfaMin = 0;
  ptr->AlfaMax = 0;
  ptr->GammaMin = 0;
  ptr->Tap = 0;
  ptr->TapMin = 0;
  ptr->TapMax = 0;
  ptr->Vn = 0;
  ptr->Rd = 0;
  ptr->AC = nullptr;
  ptr->To = nullptr;
  ptr->Next = nullptr;
  return (ptr);
}
