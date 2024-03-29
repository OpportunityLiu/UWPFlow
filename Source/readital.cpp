/* Read AC data in ITALIAN format  */

#include "readdata.h"

bool ReadADDfile();
ACbusData *ACbusInList2(INDEX BusN, char *BusCode, INDEX N1, INDEX N2, bool flagAddBus);
ElementData *ElemInList2(ACbusData *From, ACbusData *To, INDEX N, char *Type, char *Ckt, bool flagAddElement);
VALUETYPE GetVelevels(char *BusName);
ElementList *TakeElemFromList(ElementList *ELptr, ElementData *Eptr);

/* ------- Global Variables ------ */
extern FILE *InputDataFile;

/* --------------- ACbusInList2  ------------------- */
ACbusData *ACbusInList2(INDEX BusN, char *BusCode, INDEX N1, INDEX N2, bool flagAddBus) {
  ACbusData *ptr, *ptrp, *ptrn;
  INDEX i;

  if ((N1 == 0 || N2 == 0) && flagAddBus) {
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
      if (BusN == ptr->Num || (BusN == 0 && !strncmp(ptr->Name, BusCode, strlen(BusCode))))
        return (ptr);
      ptrp = ptr;
      ptr = ptr->Next;
    }
    if (!flagAddBus)
      return (ptr);
    ptr = new ACbusData;
    ptrp->Next = ptr;
    ptrn = nullptr;
  }
  ptr->Num = 0;
  strcpy(ptr->Name, BusCode);
  ptr->N = 0;
  ptr->KV = 0;
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
  ptr->PgMax = 0;
  ptr->Smax = 0;
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
  ptr->Max = 0;
  ptr->Min = 0;
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

/* --------------- ElemInList2  ------------------- */
ElementData *ElemInList2(ACbusData *From, ACbusData *To, INDEX N, char *Type, char *Ckt, bool flagAddElement) {
  ElementData *ptr, *ptrp;
  ElementList *ELptr;

  if (N == 0) {
    dataPtr->Element = new ElementData;
    ptr = dataPtr->Element;
  } else {
    for (ELptr = From->Elem; ELptr != nullptr; ELptr = ELptr->Next) {
      ptr = ELptr->Eptr;
      if (((ptr->From == From && ptr->To == To) || (ptr->From == To && ptr->To == From)) && strpbrk(ptr->Type, Type) && !strcmp(ptr->Ckt, Ckt)) {
        if (flagAddElement)
          return (nullptr);
        else
          return (ptr);
      }
    }
    if (!flagAddElement)
      return (nullptr);
    ptr = new ElementData;
  }
  ptrp = dataPtr->Element;
  if (ptrp != ptr) {
    ptrp->Prev = ptr;
    ptr->Next = ptrp;
    dataPtr->Element = ptr;
  } else
    ptr->Next = nullptr;
  ptr->Prev = nullptr;
  ptr->From = From;
  ptr->To = To;
  From->Elem = (ElementList *)AddElemToList(From->Elem, ptr);
  To->Elem = (ElementList *)AddElemToList(To->Elem, ptr);
  strcpy(ptr->Ckt, Ckt);
  strcpy(ptr->Type, Type);
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
  return (ptr);
}

/* ---------------- ReadITALY ----------------------------- */
void ReadITALY()
/* Read Bus and Element data in WSCC format. */
{
  ACbusData *ACptr, *ACptrp, *ACptrs, *ACptrEstimate;
  ElementData *Eptr;
  AreaData *Aptr;
  char Line[BUFLEN], Code[6], Ckt[2], Type[2], Zone[3];
  VALUETYPE KV, KVp, KVs, R, X, G, B, B1, B2, Zb, Tap, Ang, Taps, Tmax, Tmin;
  VALUETYPE Vlevels[10], Vmax[10], Vmin[10], Vlmax, Vlmin;
  INDEX N, KVl;
  int i, j;
  bool EstimatePl = false, EstimateQl = false;

  for (i = 0; i <= 2; strcpy(dataPtr->Title[i], "\n"), i++)
    ;
  for (i = 0; i < 10; Vlevels[i] = Vmax[i] = Vmin[i] = 0, i++)
    ;
  Sn = 100.0;
  RealParameter('$', &Sn, 1.0, 100000000.0);
  ACptr = nullptr;
  i = 0;
  for (;;) { /* Reading Loop */
    Line[79] = ' ';
    if (fgets(Line, BUFLEN, InputDataFile) == nullptr) {
      ErrorHalt("Missing F card.");
      break;
    }
    LineNum++;

    /* --------------- Title and Comments ----------------------------- */
    if (Line[79] == 'C' || Line[79] == '*') {
      if (i == 0) {
        strcpy(dataPtr->Title[i], Line);
        i++;
      } else
        continue;
    }

    /* --------------- Voltage Levels and Limits --------------------- */
    else if (Line[79] == 'Z' && !strncmp(Line, "VNOM", 4)) {
      for (j = 0; j <= 9; j++)
        Vlevels[j + 1] = GetValue(Line, 7 + j * 5, 5, 0);
    } else if (Line[79] == 'Z' && !strncmp(Line, "VMIN", 4)) {
      for (j = 0; j <= 9; j++)
        Vmin[j + 1] = GetValue(Line, 7 + j * 5, 5, 0);
    } else if (Line[79] == 'Z' && !strncmp(Line, "VMAX", 4)) {
      for (j = 0; j <= 9; j++)
        Vmax[j + 1] = GetValue(Line, 7 + j * 5, 5, 0);
    }

    /* --------------- AC bus data -------------------------------- */

    /* Load */
    else if (Line[79] == 'N') {
      N = GetInt(Line, 1, 3);
      GetStr(Line, 4, 1, 1, Zone);
      KVl = GetInt(Line, 5, 1);
      sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
      ACptr = ACbusInList2(0, Code, Nac, 1, true);
      if (ACptr->V > 0) {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        ErrorHalt("The AC bus was previously defined (check N cards).");
      }
      if (ACptr->N == 0) {
        Nac++;
        ACptr->Num = ACptr->N = Nac;
        KV = Vlevels[KVl];
        sprintf(ACptr->Name, "%5s %6.0lf", Code, KV);
        Vlmax = Vmax[KVl];
        Vlmin = Vmin[KVl];
        if (KV == 0) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("Base voltage is zero.");
          KV = 1;
        } else if (Vlmax > Vlmin && Vlmax >= KV && Vlmin <= KV) {
          ACptr->Vlmax = Vlmax / KV;
          ACptr->Vlmin = Vlmin / KV;
        }
        ACptr->KV = KV;
      } else
        KV = ACptr->KV;
      strcpy(ACptr->Zone, Zone);
      strcpy(ACptr->Owner, Zone);
      ACptr->V = GetValue(Line, 47, 5, 0) / KV;
      if (ACptr->V <= 0)
        ACptr->V = 1;
      ACptr->VCont = ACptr->V;
      if (Line[20] == 'T') {
        EstimatePl = true;
        ACptrEstimate = ACptr;
        fprintf(stderr, "***Warning: The P load at bus %d %s will be estimated from the\n", ACptr->N, ACptr->Name);
        fprintf(stderr, "            oncoming element data for this bus.\n");
      } else {
        EstimatePl = false;
        ACptr->Pl = GetValue(Line, 14, 6, 0) / Sn;
      }
      if (Line[42] == 'T') {
        EstimateQl = true;
        ACptrEstimate = ACptr;
        fprintf(stderr, "***Warning: The Q load at bus %d %s will be estimated from the\n", ACptr->N, ACptr->Name);
        fprintf(stderr, "            oncoming element data for this bus.\n");
      } else {
        EstimateQl = false;
        ACptr->Ql = GetValue(Line, 36, 6, 0) / Sn;
      }
      ACptr->Vmax = GetValue(Line, 57, 3, 0) / KV;
      if (ACptr->Vmax <= 0)
        ACptr->Vmax = 10.;
      ACptr->Vmin = GetValue(Line, 53, 3, 0) / KV;
      if (ACptr->Vmin <= 0)
        ACptr->Vmin = 0.001;
      if (ACptr->Vmax <= ACptr->Vmin) {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        ErrorHalt("AC bus V limits are wrong: Vmin >= Vmax.");
      }
      if (Line[22] == 'T') {
        Nslack++;
        strcpy(ACptr->Type, "BS");
        strcpy(ACptr->cont, "V");
        AngSlack = 0;
      } else if (Line[44] == 'T' && strcmp(ACptr->Type, "BQ")) {
        Nvolt++;
        strcpy(ACptr->Type, "BQ");
        strcpy(ACptr->cont, "V");
      } else {
        ACptr->Cont = ACptr;
      }
    }

    /* Generator */
    else if (Line[79] == 'I' || Line[79] == 'E') {
      N = GetInt(Line, 1, 3);
      GetStr(Line, 4, 1, 1, Zone);
      KVl = GetInt(Line, 5, 1);
      sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
      ACptr = ACbusInList2(0, Code, Nac, 1, true);
      if (ACptr->N == 0) {
        Nac++;
        ACptr->Num = ACptr->N = Nac;
        KV = Vlevels[KVl];
        sprintf(ACptr->Name, "%5s %6.0lf", Code, KV);
        Vlmax = Vmax[KVl];
        Vlmin = Vmin[KVl];
        if (KV == 0) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("Base voltage is zero.");
          KV = 1;
        } else if (Vlmax > Vlmin && Vlmax >= KV && Vlmin <= KV) {
          ACptr->Vlmax = Vlmax / KV;
          ACptr->Vlmin = Vlmin / KV;
        }
        ACptr->KV = KV;
      }
      ACptr->Pg += -GetValue(Line, 14, 6, 0) / Sn;
      ACptr->PgMax += -GetValue(Line, 25, 5, 0) / Sn;
      if (ACptr->PgMax == 0)
        ACptr->PgMax = 99999999.;
      ACptr->Smax += -GetValue(Line, 30, 5, 0) / Sn;
      if (ACptr->Smax == 0)
        ACptr->Smax = 99999999.;
      if (ACptr->PgMax > ACptr->Smax) {
        ACptr->PgMax = ACptr->Smax;
        fprintf(stderr, "***Warning: Bus %d %s has PgMax > Smax.\n", ACptr->N, ACptr->Name);
        fprintf(stderr, "            PgMax will be set to Smax.\n");
      }
      if (ACptr->Pg > ACptr->PgMax) {
        ACptr->Pg = ACptr->PgMax;
        fprintf(stderr, "***Warning: Bus %d %s has Pg > PgMax.\n", ACptr->N, ACptr->Name);
        fprintf(stderr, "            Pg will be set to PgMax.\n");
      }
      ACptr->Qg += -GetValue(Line, 36, 6, 0) / Sn;
      ACptr->Qmax += -GetValue(Line, 47, 5, 0) / Sn;
      ACptr->Qmin += -GetValue(Line, 43, 4, 0) / Sn;
      if (ACptr->Qmax <= ACptr->Qmin) {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        ErrorHalt("AC bus Q limits are wrong: Qmin >= Qmax.");
      }
    }

    /* Synchronous Condenser */
    else if (Line[79] == 'S') {
      N = GetInt(Line, 1, 3);
      GetStr(Line, 4, 1, 1, Zone);
      KVl = GetInt(Line, 5, 1);
      sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
      ACptr = ACbusInList2(0, Code, Nac, 1, true);
      if (ACptr->N == 0) {
        Nac++;
        ACptr->Num = ACptr->N = Nac;
        KV = Vlevels[KVl];
        sprintf(ACptr->Name, "%5s %6.0lf", Code, KV);
        Vlmax = Vmax[KVl];
        Vlmin = Vmin[KVl];
        if (KV == 0) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("Base voltage is zero.");
          KV = 1;
        } else if (Vlmax > Vlmin && Vlmax >= KV && Vlmin <= KV) {
          ACptr->Vlmax = Vlmax / KV;
          ACptr->Vlmin = Vlmin / KV;
        }
        ACptr->KV = KV;
      }
      ACptr->Qg += -GetValue(Line, 36, 6, 0) / Sn;
      ACptr->Qmax += -GetValue(Line, 47, 5, 0) / Sn;
      ACptr->Qmin += -GetValue(Line, 43, 4, 0) / Sn;
      if (ACptr->Qmax <= ACptr->Qmin) {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        ErrorHalt("AC bus Q limits are wrong: Qmin >= Qmax.");
      }
    }

    /* Shunt Capacitor */
    else if (Line[79] == 'Q') {
      N = GetInt(Line, 1, 3);
      GetStr(Line, 4, 1, 1, Zone);
      KVl = GetInt(Line, 5, 1);
      sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
      ACptr = ACbusInList2(0, Code, Nac, 1, true);
      if (ACptr->N == 0) {
        Nac++;
        ACptr->Num = ACptr->N = Nac;
        KV = Vlevels[KVl];
        sprintf(ACptr->Name, "%5s %6.0lf", Code, KV);
        Vlmax = Vmax[KVl];
        Vlmin = Vmin[KVl];
        if (KV == 0) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("Base voltage is zero.");
          KV = 1;
        } else if (Vlmax > Vlmin && Vlmax >= KV && Vlmin <= KV) {
          ACptr->Vlmax = Vlmax / KV;
          ACptr->Vlmin = Vlmin / KV;
        }
        ACptr->KV = KV;
      } else
        KV = ACptr->KV;
      KV = GetValue(Line, 43, 4, 0) / KV;
      if (KV != 0)
        ACptr->B += (-GetValue(Line, 47, 5, 0) / Sn) / (KV * KV);
    }

    /* --------------- AC element data -------------------------------- */
    else if (Line[79] == 'L' || Line[79] == 'T') {
      N = GetInt(Line, 1, 3);
      GetStr(Line, 4, 1, 1, Zone);
      KVl = GetInt(Line, 5, 1);
      sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
      ACptr = ACbusInList2(0, Code, Nac, 1, true);
      if (ACptr->N == 0) {
        Nac++;
        ACptr->Num = ACptr->N = Nac;
        KV = Vlevels[KVl];
        sprintf(ACptr->Name, "%5s %6.0lf", Code, KV);
        Vlmax = Vmax[KVl];
        Vlmin = Vmin[KVl];
        if (KV == 0) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("Base voltage is zero.");
          KV = 1;
        } else if (Vlmax > Vlmin && Vlmax >= KV && Vlmin <= KV) {
          ACptr->Vlmax = Vlmax / KV;
          ACptr->Vlmin = Vlmin / KV;
        }
        ACptr->KV = KV;
      } else
        KV = ACptr->KV;
      N = GetInt(Line, 7, 3);
      GetStr(Line, 10, 1, 1, Zone);
      KVl = GetInt(Line, 11, 1);
      sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
      ACptrp = ACbusInList2(0, Code, Nac, 1, true);
      if (ACptrp->N == 0) {
        Nac++;
        ACptrp->Num = ACptrp->N = Nac;
        KVp = Vlevels[KVl];
        sprintf(ACptrp->Name, "%5s %6.0lf", Code, KVp);
        Vlmax = Vmax[KVl];
        Vlmin = Vmin[KVl];
        if (KVp == 0) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("Base voltage is zero.");
          KVp = 1;
        } else if (Vlmax > Vlmin && Vlmax >= KVp && Vlmin <= KVp) {
          ACptrp->Vlmax = Vlmax / KVp;
          ACptrp->Vlmin = Vlmin / KVp;
        }
        ACptrp->KV = KVp;
      } else
        KVp = ACptrp->KV;
      if (ACptr == ACptrp) {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        ErrorHalt("Both AC element buses are the same.");
      }
      if (EstimatePl) {
        if (Line[20] == 'T') {
          fprintf(stderr, "***Warning: The P load estimate at bus %d %s may be incorrect\n", ACptrEstimate->N, ACptrEstimate->Name);
          fprintf(stderr, "            due to lack of P flow information in "
                          "the corresponding element data:\n");
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        } else {
          if (ACptrEstimate == ACptr)
            ACptrEstimate->Pl += -GetValue(Line, 14, 6, 0) / Sn;
          else
            ACptrEstimate->Pl += GetValue(Line, 14, 6, 0) / Sn;
        }
      }
      if (EstimateQl) {
        if (Line[42] == 'T') {
          fprintf(stderr, "***Warning: The Q load estimate at bus %d %s may be incorrect\n", ACptrEstimate->N, ACptrEstimate->Name);
          fprintf(stderr, "            due to lack of Q flow information in "
                          "the corresponding element data:\n");
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        } else {
          if (ACptrEstimate == ACptr)
            ACptrEstimate->Ql += -GetValue(Line, 36, 6, 0) / Sn;
          else
            ACptrEstimate->Ql += GetValue(Line, 36, 6, 0) / Sn;
        }
      }
      if (KV > KVp) {
        KVs = KV;
        KV = KVp;
        KVp = KVs;
        ACptrs = ACptr;
        ACptr = ACptrp;
        ACptrp = ACptrs;
      }
      if (Line[79] == 'L')
        strcpy(Type, "L");
      else {
        Tmax = GetValue(Line, 67, 5, 2);
        Tmin = GetValue(Line, 72, 5, 2);
        /*
        if(Tmax>Tmin) strcpy(Type,"R");
        else strcpy(Type,"T");
        */
        strcpy(Type, "T");
      }
      GetStr(Line, 13, 1, 1, Ckt);
      if (!strcmp(Ckt, " "))
        strcpy(Ckt, "0");
      Eptr = ElemInList2(ACptr, ACptrp, NacEl, Type, Ckt, true);
      if (Eptr != nullptr) {
        Zb = KVp * KVp / Sn;
        R = GetValue(Line, 45, 7, 5) / Zb;
        X = GetValue(Line, 52, 8, 5) / Zb;
        if (fabs(R) < 1e-10 && fabs(X) < 1e-10) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("AC element is a short circuit. Try eliminating it.");
          G = B = 0;
        } else {
          G = R / (R * R + X * X);
          B = -X / (R * R + X * X);
        }
        B2 = B1 = GetValue(Line, 60, 7, 3) * Zb / 1000000. / 2;
        if (Line[79] == 'T') {
          Taps = KV / KVp;
          Tap = GetValue(Line, 29, 5, 0) * Taps;
        } else
          Tap = Taps = 1;
        Ang = 0;
        Eptr->Imax = GetValue(Line, 23, 5, 0) / (Sn * 1000 / (sqrt(3.0) * KVp));
        Eptr->G = G;
        Eptr->B = B;
        Eptr->B1 = B1;
        Eptr->B2 = B2;
        Eptr->Tap = Tap;
        Eptr->Taps = Taps;
        Eptr->Ang = Ang;
        GetStr(Line, 4, 1, 1, Eptr->Zone);
        strcpy(Eptr->Owner, Eptr->Zone);
        NacEl++;
        Eptr->Meter = ACptr;
        strcpy(Eptr->Type, Type);
        /*
        if (!strcmp(Type,"R")) {
          Eptr->Tmax=Tmax*Taps;
          Eptr->Tmin=Tmin*Taps;
          NregV++;
          if (!strcmp(ACptrp->Type,"B")) strcpy(ACptrp->Type,"BT");
          ACptrp->Reg=AddElemToList(ACptrp->Reg,Eptr);
          Eptr->Cont=ACptrp;
        }
        */
      }
    }

    /* -------------------- DC data -------------------------------- */
    else if (!strncmp(Line, "BD ", 3) || !strncmp(Line, "BZ ", 3) || !strncmp(Line, "LD ", 3))
      ReadEPRIdc(Line);

    /* FACTS */

    /* ---------------------- SVC data ------------------------ */
    else if (!strncmp(Line, "FS ", 3))
      ReadSVC(Line);

    /* ---------------------- TCSC data ------------------------ */
    else if (!strncmp(Line, "FC ", 3))
      ReadTCSC(Line);

    /* ---------------------- STATCOM data ------------------------ */
    else if (!strncmp(Line, "FT ", 3))
      ReadSTATCOM(Line);

    /* END FACTS */

    else if (Line[79] == 'F')
      break;
    else {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      fprintf(stderr, "***Warning: The program will ignore this line.\n");
    }
  }
  fclose(InputDataFile);
  MaxIter = 50;

  if (!ReadADDfile())
    for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
      Aptr = (AreaData *)AreaInList(0, ACptr->Zone, Narea);
      if (Aptr->N == 0) {
        Narea++;
        Aptr->N = Narea;
        strcpy(Aptr->Zone[1], ACptr->Zone);
      }
      ACptr->Area = Aptr;
      if (strpbrk(ACptr->Type, "S"))
        Aptr->Slack = Aptr->BSptr = ACptr;
      else if (strpbrk(ACptr->Type, "Q") && Aptr->Slack == nullptr)
        Aptr->Slack = Aptr->BSptr = ACptr;
    }
}

/* --------------- ReadADDfile ------------------- */
bool ReadADDfile()
/* Read ADD COLAS file with new voltage information, areas, and
   generator and load data for voltage collapse studies. */
{
  FILE *InputFile;
  char *Name, Line[BUFLEN], Code[20], Ckt[2], Zone[5];
  ACbusData *ACptr, *ACptrp, *ACptrs;
  AreaData *Aptr;
  ElementList *ELptr;
  ElementData *Eptr;
  VALUETYPE KV, KVp, KVs, KVmax, KVmin, val, Tap, Taps, Tmax, Tmin, Q, Qmax, Qmin, Pn, Qn, Sum = 0;
  bool flagAreas = false, flagPrint = true, flagScards = false;
  INDEX N, NJcard = 0, N2SVCarea = 0, KVl;

  Name = NameParameter('6');
  if (!NullName(Name) && (InputFile = OpenInput(Name)) != nullptr) {
    LineNum = 0;
    for (;;) {
      if (fgets(Line, BUFLEN, InputFile) == nullptr)
        break;
      LineNum++;

      /* --------------- Comments ----------------------------- */
      if (Line[79] == 'C')
        continue;

      /* --------------- AC buses ----------------------------- */
      else if (Line[79] == 'K') {
        N = GetInt(Line, 1, 3);
        GetStr(Line, 4, 1, 1, Zone);
        KVl = GetInt(Line, 5, 1);
        sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
        ACptr = ACbusInList2(0, Code, Nac, 1, false);
        if (ACptr == nullptr) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          fprintf(stderr, "***Warning: This bus has not been defined on the "
                          "main input data file.\n");
          fprintf(stderr, "            This line in the ADD file will be ignored.\n");
        } else {
          KV = GetValue(Line, 6, 7, 0);
          if (KV == 0) {
            fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
            fprintf(stderr, "***Warning: This bus has a zero base bus voltage.\n");
            fprintf(stderr, "            This line in the ADD file will be ignored.\n");
          } else {
            sprintf(ACptr->Name, "%5s %6.0lf", Code, KV);
            KVp = ACptr->KV;
            val = KVp / KV;
            ACptr->VCont = ACptr->V = ACptr->V * val;
            ACptr->KV = KV;
            KVmin = GetValue(Line, 13, 7, 0);
            KVmax = GetValue(Line, 20, 7, 0);
            if (KVmax > KVmin && KVmax >= KV && KVmin <= KV) {
              ACptr->Vlmax = KVmax / KV;
              ACptr->Vlmin = KVmin / KV;
            }
            val = KV * KV / (KVp * KVp);
            ACptr->B = ACptr->B * val;
            for (ELptr = ACptr->Elem; ELptr != nullptr; ELptr = ELptr->Next) {
              Eptr = ELptr->Eptr;
              if (Eptr->To == ACptr) {
                val = KV * KV / (KVp * KVp);
                Eptr->G = Eptr->G * val;
                Eptr->B = Eptr->B * val;
                Eptr->B1 = Eptr->B1 * val;
                Eptr->B2 = Eptr->B2 * val;
                val = KVp / KV;
                Eptr->Imax = Eptr->Imax / val;
                Eptr->Tap = Eptr->Tap * val;
                Eptr->Taps = Eptr->Taps * val;
                Eptr->Tmax = Eptr->Tmax * val;
                Eptr->Tmin = Eptr->Tmin * val;
              } else {
                val = KV / KVp;
                Eptr->Tap = Eptr->Tap * val;
                Eptr->Taps = Eptr->Taps * val;
                Eptr->Tmax = Eptr->Tmax * val;
                Eptr->Tmin = Eptr->Tmin * val;
              }
            }
          }
          GetStr(Line, 50, 1, 1, Code);
          if (!strcmp(Code, "T")) {
            if (!strcmp(ACptr->Type, "BQ")) {
              ACptr->Qmax = 99999999.;
              ACptr->Qmin = -99999999.;
            } else {
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              fprintf(stderr, "***Warning: This bus was not defined as a PV "
                              "bus in the main input data file; \n");
              fprintf(stderr, "            hence, its Q-limits cannot be "
                              "released (check field 50).\n");
              fprintf(stderr, "            This line in the ADD file will be ignored.\n");
            }
          }
          GetStr(Line, 58, 1, 1, Code);
          if (strcmp(Code, "T"))
            ACptr->DPG = ACptr->DPg = ACptr->Pg;
          if (strpbrk(ACptr->Type, "Q,S")) {
            GetStr(Line, 66, 2, 2, Zone);
            strcpy(ACptr->Zone, Zone);
            strcpy(ACptr->Owner, Zone);
          } else {
            strcpy(ACptr->Zone, " 0");
            strcpy(ACptr->Owner, " 0");
          }
          ACptr->Nc = GetInt(Line, 68, 2);
          GetStr(Line, 70, 1, 1, Code);
          if (strpbrk(Code, "P") && !strcmp(ACptr->Type, "B")) {
            strcpy(ACptr->Type, "BC");
            ACptr->VCont = ACptr->V;
            flag2Vcontrol = true;
          } else if (strpbrk(Code, "V,S") && strpbrk(ACptr->Type, "Q,S")) {
            if (!strcmp(ACptr->Type, "BQ"))
              strcpy(ACptr->Type, "BG");
            else
              strcpy(ACptr->Type, "BGS");
            flag2Vcontrol = true;
          }
          ACptr->Qn = ACptr->Ql;
          ACptr->b = GetValue(Line, 71, 7, 0);
          ACptr->Pn = ACptr->Pl;
          ACptr->a = GetValue(Line, 81, 7, 0);
          if (ACptr->a != 0 || ACptr->b != 0)
            flagVloads = true;
        }
      }
      /* --------------- Generators --------------------------- */
      else if (Line[79] == 'G' || Line[79] == 'H') {
        if (flagPrint) {
          fprintf(stderr, "***Warning: The program ignores the G and H cards, "
                          "as capability curves\n");
          fprintf(stderr, "            are implemented here through direct "
                          "limits on Sg, Pg, Qg,\n");
          fprintf(stderr, "            Ef and Ia (see -3 option).\n");
          flagPrint = false;
        }
      }

      /* --------------- Transformers  ------------------------ */
      else if (Line[79] == 'R') {
        N = GetInt(Line, 1, 3);
        GetStr(Line, 4, 1, 1, Zone);
        KVl = GetInt(Line, 5, 1);
        sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
        ACptr = ACbusInList2(0, Code, Nac, 1, false);
        N = GetInt(Line, 7, 3);
        GetStr(Line, 10, 1, 1, Zone);
        KVl = GetInt(Line, 11, 1);
        sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
        ACptrp = ACbusInList2(0, Code, Nac, 1, false);
        GetStr(Line, 13, 1, 1, Ckt);
        if (!strcmp(Ckt, " "))
          strcpy(Ckt, "0");
        Eptr = ElemInList2(ACptr, ACptrp, NacEl, "R,T", Ckt, false);
        if (Eptr != nullptr) {
          KV = ACptr->KV;
          KVp = ACptrp->KV;
          if (KV > KVp)
            Taps = KVp / KV;
          else
            Taps = KV / KVp;
          Eptr->Taps = Taps;
          Tap = GetValue(Line, 14, 8, 0) * Taps;
          if (Tap > 0)
            Eptr->Tap = Tap;
          else {
            fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
            fprintf(stderr, "***Warning: This transformer data has a negative "
                            "tap setting.\n");
            fprintf(stderr, "            The tap will be assumed to be equal "
                            "to the nominal\n");
            fprintf(stderr, "            voltage ratio.\n");
            Eptr->Tap = Tap = 1;
          }
          Tmin = GetValue(Line, 24, 8, 0) * Taps;
          Tmax = GetValue(Line, 34, 8, 0) * Taps;
          if (Tmax > Tmin && Tmax > Tap && Tap > Tmin) {
            Eptr->Tmax = Tmax;
            Eptr->Tmin = Tmin;
          } else {
            fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
            fprintf(stderr, "***Warning: The max. and min. tap settings of "
                            "this transformer are incorrect.\n");
            fprintf(stderr, "            These limits will be assumed to be "
                            "equal to 1.1 and 0.9 p.u.,\n");
            fprintf(stderr, "            respectively.\n");
            Eptr->Tmax = 1.1;
            Eptr->Tmin = 0.9;
          }
          N = GetInt(Line, 42, 2);
          KVs = GetValue(Line, 45, 8, 0);
          if (!strcmp(Eptr->Type, "R")) {
            ACptrs = Eptr->Cont;
            if (N == 0) {
              NregV--;
              strcpy(Eptr->Type, "T");
              ACptrs->Reg = TakeElemFromList(ACptrs->Reg, Eptr);
              Eptr->Cont = nullptr;
            } else {
              if (N > 0) {
                if (KV > KVp) {
                  ACptrs = ACptr;
                  ACptrs->V = KVs / KV;
                } else {
                  ACptrs = ACptrp;
                  ACptrs->V = KVs / KVp;
                }
              } else {
                if (KV < KVp) {
                  ACptrs = ACptr;
                  ACptrs->V = KVs / KV;
                } else {
                  ACptrs = ACptrp;
                  ACptrs->V = KVs / KVp;
                }
              }
              if (ACptrs != Eptr->Cont) {
                Eptr->Cont->Reg = TakeElemFromList(Eptr->Cont->Reg, Eptr);
                ACptrs->Reg = AddElemToList(ACptrs->Reg, Eptr);
                Eptr->Cont = ACptrs;
              }
            }
          } else if (N != 0) {
            NregV++;
            strcpy(Eptr->Type, "R");
            if (N > 0) {
              if (KV > KVp)
                ACptrs = ACptr;
              else
                ACptrs = ACptrp;
            } else {
              if (KV < KVp)
                ACptrs = ACptr;
              else
                ACptrs = ACptrp;
            }
            if (!strcmp(ACptrs->Type, "B"))
              strcpy(ACptrs->Type, "BT");
            ACptrs->Reg = AddElemToList(ACptrs->Reg, Eptr);
            Eptr->Cont = ACptrs;
          }
        } else {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          fprintf(stderr, "***Warning: This transformer has not been defined "
                          "in the main data file.\n");
          fprintf(stderr, "            This line in the ADD file will be ignored.\n");
        }
      }

      /* --------------- Shunt Compensation  ------------------ */
      else if (Line[79] == 'B') {
        N = GetInt(Line, 1, 3);
        GetStr(Line, 4, 1, 1, Zone);
        KVl = GetInt(Line, 5, 1);
        sprintf(Code, "%3d%1s%1d", N, Zone, KVl);
        ACptr = ACbusInList2(0, Code, Nac, 1, false);
        if (ACptr == nullptr) {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          fprintf(stderr, "***Warning: This bus has not been defined on the "
                          "main input data file.\n");
          fprintf(stderr, "            This line in the ADD file will be ignored.\n");
        } else {
          Q = GetValue(Line, 15, 8, 0) / Sn;
          GetStr(Line, 31, 1, 1, Code);
          if (!strcmp(Code, "T")) {
            if (!strcmp(ACptr->Type, "B")) {
              ACptr->Qg = ACptr->Qg + Q;
              Qmin = GetValue(Line, 7, 8, 0) / Sn;
              Qmax = GetValue(Line, 23, 8, 0) / Sn;
              if (Qmax <= Qmin || Qmax < Q || Q < Qmin) {
                fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
                fprintf(stderr, "***Warning: The Q limits are inconsistent on "
                                "this bus and hence will\n");
                fprintf(stderr, "            be given arbitrarily large values "
                                "(~+/-10^8 p.u.).\n");
                Qmax = 99999999.;
                Qmin = -99999999.;
              }
              ACptr->Qmax = Qmax;
              ACptr->Qmin = Qmin;
              GetStr(Line, 33, 1, 1, Code);
              if (!strcmp(Code, "T")) {
                strcpy(ACptr->Type, "BZ");
                NZvolt++;
              } else {
                strcpy(ACptr->Type, "BQ");
                Nvolt++;
              }
              ACptr->Cont = nullptr;
              strcpy(ACptr->cont, "V");
            } else {
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              fprintf(stderr, "***Warning: This bus cannot be defined as a "
                              "voltage controlled bus, as it\n");
              fprintf(stderr, "            has already been defined as a "
                              "controlled bus in the main file.\n");
            }
          } else {
            GetStr(Line, 33, 1, 1, Code);
            if (!strcmp(Code, "T"))
              ACptr->B = Q / (ACptr->V * ACptr->V);
            else
              ACptr->Qg = ACptr->Qg + Q;
          }
          GetStr(Line, 35, 2, 2, Zone);
          strcpy(ACptr->Zone, Zone);
          strcpy(ACptr->Owner, Zone);
        }
      }

      /* --------------- Areas  ------------------------------- */
      else if (Line[79] == 'J') {
        NJcard++;
        sprintf(Zone, "%2d", NJcard);
        GetStr(Line, 1, 16, 16, Code);
        val = GetValue(Line, 34, 4, 0) / 100.0;
        for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next)
          if (ACptr->Area == nullptr) {
            if (!strcmp(ACptr->Zone, Zone)) {
              Aptr = (AreaData *)AreaInList(0, Code, Narea);
              if (Aptr->N == 0) {
                Narea++;
                Aptr->N = Narea;
                strcpy(Aptr->Zone[1], Zone);
                flagAreas = true;
              }
              ACptr->Area = Aptr;
              if (strpbrk(ACptr->Type, "S"))
                Aptr->Slack = Aptr->BSptr = ACptr;
              else if (strpbrk(ACptr->Type, "Q") && Aptr->Slack == nullptr)
                Aptr->Slack = Aptr->BSptr = ACptr;
              Pn = val * ACptr->Pl;
              Qn = val * ACptr->Ql;
              ACptr->Pnl = Pn;
              ACptr->Qnl = Qn;
              Sum += Pn + Qn;
            }
          }
      }

      /* --------------- Secondary Voltage Control   ----------- */
      else if (Line[79] == 'S') {
        flagScards = true;
        N2SVCarea++;
        ACptrp = nullptr;
        GetStr(Line, 1, 16, 16, Code);
        for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
          if (ACptr->Nc == N2SVCarea && strpbrk(ACptr->Type, "C")) {
            if (ACptrp == nullptr) {
              ACptrp = ACptr;
              ACptrp->Cont = nullptr;
              ACptrp->VCont = ACptr->V;
            } else {
              fprintf(stderr,
                      "***Warning: Secondary voltage control area %s has more "
                      "than 1 pilot node.\n",
                      Code);
              fprintf(stderr, "            The pilot node %d %s will be ignored.\n", ACptr->N, ACptr->Name);
              strcpy(ACptr->Type, "B");
            }
          }
        }
        if (ACptrp != nullptr) {
          for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
            if (ACptr->Nc == N2SVCarea && strpbrk(ACptr->Type, "G")) {
              ACptr->Cont = ACptrp;
              ACptrp->Kbg++;
              ACptrp->Kbg1 = ACptrp->Kbg1 + ACptr->Qmax;
              ACptrp->Kbg2 = ACptrp->Kbg2 + ACptr->Qmin;
            }
          }
        } else {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          fprintf(stderr,
                  "***Warning: The secondary voltage control area does not "
                  "have a pilot node.\n",
                  Code);
          fprintf(stderr, "            This area in the ADD file will be ignored.\n");
        }
      }

      /* --------------- Ignore data  ------------------------- */
      else {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        fprintf(stderr, "***Warning: This line in the ADD file will be ignored.\n");
      }
    }

    fclose(InputFile);
  }

  if (flag2Vcontrol && !flagScards) {
    for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
      if (!strcmp(ACptr->Type, "BC"))
        strcpy(ACptr->Type, "B");
      else if (!strcmp(ACptr->Type, "BG"))
        strcpy(ACptr->Type, "BQ");
      else if (!strcmp(ACptr->Type, "BGS"))
        strcpy(ACptr->Type, "BS");
    }
    flag2Vcontrol = false;
  }

  if (Sum != 0)
    flagKdirection = true;
  else {
    if (ExistParameter('v')) {
      fprintf(stderr, "ERROR: The -v option will yield a singular Jacobian in "
                      "voltage collapse\n");
      fprintf(stderr, "       studies since Pnl, Qnl, Pzl, and Qzl are zero in "
                      "all load buses.\n");
      InputError = true;
    } else if (ExistParameter('L')) {
      fprintf(stderr, "***Warning: The loading factor lambda will not yield "
                      "different results\n");
      fprintf(stderr, "            from the base case since Pnl, Qnl, Pzl, and "
                      "Qzl are zero\n");
      fprintf(stderr, "            in all load buses.\n");
    } else if ((ExistParameter('H') || ExistParameter('c'))) {
      fprintf(stderr, "ERROR: The Homotopy Continuation Method will not yield "
                      "different results\n");
      fprintf(stderr, "       from the base case since Pnl, Qnl, Pzl, and Qzl "
                      "are zero in all\n");
      fprintf(stderr, "       load buses.\n");
      InputError = true;
    } else if (ExistParameter('C')) {
      fprintf(stderr, "ERROR: The Point of Collapse Method will not yield "
                      "different results\n");
      fprintf(stderr, "       from the base case since Pnl, Qnl, Pzl, and Qzl "
                      "are zero in\n");
      fprintf(stderr, "       all load buses.\n");
      InputError = true;
    }
  }

  if (flagAreas) {
    for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
      if (ACptr->Area == nullptr && !strcmp(ACptr->Zone, " 0")) {
        Aptr = (AreaData *)AreaInList(0, "REST", Narea);
        if (Aptr->N == 0) {
          Narea++;
          Aptr->N = Narea;
          strcpy(Aptr->Zone[1], " 0");
        }
        ACptr->Area = Aptr;
        if (strpbrk(ACptr->Type, "S"))
          Aptr->Slack = Aptr->BSptr = ACptr;
        else if (strpbrk(ACptr->Type, "Q") && Aptr->Slack == nullptr)
          Aptr->Slack = Aptr->BSptr = ACptr;
      } else if (ACptr->Area == nullptr && ACptr->SVC != nullptr) {
        Eptr = ACptr->Elem->Eptr;
        if (Eptr->From == ACptr)
          ACptrp = Eptr->To;
        else
          ACptrp = Eptr->From;
        ACptr->Area = ACptrp->Area;
      }
    }
  }

  return (flagAreas);
}

/* --------------- TakeElemFromList  ------------------- */
ElementList *TakeElemFromList(ElementList *ELptr, ElementData *Eptr) {
  ElementList *ptr, *prevptr, *firstptr, *nextptr;
  bool flag;

  firstptr = prevptr = ptr = ELptr;
  while (ptr != nullptr) {
    nextptr = ptr->Next;
    if (ptr->Eptr == Eptr) {
      if (ptr == firstptr)
        flag = true;
      else
        flag = false;
      delete ptr;
      if (flag)
        return (nextptr);
      else {
        prevptr->Next = nextptr;
        return (firstptr);
      }
    }
    prevptr = ptr;
    ptr = nextptr;
  }
  return (firstptr);
}
