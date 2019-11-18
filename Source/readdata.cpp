/*  Read input data.
    Transform the input data into the data structures needed by
    the power flow, i.e., AC, DC and FACTS bus data and interconnecting
    elements (lines, transformers with fixed taps, and PI equivalents).

    Main. */

#include "readdata.h"

/* ------- Global Variables ------ */
VALUETYPE AngSlack;
INDEX NdcEl, LineNum;
FILE *InputDataFile;

/* ------------------- ErrorDetect -------------------------- */
void ErrorDetect(void)
/* Detect inconsistencies in input data. */
{
  ACbusData *ACptr, *ACptrp, *ACptrs;
  AClist *ptrp, *ptrs;
  DClist *ptr;
  DCbusData *DCptr, *DCptrp;
  ElementData *Eptr, *Eptrp;
  ElementList *ELptr, *ELptrp;
  AreaData *Aptr, *Aptrp;
  SVCbusData *SVCptr;         /* FACTS */
  STATCOMbusData *STATCOMptr; /* FACTS */
  VALUETYPE KVs;
  bool flag, ieee;
  int i, j;

  if (Nac == 0 || (NacEl + NdcEl) == 0)
    ErrorHalt("No AC buses and/or elements in input data.");
  if (Ndc == 0)
    dataPtr->DCbus = nullptr;
  if (NacEl == 0)
    dataPtr->Element = nullptr;
  if (Nsvc == 0)
    dataPtr->SVCbus = nullptr; /* FACTS */
  if (Ntcsc == 0)
    dataPtr->TCSCbus = nullptr; /* FACTS */
  if (Nstatcom == 0)
    dataPtr->STATCOMbus = nullptr; /* FACTS */
  if (Narea < 2) {
    dataPtr->Area = nullptr;
    Narea = 0;
  }
  if (Nslack == 0)
    ErrorHalt("No angle reference Bus for AC system.");

  /* ------------------- AC buses ------------------------ */
  ACptr = dataPtr->ACbus;
  ieee = ExistParameter('I');
  while (ACptr != nullptr) {
    if (ACptr->V == 0) {
      fprintf(stderr, "ERROR: AC/DC bus %d %s has not been defined.\n",
              ACptr->Num, ACptr->Name);
      fprintf(stderr, "       Check the bus input cards.\n");
      InputError = true;
    }
    /* FACTS */
    if (ACptr->Elem == nullptr && ACptr->DC == nullptr &&
        ACptr->TCSC == nullptr) {
      fprintf(stderr, "ERROR: AC bus %d %s is isolated.\n", ACptr->Num,
              ACptr->Name);
      fprintf(stderr, "       Check the AC/DC/FACTS input cards.\n");
      InputError = true;
    }
    /* END OF FACTS */
    if (strpbrk(ACptr->Type, "G")) {
      ptrp = new AClist;
      ACptr->ContBus = ptrp;
      ACptr->ContBus->AC = ACptr->Cont;
      ACptr->ContBus->Next = ACptr->ContBus->Prev = nullptr;
      if (ACptr->Cont != nullptr) {
        if (!strpbrk(ACptr->Cont->Type, "C")) {
          fprintf(stderr, "ERROR: The voltage controlled bus %d %s of PV bus\n",
                  ACptr->Cont->Num, ACptr->Cont->Name);
          fprintf(
              stderr,
              "       %d %s, is not a PQ bus.  Check the AC bus input data.\n",
              ACptr->Num, ACptr->Name);
          InputError = true;
        } else {
          ptrp = new AClist;
          ptrp->AC = ACptr;
          ptrp->Prev = nullptr;
          if (ACptr->Cont->ContBus == nullptr) {
            ACptr->Cont->ContBus = ptrp;
            ACptr->Cont->ContBus->Next = ACptr->Cont->ContBus->Prev = nullptr;
          } else {
            ptrs = ACptr->Cont->ContBus;
            ACptr->Cont->ContBus = ptrp;
            ptrs->Prev = ptrp;
            ptrp->Next = ptrs;
          }
        }
        if (!QRcont) {
          ACptr->Cont->Cont = ACptr->Cont;
          ACptr->Cont = nullptr;
        } else {
          if (flag2Vcontrol) {
            if (ACptr->Cont->Kbg1 > 0)
              ACptr->Kbg1 = ACptr->Qmax / ACptr->Cont->Kbg1;
            else
              ACptr->Kbg1 = 1;
            if (ACptr->Cont->Kbg2 < 0)
              ACptr->Kbg2 = ACptr->Qmin / ACptr->Cont->Kbg2;
            else
              ACptr->Kbg2 = 1;
            ACptr->Kbg = ACptr->Kbg1;
          } else {
            ACptr->Kbg1 = ACptr->Kbg;
            ACptr->Kbg2 = ACptr->Kbg;
          }
          ACptr->Cont->Qr = ACptr->Cont->Qr + ACptr->Qg;
        }
      } else {
        fprintf(stderr, "ERROR: The remote controlled bus of PV bus %d %s\n",
                ACptr->Num, ACptr->Name);
        fprintf(stderr,
                "       has not been defined.  Check the AC bus input data.\n");
        InputError = true;
      }
    }
    if (strpbrk(ACptr->Type, "C") && ACptr->Kbg < 1) {
      fprintf(stderr, "ERROR: The voltage controlled bus %d %s does not have\n",
              ACptr->Num, ACptr->Name);
      fprintf(stderr, "       any generator controlling the voltage. Check the "
                      "AC bus input data.\n");
      InputError = true;
    }
    if (Rcont && strpbrk(ACptr->Type, "T"))
      ACptr->Cont = nullptr;
    if (Narea > 1) {
      if (ACptr->Area == nullptr)
        for (Aptr = dataPtr->Area; Aptr != nullptr; Aptr = Aptr->Next) {
          for (i = 1; i <= 10; i++) {
            if (!strcmp(ACptr->Zone, Aptr->Zone[i])) {
              ACptr->Area = Aptr;
              if (strpbrk(ACptr->Type, "S"))
                Aptr->i++;
              if (Aptr->i > 1) {
                fprintf(stderr, "ERROR: Area %d %s has 2 slack buses.\n",
                        Aptr->N, Aptr->Name);
                fprintf(stderr, "       Check AC area and bus input data.\n");
                InputError = true;
              }
              break;
            }
          }
          if (ACptr->Area != nullptr)
            break;
        }
      if (ACptr->Area == nullptr) {
        fprintf(stderr, "ERROR: AC/DC bus %d %s is not in any area.\n",
                ACptr->Num, ACptr->Name);
        fprintf(stderr, "       Check area and bus input data.\n");
        InputError = true;
      } else {
        if (!strcmp(ACptr->Area->Name, "")) {
          fprintf(stderr, "ERROR: Area %d has not been defined.\n",
                  ACptr->Area->N);
          fprintf(stderr, "       Check area and bus input data.\n");
          InputError = true;
        }
        if (ACptr == ACptr->Area->Slack && !strpbrk(ACptr->Type, "S")) {
          strcat(ACptr->Type, "A");
          ACptr->Area->i++;
          if (ACptr->Area->i > 1) {
            fprintf(stderr, "ERROR: Area %d %s has 2 slack buses.\n",
                    ACptr->Area->N, ACptr->Area->Name);
            fprintf(stderr, "       Check AC area and bus input data.\n");
            InputError = true;
          }
        }
        Aptr = ACptr->Area;
        ptrp = Aptr->AC;
        Aptr->AC = new AClist;
        Aptr->AC->AC = ACptr;
        Aptr->AC->Next = ptrp;
        Aptr->AC->Prev = nullptr;
        if (ptrp != nullptr)
          ptrp->Prev = Aptr->AC;
        if (ExistParameter('6') && Aptr->BSptr == nullptr)
          Aptr->Slack = Aptr->BSptr = ACptr;
      }
    } else
      ACptr->Area = nullptr;
    if (ACptr->Ang >= 1000)
      ACptr->Ang = AngSlack;
    if (ACptr->DPg == 0 && strpbrk(ACptr->Type, "AS")) {
      if (ExistParameter('6')) {
        if (strpbrk(ACptr->Type, "S"))
          ACptr->DPg = 1;
      } else
        ACptr->DPg = 1;
    }
    if (ACptr->PgMax <= 0)
      ACptr->PgMax = 99999999.;
    if (ACptr->Smax <= 0)
      ACptr->Smax = 99999999.;
    ACptr->DPG = ACptr->DPg;
    if (strpbrk(ACptr->Type, "X")) {
      if (ACptr->Cont == nullptr) {
        fprintf(stderr,
                "ERROR: The reactance controlled bus %d %s does not have\n",
                ACptr->Num, ACptr->Name);
        fprintf(stderr, "       a controlling bus. Check the AC input data.\n");
        InputError = true;
      } else {
        ACptrp = ACptr->Cont;
        ACptr->Vmax = ACptrp->Vmax;
        ACptr->Vmin = ACptrp->Vmin;
        if (ACptr->Vmin >= ACptr->Vmax) {
          fprintf(
              stderr,
              "ERROR: The reactance controlling bus %d %s has inconsistent\n",
              ACptrp->Num, ACptrp->Name);
          fprintf(stderr, "       voltage limits. Check the AC input data.\n");
          InputError = true;
        } else if (ACptrp->steps == 0) {
          fprintf(stderr,
                  "ERROR: The reactance controlling bus %d %s has zero MVAr "
                  "steps.\n",
                  ACptrp->Num, ACptrp->Name);
          fprintf(stderr, "       Check the AC input data.\n");
          InputError = true;
        }
      }
    }

    /* -------------------------- AC elements --------------------------- */
    i = 0;
    for (ELptr = ACptr->Reg; ELptr != nullptr; ELptr = ELptr->Next) {
      Eptr = ELptr->Eptr;
      if ((!strcmp(Eptr->Type, "R") && !strpbrk(ACptr->Type, "T")) ||
          (!strcmp(Eptr->Type, "RV") && !strpbrk(ACptr->Type, "R"))) {
        fprintf(stderr, "ERROR: LTC volt. controlled bus %d %s is not PQ.\n",
                ACptr->Num, ACptr->Name);
        fprintf(stderr, "       Check the AC/DC bus input cards.\n");
        InputError = true;
      }
      if (!strcmp(Eptr->Type, "R") || !strcmp(Eptr->Type, "RV"))
        i++;
    }
    if ((ACptr->Reg == nullptr || (ACptr->Reg != nullptr && i == 0)) &&
        strpbrk(ACptr->Type, "T")) {
      fprintf(stderr,
              "ERROR: LTC volt. controlled bus %d %s does not have LTCs.\n",
              ACptr->Num, ACptr->Name);
      fprintf(stderr, "       Check the AC/DC bus and element input cards.\n");
      InputError = true;
    }
    if (strpbrk(ACptr->Type, "S,A")) {
      ptrp = dataPtr->KGbus;
      ptrs = new AClist;
      ptrs->AC = ACptr;
      ptrs->Next = ptrp;
      ptrs->Prev = nullptr;
      dataPtr->KGbus = ptrs;
      if (ptrp != nullptr)
        ptrp->Prev = ptrs;
    }
    ACptr->PG = ACptr->Pg;
    ACptr->PL = ACptr->Pl;
    ACptr->QL = ACptr->Ql;
    ACptrp = ACptr->Next;
    if (ieee) {
      ACptr->Next = ACptr->Prev;
      ACptr->Prev = ACptrp;
      dataPtr->ACbus = ACptr;
    }
    ACptr = ACptrp;
  }
  Eptr = dataPtr->Element;
  while (Eptr != nullptr) {
    if (strpbrk(Eptr->Type, "R") && !strcmp(Eptr->Zone, "")) {
      fprintf(stderr, "ERROR: Reg. transf. from %d %s to %d %s\n",
              Eptr->From->Num, Eptr->From->Name, Eptr->To->Num, Eptr->To->Name);
      fprintf(stderr, "       has not been completely defined.  Check T cards "
                      "on WSCC format.\n");
      InputError = true;
    }
    if (!strcmp(Eptr->Type, "R") || !strcmp(Eptr->Type, "RV"))
      strcpy(Eptr->Ctype, "V");
    else if (strpbrk(Eptr->Type, "PM"))
      strcpy(Eptr->Ctype, "P");
    else if (strpbrk(Eptr->Type, "QN"))
      strcpy(Eptr->Ctype, "Q");
    if (strcmp(Eptr->Owner, "")) {
      if (!strcmp(Eptr->Owner, Eptr->From->Owner))
        strcpy(Eptr->Zone, Eptr->From->Zone);
      else
        strcpy(Eptr->Zone, Eptr->To->Zone);
    }
    if (Narea > 1) {
      Aptr = Eptr->From->Area;
      Aptrp = Eptr->To->Area;
      if (Aptr != Aptrp && Aptr != nullptr && Aptrp != nullptr) {
        Aptr->Elem = (ElementList *)AddElemToList(Aptr->Elem, Eptr);
        Aptrp->Elem = (ElementList *)AddElemToList(Aptrp->Elem, Eptr);
      }
      if (Eptr->Meter == nullptr) {
        if (strcmp(Eptr->Owner, "")) {
          if (strcmp(Eptr->Owner, Eptr->From->Owner))
            Eptr->Meter = Eptr->From;
          else
            Eptr->Meter = Eptr->To;
        } else if (Eptr->Area == Aptr)
          Eptr->Meter = Eptr->To;
        else
          Eptr->Meter = Eptr->From;
      }
    }
    Eptrp = Eptr->Next;
    Eptr->Next = Eptr->Prev;
    Eptr->Prev = Eptrp;
    dataPtr->Element = Eptr;
    Eptr = Eptrp;
  }

  /* ------------------- DC buses ------------------------ */
  for (DCptr = dataPtr->DCbus; DCptr != nullptr; DCptr = DCptr->Next) {
    ACptr = DCptr->AC;
    if (ACptr == nullptr || !strcmp(DCptr->Type, "")) {
      fprintf(
          stderr,
          "ERROR: DC bus %8s has not been fully defined in the input data.\n",
          DCptr->Name);
      fprintf(stderr, "       Check the BD, LD, and/or BZ cards.\n");
      InputError = true;
    }

    /* --------- Fix data read in WSCC/BPA Format --------- */
    if (DCptr->Xc == 0) {
      for (ELptrp = ELptr = ACptr->Elem; ELptr != nullptr;
           ELptrp = ELptr, ELptr = ELptr->Next) {
        Eptr = ELptr->Eptr;
        if (!strcmp(Eptr->From->Name, DCptr->Name)) {
          ACptrp = Eptr->From;
          break;
        }
        if (!strcmp(Eptr->To->Name, DCptr->Name)) {
          ACptrp = Eptr->To;
          break;
        }
      }
      DCptr->Name[8] = '\0';
      if (ELptr == nullptr) {
        fprintf(stderr,
                "ERROR: The tranformer for DC bus %8s has not been defined in "
                "the\n",
                DCptr->Name);
        fprintf(stderr, "       input data. Check the T data cards.\n");
        InputError = true;
      } else {
        DCptr->Xc = -Eptr->B / (Eptr->G * Eptr->G + Eptr->B * Eptr->B);
        DCptr->Xc = (DCptr->Xc * ACptr->KV * ACptr->KV / Sn) * DCptr->Nbr;
        DCptr->Tap = Eptr->Tap;
        if (strpbrk(Eptr->Type, "P,Q,M,N")) {
          fprintf(stderr,
                  "ERROR: DC bus %8s regulating transfomer may only control "
                  "voltage.\n",
                  DCptr->Name);
          fprintf(stderr,
                  "       Check the related R transfomer input cards.\n");
          InputError = true;
        }
        if (!strcmp(Eptr->Type, "R")) {
          if (Eptr->Cont != ACptrp) {
            fprintf(stderr,
                    "ERROR: DC bus %8s regulating transfomer must control "
                    "voltage\n",
                    DCptr->Name);
            fprintf(stderr, "       of DC bus.  Check the related R transfomer "
                            "input cards.\n");
            InputError = true;
          }
          NregV--;
          DCptr->TapMax = Eptr->Tmax;
          DCptr->TapMin = Eptr->Tmin;
        } else {
          DCptr->TapMax = 1.1;
          DCptr->TapMin = 0.9;
        }
      }
      /* Remove extra element from AC data base */
      if (ELptrp == ELptr) {
        ACptr->Elem = ELptr->Next;
        if (ACptr->Elem == nullptr) {
          fprintf(
              stderr,
              "ERROR: AC bus %d %s is isolated from the rest of the system.\n",
              ACptr->Num, ACptr->Name);
          fprintf(stderr, "       Check the related AC element input cards.\n");
          InputError = true;
        }
      } else
        ELptrp->Next = ELptr->Next;
      delete ELptr;
      Aptr = ACptr->Area;
      if (Aptr != nullptr) {
        for (ELptrp = ELptr = Aptr->Elem; ELptr != nullptr;
             ELptrp = ELptr, ELptr = ELptr->Next)
          if (Eptr == ELptr->Eptr)
            break;
        if (ELptr != nullptr) {
          if (ELptrp == ELptr)
            Aptr->Elem = ELptr->Next;
          else
            ELptrp->Next = ELptr->Next;
          delete ELptr;
        }
      }
      NacEl--;
      if (Eptr->Prev == nullptr) {
        Eptrp = Eptr->Next;
        dataPtr->Element = Eptrp;
        Eptrp->Prev = nullptr;
      } else {
        Eptrp = Eptr->Prev;
        Eptrp->Next = Eptr->Next;
        if (Eptr->Next != nullptr)
          Eptr->Next->Prev = Eptrp;
      }
      delete Eptr;
      /* Remove DC bus from AC data base */
      Aptr = ACptr->Area;
      if (Aptr != nullptr) {
        for (ptrp = Aptr->AC; ptrp != nullptr; ptrp = ptrp->Next)
          if (ACptrp == ptrp->AC)
            break;
        if (ptrp != nullptr) {
          if (ptrp->Prev == nullptr) {
            ptrs = ptrp->Next;
            Aptr->AC = ptrs;
            ptrs->Prev = nullptr;
          } else {
            ptrs = ptrp->Prev;
            ptrs->Next = ptrp->Next;
            if (ptrp->Next != nullptr)
              ptrp->Next->Prev = ptrs;
          }
          delete ptrp;
        }
      }
      Nac--;
      for (ACptrs = ACptrp->Next; ACptrs != nullptr;
           ACptrs->N--, ACptrs = ACptrs->Next)
        ;
      if (ACptrp->Prev == nullptr) {
        ACptrs = ACptrp->Next;
        dataPtr->ACbus = ACptrs;
        ACptrs->Prev = nullptr;
      } else {
        ACptrs = ACptrp->Prev;
        ACptrs->Next = ACptrp->Next;
        if (ACptrp->Next != nullptr)
          ACptrp->Next->Prev = ACptrs;
      }
      delete ACptrp;
    }

    if (DCptr->MVA > 0) {
      KVs = ACptr->KV * DCptr->Ntrf / DCptr->Nbr;
      DCptr->Xc = DCptr->Xc / DCptr->Nbr;
      DCptr->Xc = (DCptr->Xc * KVs * KVs / DCptr->MVA) * DCptr->Nbr;
    }
    if (DCptr->Xc <= 0.1) {
      fprintf(
          stderr,
          "ERROR: DC bus %8s has a negative or zero commutation reactance ->\n",
          DCptr->Name);
      fprintf(stderr, "       Xc=%6.3lf (Ohms)\n", DCptr->Xc);
      fprintf(stderr, "       Check the BD and/or BZ cards.\n");
      InputError = true;
    }
    if (DCptr->To == nullptr) {
      fprintf(stderr, "ERROR: DC bus %8s is isolated.\n", DCptr->Name);
      fprintf(stderr, "       Check the DC line input cards.\n");
      InputError = true;
    }
    DCptr->Area = ACptr->Area;
    if (!strcmp(DCptr->Cont1, "AL") || !strcmp(DCptr->Cont2, "AL")) {
      if (DCptr->Alfa > DCptr->AlfaMax || DCptr->Alfa < DCptr->AlfaMin) {
        fprintf(stderr, "ERROR: DC bus %8s has ALPHA outside its limits.\n",
                DCptr->Name);
        fprintf(stderr, "       Check the BD (limits) and/or BZ card.\n");
        InputError = true;
      }
      DCptr->Gamma = 180 - DCptr->Alfa;
    } else if (!strcmp(DCptr->Cont1, "GA") || !strcmp(DCptr->Cont2, "GA")) {
      if (DCptr->Gamma < DCptr->GammaMin) {
        fprintf(stderr, "ERROR: DC bus %8s has GAMMA outside its limits.\n",
                DCptr->Name);
        fprintf(stderr, "       Check the BD (limits) and/or BZ card.\n");
        InputError = true;
      }
      DCptr->Alfa = 180 - DCptr->Gamma;
    } else if (!strcmp(DCptr->Type, "R")) {
      DCptr->Alfa = DCptr->AlfaMin;
      DCptr->Gamma = 180 - DCptr->Alfa;
    } else {
      DCptr->Gamma = DCptr->GammaMin;
      DCptr->Alfa = 180 - DCptr->Gamma;
    }
    if (DCptr->Alfa > DCptr->AlfaMax || DCptr->Alfa < DCptr->AlfaMin ||
        DCptr->Gamma < DCptr->GammaMin)
      fprintf(
          stderr,
          "***Warning: DC bus %8s could have wrong ALPHA or GAMMA limits.\n",
          DCptr->Name);
    if (DCptr->Tap <= 0)
      DCptr->Tap = 1;
    if (DCptr->Tap < DCptr->TapMin)
      DCptr->Tap = DCptr->TapMin;
    if (DCptr->Tap > DCptr->TapMax)
      DCptr->Tap = DCptr->TapMax;
  }

  /* ------------------- DC elements ------------------------ */
  if (Ndc != (2 * NdcEl)) {
    fprintf(stderr, "ERROR: There are inconsistencies between the DC bus and "
                    "DC line input data.\n");
    fprintf(stderr, "       Check DC input data and remember that the program "
                    "just allows for \n");
    fprintf(stderr, "       two-terminal HVDC links.\n");
    InputError = true;
  }
  for (DCptr = dataPtr->DCbus; DCptr != nullptr; DCptr = DCptr->Next) {
    DCptrp = DCptr->To;
    if (DCptr->N != 0 && DCptrp->N != 0) {
      if ((strcmp(DCptr->Type, "R") || strcmp(DCptrp->Type, "I")) &&
          (strcmp(DCptr->Type, "I") || strcmp(DCptrp->Type, "R"))) {
        fprintf(stderr,
                "ERROR: Both converters for the DC link between %8s and %8s\n",
                DCptr->Name, DCptrp->Name);
        fprintf(stderr, "       are either rectifiers or inverters.\n");
        InputError = true;
      }
      if ((!strcmp(DCptr->Cont1, "ID") || !strcmp(DCptr->Cont2, "ID")) &&
          (!strcmp(DCptrp->Cont1, "ID") || !strcmp(DCptrp->Cont2, "ID"))) {
        fprintf(stderr,
                "ERROR: Both converters for the DC link between %8s and %8s\n",
                DCptr->Name, DCptrp->Name);
        fprintf(stderr, "       are controlling the current.\n");
        InputError = true;
      }
      if (DCptr->Area == DCptrp->Area || !strcmp(DCptr->Zone, DCptrp->Zone))
        DCptr->Meter = DCptrp->Meter = nullptr;
      if (DCptr->Area != DCptrp->Area)
        for (i = 1; i <= 2; i++) {
          if (i == 1)
            Aptr = DCptr->Area;
          else
            Aptr = DCptrp->Area;
          if (Aptr != nullptr) {
            ptr = Aptr->DC;
            Aptr->DC = new DClist;
            if (i == 1)
              Aptr->DC->DC = DCptr;
            else
              Aptr->DC->DC = DCptrp;
            Aptr->DC->Next = ptr;
          }
        }
      if (DCptr->Id > 0)
        DCptrp->Id = DCptr->Id;
      else if (DCptrp->Id > 0)
        DCptr->Id = DCptrp->Id;
      DCptr->N = 0;
      DCptrp->N = 0;
    }
  }

  /* -------------------------- Areas --------------------------- */
  for (Aptr = dataPtr->Area; Aptr != nullptr; Aptr = Aptr->Next) {
    flag = true;
    i = 0;
    j = 0;
    for (ELptr = Aptr->Elem; ELptr != nullptr; ELptr = ELptr->Next) {
      Eptr = ELptr->Eptr;
      flag = false;
      i++;
      if (!strcmp(Eptr->Type, "RP"))
        j++;
      if (Eptr->From->Area != Aptr)
        Aptrp = Eptr->From->Area;
      else
        Aptrp = Eptr->To->Area;
      if (!Acont && strpbrk(Aptr->Slack->Type, "S") &&
          !strpbrk(Aptrp->Slack->Type, "S"))
        ExpandSlack(Aptr->Slack, Aptrp);
    }
    if (!flag && i == j) {
      fprintf(stderr,
              "ERROR: All tie lines for area %d %s are P reg. transf.\n",
              Aptr->N, Aptr->Name);
      fprintf(stderr,
              "       Change at least one reg. transf. to a standard one.\n");
      InputError = true;
    }
    if (flag && !strpbrk(Aptr->Slack->Type, "S")) {
      strcat(Aptr->Slack->Type, "S");
      Nslack++;
    }
    ptrp = Aptr->AC;
    while (ptrp != nullptr) {
      ptrs = ptrp->Next;
      ptrp->Next = ptrp->Prev;
      ptrp->Prev = ptrs;
      if (ptrs == nullptr)
        Aptr->AC = ptrp;
      ptrp = ptrs;
    }
  }

  /* -------------------------- FACTS --------------------------- */
  for (SVCptr = dataPtr->SVCbus; SVCptr != nullptr; SVCptr = SVCptr->Next) {
    ACptr = SVCptr->Ctrl;
    if (ACptr->Cont == nullptr) {
      fprintf(stderr,
              "ERROR: The SVC controlled bus %d %s is already controlled.\n",
              ACptr->N, ACptr->Name);
      fprintf(stderr, "       Check the AC bus cards.\n");
      InputError = true;
    }
    ACptr = SVCptr->From;
    if (ACptr->Cont == nullptr) {
      fprintf(stderr, "ERROR: The SVC bus %d %s is a voltage controlled bus.\n",
              ACptr->N, ACptr->Name);
      fprintf(stderr, "       Check the AC bus cards.\n");
      InputError = true;
    }
    for (ELptr = ACptr->Elem, i = 0; ELptr != nullptr; ELptr = ELptr->Next, i++)
      ;
    if (i > 1) {
      fprintf(stderr,
              "ERROR: The SVC bus %d %s has more than one AC element connected "
              "to it.\n",
              ACptr->N, ACptr->Name);
      fprintf(stderr, "       Check the AC element cards.\n");
      InputError = true;
    }
  }
  for (STATCOMptr = dataPtr->STATCOMbus; STATCOMptr != nullptr;
       STATCOMptr = STATCOMptr->Next) {
    ACptr = STATCOMptr->Ctrl;
    if (ACptr->Cont == nullptr) {
      fprintf(
          stderr,
          "ERROR: The STATCOM controlled bus %d %s is already controlled.\n",
          ACptr->N, ACptr->Name);
      fprintf(stderr, "       Check the AC bus cards.\n");
      InputError = true;
    }
    ACptr = STATCOMptr->From;
    if (ACptr->Cont == nullptr) {
      fprintf(stderr,
              "ERROR: The STATCOM bus %d %s is a voltage controlled bus.\n",
              ACptr->N, ACptr->Name);
      fprintf(stderr, "       Check the AC bus cards.\n");
      InputError = true;
    }
  }
}

/* --------- ExpandSlack --------- */
void ExpandSlack(ACbusData *BSptr, AreaData *Aptr) {
  AreaData *Aptrp;
  ElementList *ELptr;
  ElementData *Eptr;

  Aptr->Slack = BSptr;
  for (ELptr = Aptr->Elem; ELptr != nullptr; ELptr = ELptr->Next) {
    Eptr = ELptr->Eptr;
    if (Eptr->From->Area != Aptr)
      Aptrp = Eptr->From->Area;
    else
      Aptrp = Eptr->To->Area;
    if (!strpbrk(Aptrp->Slack->Type, "S"))
      ExpandSlack(Aptr->Slack, Aptrp);
  }
}

/* --------- WriteSummary --------- */
void WriteSummary(void) {
  ACbusData *ACptr;
  int i;

  fprintf(stderr, "Summary of input data for case:\n");
  i = 0;
  while (i <= 2 && dataPtr->Title[0][0] != '\0') {
    fprintf(stderr, "%s", dataPtr->Title[i]);
    i++;
  }
  fprintf(stderr, "            AC buses -> %d\n", Nac);
  fprintf(stderr, "            PV buses -> %d\n", Nvolt);
  fprintf(stderr, "            X buses  -> %d\n", NXvolt);
  fprintf(stderr, "            Z buses  -> %d\n", NZvolt);
  fprintf(stderr, "            AC elem. -> %d\n", NacEl);
  fprintf(stderr, "         V Reg. Trf. -> %d\n", NregV);
  fprintf(stderr, "        PQ Reg. Trf. -> %d\n", NregPQ);
  fprintf(stderr, "            DC buses -> %d\n", Ndc);
  fprintf(stderr, "            DC lines -> %d\n", NdcEl);
  fprintf(stderr, "                SVCs -> %d\n", Nsvc);     /* FACTS */
  fprintf(stderr, "               TCSCs -> %d\n", Ntcsc);    /* FACTS */
  fprintf(stderr, "            STATCOMs -> %d\n", Nstatcom); /* FACTS */
  fprintf(stderr, "           No. Areas -> %d\n", Narea);
  fprintf(stderr, "   Reference Bus(es) -> ");
  i = 0;
  for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next)
    if (strpbrk(ACptr->Type, "S")) {
      if (i > 0)
        fprintf(stderr, "                        ");
      fprintf(stderr, "%d %s (Angle=%6.2lf deg.)\n", ACptr->Num, ACptr->Name,
              ACptr->Ang * 180 / PI);
      i++;
    }
  if (i == 0)
    fprintf(stderr, "\n");
  return;
}

/* ------------------------- ReadData ------------------------------- */
void ReadData(char *Name)
/* Main routine. */
{
  InputDataFile = (FILE *)OpenInput(Name);
  dataPtr = new Data;
  dataPtr->Title[0][0] = '\0';
  dataPtr->ACbus = nullptr;
  dataPtr->DCbus = nullptr;
  dataPtr->Element = nullptr;
  dataPtr->Area = nullptr;
  dataPtr->KGbus = nullptr;
  dataPtr->SVCbus = nullptr;     /* FACTS */
  dataPtr->TCSCbus = nullptr;    /* FACTS */
  dataPtr->STATCOMbus = nullptr; /* FACTS */
  Nac = 0;
  Ndc = 0;
  NacEl = 0;
  NdcEl = 0;
  Nsvc = Nstatcom = 0;
  Ntcsc = NtcscVar = 0; /* FACTS */
  LineNum = 0;
  Nvolt = 0;
  Nslack = 0;
  Narea = 0;
  NregPQ = NregV = 0;
  NZvolt = NXvolt = 0;
  InputError = false;
  flag2Vcontrol = ExistParameter('#');
  if (ExistParameter('I'))
    ReadIEEE();
  else if (ExistParameter('6'))
    ReadITALY();
  else
    ReadWSCC();
  ErrorDetect();
  WriteSummary();
  if (InputError == true) {
    fprintf(stderr,
            "*** The data has errors! Please review the input file. ***\n");
    exit(ERROREXIT);
  } else
    fprintf(stderr, "*** The data has been read successfully ***\n");
  delete[] Name;
}
