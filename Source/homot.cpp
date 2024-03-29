/* Homotopy Continuation Method: Main. */

#include "homot.h"

void AddLoad(void);
bool CheckVIlimits(bool flagVoltage, bool flagCurrent);

/* ------- Global Variables ------ */
VALUETYPE *Dx, Dparam, param0, *x0, *x0p, Kh, Htol, AngTr, VoltageSum, VoltageSum0, DxiMax, VSFone, VSFinf, SF, TVI, ZeroDx, lambda_o, TotalPl, TotalQl,
    TotalPg, TotalQg;
INDEX RedNumEq, NewNumEq, CountSteps, NumSteps, TVIbus, TVIrank, VSFbus;
bool flagReduceSystem, *DxZero, flagPrintTotalPl, flagPrintTotalQl, flagPrintTotalPg, flagPrintTotalQg;
AClist *Vlist = nullptr, *Vlistp = nullptr;
FILE *OutputHomot;
int field, SD0;
extern bool flagBS, flagVloads;

/* --------------------------- DCsetup --------------------------------- */
bool DCsetup(void)
/* Setup inital DC control mode, i.e., rectifier -> alfa,Id (or P),
                                       inverter  -> gamma,Vd. */
{
  DCbusData *DCptrR, *DCptrI;
  bool flag = false, flagP = false, flagT = false;

  for (DCptrR = dataPtr->DCbus; DCptrR != nullptr; DCptrR = DCptrR->Next)
    if (!strcmp(DCptrR->Type, "R")) {
      DCptrI = DCptrR->To;
      DCptrI->VdN = DCptrI->Vd;
      DCptrR->AlfaN = DCptrR->Alfa;
      if (strcmp(DCptrR->Cont1, "AL") && strcmp(DCptrR->Cont2, "AL"))
        flag = true;
      else if (strcmp(DCptrR->Cont1, "AT") && strcmp(DCptrR->Cont2, "AT"))
        flag = true;
      else if (!strpbrk(DCptrR->Cont1, "IP") && !strpbrk(DCptrR->Cont2, "IP"))
        flag = true;
      else if (strcmp(DCptrI->Cont1, "GA") && strcmp(DCptrI->Cont2, "GA"))
        flag = true;
      else if (strcmp(DCptrI->Cont1, "AT") && strcmp(DCptrI->Cont2, "AT"))
        flag = true;
      else if (strcmp(DCptrI->Cont1, "VD") && strcmp(DCptrI->Cont2, "VD"))
        flag = true;
      if (!strcmp(DCptrR->Cont1, "PA") || !strcmp(DCptrR->Cont2, "PA") || !strcmp(DCptrI->Cont1, "PA") || !strcmp(DCptrI->Cont2, "PA"))
        flagP = true;
      if (!strcmp(DCptrR->Cont1, "AL") || !strcmp(DCptrR->Cont2, "AL"))
        strcpy(DCptrR->Cont1, "AL");
      else if (!strcmp(DCptrR->Cont1, "AT") || !strcmp(DCptrR->Cont2, "AT")) {
        strcpy(DCptrR->Cont1, "AT");
        flagT = true;
      } else
        strcpy(DCptrR->Cont1, "AL");
      if (flagP)
        strcpy(DCptrR->Cont2, "PA");
      else
        strcpy(DCptrR->Cont2, "ID");
      strcpy(DCptrI->Cont1, "GA");
      if (flagT)
        strcpy(DCptrI->Cont2, "AT");
      else
        strcpy(DCptrI->Cont2, "VD");
      if (DCptrR->Tap > DCptrR->TapMax) {
        flag = true;
        DCptrR->TapMax = DCptrR->Tap;
        fprintf(stderr,
                "***Warning: The program will change the maximum tap limit for "
                "rect. %s\n",
                DCptrR->Name);
      }
      if (DCptrR->Tap < DCptrR->TapMin) {
        flag = true;
        DCptrR->TapMin = DCptrR->Tap;
        fprintf(stderr,
                "***Warning: The program will change the minimum tap limit for "
                "rect. %s\n",
                DCptrR->Name);
      }
      if (DCptrI->Tap > DCptrI->TapMax) {
        flag = true;
        DCptrI->TapMax = DCptrI->Tap;
        fprintf(stderr,
                "***Warning: The program will change the maximum tap limit for "
                "inv. %s\n",
                DCptrI->Name);
      }
      if (DCptrI->Tap < DCptrI->TapMin) {
        flag = true;
        DCptrI->TapMin = DCptrI->Tap;
        fprintf(stderr,
                "***Warning: The program will change the minimum tap limit for "
                "inv. %s\n",
                DCptrI->Name);
      }
      if (DCptrI->Gamma != DCptrI->GammaMin) {
        DCptrI->GammaMin = DCptrI->Gamma;
        fprintf(stderr,
                "***Warning: The program will change the minimum gamma limit "
                "for inv. %s\n",
                DCptrI->Name);
      }
      fprintf(stderr, "***Warning: The startup control mode for HVDC link from %s to %s\n", DCptrR->Name, DCptrI->Name);
      fprintf(stderr, "            is: rectifier -> %s-%s\n", DCptrR->Cont1, DCptrR->Cont2);
      fprintf(stderr, "                inverter  -> %s-%s\n", DCptrI->Cont1, DCptrI->Cont2);
    }
  return (flag);
}

/* -------------------- ChangeDCmode ----------------------- */
bool ChangeDCmode(void)
/* Change DC control mode when DC variable at a limit. */
{
  DCbusData *DCptrR, *DCptrI;
  bool flag = false, flagp = false;
  INDEX i;

  i = NacVar;
  for (DCptrR = dataPtr->DCbus; DCptrR != nullptr; DCptrR = DCptrR->Next)
    if (!strcmp(DCptrR->Type, "R")) {
      DCptrI = DCptrR->To;
      if (!strcmp(DCptrR->Cont1, "AL") && strpbrk(DCptrR->Cont2, "IP") && (DCptrR->Tap >= DCptrR->TapMax || DCptrR->Tap <= DCptrR->TapMin)) {
        flag = flagp = true;
        strcpy(DCptrR->Cont1, "AT");
      } else if (!strcmp(DCptrR->Cont1, "AT") && strpbrk(DCptrR->Cont2, "IP") && DCptrR->Alfa <= DCptrR->AlfaMin) {
        flag = flagp = true;
        if (!strcmp(DCptrR->Cont2, "ID"))
          strcpy(DCptrI->Cont2, "ID");
        else
          strcpy(DCptrI->Cont2, "PA");
        strcpy(DCptrR->Cont1, "AL");
        strcpy(DCptrR->Cont2, "AT");
        strcpy(DCptrI->Cont1, "AT");
      } else if (!strcmp(DCptrR->Cont1, "AT") && strpbrk(DCptrR->Cont2, "IP") && DCptrR->Alfa == DCptrR->AlfaN) {
        flag = flagp = true;
        strcpy(DCptrR->Cont1, "AL");
      } else if (!strcmp(DCptrI->Cont1, "AT") && strpbrk(DCptrI->Cont2, "IP") && DCptrI->Tap <= DCptrI->TapMin && DCptrI->Gamma <= DCptrI->GammaMin) {
        flag = flagp = true;
        if (!strcmp(DCptrI->Cont2, "ID"))
          strcpy(DCptrR->Cont2, "ID");
        else
          strcpy(DCptrR->Cont2, "PA");
        strcpy(DCptrR->Cont1, "AT");
        strcpy(DCptrI->Cont1, "GA");
        strcpy(DCptrI->Cont2, "AT");
      } else if (!strcmp(DCptrI->Cont1, "AT") && strpbrk(DCptrI->Cont2, "IP") && DCptrI->Tap > DCptrI->TapMin && DCptrI->Gamma <= DCptrI->GammaMin) {
        flag = flagp = true;
        if (!strcmp(DCptrI->Cont2, "ID"))
          strcpy(DCptrR->Cont2, "ID");
        else
          strcpy(DCptrR->Cont2, "PA");
        strcpy(DCptrR->Cont1, "AT");
        strcpy(DCptrI->Cont1, "GA");
        strcpy(DCptrI->Cont2, "VD");
      } else if (!strcmp(DCptrI->Cont1, "GA") && !strcmp(DCptrI->Cont2, "VD") && (DCptrI->Tap <= DCptrI->TapMin || DCptrI->Tap >= DCptrI->TapMax)) {
        flag = flagp = true;
        strcpy(DCptrI->Cont2, "AT");
      } else if (!strcmp(DCptrI->Cont1, "GA") && !strcmp(DCptrI->Cont2, "AT") && DCptrI->Vd == DCptrI->VdN) {
        flag = flagp = true;
        strcpy(DCptrI->Cont2, "VD");
      }
      if (flagH && flag) {
        if (strcmp(DCptrR->Cont1, "VD") && strcmp(DCptrR->Cont2, "VD")) {
          i++;
          x0[i] = DCptrR->Vd;
        }
        if (strcmp(DCptrR->Cont1, "AT") && strcmp(DCptrR->Cont2, "AT")) {
          i++;
          x0[i] = DCptrR->Tap * DCptrR->Ntrf;
        }
        if (strcmp(DCptrR->Cont1, "AL") && strcmp(DCptrR->Cont2, "AL")) {
          i++;
          x0[i] = cos(DCptrR->Alfa);
        }
        if (strcmp(DCptrR->Cont1, "GA") && strcmp(DCptrR->Cont2, "GA")) {
          i++;
          x0[i] = cos(DCptrR->Gamma);
        }
        i++;
        x0[i] = DCptrR->MVA;
        if (strcmp(DCptrR->Cont1, "PA") && strcmp(DCptrR->Cont2, "PA")) {
          i++;
          x0[i] = DCptrR->P;
        }
        if (strcmp(DCptrR->Cont1, "QA") && strcmp(DCptrR->Cont2, "QA")) {
          i++;
          x0[i] = DCptrR->Q;
        }
        if (strcmp(DCptrI->Cont1, "VD") && strcmp(DCptrI->Cont2, "VD")) {
          i++;
          x0[i] = DCptrI->Vd;
        }
        if (strcmp(DCptrI->Cont1, "AT") && strcmp(DCptrI->Cont2, "AT")) {
          i++;
          x0[i] = DCptrI->Tap * DCptrI->Ntrf;
        }
        if (strcmp(DCptrI->Cont1, "AL") && strcmp(DCptrI->Cont2, "AL")) {
          i++;
          x0[i] = cos(DCptrI->Alfa);
        }
        if (strcmp(DCptrI->Cont1, "GA") && strcmp(DCptrI->Cont2, "GA")) {
          i++;
          x0[i] = cos(DCptrI->Gamma);
        }
        i++;
        x0[i] = DCptrI->MVA;
        if (strcmp(DCptrI->Cont1, "PA") && strcmp(DCptrI->Cont2, "PA")) {
          i++;
          x0[i] = DCptrI->P;
        }
        if (strcmp(DCptrI->Cont1, "QA") && strcmp(DCptrI->Cont2, "QA")) {
          i++;
          x0[i] = DCptrI->Q;
        }
        i++;
        x0[i] = DCptrI->MVA;
        if (strcmp(DCptrR->Cont1, "ID") && strcmp(DCptrR->Cont2, "ID") && strcmp(DCptrI->Cont1, "ID") && strcmp(DCptrI->Cont2, "ID")) {
          i++;
          x0[i] = DCptrR->Id;
        }
      }
      if (flagp) {
        flagp = false;
        fprintf(stderr, "***Warning: HVDC link from %s to %s will switch control mode\n", DCptrR->Name, DCptrI->Name);
        fprintf(stderr, "            to: rectifier -> %s-%s\n", DCptrR->Cont1, DCptrR->Cont2);
        fprintf(stderr, "                inverter  -> %s-%s\n", DCptrI->Cont1, DCptrI->Cont2);
      }
    }
  return (flag);
}

/* ----------------- Direction ---------------------------- */
int Direction(SparseMatrix *Mptr, VALUETYPE *vec, bool flag)
/* Find derivative -> dx/dlambda. */
{
  SparseMatrixElement *Jptr;
  AreaData *Aptr;
  ACbusData *ACptr;
  INDEX i, N;
  int m = 0, PgMax, PgMaxH;

  N = Mptr->n1;
  if (!flag) {
    if (ExistParameter('d'))
      fprintf(stderr, "Direction Factor Jacobian.\n");
    for (i = 1; i <= N; i++)
      for (vec[i] = dx[i] = 0, Jptr = Mptr->RowHead[i]; Jptr != nullptr; Jptr->Value = 0, Jptr = Jptr->RowNext)
        ;
    Aptr = ACFunJac(Mptr, &PgMax, false, true, false);
    if (DCFunJac(Mptr, false, true))
      return (-1);
    SVCFunJac(Mptr, false, true);     /* FACTS */
    TCSCFunJac(Mptr, false, true);    /* FACTS */
    STATCOMFunJac(Mptr, false, true); /* FACTS */
    if (flagH) {
      PgMaxH = HFunJac(false, true, Aptr, dx);
      if (NewCol->p[N] == 0)
        Jptr = Mptr->ColHead[N];
      else
        Jptr = Mptr->ColHead[NewCol->p[N]];
      while (Jptr != nullptr) {
        i = OldRow->p[Jptr->Row];
        if (i == 0)
          i = Jptr->Row;
        if (i < N)
          vec[i] = -Jptr->Value;
        Jptr = Jptr->ColNext;
      }
    }
    if (PgMax < 0 && (!flagH || (flagH && PgMaxH < 0))) {
      if (Aptr != nullptr) {
        fprintf(stderr, "\nError: Area %d %s does not have any spinning reserves.\n", Aptr->N, Aptr->Name);
        fprintf(stderr, "       Increase the maximum P generation in this "
                        "area, otherwise\n");
      } else {
        fprintf(stderr, "\nError: The system does not have any spinning reserves.\n");
        fprintf(stderr, "       Increase the maximum P generation in this "
                        "system, otherwise\n");
      }
      fprintf(stderr, "       the Jacobian matrix becomes singular.\n");
      fprintf(stderr, "Loading Factor -> %-6.3lf\n", lambda);
      WriteSolution(0, TrueParamStr(2), "Pg Max. Problems:");
      exit(1);
    }
    if (NewCol->p[N] != 0)
      m = factor(Mptr);
    else
      m = WARNINGEXIT;
  }
  if (m == WARNINGEXIT || flag) {
    if (ExistParameter('d'))
      fprintf(stderr, "Direction Order and Factor Jacobian.\n");
    DeleteJac(Mptr, NewRow, NewCol, OldRow, OldCol);
    Aptr = ACFunJac(Mptr, &PgMax, false, true, false);
    if (DCFunJac(Mptr, false, true))
      return (-1);
    SVCFunJac(Mptr, false, true);     /* FACTS */
    TCSCFunJac(Mptr, false, true);    /* FACTS */
    STATCOMFunJac(Mptr, false, true); /* FACTS */
    if (flagH) {
      for (i = 1; i <= N; vec[i] = dx[i] = 0, i++)
        ;
      PgMaxH = HFunJac(false, true, Aptr, dx);
      for (Jptr = Mptr->ColHead[N]; Jptr != nullptr; Jptr = Jptr->ColNext) {
        i = Jptr->Row;
        if (i < N)
          vec[i] = -Jptr->Value;
      }
    }
    if (PgMax < 0 && (!flagH || (flagH && PgMaxH < 0))) {
      if (Aptr != nullptr) {
        fprintf(stderr, "\nError: Area %d %s does not have any spinning reserves.\n", Aptr->N, Aptr->Name);
        fprintf(stderr, "       Increase the maximum P generation in this "
                        "area, otherwise\n");
      } else {
        fprintf(stderr, "\nError: The system does not have any spinning reserves.\n");
        fprintf(stderr, "       Increase the maximum P generation in this "
                        "system, otherwise\n");
      }
      fprintf(stderr, "       the Jacobian matrix becomes singular.\n");
      fprintf(stderr, "Loading Factor -> %-6.3lf\n", lambda);
      WriteSolution(0, TrueParamStr(2), "Pg Max. Problems:");
      exit(1);
    }
    SortRowsColumns(Mptr);
    if (factorns(Mptr, alpha, RowPartition, ColPartition, NewRow, NewCol, OldRow, OldCol)) {
      fprintf(stderr, "*** Singular Jacobian (possible voltage collapse, "
                      "contol or limit problems).\n");
      fprintf(stderr, "    Try changing the load levels, controls or limits, "
                      "or use the -F option.\n");
      fprintf(stderr, "Loading Factor -> %-6.3lf\n", lambda);
      WriteSolution(0, TrueParamStr(2), "Singular Jacobian:");
      exit(1);
    }
    SortRowsColumns(Mptr);
  }
  if (!flagH) {
    for (i = 1; i <= N; vec[i] = 0, i++)
      ;
    for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
      i = ACvar[ACptr->N];
      vec[i] = ACptr->Pnl * pow(ACptr->V, ACptr->a) + ACptr->Pzl * ACptr->V * ACptr->V;
      vec[i + 1] = ACptr->Qnl * pow(ACptr->V, ACptr->b) + ACptr->Qzl * ACptr->V * ACptr->V;
    }
  }
  repsolp(Mptr, vec, OldRow, NewCol);
  if (m == WARNINGEXIT)
    SD0 = 0;
  return (0);
}

/* --------------------------- ChangeParam --------------------------------- */
bool ChangeParam(void)
/* Change parameter for Homotopy method. */
{
  INDEX i;

  if (DxiMax > 1) {
    Dparam = 0;
    for (i = 1; i < Jac->n1; i++)
      Dx[i] = 0;
    LoadX0(false, true, true);
    if (BlPtr->N == Bl) {
      strcpy(BlPtr->Type, "B");
      if (BlPtr->Area != nullptr && BlPtr->Area->Slack == BlPtr)
        strcat(BlPtr->Type, "A");
      Bl = 0;
      BlPtr->Cont = BlPtr;
      if (ExistParameter('d'))
        fprintf(stderr, "Change param. to lambda.\n");
      return (true);
    } else {
      strcpy(BlPtr->Type, "BL");
      if (BlPtr->Area != nullptr && BlPtr->Area->Slack == BlPtr)
        strcat(BlPtr->Type, "A");
      Bl = BlPtr->N;
      BlPtr->Cont = nullptr;
      if (ExistParameter('d'))
        fprintf(stderr, "Change param. to V at %s %d %s.\n", BlPtr->Type, BlPtr->Num, BlPtr->Name);
      return (true);
    }
  } else
    return (false);
}

/* ------------------------- VoltageSensFactor --------------------------- */
void VoltageSensFactor(VALUETYPE *vec, bool first)
/* Calculate Voltage Sensitivity Factor and
   Tangent Vector Index and rank for a bus */
{
  ACbusData *ACptr;
  VALUETYPE val;
  INDEX i, l, N = 0, MaxBus = 0;
  ACranklist *RankList, *Ptr, *PtrMax, *PrevPtr, *PrevPtrMax;

  if (first)
    TVIbus = IntegerParameter('f', 0, 1, 9999);
  RankList = nullptr;
  for (VSFinf = VSFone = TVI = 0, ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    if (ACptr->Cont != nullptr) {
      i = ACvar[ACptr->N];
      val = fabs(vec[i + 1]);
      Ptr = RankList;
      RankList = new ACranklist;
      RankList->AC = ACptr;
      RankList->val = val;
      RankList->Next = Ptr;
      N++;
      if (val > VSFone) {
        VSFone = val;
        VSFbus = MaxBus = ACptr->Num;
      }
      VSFinf = VSFinf + val;
      if (TVIbus == ACptr->Num)
        TVI = val;
    }
  }
  if (first && TVI == 0) {
    TVIbus = MaxBus;
    TVI = VSFone;
  }
  /*    Rank buses based on Tangent Vector   */
  if (TVI != 0)
    for (l = 1; l <= N; l++) {
      if (RankList == nullptr || RankList->Next == nullptr)
        break;
      for (val = -0.1, PrevPtr = Ptr = RankList; Ptr != nullptr; PrevPtr = Ptr, Ptr = Ptr->Next) {
        if (Ptr->val > val) {
          val = Ptr->val;
          PtrMax = Ptr;
          PrevPtrMax = PrevPtr;
        }
      }
      if (PtrMax != nullptr) {
        if (PrevPtrMax != PtrMax)
          PrevPtrMax->Next = PtrMax->Next;
        else
          RankList = PtrMax->Next;
        if (ExistParameter('d')) {
          fprintf(stderr, "Ranked List: ");
          fprintf(stderr, "%4d %12s %-11.5g\n", PtrMax->AC->Num, PtrMax->AC->Name, PtrMax->val);
        }
        if (TVIbus == PtrMax->AC->Num) {
          TVIrank = l;
          break;
        }
        delete PtrMax;
      }
    }
  for (Ptr = RankList; Ptr != nullptr;) {
    PrevPtr = Ptr->Next;
    delete Ptr;
    Ptr = PrevPtr;
  }
}

/* --------------------------  CheckVIlimits ------------------------ */
bool CheckVIlimits(bool flagVoltage, bool flagCurrent) {
  ACbusData *ACptr;
  ElementData *Eptr;
  ElementList *ELptr;
  VALUETYPE Pij, Qij, Vi, di, Vj, dj, Iij, G, B, Gi, Bi, Pl, Ql, Gp, Bp, Gj, Bj, Pji, Qji, Iji;
  bool FlagStopContinuation = false;

  TotalPl = TotalQl = TotalPg = TotalQg = 0;
  for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    TotalPg = TotalPg + ACptr->PG;
    Pl = (ACptr->Pn + lambda * ACptr->Pnl) * pow(ACptr->V, ACptr->a) + (ACptr->Pz + lambda * ACptr->Pzl) * ACptr->V * ACptr->V;
    ACptr->PL = Pl;
    TotalPl = TotalPl + Pl;
    Ql = (ACptr->Qn + lambda * ACptr->Qnl) * pow(ACptr->V, ACptr->b) + (ACptr->Qz + lambda * ACptr->Qzl) * ACptr->V * ACptr->V;
    ACptr->QL = Ql;
    TotalQl = TotalQl + Ql;
    TotalQg = TotalQg + ACptr->Qg;
    if (ACptr->CheckVlimits && ACptr->Vlmax > ACptr->Vlmin && flagVoltage) {
      if (ACptr->V <= ACptr->Vlmin) {
        fprintf(stderr,
                "***Warning: Bus %d %s has reached or exceeded its minimum V "
                "limit.\n",
                ACptr->Num, ACptr->Name);
        ACptr->CheckVlimits = false;
        FlagStopContinuation = true;
      } else if (ACptr->V >= ACptr->Vlmax) {
        fprintf(stderr,
                "***Warning: Bus %d %s has reached or exceeded its maximum V "
                "limit.\n",
                ACptr->Num, ACptr->Name);
        ACptr->CheckVlimits = false;
        FlagStopContinuation = true;
      }
    }
    if (flagCurrent)
      for (ELptr = ACptr->Elem; ELptr != nullptr; ELptr = ELptr->Next) {
        Eptr = ELptr->Eptr;
        if (Eptr->CheckIlimits && Eptr->Imax > 0) {
          Vi = Eptr->From->V;
          di = Eptr->From->Ang;
          Vj = Eptr->To->V;
          dj = Eptr->To->Ang;
          G = (Eptr->G * cos(Eptr->Ang) - Eptr->B * sin(Eptr->Ang)) * Eptr->Tap;
          B = (Eptr->G * sin(Eptr->Ang) + Eptr->B * cos(Eptr->Ang)) * Eptr->Tap;
          Gi = (Eptr->G1 + Eptr->G) * pow(Eptr->Tap, 2.0) - G;
          Bi = (Eptr->B1 + Eptr->B) * pow(Eptr->Tap, 2.0) - B;
          Gp = (Eptr->G * cos(Eptr->Ang) + Eptr->B * sin(Eptr->Ang)) * Eptr->Tap;
          Bp = (-Eptr->G * sin(Eptr->Ang) + Eptr->B * cos(Eptr->Ang)) * Eptr->Tap;
          Gj = Eptr->G + Eptr->G2 - Gp;
          Bj = Eptr->B + Eptr->B2 - Bp;
          Pij = Vi * Vi * (Gi + G) - Vi * Vj * (G * cos(di - dj) + B * sin(di - dj));
          Qij = -Vi * Vi * (Bi + B) - Vi * Vj * (G * sin(di - dj) - B * cos(di - dj));
          Iij = sqrt(Pij * Pij + Qij * Qij) / Vi;
          Pji = Vj * Vj * (Gj + Gp) - Vi * Vj * (Gp * cos(dj - di) + Bp * sin(dj - di));
          Qji = -Vj * Vj * (Bj + Bp) - Vi * Vj * (Gp * sin(dj - di) - Bp * cos(dj - di));
          Iji = sqrt(Pji * Pji + Qji * Qji) / Vj;
          if (Iij >= Eptr->Imax) {
            fprintf(stderr, "***Warning: Element %d %s to %d %s has reached\n", Eptr->From->Num, Eptr->From->Name, Eptr->To->Num, Eptr->To->Name);
            fprintf(stderr, "            or exceeded its maximum I limit at "
                            "its first bus.\n");
            Eptr->CheckIlimits = false;
            FlagStopContinuation = true;
          } else if (Iji >= Eptr->Imax) {
            fprintf(stderr, "***Warning: Element %d %s to %d %s has reached\n", Eptr->From->Num, Eptr->From->Name, Eptr->To->Num, Eptr->To->Name);
            fprintf(stderr, "            or exceeded its maximum I limit at "
                            "its second bus.\n");
            Eptr->CheckIlimits = false;
            FlagStopContinuation = true;
          }
        }
      }
  }
  return (FlagStopContinuation);
}

/* ------------------------- AddLoad --------------------------- */
void AddLoad(void)
/* Add initial loading to load buses  */
{
  ACbusData *ACptr;

  for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    ACptr->Pn = ACptr->Pn + lambda * ACptr->Pnl;
    ACptr->Pz = ACptr->Pz + lambda * ACptr->Pzl;
    ACptr->Qn = ACptr->Qn + lambda * ACptr->Qnl;
    ACptr->Qz = ACptr->Qz + lambda * ACptr->Qzl;
  }
}

/* --------------------------- Homot --------------------------------- */
int Homot(void)
/* Build voltage profiles. */
{
  VALUETYPE NDx, cons, PreviousLambda, StoplambdaVal, FirstVoltageSumUp = Tol;
  int iter, iterp, CountCyclingSteps, CountVupSteps, NumHVsolutions, Steps, IndicesSteps, MaxSteps, i, N, n, CountStepsAfterVup = 0;
  bool flagLimits = false, flagSTOP = false, flagFACTSlimits = false, flagt = false, stopValueSet = false, flagMakeVlist = true,
       flagCountVupAfterFirstLimit = false, flagD = false, ForceflagCountVupAfterFirstLimit = true, flagFirstStep = true, flagCycling = false,
       flagRlimits = false, flagQlimits = false, flagVlimits = false, flagVIlimits = false, flagDClimits = false, flagVoltage = false, flagCurrent = false,
       flagSVClimits = false, flagTCSClimits = false, flagSTATCOMlimits = false; /* FACTS */

  BlPtr = nullptr;
  Dparam = param0 = 0;
  AngTr = 180 / PI;
  Kh = 1;
  RealParameter('k', &Kh, 1e-4, 1e4);
  Htol = 0.000001;
  RealParameter('o', &Htol, 1e-10, 1e-1);
  if (ExistParameter('H'))
    OutputHomot = (FILE *)OpenOutput(NameParameter('H'));
  else
    OutputHomot = (FILE *)OpenOutput(NameParameter('c'));
  N = Jac->n1;
  RedNumEq = N - 1;
  CountSteps = 0;
  NumSteps = IntegerParameter('U', 10, 2, 100);
  iter = Pflow(1, false, true, true);
  if (iter < 0) {
    fclose(OutputHomot);
    return (iter);
  } else
    fprintf(stderr, "**** Initial Power Flow Solved ****    \n\n");
  lambda = SD0 = 0;
  RealParameter('L', &lambda, -1e6, 1e6);
  lambda_o = lambda;
  NumHVsolutions = IntegerParameter('2', 5, 2, 100);
  flagVoltage = ExistParameter('7');
  flagCurrent = ExistParameter('8');
  if (lambda != 0) {
    if ((ExistParameter('D') && (!NullName(NameParameter('D')))) || flagVloads)
      flagD = true;
    iter = Pflow(iter, flagD, true, false);
    if (iter > 0) {
      fprintf(stderr, "**** Initial Lambda Case Solved ****    ");
      fprintf(stderr, "Starting loading factor -> %-10.6lg\n\n", lambda);
      AddLoad();
    } else {
      fclose(OutputHomot);
      return (iter);
    }
    lambda = 0;
  }
  flagt = DCsetup();
  field = IntegerParameter('O', 6, 6, 10);
  MaxSteps = IntegerParameter('z', 1000, 1, 9999);
  StoplambdaVal = lambda;
  if (ExistParameter('u')) {
    ZeroDx = 0.001;
    flagReduceSystem = true;
    DxZero = new bool[N];
    RealParameter('u', &ZeroDx, 0., 0.2);
    if (Kh < 1)
      ZeroDx = Kh * ZeroDx;
  } else {
    flagReduceSystem = false;
    ZeroDx = 0;
    DxZero = nullptr;
  }
  for (i = 1; i <= N; i++) {
    if (DxZero != nullptr)
      DxZero[i] = false;
  }
  for (flagL = flagR = true, CountCyclingSteps = CountVupSteps = Steps = 1;;) {
  Homotopy:
    if (ExistParameter('d'))
      fprintf(stderr, "\nContinuation Step=%d\n", Steps);
    Dparam = 0;
    PreviousLambda = lambda;
    /* Predictor Step */
    if (Direction(Jac, Dx, flagt) < 0) {
      fclose(OutputHomot);
      return (-iter);
    }
    if (SD0 == 0)
      SD0 = DetSign;
    for (NDx = SF = 0, i = 1; i <= N - 1; i++) {
      if (fabs(Dx[i]) > NDx)
        NDx = fabs(Dx[i]);
      SF = SF + fabs(Dx[i]);
    }
    if (flagFirstStep) {
      if (ExistParameter('u') && ZeroDx > 0)
        flagReducedContinuation = true;
      flagFirstStep = false;
    }
    if (ExistParameter('f'))
      VoltageSensFactor(Dx, flagMakeVlist);
    if (flagMakeVlist) {
      MakeVlist(OutputHomot);
      CheckVIlimits(flagVoltage, flagCurrent);
      VoltProf(true, OutputHomot);
      if (ExistParameter('0')) {
        IndicesSteps = Steps;
        WriteQmatrix(IndicesSteps, Dx);
      }
      flagMakeVlist = false;
    }
    Steps++;
    if (ExistParameter('f')) {
      if (Kh < 0) {
        VSFone = VSFinf = SF = TVI = 0;
        VSFbus = TVIrank = 0;
      }
      Print(OutputHomot, 0, 6, 4, VSFone);
      fprintf(OutputHomot, "    %6d    ", VSFbus);
      Print(OutputHomot, 0, 6, 4, VSFinf);
      fprintf(OutputHomot, "    ");
      Print(OutputHomot, 0, 6, 4, SF);
      fprintf(OutputHomot, "    ");
      if (TVI != 0)
        Print(OutputHomot, 0, 6, 4, 1. / TVI);
      else
        fprintf(OutputHomot, "%6d", 0);
      fprintf(OutputHomot, "    %6d", TVIrank);
    }
    fprintf(OutputHomot, "\n");
    if (NDx) {
      if (ExistParameter('d'))
        fprintf(stderr, "Lambda->%lf  SD0=%2d  Det.Sign=%2d  Kh=%lf\n", lambda, SD0, DetSign, Kh);
      if (Bl) {
        cons = -1;
        Kh = -fabs(Kh);
        SD0 = 0;
      } else {
        if (SD0 * DetSign < 0) {
          SD0 = DetSign;
          Kh = -Kh;
        }
        cons = Kh;
      }
      if (ExistParameter('d'))
        fprintf(stderr, "%38s Kh=%lf\n", "", Kh);
      Dparam = cons / NDx;
      for (i = 1; i <= N - 1; i++)
        Dx[i] = Dparam * Dx[i];
      if (flagSTOP) {
        if (ExistParameter('Z') && !NullName(NameParameter('Z')))
          PrintDirection('Z', Dx, 1.);
        break;
      }

      if (flagReducedContinuation) {
        /* Stop system reduction when V goes up, and start again
           after V goes down for 5 steps. */
        CountSteps++;
        if (VoltageSum > FirstVoltageSumUp && Kh < 0)
          flagReduceSystem = false;
        if (!flagReduceSystem) {
          if (VoltageSum <= FirstVoltageSumUp && Kh < 0)
            CountStepsAfterVup++;
          if (CountStepsAfterVup >= 5)
            flagReduceSystem = true;
        }
      }

      /* Update variables from predictor */
      cons = LoadX0(true, true, true);
      if (!CountVupSteps)
        FirstVoltageSumUp = VoltageSum0;
      if (flagReducedContinuation && flagReduceSystem && CountSteps >= NumSteps) {
        CountSteps = 0;
        if (NewNumEq < RedNumEq) {
          RedNumEq = NewNumEq;
          SD0 = 0;
        }
      }
      if (ExistParameter('d'))
        fprintf(stderr, "cons=%lf Dparam=%lf\n", cons, Dparam);
      if (ExistParameter('d'))
        fprintf(stderr, "Num.Steps=%d  Num.Eqs.System=%d   Num.Eqs.Reduced=%d\n", CountSteps, N - 1, RedNumEq);
      if (!flagt) {
        if (ExistParameter('H'))
          flagt = ChangeParam();
      } else
        flagt = false;
      if (flagt)
        goto Homotopy;

      /*  Enforce Limits */
      if (cons > 0) {
        Dparam = fabs(1 - cons + Htol * 0.1) * Dparam;
        if (ExistParameter('d'))
          fprintf(stderr, "(1-%lf)Dparam->%lf\n", cons, Dparam);
        for (i = 1; i <= N - 1; i++)
          Dx[i] = fabs(1 - cons + Htol * 0.1) * Dx[i];
        LoadX0(false, true, true);
        if ((!ExistParameter('G') || !flagCountVupAfterFirstLimit || !CountVupSteps) && fabs(Dparam) <= Htol) {
          flagRlimits = CheckRlimits();
          flagVlimits = CheckVlimits();
          flagQlimits = CheckQlimits();
          flagDClimits = ChangeDCmode();
          flagSVClimits = ChangeSVCmode();         /* FACTS */
          flagTCSClimits = ChangeTCSCmode();       /* FACTS */
          flagSTATCOMlimits = ChangeSTATCOMmode(); /* FACTS */
          if (flagRlimits || flagVlimits || flagQlimits || flagDClimits || flagSVClimits || flagTCSClimits || flagSTATCOMlimits)
            flagLimits = true; /* FACTS */
          else
            flagLimits = false;
          if (flagDClimits || flagSVClimits || flagTCSClimits || flagSTATCOMlimits)
            flagFACTSlimits = true; /* FACTS */
          else
            flagFACTSlimits = false; /* FACTS */
          if (flagLimits) {
            SD0 = 0;
            if (ForceflagCountVupAfterFirstLimit)
              flagCountVupAfterFirstLimit = true;
          }
        } else
          flagLimits = flagFACTSlimits = false; /* FACTS */
      } else
        flagLimits = flagFACTSlimits = false; /* FACTS */

      /* Corrector step */
      for (n = 1; n <= 10; n++) {
        if (n == 10)
          flagL = flagR = false;
        iterp = Pflow(iter, flagFACTSlimits, flagLimits, false);
        if (flagL == false)
          flagL = flagR = true;
        /* if ((iterp>0 && lambda>-1e-4) || fabs(Dparam)<=1e-6)  break;*/
        if (iterp > 0 && lambda > -1e-4)
          break;
        Dparam = Dparam / 2;
        for (i = 1; i <= N - 1; i++)
          Dx[i] = Dx[i] / 2;
        LoadX0(false, true, true);
        if (ExistParameter('d'))
          fprintf(stderr, "Dparam/2^%d->%lf\n", n, Dparam);
      }

      /* Turning point computations */
      if (flagCountVupAfterFirstLimit && VoltageSum > FirstVoltageSumUp && Kh < 0)
        CountVupSteps++;
      else
        CountVupSteps = 0;
      if (ExistParameter('d'))
        fprintf(stderr, "V summation->%lf     First V sum Up->%lf    Vup Steps->%d\n", VoltageSum, FirstVoltageSumUp, CountVupSteps);
      if (ExistParameter('G') && flagCountVupAfterFirstLimit && (CountVupSteps >= NumHVsolutions || (CountCyclingSteps > 0 && iterp < 0))) {
        Kh = -Kh;
        ForceflagCountVupAfterFirstLimit = flagCountVupAfterFirstLimit = false;
        if (iterp < 0)
          iterp = -iterp;
      }

      /* Stop and output commands */
      // when lamda hits its max, set the stop value
      if (lambda < PreviousLambda && !stopValueSet) {
        stopValueSet = true;
        RealParameter('S', &StoplambdaVal, 0., 1.);
        StoplambdaVal = PreviousLambda * StoplambdaVal;
        if (ExistParameter('E') && !NullName(NameParameter('E'))) {
          if (ExistParameter('d'))
            fprintf(stderr, "Write right e-vector at lambda->%lf\n", PreviousLambda);
          PrintDirection('E', Dx, 1.0);
        }
      }
      // Sometimes there is a small decrease in lamda before it actually hits
      // its max the following will reset the stopValueSet flag incase that
      // happened
      if (stopValueSet && lambda > PreviousLambda) {
        stopValueSet = false;
      }

      if ((MaxSteps > 0 && Steps >= MaxSteps) || (lambda - StoplambdaVal) <= 1e-5)
        flagSTOP = true;
      if (iterp < 0) {
        iter = -iterp;
        fprintf(stderr, "\n *** The case diverges (possible AC/DC control problems).\n");
        fprintf(stderr, "     Try reducing the step size (-k option).\n");
        fprintf(stderr, "Loading Factor -> %-6.3lf\n", lambda);
        WriteSolution(--iter, TrueParamStr(2), "Divergence:");
        exit(1);
      } else
        iter = iterp;
      flagVIlimits = CheckVIlimits(flagVoltage, flagCurrent);
      VoltProf(false, OutputHomot);
      if (flagVIlimits)
        break;
      if (ExistParameter('0') && Kh > 0) {
        IndicesSteps = Steps;
        WriteQmatrix(IndicesSteps, Dx);
      }
      if (fabs(lambda - param0) <= Htol)
        CountCyclingSteps++;
      else
        CountCyclingSteps = 0;
      if (CountCyclingSteps >= 5 && !CountVupSteps) {
        if (flagCycling) {
          fprintf(stderr, "\n *** Cycling due to possible AC/DC/FACTS control "
                          "problems.  Try using the debug (-d)\n");
          fprintf(stderr, "     option to detect the limit that is creating "
                          "the problem.\n");
          fprintf(stderr, "Loading Factor -> %-6.3lf\n", lambda);
          WriteSolution(iter, TrueParamStr(2), "Divergence:");
          exit(1);
        } else
          Kh = -Kh;
        flagCycling = true;
      }
    } else
      break;
  }
  if (ExistParameter('m'))
    MatlabV(OutputHomot);
  if (OutputHomot != nullptr)
    fclose(OutputHomot);
  if (ExistParameter('O'))
    TEFMatlabFiles();
  if (ExistParameter('0'))
    IndicesMatlab(IndicesSteps);
  return (iter);
}
