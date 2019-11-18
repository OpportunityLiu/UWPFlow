/* Write Jacobian. */

#include "write.h"

/* ---------- Global Variables --------- */
INDEX TFnum, TFbus;
char TFname[13];
extern VALUETYPE K;
extern INDEX *ACvar;

/* --------------- WriteJac ---------------------- */
void WriteJac(void) {
  SparseMatrixElement *Jptr;
  ACbusData *ACptr, *ACptrp;
  AClist *Lptr;
  DCbusData *DCptrR, *DCptrI, *DCptr;
  ElementData *Eptr;
  ElementList *ELptr;
  char Namebase[80], Name[80], str[80], type[2];
  FILE *OutFile, *OutFilep;
  int i, j, k, l;
  INDEX I, J, N, Nvar;
  bool flag = false;
  SVCbusData *SVCptr;         /* FACTS */
  TCSCbusData *TCSCptr;       /* FACTS */
  STATCOMbusData *STATCOMptr; /* FACTS */

  if (ExistParameter('J'))
    strcpy(Namebase, NameParameter('J'));
  else
    strcpy(Namebase, NameParameter('j'));
  if (NullName(Namebase))
    return;
  if (ExistParameter('J') && Bl)
    for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
      if (ACptr->Num == Bl) {
        strcpy(ACptr->Type, "B");
        if (ACptr->Area != nullptr && ACptr->Area->Slack == ACptr)
          strcat(ACptr->Type, "A");
        Bl = 0;
        ACptr->Cont = ACptr;
        break;
      }
    }
  flag = (!ExistParameter('J')) && flagPoC;
  DeleteJac(Jac, NewRow, NewCol, OldRow, OldCol);
  RowPer = NewRow;
  ColPer = NewCol;
  Nvar =
      NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar + 7 * Nstatcom; /*  FACTS  */
  ACFunJac(Jac, nullptr, true, true, false);
  DCFunJac(Jac, true, true);
  SVCFunJac(Jac, true, true);     /*  FACTS  */
  TCSCFunJac(Jac, true, true);    /*  FACTS  */
  STATCOMFunJac(Jac, true, true); /*  FACTS  */
  if (flagH) {
    Nvar++;
    HFunJac(true, true, nullptr, Dx);
  } else if (flag) {
    Nvar = 2 * Nvar + 1;
    ACFunHes(true, true);
    DCFunHes(true, true);
    SVCFunHes(true, true);     /* FACTS  */
    TCSCFunHes(true, true);    /* FACTS  */
    STATCOMFunHes(true, true); /* FACTS  */
  }
  SortRowsColumns(Jac);
  strcpy(Name, Namebase);
  strcat(Name, ".jac");
  OutFile = OpenOutput(Name);
  fprintf(OutFile, "%d %d\n", Nvar, Nvar);
  for (i = 1; i <= Nvar; i++) {
    for (Jptr = Jac->RowHead[i]; Jptr != nullptr; Jptr = Jptr->RowNext)
      fprintf(OutFile, "%4d %4d %-11.5g\n", Jptr->Row, Jptr->Col, Jptr->Value);
  }
  fprintf(OutFile, "%4d %4d %-11.5g\n", 0, 0, 0.);
  fclose(OutFile);
  strcpy(Name, Namebase);
  strcat(Name, ".var");
  OutFile = OpenOutput(Name);
  strcpy(Name, Namebase);
  strcat(Name, ".mis");
  OutFilep = OpenOutput(Name);
  fprintf(OutFile, "%d 1\n", Nvar);
  fprintf(OutFilep, "%d 1\n", Nvar);
  for (i = 0, ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    if (ACptr->Cont != nullptr) {
      if (strpbrk(ACptr->Type, "S")) {
        sprintf(str, "kg%-d", ACptr->Num);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Kg);
      } else {
        sprintf(str, "d%-d", ACptr->Num);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Ang);
      }
      sprintf(str, "dP%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      sprintf(str, "V%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->V);
      sprintf(str, "dQ%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    } else if (QRcont && strpbrk(ACptr->Type, "C")) {
      sprintf(str, "d%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Ang);
      sprintf(str, "dP%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      for (Lptr = ACptr->ContBus; Lptr != nullptr; Lptr = Lptr->Next) {
        ACptrp = Lptr->AC;
        if (strpbrk(ACptrp->cont, "V"))
          break;
      }
      sprintf(str, "Qr%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Qr);
      sprintf(str, "dQ%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    } else if (Rcont && strpbrk(ACptr->Type, "T")) {
      sprintf(str, "d%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Ang);
      sprintf(str, "dP%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      for (ELptr = ACptr->Reg; ELptr != nullptr; ELptr = ELptr->Next) {
        Eptr = ELptr->Eptr;
        I = Eptr->From->Num;
        J = Eptr->To->Num;
        if (!strcmp(Eptr->Type, "R"))
          break;
      }
      sprintf(str, "1/t%-d_%-d", I, J);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, Eptr->Tap);
      sprintf(str, "dQ%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    } else if (strpbrk(ACptr->Type, "L")) {
      sprintf(str, "d%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Ang);
      sprintf(str, "dP%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, "l", lambda);
      sprintf(str, "dQ%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    } else if (strpbrk(ACptr->Type, "Q") || strpbrk(ACptr->Type, "V") ||
               (!QRcont && strpbrk(ACptr->Type, "G"))) {
      if (strpbrk(ACptr->Type, "S")) {
        sprintf(str, "kg%-d", ACptr->Num);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Kg);
      } else {
        sprintf(str, "d%-d", ACptr->Num);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Ang);
      }
      sprintf(str, "dP%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      sprintf(str, "Qg%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Qg);
      sprintf(str, "dQ%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    } else if (strpbrk(ACptr->Type, "Z")) {
      sprintf(str, "d%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Ang);
      sprintf(str, "dP%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      sprintf(str, "Qz%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Qg);
      sprintf(str, "dQ%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    } else if (strpbrk(ACptr->Type, "S")) {
      sprintf(str, "kg%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Kg);
      sprintf(str, "dP%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      sprintf(str, "Qg%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Qg);
      sprintf(str, "dQ%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    }
    if (Acont && strpbrk(ACptr->Type, "A")) {
      sprintf(str, "kg%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Kg);
      sprintf(str, "dPA%-d", ACptr->Area->N);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    }
    if (PQcont)
      for (ELptr = ACptr->Reg; ELptr != nullptr; ELptr = ELptr->Next) {
        Eptr = ELptr->Eptr;
        if (strpbrk(Eptr->Type, "PQNM")) {
          if (Eptr->From == ACptr) {
            I = Eptr->From->Num;
            J = Eptr->To->Num;
          } else {
            J = Eptr->From->Num;
            I = Eptr->To->Num;
          }
          if (!strcmp(Eptr->Type, "RP")) {
            sprintf(str, "a%-d_%-d", I, J);
            fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, Eptr->Ang);
            sprintf(str, "dP%-d_%-d", I, J);
            fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
          } else if (strpbrk(Eptr->Type, "PM")) {
            sprintf(str, "P%-d_%-d", I, J);
            fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, Eptr->Cvar);
            sprintf(str, "dP%-d_%-d", I, J);
            fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
          } else if (!strcmp(Eptr->Type, "RQ")) {
            sprintf(str, "1/t%-d_%-d", I, J);
            fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, Eptr->Tap);
            sprintf(str, "dQ%-d_%-d", I, J);
            fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
          } else {
            sprintf(str, "Q%-d_%-d", I, J);
            fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, Eptr->Cvar);
            sprintf(str, "dQ%-d_%-d", I, J);
            fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
          }
        }
      }
    if (ACptr->Gen != nullptr) {
      i = ACptr->Gen->Nvar;
      sprintf(str, "dPg%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dQg%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dEq%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dEd%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dVd%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dVq%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dId%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dIq%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dVr%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dVi%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      sprintf(str, "dIa%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", ++i, str, dF[i]);
      i = ACptr->Gen->Nvar;
      if (strpbrk(ACptr->cont, "E")) {
        sprintf(str, "Qg%-d", ACptr->Num);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Qg);
      } else {
        sprintf(str, "Eq%-d", ACptr->Num);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Eq);
      }
      sprintf(str, "dg%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->dg);
      sprintf(str, "Vr%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Vr);
      sprintf(str, "Vi%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Vi);
      sprintf(str, "Ir%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Ir);
      sprintf(str, "Ii%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Ii);
      sprintf(str, "Vq%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Vq);
      sprintf(str, "Vd%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Vd);
      sprintf(str, "Iq%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Iq);
      sprintf(str, "Id%-d", ACptr->Num);
      fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Id);
      if (strpbrk(ACptr->cont, "I")) {
        sprintf(str, "Qg%-d", ACptr->Num);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Qg);
      } else {
        sprintf(str, "Ia%-d", ACptr->Num);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++i, str, ACptr->Gen->Ia);
      }
    }
  }
  for (k = 0, DCptrR = dataPtr->DCbus; DCptrR != nullptr;
       DCptrR = DCptrR->Next) {
    DCptrI = DCptrR->To;
    if (!strcmp(DCptrR->Type, "R")) {
      for (k++, l = 1; l <= 11; l++) {
        sprintf(str, "Fdc%-d_%-d", k, l);
        i++;
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      }
      for (l = i - 11, j = 1; j <= 2; j++) {
        if (j == 1) {
          DCptr = DCptrR;
          strcpy(type, "r");
        } else {
          DCptr = DCptrI;
          strcpy(type, "i");
        }
        if (strcmp(DCptr->Cont1, "VD") && strcmp(DCptr->Cont2, "VD")) {
          sprintf(str, "Vd%1s%-d", type, k);
          fprintf(OutFile, "%4d %8s %-11.5g\n", ++l, str, DCptr->Vd);
        }
        if (strcmp(DCptr->Cont1, "AT") && strcmp(DCptr->Cont2, "AT")) {
          sprintf(str, "t%1s%-d", type, k);
          fprintf(OutFile, "%4d %8s %-11.5g\n", ++l, str, DCptr->Tap);
        }
        if (strcmp(DCptr->Cont1, "AL") && strcmp(DCptr->Cont2, "AL")) {
          sprintf(str, "al%1s%-d", type, k);
          fprintf(OutFile, "%4d %8s %-11.5g\n", ++l, str, DCptr->Alfa);
        }
        if (strcmp(DCptr->Cont1, "GA") && strcmp(DCptr->Cont2, "GA")) {
          sprintf(str, "ga%1s%-d", type, k);
          fprintf(OutFile, "%4d %8s %-11.5g\n", ++l, str, DCptr->Gamma);
        }
        sprintf(str, "s%1s%-d", type, k);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++l, str, DCptr->MVA);
        if (strcmp(DCptr->Cont1, "PA") && strcmp(DCptr->Cont2, "PA")) {
          sprintf(str, "P%1s%-d", type, k);
          fprintf(OutFile, "%4d %8s %-11.5g\n", ++l, str, DCptr->P);
        }
        if (strcmp(DCptr->Cont1, "QA") && strcmp(DCptr->Cont2, "QA")) {
          sprintf(str, "Q%1s%-d", type, k);
          fprintf(OutFile, "%4d %8s %-11.5g\n", ++l, str, DCptr->Q);
        }
      }
      if (strcmp(DCptrR->Cont1, "ID") && strcmp(DCptrR->Cont2, "ID") &&
          strcmp(DCptrI->Cont1, "ID") && strcmp(DCptrI->Cont2, "ID")) {
        sprintf(str, "Id%-d", k);
        fprintf(OutFile, "%4d %8s %-11.5g\n", ++l, str, DCptrR->Id);
      }
    }
  }

  /*   FACTS   */
  for (k = 0, SVCptr = dataPtr->SVCbus; SVCptr != nullptr;
       SVCptr = SVCptr->Next) {
    k++;
    l = 0;
    sprintf(str, "Qsvc%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, SVCptr->Qsvc);
    sprintf(str, "Fsvc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "Bv%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, SVCptr->Bv);
    sprintf(str, "Fsvc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    if (!strcmp(SVCptr->Cont, "AL")) {
      sprintf(str, "alpha%-d", k);
      i++;
      fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, SVCptr->alpha_svc);
    } else {
      sprintf(str, "Vrefc%-d", k);
      i++;
      fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, SVCptr->Vvar);
    }
    sprintf(str, "Fsvc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
  }

  for (k = 0, TCSCptr = dataPtr->TCSCbus; TCSCptr != nullptr;
       TCSCptr = TCSCptr->Next) {
    k++;
    l = 0;
    sprintf(str, "Ptcsc%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, TCSCptr->Ptcsc);
    sprintf(str, "Ftcsc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "Qtcsck%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, TCSCptr->Qtcsck);
    sprintf(str, "Ftcsc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "Qtcscm%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, TCSCptr->Qtcscm);
    sprintf(str, "Ftcsc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "Be%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, TCSCptr->Be);
    sprintf(str, "Ftcsc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "alpha%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, TCSCptr->alpha_tcsc);
    sprintf(str, "Ftcsc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "Itcsc%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, TCSCptr->Itcsc);
    sprintf(str, "Ftcsc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "delta%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, TCSCptr->delta_t);
    sprintf(str, "Ftcsc%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
  }
  for (k = 0, STATCOMptr = dataPtr->STATCOMbus; STATCOMptr != nullptr;
       STATCOMptr = STATCOMptr->Next) {
    k++;
    l = 0;
    if (!strcmp(STATCOMptr->Cont, "PW") || !strcmp(STATCOMptr->Cont, "AL")) {
      sprintf(str, "Istat%-d", k);
      i++;
      fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, STATCOMptr->I);
    } else {
      sprintf(str, "Vrefc%-d", k);
      i++;
      fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, STATCOMptr->Vvar);
    }
    sprintf(str, "Vcont%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "theta%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, STATCOMptr->theta);
    sprintf(str, "Cont%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "Vdc%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, STATCOMptr->Vdc);
    sprintf(str, "Loss%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "k%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, STATCOMptr->k);
    sprintf(str, "ReS%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "alpha%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, STATCOMptr->alpha);
    sprintf(str, "ImS%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "Pstat%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, STATCOMptr->P);
    sprintf(str, "Pstat%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
    sprintf(str, "Qstat%-d", k);
    i++;
    fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, STATCOMptr->Q);
    sprintf(str, "Qstat%-d_%-d", k, ++l);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
  }
  /* END FACTS */

  if (flagH) {
    fprintf(OutFile, "%4d %8s %-11.5g\n", Jac->n1, "l", lambda);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", Jac->n1, "dH", dF[Jac->n1]);
  } else if (flag) {
    N = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar + 7 * Nstatcom; /* FACTS */
    for (i = N, ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
      i++;
      sprintf(str, "w%-d", i - N);
      fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
      sprintf(str, "gP%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      i++;
      sprintf(str, "w%-d", i - N);
      fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
      sprintf(str, "gQ%-d", ACptr->Num);
      fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      if (Acont && strpbrk(ACptr->Type, "A")) {
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gPA%-d", ACptr->Area->N);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      }
      if (PQcont)
        for (ELptr = ACptr->Reg; ELptr != nullptr; ELptr = ELptr->Next) {
          Eptr = ELptr->Eptr;
          if (strpbrk(Eptr->Type, "PQNM")) {
            if (Eptr->From == ACptr) {
              I = Eptr->From->Num;
              J = Eptr->To->Num;
            } else {
              J = Eptr->From->Num;
              I = Eptr->To->Num;
            }
            if (!strcmp(Eptr->Type, "PM")) {
              i++;
              sprintf(str, "w%-d", i - N);
              fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
              sprintf(str, "gP%-d_%-d", I, J);
              fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
            } else {
              i++;
              sprintf(str, "w%-d", i - N);
              fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
              sprintf(str, "gQ%-d_%-d", I, J);
              fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
            }
          }
        }
      if (ACptr->Gen != nullptr) {
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gPg%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gQg%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gEq%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gEd%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gVd%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gVq%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gId%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gIq%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gVr%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gVi%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        i++;
        sprintf(str, "w%-d", i - N);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gIa%-d", ACptr->Num);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      }
    }
    for (k = 0, DCptrR = dataPtr->DCbus; DCptrR != nullptr;
         DCptrR = DCptrR->Next) {
      if (!strcmp(DCptrR->Type, "R")) {
        for (k++, l = 1; l <= 11; l++) {
          i++;
          sprintf(str, "w%-d", i - N);
          fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
          sprintf(str, "gdc%-d_%-d", k, l);
          fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
        }
      }
    }

    /* FACTS */
    for (k = 0, SVCptr = dataPtr->SVCbus; SVCptr != nullptr;
         SVCptr = SVCptr->Next) {
      for (k++, l = 1; l <= 3; l++) {
        i++;
        sprintf(str, "wsvc%-d_%-d", k, l);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gsvc%-d_%-d", k, l);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      }
    }
    for (k = 0, TCSCptr = dataPtr->TCSCbus; TCSCptr != nullptr;
         TCSCptr = TCSCptr->Next) {
      for (k++, l = 1; l <= 7; l++) {
        i++;
        sprintf(str, "wtcsc%-d_%-d", k, l);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gtcsc%-d_%-d", k, l);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      }
    }
    for (k = 0, STATCOMptr = dataPtr->STATCOMbus; STATCOMptr != nullptr;
         STATCOMptr = STATCOMptr->Next) {
      for (k++, l = 1; l <= 7; l++) {
        i++;
        sprintf(str, "wstat%-d_%-d", k, l);
        fprintf(OutFile, "%4d %8s %-11.5g\n", i, str, x0[i - N]);
        sprintf(str, "gstat%-d_%-d", k, l);
        fprintf(OutFilep, "%4d %8s %-11.5g\n", i, str, dF[i]);
      }
    }
    /* END OF FACTS */

    fprintf(OutFile, "%4d %8s %-11.5g\n", Jac->n1, "l", lambda);
    fprintf(OutFilep, "%4d %8s %-11.5g\n", Jac->n1, "gl", dF[Jac->n1]);
  }
  fprintf(OutFile, "%4d %8s %-11.5g\n", 0, "0", 0.);
  fprintf(OutFilep, "%4d %8s %-11.5g\n", 0, "0", 0.);
  fclose(OutFile);
  fclose(OutFilep);
}

/* --------------- WriteQmatrix ---------------------- */
void WriteQmatrix(INDEX count, VALUETYPE *vec) {
  SparseMatrixElement *Jptr;
  ACbusData *ACptr;
  char Namebase[80], Name[80], MaxName[13];
  FILE *OutFile;
  INDEX i, j, I, J, m, MaxBus, MaxNum, Nvar;
  VALUETYPE val = 0, valMax = -0.1;
  bool flagTFbus = false;

  strcpy(Namebase, NameParameter('0'));
  if (NullName(Namebase))
    return;
  sprintf(Name, "%s%d.m", Namebase, count);
  OutFile = OpenOutput(Name);
  if (count == 1)
    TFbus = IntegerParameter('1', 0, 1, 9999);
  Nvar = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar + 7 * Nstatcom; /* FACTS */
  for (i = 1; i <= Nvar + 1; i++)
    for (Jptr = Jac->RowHead[i]; Jptr != nullptr;
         Jptr->Value = 0, Jptr = Jptr->RowNext)
      ;
  ACFunJac(Jac, nullptr, false, true, false);
  DCFunJac(Jac, false, true);
  SVCFunJac(Jac, false, true);     /*  FACTS  */
  TCSCFunJac(Jac, false, true);    /*  FACTS  */
  STATCOMFunJac(Jac, false, true); /*  FACTS  */
  fprintf(OutFile, "lambda(%d)=%-10.6g;\n", count, lambda + lambda_o);
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Full System Jacobian:\n", "%%");
  fprintf(OutFile, "N=%d;\nJ=zeros(N);\n", Nvar);
  for (i = 1; i <= Nvar + 1; i++) {
    for (Jptr = Jac->RowHead[i]; Jptr != nullptr; Jptr = Jptr->RowNext) {
      I = Jptr->Row;
      if (OldRow->p[I] != 0)
        I = OldRow->p[I];
      J = Jptr->Col;
      if (OldCol->p[J] != 0)
        J = OldCol->p[J];
      if (I <= Nvar && J <= Nvar)
        fprintf(OutFile, "J(%d,%d)=%-10.10g;\n", I, J, Jptr->Value);
    }
  }
  fprintf(OutFile, "J=sparse(J);\n");
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Voltage (load) Buses:\n", "%%");
  fprintf(OutFile, "Vars=[(1:N)' zeros(N,1) zeros(N,1) ones(N,1)];\n");
  fprintf(OutFile, "clear Buses Vbuses\n");
  for (m = 0, ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    i = ACvar[ACptr->N];
    j = i + 1;
    fprintf(OutFile, "Vars(%d,2)=%d; Vars(%d,2)=%d;\n", i, -ACptr->Num, j,
            ACptr->Num);
    if (!strpbrk(ACptr->Type, "S"))
      fprintf(OutFile, "Vars(%d,3)=1;\n", i);
    if ((ACptr->Cont != nullptr && (QRcont || !strpbrk(ACptr->Type, "G"))) ||
        (!Rcont && strpbrk(ACptr->Type, "T")) ||
        (!QRcont && strpbrk(ACptr->Type, "C"))) {
      i++;
      m++;
      fprintf(OutFile, "Vbuses(%d,1)=%d; Vbuses(%d,2)=%d;\n", m, i, m,
              ACptr->Num);
      fprintf(OutFile, "Vars(%d,3)=1; Vars(%d,4)=0;\n", i, i);
      if (count == 1) {
        val = fabs(vec[i]);
        if (TFbus == ACptr->Num && !flagTFbus) {
          flagTFbus = true;
          TFnum = i;
          strcpy(TFname, ACptr->Name);
        }
        if (val > valMax) {
          valMax = val;
          MaxBus = ACptr->Num;
          MaxNum = i;
          strcpy(MaxName, ACptr->Name);
        }
      }
    }
  }
  if (count == 1 && !flagTFbus) {
    TFnum = MaxNum;
    TFbus = MaxBus;
    strcpy(TFname, MaxName);
  }
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Bus for test functions and tangent vector:\n", "%%");
  fprintf(OutFile, "l=%d;\n", TFbus);
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Minimum singular value and real |e-value| for full J:\n",
          "%%");
  fprintf(OutFile, "[val,svecJ]=inviter(J'*J); svJ(%d)=sqrt(val);\n", count);
  fprintf(OutFile, "[val,evecJ]=inviter(J); evJ(%d)=abs(val);\n", count);
  fprintf(OutFile, "  %s Critical bus numbers and ranking of bus l:\n", "%%");
  fprintf(OutFile, "  [val,max_sv]=max(abs(svecJ));\n");
  fprintf(OutFile, "  crsvJ(%d,1)=Vars(max_sv,2);\n", count);
  fprintf(OutFile, "  [val,max_ev]=max(abs(evecJ));\n");
  fprintf(OutFile, "  crevJ(%d,1)=Vars(max_ev,2);\n", count);
  fprintf(OutFile, "  for i=1:N,\n");
  fprintf(OutFile, "    if (Vars(i,2)==l), lnumJ=i; break; end\n");
  fprintf(OutFile, "  end\n");
  fprintf(OutFile, "  crsvJ(%d,2)=rankbus(svecJ,lnumJ);\n", count);
  fprintf(OutFile, "  crevJ(%d,2)=rankbus(evecJ,lnumJ);\n", count);
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Compute J_PF matrix:\n", "%%");
  fprintf(OutFile, "P0=zeros(N);\n");
  fprintf(OutFile, "j=0;\n");
  fprintf(OutFile, "for i=1:N,\n");
  fprintf(OutFile, "  if(Vars(i,3)>0),\n");
  fprintf(OutFile, "    j=j+1;\n");
  fprintf(OutFile, "    Buses(j)=Vars(i,2);\n");
  fprintf(OutFile, "    if (Vars(i,2)==l), lnumPF=j; end\n");
  fprintf(OutFile, "    P0(Vars(i,1),j)=1;\n");
  fprintf(OutFile, "  end\n");
  fprintf(OutFile, "end\n");
  fprintf(OutFile, "M=j;\n");
  fprintf(OutFile, "for i=1:N,\n");
  fprintf(OutFile, "  if(Vars(i,3)==0), j=j+1; P0(Vars(i,1),j)=1; end\n");
  fprintf(OutFile, "end\n");
  fprintf(OutFile, "P0=sparse(P0); Jp=P0'*J*P0;\n");
  fprintf(OutFile, "J_PF=Jp(1:M,1:M);\n");
  fprintf(
      OutFile,
      "%s Minimum singular value and real |e-value| for standard P.F. J_PF:\n",
      "%%");
  fprintf(OutFile, "[val,svecPF]=inviter(J_PF'*J_PF); svPF(%d)=sqrt(val);\n",
          count);
  fprintf(OutFile, "[val,evecPF]=inviter(J_PF); evPF(%d)=abs(val);\n", count);
  fprintf(OutFile, "  %s Critical bus numbers and ranking of bus l:\n", "%%");
  fprintf(OutFile, "  [val,max_sv]=max(abs(svecPF));\n");
  fprintf(OutFile, "  crsvPF(%d,1)=Buses(max_sv);\n", count);
  fprintf(OutFile, "  [val,max_ev]=max(abs(evecPF));\n");
  fprintf(OutFile, "  crevPF(%d,1)=Buses(max_ev);\n", count);
  fprintf(OutFile, "  crsvPF(%d,2)=rankbus(svecPF,lnumPF);\n", count);
  fprintf(OutFile, "  crevPF(%d,2)=rankbus(evecPF,lnumPF);\n", count);
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Compute J_QV matrix:\n", "%%");
  fprintf(OutFile, "P1=zeros(N);\n");
  fprintf(OutFile, "Mp=%d; j=Mp;\n", m);
  fprintf(OutFile, "for i=1:Mp,\n");
  fprintf(OutFile, "  if (Vbuses(i,2)==l), lnumQV=i; end\n");
  fprintf(OutFile, "  P1(Vbuses(i,1),i)=1;\n");
  fprintf(OutFile, "end;\n");
  fprintf(OutFile, "for i=1:N,\n");
  fprintf(OutFile, "  if(Vars(i,4)>0), j=j+1; P1(Vars(i,1),j)=1; end\n");
  fprintf(OutFile, "end\n");
  fprintf(OutFile, "P1=sparse(P1); Jp=P1'*J*P1;\n");
  fprintf(OutFile, "A1=Jp(1:Mp,1:Mp);  B1=Jp(1:Mp,Mp+1:N); C1=Jp(Mp+1:N,1:Mp); "
                   "D1=Jp(Mp+1:N,Mp+1:N);\n");
  fprintf(OutFile, "J_QV=A1-B1*(D1\\C1);\n");
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Minimum singular value and real |e-value| for J_QV:\n",
          "%%");
  fprintf(OutFile, "[val,svecQV]=inviter(J_QV'*J_QV); svQV(%d)=sqrt(val);\n",
          count);
  fprintf(OutFile, "[val,evecQV]=inviter(J_QV,-0.001); evQV(%d)=abs(val);\n",
          count);
  fprintf(OutFile, "  %s Critical V bus numbers and ranking of bus l:\n", "%%");
  fprintf(OutFile, "  [val,max_sv]=max(abs(svecQV));\n");
  fprintf(OutFile, "  crsvQV(%d,1)=Vbuses(max_sv,2);\n", count);
  fprintf(OutFile, "  [val,max_ev]=max(abs(evecQV));\n");
  fprintf(OutFile, "  crevQV(%d,1)=Vbuses(max_ev,2);\n", count);
  fprintf(OutFile, "  crsvQV(%d,2)=rankbus(svecQV,lnumQV);\n", count);
  fprintf(OutFile, "  crevQV(%d,2)=rankbus(evecQV,lnumQV);\n", count);
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Reduced determinant detD_ll, for bus l=%d, '%12s':\n",
          "%%", TFbus, TFname);
  fprintf(OutFile, "P2=speye(N);\n");
  fprintf(OutFile, "v=P2(1,:); P2(1,:)=P2(%d,:); P2(%d,:)=v; \n", TFnum - 1,
          TFnum - 1);
  fprintf(OutFile, "v=P2(2,:); P2(2,:)=P2(%d,:); P2(%d,:)=v; \n", TFnum, TFnum);
  fprintf(OutFile, "Jp=P2'*J*P2;\n");
  fprintf(OutFile,
          "A2=Jp(1:2,1:2);  B2=Jp(1:2,3:N); C2=Jp(3:N,1:2); D2=Jp(3:N,3:N);\n");
  fprintf(OutFile, "detD_ll(%d)=det(A2-B2*(D2\\C2));\n", count);
  fprintf(OutFile, "\n");
  fprintf(OutFile, "%s Test Function t_ll, for bus l=%d, '%12s':\n", "%%",
          TFbus, TFname);
  fprintf(OutFile, "el=zeros(N,1); el(%d)=1; el=sparse(el);\n", TFnum);
  fprintf(OutFile, "J_ll=(speye(N)-el*el')*J+el*el';\n");
  fprintf(OutFile, "t_ll(%d)=el'*J*(J_ll\\el);\n", count);
  fclose(OutFile);
}
