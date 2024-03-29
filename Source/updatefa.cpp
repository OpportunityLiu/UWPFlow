/* Update FACTS variables. */

#include "constant.h"
#include "param.h"
#include "pflow.h"
#include "sparse.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void UpdateSVCvar(VALUETYPE cons, INDEX j);
void UpdateTCSCvar(VALUETYPE cons, INDEX j);
void UpdateSTATCOMvar(VALUETYPE cons, INDEX j);
bool ChangeSVCmode(void);
bool ChangeTCSCmode(void);
bool ChangeSTATCOMmode(void);

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX MaxIter, Nac, NacEl, NregPQ, NregV, Ndc, Nslack, Nvolt, Narea, NacVar, Bl, Nsvc, Ntcsc, NtcscVar, Nstatcom; /* FACTS */
extern INDEX *ACvar;
extern VALUETYPE *dx, *dF, tol, Tol, Sn, lambda, *x0;
extern VALUETYPE K1, K2, MaxdFi, alpha;
extern IntegerVector *NewRow, *OldRow, *NewCol, *OldCol, *RowPartition, *ColPartition;
extern IntegerVector *RowPer, *ColPer;
extern bool Acont, PQcont, QRcont, Rcont, PQlim, Tlim, Qlim, Vlim, flagH, flagPoC, flagL, flagR, flagBS;

void UpdateSVCvar(VALUETYPE cons, INDEX j) {
  SVCbusData *SVCptr;
  INDEX i;

  i = NacVar + 11 * Ndc / 2;
  for (SVCptr = dataPtr->SVCbus; SVCptr != nullptr; SVCptr = SVCptr->Next) {
    SVCptr->Qsvc = SVCptr->Qsvc + cons * dx[i + 1];
    SVCptr->Bv = SVCptr->Bv + cons * dx[i + 2];
    if (!strcmp(SVCptr->Cont, "AL")) {
      if (j == 0)
        SVCptr->val = SVCptr->alpha_svc;
      SVCptr->val = SVCptr->val + cons * dx[i + 3];
      if (SVCptr->val >= SVCptr->AlphaMax)
        SVCptr->alpha_svc = SVCptr->AlphaMax;
      else if (SVCptr->val <= SVCptr->AlphaMin)
        SVCptr->alpha_svc = SVCptr->AlphaMin;
      else
        SVCptr->alpha_svc = SVCptr->val;
    } else if (!strcmp(SVCptr->Cont, "MN")) {
      if (j == 0)
        SVCptr->val = SVCptr->Vvar;
      SVCptr->val = SVCptr->val + cons * dx[i + 3];
      if (SVCptr->val <= SVCptr->Vref)
        SVCptr->Vvar = SVCptr->Vref;
      else
        SVCptr->Vvar = SVCptr->val;
    } else {
      if (j == 0)
        SVCptr->val = SVCptr->Vvar;
      SVCptr->val = SVCptr->val + cons * dx[i + 3];
      if (SVCptr->val >= SVCptr->Vref)
        SVCptr->Vvar = SVCptr->Vref;
      else
        SVCptr->Vvar = SVCptr->val;
    }
    i = i + 3;
  }
}

bool ChangeSVCmode() {
  SVCbusData *SVCptr;
  bool flag = false;
  INDEX i;

  i = NacVar + 11 * Ndc / 2;
  for (SVCptr = dataPtr->SVCbus; SVCptr != nullptr; SVCptr = SVCptr->Next) {
    if (!strcmp(SVCptr->Cont, "AL") && SVCptr->alpha_svc >= SVCptr->AlphaMax) {
      strcpy(SVCptr->Cont, "MX");
      SVCptr->alpha_svc = SVCptr->AlphaMax;
      SVCptr->Vvar = SVCptr->Vref;
      if (flagH)
        x0[i + 3] = SVCptr->Vvar;
      fprintf(stderr, "***Warning: %s SVC firing angle is at its maximum limit.\n", SVCptr->Name);
      fprintf(stderr, "            The SVC will be treated as a fixed "
                      "capacitive reactance.\n");
      flag = true;
    } else if (!strcmp(SVCptr->Cont, "AL") && SVCptr->alpha_svc <= SVCptr->AlphaMin) {
      strcpy(SVCptr->Cont, "MN");
      SVCptr->alpha_svc = SVCptr->AlphaMin;
      SVCptr->Vvar = SVCptr->Vref;
      if (flagH)
        x0[i + 3] = SVCptr->Vvar;
      fprintf(stderr, "***Warning: %s SVC firing angle is at its minimum limit.\n", SVCptr->Name);
      fprintf(stderr, "            The SVC will be treated as a fixed "
                      "inductive reactance.\n");
      flag = true;
    } else if (!strcmp(SVCptr->Cont, "MX") && SVCptr->Vvar >= SVCptr->Vref) {
      strcpy(SVCptr->Cont, "AL");
      SVCptr->alpha_svc = SVCptr->AlphaMax;
      SVCptr->Vvar = SVCptr->Vref;
      if (flagH)
        x0[i + 3] = SVCptr->alpha_svc;
      fprintf(stderr, "***Warning: %s SVC firing angle is within limits.\n", SVCptr->Name);
      fprintf(stderr, "            The SVC voltage is now within controllable range.\n");
      flag = true;
    } else if (!strcmp(SVCptr->Cont, "MN") && SVCptr->Vvar <= SVCptr->Vref) {
      strcpy(SVCptr->Cont, "AL");
      SVCptr->alpha_svc = SVCptr->AlphaMin;
      SVCptr->Vvar = SVCptr->Vref;
      if (flagH)
        x0[i + 3] = SVCptr->alpha_svc;
      fprintf(stderr, "***Warning: %s SVC firing angle is within limits.\n", SVCptr->Name);
      fprintf(stderr, "            The SVC voltage is now within controllable range.\n");
      flag = true;
    }
    i = i + 3;
  }
  return (flag);
}

void UpdateTCSCvar(VALUETYPE cons, INDEX j) {
  TCSCbusData *TCSCptr;
  INDEX i;

  i = NacVar + 11 * Ndc / 2 + 3 * Nsvc;
  for (TCSCptr = dataPtr->TCSCbus; TCSCptr != nullptr; TCSCptr = TCSCptr->Next) {
    TCSCptr->Ptcsc = TCSCptr->Ptcsc + cons * dx[i + 1];
    TCSCptr->Qtcsck = TCSCptr->Qtcsck + cons * dx[i + 2];
    TCSCptr->Qtcscm = TCSCptr->Qtcscm + cons * dx[i + 3];
    TCSCptr->Be = TCSCptr->Be + cons * dx[i + 4];
    if (j == 0)
      TCSCptr->val = TCSCptr->alpha_tcsc;
    TCSCptr->val = TCSCptr->val + cons * dx[i + 5];
    if (TCSCptr->val >= TCSCptr->AlphaMax)
      TCSCptr->alpha_tcsc = TCSCptr->AlphaMax;
    else if (TCSCptr->val <= TCSCptr->AlphaMin)
      TCSCptr->alpha_tcsc = TCSCptr->AlphaMin;
    else
      TCSCptr->alpha_tcsc = TCSCptr->val;
    TCSCptr->Itcsc = TCSCptr->Itcsc + cons * dx[i + 6];
    TCSCptr->delta_t = TCSCptr->delta_t + cons * dx[i + 7];
    i = i + 7;
  }
}

bool ChangeTCSCmode() {
  TCSCbusData *TCSCptr;
  VALUETYPE alpha, Xc, Xl, Be;
  bool flag = false;

  for (TCSCptr = dataPtr->TCSCbus; TCSCptr != nullptr; TCSCptr = TCSCptr->Next) {
    alpha = TCSCptr->alpha_tcsc;
    Xc = TCSCptr->Xc;
    Xl = TCSCptr->Xl;
    if (!strpbrk(TCSCptr->Cont, "X") && TCSCptr->alpha_tcsc >= TCSCptr->AlphaMax) {
      strcpy(TCSCptr->Cont, "X");
      Be = 1.0 / Xc - (2.0 * PI - 2.0 * alpha + sin(2.0 * alpha)) / (PI * Xl);
      TCSCptr->Bset = Be;
      fprintf(stderr, "***Warning: %s TCSC firing angle is at its maximum limit.\n", TCSCptr->Name);
      fprintf(stderr, "            The TCSC will be treated as a fixed series "
                      "capacitor.\n");
      flag = true;
    } else if (!strpbrk(TCSCptr->Cont, "X") && TCSCptr->alpha_tcsc <= TCSCptr->AlphaMin) {
      strcpy(TCSCptr->Cont, "X");
      Be = 1.0 / Xc - (2.0 * PI - 2.0 * alpha + sin(2.0 * alpha)) / (PI * Xl);
      TCSCptr->Bset = Be;
      fprintf(stderr, "***Warning: %s TCSC firing angle is at its minimum limit.\n", TCSCptr->Name);
      fprintf(stderr, "            The TCSC will be treated as a fixed series "
                      "capacitor.\n");
      flag = true;
    }
  }
  return (flag);
}

void UpdateSTATCOMvar(VALUETYPE cons, INDEX j) {
  STATCOMbusData *STATCOMptr;
  VALUETYPE Q, theta, vals;
  INDEX i;

  i = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar;
  for (STATCOMptr = dataPtr->STATCOMbus; STATCOMptr != nullptr; STATCOMptr = STATCOMptr->Next) {
    Q = STATCOMptr->Q;
    if (!strcmp(STATCOMptr->Cont, "PW") || !strcmp(STATCOMptr->Cont, "AL")) {
      if (j == 0)
        STATCOMptr->val = STATCOMptr->I;
      STATCOMptr->val = STATCOMptr->val + cons * dx[i + 1];
      if (STATCOMptr->val <= 0)
        STATCOMptr->I = 0.0001;
      else if (STATCOMptr->val >= STATCOMptr->Imax && Q > 0)
        STATCOMptr->I = STATCOMptr->Imax;
      else if (STATCOMptr->val >= STATCOMptr->Imin && Q < 0)
        STATCOMptr->I = STATCOMptr->Imin;
      else
        STATCOMptr->I = STATCOMptr->val;
    } else if (!strcmp(STATCOMptr->Cont, "MX")) {
      if (j == 0)
        STATCOMptr->val = STATCOMptr->Vvar;
      STATCOMptr->val = STATCOMptr->val + cons * dx[i + 1];
      if (STATCOMptr->val <= STATCOMptr->Vref)
        STATCOMptr->Vvar = STATCOMptr->Vref;
      else
        STATCOMptr->Vvar = STATCOMptr->val;
    } else {
      if (j == 0)
        STATCOMptr->val = STATCOMptr->Vvar;
      STATCOMptr->val = STATCOMptr->val + cons * dx[i + 1];
      if (STATCOMptr->val >= STATCOMptr->Vref)
        STATCOMptr->Vvar = STATCOMptr->Vref;
      else
        STATCOMptr->Vvar = STATCOMptr->val;
    }
    theta = STATCOMptr->theta + cons * dx[i + 2];
    if (theta >= 0)
      vals = 1.00;
    else
      vals = -1.00;
    if (fabs(theta) > 2 * PI)
      theta = theta - vals * floor(fabs(theta) / (2 * PI)) * 2 * PI;
    if (fabs(theta) > PI)
      theta = theta - vals * 2 * PI;
    STATCOMptr->theta = theta;
    STATCOMptr->Vdc = STATCOMptr->Vdc + cons * dx[i + 3];
    STATCOMptr->k = STATCOMptr->k + cons * dx[i + 4];
    STATCOMptr->alpha = STATCOMptr->alpha + cons * dx[i + 5];
    STATCOMptr->P = STATCOMptr->P + cons * dx[i + 6];
    STATCOMptr->Q = STATCOMptr->Q + cons * dx[i + 7];
    i = i + 7;
  }
}

bool ChangeSTATCOMmode() {
  STATCOMbusData *STATCOMptr;
  bool flag = false;
  VALUETYPE Q;
  INDEX i;

  i = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar;
  for (STATCOMptr = dataPtr->STATCOMbus; STATCOMptr != nullptr; STATCOMptr = STATCOMptr->Next) {
    Q = STATCOMptr->Q;
    if ((!strcmp(STATCOMptr->Cont, "PW") || !strcmp(STATCOMptr->Cont, "AL")) && STATCOMptr->I >= STATCOMptr->Imax && Q > 0) {
      strcpy(STATCOMptr->Cont, "MX");
      STATCOMptr->I = STATCOMptr->Imax;
      STATCOMptr->Vvar = STATCOMptr->Vref;
      if (flagH)
        x0[i + 1] = STATCOMptr->Vvar;
      fprintf(stderr, "***Warning: %s STATCOM current is at its maximum limit.\n", STATCOMptr->Name);
      fprintf(stderr, "            The STATCOM will be treated as a fixed "
                      "inductive current source.\n");
      flag = true;
    } else if ((!strcmp(STATCOMptr->Cont, "PW") || !strcmp(STATCOMptr->Cont, "AL")) && STATCOMptr->I >= STATCOMptr->Imin && Q < 0) {
      strcpy(STATCOMptr->Cont, "MN");
      STATCOMptr->I = STATCOMptr->Imin;
      STATCOMptr->Vvar = STATCOMptr->Vref;
      if (flagH)
        x0[i + 1] = STATCOMptr->Vvar;
      fprintf(stderr, "***Warning: %s STATCOM current is at its minimum limit.\n", STATCOMptr->Name);
      fprintf(stderr, "            The STATCOM will be treated as a fixed "
                      "capacitive current source.\n");
      flag = true;
    } else if (!strcmp(STATCOMptr->Cont, "MN") && STATCOMptr->Vvar >= STATCOMptr->Vref) {
      strcpy(STATCOMptr->Cont, STATCOMptr->Cont1);
      STATCOMptr->I = STATCOMptr->Imax;
      STATCOMptr->Vvar = STATCOMptr->Vref;
      if (flagH)
        x0[i + 1] = STATCOMptr->I;
      fprintf(stderr, "***Warning: %s STATCOM current is within limits.\n", STATCOMptr->Name);
      fprintf(stderr, "            The STATCOM voltage is now within "
                      "controllable range.\n");
      flag = true;
    } else if (!strcmp(STATCOMptr->Cont, "MX") && STATCOMptr->Vvar <= STATCOMptr->Vref) {
      strcpy(STATCOMptr->Cont, STATCOMptr->Cont1);
      STATCOMptr->I = STATCOMptr->Imin;
      STATCOMptr->Vvar = STATCOMptr->Vref;
      if (flagH)
        x0[i + 1] = STATCOMptr->I;
      fprintf(stderr, "***Warning: %s STATCOM current is within limits.\n", STATCOMptr->Name);
      fprintf(stderr, "            The STATCOM voltage is now within "
                      "controllable range.\n");
      flag = true;
    }
    i = i + 7;
  }
  return (flag);
}
