#include "constant.h"
#include "param.h"
#include "pflow.h"
#include "sparse.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void SVCinit(void);
void TCSCinit(void);
void STATCOMinit(void);

extern Data *dataPtr;
extern VALUETYPE Sn, K3;

/* --------------------------- SVCinit ------------------------------------ */
/* Initialize the SVC variables and transform them to the system p.u. base. */
void SVCinit(void) {
  SVCbusData *SVCptr;
  INDEX i, j;
  VALUETYPE Vc, V, Vref, Xsl, Bv, Xl, Xc, Be, Ssvc, I, alpha, x, dx, dxp, k, f,
      df;

  for (SVCptr = dataPtr->SVCbus; SVCptr != nullptr; SVCptr = SVCptr->Next) {
    Vc = SVCptr->Ctrl->V;
    V = SVCptr->From->V;
    Vref = SVCptr->Vref;
    Ssvc = SVCptr->SVC_base;
    Xsl = SVCptr->slope / 100.0 * Sn / Ssvc;
    if (Xsl != 0)
      I = (Vc - Vref) / Xsl;
    else
      I = 0;
    Be = -I / V;
    Xl = SVCptr->Xl * Sn / Ssvc;
    Xc = SVCptr->Xc * Sn / Ssvc;
    Bv = 1.0 / Xc - Be;
    alpha = SVCptr->alpha_svc = SVCptr->alpha_svc * K3;
    SVCptr->AlphaMin = SVCptr->AlphaMin * K3;
    SVCptr->AlphaMax = SVCptr->AlphaMax * K3;
    if (alpha == 0.) {
      x = PI;
      for (i = 0; i <= 10; i = i + 1) {
        f = sin(x) - x + 2.0 * PI - PI * Xl * Bv;
        if (fabs(f) <= 0.001)
          break;
        df = cos(x) - 1.0;
        if (df != 0.) {
          dx = -f / df;
          k = 1 / 2.0;
          if (i == 0)
            dxp = dx;
          else {
            if (fabs(dx) > fabs(dxp))
              for (j = 0; j <= 10; j = j + 1) {
                dx = k * dx;
                if (fabs(dx) < fabs(dxp))
                  break;
                k = k / 2.0;
              }
          }
          x = x + dx;
          dxp = dx;
          if (k < 0.001)
            break;
        } else
          break;
      }
      alpha = x / 2.0;
    }
    if (alpha <= SVCptr->AlphaMin) {
      alpha = SVCptr->alpha_svc = SVCptr->AlphaMin;
      strcpy(SVCptr->Cont, "MN");
      fprintf(
          stderr,
          "***Warning: SVC %s initial firing angle is at its minimum limit.\n",
          SVCptr->Name);
      fprintf(stderr, "            The SVC will be treated as a fixed "
                      "inductive reactance.\n");
    } else if (alpha >= SVCptr->AlphaMax) {
      alpha = SVCptr->alpha_svc = SVCptr->AlphaMax;
      strcpy(SVCptr->Cont, "MX");
      fprintf(
          stderr,
          "***Warning: SVC %s initial firing angle is at its maximum limit.\n",
          SVCptr->Name);
      fprintf(stderr, "            The SVC will be treated as a fixed "
                      "capacitive reactance.\n");
    }
    Bv = (2.0 * (PI - alpha) + sin(2.0 * alpha)) / (PI * Xl);
    SVCptr->alpha_svc = alpha;
    SVCptr->Bv = Bv;
    SVCptr->Qsvc = V * V / Xc - V * V * Bv;
    SVCptr->slope = Xsl;
    SVCptr->Xc = Xc;
    SVCptr->Xl = Xl;
    SVCptr->Vvar = SVCptr->Vref;
  }
}

/* --------------------------- TCSCinit ------------------------------------ */
void TCSCinit(void) {
  TCSCbusData *TCSCptr;
  VALUETYPE Stcsc, Vn, Max, Ctrl, Kf;
  VALUETYPE Vk, thk, Vm, thm, Xc, Xl, Ptcsc, Qtcsck, Qtcscm, Be, alpha, Itcsc,
      delta_t;

  for (TCSCptr = dataPtr->TCSCbus; TCSCptr != nullptr;
       TCSCptr = TCSCptr->Next) {
    /* --------------- Initialize the TCSC system variables
     * ---------------------- */
    if (TCSCptr->To->Ang == 0.0)
      TCSCptr->To->Ang = 0.01;
    Vk = TCSCptr->From->V;
    thk = TCSCptr->From->Ang;
    Vm = TCSCptr->To->V;
    thm = TCSCptr->To->Ang;
    /* ----------------- Per unit conversion
     * ------------------------------------ */
    Stcsc = TCSCptr->TCSC_base;
    if (Stcsc == 0)
      Stcsc = TCSCptr->TCSC_base = Sn;
    Vn = TCSCptr->From->KV;
    TCSCptr->Xc = TCSCptr->Xc * Sn / Stcsc;
    TCSCptr->Xl = TCSCptr->Xl * Sn / Stcsc;
    TCSCptr->AlphaMin = TCSCptr->AlphaMin * K3;
    TCSCptr->AlphaMax = TCSCptr->AlphaMax * K3;
    Xc = TCSCptr->Xc;
    Xl = TCSCptr->Xl;
    Kf = sqrt(Xc / Xl);
    Max = TCSCptr->Max;
    Ctrl = TCSCptr->Control;
    /* --------------- Initialize the TCSC system variables
     * ---------------------- */
    alpha = TCSCptr->alpha_tcsc;
    if (alpha == 0)
      alpha = 2.8;
    else if (alpha < TCSCptr->AlphaMin) {
      TCSCptr->alpha_tcsc = alpha = 2.8;
      fprintf(stderr,
              "***Warning: %s TCSC firing angle is below its minimum limit.\n",
              TCSCptr->Name);
      fprintf(stderr,
              "            The TCSC will be initialized with a flat start.\n");
    } else if (alpha > TCSCptr->AlphaMax) {
      TCSCptr->alpha_tcsc = alpha = 2.8;
      fprintf(stderr,
              "***Warning: %s TCSC firing angle is above its maximum limit.\n",
              TCSCptr->Name);
      fprintf(stderr,
              "            The TCSC will be initialized with a flat start.\n");
    }
    if (!strcmp(TCSCptr->Cont, "X")) {
      Be = 100.0 / (Ctrl * Max * Xc);
      TCSCptr->Bset = Be;
      TCSCptr->Control = Be;
      /* } else Be=1.0/Xc-(2.0*PI-2.0*alpha+sin(2.0*alpha))/(PI*Xl); */
    } else {
      Be = PI * cos(Kf * (-PI + alpha)) * (pow(Kf, 4.0) - 2.0 * Kf * Kf + 1.0) /
           Xc /
           (-PI * pow(Kf, 4.0) * cos(Kf * (-PI + alpha)) +
            PI * cos(Kf * (-PI + alpha)) +
            2.0 * pow(Kf, 4.0) * alpha * cos(Kf * (-PI + alpha)) -
            2.0 * alpha * Kf * Kf * cos(Kf * (-PI + alpha)) +
            pow(Kf, 4.0) * sin(-2.0 * PI + 2.0 * alpha) *
                cos(Kf * (-PI + alpha)) -
            sin(-2.0 * PI + 2.0 * alpha) * Kf * Kf * cos(Kf * (-PI + alpha)) -
            4.0 * Kf * Kf * Kf * pow(cos(-PI + alpha), 2.0) *
                sin(Kf * (-PI + alpha)) +
            4.0 * Kf * Kf * cos(-PI + alpha) * sin(-PI + alpha) *
                cos(Kf * (-PI + alpha)));
    }
    if (!strcmp(TCSCptr->Cont, "P")) {
      Ptcsc = Ctrl / Sn;
      TCSCptr->Control = Ptcsc;
    } else
      Ptcsc = -Vk * Vm * Be * sin(thk - thm);
    Qtcsck = -Vk * Vk * Be + Vk * Vm * Be * cos(thk - thm);
    Qtcscm = -Vm * Vm * Be + Vk * Vm * Be * cos(thk - thm);
    if (!strcmp(TCSCptr->Cont, "I")) {
      Itcsc = Ctrl * (sqrt(3.0) * Vn) / Sn;
      TCSCptr->Control = Itcsc;
    } else {
      if (Ptcsc >= 0)
        Itcsc = (sqrt(Ptcsc * Ptcsc + Qtcsck * Qtcsck)) / Vk;
      else
        Itcsc = -(sqrt(Ptcsc * Ptcsc + Qtcsck * Qtcsck)) / Vk;
    }
    if (!strcmp(TCSCptr->Cont, "D")) {
      delta_t = Ctrl * K3;
      TCSCptr->Control = delta_t;
    } else
      delta_t = thk - thm;
    TCSCptr->Ptcsc = Ptcsc;
    TCSCptr->Qtcsck = Qtcsck;
    TCSCptr->Qtcscm = Qtcscm;
    TCSCptr->Be = Be;
    TCSCptr->alpha_tcsc = alpha;
    TCSCptr->Itcsc = Itcsc;
    TCSCptr->delta_t = delta_t;
  }
}

/* ----------------------------- STATCOMinit ----------------------------------
 */
/* Initialize the STATCOM variables assuming Vc and V are known and neglecting
 */
/* losses, i.e., R and Rc are assumed to be zero to simplify the equations, and
 */
/* transform variables to the system p.u. base. */
void STATCOMinit(void) {
  STATCOMbusData *STATCOMptr;
  VALUETYPE Sb, Vco, Vo, deltao, Io, thetao, Imin, Imax, R, G, B, Gc, Xsl, Ko,
      Vdco, Vref, Cref, Qo;

  for (STATCOMptr = dataPtr->STATCOMbus; STATCOMptr != nullptr;
       STATCOMptr = STATCOMptr->Next) {
    Sb = STATCOMptr->MVA;
    Vo = STATCOMptr->From->V;
    deltao = STATCOMptr->From->Ang;
    Xsl = STATCOMptr->slope / 100.0 * Sn / Sb;
    R = STATCOMptr->R * Sn / Sb;
    G = STATCOMptr->G * Sb / Sn;
    B = STATCOMptr->B * Sb / Sn;
    Gc = STATCOMptr->Gc * Sb / Sn;
    Imin = STATCOMptr->Imin * Sb / Sn;
    Imax = STATCOMptr->Imax * Sb / Sn;
    Vref = STATCOMptr->Vref;
    Vco = STATCOMptr->Ctrl->V;
    Cref = STATCOMptr->Contref;
    if (Xsl != 0)
      Io = (Vco - Vref) / Xsl;
    else
      Io = 0;
    if (Io < -Imin) {
      Io = -Imin;
      strcpy(STATCOMptr->Cont, "MN");
      fprintf(stderr,
              "***Warning: STATCOM %s has been initialized at its minimum "
              "current limit.\n",
              STATCOMptr->Name);
      fprintf(stderr, "            It will be treated as a fixed capacitive "
                      "current source. \n");
    } else if (Io > Imax) {
      Io = Imax;
      strcpy(STATCOMptr->Cont, "MX");
      fprintf(stderr,
              "***Warning: STATCOM %s has been initialized at its maximum "
              "current limit.\n",
              STATCOMptr->Name);
      fprintf(stderr, "            It will be treated as a fixed inductive "
                      "current source. \n");
    }
    if (Io <= 0)
      thetao = STATCOMptr->theta = deltao + PI / 2.0; /* capacitive */
    else
      thetao = STATCOMptr->theta = deltao - PI / 2.0; /* inductive  */
    Io = fabs(Io);
    if (Io == 0)
      Io = 0.0001;
    STATCOMptr->I = Io;
    STATCOMptr->Imin = Imin;
    STATCOMptr->Imax = Imax;
    Qo = STATCOMptr->Q = Vo * Io * sin(deltao - thetao);
    STATCOMptr->alpha = deltao;
    if (!strcmp(STATCOMptr->Cont1, "PW")) {
      Vdco = STATCOMptr->Vdc = Cref;
      Ko = STATCOMptr->k = 1.0 / Vdco * (Qo / (Vo * B) + Vo);
    } else {
      Ko = STATCOMptr->k = Cref;
      Vdco = STATCOMptr->Vdc = 1.0 / Ko * (Qo / (Vo * B) + Vo);
    }
    STATCOMptr->P = Gc * Vdco * Vdco + R * Io * Io;
    STATCOMptr->slope = Xsl;
    STATCOMptr->R = R;
    STATCOMptr->G = G;
    STATCOMptr->Gc = Gc;
    STATCOMptr->B = B;
    STATCOMptr->Vvar = STATCOMptr->Vref;
  }
}
