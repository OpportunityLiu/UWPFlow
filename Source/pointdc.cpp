/* Point of Collapse:  DC Hessian. */

#include "pointl.h"

/* ------------------ DCFunHes ----------------------------- */
#ifdef ANSIPROTO
bool DCFunHes(bool flagF, bool flagJ)
#else
bool DCFunHes(flagF, flagJ)
    bool flagF,
    flagJ;
#endif
/* Construct the DC part of the PoC Jacobian. */
{
  INDEX i, j, k, l, kp, lp, m, n, N, DCvar[16];
  DCbusData *DCptrR, *DCptrI;
  ACbusData *BEptr;
  VALUETYPE Sa1, Sa2;
  VALUETYPE Vdr, Vr, ar, cosar, cosgr, Xcr, Sr, Pr, Dr, Id;
  VALUETYPE Vdi, Vi, ai, cosai, cosgi, Xci, Si, Pi, Di, Rd;
  bool dVr = false, dVi = false;

  i = NacVar;
  N = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar + 7 * Nstatcom; /*  FACTS  */
  for (DCptrR = dataPtr->DCbus; DCptrR != nullptr; DCptrR = DCptrR->Next)
  {
    DCptrI = DCptrR->To;
    if (!strcmp(DCptrR->Type, "R"))
    {
      Id = DCptrR->Id;
      Rd = DCptrR->Rd;
      Vdr = DCptrR->Vd;
      Vdi = DCptrI->Vd;
      k = ACvar[DCptrR->AC->N];
      l = ACvar[DCptrI->AC->N];
      kp = ACvar[DCptrR->AC->N] + 1;
      lp = ACvar[DCptrI->AC->N] + 1;
      if (DCptrR->AC->Cont != nullptr)
        dVr = true;
      if (DCptrI->AC->Cont != nullptr)
        dVi = true;
      Vr = DCptrR->AC->V;
      Vi = DCptrI->AC->V;
      ar = DCptrR->Tap * DCptrR->Ntrf;
      ai = DCptrI->Tap * DCptrI->Ntrf;
      Xcr = DCptrR->Xc;
      Xci = DCptrI->Xc;
      cosar = cos(DCptrR->Alfa);
      cosai = cos(DCptrI->Alfa);
      cosgr = cos(DCptrR->Gamma);
      cosgi = cos(DCptrI->Gamma);
      Pr = DCptrR->P;
      Pi = DCptrI->P;
      Sr = DCptrR->MVA;
      Si = DCptrI->MVA;
      Dr = Sr * Sr - Pr * Pr;
      Di = Si * Si - Pi * Pi;
      if (Dr <= 0 || Di <= 0)
        return (true);
      if (Acont && DCptrR->Meter != nullptr)
      {
        if (DCptrR->Meter == DCptrR)
        {
          Sa1 = -1;
          Sa2 = 1;
        }
        else
        {
          Sa1 = 1;
          Sa2 = -1;
        }
        BEptr = DCptrR->Area->Slack;
        if (!strpbrk(BEptr->Type, "S"))
          m = ACvar[BEptr->N] + 2;
        else
          m = 0;
        BEptr = DCptrI->Area->Slack;
        if (!strpbrk(BEptr->Type, "S"))
          n = ACvar[BEptr->N] + 2;
        else
          n = 0;
      }
      else
        m = n = 0;
      if (flagF)
      {
        j = i + N;
        if (dVr)
          dF[kp + N] = dF[kp + N] - K1 * ar * cosar * x0[i + 1] - K1 * ar * Id * x0[i + 2] - ar * (cosar + cosgr) * x0[i + 9];
        if (strcmp(DCptrR->Cont1, "VD") && strcmp(DCptrR->Cont2, "VD"))
          dF[++j] = x0[i + 1] + Id * x0[i + 3] + x0[i + 11];
        if (strcmp(DCptrR->Cont1, "AT") && strcmp(DCptrR->Cont2, "AT"))
          dF[++j] = -K1 * Vr * cosar * x0[i + 1] - K1 * Vr * Id * x0[i + 2] - Vr * (cosar + cosgr) * x0[i + 9];
        if (strcmp(DCptrR->Cont1, "AL") && strcmp(DCptrR->Cont2, "AL"))
          dF[++j] = -K1 * ar * Vr * x0[i + 1] - ar * Vr * x0[i + 9];
        if (strcmp(DCptrR->Cont1, "GA") && strcmp(DCptrR->Cont2, "GA"))
          dF[++j] = -ar * Vr * x0[i + 9];
        dF[++j] = x0[i + 2] + Sr / sqrt(Dr) * x0[i + 4];
        if (strcmp(DCptrR->Cont1, "PA") && strcmp(DCptrR->Cont2, "PA"))
        {
          dF[++j] = x0[i + 3] - Pr / sqrt(Dr) * x0[i + 4] + x0[k];
          if (m != 0 && DCptrR->Meter == DCptrR)
            dF[j] = dF[j] - Sa1 * x0[m];
          if (n != 0 && DCptrR->Meter == DCptrR)
            dF[j] = dF[j] - Sa2 * x0[n];
        }
        if (strcmp(DCptrR->Cont1, "QA") && strcmp(DCptrR->Cont2, "QA"))
          dF[++j] = x0[i + 4] + x0[kp];
        if (dVi)
          dF[lp + N] = dF[lp + N] - K1 * ai * cosgi * x0[i + 5] - K1 * ai * Id * x0[i + 6] - ai * (cosai + cosgi) * x0[i + 10];
        if (strcmp(DCptrI->Cont1, "VD") && strcmp(DCptrI->Cont2, "VD"))
          dF[++j] = x0[i + 5] - Id * x0[i + 7] - x0[i + 11];
        if (strcmp(DCptrI->Cont1, "AT") && strcmp(DCptrI->Cont2, "AT"))
          dF[++j] = -K1 * Vi * cosgi * x0[i + 5] - K1 * Vi * Id * x0[i + 6] - Vi * (cosai + cosgi) * x0[i + 10];
        if (strcmp(DCptrI->Cont1, "AL") && strcmp(DCptrI->Cont2, "AL"))
          dF[++j] = -ai * Vi * x0[i + 10];
        if (strcmp(DCptrI->Cont1, "GA") && strcmp(DCptrI->Cont2, "GA"))
          dF[++j] = -K1 * ai * Vi * x0[i + 5] - ai * Vi * x0[i + 10];
        dF[++j] = x0[i + 6] + Si / sqrt(Di) * x0[i + 8];
        if (strcmp(DCptrI->Cont1, "PA") && strcmp(DCptrI->Cont2, "PA"))
        {
          dF[++j] = x0[i + 7] - Pi / sqrt(Di) * x0[i + 8] + x0[l];
          if (m != 0 && DCptrI->Meter == DCptrI)
            dF[j] = dF[j] - Sa1 * x0[m];
          if (n != 0 && DCptrI->Meter == DCptrI)
            dF[j] = dF[j] - Sa2 * x0[n];
        }
        if (strcmp(DCptrI->Cont1, "QA") && strcmp(DCptrI->Cont2, "QA"))
          dF[++j] = x0[i + 8] + x0[lp];
        if (strcmp(DCptrR->Cont1, "ID") && strcmp(DCptrR->Cont2, "ID") &&
            strcmp(DCptrI->Cont1, "ID") && strcmp(DCptrI->Cont2, "ID"))
          dF[++j] = K2 * Xcr * x0[i + 1] - K1 * ar * Vr * x0[i + 2] + Vdr * x0[i + 3] + K2 * Xci * x0[i + 5] - K1 * ai * Vi * x0[i + 6] - Vdi * x0[i + 7] + sqrt(2.0) * (Xcr * x0[i + 9] + Xci * x0[i + 10]) - Rd * x0[i + 11];
      }
      if (flagJ)
      {
        j = i;
        if (strcmp(DCptrR->Cont1, "VD") && strcmp(DCptrR->Cont2, "VD"))
          DCvar[1] = ++j;
        else
          DCvar[1] = 0;
        if (strcmp(DCptrR->Cont1, "AT") && strcmp(DCptrR->Cont2, "AT"))
          DCvar[2] = ++j;
        else
          DCvar[2] = 0;
        if (strcmp(DCptrR->Cont1, "AL") && strcmp(DCptrR->Cont2, "AL"))
          DCvar[3] = ++j;
        else
          DCvar[3] = 0;
        if (strcmp(DCptrR->Cont1, "GA") && strcmp(DCptrR->Cont2, "GA"))
          DCvar[4] = ++j;
        else
          DCvar[4] = 0;
        DCvar[5] = ++j;
        if (strcmp(DCptrR->Cont1, "PA") && strcmp(DCptrR->Cont2, "PA"))
          DCvar[6] = ++j;
        else
          DCvar[6] = 0;
        if (strcmp(DCptrR->Cont1, "QA") && strcmp(DCptrR->Cont2, "QA"))
          DCvar[7] = ++j;
        else
          DCvar[7] = 0;
        if (strcmp(DCptrI->Cont1, "VD") && strcmp(DCptrI->Cont2, "VD"))
          DCvar[8] = ++j;
        else
          DCvar[8] = 0;
        if (strcmp(DCptrI->Cont1, "AT") && strcmp(DCptrI->Cont2, "AT"))
          DCvar[9] = ++j;
        else
          DCvar[9] = 0;
        if (strcmp(DCptrI->Cont1, "AL") && strcmp(DCptrI->Cont2, "AL"))
          DCvar[10] = ++j;
        else
          DCvar[10] = 0;
        if (strcmp(DCptrI->Cont1, "GA") && strcmp(DCptrI->Cont2, "GA"))
          DCvar[11] = ++j;
        else
          DCvar[11] = 0;
        DCvar[12] = ++j;
        if (strcmp(DCptrI->Cont1, "PA") && strcmp(DCptrI->Cont2, "PA"))
          DCvar[13] = ++j;
        else
          DCvar[13] = 0;
        if (strcmp(DCptrI->Cont1, "QA") && strcmp(DCptrI->Cont2, "QA"))
          DCvar[14] = ++j;
        else
          DCvar[14] = 0;
        if (strcmp(DCptrR->Cont1, "ID") && strcmp(DCptrR->Cont2, "ID") &&
            strcmp(DCptrI->Cont1, "ID") && strcmp(DCptrI->Cont2, "ID"))
          DCvar[15] = ++j;
        else
          DCvar[15] = 0;
        if (dVr)
        {
          if (DCvar[2])
            JacElement(Jac, kp + N, DCvar[2], -K1 * cosar * x0[i + 1] - K1 * Id * x0[i + 2] - (cosar + cosgr) * x0[i + 9]);
          if (DCvar[3])
            JacElement(Jac, kp + N, DCvar[3], -K1 * ar * x0[i + 1] - ar * x0[i + 9]);
          if (DCvar[4])
            JacElement(Jac, kp + N, DCvar[4], -ar * x0[i + 9]);
          if (DCvar[15])
            JacElement(Jac, kp + N, DCvar[15], -K1 * ar * x0[i + 2]);
        }
        else
        {
          if (DCvar[2])
            JacElement(Jac, kp + N, DCvar[2], 0.);
          if (DCvar[3])
            JacElement(Jac, kp + N, DCvar[3], 0.);
          if (DCvar[4])
            JacElement(Jac, kp + N, DCvar[4], 0.);
          if (DCvar[15])
            JacElement(Jac, kp + N, DCvar[15], 0.);
        }
        if (strcmp(DCptrR->Cont1, "VD") && strcmp(DCptrR->Cont2, "VD"))
        {
          if (DCvar[15])
            JacElement(Jac, DCvar[1] + N, DCvar[15], x0[i + 3]);
        }
        if (strcmp(DCptrR->Cont1, "AT") && strcmp(DCptrR->Cont2, "AT"))
        {
          if (dVr)
            JacElement(Jac, DCvar[2] + N, kp, -K1 * cosar * x0[1 + 1] - K1 * Id * x0[i + 2] - (cosar * cosgr) * x0[i + 9]);
          else
            JacElement(Jac, DCvar[2] + N, kp, 0.);
          if (DCvar[3])
            JacElement(Jac, DCvar[2] + N, DCvar[3], -K1 * Vr * x0[i + 1] - Vr * x0[i + 9]);
          if (DCvar[4])
            JacElement(Jac, DCvar[2] + N, DCvar[4], -Vr * x0[i + 9]);
          if (DCvar[15])
            JacElement(Jac, DCvar[2] + N, DCvar[15], -K1 * Vr * x0[i + 2]);
        }
        if (strcmp(DCptrR->Cont1, "AL") && strcmp(DCptrR->Cont2, "AL"))
        {
          if (dVr)
            JacElement(Jac, DCvar[3] + N, kp, -K1 * ar * x0[1 + 1] - ar * x0[i + 9]);
          else
            JacElement(Jac, DCvar[3] + N, kp, 0.);
          if (DCvar[2])
            JacElement(Jac, DCvar[3] + N, DCvar[2], -K1 * Vr * x0[1 + 1] - Vr * x0[i + 9]);
        }
        if (strcmp(DCptrR->Cont1, "GA") && strcmp(DCptrR->Cont2, "GA"))
        {
          if (dVr)
            JacElement(Jac, DCvar[4] + N, kp, -ar * x0[i + 9]);
          else
            JacElement(Jac, DCvar[4] + N, kp, 0.);
          if (DCvar[2])
            JacElement(Jac, DCvar[4] + N, DCvar[2], -Vr * x0[i + 9]);
        }
        JacElement(Jac, DCvar[5] + N, DCvar[5], -Pr * Pr * x0[i + 4] / (Dr * sqrt(Dr)));
        if (DCvar[6])
          JacElement(Jac, DCvar[5] + N, DCvar[6], Sr * Pr * x0[i + 4] / (Dr * sqrt(Dr)));
        if (strcmp(DCptrR->Cont1, "PA") && strcmp(DCptrR->Cont2, "PA"))
        {
          JacElement(Jac, DCvar[6] + N, DCvar[5], Sr * Pr * x0[i + 4] / (Dr * sqrt(Dr)));
          JacElement(Jac, DCvar[6] + N, DCvar[6], -Sr * Sr * x0[i + 4] / (Dr * sqrt(Dr)));
        }
        if (dVi)
        {
          if (DCvar[9])
            JacElement(Jac, lp + N, DCvar[9], -K1 * cosgi * x0[i + 5] - K1 * Id * x0[i + 6] - (cosai + cosgi) * x0[i + 10]);
          if (DCvar[10])
            JacElement(Jac, lp + N, DCvar[10], -ai * x0[i + 10]);
          if (DCvar[11])
            JacElement(Jac, lp + N, DCvar[11], -K1 * ai * x0[i + 5] - ai * x0[i + 10]);
          if (DCvar[15])
            JacElement(Jac, lp + N, DCvar[15], -K1 * ai * x0[i + 6]);
        }
        else
        {
          if (DCvar[9])
            JacElement(Jac, lp + N, DCvar[9], 0.);
          if (DCvar[10])
            JacElement(Jac, lp + N, DCvar[10], 0.);
          if (DCvar[11])
            JacElement(Jac, lp + N, DCvar[11], 0.);
          if (DCvar[15])
            JacElement(Jac, lp + N, DCvar[15], 0.);
        }
        if (strcmp(DCptrI->Cont1, "VD") && strcmp(DCptrI->Cont2, "VD"))
        {
          if (DCvar[15])
            JacElement(Jac, DCvar[8] + N, DCvar[15], -x0[i + 7]);
        }
        if (strcmp(DCptrI->Cont1, "AT") && strcmp(DCptrI->Cont2, "AT"))
        {
          if (dVi)
            JacElement(Jac, DCvar[9] + N, lp, -K1 * cosgi * x0[i + 5] - K1 * Id * x0[i + 6] - (cosai + cosgi) * x0[i + 10]);
          else
            JacElement(Jac, DCvar[9] + N, lp, 0.);
          if (DCvar[10])
            JacElement(Jac, DCvar[9] + N, DCvar[10], -Vi * x0[i + 10]);
          if (DCvar[11])
            JacElement(Jac, DCvar[9] + N, DCvar[11], -K1 * Vi * x0[i + 5] - Vi * x0[i + 10]);
          if (DCvar[15])
            JacElement(Jac, DCvar[9] + N, DCvar[15], -K1 * Vi * x0[i + 6]);
        }
        if (strcmp(DCptrI->Cont1, "AL") && strcmp(DCptrI->Cont2, "AL"))
        {
          if (dVi)
            JacElement(Jac, DCvar[10] + N, lp, -ai * x0[i + 10]);
          else
            JacElement(Jac, DCvar[10] + N, lp, 0.);
          if (DCvar[9])
            JacElement(Jac, DCvar[10] + N, DCvar[9], -Vi * x0[i + 10]);
        }
        if (strcmp(DCptrI->Cont1, "GA") && strcmp(DCptrI->Cont2, "GA"))
        {
          if (dVi)
            JacElement(Jac, DCvar[11] + N, lp, -K1 * ai * x0[i + 5] - ai * x0[i + 10]);
          else
            JacElement(Jac, DCvar[11] + N, lp, 0.);
          if (DCvar[9])
            JacElement(Jac, DCvar[11] + N, DCvar[9], -K1 * Vi * x0[i + 5] - Vi * x0[i + 10]);
        }
        JacElement(Jac, DCvar[12] + N, DCvar[12], -Pi * Pi * x0[i + 8] / (Di * sqrt(Di)));
        if (DCvar[13])
          JacElement(Jac, DCvar[12] + N, DCvar[13], Si * Pi * x0[i + 8] / (Di * sqrt(Di)));
        if (strcmp(DCptrI->Cont1, "PA") && strcmp(DCptrI->Cont2, "PA"))
        {
          JacElement(Jac, DCvar[13] + N, DCvar[12], Si * Pi * x0[i + 8] / (Di * sqrt(Di)));
          JacElement(Jac, DCvar[13] + N, DCvar[13], -Si * Si * x0[i + 8] / (Di * sqrt(Di)));
        }
        if (strcmp(DCptrR->Cont1, "ID") && strcmp(DCptrR->Cont2, "ID") &&
            strcmp(DCptrI->Cont1, "ID") && strcmp(DCptrI->Cont2, "ID"))
        {
          if (dVr)
            JacElement(Jac, DCvar[15] + N, kp, -K1 * ar * x0[i + 2]);
          else
            JacElement(Jac, DCvar[15] + N, kp, 0.);
          if (DCvar[1])
            JacElement(Jac, DCvar[15] + N, DCvar[1], x0[i + 3]);
          if (DCvar[2])
            JacElement(Jac, DCvar[15] + N, DCvar[2], -K1 * Vr * x0[i + 2]);
          if (dVi)
            JacElement(Jac, DCvar[15] + N, lp, -K1 * ai * x0[i + 6]);
          else
            JacElement(Jac, DCvar[15] + N, lp, 0.);
          if (DCvar[8])
            JacElement(Jac, DCvar[15] + N, DCvar[8], -x0[i + 7]);
          if (DCvar[9])
            JacElement(Jac, DCvar[15] + N, DCvar[9], -K1 * Vi * x0[i + 6]);
        }
      }
      i = i + 11;
    }
  }
  return (false);
}
