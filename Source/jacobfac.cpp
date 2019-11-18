/* FACTS mismatch vectors and Jacobian */

#include "jacob.h"

void SVCFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ);
void TCSCFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ);
void STATCOMFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ);

/*  Construct the SVC part of the mismatch vector and Jacobian  */
void SVCFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ)
{
  INDEX i, k, l;
  SVCbusData *SVCptr;
  VALUETYPE Vref, Vk, Vl, Xl, Xc, Bv, Xsl, Qsvc, alpha;
  bool flag1;

  i = NacVar + 11 * Ndc / 2;
  for (SVCptr = dataPtr->SVCbus; SVCptr != nullptr; SVCptr = SVCptr->Next)
  {
    k = ACvar[SVCptr->From->N];
    Vk = SVCptr->From->V;
    l = ACvar[SVCptr->Ctrl->N];
    Vl = SVCptr->Ctrl->V;
    if (!strcmp(SVCptr->Cont, "AL"))
      flag1 = false;
    else
      flag1 = true;
    if (!flag1)
      Vref = SVCptr->Vref;
    else
      Vref = SVCptr->Vvar;
    Xl = SVCptr->Xl;
    Xc = SVCptr->Xc;
    Xsl = SVCptr->slope;
    Bv = SVCptr->Bv;
    Qsvc = SVCptr->Qsvc;
    alpha = SVCptr->alpha_svc;
    if (flagF)
    {
      dF[i + 1] = Vl - Vref - Xsl * Vk * (Bv - 1.0 / Xc);
      dF[i + 2] = Vk * Vk * (Bv - 1.0 / Xc) + Qsvc;
      dF[i + 3] = -PI * Xl * Bv + 2.0 * (PI - alpha) + sin(2.0 * alpha);
      dF[k + 1] = dF[k + 1] + Qsvc;
    }
    if (flagJ)
    {
      JacElement(Mptr, k + 1, i + 1, 1.0);
      JacElement(Mptr, i + 1, k + 1, -Xsl * (Bv - 1.0 / Xc));
      JacElement(Mptr, i + 1, l + 1, 1.0);
      JacElement(Mptr, i + 1, i + 2, -Xsl * Vk);
      JacElement(Mptr, i + 2, k + 1, -2.0 * Vk * (1.0 / Xc - Bv));
      JacElement(Mptr, i + 2, i + 1, 1.0);
      JacElement(Mptr, i + 2, i + 2, Vk * Vk);
      JacElement(Mptr, i + 3, i + 2, -PI * Xl);
      if (!flag1)
      {
        JacElement(Mptr, i + 1, i + 3, 0.0);
        JacElement(Mptr, i + 3, i + 3, -2.0 * (1.0 - cos(2.0 * alpha)));
      }
      else
      {
        JacElement(Mptr, i + 1, i + 3, -1.0);
        JacElement(Mptr, i + 3, i + 3, 0.0);
      }
      if (flagH)
      {
        if (strpbrk(SVCptr->From->Type, "L"))
        {
          JacElement(Mptr, i + 1, Mptr->n1, -Xsl * (Bv - 1.0 / Xc));
          JacElement(Mptr, i + 2, Mptr->n1, 2.0 * Vk * (Bv - 1.0 / Xc));
        }
        if (strpbrk(SVCptr->Ctrl->Type, "L"))
        {
          JacElement(Mptr, i + 1, Mptr->n1, 1.0);
        }
      }
    }
    i = i + 3;
  }
}

/*  Construct the TCSC part of the mismatch vector and Jacobian  */
void TCSCFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ)
{
  INDEX i, k, kp, m, mp;
  TCSCbusData *TCSCptr;
  VALUETYPE Vk, Vm, thk, thm, Xl, Xc, Ptcsc, Qtcsck, Qtcscm, Be, alpha, Itcsc, delta_t;
  VALUETYPE Bset, Pset, Iset, delta_set, Ctrl, Kf;
  VALUETYPE s1, s2, s3, s4, s5, s6, s7, s8;
  bool dVk = false, dVm = false;

  i = NacVar + 11 * Ndc / 2 + 3 * Nsvc;
  for (TCSCptr = dataPtr->TCSCbus; TCSCptr != nullptr; TCSCptr = TCSCptr->Next)
  {
    k = ACvar[TCSCptr->From->N];
    kp = k + 1;
    Vk = TCSCptr->From->V;
    thk = TCSCptr->From->Ang;
    m = ACvar[TCSCptr->To->N];
    mp = m + 1;
    Vm = TCSCptr->To->V;
    thm = TCSCptr->To->Ang;
    Xc = TCSCptr->Xc;
    Xl = TCSCptr->Xl;
    Kf = sqrt(Xc / Xl);
    Ctrl = TCSCptr->Control;
    if (!strcmp(TCSCptr->Cont, "X"))
      Bset = TCSCptr->Bset;
    else if (!strcmp(TCSCptr->Cont, "P"))
      Pset = Ctrl;
    else if (!strcmp(TCSCptr->Cont, "I"))
      Iset = Ctrl;
    else if (!strcmp(TCSCptr->Cont, "D"))
      delta_set = Ctrl;
    Ptcsc = TCSCptr->Ptcsc;
    Qtcsck = TCSCptr->Qtcsck;
    Qtcscm = TCSCptr->Qtcscm;
    Be = TCSCptr->Be;
    alpha = TCSCptr->alpha_tcsc;
    Itcsc = TCSCptr->Itcsc;
    delta_t = TCSCptr->delta_t;
    if (TCSCptr->From->Cont != nullptr)
      dVk = true;
    else
      dVk = false;
    if (TCSCptr->To->Cont != nullptr)
      dVm = true;
    else
      dVm = false;
    if (flagF)
    {
      dF[i + 1] = -Vk * Vm * Be * sin(thk - thm) - Ptcsc;
      dF[i + 2] = -Vk * Vk * Be + Vk * Vm * Be * cos(thk - thm) - Qtcsck;
      dF[i + 3] = -Vm * Vm * Be + Vk * Vm * Be * cos(thk - thm) - Qtcscm;
      /*
     dF[i+4]=1.0/Xc-(2.0*(PI-alpha)+sin(2.0*alpha))/(PI*Xl)-Be;
     */
      dF[i + 4] = PI * cos(Kf * (-PI + alpha)) * (pow(Kf, 4.0) - 2.0 * Kf * Kf + 1.0) / Xc / (-PI * pow(Kf, 4.0) * cos(Kf * (-PI + alpha)) + PI * cos(Kf * (-PI + alpha)) + 2.0 * pow(Kf, 4.0) * alpha * cos(Kf * (-PI + alpha)) - 2.0 * alpha * Kf * Kf * cos(Kf * (-PI + alpha)) + pow(Kf, 4.0) * sin(-2.0 * PI + 2.0 * alpha) * cos(Kf * (-PI + alpha)) - sin(-2.0 * PI + 2.0 * alpha) * Kf * Kf * cos(Kf * (-PI + alpha)) - 4.0 * Kf * Kf * Kf * pow(cos(-PI + alpha), 2.0) * sin(Kf * (-PI + alpha)) + 4.0 * Kf * Kf * cos(-PI + alpha) * sin(-PI + alpha) * cos(Kf * (-PI + alpha))) - Be;
      if (Itcsc >= 0)
        dF[i + 5] = sqrt(Ptcsc * Ptcsc + Qtcsck * Qtcsck) - Itcsc * Vk;
      else
        dF[i + 5] = -sqrt(Ptcsc * Ptcsc + Qtcsck * Qtcsck) - Itcsc * Vk;
      dF[i + 6] = thk - thm - delta_t;
      if (!strcmp(TCSCptr->Cont, "X"))
        dF[i + 7] = Bset - Be;
      else if (!strcmp(TCSCptr->Cont, "P"))
        dF[i + 7] = Pset - Ptcsc;
      else if (!strcmp(TCSCptr->Cont, "I"))
        dF[i + 7] = Iset - Itcsc;
      else if (!strcmp(TCSCptr->Cont, "D"))
        dF[i + 7] = delta_set - delta_t;
      dF[k] = dF[k] - Ptcsc;
      dF[k + 1] = dF[k + 1] - Qtcsck;
      dF[m] = dF[m] + Ptcsc;
      dF[m + 1] = dF[m + 1] - Qtcscm;
    }
    if (flagJ)
    {
      JacElement(Mptr, k, i + 1, -1.0);
      JacElement(Mptr, k + 1, i + 2, -1.0);
      JacElement(Mptr, m, i + 1, 1.0);
      JacElement(Mptr, m + 1, i + 3, -1.0);
      if (!strpbrk(TCSCptr->From->Type, "S"))
        JacElement(Mptr, i + 1, k, -Vk * Vm * Be * cos(thk - thm));
      if (dVk)
      {
        JacElement(Mptr, i + 1, kp, -Vm * Be * sin(thk - thm));
        JacElement(Mptr, i + 2, kp, -2.0 * Vk * Be + Vm * Be * cos(thk - thm));
        JacElement(Mptr, i + 3, kp, Vm * Be * cos(thk - thm));
        JacElement(Mptr, i + 5, kp, -Itcsc);
      }
      else
      {
        JacElement(Mptr, i + 1, kp, 0.0);
        JacElement(Mptr, i + 2, kp, 0.0);
        JacElement(Mptr, i + 3, kp, 0.0);
        JacElement(Mptr, i + 5, kp, 0.0);
      }
      if (!strpbrk(TCSCptr->To->Type, "S"))
        JacElement(Mptr, i + 1, m, Vk * Vm * Be * cos(thk - thm));
      if (dVm)
      {
        JacElement(Mptr, i + 1, mp, -Vk * Be * sin(thk - thm));
        JacElement(Mptr, i + 2, mp, Vk * Be * cos(thk - thm));
        JacElement(Mptr, i + 3, mp, -2 * Vm * Be + Vk * Be * cos(thk - thm));
      }
      else
      {
        JacElement(Mptr, i + 1, mp, 0.0);
        JacElement(Mptr, i + 2, mp, 0.0);
        JacElement(Mptr, i + 3, mp, 0.0);
      }
      JacElement(Mptr, i + 1, i + 1, -1.0);
      JacElement(Mptr, i + 1, i + 4, -Vk * Vm * sin(thk - thm));
      if (!strpbrk(TCSCptr->From->Type, "S"))
        JacElement(Mptr, i + 2, k, -Vk * Vm * Be * sin(thk - thm));
      if (!strpbrk(TCSCptr->To->Type, "S"))
        JacElement(Mptr, i + 2, m, Vk * Vm * Be * sin(thk - thm));
      JacElement(Mptr, i + 2, i + 2, -1.0);
      JacElement(Mptr, i + 2, i + 4, -Vk * Vk + Vk * Vm * cos(thk - thm));
      if (!strpbrk(TCSCptr->From->Type, "S"))
        JacElement(Mptr, i + 3, k, -Vk * Vm * Be * sin(thk - thm));
      if (!strpbrk(TCSCptr->To->Type, "S"))
        JacElement(Mptr, i + 3, m, Vk * Vm * Be * sin(thk - thm));
      JacElement(Mptr, i + 3, i + 3, -1.0);
      JacElement(Mptr, i + 3, i + 4, -Vm * Vm + Vk * Vm * cos(thk - thm));
      JacElement(Mptr, i + 4, i + 4, -1.0);
      /*  JacElement(Mptr,i+4,i+5,(2.0-2.0*cos(2*alpha))/(PI*Xl));  */
      s1 = -PI * sin(Kf * (-PI + alpha)) * Kf * (pow(Kf, 4.0) - 2.0 * Kf * Kf + 1.0) / Xc / (-PI * pow(Kf, 4.0) * cos(Kf * (-PI + alpha)) + PI * cos(Kf * (-PI + alpha)) + 2.0 * pow(Kf, 4.0) * alpha * cos(Kf * (-PI + alpha)) - 2.0 * alpha * Kf * Kf * cos(Kf * (-PI + alpha)) + pow(Kf, 4.0) * sin(-2.0 * PI + 2.0 * alpha) * cos(Kf * (-PI + alpha)) - sin(-2.0 * PI + 2.0 * alpha) * Kf * Kf * cos(Kf * (-PI + alpha)) - 4.0 * Kf * Kf * Kf * pow(cos(-PI + alpha), 2.0) * sin(Kf * (-PI + alpha)) + 4.0 * Kf * Kf * cos(-PI + alpha) * sin(-PI + alpha) * cos(Kf * (-PI + alpha)));
      s3 = -PI * cos(Kf * (-PI + alpha));
      s5 = (pow(Kf, 4.0) - 2.0 * Kf * Kf + 1.0) / Xc;
      s7 = 1 / (pow(-PI * pow(Kf, 4.0) * cos(Kf * (-PI + alpha)) + PI * cos(Kf * (-PI + alpha)) + 2.0 * pow(Kf, 4.0) * alpha * cos(Kf * (-PI + alpha)) - 2.0 * alpha * Kf * Kf * cos(Kf * (-PI + alpha)) + pow(Kf, 4.0) * sin(-2.0 * PI + 2.0 * alpha) * cos(Kf * (-PI + alpha)) - sin(-2.0 * PI + 2.0 * alpha) * Kf * Kf * cos(Kf * (-PI + alpha)) - 4.0 * Kf * Kf * Kf * pow(cos(-PI + alpha), 2.0) * sin(Kf * (-PI + alpha)) + 4.0 * Kf * Kf * cos(-PI + alpha) * sin(-PI + alpha) * cos(Kf * (-PI + alpha)), 2.0));
      s8 = PI * pow(Kf, 5.0) * sin(Kf * (-PI + alpha)) - PI * sin(Kf * (-PI + alpha)) * Kf + 2.0 * pow(Kf, 4.0) * cos(Kf * (-PI + alpha)) - 2.0 * pow(Kf, 5.0) * alpha * sin(Kf * (-PI + alpha)) - 2.0 * Kf * Kf * cos(Kf * (-PI + alpha)) + 2.0 * alpha * Kf * Kf * Kf * sin(Kf * (-PI + alpha)) + 2.0 * pow(Kf, 4.0) * cos(-2.0 * PI + 2.0 * alpha) * cos(Kf * (-PI + alpha)) - pow(Kf, 5.0) * sin(-2.0 * PI + 2.0 * alpha) * sin(Kf * (-PI + alpha)) - 2.0 * cos(-2.0 * PI + 2.0 * alpha) * Kf * Kf * cos(Kf * (-PI + alpha)) + sin(-2.0 * PI + 2.0 * alpha) * Kf * Kf * Kf * sin(Kf * (-PI + alpha)) + 4.0 * Kf * Kf * Kf * cos(-PI + alpha) * sin(Kf * (-PI + alpha)) * sin(-PI + alpha) - 4.0 * pow(Kf, 4.0) * pow(cos(-PI + alpha), 2.0) * cos(Kf * (-PI + alpha)) - 4.0 * Kf * Kf * pow(sin(-PI + alpha), 2.0) * cos(Kf * (-PI + alpha)) + 4.0 * Kf * Kf * pow(cos(-PI + alpha), 2.0) * cos(Kf * (-PI + alpha));
      s6 = s7 * s8;
      s4 = s5 * s6;
      s2 = s3 * s4;
      JacElement(Mptr, i + 4, i + 5, s1 + s2);
      if (Itcsc >= 0)
      {
        if (Ptcsc == 0)
          JacElement(Mptr, i + 5, i + 1, 0.);
        else
          JacElement(Mptr, i + 5, i + 1, Ptcsc / (sqrt(Ptcsc * Ptcsc + Qtcsck * Qtcsck)));
        if (Qtcsck == 0)
          JacElement(Mptr, i + 5, i + 2, 0.);
        else
          JacElement(Mptr, i + 5, i + 2, Qtcsck / (sqrt(Ptcsc * Ptcsc + Qtcsck * Qtcsck)));
      }
      else
      {
        if (Ptcsc == 0)
          JacElement(Mptr, i + 5, i + 1, 0.);
        else
          JacElement(Mptr, i + 5, i + 1, -Ptcsc / (sqrt(Ptcsc * Ptcsc + Qtcsck * Qtcsck)));
        if (Qtcsck == 0)
          JacElement(Mptr, i + 5, i + 2, 0.);
        else
          JacElement(Mptr, i + 5, i + 2, -Qtcsck / (sqrt(Ptcsc * Ptcsc + Qtcsck * Qtcsck)));
      }
      JacElement(Mptr, i + 5, i + 6, -Vk);
      if (!strpbrk(TCSCptr->From->Type, "S"))
        JacElement(Mptr, i + 6, k, 1.0);
      if (!strpbrk(TCSCptr->To->Type, "S"))
        JacElement(Mptr, i + 6, m, -1.0);
      JacElement(Mptr, i + 6, i + 7, -1.0);
      if (!strcmp(TCSCptr->Cont, "X"))
        JacElement(Mptr, i + 7, i + 4, -1.0);
      else
      {
        JacElement(Mptr, i + 7, i + 4, 0.0);
        if (!strcmp(TCSCptr->Cont, "P"))
          JacElement(Mptr, i + 7, i + 1, -1.0);
        else if (!strcmp(TCSCptr->Cont, "I"))
          JacElement(Mptr, i + 7, i + 6, -1.0);
        else if (!strcmp(TCSCptr->Cont, "D"))
          JacElement(Mptr, i + 7, i + 7, -1.0);
      }
      if (flagH)
      {
        if (strpbrk(TCSCptr->From->Type, "L"))
        {
          JacElement(Mptr, i + 1, Mptr->n1, -Vm * Be * sin(thk - thm));
          JacElement(Mptr, i + 2, Mptr->n1, -2.0 * Vk * Be + Vm * Be * cos(thk - thm));
          JacElement(Mptr, i + 3, Mptr->n1, Vm * Be * cos(thk - thm));
          JacElement(Mptr, i + 5, Mptr->n1, -Itcsc);
        }
        else if (strpbrk(TCSCptr->To->Type, "L"))
        {
          JacElement(Mptr, i + 1, Mptr->n1, -Vk * Be * sin(thk - thm));
          JacElement(Mptr, i + 2, Mptr->n1, Vk * Be * cos(thk - thm));
          JacElement(Mptr, i + 3, Mptr->n1, -2.0 * Vm * Be + Vk * Be * cos(thk - thm));
        }
      }
    }
    i = i + 7;
  }
}

/*  Construct the STATCOM part of the mismatch vector and Jacobian  */
void STATCOMFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ)
{
  INDEX i, k, l;
  STATCOMbusData *STATCOMptr;
  VALUETYPE Vk, delta, Vl, R, G, B, Gc, Xsl, Vref, P, Q, I, theta, Vdc, alpha, K, Cref;
  bool flagLimits, flagPWM;

  i = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar;
  for (STATCOMptr = dataPtr->STATCOMbus; STATCOMptr != nullptr; STATCOMptr = STATCOMptr->Next)
  {
    k = ACvar[STATCOMptr->From->N];
    Vk = STATCOMptr->From->V;
    delta = STATCOMptr->From->Ang;
    l = ACvar[STATCOMptr->Ctrl->N];
    Vl = STATCOMptr->Ctrl->V;
    if (!strcmp(STATCOMptr->Cont, "PW") || !strcmp(STATCOMptr->Cont, "AL"))
      flagLimits = false;
    else
      flagLimits = true;
    if (!flagLimits)
      Vref = STATCOMptr->Vref;
    else
      Vref = STATCOMptr->Vvar;
    if (!strcmp(STATCOMptr->Cont1, "PW"))
      flagPWM = true;
    else
      flagPWM = false;
    R = STATCOMptr->R;
    G = STATCOMptr->G;
    B = STATCOMptr->B;
    Gc = STATCOMptr->Gc;
    Xsl = STATCOMptr->slope;
    I = STATCOMptr->I;
    theta = STATCOMptr->theta;
    Vdc = STATCOMptr->Vdc;
    K = STATCOMptr->k;
    alpha = STATCOMptr->alpha;
    P = STATCOMptr->P;
    Q = STATCOMptr->Q;
    Cref = STATCOMptr->Contref;
    if (flagF)
    {
      if (Q > 0)
        dF[i + 1] = Vl - Vref - Xsl * I;
      else
        dF[i + 1] = Vl - Vref + Xsl * I;
      if (flagPWM)
        dF[i + 2] = Vdc - Cref;
      else
        dF[i + 2] = K - Cref;
      dF[i + 3] = P - Gc * Vdc * Vdc - R * I * I;
      dF[i + 4] = P - Vk * I * cos(delta - theta);
      dF[i + 5] = Q - Vk * I * sin(delta - theta);
      dF[i + 6] = P - G * Vk * Vk + G * K * Vdc * Vk * cos(delta - alpha) + B * K * Vdc * Vk * sin(delta - alpha);
      dF[i + 7] = Q + B * Vk * Vk - B * K * Vdc * Vk * cos(delta - alpha) + G * K * Vdc * Vk * sin(delta - alpha);
      dF[k] = dF[k] - P;
      dF[k + 1] = dF[k + 1] - Q;
    }
    if (flagJ)
    {
      if (!flagLimits)
      {
        if (Q > 0)
          JacElement(Mptr, i + 1, i + 1, -Xsl);
        else
          JacElement(Mptr, i + 1, i + 1, +Xsl);
      }
      else
        JacElement(Mptr, i + 1, i + 1, -1.0);
      JacElement(Mptr, i + 1, l + 1, 1.0);

      if (flagPWM)
        JacElement(Mptr, i + 2, i + 3, 1.0);
      else
        JacElement(Mptr, i + 2, i + 4, 1.0);

      if (!flagLimits)
        JacElement(Mptr, i + 3, i + 1, -2.0 * R * I);
      else
        JacElement(Mptr, i + 3, i + 1, 0.0);
      JacElement(Mptr, i + 3, i + 3, -2.0 * Gc * Vdc);
      JacElement(Mptr, i + 3, i + 6, 1.0);

      if (!flagLimits)
        JacElement(Mptr, i + 4, i + 1, -Vk * cos(delta - theta));
      else
        JacElement(Mptr, i + 4, i + 1, 0.0);
      JacElement(Mptr, i + 4, i + 2, -Vk * I * sin(delta - theta));
      JacElement(Mptr, i + 4, i + 6, 1.0);
      JacElement(Mptr, i + 4, k + 1, -I * cos(delta - theta));
      if (!strpbrk(STATCOMptr->From->Type, "S"))
        JacElement(Mptr, i + 4, k, Vk * I * sin(delta - theta));

      if (!flagLimits)
        JacElement(Mptr, i + 5, i + 1, -Vk * sin(delta - theta));
      else
        JacElement(Mptr, i + 5, i + 1, 0.0);
      JacElement(Mptr, i + 5, i + 2, Vk * I * cos(delta - theta));
      JacElement(Mptr, i + 5, i + 7, 1.0);
      JacElement(Mptr, i + 5, k + 1, -I * sin(delta - theta));
      if (!strpbrk(STATCOMptr->From->Type, "S"))
        JacElement(Mptr, i + 5, k, -Vk * I * cos(delta - theta));

      JacElement(Mptr, i + 6, i + 3, G * K * Vk * cos(delta - alpha) + B * K * Vk * sin(delta - alpha));
      JacElement(Mptr, i + 6, i + 4, G * Vdc * Vk * cos(delta - alpha) + B * Vdc * Vk * sin(delta - alpha));
      JacElement(Mptr, i + 6, i + 5, G * K * Vdc * Vk * sin(delta - alpha) - B * K * Vdc * Vk * cos(delta - alpha));
      JacElement(Mptr, i + 6, i + 6, 1.0);
      JacElement(Mptr, i + 6, k + 1, -2.0 * G * Vk + G * K * Vdc * cos(delta - alpha) + B * K * Vdc * sin(delta - alpha));
      if (!strpbrk(STATCOMptr->From->Type, "S"))
        JacElement(Mptr, i + 6, k, -G * K * Vdc * Vk * sin(delta - alpha) + B * K * Vdc * Vk * cos(delta - alpha));

      JacElement(Mptr, i + 7, i + 3, -B * K * Vk * cos(delta - alpha) + G * K * Vk * sin(delta - alpha));
      JacElement(Mptr, i + 7, i + 4, -B * Vdc * Vk * cos(delta - alpha) + G * Vdc * Vk * sin(delta - alpha));
      JacElement(Mptr, i + 7, i + 5, -B * K * Vdc * Vk * sin(delta - alpha) - G * K * Vdc * Vk * cos(delta - alpha));
      JacElement(Mptr, i + 7, i + 7, 1.0);
      JacElement(Mptr, i + 7, k + 1, 2.0 * B * Vk - B * K * Vdc * cos(delta - alpha) - G * K * Vdc * sin(delta - alpha));
      if (!strpbrk(STATCOMptr->From->Type, "S"))
        JacElement(Mptr, i + 7, k, B * K * Vdc * Vk * sin(delta - alpha) + G * K * Vdc * Vk * cos(delta - alpha));

      JacElement(Mptr, k, i + 6, -1.0);

      JacElement(Mptr, k + 1, i + 7, -1.0);

      if (flagH)
      {
        if (strpbrk(STATCOMptr->From->Type, "L"))
        {
          JacElement(Mptr, i + 4, Mptr->n1, -I * cos(delta - theta));
          JacElement(Mptr, i + 5, Mptr->n1, -I * sin(delta - theta));
          JacElement(Mptr, i + 6, Mptr->n1, -2.0 * G * Vk + G * K * Vdc * cos(delta - alpha) + B * K * Vdc * sin(delta - alpha));
          JacElement(Mptr, i + 7, Mptr->n1, 2.0 * B * Vk - B * K * Vdc * cos(delta - alpha) - G * K * Vdc * sin(delta - alpha));
        }
        if (strpbrk(STATCOMptr->Ctrl->Type, "L"))
        {
          JacElement(Mptr, i + 1, Mptr->n1, 1.0);
        }
      }
    }
    i = i + 7;
  }
}
