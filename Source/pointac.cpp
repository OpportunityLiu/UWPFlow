/* Point of Collapse:  AC Hessian. */

#include "pointl.h"

/* ------------------ ACFunHes ----------------------------- */
void ACFunHes(bool flagF, bool flagJ)
/* Construct the AC part of the PoC Jacobian. */
{
  SparseMatrixElement *Jptr;
  ACbusData *ACptr, *To, *From, *BEptr, *BSptr;
  AClist *ALptr;
  ElementList *ELptr;
  ElementData *Eptr;
  VALUETYPE Vi, Vj, di, dj, gij, bij, gsij, bsij, gji, bji;
  VALUETYPE v1, dv1, v2, dv2, v3, dv3, v4, dv4, v5, v6, dv6, v7, v8, dv8, v9,
      v10, dv10, v11;
  VALUETYPE Ddi, Dvi, Ddip, Dvip, Ddj, Dvj, Dval, a, b, Dlambda, DPg;
  VALUETYPE Ra, Xd, Xq, dg, Vd, Vq, Vr, Vim, Iq, Id, Ir, Iim, Ia;
  INDEX N, i, j, k, l;

  N = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar + 7 * Nstatcom; /*  FACTS  */
  for (v1 = 0, i = 1 + N; i < Jac->n1; i++) {
    if (fabs(x0[i - N]) > v1) {
      k = i - N;
      v1 = fabs(x0[k]);
      if (x0[k] > 0)
        v2 = 1.0;
      else
        v2 = -1.0;
    }
    if (flagJ)
      JacElement(Jac, Jac->n1, i, 0.);
    if (flagF)
      dF[i] = 0;
  }
  if (flagF)
    dF[Jac->n1] = v1 - 1;
  if (flagJ)
    JacElement(Jac, Jac->n1, k + N, v2);
  if (flagJ)
    for (i = 1; i <= Jac->n1; i++) {
      Jptr = Jac->ColHead[i];
      while (Jptr != nullptr) {
        j = Jptr->Col;
        if (OldCol->p[j])
          j = OldCol->p[j];
        l = Jptr->Row;
        if (OldRow->p[l])
          l = OldRow->p[l];
        if (flagJ && j <= N && l <= N)
          JacElement(Jac, j + N, l + N, Jptr->Value);
        Jptr = Jptr->ColNext;
      }
    }
  for (ALptr = dataPtr->KGbus; ALptr != nullptr; ALptr = ALptr->Next) {
    BSptr = ALptr->AC;
    if (strpbrk(ALptr->AC->Type, "S"))
      break;
  }
  for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    if (ACptr->Area != nullptr)
      BSptr = ACptr->Area->Slack;
    i = ACvar[ACptr->N];
    Vi = ACptr->V;
    di = ACptr->Ang;
    a = ACptr->a;
    b = ACptr->b;
    From = ACptr;
    Ddi = Dvi = Ddip = Dvip = 0;
    /*if (flagF) dF[Jac->n1]+=-ACptr->DPl*x0[i]-ACptr->DQl*x0[i+1];*/
    for (ELptr = ACptr->Elem; ELptr != nullptr; ELptr = ELptr->Next) {
      Eptr = ELptr->Eptr;
      if (Eptr->From == ACptr)
        To = Eptr->To;
      else
        To = Eptr->From;
      j = ACvar[To->N];
      Vj = To->V;
      dj = To->Ang;
      if (Eptr->From == From) {
        gij = (Eptr->G * cos(Eptr->Ang) - Eptr->B * sin(Eptr->Ang)) * Eptr->Tap;
        bij = (Eptr->G * sin(Eptr->Ang) + Eptr->B * cos(Eptr->Ang)) * Eptr->Tap;
        gsij = (Eptr->G1 + Eptr->G) * Eptr->Tap * Eptr->Tap - gij;
        bsij = (Eptr->B1 + Eptr->B) * Eptr->Tap * Eptr->Tap - bij;
        gji = (Eptr->G * cos(Eptr->Ang) + Eptr->B * sin(Eptr->Ang)) * Eptr->Tap;
        bji =
            (-Eptr->G * sin(Eptr->Ang) + Eptr->B * cos(Eptr->Ang)) * Eptr->Tap;
      } else {
        gij = (Eptr->G * cos(Eptr->Ang) + Eptr->B * sin(Eptr->Ang)) * Eptr->Tap;
        bij =
            (-Eptr->G * sin(Eptr->Ang) + Eptr->B * cos(Eptr->Ang)) * Eptr->Tap;
        gsij = Eptr->G + Eptr->G2 - gij;
        bsij = Eptr->B + Eptr->B2 - bij;
        gji = (Eptr->G * cos(Eptr->Ang) - Eptr->B * sin(Eptr->Ang)) * Eptr->Tap;
        bji = (Eptr->G * sin(Eptr->Ang) + Eptr->B * cos(Eptr->Ang)) * Eptr->Tap;
      }
      v1 = sin(di - dj) * (gij * x0[i] + gji * x0[j]);
      dv1 = cos(di - dj) * (gij * x0[i] + gji * x0[j]);
      v2 = cos(di - dj) * (bij * x0[i] - bji * x0[j]);
      dv2 = -sin(di - dj) * (bij * x0[i] - bji * x0[j]);
      v3 = cos(di - dj) * (gij * x0[i + 1] - gji * x0[j + 1]);
      dv3 = -sin(di - dj) * (gij * x0[i + 1] - gji * x0[j + 1]);
      v4 = sin(di - dj) * (bij * x0[i + 1] + bji * x0[j + 1]);
      dv4 = cos(di - dj) * (bij * x0[i + 1] + bji * x0[j + 1]);
      v5 = -2 * ((gij + gsij) * x0[i] - (bij + bsij) * x0[i + 1]);
      if (Acont && From->Area != To->Area) {
        BEptr = From->Area->Slack;
        if (!strpbrk(BEptr->Type, "S")) {
          k = ACvar[BEptr->N] + 2;
          if (Eptr->Meter == From) {
            v6 = -(gij * sin(di - dj) - bij * cos(di - dj)) * x0[k];
            dv6 = -(gij * cos(di - dj) + bij * sin(di - dj)) * x0[k];
            v7 = -2 * (gij + gsij) * x0[k];
          } else {
            v6 = -(gji * sin(dj - di) - bji * cos(dj - di)) * x0[k];
            dv6 = (gji * cos(dj - di) + bji * sin(dj - di)) * x0[k];
            v7 = 0;
          }
        } else
          v6 = dv6 = v7 = 0;
        BEptr = To->Area->Slack;
        if (!strpbrk(BEptr->Type, "S")) {
          k = ACvar[BEptr->N] + 2;
          if (Eptr->Meter == From) {
            v8 = (gij * sin(di - dj) - bij * cos(di - dj)) * x0[k];
            dv8 = (gij * cos(di - dj) + bij * sin(di - dj)) * x0[k];
            v9 = 2 * (gij + gsij) * x0[k];
          } else {
            v8 = (gji * sin(dj - di) - bji * cos(dj - di)) * x0[k];
            dv8 = -(gji * cos(dj - di) + bji * sin(dj - di)) * x0[k];
            v9 = 0;
          }
        } else
          v8 = dv8 = v9 = 0;
      } else
        v6 = dv6 = v7 = v8 = dv8 = v9 = 0;
      if (PQcont && strpbrk(Eptr->Type, "PM")) {
        k = ACvar[Eptr->Cont->N] + 1 + Eptr->Cont->Ncont - Eptr->Ncont;
        if (Acont && strpbrk(Eptr->Cont->Type, "A"))
          k++;
        if (Eptr->Cont == From) {
          v10 = -(gij * sin(di - dj) - bij * cos(di - dj)) * x0[k];
          dv10 = -(gij * cos(di - dj) + bij * sin(di - dj)) * x0[k];
          v11 = -2 * (gij + gsij) * x0[k];
        } else {
          v10 = (gji * sin(dj - di) - bji * cos(dj - di)) * x0[k];
          dv10 = -(gji * cos(dj - di) + bji * sin(dj - di)) * x0[k];
          v11 = 0;
        }
      } else if (PQcont && strpbrk(Eptr->Type, "QN")) {
        k = ACvar[Eptr->Cont->N] + 1 + Eptr->Cont->Ncont - Eptr->Ncont;
        if (Acont && strpbrk(Eptr->Cont->Type, "A"))
          k++;
        if (Eptr->Cont == From) {
          v10 = (gij * cos(di - dj) + bij * sin(di - dj)) * x0[k];
          dv10 = -(gij * sin(di - dj) - bij * cos(di - dj)) * x0[k];
          v11 = 2 * (bij + bsij) * x0[k];
        } else {
          v10 = -(gji * cos(dj - di) + bji * sin(dj - di)) * x0[k];
          dv10 = -(gji * sin(dj - di) - bji * cos(dj - di)) * x0[k];
          v11 = 0;
        }
      } else
        v10 = dv10 = v11 = 0;
      if (flagF) {
        if (!strpbrk(From->Type, "S"))
          dF[i + N] += Vi * Vj * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
        if (From->Cont != nullptr)
          dF[i + N + 1] += Vj * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10) +
                           Vi * (v5 + v7 + v9 + v11);
      }
      if (flagJ) {
        if (!strpbrk(From->Type, "S")) {
          Ddj = Vi * Vj * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10);
          Ddi = Ddi - Ddj;
          Dvj = Vi * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
          Dvi = Dvi + Vj * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
          if (!strpbrk(To->Type, "S"))
            JacElement(Jac, i + N, j, Ddj);
          if (To->Cont != nullptr)
            JacElement(Jac, i + N, j + 1, Dvj);
          else
            JacElement(Jac, i + N, j + 1, 0.);
          if (PQcont && strpbrk(Eptr->Type, "PQMN")) {
            k = ACvar[Eptr->Cont->N] + 1 + Eptr->Cont->Ncont - Eptr->Ncont;
            if (Acont && strpbrk(Eptr->Cont->Type, "A"))
              k++;
            if (!strcmp(Eptr->Type, "RQ"))
              Dval = Vi * Vj / Eptr->Tap * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
            else if (!strcmp(Eptr->Type, "RP"))
              Dval = Vi * Vj * (-dv2 + dv1 - dv4 - dv3 - dv6 - dv8 - dv10);
            else
              Dval = 0;
            JacElement(Jac, i + N, k, Dval);
          } else if (Rcont && Eptr->Cont != nullptr) {
            k = ACvar[Eptr->Cont->N] + 1;
            if (!strcmp(Eptr->Type, "R"))
              Dval = Vi * Vj / Eptr->Tap * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
            else
              Dval = 0;
            JacElement(Jac, i + N, k, Dval);
          }
        }
        if (From->Cont != nullptr) {
          Ddj = -Vj * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
          Ddip = Ddip - Ddj;
          Dvj = dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10;
          Dvip = Dvip + v5 + v7 + v9 + v11;
          if (!strpbrk(To->Type, "S"))
            JacElement(Jac, i + 1 + N, j, Ddj);
          if (To->Cont != nullptr)
            JacElement(Jac, i + 1 + N, j + 1, Dvj);
          else
            JacElement(Jac, i + 1 + N, j + 1, 0.);
          if (PQcont && strpbrk(Eptr->Type, "PQMN")) {
            k = ACvar[Eptr->Cont->N] + 1 + Eptr->Cont->Ncont - Eptr->Ncont;
            if (Acont && strpbrk(Eptr->Cont->Type, "A"))
              k++;
            if (!strcmp(Eptr->Type, "RQ")) {
              Dval =
                  Vj / Eptr->Tap * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10);
              if (Eptr->From == From)
                Dval = Dval + 2 * Vi / Eptr->Tap * (v5 + v7 + v9 + v11);
            } else if (!strcmp(Eptr->Type, "RP"))
              Dval = Vj * (-v2 + v1 - v4 - v3 - v6 - v8 - v10);
            else
              Dval = 0;
            JacElement(Jac, i + 1 + N, k, Dval);
          } else if (Rcont && Eptr->Cont != nullptr) {
            k = ACvar[Eptr->Cont->N] + 1;
            if (!strcmp(Eptr->Type, "R")) {
              Dval =
                  Vj / Eptr->Tap * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10);
              if (Eptr->From == From)
                Dval = Dval + 2 * Vi / Eptr->Tap * (v5 + v7 + v9 + v11);
            } else
              Dval = 0;
            JacElement(Jac, i + 1 + N, k, Dval);
          }
        } else {
          JacElement(Jac, i + 1 + N, j, 0.);
          JacElement(Jac, i + 1 + N, j + 1, 0.);
          if (PQcont && strpbrk(Eptr->Type, "PQMN")) {
            k = ACvar[Eptr->Cont->N] + 1 + Eptr->Cont->Ncont - Eptr->Ncont;
            if (Acont && strpbrk(Eptr->Cont->Type, "A"))
              k++;
            JacElement(Jac, i + 1 + N, k, 0.);
          } else if (Rcont && Eptr->Cont != nullptr) {
            k = ACvar[Eptr->Cont->N] + 1;
            JacElement(Jac, i + 1 + N, k, 0.);
          }
        }
      }
    }
    if (flagF) {
      DPg = ACptr->DPG;
      if (DPg) {
        k = ACvar[BSptr->N];
        if (strpbrk(BSptr->Type, "S"))
          dF[k + N] += DPg * x0[i];
        else if (Acont)
          dF[k + N + 2] += DPg * x0[i];
      }
      if (ACptr->Cont != nullptr)
        dF[i + N + 1] +=
            -2 * Vi * (ACptr->G * x0[i] - ACptr->B * x0[i + 1]) -
            2 * Vi * (ACptr->Pz + ACptr->Pzl * lambda) * x0[i] -
            a * pow(Vi, a - 1) * (ACptr->Pn + ACptr->Pnl * lambda) * x0[i] -
            2 * Vi * (ACptr->Qz + ACptr->Qzl * lambda) * x0[i + 1] -
            b * pow(Vi, b - 1) * (ACptr->Qn + ACptr->Qnl * lambda) * x0[i + 1];
      if (QRcont && strpbrk(ACptr->Type, "G")) {
        j = ACvar[ACptr->Cont->N];
        if (strpbrk(ACptr->cont, "V"))
          dF[j + N + 1] += ACptr->Kbg * x0[i + 1];
        else if (ACptr->Gen != nullptr) {
          j = ACptr->Gen->Nvar;
          if (strpbrk(ACptr->cont, "E"))
            dF[j + N + 1] += x0[i + 1];
          else
            dF[j + N + 11] += x0[i + 1];
        }
      } else if (strpbrk(ACptr->Type, "V,Q,S,Z") ||
                 (!QRcont && strpbrk(ACptr->Type, "G"))) {
        if (strpbrk(ACptr->cont, "V"))
          dF[i + N + 1] += x0[i + 1];
        else if (ACptr->Gen != nullptr) {
          j = ACptr->Gen->Nvar;
          if (strpbrk(ACptr->cont, "E"))
            dF[j + N + 1] += x0[i + 1];
          else
            dF[j + N + 11] += x0[i + 1];
        }
      }
    }
    if (flagJ) {
      if (!strpbrk(From->Type, "S"))
        JacElement(Jac, i + N, i, Ddi);
      if (From->Cont != nullptr) {
        Dvip = Dvip - 2 * (ACptr->G * x0[i] - ACptr->B * x0[i + 1]) -
               2 * (ACptr->Pz + ACptr->Pzl * lambda) * x0[i] -
               a * (a - 1) * pow(Vi, a - 2) *
                   (ACptr->Pn + ACptr->Pnl * lambda) * x0[i] -
               2 * (ACptr->Qz + ACptr->Qzl * lambda) * x0[i + 1] -
               b * (b - 1) * pow(Vi, b - 2) *
                   (ACptr->Qn + ACptr->Qnl * lambda) * x0[i + 1];
        JacElement(Jac, i + N, i + 1, Dvi);
        JacElement(Jac, i + 1 + N, i, Ddip);
        JacElement(Jac, i + 1 + N, i + 1, Dvip);
        Dlambda =
            -(2 * Vi * ACptr->Pzl + a * pow(Vi, a - 1) * ACptr->Pnl) * x0[i] -
            (2 * Vi * ACptr->Qzl + b * pow(Vi, b - 1) * ACptr->Qnl) * x0[i + 1];
        if (Dlambda) {
          JacElement(Jac, i + 1 + N, Jac->n1, Dlambda);
          /*JacElement(Jac,Jac->n1,i+1,Dlambda);*/
        }
      } else {
        if (!strpbrk(From->Type, "S")) {
          JacElement(Jac, i + N, i + 1, 0.);
          JacElement(Jac, i + 1 + N, i, 0.);
        }
        JacElement(Jac, i + 1 + N, i + 1, 0.);
      }
      Dlambda = -Vi * Vi * ACptr->Pzl - pow(Vi, a) * ACptr->Pnl;
      if (Dlambda) {
        JacElement(Jac, i, Jac->n1, Dlambda);
        /*JacElement(Jac,Jac->n1,i+N,Dlambda);*/
      }
      Dlambda = -Vi * Vi * ACptr->Qzl - pow(Vi, b) * ACptr->Qnl;
      if (Dlambda) {
        JacElement(Jac, i + 1, Jac->n1, Dlambda);
        /*JacElement(Jac,Jac->n1,i+1+N,Dlambda);*/
      }
    }

    /* -------------- Regulating Transf. ----------------------- */
    if (PQcont || Rcont)
      for (ELptr = ACptr->Reg; ELptr != nullptr; ELptr = ELptr->Next) {
        Eptr = ELptr->Eptr;
        From = Eptr->From;
        i = ACvar[From->N];
        Vi = From->V;
        di = From->Ang;
        To = Eptr->To;
        j = ACvar[To->N];
        Vj = To->V;
        dj = To->Ang;
        gij = (Eptr->G * cos(Eptr->Ang) - Eptr->B * sin(Eptr->Ang)) * Eptr->Tap;
        bij = (Eptr->G * sin(Eptr->Ang) + Eptr->B * cos(Eptr->Ang)) * Eptr->Tap;
        gsij = (Eptr->G1 + Eptr->G) * Eptr->Tap * Eptr->Tap - gij;
        bsij = (Eptr->B1 + Eptr->B) * Eptr->Tap * Eptr->Tap - bij;
        gji = (Eptr->G * cos(Eptr->Ang) + Eptr->B * sin(Eptr->Ang)) * Eptr->Tap;
        bji =
            (-Eptr->G * sin(Eptr->Ang) + Eptr->B * cos(Eptr->Ang)) * Eptr->Tap;
        v1 = sin(di - dj) * (gij * x0[i] + gji * x0[j]);
        dv1 = cos(di - dj) * (gij * x0[i] + gji * x0[j]);
        v2 = cos(di - dj) * (bij * x0[i] - bji * x0[j]);
        dv2 = -sin(di - dj) * (bij * x0[i] - bji * x0[j]);
        v3 = cos(di - dj) * (gij * x0[i + 1] - gji * x0[j + 1]);
        dv3 = -sin(di - dj) * (gij * x0[i + 1] - gji * x0[j + 1]);
        v4 = sin(di - dj) * (bij * x0[i + 1] + bji * x0[j + 1]);
        dv4 = cos(di - dj) * (bij * x0[i + 1] + bji * x0[j + 1]);
        v5 = -2 * ((gij + gsij) * x0[i] - (bij + bsij) * x0[i + 1]);
        if (Acont && From->Area != To->Area) {
          BEptr = From->Area->Slack;
          if (!strpbrk(BEptr->Type, "S")) {
            k = ACvar[BEptr->N] + 2;
            if (Eptr->Meter == From) {
              v6 = -(gij * sin(di - dj) - bij * cos(di - dj)) * x0[k];
              dv6 = -(gij * cos(di - dj) + bij * sin(di - dj)) * x0[k];
              v7 = -2 * (gij + gsij) * x0[k];
            } else {
              v6 = -(gji * sin(dj - di) - bji * cos(dj - di)) * x0[k];
              dv6 = (gji * cos(dj - di) + bji * sin(dj - di)) * x0[k];
              v7 = 0;
            }
          } else
            v6 = dv6 = v7 = 0;
          BEptr = To->Area->Slack;
          if (!strpbrk(BEptr->Type, "S")) {
            k = ACvar[BEptr->N] + 2;
            if (Eptr->Meter == From) {
              v8 = (gij * sin(di - dj) - bij * cos(di - dj)) * x0[k];
              dv8 = (gij * cos(di - dj) + bij * sin(di - dj)) * x0[k];
              v9 = 2 * (gij + gsij) * x0[k];
            } else {
              v8 = (gji * sin(dj - di) - bji * cos(dj - di)) * x0[k];
              dv8 = -(gji * cos(dj - di) + bji * sin(dj - di)) * x0[k];
              v9 = 0;
            }
          } else
            v8 = dv8 = v9 = 0;
        } else
          v6 = dv6 = v7 = v8 = dv8 = v9 = 0;
        if (PQcont && !strcmp(Eptr->Type, "RP")) {
          k = ACvar[Eptr->Cont->N] + 1 + Eptr->Cont->Ncont - Eptr->Ncont;
          if (Acont && strpbrk(Eptr->Cont->Type, "A"))
            k++;
          if (Eptr->Cont == From) {
            v10 = -(gij * sin(di - dj) - bij * cos(di - dj)) * x0[k];
            dv10 = -(gij * cos(di - dj) + bij * sin(di - dj)) * x0[k];
          } else {
            v10 = (gji * sin(dj - di) - bji * cos(dj - di)) * x0[k];
            dv10 = -(gji * cos(dj - di) + bji * sin(dj - di)) * x0[k];
          }
          if (flagF)
            dF[k + N] = Vi * Vj * (v1 - v2 - v3 - v4 - v6 - v8 - v10);
          if (flagJ) {
            Dval = Vi * Vj * (dv2 - dv1 + dv4 + dv3 + dv6 + dv8 + dv10);
            JacElement(Jac, k + N, k, Dval);
            if (!strpbrk(From->Type, "S")) {
              Ddi = Vi * Vj * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10);
              JacElement(Jac, k + N, i, Ddi);
            }
            if (From->Cont != nullptr) {
              Dvi = -Vj * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
              JacElement(Jac, k + N, i + 1, Dvi);
            } else
              JacElement(Jac, k + N, i + 1, 0.);
            if (!strpbrk(To->Type, "S")) {
              Ddj = -Vi * Vj * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10);
              JacElement(Jac, k + N, j, Ddj);
            }
            if (To->Cont != nullptr) {
              Dvj = -Vi * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
              JacElement(Jac, k + N, j + 1, Dvj);
            } else
              JacElement(Jac, k + N, j + 1, 0.);
          }
        } else if (PQcont && !strcmp(Eptr->Type, "RQ")) {
          k = ACvar[Eptr->Cont->N] + 1 + Eptr->Cont->Ncont - Eptr->Ncont;
          if (Acont && strpbrk(Eptr->Cont->Type, "A"))
            k++;
          if (Eptr->Cont == From) {
            v10 = (gij * cos(di - dj) + bij * sin(di - dj)) * x0[k];
            dv10 = -(gij * sin(di - dj) - bij * cos(di - dj)) * x0[k];
            v11 = 2 * (bij + bsij) * x0[k];
          } else {
            v10 = -(gji * cos(dj - di) + bji * sin(dj - di)) * x0[k];
            dv10 = -(gji * sin(dj - di) - bji * cos(dj - di)) * x0[k];
            v11 = 0;
          }
          if (flagF)
            dF[k + N] += Vi * Vj / Eptr->Tap *
                             (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10) +
                         Vi * Vi / Eptr->Tap * (v5 + v7 + v9 + v11);
          if (flagJ) {
            Dval = Vi * Vi / (Eptr->Tap * Eptr->Tap) * (v5 + v7 + v9 + v11);
            JacElement(Jac, k + N, k, Dval);
            if (!strpbrk(From->Type, "S")) {
              Ddi = Vi * Vj / Eptr->Tap * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
              JacElement(Jac, k + N, i, Ddi);
            }
            if (From->Cont != nullptr) {
              Dvi =
                  Vj / Eptr->Tap * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10) +
                  2 * Vi / Eptr->Tap * (v5 + v7 + v9 + v11);
              JacElement(Jac, k + N, i + 1, Dvi);
            } else
              JacElement(Jac, k + N, i + 1, 0.);
            if (!strpbrk(To->Type, "S")) {
              Ddj = -Vi * Vj / Eptr->Tap * (-v1 + v2 + v3 + v4 + v6 + v8 + v10);
              JacElement(Jac, k + N, j, Ddj);
            }
            if (To->Cont != nullptr) {
              Dvj = Vi / Eptr->Tap * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8 - dv10);
              JacElement(Jac, k + N, j + 1, Dvj);
            } else
              JacElement(Jac, k + N, j + 1, 0.);
          }
        } else if (Rcont && !strcmp(Eptr->Type, "R")) {
          k = ACvar[Eptr->Cont->N] + 1;
          if (flagF)
            dF[k + N] +=
                Vi * Vj / Eptr->Tap * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8) +
                Vi * Vi / Eptr->Tap * (v5 + v7 + v9);
          if (flagJ) {
            Dval = Vi * Vi / (Eptr->Tap * Eptr->Tap) * (v5 + v7 + v9);
            JacElement(Jac, k + N, k, Dval);
            if (!strpbrk(From->Type, "S")) {
              Ddi = Vi * Vj / Eptr->Tap * (-v1 + v2 + v3 + v4 + v6 + v8);
              JacElement(Jac, k + N, i, Ddi);
            }
            if (From->Cont != nullptr) {
              Dvi = Vj / Eptr->Tap * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8) +
                    2 * Vi / Eptr->Tap * (v5 + v7 + v9);
              JacElement(Jac, k + N, i + 1, Dvi);
            } else
              JacElement(Jac, k + N, i + 1, 0.);
            if (!strpbrk(To->Type, "S")) {
              Ddj = -Vi * Vj / Eptr->Tap * (-v1 + v2 + v3 + v4 + v6 + v8);
              JacElement(Jac, k + N, j, Ddj);
            }
            if (To->Cont != nullptr) {
              Dvj = Vi / Eptr->Tap * (dv1 - dv2 - dv3 - dv4 - dv6 - dv8);
              JacElement(Jac, k + N, j + 1, Dvj);
            } else
              JacElement(Jac, k + N, j + 1, 0.);
          }
        } else if (flagF && PQcont && strpbrk(Eptr->Type, "PQMN")) {
          k = ACvar[Eptr->Cont->N] + 1 + Eptr->Cont->Ncont - Eptr->Ncont;
          if (Acont && strpbrk(Eptr->Cont->Type, "A"))
            k++;
          dF[k + N] = x0[k];
        }
      }

    /* -------------- Generator Model ----------------------- */
    if (ACptr->Gen != nullptr) {
      i = ACptr->Gen->Nvar;
      Ra = ACptr->Gen->Ra;
      Xd = ACptr->Gen->Xd;
      Xq = ACptr->Gen->Xq;
      dg = ACptr->Gen->dg;
      Vr = ACptr->Gen->Vr;
      Vim = ACptr->Gen->Vi;
      Ir = ACptr->Gen->Ir;
      Iim = ACptr->Gen->Ii;
      Vq = ACptr->Gen->Vq;
      Vd = ACptr->Gen->Vd;
      Iq = ACptr->Gen->Iq;
      Id = ACptr->Gen->Id;
      Ia = ACptr->Gen->Ia;
      if (flagF) {
        if (strpbrk(ACptr->cont, "V")) {
          if (QRcont && strpbrk(ACptr->Type, "G")) {
            j = ACvar[ACptr->Cont->N];
            dF[j + N + 1] += ACptr->Kbg * x0[i + 2];
          } else {
            j = ACvar[ACptr->N];
            dF[j + N + 1] += x0[i + 2];
          }
        }
        if (strpbrk(ACptr->cont, "E"))
          dF[i + N + 1] += x0[i + 2];
        else
          dF[i + N + 1] += x0[i + 3];
        dF[i + N + 2] += (sin(dg) * Vq + cos(dg) * Vd) * x0[i + 5] +
                         (-cos(dg) * Vq + sin(dg) * Vd) * x0[i + 6] +
                         (sin(dg) * Iq + cos(dg) * Id) * x0[i + 7] +
                         (-cos(dg) * Iq + sin(dg) * Id) * x0[i + 8];
        dF[i + N + 3] +=
            -Ir * x0[i + 1] + Iim * x0[i + 2] + x0[i + 5] + x0[i + 9];
        dF[i + N + 4] +=
            -Iim * x0[i + 1] - Ir * x0[i + 2] + x0[i + 6] + x0[i + 10];
        dF[i + N + 5] +=
            -Vr * x0[i + 1] - Vim * x0[i + 2] + x0[i + 7] - 2 * Ir * x0[i + 11];
        dF[i + N + 6] += -Vim * x0[i + 1] + Vr * x0[i + 2] + x0[i + 8] -
                         2 * Iim * x0[i + 11];
        dF[i + N + 7] += -x0[i + 3] - cos(dg) * x0[i + 5] - sin(dg) * x0[i + 6];
        dF[i + N + 8] += x0[i + 4] + sin(dg) * x0[i + 5] - cos(dg) * x0[i + 6];
        dF[i + N + 9] += -Ra * x0[i + 3] + Xq * x0[i + 4] -
                         cos(dg) * x0[i + 7] - sin(dg) * x0[i + 8];
        dF[i + N + 10] += Xd * x0[i + 3] + Ra * x0[i + 4] +
                          sin(dg) * x0[i + 7] - cos(dg) * x0[i + 8];
        if (strpbrk(ACptr->cont, "I"))
          dF[i + N + 11] += x0[i + 2];
        else
          dF[i + N + 11] += 2 * Ia * x0[i + 11];
        j = ACvar[BSptr->N];
        if (DPg) {
          if (strpbrk(BSptr->Type, "S"))
            dF[j + N] += DPg * x0[i + 1];
          else if (Acont)
            dF[j + N + 2] += DPg * x0[i + 1];
        }
        j = ACvar[ACptr->N];
        if (!strpbrk(ACptr->Type, "S"))
          dF[j + N] += Vi * sin(di) * x0[i + 9] - Vi * cos(di) * x0[i + 10];
        if (ACptr->Cont != nullptr)
          dF[j + N + 1] += -cos(di) * x0[i + 9] - sin(di) * x0[i + 10];
      }
      if (flagJ) {
        JacElement(Jac, i + N + 2, i + 2,
                   (cos(dg) * Vq - sin(dg) * Vd) * x0[i + 5] +
                       (sin(dg) * Vq + cos(dg) * Vd) * x0[i + 6] +
                       (cos(dg) * Iq - sin(dg) * Id) * x0[i + 7] +
                       (sin(dg) * Iq + cos(dg) * Id) * x0[i + 8]);
        JacElement(Jac, i + N + 2, i + 7,
                   sin(dg) * x0[i + 5] - cos(dg) * x0[i + 6]);
        JacElement(Jac, i + N + 2, i + 8,
                   cos(dg) * x0[i + 5] + sin(dg) * x0[i + 6]);
        JacElement(Jac, i + N + 2, i + 9,
                   sin(dg) * x0[i + 7] - cos(dg) * x0[i + 8]);
        JacElement(Jac, i + N + 2, i + 10,
                   cos(dg) * x0[i + 7] + sin(dg) * x0[i + 8]);
        JacElement(Jac, i + N + 3, i + 5, -x0[i + 1]);
        JacElement(Jac, i + N + 3, i + 6, x0[i + 2]);
        JacElement(Jac, i + N + 4, i + 5, -x0[i + 2]);
        JacElement(Jac, i + N + 4, i + 6, -x0[i + 1]);
        JacElement(Jac, i + N + 5, i + 3, -x0[i + 1]);
        JacElement(Jac, i + N + 5, i + 4, -x0[i + 2]);
        JacElement(Jac, i + N + 5, i + 5, -2 * x0[i + 11]);
        JacElement(Jac, i + N + 6, i + 3, x0[i + 2]);
        JacElement(Jac, i + N + 6, i + 4, -x0[i + 1]);
        JacElement(Jac, i + N + 6, i + 6, -2 * x0[i + 11]);
        JacElement(Jac, i + N + 7, i + 2,
                   sin(dg) * x0[i + 5] - cos(dg) * x0[i + 6]);
        JacElement(Jac, i + N + 8, i + 2,
                   cos(dg) * x0[i + 5] + sin(dg) * x0[i + 6]);
        JacElement(Jac, i + N + 9, i + 2,
                   sin(dg) * x0[i + 7] - cos(dg) * x0[i + 8]);
        JacElement(Jac, i + N + 10, i + 2,
                   cos(dg) * x0[i + 7] + sin(dg) * x0[i + 8]);
        if (!strpbrk(ACptr->cont, "I"))
          JacElement(Jac, i + N + 11, i + 11, 2 * x0[i + 11]);
        else
          JacElement(Jac, i + N + 11, i + 11, 0.);
        j = ACvar[ACptr->N];
        if (!strpbrk(ACptr->Type, "S")) {
          JacElement(Jac, j + N, j,
                     Vi * cos(di) * x0[i + 9] + Vi * sin(di) * x0[i + 10]);
          if (ACptr->Cont != nullptr) {
            JacElement(Jac, j + N, j + 1,
                       sin(di) * x0[i + 9] - cos(di) * x0[i + 10]);
            JacElement(Jac, j + N + 1, j,
                       sin(di) * x0[i + 9] - cos(di) * x0[i + 10]);
          } else {
            JacElement(Jac, j + N, j + 1, 0.);
            JacElement(Jac, j + N + 1, j, 0.);
          }
        }
      }
    }
  }
}
