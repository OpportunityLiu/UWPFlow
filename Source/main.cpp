/* AC/DC Power Flow Program.

*/

#include "constant.h"
#include "param.h"
#include "pflow.h"
#include "sparse.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ErrorStop(char *Msg);
void DCinit(void);
void SVCinit(void);     /* FACTS */
void TCSCinit(void);    /* FACTS */
void STATCOMinit(void); /* FACTS */
void Jacobian(void);
void ReadData(char *Name);
bool ReadInit(void);
void WriteSolution(INDEX Iter, char *File1, char *str);
int Pflow(int iter, bool flagF, bool flagD, bool flagP);
int Homot(void);
int PoCPoint(void);
int main(int argc, const char **argv);
void InitializeLoad(void);
int Evector(int M, int iter, VALUETYPE tol, bool RightEvector,
            VALUETYPE *EigenValue);
void PrintLeftEvector(INDEX N, FILE *Out);
void PrintDirection(char Option, VALUETYPE *vector, VALUETYPE Max);
void CleanUp(void);
void DeleteJac(SparseMatrix *Mptr, IntegerVector *P1Row, IntegerVector *P1Col,
               IntegerVector *P2Row, IntegerVector *P2Col);

/* ------- Global Variables ------ */
Data *dataPtr = nullptr;
SparseMatrix *Jac = nullptr;
INDEX MaxIter, Nac, NacEl, NregPQ, NregV, Ngen, Ndc, Nslack, Nvolt, Narea,
    NacVar, Bl, *ACvar, NZvolt, NXvolt, Nsvc, Ntcsc, NtcscVar,
    Nstatcom; /* FACTS */
ACbusData *BlPtr;
VALUETYPE *dx, *dF, tol, Tol, Sn, lambda;
VALUETYPE K1, K2, K3 = PI / 180.0, MaxdFi, alpha;
IntegerVector *NewRow = nullptr, *OldRow = nullptr, *NewCol = nullptr,
              *OldCol = nullptr, *RowPartition = nullptr,
              *ColPartition = nullptr;
bool Acont, PQcont, QRcont, Rcont, Xcont, PQlim, Tlim, Qlim, Vlim, Elim, Ilim,
    Zlim, Xlim, flag2Vcontrol, flagH, flagPoC, flagL, flagR, flagPgMax,
    flagSmax, flagReducedContinuation, flagVloads, flagKdirection, flagBS;
extern VALUETYPE *x0, *x0p, *Dx;
extern bool *DxZero;
extern AClist *Vlist, *Vlistp;

/* ------------------- DCinit -------------------- */
void DCinit(void)
/*    Intialize the DC system variables. */
{
  ACbusData *ACptrR, *ACptrI;
  DCbusData *DCptrR, *DCptrI;
  VALUETYPE KVr, KVi, Vr, Vi, ar, ai, cosar, cosgi, Xcr, Xci, Id, Rd, Ld, Vn,
      In;

  for (DCptrR = dataPtr->DCbus; DCptrR != nullptr; DCptrR = DCptrR->Next) {
    DCptrI = DCptrR->To;
    if (!strcmp(DCptrR->Type, "R")) {
      KVr = DCptrR->Vn;
      KVi = DCptrI->Vn;
      ACptrR = DCptrR->AC;
      ACptrI = DCptrI->AC;
      Vr = KVr * ACptrR->V;
      Vi = KVi * ACptrI->V;
      ar = DCptrR->Tap * DCptrR->Ntrf;
      ai = DCptrI->Tap * DCptrI->Ntrf;
      Xcr = DCptrR->Xc;
      Xci = DCptrI->Xc;
      cosar = cos(DCptrR->Alfa * K3);
      cosgi = cos(DCptrI->Gamma * K3);
      Id = DCptrR->Id / 1000;
      Rd = DCptrR->Rd;
      Ld = DCptrR->Ld;
      if (DCptrR->Vd <= 0)
        DCptrR->Vd = fabs(K1 * ar * Vr * cosar - K2 * Xcr * Id);
      if (DCptrI->Vd <= 0)
        DCptrI->Vd = fabs(K1 * ai * Vi * cosgi - K2 * Xci * Id);
      if (Id <= 0) {
        Id = fabs(DCptrR->Vd - DCptrI->Vd) / Rd;
        DCptrR->Id = DCptrI->Id = Id;
      }
      if (DCptrR->P <= 0)
        DCptrR->P = DCptrR->Vd * Id;
      if (DCptrI->P <= 0)
        DCptrI->P = DCptrI->Vd * Id;
      if (DCptrR->Q <= 0)
        DCptrR->Q = DCptrR->Vd * Id * tan(acos(cosar));
      if (DCptrI->Q <= 0)
        DCptrI->Q = DCptrI->Vd * Id * tan(acos(cosgi));
      /* ------------- Per unit conversion ----------------------------- */
      Vn = KVr;
      In = Sn / Vn;
      DCptrR->Vn = DCptrI->Vn = Vn;
      DCptrR->Rd = DCptrI->Rd = Rd * In / Vn;
      DCptrR->Ld = DCptrI->Ld = Ld * In / Vn;
      DCptrR->Id = DCptrI->Id = Id / In;
      DCptrR->Vd = DCptrR->Vd / Vn;
      DCptrI->Vd = DCptrI->Vd / Vn;
      DCptrR->Xc = DCptrR->Xc * In / Vn;
      DCptrI->Xc = DCptrI->Xc * In / Vn;
      DCptrR->Alfa = DCptrR->Alfa * K3;
      DCptrI->Alfa = DCptrI->Alfa * K3;
      DCptrR->Gamma = DCptrR->Gamma * K3;
      DCptrI->Gamma = DCptrI->Gamma * K3;
      DCptrR->AlfaMin = DCptrR->AlfaMin * K3;
      DCptrI->AlfaMin = DCptrI->AlfaMin * K3;
      DCptrR->AlfaMax = DCptrR->AlfaMax * K3;
      DCptrI->AlfaMax = DCptrI->AlfaMax * K3;
      DCptrR->GammaMin = DCptrR->GammaMin * K3;
      DCptrI->GammaMin = DCptrI->GammaMin * K3;
      DCptrI->Ntrf = DCptrI->Ntrf * KVi / KVr;
      DCptrR->P = -DCptrR->P / Sn;
      DCptrI->P = DCptrI->P / Sn;
      DCptrR->Q = -DCptrR->Q / Sn;
      DCptrI->Q = -DCptrI->Q / Sn;
      DCptrR->MVA = sqrt(DCptrR->P * DCptrR->P + DCptrR->Q * DCptrR->Q);
      DCptrI->MVA = sqrt(DCptrI->P * DCptrI->P + DCptrI->Q * DCptrI->Q);
    }
  }
}

/* --------------------------- InitializeLoad ------------------------------ */
void InitializeLoad(void)
/* Initialize load parameters for OH load model. */
{
  ACbusData *ACptr;
  VALUETYPE Pn, Qn, Pz, Qz, val;

  for (ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    if (ACptr->Pz == 0 && ACptr->a == 0)
      ACptr->Pn = ACptr->Pl;
    else {
      Pn = ACptr->Pn;
      Pz = 1 - Pn;
      val =
          ACptr->Pl / (Pn * pow(ACptr->V, ACptr->a) + Pz * ACptr->V * ACptr->V);
      ACptr->Pn = Pn * val;
      ACptr->Pz = Pz * val;
    }
    if (ACptr->Qz == 0 && ACptr->b == 0)
      ACptr->Qn = ACptr->Ql;
    else {
      Qn = ACptr->Qn;
      Qz = 1 - Qn;
      val =
          ACptr->Ql / (Qn * pow(ACptr->V, ACptr->b) + Qz * ACptr->V * ACptr->V);
      ACptr->Qn = Qn * val;
      ACptr->Qz = Qz * val;
    }
  }
}

/* --------------------------- Main Program  ------------------------------ */
int main(int argc, const char **argv)
/* Main Routine. */
{
  char *Name;
  bool flagv = false, flagD = false, PrintEvalue = false;
  VALUETYPE V, EigenValue;
  INDEX N1;
  FILE *Out;
  int i;

  SetArguments(argc, argv);

  if (HelpRequested()) {
    // display error message
    ErrorStop("");
    // exit from main
    return (1);
  }
  if (ExistParameter('l')) {
    Name = NameParameter('l');
    if (!NullName(Name))
      freopen(Name, "wt", stderr);
  }

  fprintf(stderr, "UW Continuation Power Flow (c)1992,1996,1999, 2006 C. "
                  "Canizares, F. Alvarado and S. Zhang.\n");
  if (ExistParameter('A') || ExistParameter('6'))
    Acont = false;
  else
    Acont = true;
  if (ExistParameter('P'))
    PQcont = false;
  else
    PQcont = true;
  if (ExistParameter('Q') && NullName(NameParameter('Q')))
    QRcont = false;
  else
    QRcont = true;
  if (ExistParameter('Q') && !strcmp(NameParameter('Q'), "X"))
    Xcont = false;
  else
    Xcont = true;
  if (ExistParameter('R'))
    Rcont = false;
  else
    Rcont = true;
  if (ExistParameter('N'))
    Acont = PQcont = QRcont = Rcont = Xcont = false;
  if (ExistParameter('a'))
    Tlim = false;
  else
    Tlim = true;
  if (ExistParameter('P') || ExistParameter('p'))
    PQlim = false;
  else
    PQlim = true;
  if (ExistParameter('q') && NullName(NameParameter('q')))
    Qlim = false;
  else
    Qlim = true;
  if (ExistParameter('q') && !strcmp(NameParameter('q'), "z"))
    Zlim = false;
  else
    Zlim = true;
  if (ExistParameter('q') && !strcmp(NameParameter('q'), "x"))
    Xlim = false;
  else
    Xlim = true;
  if (ExistParameter('r'))
    Vlim = false;
  else
    Vlim = true;
  if (ExistParameter('3')) {
    if (ExistParameter('4'))
      Elim = false;
    else
      Elim = true;
    if (ExistParameter('5'))
      Ilim = false;
    else
      Ilim = true;
  } else
    Elim = Ilim = false;
  flagPgMax = ExistParameter('X');
  flagSmax = ExistParameter('9');
  if (ExistParameter('n')) {
    Tlim = PQlim = Qlim = Vlim = Elim = Ilim = Xlim = Zlim = false;
    flagPgMax = flagSmax = true;
  }
  flagVloads = flagKdirection = false;

  // Call ReadData with the name of the inpute file
  ReadData(TrueParamStr(1));
  Tol = 1e-4;
  RealParameter('T', &Tol, 1e-10, 1e-3);
  tol = 0.0001;
  RealParameter('t', &tol, Tol, 0.1);
  alpha = 0.01;
  RealParameter('F', &alpha, 0.001, 1.);
  MaxIter = (INDEX)IntegerParameter('M', MaxIter, 0, 1000);
  if (ExistParameter('d')) {
    fprintf(stderr, "Tol %lf   tol %lf   alpha %lf\n", Tol, tol, alpha);
  }
  flagL = flagR = false;
  if ((ExistParameter('H') || ExistParameter('c')) &&
      (!NullName(NameParameter('K')) || flagKdirection))
    flagH = true;
  else
    flagH = false;
  if (ExistParameter('C') && (!NullName(NameParameter('K')) || flagKdirection))
    flagPoC = true;
  else
    flagPoC = false;
  if ((ExistParameter('D') && (!NullName(NameParameter('D')))) || flagVloads)
    flagD = true;
  else
    flagD = false;
  flagv = ReadInit();
  flagReducedContinuation = false;
  flagBS = ExistParameter('x');
  K1 = sqrt(2.0) * 3 / PI;
  K2 = 3 / PI;
  lambda = 0;
  DCinit();
  SVCinit();     /* FACTS */
  TCSCinit();    /* FACTS */
  STATCOMinit(); /* FACTS */
  Jacobian();
  if (flagH) {
    i = Homot();
    if (i < 0) {
      i = -i;
      fprintf(
          stderr,
          "\n *** The DC equations have a square root of a negative number.\n");
      fprintf(stderr, "     Try changing the DC controls.\n");
      fprintf(stderr, "Loading factor -> %-10.6lg\n\n", lambda);
      WriteSolution(--i, TrueParamStr(2), "DC problems:");
      exit(1);
    }
    fprintf(stderr, "**** Voltage Profile Case Solved **** \n\n");
    WriteSolution(--i, TrueParamStr(2), "U.E.P. Solution:");
    if (ExistParameter('y') && !NullName(NameParameter('y'))) {
      if (ExistParameter('d'))
        fprintf(stderr,
                "Write left e-vector for base case (see file 'evect.dat').\n");
      Evector(40, 0, 0.00001, false, &EigenValue);
      fprintf(stderr, "Minimum |e_value| -> %-10.6lg\n\n", EigenValue);
      PrintEvalue = true;
      Out = OpenOutput(NameParameter('y'));
      N1 = NacVar + 11 * Ndc / 2 + 3 * Nsvc + NtcscVar +
           7 * Nstatcom; /* FACTS */
      PrintLeftEvector(N1, Out);
    }
    if (ExistParameter('Y') && !NullName(NameParameter('Y'))) {
      if (ExistParameter('d'))
        fprintf(stderr,
                "Write right e-vector for base case (see file 'evect.dat').\n");
      Evector(40, 0, 0.00001, true, &EigenValue);
      if (!PrintEvalue)
        fprintf(stderr, "Minimum |e_value| -> %-10.6lg\n\n", EigenValue);
      PrintDirection('Y', x0, 1.0);
    }
  } else if (flagPoC) {
    i = PoCPoint();
    if (i < 0) {
      i = -i;
      fprintf(
          stderr,
          "\n *** The DC equations have a square root of a negative number.\n");
      fprintf(stderr, "     Try changing the DC controls.\n");
      fprintf(stderr, "Loading factor -> %-10.6lg\n\n", lambda);
      WriteSolution(--i, TrueParamStr(2), "DC problems:");
      exit(1);
    }
    fprintf(stderr, "**** Point of Collapse Case Solved ****   ");
    fprintf(stderr, "Loading factor -> %-10.6lg\n\n", lambda);
    WriteSolution(--i, TrueParamStr(2), "PoC Solution:");
  } else if (flagv) {
    i = 1;
    V = BlPtr->V;
    if (flagD) {
      i = Pflow(i, false, true, true);
      if (i < 0) {
        i = -i;
        fprintf(stderr, "\n *** The DC equations have a square root of a "
                        "negative number.\n");
        fprintf(stderr, "     Try changing the DC controls.\n");
        WriteSolution(--i, TrueParamStr(2), "DC problems:");
        exit(1);
      }
      fprintf(
          stderr,
          "**** Base Case Solved (to calculate OH load parameters) ****\n\n");
    } else
      InitializeLoad();
    RealParameter('L', &lambda, -1e6, 1e6);
    strcpy(BlPtr->Type, "BL");
    if (BlPtr->Area != nullptr && BlPtr->Area->Slack == BlPtr)
      strcat(BlPtr->Type, "A");
    Bl = BlPtr->N;
    BlPtr->V = V;
    BlPtr->Cont = nullptr;
    i = Pflow(i, flagD, true, false);
    if (i < 0) {
      i = -i;
      fprintf(
          stderr,
          "\n *** The DC equations have a square root of a negative number.\n");
      fprintf(stderr, "     Try changing the DC controls.\n");
      fprintf(stderr, "Loading factor -> %-10.6lg\n\n", lambda);
      WriteSolution(--i, TrueParamStr(2), "DC problems:");
      exit(1);
    }
    fprintf(stderr, "**** Voltage/Lambda Case Solved ****    ");
    fprintf(stderr, "Loading factor -> %-10.6lg\n\n", lambda);
    WriteSolution(--i, TrueParamStr(2), "Voltage/Lambda Solution:");
  } else if (ExistParameter('L') && (ExistParameter('K') || flagKdirection)) {
    i = 1;
    if (ExistParameter('b') || flagD) {
      i = Pflow(i, false, true, true);
      if (i < 0) {
        i = -i;
        fprintf(stderr, "\n *** The DC equations have a square root of a "
                        "negative number.\n");
        fprintf(stderr, "     Try changing the DC controls.\n");
        WriteSolution(--i, TrueParamStr(2), "DC problems:");
        exit(1);
      }
      fprintf(stderr,
              "**** Base Case Solved (to initialize power flow) ****\n\n");
    } else
      InitializeLoad();
    RealParameter('L', &lambda, -1e6, 1e6);
    if (lambda == 0 || (NullName(NameParameter('K')) && !flagKdirection)) {
      fprintf(stderr, "***Warning: The program has detected the -L option but "
                      "either lambda is zero\n");
      fprintf(stderr,
              "            or there is no gen./load variations defined.\n");
      if (!ExistParameter('b') && !flagD) {
        fprintf(stderr,
                "            The program will just solve the base case.\n");
        i = Pflow(i, false, true, true);
        if (i < 0) {
          i = -i;
          fprintf(stderr, "\n *** The DC equations have a square root of a "
                          "negative number.\n");
          fprintf(stderr, "     Try changing the DC controls.\n");
          WriteSolution(--i, TrueParamStr(2), "DC problems:");
          exit(1);
        }
        fprintf(stderr, "**** Base Case Solved ****\n\n");
      }
    } else {
      i = Pflow(i, flagD, true, false);
      if (i < 0) {
        i = -i;
        fprintf(stderr, "\n *** The DC equations have a square root of a "
                        "negative number.\n");
        fprintf(stderr, "     Try changing the DC controls.\n");
        fprintf(stderr, "Loading factor -> %-10.6lg\n\n", lambda);
        WriteSolution(--i, TrueParamStr(2), "DC problems:");
        exit(1);
      }
      fprintf(stderr, "**** Lambda Case Solved ****    ");
      fprintf(stderr, "Loading factor -> %-10.6lg\n\n", lambda);
    }
    WriteSolution(--i, TrueParamStr(2), "Lambda Solution:");
    if (ExistParameter('y') && !NullName(NameParameter('y'))) {
      N1 = NacVar + 11 * Ndc / 2 + 1 + 3 * Nsvc + NtcscVar +
           7 * Nstatcom; /* FACTS */
      x0 = new VALUETYPE[N1];
      if (ExistParameter('d'))
        fprintf(stderr,
                "Write left e-vector for base case (see file 'evect.dat').\n");
      Evector(40, 0, 0.00001, false, &EigenValue);
      fprintf(stderr, "Minimum |e_value| -> %-10.6lg\n\n", EigenValue);
      PrintEvalue = true;
      Out = OpenOutput(NameParameter('y'));
      PrintLeftEvector(N1, Out);
    }
    if (ExistParameter('Y') && !NullName(NameParameter('Y'))) {
      N1 = NacVar + 11 * Ndc / 2 + 1 + 3 * Nsvc + NtcscVar +
           7 * Nstatcom; /* FACTS */
      if (x0 == nullptr) {
        x0 = new VALUETYPE[N1];
      }
      if (ExistParameter('d'))
        fprintf(stderr,
                "Write right e-vector for base case (see file 'evect.dat').\n");
      Evector(40, 0, 0.00001, true, &EigenValue);
      if (!PrintEvalue)
        fprintf(stderr, "Minimum |e_value| -> %-10.6lg\n\n", EigenValue);
      PrintDirection('Y', x0, 1.0);
    }
  } else {
    i = Pflow(1, false, true, true);
    if (i < 0) {
      i = -i;
      fprintf(
          stderr,
          "\n *** The DC equations have a square root of a negative number.\n");
      fprintf(stderr, "     Try changing the DC controls.\n");
      WriteSolution(--i, TrueParamStr(2), "DC problems:");
      exit(1);
    }
    fprintf(stderr, "**** Base Case Solved ****\n\n");
    WriteSolution(--i, TrueParamStr(2), "Base Solution:");
    if (ExistParameter('y') && !NullName(NameParameter('y'))) {
      N1 = NacVar + 11 * Ndc / 2 + 1 + 3 * Nsvc + NtcscVar +
           7 * Nstatcom; /* FACTS */
      x0 = new VALUETYPE[N1];
      if (ExistParameter('d'))
        fprintf(stderr,
                "Write left e-vector for base case (see file 'evect.dat').\n");
      Evector(40, 0, 0.00001, false, &EigenValue);
      fprintf(stderr, "Minimum |e_value| -> %-10.6lg\n\n", EigenValue);
      PrintEvalue = true;
      Out = OpenOutput(NameParameter('y'));
      PrintLeftEvector(N1, Out);
    }
    if (ExistParameter('Y') && !NullName(NameParameter('Y'))) {
      N1 = NacVar + 11 * Ndc / 2 + 1 + 3 * Nsvc + NtcscVar +
           7 * Nstatcom; /* FACTS */
      if (x0 == nullptr) {
        x0 = new VALUETYPE[N1];
      }
      if (ExistParameter('d'))
        fprintf(stderr,
                "Write right e-vector for base case (see file 'evect.dat').\n");
      Evector(40, 0, 0.00001, true, &EigenValue);
      if (!PrintEvalue)
        fprintf(stderr, "Minimum |e_value| -> %-10.6lg\n\n", EigenValue);
      PrintDirection('Y', x0, 1.0);
    }
  }

  return (0);
}
