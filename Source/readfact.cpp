/* Read SVC data in WSCC type format  */

#include "readdata.h"

/* --------- Global Input File --------- */
extern FILE *InputDataFile;

/* --------------------- Read SVC ----------------------------- */
void ReadSVC(const char *Line)
{
  SVCbusData *SVCptr;
  SVClist *ptr;
  ACbusData *ACptr, *ptrac, *ptrac1;
  ElementData *Eptr;
  char Name[31], str[31], Name1[31], Name2[31], Name3[31];
  VALUETYPE KV, KV1, KV2, Xth_l, X, Taps;
  VALUETYPE G, B, G1, B1, Tap, Ang, G2, B2, Ssvc;
  INDEX Sec;

  /* ---------------------- SVC data ------------------------ */
  if (!strncmp(Line, "FS ", 3))
  {
    GetStr(Line, 7, 12, 12, Name);
    KV = GetValue(Line, 15, 4, 0);
    GetStr(Line, 20, 12, 12, Name1);
    KV1 = GetValue(Line, 28, 4, 0);
    ptrac1 = ACbusInList(0, Name1, KV1, Nac, 1);
    KV2 = GetValue(Line, 64, 4, 0);
    Xth_l = GetValue(Line, 68, 8, 7);
    strcpy(Name2, Name);
    if (Xth_l != 0.0)
    {
      Name2[7] = 'L';
      Name2[8] = '\0';
      GetStr(Line, 64, 4, 4, Name3);
      strcat(Name2, Name3);
      ptrac = ACbusInList(0, Name2, KV2, Nac, 1);
      if (ptrac->N == 0)
      {
        Nac++;
        ptrac->Num = ptrac->N = Nac;
      }
      ACptr = ACbusInList(0, Name, KV, Nac, 1);
      if (ACptr->N == 0)
      {
        Nac++;
        ACptr->Num = ACptr->N = Nac;
      }
      strcpy(ptrac->Type, "B");
      strcpy(ptrac->Zone, ACptr->Zone);
      strcpy(ptrac->Owner, ACptr->Owner);
      ptrac->V = 1.00;
      /*  Create new element between buses Name2 & Name */
      if (KV == 0)
      {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        ErrorHalt("Base voltage at bus 1 is zero.");
      }
      if (KV2 == 0)
      {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        ErrorHalt("Base voltage at bus 2 is zero.");
      }
      Ssvc = GetValue(Line, 56, 4, 0);
      if (Ssvc == 0.0)
        Ssvc = Sn;
      X = Xth_l * Sn / Ssvc;
      if (fabs(X) < 0.0001)
      {
        fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
        ErrorHalt("AC element is a short circuit. Try eliminating it.");
        B = 0;
      }
      else
      {
        B = -1.0 / X;
      }
      G = 0.0;
      G1 = 0.0;
      B1 = 0.0;
      Tap = Taps = 1;
      Ang = 0;
      Sec = 0;
      G2 = G1;
      B2 = B1;
      strcpy(str, " ");
      Eptr = ElemInList(ACptr, ptrac, NacEl, 0, "", str);
      Eptr->Sec = Sec;
      Eptr->G = G;
      Eptr->B = B;
      Eptr->G1 = G1;
      Eptr->B1 = B1;
      Eptr->G2 = G2;
      Eptr->B2 = B2;
      Eptr->Tap = Tap;
      Eptr->Taps = Taps;
      Eptr->Ang = Ang;
      strcpy(Eptr->Owner, ACptr->Owner);
      NacEl++;
      strcpy(Eptr->Zone, " ");
      strcpy(Eptr->Type, "L");
      ptrac->Area = ACptr->Area;
      ptrac->Cont = ptrac;
    }
    else
      ptrac = ACbusInList(0, Name, KV, Nac, 1);
    SVCptr = SVCbusInList(Name2, 0, ptrac, ptrac1);
    if (SVCptr->Xc > 0)
    {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      ErrorHalt("The SVC bus was previously defined (check FS cards).");
    }
    if (SVCptr->N == 0)
    {
      Nsvc++;
      SVCptr->N = Nsvc;
    }
    ptr = ptrac->SVC;
#ifdef WINDOWS
    ptrac->SVC = new SVClist;
#else
    ptrac->SVC = (SVClist *)malloc(sizeof(SVClist));
    if (ptrac->SVC == nullptr)
    {
      fclose(InputDataFile);
      ErrorHalt("Insufficient memory to allocate SVC element data.");
      exit(ERROREXIT);
    }
#endif
    ptrac->SVC->SVC = SVCptr;
    ptrac->SVC->Next = ptr;
    SVCptr->Xth_l = Xth_l;
    if (KV2 == 0.0)
      KV2 = KV;
    SVCptr->Vsvc = KV2;
    SVCptr->Xc = GetValue(Line, 32, 8, 7);
    SVCptr->Xl = GetValue(Line, 40, 8, 7);
    SVCptr->AlphaMin = GetValue(Line, 48, 3, 0);
    if (SVCptr->AlphaMin == 0)
      SVCptr->AlphaMin = 90;
    SVCptr->AlphaMax = GetValue(Line, 51, 3, 0);
    if (SVCptr->AlphaMax == 0)
      SVCptr->AlphaMax = 178;
    SVCptr->slope = GetValue(Line, 54, 2, 0);
    SVCptr->SVC_base = GetValue(Line, 56, 4, 0);
    SVCptr->Vref = GetValue(Line, 60, 4, 3);
    if (SVCptr->Vref == 0)
      SVCptr->Vref = 1.0;
    SVCptr->alpha_svc = GetValue(Line, 76, 5, 1);
  }
}

/* --------------------- Read TCSC ----------------------------- */
#ifdef ANSIPROTO
void ReadTCSC(const char *Line)
#else
void ReadTCSC(Line) char *Line;
#endif
{
  TCSCbusData *TCSCptr;
  TCSClist *ptr;
  ACbusData *ptrac, *ptrac1;
  char Name[31], Name1[31];
  VALUETYPE KV, KV1;

  /* ---------------------- TCSC data ------------------------ */
  if (!strncmp(Line, "FC ", 3))
  {
    GetStr(Line, 7, 12, 12, Name);
    KV = GetValue(Line, 15, 4, 0);
    GetStr(Line, 20, 12, 12, Name1);
    KV1 = GetValue(Line, 28, 4, 0);
    ptrac = ACbusInList(0, Name, KV, Nac, 1);
    if (ptrac->N == 0)
    {
      Nac++;
      ptrac->Num = ptrac->N = Nac;
    }
    ptrac1 = ACbusInList(0, Name1, KV1, Nac, 1);
    if (ptrac1->N == 0)
    {
      Nac++;
      ptrac1->Num = ptrac1->N = Nac;
    }
    if (KV == 0)
    {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      ErrorHalt("Base voltage at bus From is zero.");
    }
    if (KV1 == 0)
    {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      ErrorHalt("Base voltage at bus To is zero.");
    }
    if (ptrac == ptrac1)
    {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      ErrorHalt("Both TCSC buses are the same.");
    }
    TCSCptr = TCSCbusInList(Name, 0, ptrac, ptrac1);
    if (TCSCptr->Xc > 0)
    {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      ErrorHalt("The TCSC bus was previously defined (check FC cards).");
    }
    if (TCSCptr->N == 0)
    {
      Ntcsc++;
      TCSCptr->N = Ntcsc;
    }
    ptr = ptrac->TCSC;
#ifdef WINDOWS
    ptrac->TCSC = new TCSClist;
#else
    ptrac->TCSC = (TCSClist *)malloc(sizeof(TCSClist));
    if (ptrac->TCSC == nullptr)
    {
      fclose(InputDataFile);
      ErrorHalt("Insufficient memory to allocate TCSC element data.");
      exit(ERROREXIT);
    }
#endif
    ptrac->TCSC->TCSC = TCSCptr;
    ptrac->TCSC->Next = ptr;
    TCSCptr->From = ptrac;
    ptr = ptrac1->TCSC;
#ifdef WINDOWS
    ptrac1->TCSC = new TCSClist;
#else
    ptrac1->TCSC = (TCSClist *)malloc(sizeof(TCSClist));
    if (ptrac1->TCSC == nullptr)
    {
      fclose(InputDataFile);
      ErrorHalt("Insufficient memory to allocate TCSC element data.");
      exit(ERROREXIT);
    }
#endif
    ptrac1->TCSC->TCSC = TCSCptr;
    ptrac1->TCSC->Next = ptr;
    TCSCptr->To = ptrac1;
    TCSCptr->Xc = GetValue(Line, 32, 8, 7);
    TCSCptr->Xl = GetValue(Line, 40, 8, 7);
    TCSCptr->AlphaMin = GetValue(Line, 48, 3, 0);
    if (TCSCptr->AlphaMin <= 0)
      TCSCptr->AlphaMin = 143;
    TCSCptr->AlphaMax = GetValue(Line, 51, 3, 0);
    if (TCSCptr->AlphaMax <= 0 || TCSCptr->AlphaMax > 175)
      TCSCptr->AlphaMax = 175;
    TCSCptr->Control = GetValue(Line, 54, 8, 3);
    if (Line[61] == 'P' || Line[61] == 'I' || Line[61] == 'D')
      TCSCptr->Cont[0] = Line[61];
    else
      TCSCptr->Cont[0] = 'X';
    TCSCptr->Cont[1] = '\0';
    NtcscVar = 7 * Ntcsc;
    TCSCptr->TCSC_base = GetValue(Line, 63, 4, 0);
    TCSCptr->alpha_tcsc = GetValue(Line, 67, 5, 1) * K3;
    TCSCptr->Max = GetValue(Line, 72, 3, 0);
  }
}

/* --------------------- Read STATCOM ----------------------------- */
#ifdef ANSIPROTO
void ReadSTATCOM(const char *Line)
#else
void ReadSTATCOM(Line) char *Line;
#endif
{
  STATCOMbusData *STATCOMptr;
  STATCOMlist *ptr;
  ACbusData *ptrac, *ptrac1;
  char Name[31], Name1[31], Cont[2];
  VALUETYPE KV, KV1, R, X, Gc, MVA, val;

  /* ---------------------- STATCOM data ------------------------ */
  if (!strncmp(Line, "FT ", 3))
  {
    GetStr(Line, 7, 12, 12, Name);
    KV = GetValue(Line, 15, 4, 0);
    ptrac = ACbusInList(0, Name, KV, Nac, 1);
    GetStr(Line, 20, 12, 12, Name1);
    KV1 = GetValue(Line, 28, 4, 0);
    ptrac1 = ACbusInList(0, Name1, KV1, Nac, 1);
    STATCOMptr = STATCOMbusInList(Name, 0, ptrac, ptrac1);
    if (STATCOMptr->B > 0)
    {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      ErrorHalt("The STATCOM bus was previously defined (check FT cards).");
    }
    if (STATCOMptr->N == 0)
    {
      Nstatcom++;
      STATCOMptr->N = Nstatcom;
    }
    ptr = ptrac->STATCOM;
#ifdef WINDOWS
    ptrac->STATCOM = new STATCOMlist;
#else
    ptrac->STATCOM = (STATCOMlist *)malloc(sizeof(STATCOMlist));
    if (ptrac->STATCOM == nullptr)
    {
      fclose(InputDataFile);
      ErrorHalt("Insufficient memory to allocate STATCOM element data.");
      exit(ERROREXIT);
    }
#endif
    ptrac->STATCOM->STATCOM = STATCOMptr;
    ptrac->STATCOM->Next = ptr;
    R = GetValue(Line, 32, 6, 5);
    if (R < 0)
    {
      fprintf(stderr, "***Warning: The STATCOM ac resistance R is negative. The corresponding \n");
      fprintf(stderr, "            positive value will be used.");
      R = -R;
    }
    X = GetValue(Line, 38, 6, 5);
    if (X <= 0)
    {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      ErrorHalt("The STATCOM reactance must be greater than zero.");
    }
    Gc = GetValue(Line, 44, 6, 5);
    if (Gc < 0)
    {
      fprintf(stderr, "***Warning: The STATCOM dc conductance Gc is negative. The corresponding \n");
      fprintf(stderr, "            positive value will be used.");
      Gc = -Gc;
    }
    MVA = GetValue(Line, 65, 5, 0);
    if (MVA <= 0)
      MVA = Sn;
    STATCOMptr->MVA = MVA;
    STATCOMptr->R = R;
    if (R > 0 || X > 0)
    {
      STATCOMptr->G = R / (R * R + X * X);
      STATCOMptr->B = -X / (R * R + X * X);
    }
    STATCOMptr->Gc = Gc;
    STATCOMptr->Imin = fabs(GetValue(Line, 50, 6, 4));
    STATCOMptr->Imax = fabs(GetValue(Line, 56, 6, 4));
    if (STATCOMptr->Imin == 0)
      STATCOMptr->Imin = 99999999.;
    if (STATCOMptr->Imax == 0)
      STATCOMptr->Imax = 99999999.;
    STATCOMptr->slope = GetValue(Line, 62, 3, 0);
    STATCOMptr->Vref = GetValue(Line, 70, 4, 3);
    if (STATCOMptr->Vref <= 0)
      STATCOMptr->Vref = 1.0;
    val = GetValue(Line, 75, 6, 5);
    GetStr(Line, 74, 1, 1, Cont);
    if (!strcmp(Cont, "P"))
    {
      strcpy(STATCOMptr->Cont, "AL");
      strcpy(STATCOMptr->Cont1, "AL");
      if (val <= 0)
        STATCOMptr->Contref = 0.9;
      else
        STATCOMptr->Contref = val;
    }
    else
    {
      strcpy(STATCOMptr->Cont, "PW");
      strcpy(STATCOMptr->Cont1, "PW");
      if (val <= 0)
        STATCOMptr->Contref = 1;
      else
        STATCOMptr->Contref = val;
    }
  }
}
