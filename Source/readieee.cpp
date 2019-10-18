/* Read AC data in IEEE common fotmat. */

#include "readdata.h"

/* --------- Global Input File --------- */
extern FILE *InputDataFile;

/* ---------------- ReadIEEE ----------------------------- */
void ReadIEEE()
/* Read Bus and Element data in WSCC format. */
{
  ACbusData *ACptr, *ACptrp, *ACptrs;
  ElementData *Eptr;
  ElementList *ELptr;
  AreaData *Aptr;
  char Line[75 + BUFLEN], Linep[BUFLEN];
  char Name[31], str[6];
  VALUETYPE KV, R, X, Sb;
  BOOLEAN flag = FALSE, card = FALSE, empty = FALSE, flagp = FALSE;
  INDEX i, j, k, NumBusDigits = 4, NumAreaDigits = 2, CircData = 21;

  Line[0] = '\0';
  for (i = 0; i <= 2; strcpy(dataPtr->Title[i], "\n"), i++)
    ;
  Sn = 100.0;
  RealParameter('$', &Sn, 1.0, 100000000.0);
  if (!strcmp(NameParameter('I'), "p"))
  {
    NumBusDigits = 5;
    NumAreaDigits = 3;
    CircData = 23;
  }
  for (;;)
  { /* Reading Loop */
    if (LineNum && Line[0] != '%')
      strcpy(Linep, Line);
    if (fgets(Line, BUFLEN, InputDataFile) == NULL)
    {
      fprintf(stderr, "***Warning: END OF DATA card missing in IEEE input data file.");
      break;
    }
    LineNum++;
    if (!strncmp(Line, "CARD", 4))
      card = TRUE;

    /* --------- Comment cards (Electrocon), Sn and Title ----------------- */
    else if (!strncmp(Line, "COMMENT DATA", 7))
      for (;;)
      {
        if (!flagp)
        {
          flagp = TRUE;
          GetStr(Linep, 39, 35, 35, dataPtr->Title[0]);
          Sb = GetValue(Linep, 32, 6, 0);
          if (Sb > 0)
            Sn = Sb;
        }
        if (fgets(Line, BUFLEN, InputDataFile) == NULL)
        {
          ErrorHalt("Missing $$$ card in COMMENT DATA.");
          break;
        }
        LineNum++;
        if (!strncmp(Line, "$$$", 3))
          break;
      }

    /* -------------------- AC bus data, Sn and Title --------------------- */
    else if (!strncmp(Line, "BUS DATA", 3))
      for (;;)
      {
        if (!flagp)
        {
          flagp = TRUE;
          GetStr(Linep, 46, 28, 28, dataPtr->Title[0]);
          Sb = GetValue(Linep, 32, 6, 0);
          if (Sb > 0)
            Sn = Sb;
        }
      BusData:
        if (fgets(Line, BUFLEN, InputDataFile) == NULL)
        {
          ErrorHalt("Missing -999 card in BUS DATA.");
          break;
        }
        LineNum++;
        if (Line[0] == '%')
          goto BusData;
        if (!strncmp(Line, "-999", 4))
          break;
        if (card)
        {
          if (fgets(Linep, BUFLEN, InputDataFile) == NULL)
          {
            ErrorHalt("Missing -999 card in BUS DATA.");
            break;
          }
          LineNum++;
          for (i = strlen(Line) - 1; i <= 74; Line[i] = ' ', i++)
            ;
          Line[75] = '\0';
          strcat(Line, Linep);
        }
        i = GetInt(Line, 25, 2);
        if (i != 4)
        {
          i = GetInt(Line, 1, NumBusDigits);
          GetStr(Line, 6, 12, 12, Name);
          KV = GetValue(Line, 77, 7, 2);
          for (j = 0; j <= 11; j++)
          {
            if (Name[j] != ' ')
              break;
            if (j == 11)
              empty = TRUE;
          }
          if (empty || (strpbrk(Name, "0") && strlen(Name) == 1))
          {
            empty = FALSE;
            strcpy(Name, "BUS_");
            GetStr(Line, 1, 5, 5, str);
            for (j = 0; j <= 3; j++)
              if (str[j] != ' ')
                break;
            strcat(Name, &str[j]);
            for (j = strlen(Name); j <= 11; Name[j] = ' ', j++)
              ;
            Name[12] = '\0';
          }
          ACptr = (ACbusData *)ACbusInList(i, Name, KV, Nac, 0);
          if (ACptr->N == 0)
          {
            Nac++;
            ACptr->Num = i;
            ACptr->N = Nac;
          }
          i = GetInt(Line, 19, 2);
          if (i)
          {
            ACptr->Area = (AreaData *)AreaInList(i, "", Narea);
            if (ACptr->Area->N == 0)
            {
              Narea++;
              ACptr->Area->N = i;
            }
          }
          GetStr(Line, 21, 3, 3, ACptr->Zone);
          i = GetInt(Line, 25, 2);
          ACptr->V = GetValue(Line, 28, 6, 4);
          if (ACptr->V <= 0)
            ACptr->V = 1;
          ACptr->Ang = GetValue(Line, 34, 7, 2) * K3;
          ACptr->Pl = GetValue(Line, 41, 9, 2) / Sn;
          ACptr->Ql = GetValue(Line, 50, 9, 2) / Sn;
          ACptr->Pg = GetValue(Line, 59, 9, 2) / Sn;
          ACptr->Qg = GetValue(Line, 68, 8, 2) / Sn;
          KV = GetValue(Line, 85, 6, 4);
          ACptr->G = GetValue(Line, 107, 8, 4);
          ACptr->B = GetValue(Line, 115, 8, 4);
          if (i == 0 && strcmp(ACptr->Type, "BC"))
            ACptr->Cont = ACptr;
          else if (i == 1)
          {
            strcpy(ACptr->Type, "BV");
            strcpy(ACptr->cont, "Q");
            ACptr->VCont = ACptr->Qg;
            ACptr->Vmax = GetValue(Line, 91, 8, 2);
            if (ACptr->Vmax <= 0)
              ACptr->Vmax = 10.;
            ACptr->Vmin = GetValue(Line, 99, 8, 2);
            if (ACptr->Vmin <= 0)
              ACptr->Vmin = 0.001;
            if (ACptr->Vmax <= ACptr->Vmin)
            {
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              ErrorHalt("AC bus V limits are wrong: Vmin >= Vmax.");
            }
            ACptr->Cont = ACptr;
            ACptr->Qmax = 99999999.;
            ACptr->Qmin = -99999999.;
          }
          else if (i == 2 || i == 3)
          {
            if (ExistParameter('g'))
              ACptr->Qg = 0;
            ACptr->Qmax = GetValue(Line, 91, 8, 2) / Sn;
            ACptr->Qmin = GetValue(Line, 99, 8, 2) / Sn;
            if (ACptr->Qmax < ACptr->Qmin)
            {
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              ErrorHalt("AC bus Q limits are wrong: Qmin > Qmax.");
            }
            /*if (ACptr->Qmax==ACptr->Qmin && ACptr->Qmax==0.) {
            ACptr->Qmax=9999999.;
            ACptr->Qmin= -999999.;
          }*/
            if (ACptr->Qmax == ACptr->Qmin)
              ACptr->Cont = ACptr;
            else
            {
              Nvolt++;
              ACptr->val = GetValue(Line, 85, 6, 4);
              if (ACptr->val <= 0)
                ACptr->val = 1;
              j = GetInt(Line, 124, NumBusDigits);
              if (j > 0 && j != ACptr->Num)
              {
                ACptr->Nc = j;
                strcpy(ACptr->Type, "BG");
                strcpy(ACptr->cont, "V");
                ACptr->Kbg = 1;
                ACptr->VCont = ACptr->val;
              }
              else
              {
                strcpy(ACptr->Type, "BQ");
                strcpy(ACptr->cont, "V");
                ACptr->VCont = ACptr->V = ACptr->val;
              }
            }
            if (i == 3)
            {
              Nslack++;
              strcat(ACptr->Type, "S");
              if (ACptr->Area != NULL)
              {
                ACptr->Area->i++;
                if (ACptr->Area->i > 1)
                {
                  if (card)
                    LineNum--;
                  fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
                  ErrorHalt("The Area has more than one slack bus.");
                }
              }
            }
          }
        }
      }

    /* -------------------- AC elemet data ------------------------- */
    else if (!strncmp(Line, "BRANCH DATA", 6))
      for (;;)
      {
      BranchData:
        if (fgets(Line, BUFLEN, InputDataFile) == NULL)
        {
          ErrorHalt("Missing -999 card in BRANCH DATA.");
          break;
        }
        LineNum++;
        if (Line[0] == '%')
          goto BranchData;
        if (!strcmp(NameParameter('I'), "p") && Line[strlen(Line) - 2] == '0')
          goto BranchData;
        if (!strncmp(Line, "-999", 4))
          break;
        k = GetInt(Line, 19, 1);
        if (card && k)
        {
          if (fgets(Linep, BUFLEN, InputDataFile) == NULL)
          {
            ErrorHalt("Missing -999 card in BRANCH DATA.");
            break;
          }
          LineNum++;
          for (i = strlen(Line) - 1; i <= 74; Line[i] = ' ', i++)
            ;
          Line[75] = '\0';
          strcat(Line, Linep);
        }
        i = GetInt(Line, 1, NumBusDigits);
        ACptr = (ACbusData *)ACbusInList(i, "", 0., Nac, 1);
        if (ACptr->N == 0)
        {
          if (card)
            LineNum--;
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("The tap bus has not been defined (check BUS DATA cards).");
        }
        j = GetInt(Line, 6, NumBusDigits);
        ACptrp = (ACbusData *)ACbusInList(j, "", 0., Nac, 1);
        if (ACptrp->N == 0)
        {
          if (card)
            LineNum--;
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("The Z bus has not been defined (check BUS DATA cards).");
        }
        Eptr = (ElementData *)ElemInList(ACptr, ACptrp, NacEl++, 0, "", "");
        R = GetValue(Line, 20, 10, 6);
        X = GetValue(Line, 30, 10, 6);
        if (fabs(R) < 0.0000001 && fabs(X) < 0.0000001)
        {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("AC element is a short circuit. Try eliminating it.");
        }
        else
        {
          Eptr->G = R / (R * R + X * X);
          Eptr->B = -X / (R * R + X * X);
        }
        Eptr->B1 = Eptr->B2 = GetValue(Line, 41, 9, 5) / 2;
        Eptr->Imax = GetInt(Line, 51, 5) / Sn;
        i = GetInt(Line, 11, 2);
        if (i)
        {
          Eptr->Area = (AreaData *)AreaInList(i, "", Narea);
          if (Eptr->Area->N == 0)
          {
            if (card)
              LineNum--;
            fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
            ErrorHalt("The Area has not been defined (check BUS DATA cards).");
          }
        }
        GetStr(Line, 13, 3, 3, Eptr->Zone);
        GetStr(Line, 17, 1, 1, Eptr->Ckt);
        if (k == 0)
          strcpy(Eptr->Type, "L");
        else if (k >= 1)
        {
          Eptr->Tap = 1 / GetValue(Line, 77, 6, 4);
          Eptr->Ang = GetValue(Line, 84, 7, 2) * K3;
          j = GetInt(Line, 69, NumBusDigits);
          if (j > 0)
          {
            ACptrs = (ACbusData *)ACbusInList(j, "", 0., Nac, 1);
            if (ACptrs->N == 0)
            {
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              ErrorHalt("The controlled bus has not been defined (check BUS DATA cards).");
            }
          }
          else
            ACptrs = ACptr;
          if (k == 1)
            strcpy(Eptr->Type, "T");
          else if (k == 2)
          {
            Eptr->Tmin = GetValue(Line, 91, 7, 2);
            Eptr->Tmax = GetValue(Line, 98, 7, 2);
            if (Eptr->Tmax < 0)
              Eptr->Tmax = 1.1;
            if (Eptr->Tmin < 0)
              Eptr->Tmin = 0.9;
            if (Eptr->Tmax < Eptr->Tmin)
            {
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              ErrorHalt("LTC limits are wrong: Tmin > Tmax.");
            }
            if (Eptr->Tmax == Eptr->Tmin)
            {
              strcpy(Eptr->Type, "T");
              Eptr->Tmax = Eptr->Tmin = 0;
            }
            else if (strcmp(ACptrs->Type, "BG") && strcmp(ACptrs->Type, "BQ"))
            {
              NregV++;
              ACptrs->Reg = (ElementList *)AddElemToList(ACptrs->Reg, Eptr);
              Eptr->Cont = ACptrs;
              Eptr->Min = GetValue(Line, 113, 7, 5);
              Eptr->Max = GetValue(Line, 120, 7, 5);
              if (Eptr->Max < 0)
                Eptr->Max = 1.;
              if (Eptr->Min < 0)
                Eptr->Min = 1.;
              if (Eptr->Max < Eptr->Min)
              {
                fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
                ErrorHalt("LTC limits are wrong: Vmin > Vmax.");
              }
              if (Eptr->Max == Eptr->Min)
              {
                strcpy(Eptr->Type, "R");
                strcpy(ACptrs->Type, "BT");
                ACptrs->V = Eptr->Max;
                if (ACptrs->V == 0)
                  ACptrs->V = 1.;
                Eptr->Max = Eptr->Min = 0;
              }
              else
              {
                strcpy(Eptr->Type, "RV");
                if (strcmp(ACptrs->Type, "BT"))
                  strcpy(ACptrs->Type, "BR");
                if (ACptrs->Vmax > Eptr->Max || ACptrs->Vmax == 0)
                  ACptrs->Vmax = Eptr->Max;
                if (ACptrs->Vmin > Eptr->Min || ACptrs->Vmin == 0)
                {
                  ACptrs->Vmin = Eptr->Min;
                  if (ACptrs->Vmin == 0)
                    ACptr->Vmin = 0.00001;
                }
              }
            }
            else
            {
              strcpy(Eptr->Type, "T");
              Eptr->Tmax = Eptr->Tmin = 0;
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              fprintf(stderr, "***Warning: The controlled bus of this LTC is not PQ (check BUS DATA cards).\n");
              fprintf(stderr, "            It will be treated as a fixed tap transformer.\n");
            }
          }
          else if (k == 3)
          {
            Eptr->Tmin = GetValue(Line, 91, 7, 2);
            Eptr->Tmax = GetValue(Line, 98, 7, 2);
            if (Eptr->Tmax < 0)
              Eptr->Tmax = 1.1;
            if (Eptr->Tmin < 0)
              Eptr->Tmin = 0.9;
            if (Eptr->Tmax < Eptr->Tmin)
            {
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              ErrorHalt("LTC limits are wrong: Tmin > Tmax.");
            }
            if (Eptr->Tmax == Eptr->Tmin)
            {
              strcpy(Eptr->Type, "T");
              Eptr->Tmin = Eptr->Tmax = 0;
            }
            else
            {
              NregPQ++;
              if (ACptrs != ACptr && ACptrs != ACptrp)
                ACptrs = ACptr;
              ACptrs->Reg = (ElementList *)AddElemToList(ACptrs->Reg, Eptr);
              Eptr->Cont = ACptrs;
              Eptr->Min = GetValue(Line, 113, 7, 5) / Sn;
              Eptr->Max = GetValue(Line, 120, 7, 5) / Sn;
              if (Eptr->Max < Eptr->Min)
              {
                fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
                ErrorHalt("LTC limits are wrong: Qmin > Qmax.");
              }
              if (Eptr->Max == Eptr->Min)
              {
                strcpy(Eptr->Type, "RQ");
                Eptr->Cvar = Eptr->Max;
                Eptr->Max = Eptr->Min = 0;
              }
              else
                strcpy(Eptr->Type, "RN");
              Eptr->Ncont = ACptrs->Ncont;
              ACptrs->Ncont++;
            }
          }
          else if (k == 4)
          {
            Eptr->Tmin = GetValue(Line, 91, 7, 2) * K3;
            Eptr->Tmax = GetValue(Line, 98, 7, 2) * K3;
            if (Eptr->Tmax < Eptr->Tmin)
            {
              fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
              ErrorHalt("LTC limits are wrong: Tmin > Tmax.");
            }
            if (Eptr->Tmax == Eptr->Tmin)
            {
              strcpy(Eptr->Type, "T");
              Eptr->Tmin = Eptr->Tmax = 0;
            }
            else
            {
              NregPQ++;
              if (ACptrs != ACptr && ACptrs != ACptrp)
                ACptrs = ACptr;
              ACptrs->Reg = (ElementList *)AddElemToList(ACptrs->Reg, Eptr);
              Eptr->Cont = ACptrs;
              Eptr->Min = GetValue(Line, 113, 7, 5) / Sn;
              Eptr->Max = GetValue(Line, 120, 7, 5) / Sn;
              if (Eptr->Max < Eptr->Min)
              {
                fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
                ErrorHalt("LTC limits are wrong: Pmin > Pmax.");
              }
              if (Eptr->Max == Eptr->Min)
              {
                strcpy(Eptr->Type, "RP");
                Eptr->Cvar = Eptr->Max;
                Eptr->Max = Eptr->Min = 0;
              }
              else
                strcpy(Eptr->Type, "RM");
              Eptr->Ncont = ACptrs->Ncont;
              ACptrs->Ncont++;
            }
          }
        }
      }

    /* -------------------- Area data -------------------------------- */
    else if (!strncmp(Line, "INTERCHANGE DATA", 11))
      for (;;)
      {
      AreaIntData:
        if (fgets(Line, BUFLEN, InputDataFile) == NULL)
        {
          ErrorHalt("Missing -9 card in INTERCHANGE DATA.");
          break;
        }
        LineNum++;
        if (Line[0] == '%')
          goto AreaIntData;
        if (!strncmp(Line, "-9", 2))
          break;
        flag = TRUE;
        i = GetInt(Line, 1, NumAreaDigits);
        Aptr = (AreaData *)AreaInList(i, "", Narea);
        if (Aptr->N == 0)
        {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("The Area does not contain any buses (check BUS DATA cards).");
        }
        GetStr(Line, 46, 30, 30, Name);
        for (j = 0; j <= 29; j++)
        {
          if (Name[j] != ' ')
            break;
          if (j == 29)
            empty = TRUE;
        }
        if (empty || (strpbrk(Name, "0") && strlen(Name) == 1))
        {
          empty = FALSE;
          strcpy(Name, "AREA_");
          GetStr(Line, 1, 4, 4, str);
          for (j = 0; j <= 3; j++)
            if (str[j] != ' ')
              break;
          strcat(Name, &str[j]);
          for (j = strlen(Name); j <= 29; Name[j] = ' ', j++)
            ;
          Name[30] = '\0';
        }
        strcpy(Aptr->Name, Name);
        i = GetInt(Line, 4, NumBusDigits);
        ACptr = (ACbusData *)ACbusInList(i, "", 0., Nac, 1);
        if (ACptr->N == 0)
        {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("The area slack bus has not been defined (check BUS DATA cards).");
        }
        if (ACptr->Area != Aptr)
        {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("The slack bus is not in the area (check BUS DATA cards).");
        }
        Aptr->Slack = Aptr->BSptr = ACptr;
        Aptr->P = GetValue(Line, 21, 8, 1) / Sn;
      }

    /* -------------------- Tie Lines data -------------------------------- */
    else if (!strncmp(Line, "TIE LINE DATA", 3))
      for (;;)
      {
      TieLineData:
        if (fgets(Line, BUFLEN, InputDataFile) == NULL)
        {
          ErrorHalt("Missing -999 card in TIE LINE DATA.");
          break;
        }
        LineNum++;
        if (Line[0] == '%')
          goto TieLineData;
        if (!strncmp(Line, "-999", 4))
          break;
        i = GetInt(Line, 1, NumBusDigits);
        ACptr = (ACbusData *)ACbusInList(i, "", 0., Nac, 1);
        if (ACptr->N == 0)
        {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("The metered bus has not been defined (check BUS DATA cards).");
        }
        i = GetInt(Line, 11, NumBusDigits);
        ACptrp = (ACbusData *)ACbusInList(i, "", 0., Nac, 1);
        if (ACptrp->N == 0)
        {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          ErrorHalt("The nonmetered bus has not been defined (check BUS DATA cards).");
        }
        GetStr(Line, CircData, 1, 1, str);
        for (ELptr = ACptr->Elem; ELptr != NULL; ELptr = ELptr->Next)
        {
          Eptr = ELptr->Eptr;
          if (((Eptr->From == ACptr && Eptr->To == ACptrp) ||
               (Eptr->From == ACptrp && Eptr->To == ACptr)) &&
              !strcmp(Eptr->Ckt, str))
            break;
        }
        if (ELptr != NULL)
        {
          if (Eptr->Meter == NULL)
            Eptr->Meter = ACptr;
          else
          {
            fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
            fprintf(stderr, "***Warning: This tie line was already defined (check circuit code).\n");
          }
        }
        else
        {
          fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
          fprintf(stderr, "***Warning: This element has not been defined (check BRANCH DATA cards).\n");
        }
      }

    /* -------------------- DC data -------------------------------- */
    else if (!strncmp(Line, "BD ", 3) || !strncmp(Line, "BZ ", 3) ||
             !strncmp(Line, "LD ", 3))
      ReadEPRIdc(Line);

    /* FACTS */

    /* ---------------------- SVC data ------------------------ */
    else if (!strncmp(Line, "FS ", 3))
      ReadSVC(Line);

    /* ---------------------- TCSC data ------------------------ */
    else if (!strncmp(Line, "FC ", 3))
      ReadTCSC(Line);

    /* ---------------------- STATCOM data ------------------------ */
    else if (!strncmp(Line, "FT ", 3))
      ReadSTATCOM(Line);

    /* END FACTS */

    else if (!strncmp(Line, "END OF DATA", 11))
      break;
    else if (flagp && Line[0] != 'C')
    {
      fprintf(stderr, "Input Line-> %d\n%s", LineNum, Line);
      fprintf(stderr, "***Warning: The program will ignore this line.\n");
    }
  }
  fclose(InputDataFile);

  if (!flag)
  {
    Narea = 0;
    for (ACptr = dataPtr->ACbus; ACptr != NULL; ACptr->Area = NULL, ACptr = ACptr->Next)
      ;
  }
  MaxIter = 50;
  for (ACptr = dataPtr->ACbus; ACptr != NULL; ACptr = ACptr->Next)
    if (strpbrk(ACptr->Type, "G"))
    {
      ACptrp = (ACbusData *)ACbusInList(ACptr->Nc, "", 0., Nac, 1);
      if (ACptrp->N == 0)
      {
        fprintf(stderr, "Error: The voltage controlled bus %d %s\n", ACptrp->Num, ACptrp->Name);
        fprintf(stderr, "       has not been defined.  Check AC bus data.\n");
        /*WriteSummary();
         exit(ERROREXIT);*/
        InputError = TRUE;
      }
      if (!strcmp(ACptrp->Type, "B"))
        strcpy(ACptrp->Type, "BC");
      else if (strcmp(ACptrp->Type, "BC"))
      {
        fprintf(stderr, "Error: The voltage controlled bus %d %s\n", ACptrp->Num, ACptrp->Name);
        fprintf(stderr, "       is not a PQ bus.  Check AC bus data.\n");
        /*WriteSummary();
         exit(ERROREXIT);*/
        InputError = TRUE;
      }
      ACptr->Cont = ACptrp;
      ACptrp->VCont = ACptrp->V = ACptr->VCont;
      ACptrp->Cont = NULL;
      ACptrp->Kbg++;
      if (flag2Vcontrol)
      {
        ACptrp->Kbg1 = ACptrp->Kbg1 + ACptr->Qmax;
        ACptrp->Kbg2 = ACptrp->Kbg2 + ACptr->Qmin;
      }
      if (ACptrp->Kbg > 1)
      {
        if (!flag2Vcontrol)
        {
          fprintf(stderr, "***Warning: The IEEE common format assumes that multiple generators\n");
          fprintf(stderr, "            controlling voltage at bus %d %s\n", ACptrp->Num, ACptrp->Name);
          fprintf(stderr, "            share equally the control.\n");
        }
        fprintf(stderr, "***Warning: The program will use the first generator controlling bus \n");
        fprintf(stderr, "            to define the remote voltage at bus %d %s\n", ACptrp->Num, ACptrp->Name);
      }
    }
}
