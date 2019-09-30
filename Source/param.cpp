/* param.c 030390 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#ifdef ANSIPROTO
#include <float.h>
#endif

#include "constant.h"
#include "param.h"
#include "sparse.h"

/* ==================== Global Variables ============================ */
int GlobalArgc;
const char **GlobalArgv;
#ifdef WINDOWS
extern CString dir;
#endif

/* ========================== SetArguments ========================== */
#ifdef ANSIPROTO
void SetArguments(int argc, const char **argv)
#else
void SetArguments(argc,argv)
int argc;
char **argv;
#endif
{
  GlobalArgc = argc;
  GlobalArgv = argv;
}

/* ========================== ExistParameters ======================= */
#ifdef ANSIPROTO
BOOLEAN ExistParameter(char ch)
#else
BOOLEAN ExistParameter(ch)
char ch;
#endif
{
  int i;
  /* BEGIN */
  for (i=1; i<GlobalArgc; i++) {
    if (GlobalArgv[i][0] == '-') {
      if (GlobalArgv[i][1] == ch) return(TRUE);
    }
  }
  return(FALSE);
}

/* ============================ TrueParamStr ========================== */
#ifdef ANSIPROTO
char * TrueParamStr(int Item)
#else
char * TrueParamStr(Item)
int Item;
#endif
{
  int i,k;
  char * TempPtr;
  /* BEGIN */
  k = 0;
  for (i=1; i<GlobalArgc; i++) {
    if (GlobalArgv[i][0] != '-') k++;
    if (k >= Item) {
#ifdef WINDOWS
      TempPtr = new char[strlen(GlobalArgv[i])+1];
#else
      TempPtr = (char *) malloc((strlen(GlobalArgv[i])+1)*sizeof(char));
#endif
      if (GlobalArgv[i][0]!='>') strcpy(TempPtr,GlobalArgv[i]);
      else strcpy(TempPtr,&GlobalArgv[i][1]);
      return(TempPtr);
    }
  }
  return(NULL);
}

/* =========================== HelpRequested ======================= */
BOOLEAN HelpRequested()
{
  return(ExistParameter('h'));
}

/* =========================== DetermineOutputPrecision ============= */
int DetermineOutputPrecision()
{
  int i;
  int precis;
  /* BEGIN */
  for (i=1; i<GlobalArgc; i++) {
    if (GlobalArgv[i][0] == '-') {
      if (GlobalArgv[i][1] == 'p') {
        sscanf(&GlobalArgv[i][2],"%d",&precis);
        return(precis);
      }
    }
  }
  return(7);
}

/* =========================== NameParameter ========================= */
#ifdef ANSIPROTO
char *NameParameter(char ch)
#else
char *NameParameter(ch)
char ch;
#endif
{
  int I;
  /* BEGIN */
  for (I=1; I<GlobalArgc; I++) {
    if (GlobalArgv[I][0] == '-') {
      if (GlobalArgv[I][1] == ch) {
        return((char *) &GlobalArgv[I][2]);
      }
    }
  }
  return(NULL);
}

/* =========================== IntegerParameter ======================== */
#ifdef ANSIPROTO
int IntegerParameter(char ch,int Default,int MinVal,int MaxVal)
#else
int IntegerParameter(ch,Default,MinVal,MaxVal)
char ch;
int Default,MinVal,MaxVal;
#endif
{
  int I,TempVal,val;
  /* BEGIN */
  for (I=1; I<GlobalArgc; I++) {
    if (GlobalArgv[I][0] == '-') {
      if (GlobalArgv[I][1] == ch) {
        if(!(val=sscanf(&GlobalArgv[I][2],"%d",&TempVal)) || val == EOF) return(Default);
        if (TempVal < MinVal) return(MinVal);
        if (TempVal > MaxVal) return(MaxVal);
        return(TempVal);
      }
    }
  }
  return(Default);
}

/* ============================ IntegerPosParam ======================== */
#ifdef ANSIPROTO
int IntegerPosParam(int Item,int Default)
#else
int IntegerPosParam(Item,Default)
int Item,Default;
#endif
{
  int i,k;
  int TempVal;
  /* BEGIN */
  k = 0;
  for (i=1; (i<GlobalArgc) && (k!=Item); i++) {
    if (GlobalArgv[i][0] != '-') k++;
  }
  i--;
  if (k == Item) {
    if (sscanf(GlobalArgv[i],"%d",&TempVal)) {
      return(TempVal);
    }
  }
  return(Default);
}

/* ======================== RealParameter ============================== */
#ifdef ANSIPROTO
void RealParameter(char ch,double *deflt,double MinVal,double MaxVal)
#else
void RealParameter(ch, deflt, MinVal, MaxVal)
char ch;
double *deflt, MinVal, MaxVal;
#endif
{
  int I;
  /* BEGIN */
  for (I=1; I<GlobalArgc; I++) {
    if (GlobalArgv[I][0] == '-') {
      if (GlobalArgv[I][1] == ch) {
        sscanf(&GlobalArgv[I][2],"%lf",deflt);
        if (*deflt < MinVal) *deflt = MinVal;
        if (*deflt > MaxVal) *deflt = MaxVal;
        return;
      }
    }
  }
  return;
}

/* =========================== RealPosParam =========================== */
#ifdef ANSIPROTO
void RealPosParam(int Item,double *deflt)
#else
void RealPosParam(Item, deflt)
int Item;
double *deflt;
#endif
{
  int i,k;
  /* BEGIN */
  k = 0;
  for (i=1; (i<GlobalArgc) && (k!=Item); i++) {
    if ( (GlobalArgv[i][0] != '-') ||
        ((GlobalArgv[i][0] == '-') && (isdigit(GlobalArgv[i][1])!=0)) ) k++;
  }
  i--;
  if (k == Item) sscanf(GlobalArgv[i],"%lf",deflt);
}

/* ================================ DisplayReal ======================= */
#ifdef ANSIPROTO
char * DisplayReal(double value,int precision,char *buff)
#else
char * DisplayReal(value,precision,buff)
double value;
int precision;
char *buff;
#endif
{
  char format[80];
  int decimals;
  ELEMENTVALUETYPE power10;
  /* BEGIN */
  decimals = precision - 2;
  if (value < 0.0) decimals--;
  power10 = (ELEMENTVALUETYPE) 10.0;
  while ((decimals > 0) && (fabs(value) > power10)) {
    decimals--;
    power10 = power10 * 10;
  }
  sprintf_s(format,"%%-%d.%dg",precision,decimals);
  sprintf(buff,format,value);
  return buff;
}

/* ================================ NullName ========================== */
#ifdef ANSIPROTO
BOOLEAN NullName(const char *Name)
#else
BOOLEAN NullName(Name)
char *Name;
#endif
{
  if (Name == NULL) return(TRUE);
  if (strcmp(Name,"") == 0)
     return(TRUE);
  else
     return(FALSE);
}

/* =============================== OpenInput ========================== */
//Opens a file with file name 'Name' for reading. Returns the input stream pointer
#ifdef ANSIPROTO
FILE *OpenInput(const char *Name)
#else
FILE * OpenInput(Name)
char * Name;
#endif
{
  FILE * InputFile;
  char s[80];

#ifdef WINDOWS
  char charDir[300];
  strcpy_s(charDir, dir.GetBuffer(dir.GetLength()));
#endif

  if (!NullName(Name)) {
#ifdef WINDOWS
      if ((InputFile = fopen(strcat(charDir, Name),"rt")) == NULL) {
#else
      if ((InputFile = fopen(Name,"rt")) == NULL) {
#endif
         sprintf_s(s,"No input data file -> %s",Name);
         ErrorHalt(s);
         exit(ERROREXIT);
         /* fprintf(stderr,"Input from standard input\n");
         InputFile = stdin;*/
      } 
  } else {
    sprintf_s(s,"No input data file -> %s",Name);
    ErrorHalt(s);
    exit(ERROREXIT);
    /* InputFile = stdin;*/
  }
  return(InputFile);
}

/* ============================== OpenOutput ========================== */
#ifdef ANSIPROTO
FILE *OpenOutput(const char *Name)
#else
FILE * OpenOutput(Name)
char * Name;
#endif
{
  FILE * OutputFile;

#ifdef WINDOWS
  char charDir[300];
  strcpy_s(charDir, dir.GetBuffer(dir.GetLength()));
#endif

  if (!NullName(Name)) {
#ifdef WINDOWS
    if ((OutputFile = fopen(strcat(charDir, Name),"wt")) == NULL) {
#else
    if ((OutputFile = fopen(Name,"wt")) == NULL) {
#endif
       fprintf(stderr,"Output to standard output\n");
#ifdef WINDOWS
	   OutputFile = NULL;
#else
	   OutputFile = stdout;
#endif
	   //explanation for above: a really strange bug where stdout gets initialized to something 
	   //by itself. I(that is, Shu) really couldn't figure it out so I set OutputFile to NULL instead which will
	   //print to screen anyways. 
    }
  }
  else
#ifdef WINDOWS
	   OutputFile = NULL;
#else
	   OutputFile = stdout;
#endif
  return(OutputFile);
}

#ifdef ANSIPROTO
int SizeParameter(char ch,int i,int imax)
#else
int SizeParameter(ch,i,imax)
char ch;
int i,imax;
#endif
{
  double TempSize;
  TempSize = 1.0;
  RealParameter(ch,&TempSize,(double) 0.02,(double) imax);
  if (TempSize > (double) 1.0) return(IntegerParameter(ch,i,50,imax));
  else                         return((int) (TempSize * (double) imax));
}
