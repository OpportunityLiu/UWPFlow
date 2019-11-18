#pragma once

/* param.h 030390 */

#include "constant.h"

void SetArguments(int argc, const char **argv);
bool ExistParameter(char ch);
char *TrueParamStr(int Item);
bool HelpRequested(void);
int DetermineOutputPrecision(void);
char *NameParameter(char ch);
int IntegerParameter(char ch, int Default, int MinVal, int MaxVal);
int IntegerPosParam(int Item, int Default);
void RealParameter(char ch, double *Default, double MinVal, double MaxVal);
void RealPosParam(int Item, double *Default);
char *DisplayReal(double value, int precision, char *buffer);
bool NullName(const char *Name);
FILE *OpenInput(const char *Name);
FILE *OpenOutput(const char *Name);
int SizeParameter(char ch, int i, int imax);
