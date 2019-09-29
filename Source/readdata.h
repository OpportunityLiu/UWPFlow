#include <stdlib.h>
//#ifndef WINDOWS
//#include <stdio.h>
//#else
#include "pfwstdio.h"
//#endif
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"

#ifdef ANSIPROTO
char *GetStr(const char *ptr,int Pos,int Leng,int Tot,char *str);
VALUETYPE GetValue(const char *ptr,int Pos,int Leng,int Dec);
INDEX GetInt(const char *ptr,int Pos,int Leng);
ACbusData *ACbusInList(INDEX BusN,char *BusName,VALUETYPE V,INDEX N1,INDEX N2);
SVCbusData *SVCbusInList(char *BusName,INDEX N,ACbusData *ptrac,ACbusData *ptrac1);         /* FACTS */
TCSCbusData *TCSCbusInList(char *BusName,INDEX N,ACbusData *ptrac,ACbusData *ptrac1);       /* FACTS */
STATCOMbusData *STATCOMbusInList(char *BusName,INDEX N,ACbusData *ptrac,ACbusData *ptrac1); /* FACTS */
AreaData *AreaInList(INDEX i,char *Name,INDEX N);
ElementList *AddElemToList(ElementList *ELptr,ElementData *Eptr);
ElementData *ElemInList(ACbusData *From,ACbusData *To,INDEX N1,INDEX N2,char *Type,char *Ckt);
DCbusData *DCbusInList(char *BusName,INDEX N);
void ReadWSCC(void);
void ReadIEEE(void);
void ReadEPRIdc(const char *Line);
void ReadSVC(const char *Line);               /* FACTS */
void ReadTCSC(const char *Line);              /* FACTS */
void ReadSTATCOM(const char *Line);           /* FACTS */
void Multiply(VALUETYPE *a,VALUETYPE *b,VALUETYPE c,VALUETYPE d);
void Divide(VALUETYPE *a,VALUETYPE *b,VALUETYPE c,VALUETYPE d);
BOOLEAN AddSection(ACbusData *From,ACbusData *To,char *Line,char *Ckt,INDEX Sec);
void ErrorDetect(void);
void WriteSummary(void);
void ExpandSlack(ACbusData *BSptr,AreaData *Aptr);
void ReadData(char *Name);
BOOLEAN ReadInit(void);
BOOLEAN ReadOHload(char *File);
void ReadITALY(void);
#else
char *GetStr();
VALUETYPE GetValue();
INDEX GetInt();
ACbusData *ACbusInList();
SVCbusData *SVCbusInList();          // FACTS
TCSCbusData *TCSCbusInList();        // FACTS
STATCOMbusData *STATCOMbusInList();  // FACTS
AreaData *AreaInList();
ElementList *AddElemToList();
ElementData *ElemInList();
DCbusData *DCbusInList();
void ReadWSCC();
void ReadIEEE();
void ReadEPRIdc();
void ReadSVC();       // FACTS
void ReadTCSC();      // FACTS
void ReadSTATCOM();   // FACTS
void Multiply();
void Divide();
BOOLEAN AddSection();
void ErrorDetect();
void WriteSummary();
void ExpandSlack();
void ReadData();
BOOLEAN ReadInit();
BOOLEAN ReadOHload();
void ReadITALY();
#endif   


/* ------- Global Variables (some defined in pflow.c) ------ */
extern Data *dataPtr;
extern INDEX MaxIter,Nac,NacEl,NregV,NregPQ,Ndc,Nslack,Nvolt,Narea,Ngen,NZvolt,NXvolt,
	            Nsvc,Ntcsc,NtcscVar,Nstatcom; /* FACTS */
extern VALUETYPE lambda,Sn,K3;
extern VALUETYPE AngSlack;
extern INDEX NdcEl,LineNum,Bl;
extern ACbusData *BlPtr;
extern int GlobalArgc;
extern char **GlobalArgv;
extern BOOLEAN InputError,Acont,PQcont,QRcont,Rcont,Xcont,flagVloads,flagKdirection,flag2Vcontrol;
