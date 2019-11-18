#pragma once

#include "constant.h"
#include "param.h"
#include "pflow.h"
#include "sparse.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void JacElement(SparseMatrix *Mptr, INDEX I, INDEX J, VALUETYPE val);
void DeleteJac(SparseMatrix *Mptr, IntegerVector *P1Row, IntegerVector *P1Col,
               IntegerVector *P2Row, IntegerVector *P2Col);
AreaData *ACFunJac(SparseMatrix *Mptr, int *val, bool flagF, bool flagJ,
                   bool flagP);
bool DCFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ);
void SVCFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ);     /* FACTS */
void TCSCFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ);    /* FACTS */
void STATCOMFunJac(SparseMatrix *Mptr, bool flagF, bool flagJ); /* FACTS */
int HFunJac(bool FlagFunction, bool FlagJacobian, AreaData *Aptr,
            VALUETYPE *vec);
VALUETYPE LoadX0(bool FlagLoadX0, bool FlagUpdateVar, bool FlagMakeDxZero);
int factorns(SparseMatrix *Mptr, double Param, IntegerVector *PartRow,
             IntegerVector *PartCol, IntegerVector *P1Row, IntegerVector *P1Col,
             IntegerVector *P2Row, IntegerVector *P2Col);
int factor(SparseMatrix *Mptr);
void repsolp(SparseMatrix *Mptr, VALUETYPE *Vptr, IntegerVector *PermR,
             IntegerVector *PermC);
void WriteSolution(INDEX Iter, char *File1, char *str);
int Pflow(int iter, bool flagF, bool flagD, bool flagP);
void MakeVlist(FILE *Out);
void VoltProf(bool flag, FILE *Out);
int Direction(SparseMatrix *Mptr, VALUETYPE *vec, bool flag);
bool ChangeParam(void);
int Homot(void);
bool CheckRlimits(void);
bool CheckVlimits(void);
bool CheckQlimits(void);
bool ChangeDCmode(void);
bool ChangeSVCmode(void);     /* FACTS */
bool ChangeTCSCmode(void);    /* FACTS */
bool ChangeSTATCOMmode(void); /* FACTS */
bool DCsetup(void);
void Print(FILE *File, int spaces, int width, int decimals, VALUETYPE val);
void VoltageSensFactor(VALUETYPE *vec, bool first);
void TEFac(bool flag);
void TEFdc(FILE *Out);
void MatlabV(FILE *Out);
void TEFMatlabFiles(void);
bool InList(ACbusData *ACptr, AClist *Vptr);
void PrintDirection(char Option, VALUETYPE *vector, VALUETYPE Max);
void WriteQmatrix(INDEX count, VALUETYPE *vec);
void IndicesMatlab(INDEX count);

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX Nac, NacEl, Ndc, Narea, NacVar, Nvolt, NregV, Nslack, Nsvc, Ntcsc,
    NtcscVar, Nstatcom; /* FACTS */
extern INDEX *ACvar, Bl;
extern VALUETYPE *dF, lambda, alpha, Vac, Sn, *dx, Tol;
extern IntegerVector *NewRow, *OldRow, *NewCol, *OldCol, *RowPartition,
    *ColPartition;
extern bool Acont, PQcont, QRcont, Rcont, PQlim, Tlim, Qlim, Vlim, Elim, Ilim,
    Xlim, Zlim, flagH, flagL, flagR, flagReducedContinuation, flagPgMax,
    flagSmax;
extern ACbusData *BlPtr;
extern INDEX TFnum, TFbus;
extern int DetSign;
extern char TFname[13];

typedef struct ACranklist {
  struct ACbusData *AC;
  VALUETYPE val;
  struct ACranklist *Next;
} ACranklist;
