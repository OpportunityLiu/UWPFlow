#pragma once

#include "constant.h"

struct AreaData {
  INDEX N;
  INDEX i;
  char Name[31];
  VALUETYPE P;
  VALUETYPE SPg;
  char Zone[11][4];
  struct ACbusData *Slack;
  struct ACbusData *BSptr;
  struct AClist *AC;
  struct DClist *DC;
  struct ElementList *Elem;
  struct AreaData *Next;
};

struct ACbusData {
  char Name[13];
  INDEX N;
  INDEX Num;
  VALUETYPE KV;
  char Type[4];
  char Zone[4];
  char Owner[4];
  char cont[2];
  struct AreaData *Area;
  int Ncont;
  struct ElementList *Reg;
  struct ElementList *Elem;
  struct GenModel *Gen; /* Generator Steady State Model */
  VALUETYPE V;
  VALUETYPE VCont;
  VALUETYPE Ang;
  VALUETYPE Pg;
  VALUETYPE PG;
  VALUETYPE Qg;
  VALUETYPE Pl;
  VALUETYPE PL;
  VALUETYPE Ql;
  VALUETYPE QL;
  VALUETYPE G;
  VALUETYPE B;
  VALUETYPE Bz;
  INDEX step;
  INDEX steps;
  VALUETYPE Bx[73];
  VALUETYPE PgMax;
  VALUETYPE Pmax;
  INDEX flagPgMax;
  VALUETYPE a;
  VALUETYPE b;
  VALUETYPE Pz;
  VALUETYPE Pn;
  VALUETYPE Qz;
  VALUETYPE Qn;
  VALUETYPE Pzl;
  VALUETYPE Pnl;
  VALUETYPE Qzl;
  VALUETYPE Qnl;
  VALUETYPE Kg;
  VALUETYPE Qmax;
  VALUETYPE Qmin;
  VALUETYPE Max;
  VALUETYPE Min;
  VALUETYPE Smax;
  VALUETYPE Vmax;
  VALUETYPE Vmin;
  VALUETYPE Vlmax;
  VALUETYPE Vlmin;
  bool CheckVlimits;
  VALUETYPE DPg;
  VALUETYPE DPG;
  VALUETYPE Qr;
  VALUETYPE Kbg;
  VALUETYPE Kbg1;
  VALUETYPE Kbg2;
  VALUETYPE val;
  VALUETYPE valp;
  VALUETYPE vals;
  VALUETYPE valt;
  INDEX Nc;
  struct DClist *DC;
  struct SVClist *SVC;         /* FACTS */
  struct TCSClist *TCSC;       /* FACTS */
  struct STATCOMlist *STATCOM; /* FACTS */
  struct ACbusData *Cont;
  struct AClist *ContBus;
  struct ACbusData *Next;
  struct ACbusData *Prev;
};

/* Generator Steady State Model */
struct GenModel {
  INDEX Nvar;
  VALUETYPE Ra;
  VALUETYPE Xd;
  VALUETYPE Xq;
  VALUETYPE IaMax;
  VALUETYPE EqMax;
  VALUETYPE EqMin;
  VALUETYPE Eq;
  VALUETYPE dg;
  VALUETYPE Vr;
  VALUETYPE Vi;
  VALUETYPE Ir;
  VALUETYPE Ii;
  VALUETYPE Vq;
  VALUETYPE Vd;
  VALUETYPE Iq;
  VALUETYPE Id;
  VALUETYPE Ia;
};

struct DCbusData {
  char Name[13];
  INDEX N;
  char Type[2];
  char Cont1[3];
  char Cont2[3];
  char Zone[4];
  char Lzone[4];
  struct DCbusData *Meter;
  struct AreaData *Area;
  VALUETYPE Xc;
  VALUETYPE Nbr;
  VALUETYPE Ntrf;
  VALUETYPE MVA;
  VALUETYPE Vd;
  VALUETYPE VdN;
  VALUETYPE Id;
  VALUETYPE P;
  VALUETYPE Q;
  VALUETYPE Alfa;
  VALUETYPE AlfaN;
  VALUETYPE Gamma;
  VALUETYPE AlfaMin;
  VALUETYPE AlfaMax;
  VALUETYPE GammaMin;
  VALUETYPE Tap;
  VALUETYPE TapMin;
  VALUETYPE TapMax;
  VALUETYPE Vn;
  VALUETYPE Rd;
  VALUETYPE Ld;
  VALUETYPE val[4];
  ACbusData *AC;
  DCbusData *To;
  DCbusData *Next;
};

struct ElementData {
  struct ACbusData *From;
  struct ACbusData *To;
  char Ckt[2];
  char Zone[4];
  char Owner[4];
  char Type[3];
  struct AreaData *Area;
  struct ACbusData *Meter;
  struct ACbusData *Cont;
  INDEX Sec;
  INDEX Ncont;
  VALUETYPE G;
  VALUETYPE B;
  VALUETYPE G1;
  VALUETYPE B1;
  VALUETYPE G2;
  VALUETYPE B2;
  VALUETYPE Tap;
  VALUETYPE Taps;
  VALUETYPE Ang;
  char Ctype[2];
  VALUETYPE Cvar;
  VALUETYPE Tmin;
  VALUETYPE Tmax;
  VALUETYPE Min;
  VALUETYPE Max;
  VALUETYPE Imax;
  bool CheckIlimits;
  VALUETYPE val;
  struct ElementData *Next;
  struct ElementData *Prev;
};

struct AClist {
  struct ACbusData *AC;
  struct AClist *Next;
  struct AClist *Prev;
  struct AreaData *Area;
  INDEX N;
  char Type[3];
};

struct DClist {
  struct DCbusData *DC;
  struct DClist *Next;
};

typedef struct ElementList {
  struct ElementData *Eptr;
  struct ElementList *Next;
} ElementList;

struct Data {
  char Title[3][BUFLEN + 1];
  struct ACbusData *ACbus;
  struct DCbusData *DCbus;
  struct ElementData *Element;
  struct AreaData *Area;
  struct SVCbusData *SVCbus;         /* FACTS */
  struct TCSCbusData *TCSCbus;       /* FACTS */
  struct STATCOMbusData *STATCOMbus; /* FACTS */
  struct AClist *KGbus;
};

/* FACTS */

struct SVCbusData {
  char Name[13];
  INDEX N;
  char Type[3];
  char Cont[3];
  VALUETYPE Xth_l;
  VALUETYPE Vsvc;
  VALUETYPE Xc;
  VALUETYPE Xl;
  VALUETYPE AlphaMin;
  VALUETYPE AlphaMax;
  VALUETYPE slope;
  VALUETYPE SVC_base;
  VALUETYPE Qsvc;
  VALUETYPE Bv;
  VALUETYPE alpha_svc;
  VALUETYPE Vref;
  VALUETYPE val;
  VALUETYPE Vvar;
  struct ACbusData *From;
  struct ACbusData *Ctrl;
  struct SVCbusData *Prev;
  struct SVCbusData *Next;
};

struct TCSCbusData {
  char Name[13];
  INDEX N;
  char Type[3];
  char Cont[2];
  VALUETYPE Xc;
  VALUETYPE Xl;
  VALUETYPE AlphaMin;
  VALUETYPE AlphaMax;
  VALUETYPE Control;
  VALUETYPE Bset;
  VALUETYPE TCSC_base;
  VALUETYPE Ptcsc;
  VALUETYPE Qtcsck;
  VALUETYPE Qtcscm;
  VALUETYPE Be;
  VALUETYPE alpha_tcsc;
  VALUETYPE Itcsc;
  VALUETYPE delta_t;
  VALUETYPE val;
  VALUETYPE Max;
  struct ACbusData *From;
  struct ACbusData *To;
  struct TCSCbusData *Prev;
  struct TCSCbusData *Next;
};

struct STATCOMbusData {
  char Name[13];
  INDEX N;
  char Type[3];
  char Cont[3];
  char Cont1[3];
  VALUETYPE I;
  VALUETYPE theta;
  VALUETYPE k;
  VALUETYPE Vdc;
  VALUETYPE alpha;
  VALUETYPE R;
  VALUETYPE G;
  VALUETYPE B;
  VALUETYPE Gc;
  VALUETYPE Imin;
  VALUETYPE Imax;
  VALUETYPE slope;
  VALUETYPE P;
  VALUETYPE Q;
  VALUETYPE MVA;
  VALUETYPE Vref;
  VALUETYPE Contref;
  VALUETYPE val;
  VALUETYPE Vvar;
  struct ACbusData *From;
  struct ACbusData *Ctrl;
  struct STATCOMbusData *Prev;
  struct STATCOMbusData *Next;
};

struct SVClist {
  struct SVCbusData *SVC;
  struct SVClist *Next;
};

struct TCSClist {
  struct TCSCbusData *TCSC;
  struct TCSClist *Next;
};

struct STATCOMlist {
  struct STATCOMbusData *STATCOM;
  struct STATCOMlist *Next;
};
