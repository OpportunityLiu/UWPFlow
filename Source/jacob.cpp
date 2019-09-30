/* Initialize variables. */

#include "jacob.h"

/* ------- Global Variables ------ */
IntegerVector *RowPer,*ColPer;

/* -------------------- JacElement ---------------------- */
#ifdef ANSIPROTO
void JacElement(SparseMatrix *Mptr,INDEX I,INDEX J,VALUETYPE val)
#else
void JacElement(Mptr,I,J,val)
SparseMatrix *Mptr; 
INDEX I,J;
VALUETYPE val;
#endif
/* Jacobian Element */
{
  SparseMatrixElement *Ptr; 
  INDEX i,j;

/*  if (val==0) return;*/
  i=I;  j=J;
  if (RowPer->p[I]!=0) I=RowPer->p[I];
  if (ColPer->p[J]!=0) J=ColPer->p[J];
  for(Ptr=Mptr->RowHead[I];Ptr!=NULL;Ptr=Ptr->RowNext)
    if(Ptr->Row==I && Ptr->Col==J) {
      if (flagReducedContinuation && (DxZero[i] || DxZero[j])) {
        if (i==j) Ptr->Value=1.;
        else Ptr->Value=0;
      } else Ptr->Value=Ptr->Value+val;
      break;
    }
  if (Ptr==NULL && (i==j || !(flagReducedContinuation && (DxZero[i] || DxZero[j])))) {
#ifdef WINDOWS
    Ptr=new SparseMatrixElement;
#else
    Ptr=(SparseMatrixElement *) malloc(sizeof(SparseMatrixElement));
    if(Ptr==NULL) {ErrorHalt("Insufficient memory to allocate Jacobian."); exit(ERROREXIT);}
#endif
    Ptr->Row=I;
    Ptr->Col=J;
    Ptr->RowNext=Mptr->RowHead[I];
    Mptr->RowHead[I]=Ptr;
    Ptr->ColNext=Mptr->ColHead[J];
    Mptr->ColHead[J]=Ptr;
    if (flagReducedContinuation && (DxZero[i] || DxZero[j])) Ptr->Value=1.;
    else Ptr->Value=val;
  }
}


/* ------------------- Jacobian --------------------- */
#ifdef ANSIPROTO
void Jacobian(void)
#else
void Jacobian()
#endif
/* Allocate memory for the Jacobian. */
{
  ACbusData *ACptr;
  INDEX i,N,N1,N2;

  Dx=x0=x0p=dx=dF=NULL;
  DxZero=NULL;
  Vlist=Vlistp=NULL;
#ifdef WINDOWS
  ACvar= new INDEX[Nac+1];
#else
  ACvar=(INDEX *) calloc(Nac+1,sizeof(INDEX));
  if (ACvar==NULL) {ErrorHalt("Insufficient memory to allocate AC variable vector."); exit(ERROREXIT);}
#endif
  N=1;
  for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next){
    ACvar[ACptr->N]=N;
    if (PQcont) N1=ACptr->Ncont;
    else N1=0;
    if (Acont && strpbrk(ACptr->Type,"A")) N=N+3+N1;
    else N=N+2+N1;
    if (ACptr->Gen!=NULL) {
      ACptr->Gen->Nvar=N-1;
      N=N+11;
    }
  }
  NacVar=N-1;
  N1=NacVar+11*Ndc/2 +3*Nsvc+NtcscVar+7*Nstatcom;    /* FACTS  */
#ifdef WINDOWS
  RowPartition= new IntegerVector;
  ColPartition= new IntegerVector;
#else
  RowPartition=(IntegerVector *) malloc(sizeof(IntegerVector));
  ColPartition=(IntegerVector *) malloc(sizeof(IntegerVector));
  if (ColPartition==NULL) {ErrorHalt("Insufficient memory to allocate partition vectors."); exit(ERROREXIT);}
#endif
  if (flagH) {
    N1++; N2=N1;
    RowPartition->N=2;
    ColPartition->N=1;
#ifdef WINDOWS
    RowPartition->p= new INDEX[3];
    ColPartition->p= new INDEX[3];
#else
    RowPartition->p=(INDEX *) calloc(3,sizeof(INDEX));
    ColPartition->p=(INDEX *) calloc(3,sizeof(INDEX));
    if (ColPartition->p==NULL) {ErrorHalt("Insufficient memory to allocate partition vectors."); exit(ERROREXIT);}
#endif
    RowPartition->p[1]=N1-1;
    RowPartition->p[2]=N1;
    ColPartition->p[1]=N1;
    RowPartition->p[0]=ColPartition->p[0]=ColPartition->p[2]=0;
#ifdef WINDOWS
    Dx= new VALUETYPE[N1+1];
    x0= new VALUETYPE[N1+1];
    x0p= new VALUETYPE[N1+1];
#else
    Dx=(VALUETYPE *) calloc(N1+1,sizeof(VALUETYPE));
    x0=(VALUETYPE *) calloc(N1+1,sizeof(VALUETYPE));
    x0p=(VALUETYPE *) calloc(N1+1,sizeof(VALUETYPE));
    if (x0p==NULL) {ErrorHalt("Insufficient memory to allocate Homotopy vectors."); exit(ERROREXIT);}
#endif
    for(i=0;i<N1+1;i++) Dx[i]=x0[i]=x0p[i]=0.;
  } else if (flagPoC) {
    N2=N1;
    N1=2*N1+1;
    RowPartition->N=ColPartition->N=1;
#ifdef WINDOWS
    RowPartition->p= new INDEX[2];
    ColPartition->p= new INDEX[2];
#else
    RowPartition->p=(INDEX *) calloc(2,sizeof(INDEX));
    ColPartition->p=(INDEX *) calloc(2,sizeof(INDEX));
    if (ColPartition->p==NULL) {ErrorHalt("Insufficient memory to allocate partition vectors."); exit(ERROREXIT);}
#endif
    RowPartition->p[1]=ColPartition->p[1]=N2;
    RowPartition->p[0]=ColPartition->p[0]=0;
#ifdef WINDOWS
    Dx= new VALUETYPE[N2+1];
    x0= new VALUETYPE[N2+1];
#else
    Dx=(VALUETYPE *) calloc(N2+1,sizeof(VALUETYPE));
    x0=(VALUETYPE *) calloc(N2+1,sizeof(VALUETYPE));
    if (x0==NULL) {ErrorHalt("Insufficient memory to allocate left Eigenvector."); exit(ERROREXIT);}
#endif
    for(i=0;i<N2+1;i++) Dx[i]=x0[i]=0.;
  } else {
    N2=N1;
    RowPartition->N=ColPartition->N=1;
#ifdef WINDOWS
    RowPartition->p= new INDEX[2];
    ColPartition->p= new INDEX[2];
#else
    RowPartition->p=(INDEX *) calloc(2,sizeof(INDEX));
    ColPartition->p=(INDEX *) calloc(2,sizeof(INDEX));
    if (ColPartition->p==NULL) {ErrorHalt("Insufficient memory to allocate partition vectors."); exit(ERROREXIT);}
#endif
    RowPartition->p[1]=ColPartition->p[1]=N2;
    RowPartition->p[0]=ColPartition->p[0]=0;
  }
#ifdef WINDOWS
  dx= new VALUETYPE[N1+2];
  dF= new VALUETYPE[N1+2];
  Jac= new SparseMatrix;
#else
  dx=(VALUETYPE *) calloc(N1+2,sizeof(VALUETYPE));
  dF=(VALUETYPE *) calloc(N1+2,sizeof(VALUETYPE));
  if (dF==NULL) {ErrorHalt("Insufficient memory to allocate AC/DC/FACTS mismatches."); exit(ERROREXIT);}
  Jac=(SparseMatrix *) malloc(sizeof(SparseMatrix));
  if (Jac==NULL) {ErrorHalt("Insufficient memory to allocate P.F. Jacobian."); exit(ERROREXIT);}
#endif
  for(i=0;i<N1+2;i++) dx[i]=dF[i]=0.;
  Jac->n1=N2;
  Jac->n2=N2;
#ifdef WINDOWS
  Jac->RowHead= new SparseMatrixElement*[N1+2];
  Jac->ColHead= new SparseMatrixElement*[N1+2];
#else
  Jac->RowHead=(SparseMatrixElement **) calloc(N1+2,sizeof(SparseMatrixElement *));
  Jac->ColHead=(SparseMatrixElement **) calloc(N1+2,sizeof(SparseMatrixElement *));
  if(Jac->ColHead==NULL) {ErrorHalt("Insufficient memory to allocate Jacobian."); exit(ERROREXIT);}
#endif
  for(i=0;i<N1+2;i++) Jac->RowHead[i]=Jac->ColHead[i]=NULL;
  NewRow=AllocatePermutation(N1);
  OldRow=AllocatePermutation(N1);
  NewCol=AllocatePermutation(N1);
  OldCol=AllocatePermutation(N1);
  NewRow->N=NewCol->N=N2;
  OldRow->N=OldCol->N=N2;
}

