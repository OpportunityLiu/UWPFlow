/* AC/DC Transient Energy Functions. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constant.h"
#include "param.h"
#include "sparse.h"
#include "pflow.h"

#ifdef ANSIPROTO
void TEFac(BOOLEAN flag);
void TEFdc(FILE *Out);
void MatlabV(FILE *Out);
void TEFMatlabFiles(void);
void Print(FILE *File,int spaces,int width,int decimals,VALUETYPE val);
#else
void TEFac();
void TEFdc();
void MatlabV();
void TEFMatlabFiles();
void Print();
#endif

/* ------- Global Variables ------ */
extern Data *dataPtr;
extern SparseMatrix *Jac;
extern INDEX Nac,Ngen,Ndc,NacVar,Narea;
extern INDEX *ACvar;
extern VALUETYPE *x0p,lambda,Sn;
extern AClist *Vlist;
extern BOOLEAN flagBS,flagPgMax;
VALUETYPE Vac;
extern int field;
extern BOOLEAN flagPrintTotalPl,flagPrintTotalQl,flagPrintTotalPg,flagPrintTotalQg;


/* ------------------ TEFac ----------------------------- */
#ifdef ANSIPROTO
void TEFac(BOOLEAN flag)
#else
void TEFac(flag)
BOOLEAN flag;
#endif
/* Caculate AC TEF. */
{
  ACbusData *ACptr,*BEptr,*To;
  ElementList *ELptr;
  ElementData *Eptr;
  VALUETYPE Vi,Vj,di,dj,Pi,Qi,Pl,Ql,val;
  VALUETYPE gii,bii,gij,bij,gsij,bsij,DPg,Pg;
  INDEX i,j;

  if (Narea<2){
    for(BEptr=dataPtr->ACbus;BEptr!=NULL;BEptr=BEptr->Next)
      if(strpbrk(BEptr->Type,"S")) break;
  }
  for(Vac=0,ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next) {
    if(ACptr->Area!=NULL) BEptr=ACptr->Area->Slack;
    if (flagBS) {if (ACptr!=BEptr) DPg=0; else DPg=1;}
    else DPg=ACptr->DPg;
    i=ACvar[ACptr->N];
    Vi=ACptr->V;
    di=ACptr->Ang;
    if (flag) {
      x0p[i]=di;
      x0p[i+1]=Vi;
    }
    Pl = (ACptr->Pn+lambda*ACptr->Pnl)*pow(Vi,ACptr->a);
    Ql = (ACptr->Qn-lambda*ACptr->Qnl)*pow(Vi,ACptr->b);
    Pg=ACptr->Pg+DPg*BEptr->Kg;
    if (!flagPgMax && Pg>ACptr->PgMax) Pg=ACptr->PgMax;
    Pi = Pg-Pl;
    Qi = -Ql;
    Qi=Qi+ACptr->Qg;
    gii=ACptr->G+ACptr->Pz+lambda*ACptr->Pzl;
    bii=ACptr->B+ACptr->Qz+lambda*ACptr->Qzl;
    for(val=0,ELptr=ACptr->Elem; ELptr!=NULL; ELptr=ELptr->Next) {
      Eptr=ELptr->Eptr;
      if (Eptr->From==ACptr) To=Eptr->To;
      else To=Eptr->From;
      j=ACvar[To->N];
      Vj=To->V;
      dj=To->Ang;
      if (flag) {
        x0p[j]=dj;
        x0p[j+1]=Vj;
      }
      if(Eptr->From==ACptr) {
        gij=(Eptr->G*cos(Eptr->Ang)-Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
        bij=(Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
        gsij=(Eptr->G1+Eptr->G)*pow(Eptr->Tap,2.0)-gij;
        bsij=(Eptr->B1+Eptr->B)*pow(Eptr->Tap,2.0)-bij;
      } else {
        gij=(Eptr->G*cos(Eptr->Ang)+Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
        bij=(-Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
        gsij=Eptr->G+Eptr->G2-gij;
        bsij=Eptr->B+Eptr->B2-bij;
      }
      val=val+0.5*bij*Vi*Vj*cos(di-dj);
      val=val-gij*x0p[i+1]*x0p[j+1]*cos(x0p[i]-x0p[j])*di-gij*x0p[j+1]*sin(x0p[i]-x0p[j])*Vi;
      gii=gii+gij+gsij;
      bii=bii+bij+bsij;
    }
    val=val-0.5*bii*Vi*Vi-Pi*di;
    if (Vi>0) val=val-Qi*log(Vi);
    val=val+gii*x0p[i+1]*x0p[i+1]*di;
    Vac=Vac+val;
  }
}

/* ------------------ TEFdc ----------------------------- */
#ifdef ANSIPROTO
void TEFdc(FILE *Out)
#else
void TEFdc(Out)
FILE *Out;
#endif
/* Caculate DC TEF and coupling terms. */
{
  DCbusData *DCptrR,*DCptrI,*DCptr;
  VALUETYPE a,b,d,t,cosag,Id,Xc;
  INDEX j;

  for(DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if(!strcmp(DCptrR->Type,"R")){
      for (j=1;j<=2;j++) {
        if (j==1) {
          DCptr=DCptrR;
          cosag=cos(DCptr->Alfa);
        } else {
          DCptr=DCptrI;
          cosag=cos(DCptr->Gamma);
        }
        Id=DCptr->Id;
        t=DCptr->Tap*DCptr->Ntrf;
        Xc=DCptr->Xc;
        a= -Xc*Xc*Id*Id;
        b=2*sqrt(2.0)*t*Xc*Id*cosag;
        /* c=2*t*t*(1-cosag*cosag); */
        d=8*t*t*Xc*Xc*Id*Id;
        Print(Out,0,field,4,a);  fprintf(Out,"    ");
        Print(Out,0,field,4,b);  fprintf(Out,"    ");
        Print(Out,0,field,4,d);  fprintf(Out,"    ");
      }
    }
  }
}


/* --------------------------- MatlabV --------------------------------- */
#ifdef ANSIPROTO
void MatlabV(FILE *Out)
#else
void MatlabV(Out)
FILE *Out;
#endif
/* Write Matlab commands in output file for ploting profiles. */
{
  AClist *Lptr;
  char LineType[4][5];
  DCbusData *DCptr,*DCptrp;
  INDEX count,countp,i,j,k,l,js,ks;
  VALUETYPE k1,k2,ki,kp,P11,P12,P22,P23,P33,Beta,KV;

  strcpy(LineType[0],"'-'");
  strcpy(LineType[1],"'-.'");
  strcpy(LineType[2],"':'");
  strcpy(LineType[3],"'--'");
  fprintf(Out,"];\n");
  fprintf(Out,"%s Plot profiles:\n","%%");
  fprintf(Out,"%s Change the value of K to scale L.F.\nK=1;\n","%%");
  fprintf(Out,"figure;\nplot(");
  for(count=countp=0,Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) {
    count++;
    if (countp<4) {
      if  (!strcmp(Lptr->Type,"V")) {
        KV=Lptr->AC->KV;
        if (KV<=0) KV=1;
        if (countp==0) fprintf(Out,"K*x(:,1),%5.1f*x(:,%1d)",KV,count+1);
        else           fprintf(Out,",K*x(:,1),%5.1f*x(:,%1d)",KV,count+1);
        fprintf(Out,",%s", LineType[countp]);
      } else {
        if (countp==0) fprintf(Out,"K*x(:,1),x(:,%1d)",count+1);
        else           fprintf(Out,",K*x(:,1),x(:,%1d)",count+1);
        fprintf(Out,",%s", LineType[countp]);
      }
    }
    countp++;
  }
  if (countp<4 && flagPrintTotalPl) {
    if (countp==0) fprintf(Out,"K*x(:,1),x(:,%1d)",count+1);
    else           fprintf(Out,",K*x(:,1),x(:,%1d)",count+1);
    fprintf(Out,",%s", LineType[countp]);
    countp++;
  }
  if (countp<4 && flagPrintTotalQl) {
    if (countp==0) fprintf(Out,"K*x(:,1),x(:,%1d)",count+1);
    else           fprintf(Out,",K*x(:,1),x(:,%1d)",count+1);
    fprintf(Out,",%s", LineType[countp]);
    countp++;
  }
  if (countp<4 && flagPrintTotalPg) {
    if (countp==0) fprintf(Out,"K*x(:,1),x(:,%1d)",count+1);
    else           fprintf(Out,",K*x(:,1),x(:,%1d)",count+1);
    fprintf(Out,",%s", LineType[countp]);
    countp++;
  }
  if (countp<4 && flagPrintTotalQg) {
    if (countp==0) fprintf(Out,"K*x(:,1),x(:,%1d)",count+1);
    else           fprintf(Out,",K*x(:,1),x(:,%1d)",count+1);
    fprintf(Out,",%s", LineType[countp]);
    countp++;
  }
  fprintf(Out,");\n");
  fprintf(Out,"legend(");
  for(countp=0,Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) {
    if (countp<4) {
      if (Lptr->AC!=NULL && !strcmp(Lptr->Type,"V")) {
        if (countp==0) fprintf(Out,"'kV_{%s}'",Lptr->AC->Name);
        else           fprintf(Out,",'kV_{%s}'",Lptr->AC->Name);
      }
      else if (Lptr->AC!=NULL && !strcmp(Lptr->Type,"D")) {
        if (countp==0) fprintf(Out,"'deg_{%s}'",Lptr->AC->Name);
        else           fprintf(Out,",'deg_{%s}'",Lptr->AC->Name);
      }
      else if (Lptr->AC!=NULL && !strcmp(Lptr->Type,"PL")) {
        if (countp==0) fprintf(Out,"'L MW_{%s}'",Lptr->AC->Name);
        else           fprintf(Out,",'L MW_{%s}'",Lptr->AC->Name);
      }
      else if (Lptr->AC!=NULL && !strcmp(Lptr->Type,"QL")) {
        if (countp==0) fprintf(Out,"'L MVar_{%s}'",Lptr->AC->Name);
        else           fprintf(Out,",'L MVar_{%s}'",Lptr->AC->Name);
      }
      else if (Lptr->AC!=NULL && !strcmp(Lptr->Type,"PG")) {
        if (countp==0) fprintf(Out,"'G MW_{%s}'",Lptr->AC->Name);
        else           fprintf(Out,",'G MW_{%s}'",Lptr->AC->Name);
      }
      else if (Lptr->AC!=NULL && !strcmp(Lptr->Type,"QG")) {
        if (countp==0) fprintf(Out,"'G MVar_{%s}'",Lptr->AC->Name);
        else           fprintf(Out,",'G MVar_{%s}'",Lptr->AC->Name);
      }
      else if(Lptr->Area!=NULL) {
        if (countp==0) fprintf(Out,"'A MW_{%s}'",Lptr->Type,Lptr->Area->Name);
        else           fprintf(Out,",'A MW_{%s}'",Lptr->Type,Lptr->Area->Name);
      }
    }
    countp++;
  }
  if (countp<4 && flagPrintTotalPl) {
    if (countp==0) fprintf(Out,"'L MW_{TOTAL}'");
    else           fprintf(Out,",'L MW_{TOTAL}'");
    countp++;
  }
  if (countp<4 && flagPrintTotalQl) {
    if (countp==0) fprintf(Out,"'L MVar_{TOTAL}'");
    else           fprintf(Out,",'L MVar_{TOTAL}'");
    countp++;
  }
  if (countp<4 && flagPrintTotalPg) {
    if (countp==0) fprintf(Out,"'G MW_{TOTAL}'");
    else           fprintf(Out,",'G MW_{TOTAL}'");
    countp++;
  }
  if (countp<4 && flagPrintTotalQg) {
    if (countp==0) fprintf(Out,"'G MVar_{TOTAL}'");
    else           fprintf(Out,",'G MVar_{TOTAL}'");
    countp++;
  }
  fprintf(Out,");\n");
  for(count=i=0,Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) {
    count++;
    if (Lptr->AC!=NULL) {
      if (!strcmp(Lptr->Type,"V") && Lptr->AC->Vlmax>Lptr->AC->Vlmin) {
        KV=Lptr->AC->KV;
        if (KV<=0) KV=1;
        if (i==0) {
          fprintf(Out,"hold;\nlmax=max(K*x(:,1)); lmin=min(K*x(:,1));\n");
          i=1;
        }
        fprintf(Out,"plot([lmin lmax],[%5.1lf %5.1lf],':');  ",KV*Lptr->AC->Vlmax,KV*Lptr->AC->Vlmax);
        fprintf(Out,"plot([lmin lmax],[%5.1lf %5.1lf],':');\n",KV*Lptr->AC->Vlmin,KV*Lptr->AC->Vlmin);
      }
      if (ExistParameter('e') && Lptr->AC->Gen!=NULL) count=count+3;
    }
  }
  fprintf(Out,"title('Profiles'); xlabel('L.F. [p.u.]'); \n");
  if (ExistParameter('O')) {
    fprintf(Out,"%s Define variables for AC/DC TEF profiles:\n","%%");
    fprintf(Out,"L=x(:,1); Vac=x(:,%1d);\n",count+2);
    for (count+=2,i=0,DCptr=dataPtr->DCbus;DCptr!=NULL;DCptr=DCptr->Next) if(!strcmp(DCptr->Type,"R")) {
      i++;  DCptrp=DCptr->To;
      if (i==1) fprintf(Out,"k=input('V=Vac+k*Vdc -> k=');\n");
      for(l=1,Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) if (Lptr->AC!=NULL) {
        if (Lptr->AC==DCptr->AC && !strcmp(Lptr->Type,"V")) j=l;
        if (Lptr->AC==DCptrp->AC && !strcmp(Lptr->Type,"V")) k=l;
        if (Lptr->AC==DCptr->AC && !strcmp(Lptr->Type,"D")) js=l;
        if (Lptr->AC==DCptrp->AC && !strcmp(Lptr->Type,"D")) ks=l;
        l++;
      }
      fprintf(Out,"Xc=[%lf %lf];\n",DCptr->Xc,DCptrp->Xc);
      fprintf(Out,"V=[x(:,%1d) x(:,%1d)];\n",j+1,k+1);
      fprintf(Out,"di=[x(:,%1d) x(:,%1d)];\n",js+1,ks+1);
      fprintf(Out,"a=[x(:,%1d) x(:,%1d)];\n",count,count+3);
      fprintf(Out,"b=[x(:,%1d) x(:,%1d)];\n",count+1,count+4);
      fprintf(Out,"d=[x(:,%1d) x(:,%1d)];\n",count+2,count+5);
      k1=3*sqrt(2.0)/PI/DCptr->Ld*DCptr->AC->V*DCptr->Tap*DCptr->Ntrf;
      k2=(3/PI*(DCptr->Xc-DCptrp->Xc)+DCptr->Rd)/DCptr->Ld;
      /* if (ExistParameter('O')) {
        printf("Input DC bus %8s PI controller gains:\n",DCptr->Name);
        for(;;) {
          printf("  Ki=");
          if (!scanf("%lf",&ki) || ki<=0) printf("Error in input data. Try again.\n");
          else break;
        }
        for(;;) {
          printf("  Kp=");
          if (!scanf("%lf",&kp) || kp<=0) printf("Error in input data. Try again.\n");
          else break;
        }
      } else */{ki=75.; kp=1.;}
      ki=ki*Sn/DCptr->Vn;
      kp=kp*Sn/DCptr->Vn;
      P22=0.5*(ki+k1)/(k1*(k1*kp+k2));
      if (ki) {
        P11=0.5/ki*((ki+k1)/(k1*kp+k2)+(k1*kp+k2)/k1);
        P12= -0.5/k1;
      }
      k1=3*sqrt(2.0)/PI/DCptrp->Ld*DCptrp->AC->V*DCptrp->Tap*DCptrp->Ntrf;
      /* if (ExistParameter('O')) {
       printf("Input DC bus %8s PI controller gains:\n",DCptr->Name);
        for(;;) {
          printf("  Ki=");
          if (!scanf("%lf",&ki) || ki<=0) printf("Error in input data. Try again.\n");
          else break;
        }
        for(;;) {
          printf("  Kp=");
          if (!scanf("%lf",&kp) || kp<=0) printf("Error in input data. Try again.\n");
          else break;
        }
      } else */ {ki=75.; kp=1.;}
      ki=ki*Sn/DCptr->Vn;
      kp=kp*Sn/DCptr->Vn;
      Beta=P22/(0.5*(ki+k1)/(k1*(k1*kp+k2)));
      if (ExistParameter('d')) fprintf(stderr,"Beta=%lf\n",Beta);
      if (ki) {
        P33=Beta*0.5/ki*((ki+k1)/(k1*kp+k2)+(k1*kp+k2)/k1);
        P23=Beta*0.5/k1;
      }
      fprintf(Out,"P=[%lf %lf %lf\n",P11,P12,0.);
      fprintf(Out,"   %lf %lf %lf\n",P12,P22,P23);
      fprintf(Out,"   %lf %lf %lf];\n",0.,P23,P33);
      fprintf(Out,"[Vac,Vdc%d]=addtotef(L,di,a,b,d,V,Vac,P,Xc);\n",i);
      if (i==1) fprintf(Out,"Vdc=Vdc%1d;\n",i);
      else fprintf(Out,"Vdc=Vdc+Vdc%1d;\n",i);
      count+=6;
    }
    if (i)  fprintf(Out,"[lambda,TEF]=tefprof(L,Vac+k*Vdc);\n");
    else fprintf(Out,"[lambda,TEF]=tefprof(L,Vac);\n");
  }
  if (Out!=NULL) fclose(Out);
}


/* --------------------------- TEFMatlabFiles --------------------------------- */
#ifdef ANSIPROTO
void TEFMatlabFiles(void)
#else
void TEFMatlabFiles()
#endif
/* Create the Matlab .m files needed to compute the TEF profiles. */
{
  FILE *Out;

  /* ----------------  Create 'tefprof.m' file needed for Matlab computations --------------- */
  if (ExistParameter('O')) {
    Out=OpenOutput("tefprof.m");
    fprintf(Out,"function [L1,V3]=tefprof(L,V)\n");
    fprintf(Out,"%s\n","%%");
    fprintf(Out,"%s  Plot TEF profiles using linear interpolation.\n","%%");
    fprintf(Out,"%s\n","%%");
    fprintf(Out,"[m,n]=max(L);\n");
    fprintf(Out,"l=length(L);\n");
    fprintf(Out,"L1=L(1:n);\n");
    fprintf(Out,"V1=V(1:n);\n");
    fprintf(Out,"L2=L(n:l);\n");
    fprintf(Out,"V2=V(n:l);\n");
    fprintf(Out,"L3(1)=L2(1);\n");
    fprintf(Out,"V3(1)=V2(1);\n");
    fprintf(Out,"J=1;\n");
    fprintf(Out,"for I=2:length(L2)\n");
    fprintf(Out,"  if(L2(I)<L2(I-1)),\n");
    fprintf(Out,"    J=J+1;\n");
    fprintf(Out,"    L3(J)=L2(I);\n");
    fprintf(Out,"    V3(J)=V2(I);\n");
    fprintf(Out,"  end\n");
    fprintf(Out,"end\n");
    fprintf(Out,"L2=L3';\n");
    fprintf(Out,"V2=V3';\n");
    fprintf(Out,"L3=flipud(L2);\n");
    fprintf(Out,"if (L3(1)>0), L3(1)=0; end\n");
    fprintf(Out,"V3=flipud(V2);\n");
    fprintf(Out,"V3=abs(interp1(L3,V3,L1)-V1);\n");
    fprintf(Out,"figure;\n");
    fprintf(Out,"plot(L1,V3); \n");
    fprintf(Out,"title('TEF profile'); xlabel('L.F. [p.u.]'); ylabel('TEF [p.u.]');\n");
    fclose(Out);
  /* ----------------  Create 'addtotef.m' file needed for Matlab computations --------------- */
    if (dataPtr->DCbus!=NULL){
      Out=OpenOutput("addtotef.m");
      fprintf(Out,"function [Vac,Vdc]=addtotef(L,di,a,b,d,V,Vac,P,Xc)\n");
      fprintf(Out,"%s\n","%%");
      fprintf(Out,"%s  Add coupling Qdc terms to AC Vac.\n","%%");
      fprintf(Out,"%s\n","%%");
      fprintf(Out,"\n");
      fprintf(Out,"x=b(:,1)./sqrt(d(:,1));     %c cosar\n");
      fprintf(Out,"x=[x sqrt(-a(:,1))./Xc(1)]; %c Id\n");
      fprintf(Out,"x=[x b(:,2)./sqrt(d(:,2))]; %c cosgi\n");
      fprintf(Out,"\n");
      fprintf(Out,"c=d(:,1)./(4.*Xc(1).^2.*x(:,2).^2).*(ones(d(:,1))-x(:,1).^2);\n");
      fprintf(Out,"c=[c d(:,2)./(4.*Xc(2).^2.*x(:,2).^2).*(ones(d(:,2))-x(:,3).^2)];\n");
      fprintf(Out,"\n");
      fprintf(Out,"\n");
      fprintf(Out,"[m,n]=max(L);\n");
      fprintf(Out,"l=length(L);\n");
      fprintf(Out,"Vdc=zeros(l,1);\n");
      fprintf(Out,"L1=L(1:n);\n");
      fprintf(Out,"di1=di(1:n,:);\n");
      fprintf(Out,"V1=V(1:n,:);\n");
      fprintf(Out,"x1=x(1:n,:);\n");
      fprintf(Out,"L2=L(n:l);\n");
      fprintf(Out,"di2=di(n:l,:);\n");
      fprintf(Out,"a2=a(n:l,:);\n");
      fprintf(Out,"b2=b(n:l,:);\n");
      fprintf(Out,"c2=c(n:l,:);\n");
      fprintf(Out,"d2=d(n:l,:);\n");
      fprintf(Out,"V2=V(n:l,:);\n");
      fprintf(Out,"x2=x(n:l,:);\n");
      fprintf(Out,"cons=pi/180;\n");
      fprintf(Out,"for I=1:length(L2),\n");
      fprintf(Out,"  N=n+I-1;\n");
      fprintf(Out,"  X=L2(I);\n");
      fprintf(Out,"  for J=2:length(L1),\n");
      fprintf(Out,"   X1=L1(J-1); X2=L1(J);\n");
      fprintf(Out,"   if (X>=X1).*(X<=X2),\n");
      fprintf(Out,"    for K=1:2,\n");
      fprintf(Out,"      Id=x2(I,2);\n");
      fprintf(Out,"      if (K==1), cosag=x2(I,1); s=1;\n");
      fprintf(Out,"      else cosag=x2(I,3); s=-1; end;\n");
      fprintf(Out,"      a=a2(I,K);\n");
      fprintf(Out,"      b=b2(I,K);\n");
      fprintf(Out,"      c=c2(I,K);\n");
      fprintf(Out,"      d=d2(I,K);\n");
      fprintf(Out,"      V=V2(I,K);\n");
      fprintf(Out,"      di=cons*di2(I,K);\n");
      fprintf(Out,"      t=sqrt(-d/8/a);\n");
      fprintf(Out,"      Vd=3*sqrt(2.0)/pi*t*V*cosag-3/pi*sqrt(-a);\n");
      fprintf(Out,"      Pi=-s*Vd*Id;\n");
      fprintf(Out,"      Y1=V1(J-1,K);  Y2=V1(J,K);\n");
      fprintf(Out,"      m=(Y1-Y2)/(X1-X2);\n");
      fprintf(Out,"      p=Y1-m*X1;\n");
      fprintf(Out,"      V0=m*X+p;\n");
      fprintf(Out,"      U=a+b*V+c*V*V;\n");
      fprintf(Out,"      U0=a+b*V0+c*V0*V0;\n");
      fprintf(Out,"      q=3*Id/pi*(sqrt(U)+0.5*b/sqrt(c)*log(2*sqrt(c*U)+2*c*V+b)-sqrt(-a)*asin((b*V+2*a)/V/sqrt(d)));\n");
      fprintf(Out,"      q0=3*Id/pi*(sqrt(U0)+0.5*b/sqrt(c)*log(2*sqrt(c*U0)+2*c*V0+b)-sqrt(-a)*asin((b*V0+2*a)/V0/sqrt(d)));\n");
      fprintf(Out,"      Y1=di1(J-1,K);  Y2=di1(J,K);\n");
      fprintf(Out,"      m=(Y1-Y2)/(X1-X2);\n");
      fprintf(Out,"      p=Y1-m*X1;\n");
      fprintf(Out,"      di0=cons*(m*X+p);\n");
      fprintf(Out,"      Vac(N)=Vac(N)-Pi*(di-di0)+q-q0;\n");
      fprintf(Out,"    end;\n");
      fprintf(Out,"    Y1=x1(J-1,1);  Y2=x1(J,1);\n");
      fprintf(Out,"    m=(Y1-Y2)/(X1-X2);\n");
      fprintf(Out,"    p=Y1-m*X1;\n");
      fprintf(Out,"    dx1=x2(I,1)-(m*X+p);\n");
      fprintf(Out,"    Y1=x1(J-1,2);  Y2=x1(J,2);\n");
      fprintf(Out,"    m=(Y1-Y2)/(X1-X2);\n");
      fprintf(Out,"    p=Y1-m*X1;\n");
      fprintf(Out,"    dx2=x2(I,2)-(m*X+p);\n");
      fprintf(Out,"    Y1=x1(J-1,3);  Y2=x1(J,3);\n");
      fprintf(Out,"    m=(Y1-Y2)/(X1-X2);\n");
      fprintf(Out,"    p=Y1-m*X1;\n");
      fprintf(Out,"    dx3=x2(I,3)-(m*X+p);\n");
      fprintf(Out,"    dx=[dx1 dx2 dx3]';\n");
      fprintf(Out,"    Vdc(N)=0.5*dx'*P*dx;\n");
      fprintf(Out,"    break;\n");
      fprintf(Out,"   end;\n");
      fprintf(Out,"  end;\n");
      fprintf(Out,"end;\n");
      fclose(Out);
    }
  }
}
