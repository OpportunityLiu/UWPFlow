#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

/* SVC and TCSC mismatch vector and Jacobian for Direct Method */

#include "pointl.h"


/* ------------------ SVCFunHes ----------------------------- */
#ifdef ANSIPROTO
void SVCFunHes(BOOLEAN flagF,BOOLEAN flagJ)
#else
void SVCFunHes(flagF,flagJ)
BOOLEAN flagF,flagJ;
#endif
/* Construct the SVC part of the PoC Jacobian and mismatch. */
{
 INDEX i,j,k,l,N;
 SVCbusData *SVCptr;
 VALUETYPE Vk,Xl,Xc,Bv,Xsl,alpha;
 BOOLEAN flag1;

 i=NacVar+11*Ndc/2;
 N=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;
 j=N+i;
 for(SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
   k=ACvar[SVCptr->From->N];
   Vk=SVCptr->From->V;
   l=ACvar[SVCptr->Ctrl->N];
   if (!strcmp(SVCptr->Cont,"AL")) flag1=FALSE;
   else flag1=TRUE;
   Xl=SVCptr->Xl;
   Xc=SVCptr->Xc;
   Xsl=SVCptr->slope;
   Bv=SVCptr->Bv;
   alpha=SVCptr->alpha_svc;
   if(flagF){
     dF[j+1]=x0[i+2]+x0[k+1];
     dF[j+2]=-Xsl*Vk*x0[i+1]+Vk*Vk*x0[i+2]-PI*Xl*x0[i+3];
     if(!flag1) dF[j+3]=2*(cos(2.0*alpha)-1.0)*x0[i+3];
     else       dF[j+3]=-x0[i+1];
     dF[k+1+N]=dF[k+1+N]-Xsl*Bv*x0[i+1]+2*(Bv-1.0/Xc)*Vk*x0[i+2];
     dF[l+1+N]=dF[l+1+N]+x0[i+1];
   }
   if(flagJ){
     JacElement(Jac,j+2,k+1,2.0*Vk*x0[i+2]-Xsl*x0[i+1]);
     if (!flag1) JacElement(Jac,j+3,i+3,-4.0*sin(2.0*alpha)*x0[i+3]);
     else JacElement(Jac,j+3,i+3,0.0);
     JacElement(Jac,k+1+N,k+1,2.0*(Bv-1.0/Xc)*x0[i+2]);
     JacElement(Jac,k+1+N,i+2,-Xsl*x0[i+1]+2.0*Vk*x0[i+2]);
   }
   i=i+3;
   j=j+3;
 }
}

/* ------------------ TCSCFunHes ----------------------------- */
#ifdef ANSIPROTO
void TCSCFunHes(BOOLEAN flagF,BOOLEAN flagJ)
#else
void TCSCFunHes(flagF,flagJ)
BOOLEAN flagF,flagJ;
#endif
/* Construct the TCSC part of the PoC Jacobian and mismatch. */
{
 INDEX i,j,k,m,N;
 TCSCbusData *TCSCptr;
 VALUETYPE Vk,Vm,thk,thm,Xc,Xl,Ptcsc,Qtcsck,Be,alpha;
 VALUETYPE Itcsc,Stcsc,D,sign,Kf;
 VALUETYPE s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11;
 BOOLEAN dVk=FALSE,dVm=FALSE;

 i=NacVar+11*Ndc/2+3*Nsvc;
 N=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;
 j=N+i;
 for(TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;TCSCptr=TCSCptr->Next){
   k=ACvar[TCSCptr->From->N];
   Vk=TCSCptr->From->V;
   thk=TCSCptr->From->Ang;
   m=ACvar[TCSCptr->To->N];
   Vm=TCSCptr->To->V;
   thm=TCSCptr->To->Ang;
   Xl=TCSCptr->Xl;
   Xc=TCSCptr->Xc;
   Kf=sqrt(Xc/Xl);
   Ptcsc=TCSCptr->Ptcsc;
   Qtcsck=TCSCptr->Qtcsck;
   Stcsc=sqrt(Ptcsc*Ptcsc+Qtcsck*Qtcsck);
   D=Stcsc*Stcsc*Stcsc;
   Be=TCSCptr->Be;
   alpha=TCSCptr->alpha_tcsc;
   Itcsc=TCSCptr->Itcsc;
   if (Itcsc>=0) sign=1.0; else sign=-1.0;
   if (TCSCptr->From->Cont!=NULL) dVk=TRUE;
   if (TCSCptr->To->Cont!=NULL) dVm=TRUE;
   if(flagF){
     dF[j+1]=-x0[i+1]+ sign*Ptcsc/Stcsc*x0[i+5]-x0[k]+x0[m];
     dF[j+2]=-x0[i+2]+ sign*Qtcsck/Stcsc*x0[i+5]-x0[k+1];
     dF[j+3]=-x0[i+3]-x0[m+1];
     dF[j+4]=-Vm*Vm*sin(thk-thm)*x0[i+1]+(-Vk*Vk+Vk*Vm*cos(thk-thm))*x0[i+2]
             +(-Vm*Vm+Vk*Vm*cos(thk-thm))*x0[i+3]-x0[i+4];
     /* dF[j+5]=-2.0/PI/Xl*(-1.0+cos(2.0*alpha))*x0[i+4];  */
     s1 = -PI*sin(Kf*(-PI+alpha))*Kf*(pow(Kf,4.0)-2.0*Kf*Kf+1.0)/Xc/(-PI*
          pow(Kf,4.0)*cos(Kf*(-PI+alpha))+PI*cos(Kf*(-PI+alpha))+2.0*pow(Kf,4.0)*
          alpha*cos(Kf*(-PI+alpha))-2.0*alpha*Kf*Kf*cos(Kf*(-PI+alpha))+
          pow(Kf,4.0)*sin(-2.0*PI+2.0*alpha)*cos(Kf*(-PI+alpha))-
          sin(-2.0*PI+2.0*alpha)*Kf*Kf*cos(Kf*(-PI+alpha))-4.0*Kf*Kf*Kf*
          pow(cos(-PI+alpha),2.0)*sin(Kf*(-PI+alpha))+4.0*Kf*Kf*cos(-PI+alpha)*
          sin(-PI+alpha)*cos(Kf*(-PI+alpha)));
     s3 = -PI*cos(Kf*(-PI+alpha));
     s5 = (pow(Kf,4.0)-2.0*Kf*Kf+1.0)/Xc;
     s7 = 1/(pow(-PI*pow(Kf,4.0)*cos(Kf*(-PI+alpha))+PI*cos(Kf*(-PI+alpha))
          +2.0*pow(Kf,4.0)*alpha*cos(Kf*(-PI+alpha))-2.0*alpha*Kf*Kf*
          cos(Kf*(-PI+alpha))+pow(Kf,4.0)*sin(-2.0*PI+2.0*alpha)*
          cos(Kf*(-PI+alpha))-sin(-2.0*PI+2.0*alpha)*Kf*Kf*cos(Kf*(-PI+alpha))
          -4.0*Kf*Kf*Kf*pow(cos(-PI+alpha),2.0)*sin(Kf*(-PI+alpha))
          +4.0*Kf*Kf*cos(-PI+alpha)*sin(-PI+alpha)*cos(Kf*(-PI+alpha)),2.0));
     s8 = PI*pow(Kf,5.0)*sin(Kf*(-PI+alpha))-PI*sin(Kf*(-PI+alpha))*Kf
          +2.0*pow(Kf,4.0)*cos(Kf*(-PI+alpha))-2.0*pow(Kf,5.0)*alpha*
          sin(Kf*(-PI+alpha))-2.0*Kf*Kf*cos(Kf*(-PI+alpha))+2.0*alpha*Kf*Kf*Kf*
          sin(Kf*(-PI+alpha))+2.0*pow(Kf,4.0)*cos(-2.0*PI+2.0*alpha)*
          cos(Kf*(-PI+alpha))-pow(Kf,5.0)*sin(-2.0*PI+2.0*alpha)*
          sin(Kf*(-PI+alpha))-2.0*cos(-2.0*PI+2.0*alpha)*Kf*Kf*
          cos(Kf*(-PI+alpha))+sin(-2.0*PI+2.0*alpha)*Kf*Kf*Kf*
          sin(Kf*(-PI+alpha))+4.0*Kf*Kf*Kf*cos(-PI+alpha)*sin(Kf*(-PI+alpha))*
          sin(-PI+alpha)-4.0*pow(Kf,4.0)*pow(cos(-PI+alpha),2.0)*
          cos(Kf*(-PI+alpha))-4.0*Kf*Kf*pow(sin(-PI+alpha),2.0)*
          cos(Kf*(-PI+alpha))+4.0*Kf*Kf*pow(cos(-PI+alpha),2.0)*
          cos(Kf*(-PI+alpha));
     s6 = s7*s8;
     s4 = s5*s6;
     s2 = s3*s4;
     dF[j+5]=(s1+s2)*x0[i+4];
     dF[j+6]=-Vk*x0[i+5];
     dF[j+7]=-x0[i+6];
     if(!strcmp(TCSCptr->Cont,"X"))dF[j+4]=dF[j+4]-x0[i+7];
     else if(!strcmp(TCSCptr->Cont,"P"))dF[j+1]=dF[j+1]-x0[i+7];
     else if(!strcmp(TCSCptr->Cont,"I"))dF[j+6]=dF[j+6]-x0[i+7];
     else if(!strcmp(TCSCptr->Cont,"D"))dF[j+7]=dF[j+6]-x0[i+7];
     if (!strpbrk(TCSCptr->From->Type,"S"))
       dF[k+N]=dF[k+N]-Vk*Vm*Be*cos(thk-thm)*x0[i+1]-Vk*Vm*Be*sin(thk-thm)*x0[i+2]
               -Vk*Vm*Be*sin(thk-thm)*x0[i+3]+x0[i+6];
     if (dVk)
       dF[k+N+1]=dF[k+N+1]-Vm*Be*sin(thk-thm)*x0[i+1]+(-2.0*Vk*Be+Vm*Be*cos(thk-thm))*x0[i+2]
                 +Vm*Be*cos(thk-thm)*x0[i+3]-Itcsc*x0[i+5];
     if (!strpbrk(TCSCptr->To->Type,"S"))
       dF[m+N]=dF[m+N]+Vk*Vm*Be*cos(thk-thm)*x0[i+1]+Vk*Vm*Be*sin(thk-thm)*x0[i+2]
               +Vk*Vm*Be*sin(thk-thm)*x0[i+3]-x0[i+6];
     if (dVm)
       dF[m+N+1]=dF[m+N+1]-Vk*Be*sin(thk-thm)*x0[i+1]+(-2.0*Vm*Be+Vk*Be*cos(thk-thm))*x0[i+3]
                 +Vk*Be*cos(thk-thm)*x0[i+2];
   }
   if(flagJ){
     JacElement(Jac,j+1,i+1,sign*(1.0/Stcsc-Ptcsc*Ptcsc/D)*x0[i+5]);
     JacElement(Jac,j+1,i+2,sign*(-Ptcsc*Qtcsck/D)*x0[i+5]);
     JacElement(Jac,j+2,i+1,sign*(-Qtcsck*Ptcsc/D)*x0[i+5]);
     JacElement(Jac,j+2,i+2,sign*(1.0/Stcsc-Qtcsck*Qtcsck/D)*x0[i+5]);
     if (dVk)
       JacElement(Jac,j+4,k+1,-Vm*sin(thk-thm)*x0[i+1]+(-2.0*Vk+Vm*cos(thk-thm))*x0[i+2]+Vm*cos(thk-thm)*x0[i+3]);
     else JacElement(Jac,j+4,k+1,0.0);
     if (dVm)
       JacElement(Jac,j+4,m+1,-Vk*sin(thk-thm)*x0[i+1]+(-2.0*Vm+Vk*cos(thk-thm))*x0[i+3]+Vk*cos(thk-thm)*x0[i+2]);
     else JacElement(Jac,j+4,m+1,0.0);
     if (!strpbrk(TCSCptr->From->Type,"S"))
       JacElement(Jac,j+4,k,-Vk*Vm*cos(thk-thm)*x0[i+1]-Vk*Vm*sin(thk-thm)*x0[i+2]-Vk*Vm*sin(thk-thm)*x0[i+3]);
     if (!strpbrk(TCSCptr->To->Type,"S"))
       JacElement(Jac,j+4,m,Vk*Vm*cos(thk-thm)*x0[i+1]+Vk*Vm*sin(thk-thm)*x0[i+2]+Vk*Vm*sin(thk-thm)*x0[i+3]);
     /* JacElement(Jac,j+5,i+5,4.0/PI/Xl*sin(2.0*alpha)*x0[i+4]);  */
     s2 = -PI*cos(Kf*(-PI+alpha))*Kf*Kf*(pow(Kf,4.0)-2.0*Kf*Kf+1.0)/Xc/(-PI*
          pow(Kf,4.0)*cos(Kf*(-PI+alpha))+PI*cos(Kf*(-PI+alpha))+2.0*pow(Kf,4.0)*alpha*
          cos(Kf*(-PI+alpha))-2.0*alpha*Kf*Kf*cos(Kf*(-PI+alpha))+pow(Kf,4.0)*sin(-2.0*PI
          +2.0*alpha)*cos(Kf*(-PI+alpha))-sin(-2.0*PI+2.0*alpha)*Kf*Kf*cos(Kf*(-PI+alpha)
          )-4.0*Kf*Kf*Kf*pow(cos(-PI+alpha),2.0)*sin(Kf*(-PI+alpha))+4.0*Kf*Kf*cos(-PI+
          alpha)*sin(-PI+alpha)*cos(Kf*(-PI+alpha)));
     s4 = 2.0*PI*sin(Kf*(-PI+alpha))*Kf;
     s6 = (pow(Kf,4.0)-2.0*Kf*Kf+1.0)/Xc;
     s8 = 1/(pow(-PI*pow(Kf,4.0)*cos(Kf*(-PI+alpha))+PI*cos(Kf*(-PI+alpha))+
          2.0*pow(Kf,4.0)*alpha*cos(Kf*(-PI+alpha))-2.0*alpha*Kf*Kf*cos(Kf*(-PI+alpha))+
          pow(Kf,4.0)*sin(-2.0*PI+2.0*alpha)*cos(Kf*(-PI+alpha))-sin(-2.0*PI+2.0*alpha)*
          Kf*Kf*cos(Kf*(-PI+alpha))-4.0*Kf*Kf*Kf*pow(cos(-PI+alpha),2.0)*sin(Kf*(-PI+
          alpha))+4.0*Kf*Kf*cos(-PI+alpha)*sin(-PI+alpha)*cos(Kf*(-PI+alpha)),2.0));
     s9 = PI*pow(Kf,5.0)*sin(Kf*(-PI+alpha))-PI*sin(Kf*(-PI+alpha))*Kf+2.0*
          pow(Kf,4.0)*cos(Kf*(-PI+alpha))-2.0*pow(Kf,5.0)*alpha*sin(Kf*(-PI+alpha))-2.0*
          Kf*Kf*cos(Kf*(-PI+alpha))+2.0*alpha*Kf*Kf*Kf*sin(Kf*(-PI+alpha))+2.0*pow(Kf,4.0
          )*cos(-2.0*PI+2.0*alpha)*cos(Kf*(-PI+alpha))-pow(Kf,5.0)*sin(-2.0*PI+2.0*alpha)
          *sin(Kf*(-PI+alpha))-2.0*cos(-2.0*PI+2.0*alpha)*Kf*Kf*cos(Kf*(-PI+alpha))+sin(
          -2.0*PI+2.0*alpha)*Kf*Kf*Kf*sin(Kf*(-PI+alpha))+4.0*Kf*Kf*Kf*cos(-PI+alpha)*sin
          (Kf*(-PI+alpha))*sin(-PI+alpha)-4.0*pow(Kf,4.0)*pow(cos(-PI+alpha),2.0)*cos(Kf*
          (-PI+alpha))-4.0*Kf*Kf*pow(sin(-PI+alpha),2.0)*cos(Kf*(-PI+alpha))+4.0*Kf*Kf*
          pow(cos(-PI+alpha),2.0)*cos(Kf*(-PI+alpha));
     s7 = s8*s9;
     s5 = s6*s7;
     s3 = s4*s5;
     s1 = s2+s3;
     s2 = s1;
     s5 = 2.0*PI*cos(Kf*(-PI+alpha));
     s7 = (pow(Kf,4.0)-2.0*Kf*Kf+1.0)/Xc;
     s9 = 1/(pow(-PI*pow(Kf,4.0)*cos(Kf*(-PI+alpha))+PI*cos(Kf*(-PI+alpha))+
          2.0*pow(Kf,4.0)*alpha*cos(Kf*(-PI+alpha))-2.0*alpha*Kf*Kf*cos(Kf*(-PI+alpha))+
          pow(Kf,4.0)*sin(-2.0*PI+2.0*alpha)*cos(Kf*(-PI+alpha))-sin(-2.0*PI+2.0*alpha)*
          Kf*Kf*cos(Kf*(-PI+alpha))-4.0*Kf*Kf*Kf*pow(cos(-PI+alpha),2.0)*sin(Kf*(-PI+
          alpha))+4.0*Kf*Kf*cos(-PI+alpha)*sin(-PI+alpha)*cos(Kf*(-PI+alpha)),3.0));
     s10 = pow(PI*pow(Kf,5.0)*sin(Kf*(-PI+alpha))-PI*sin(Kf*(-PI+alpha))*Kf+
           2.0*pow(Kf,4.0)*cos(Kf*(-PI+alpha))-2.0*pow(Kf,5.0)*alpha*sin(Kf*(-PI+alpha))
           -2.0*Kf*Kf*cos(Kf*(-PI+alpha))+2.0*alpha*Kf*Kf*Kf*sin(Kf*(-PI+alpha))+2.0*pow(
           Kf,4.0)*cos(-2.0*PI+2.0*alpha)*cos(Kf*(-PI+alpha))-pow(Kf,5.0)*sin(-2.0*PI+2.0*
           alpha)*sin(Kf*(-PI+alpha))-2.0*cos(-2.0*PI+2.0*alpha)*Kf*Kf*cos(Kf*(-PI+alpha))
           +sin(-2.0*PI+2.0*alpha)*Kf*Kf*Kf*sin(Kf*(-PI+alpha))+4.0*Kf*Kf*Kf*cos(-PI+alpha
           )*sin(Kf*(-PI+alpha))*sin(-PI+alpha)-4.0*pow(Kf,4.0)*pow(cos(-PI+alpha),2.0)*
           cos(Kf*(-PI+alpha))-4.0*Kf*Kf*pow(sin(-PI+alpha),2.0)*cos(Kf*(-PI+alpha))+4.0*
           Kf*Kf*pow(cos(-PI+alpha),2.0)*cos(Kf*(-PI+alpha)),2.0);
     s8 = s9*s10;
     s6 = s7*s8;
     s4 = s5*s6;
     s6 = -PI*cos(Kf*(-PI+alpha));
     s8 = (pow(Kf,4.0)-2.0*Kf*Kf+1.0)/Xc;
     s10 = 1/(pow(-PI*pow(Kf,4.0)*cos(Kf*(-PI+alpha))+PI*cos(Kf*(-PI+alpha))+
           2.0*pow(Kf,4.0)*alpha*cos(Kf*(-PI+alpha))-2.0*alpha*Kf*Kf*cos(Kf*(-PI+alpha))+
           pow(Kf,4.0)*sin(-2.0*PI+2.0*alpha)*cos(Kf*(-PI+alpha))-sin(-2.0*PI+2.0*alpha)*
           Kf*Kf*cos(Kf*(-PI+alpha))-4.0*Kf*Kf*Kf*pow(cos(-PI+alpha),2.0)*sin(Kf*(-PI+
           alpha))+4.0*Kf*Kf*cos(-PI+alpha)*sin(-PI+alpha)*cos(Kf*(-PI+alpha)),2.0));
     s11 = PI*pow(Kf,6.0)*cos(Kf*(-PI+alpha))-PI*cos(Kf*(-PI+alpha))*Kf*Kf-
           pow(Kf,6.0)*sin(-2.0*PI+2.0*alpha)*cos(Kf*(-PI+alpha))-2.0*pow(Kf,6.0)*alpha*
           cos(Kf*(-PI+alpha))+12.0*pow(Kf,4.0)*cos(-PI+alpha)*cos(Kf*(-PI+alpha))*sin(-PI
           +alpha)+4.0*sin(-2.0*PI+2.0*alpha)*Kf*Kf*cos(Kf*(-PI+alpha))-4.0*pow(Kf,5.0)*
           sin(Kf*(-PI+alpha))-4.0*pow(Kf,5.0)*cos(-2.0*PI+2.0*alpha)*sin(Kf*(-PI+alpha))+
           4.0*Kf*Kf*Kf*sin(Kf*(-PI+alpha))+4.0*cos(-2.0*PI+2.0*alpha)*Kf*Kf*Kf*sin(Kf*(-
           PI+alpha))+4.0*pow(Kf,5.0)*pow(cos(-PI+alpha),2.0)*sin(Kf*(-PI+alpha))+2.0*pow(
           Kf,4.0)*alpha*cos(Kf*(-PI+alpha))-3.0*pow(Kf,4.0)*sin(-2.0*PI+2.0*alpha)*cos(Kf
           *(-PI+alpha))-16.0*Kf*Kf*cos(-PI+alpha)*sin(-PI+alpha)*cos(Kf*(-PI+alpha));
     s9 = s10*s11;
     s7 = s8*s9;
     s5 = s6*s7;
     s3 = s4+s5;
     JacElement(Jac,j+5,i+5,(s2+s3)*x0[i+4]);
     if (dVk) JacElement(Jac,j+6,k+1,-x0[i+5]);
     else JacElement(Jac,j+6,k+1,0.0);
     if (!strpbrk(TCSCptr->From->Type,"S")) {
       JacElement(Jac,k+N,i+4,-Vk*Vm*cos(thk-thm)*x0[i+1]-Vk*Vm*sin(thk-thm)*x0[i+2]-Vk*Vm*sin(thk-thm)*x0[i+3]);
       if (dVk) JacElement(Jac,k+N,k+1,-Vm*Be*cos(thk-thm)*x0[i+1]-Vm*Be*sin(thk-thm)*x0[i+2]-Vm*Be*sin(thk-thm)*x0[i+3]);
       else JacElement(Jac,k+N,k+1,0.0);
       if (dVm) JacElement(Jac,k+N,m+1,-Vk*Be*cos(thk-thm)*x0[i+1]-Vk*Be*sin(thk-thm)*x0[i+2]-Vk*Be*sin(thk-thm)*x0[i+3]);
       else JacElement(Jac,k+N,m+1,0.0);
       JacElement(Jac,k+N,k,Vk*Vm*Be*(sin(thk-thm)*x0[i+1]-cos(thk-thm)*x0[i+2]-cos(thk-thm)*x0[i+3]));
       if (!strpbrk(TCSCptr->To->Type,"S"))
         JacElement(Jac,k+N,k,-Vk*Vm*Be*(sin(thk-thm)*x0[i+1]-cos(thk-thm)*x0[i+2]-cos(thk-thm)*x0[i+3]));
     }
     if (dVk) {
       JacElement(Jac,k+N+1,i+4,-Vm*sin(thk-thm)*x0[i+1]+(-2.0*Vk+Vm*cos(thk-thm))*x0[i+2]+Vm*cos(thk-thm)*x0[i+3]);
       JacElement(Jac,k+N+1,i+6,-x0[i+5]);
       JacElement(Jac,k+N+1,k+1,-2*Be*x0[i+2]);
       if (dVm) JacElement(Jac,k+N+1,m+1,-Be*sin(thk-thm)*x0[i+1]+Be*cos(thk-thm)*x0[i+2]+Be*cos(thk-thm)*x0[i+3]);
       else JacElement(Jac,k+N+1,m+1,0.0);
       if (!strpbrk(TCSCptr->From->Type,"S"))
         JacElement(Jac,k+N+1,k,-Vm*Be*cos(thk-thm)*x0[i+1]-Vm*Be*sin(thk-thm)*x0[i+2]-Vm*Be*sin(thk-thm)*x0[i+3]);
       if (!strpbrk(TCSCptr->To->Type,"S"))
         JacElement(Jac,k+N+1,m,Vm*Be*cos(thk-thm)*x0[i+1]+Vm*Be*sin(thk-thm)*x0[i+2]+Vm*Be*sin(thk-thm)*x0[i+3]);
     } else {
       JacElement(Jac,k+N+1,i+4,0.0);
       JacElement(Jac,k+N+1,i+6,0.0);
       JacElement(Jac,k+N+1,k+1,0.0);
       JacElement(Jac,k+N+1,m+1,0.0);
       if (!strpbrk(TCSCptr->From->Type,"S")) JacElement(Jac,k+N+1,k,0.0);
       if (!strpbrk(TCSCptr->To->Type,"S")) JacElement(Jac,k+N+1,m,0.0);
     }
     if (!strpbrk(TCSCptr->To->Type,"S")) {
       JacElement(Jac,m+N,i+4,Vk*Vm*cos(thk-thm)*x0[i+1]+Vk*Vm*sin(thk-thm)*x0[i+2]+Vk*Vm*sin(thk-thm)*x0[i+3]);
       if (dVk) JacElement(Jac,m+N,k+1,Vm*Be*cos(thk-thm)*x0[i+1]+Vm*Be*sin(thk-thm)*x0[i+2]+Vm*Be*sin(thk-thm)*x0[i+3]);
       else JacElement(Jac,m+N,k+1,0.0);
       if (dVm) JacElement(Jac,m+N,m+1,Vk*Be*cos(thk-thm)*x0[i+1]+Vk*Be*sin(thk-thm)*x0[i+2]+Vk*Be*sin(thk-thm)*x0[i+3]);
       else JacElement(Jac,m+N,m+1,0.0);
       if (!strpbrk(TCSCptr->From->Type,"S"))
         JacElement(Jac,m+N,k,-Vk*Vm*Be*(sin(thk-thm)*x0[i+1]-cos(thk-thm)*x0[i+2]-cos(thk-thm)*x0[i+3]));
       JacElement(Jac,m+N,k,Vk*Vm*Be*(sin(thk-thm)*x0[i+1]-cos(thk-thm)*x0[i+2]-cos(thk-thm)*x0[i+3]));
     }
     if (dVm) {
       JacElement(Jac,m+N+1,i+4,-Vk*sin(thk-thm)*x0[i+1]+(-2.0*Vm+Vk*cos(thk-thm))*x0[i+3]+Vk*cos(thk-thm)*x0[i+2]);
       if (dVk) JacElement(Jac,m+N+1,k+1,-Be*sin(thk-thm)*x0[i+1]+Be*cos(thk-thm)*x0[i+3]+Be*cos(thk-thm)*x0[i+2]);
       else JacElement(Jac,m+N+1,k+1,0.0);
       JacElement(Jac,m+N+1,m+1,-2.0*Be*x0[i+3]);
       if (!strpbrk(TCSCptr->From->Type,"S"))
         JacElement(Jac,m+N+1,k,-Vk*Be*cos(thk-thm)*x0[i+1]-Vk*Be*sin(thk-thm)*x0[i+2]-Vk*Be*sin(thk-thm)*x0[i+3]);
       if (!strpbrk(TCSCptr->To->Type,"S"))
         JacElement(Jac,m+N+1,k,Vk*Be*cos(thk-thm)*x0[i+1]+Vk*Be*sin(thk-thm)*x0[i+2]+Vk*Be*sin(thk-thm)*x0[i+3]);
     } else {
       JacElement(Jac,m+N+1,i+4,0.0);
       JacElement(Jac,m+N+1,k+1,0.0);
       JacElement(Jac,m+N+1,m+1,0.0);
       if (!strpbrk(TCSCptr->From->Type,"S")) JacElement(Jac,m+N+1,k,0.0);
       if (!strpbrk(TCSCptr->To->Type,"S")) JacElement(Jac,m+N+1,m,0.0);
     }
   }
   i=i+7;
   j=j+7;
 }
}


/* ------------------ STATCOMFunHes ----------------------------- */
#ifdef ANSIPROTO
void STATCOMFunHes(BOOLEAN flagF,BOOLEAN flagJ)
#else
void STATCOMFunHes(flagF,flagJ)
BOOLEAN flagF,flagJ;
#endif
/* Construct the STATCOM part of the PoC Jacobian and mismatch. */
{
 INDEX i,j,k,l,N;
 STATCOMbusData *STATCOMptr;
 VALUETYPE Vk,Xsl,delta,R,G,B,Gc,I,theta,Vdc,K,alpha,Q;
 BOOLEAN flagLimits,flagPWM;

 i=NacVar+11*Ndc/2+3*Nsvc+NtcscVar;
 N=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;
 j=N+i;
 for(STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
   k=ACvar[STATCOMptr->From->N];
   Vk=STATCOMptr->From->V;
   delta=STATCOMptr->From->Ang;
   l=ACvar[STATCOMptr->Ctrl->N];
   if (!strcmp(STATCOMptr->Cont,"PW") || !strcmp(STATCOMptr->Cont,"AL")) flagLimits=FALSE;
   else                                                                  flagLimits=TRUE;
   if (!strcmp(STATCOMptr->Cont1,"PW")) flagPWM=TRUE;
   else                                 flagPWM=FALSE;
   R=STATCOMptr->R;
   G=STATCOMptr->G;
   B=STATCOMptr->B;
   Gc=STATCOMptr->Gc;
   Xsl=STATCOMptr->slope;
   I=STATCOMptr->I;
   theta=STATCOMptr->theta;
   Vdc=STATCOMptr->Vdc;
   K=STATCOMptr->k;
   alpha=STATCOMptr->alpha;
   Q=STATCOMptr->Q;
   if(flagF){
     if (!flagLimits) {
        dF[j+1]= x0[i+3]*(-2.0*R*I)+ x0[i+4]*(-Vk*cos(delta-theta))+x0[i+5]*(-Vk*sin(delta-theta));
        if (Q>0) dF[j+1] += x0[i+1]*(-Xsl);
        else     dF[j+1] += x0[i+1]*Xsl;
     } else {
        dF[j+1] = -x0[i+1];
     }

     dF[j+2]= x0[i+4]*(-Vk*I*sin(delta-theta))+x0[i+5]*(Vk*I*cos(delta-theta));

     dF[j+3]= x0[i+3]*(-2.0*Gc*Vdc)
             +x0[i+6]*(G*K*Vk*cos(delta-alpha)+B*K*Vk*sin(delta-alpha))
             +x0[i+7]*(-B*K*Vk*cos(delta-alpha)+G*K*Vk*sin(delta-alpha));
     if (flagPWM) dF[j+3]+= x0[i+2];

     dF[j+4]= x0[i+6]*(G*Vdc*Vk*cos(delta-alpha)+B*Vdc*Vk*sin(delta-alpha))
             +x0[i+7]*(-B*Vdc*Vk*cos(delta-alpha)+G*Vdc*Vk*sin(delta-alpha));
     if (!flagPWM) dF[j+4]+= x0[i+2];

     dF[j+5]= x0[i+6]*(G*K*Vdc*Vk*sin(delta-alpha)-B*K*Vdc*Vk*cos(delta-alpha))
             +x0[i+7]*(-B*K*Vdc*Vk*sin(delta-alpha)-G*K*Vdc*Vk*cos(delta-alpha));

     dF[j+6]= x0[i+3]+x0[i+4]+x0[i+6]-x0[k];

     dF[j+7]= x0[i+5]+x0[i+7]-x0[k+1];

     dF[l+1+N]+= x0[i+1];

     dF[k+1+N]+= x0[i+4]*(-I*cos(delta-theta))+x0[i+5]*(-I*sin(delta-theta))
                +x0[i+6]*(-2.0*G*Vk+G*K*Vdc*cos(delta-alpha)+B*K*Vdc*sin(delta-alpha))
                +x0[i+7]*(2.0*B*Vk-B*K*Vdc*cos(delta-alpha)-G*K*Vdc*sin(delta-alpha));

     if (!strpbrk(STATCOMptr->From->Type,"S"))
        dF[k+N]+= x0[i+4]*(Vk*I*sin(delta-theta))+x0[i+5]*(-Vk*I*cos(delta-theta))
                 +x0[i+6]*(-G*K*Vdc*Vk*sin(delta-alpha)+B*K*Vdc*Vk*cos(delta-alpha))
                 +x0[i+7]*(B*K*Vdc*Vk*sin(delta-alpha)+G*K*Vdc*Vk*cos(delta-alpha));
   }
   if(flagJ){
     if (!flagLimits) {
        JacElement(Jac,j+1,i+1,x0[i+3]*(-2.0*R));
        JacElement(Jac,j+1,i+2,x0[i+4]*(-Vk*sin(delta-theta))+x0[i+5]*(Vk*cos(delta-theta)));
        JacElement(Jac,j+1,k+1,x0[i+4]*(-cos(delta-theta))+x0[i+5]*(-sin(delta-theta)));
        if (!strpbrk(STATCOMptr->From->Type,"S"))
           JacElement(Jac,j+1,k,x0[i+4]*(Vk*sin(delta-theta))+x0[i+5]*(-Vk*cos(delta-theta)));
     } else {
        JacElement(Jac,j+1,i+1,0.);
        JacElement(Jac,j+1,i+2,0.);
        JacElement(Jac,j+1,k+1,0.);
        if (!strpbrk(STATCOMptr->From->Type,"S")) JacElement(Jac,j+1,k,0.);
     }

     if (!flagLimits)
        JacElement(Jac,j+2,i+1,x0[i+4]*(-Vk*sin(delta-theta))+x0[i+5]*(Vk*cos(delta-theta)));
     else JacElement(Jac,j+2,i+1,0.);
     JacElement(Jac,j+2,i+2,x0[i+4]*(Vk*I*cos(delta-theta))+x0[i+5]*(Vk*I*sin(delta-theta)));
     JacElement(Jac,j+2,k+1,x0[i+4]*(-I*sin(delta-theta))+x0[i+5]*(I*cos(delta-theta)));
     if (!strpbrk(STATCOMptr->From->Type,"S"))
        JacElement(Jac,j+2,k,x0[i+4]*(-Vk*I*cos(delta-theta))+x0[i+5]*(-Vk*I*sin(delta-theta)));

     JacElement(Jac,j+3,i+3,x0[i+3]*(-2.0*Gc));
     JacElement(Jac,j+3,i+4,x0[i+6]*(G*Vk*cos(delta-alpha)+B*Vk*sin(delta-alpha))
                            +x0[i+7]*(-B*Vk*cos(delta-alpha)+G*Vk*sin(delta-alpha)));
     JacElement(Jac,j+3,i+5,x0[i+6]*(G*K*Vk*sin(delta-alpha)-B*K*Vk*cos(delta-alpha))
                            +x0[i+7]*(-B*K*Vk*sin(delta-alpha)-G*K*Vk*cos(delta-alpha)));
     JacElement(Jac,j+3,k+1,x0[i+6]*(G*K*cos(delta-alpha)+B*K*sin(delta-alpha))
                            +x0[i+7]*(-B*K*cos(delta-alpha)+G*K*sin(delta-alpha)));
     if (!strpbrk(STATCOMptr->From->Type,"S"))
        JacElement(Jac,j+3,k,x0[i+6]*(-G*K*Vk*sin(delta-alpha)+B*K*Vk*cos(delta-alpha))
                             +x0[i+7]*(B*K*Vk*sin(delta-alpha)+G*K*Vk*cos(delta-alpha)));

     JacElement(Jac,j+4,i+3,x0[i+6]*(G*Vk*cos(delta-alpha)+B*Vk*sin(delta-alpha))
                            +x0[i+7]*(-B*Vk*cos(delta-alpha)+G*Vk*sin(delta-alpha)));
     JacElement(Jac,j+4,i+5,x0[i+6]*(G*Vdc*Vk*sin(delta-alpha)-B*Vdc*Vk*cos(delta-alpha))
                            +x0[i+7]*(-B*Vdc*Vk*sin(delta-alpha)-G*Vdc*Vk*cos(delta-alpha)));
     JacElement(Jac,j+4,k+1,x0[i+6]*(G*Vdc*cos(delta-alpha)+B*Vdc*sin(delta-alpha))
                            +x0[i+7]*(-B*Vdc*cos(delta-alpha)+G*Vdc*sin(delta-alpha)));
     if (!strpbrk(STATCOMptr->From->Type,"S"))
        JacElement(Jac,j+4,k,x0[i+6]*(-G*Vdc*Vk*sin(delta-alpha)+B*Vdc*Vk*cos(delta-alpha))
                             +x0[i+7]*(B*Vdc*Vk*sin(delta-alpha)+G*Vdc*Vk*cos(delta-alpha)));

     JacElement(Jac,j+5,i+3,x0[i+6]*(G*K*Vk*sin(delta-alpha)-B*K*Vk*cos(delta-alpha))
                            +x0[i+7]*(-B*K*Vk*sin(delta-alpha)-G*K*Vk*cos(delta-alpha)));
     JacElement(Jac,j+5,i+4,x0[i+6]*(G*Vdc*Vk*sin(delta-alpha)-B*Vdc*Vk*cos(delta-alpha))
                            +x0[i+7]*(-B*Vdc*Vk*sin(delta-alpha)-G*Vdc*Vk*cos(delta-alpha)));
     JacElement(Jac,j+5,i+5,x0[i+6]*(-G*K*Vdc*Vk*cos(delta-alpha)-B*K*Vdc*Vk*sin(delta-alpha))
                            +x0[i+7]*(B*K*Vdc*Vk*cos(delta-alpha)-G*K*Vdc*Vk*sin(delta-alpha)));
     JacElement(Jac,j+5,k+1,x0[i+6]*(G*K*Vdc*sin(delta-alpha)-B*K*Vdc*cos(delta-alpha))
                            +x0[i+7]*(-B*K*Vdc*sin(delta-alpha)-G*K*Vdc*cos(delta-alpha)));
     if (!strpbrk(STATCOMptr->From->Type,"S"))
       JacElement(Jac,j+5,k,x0[i+6]*(G*K*Vdc*Vk*cos(delta-alpha)+B*K*Vdc*Vk*sin(delta-alpha))
                            +x0[i+7]*(-B*K*Vdc*Vk*cos(delta-alpha)+G*K*Vdc*Vk*sin(delta-alpha)));

     if (!flagLimits)
       JacElement(Jac,k+1+N,i+1,x0[i+4]*(-cos(delta-theta))+x0[i+5]*(-sin(delta-theta)));
     else JacElement(Jac,k+1+N,i+1,0.);
     JacElement(Jac,k+1+N,i+2,x0[i+4]*(-I*sin(delta-theta))+x0[i+5]*(I*cos(delta-theta)));
     JacElement(Jac,k+1+N,i+3,x0[i+6]*(G*K*cos(delta-alpha)+B*K*sin(delta-alpha))
                              +x0[i+7]*(-B*K*cos(delta-alpha)-G*K*sin(delta-alpha)));
     JacElement(Jac,k+1+N,i+4,x0[i+6]*(G*Vdc*cos(delta-alpha)+B*Vdc*sin(delta-alpha))
                              +x0[i+7]*(-B*Vdc*cos(delta-alpha)-G*Vdc*sin(delta-alpha)));
     JacElement(Jac,k+1+N,i+5,x0[i+6]*(G*K*Vdc*sin(delta-alpha)-B*K*Vdc*cos(delta-alpha))
                              +x0[i+7]*(-B*K*Vdc*sin(delta-alpha)+G*K*Vdc*cos(delta-alpha)));
     JacElement(Jac,k+1+N,k+1,x0[i+6]*(-2.0*G)+x0[i+7]*(2.0*B));
     if (!strpbrk(STATCOMptr->From->Type,"S"))
       JacElement(Jac,k+1+N,k,x0[i+4]*(I*sin(delta-theta))+x0[i+5]*(-I*cos(delta-theta))
                              +x0[i+6]*(-G*K*Vdc*sin(delta-alpha)+B*K*Vdc*cos(delta-alpha))
                              +x0[i+7]*(B*K*Vdc*sin(delta-alpha)-G*K*Vdc*cos(delta-alpha)));


     if (!strpbrk(STATCOMptr->From->Type,"S")) {
        if (!flagLimits)
           JacElement(Jac,k+N,i+1,x0[i+4]*(Vk*sin(delta-theta))+x0[i+5]*(-Vk*cos(delta-theta)));
        else JacElement(Jac,k+N,i+1,0.);
        JacElement(Jac,k+N,i+2,x0[i+4]*(-Vk*I*cos(delta-theta))+x0[i+5]*(-Vk*I*sin(delta-theta)));
        JacElement(Jac,k+N,i+3,x0[i+6]*(-G*K*Vk*sin(delta-alpha)+B*K*Vk*cos(delta-alpha))
                               +x0[i+7]*(B*K*Vk*sin(delta-alpha)+G*K*Vk*cos(delta-alpha)));
        JacElement(Jac,k+N,i+4,x0[i+6]*(-G*Vdc*Vk*sin(delta-alpha)+B*Vdc*Vk*cos(delta-alpha))
                               +x0[i+7]*(B*Vdc*Vk*sin(delta-alpha)+G*Vdc*Vk*cos(delta-alpha)));
        JacElement(Jac,k+N,i+5,x0[i+6]*(G*K*Vdc*Vk*cos(delta-alpha)+B*K*Vdc*Vk*sin(delta-alpha))
                               +x0[i+7]*(-B*K*Vdc*Vk*cos(delta-alpha)+G*K*Vdc*Vk*sin(delta-alpha)));
        JacElement(Jac,k+N,k+1,x0[i+4]*(I*sin(delta-theta))+x0[i+5]*(-I*cos(delta-theta))
                               +x0[i+6]*(-G*K*Vdc*sin(delta-alpha)+B*K*Vdc*cos(delta-alpha))
                               +x0[i+7]*(B*K*Vdc*sin(delta-alpha)+G*K*Vdc*cos(delta-alpha)));
        JacElement(Jac,k+N,k,x0[i+4]*(Vk*I*cos(delta-theta))+x0[i+5]*(Vk*I*sin(delta-theta))
                             +x0[i+6]*(-G*K*Vdc*Vk*cos(delta-alpha)-B*K*Vdc*Vk*sin(delta-alpha))
                             +x0[i+7]*(B*K*Vdc*Vk*cos(delta-alpha)-G*K*Vdc*Vk*sin(delta-alpha)));
     }
   }
   i=i+7;
   j=j+7;
 }
}

