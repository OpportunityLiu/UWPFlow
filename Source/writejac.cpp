/* Write Jacobian. */

#include "write.h"


/* ---------- Global Variables --------- */
INDEX TFnum,TFbus;
char TFname[13];
extern VALUETYPE K;
extern INDEX *ACvar;

/* --------------- WriteJac ---------------------- */
#ifdef ANSIPROTO
void WriteJac(void)
#else
void WriteJac()
#endif
{
  SparseMatrixElement *Jptr;
  ACbusData *ACptr,*ACptrp;
  AClist *Lptr;
  DCbusData *DCptrR,*DCptrI,*DCptr;
  ElementData *Eptr;
  ElementList *ELptr;
  char Namebase[80],Name[80],str[80],type[2];
  FILE *OutFile,*OutFilep;
  int i,j,k,l;
  INDEX I,J,N,Nvar;
  BOOLEAN flag=FALSE;
  SVCbusData *SVCptr;                  /* FACTS */
  TCSCbusData *TCSCptr;                /* FACTS */
  STATCOMbusData *STATCOMptr;          /* FACTS */

  if (ExistParameter('J')) strcpy_s(Namebase,NameParameter('J'));
  else strcpy_s(Namebase,NameParameter('j'));
  if(NullName(Namebase)) return;
  if (ExistParameter('J') && Bl)
    for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
      if(ACptr->Num==Bl) {
	strcpy_s(ACptr->Type,"B");
	if(ACptr->Area!=NULL && ACptr->Area->Slack==ACptr) strcat_s(ACptr->Type,"A");
	Bl=0;
	ACptr->Cont=ACptr;
	break;
      }
    }
  flag=(!ExistParameter('J')) && flagPoC;
  DeleteJac(Jac,NewRow,NewCol,OldRow,OldCol);
  RowPer=NewRow; ColPer=NewCol;
  Nvar=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;  /*  FACTS  */ 
  ACFunJac(Jac,NULL,TRUE,TRUE,FALSE);
  DCFunJac(Jac,TRUE,TRUE);
  SVCFunJac(Jac,TRUE,TRUE);                  /*  FACTS  */
  TCSCFunJac(Jac,TRUE,TRUE);                 /*  FACTS  */
  STATCOMFunJac(Jac,TRUE,TRUE);              /*  FACTS  */
  if(flagH) { Nvar++; HFunJac(TRUE,TRUE,NULL,Dx);}
  else if (flag) {
    Nvar=2*Nvar+1;
    ACFunHes(TRUE,TRUE);
    DCFunHes(TRUE,TRUE);
    SVCFunHes(TRUE,TRUE);                   /* FACTS  */
    TCSCFunHes(TRUE,TRUE);                  /* FACTS  */
    STATCOMFunHes(TRUE,TRUE);               /* FACTS  */
  }
  SortRowsColumns(Jac);
  strcpy_s(Name,Namebase);
  strcat_s(Name,".jac");
  OutFile=OpenOutput(Name);
  fprintf(OutFile,"%d %d\n",Nvar,Nvar);
  for(i=1;i<=Nvar;i++) {
    for(Jptr=Jac->RowHead[i];Jptr!=NULL;Jptr=Jptr->RowNext)
      fprintf(OutFile,"%4d %4d %-11.5g\n",Jptr->Row,Jptr->Col,Jptr->Value);
  }
  fprintf(OutFile,"%4d %4d %-11.5g\n",0,0,0.);
  fclose(OutFile);
  strcpy_s(Name,Namebase);
  strcat_s(Name,".var");
  OutFile=OpenOutput(Name);
  strcpy_s(Name,Namebase);
  strcat_s(Name,".mis");
  OutFilep=OpenOutput(Name);
  fprintf(OutFile,"%d 1\n",Nvar);
  fprintf(OutFilep,"%d 1\n",Nvar);
  for (i=0,ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
    if (ACptr->Cont!=NULL) {
      if (strpbrk(ACptr->Type,"S")) {
        sprintf_s(str,"kg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Kg);
      }
      else {
        sprintf_s(str,"d%-d",ACptr->Num);  fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Ang);
      }
      sprintf_s(str,"dP%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"V%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->V);
      sprintf_s(str,"dQ%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
    }
    else if(QRcont && strpbrk(ACptr->Type,"C")){
      sprintf_s(str,"d%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Ang);
      sprintf_s(str,"dP%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      for(Lptr=ACptr->ContBus;Lptr!=NULL;Lptr=Lptr->Next){
        ACptrp=Lptr->AC;
        if (strpbrk(ACptrp->cont,"V")) break;
      }
      sprintf_s(str,"Qr%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Qr);
      sprintf_s(str,"dQ%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
    }
    else if(Rcont && strpbrk(ACptr->Type,"T")){
      sprintf_s(str,"d%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Ang);
      sprintf_s(str,"dP%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      for(ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next){
	     Eptr=ELptr->Eptr;
	     I=Eptr->From->Num;
	     J=Eptr->To->Num;
	     if(!strcmp(Eptr->Type,"R")) break;
      }
      sprintf_s(str,"1/t%-d_%-d",I,J);
      fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,Eptr->Tap);
      sprintf_s(str,"dQ%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
    }
    else if(strpbrk(ACptr->Type,"L")) {
      sprintf_s(str,"d%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Ang);
      sprintf_s(str,"dP%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      fprintf(OutFile,"%4d %8s %-11.5g\n",++i,"l",lambda);
      sprintf_s(str,"dQ%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
    }
    else if(strpbrk(ACptr->Type,"Q") || strpbrk(ACptr->Type,"V") || (!QRcont && strpbrk(ACptr->Type,"G"))) {
      if (strpbrk(ACptr->Type,"S")) {
        sprintf_s(str,"kg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Kg);
      }
      else {
        sprintf_s(str,"d%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Ang);
      }
      sprintf_s(str,"dP%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Qg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Qg);
      sprintf_s(str,"dQ%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
    }
    else if(strpbrk(ACptr->Type,"Z")) {
      sprintf_s(str,"d%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Ang);
      sprintf_s(str,"dP%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Qz%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Qg);
      sprintf_s(str,"dQ%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
    }
    else if(strpbrk(ACptr->Type,"S")){
      sprintf_s(str,"kg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Kg);
      sprintf_s(str,"dP%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Qg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Qg);
      sprintf_s(str,"dQ%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
    }
    if(Acont && strpbrk(ACptr->Type,"A")){
      sprintf_s(str,"kg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Kg);
      sprintf_s(str,"dPA%-d",ACptr->Area->N); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
    }
    if (PQcont) for(ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next) {
      Eptr=ELptr->Eptr;
      if(strpbrk(Eptr->Type,"PQNM")) {
      	 if (Eptr->From==ACptr) {
	         I=Eptr->From->Num;
	         J=Eptr->To->Num;
      	 } else {
	         J=Eptr->From->Num;
	         I=Eptr->To->Num;
	       }
	       if(!strcmp(Eptr->Type,"RP")){
	         sprintf_s(str,"a%-d_%-d",I,J); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,Eptr->Ang);
	         sprintf_s(str,"dP%-d_%-d",I,J); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
	       }
	       else if(strpbrk(Eptr->Type,"PM")){
	         sprintf_s(str,"P%-d_%-d",I,J); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,Eptr->Cvar);
	         sprintf_s(str,"dP%-d_%-d",I,J); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      	 }
	       else if(!strcmp(Eptr->Type,"RQ")){
	         sprintf_s(str,"1/t%-d_%-d",I,J); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,Eptr->Tap);
	         sprintf_s(str,"dQ%-d_%-d",I,J); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
	       }
	       else {
	         sprintf_s(str,"Q%-d_%-d",I,J); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,Eptr->Cvar);
	         sprintf_s(str,"dQ%-d_%-d",I,J); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
	       }
      }
    }
    if (ACptr->Gen!=NULL) {
      i=ACptr->Gen->Nvar;
      sprintf_s(str,"dPg%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dQg%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dEq%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dEd%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dVd%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dVq%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dId%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dIq%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dVr%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dVi%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      sprintf_s(str,"dIa%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",++i,str,dF[i]);
      i=ACptr->Gen->Nvar;
      if (strpbrk(ACptr->cont,"E")) {
        sprintf_s(str,"Qg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Qg);
      } else {
        sprintf_s(str,"Eq%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Eq);
      }
      sprintf_s(str,"dg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->dg);
      sprintf_s(str,"Vr%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Vr);
      sprintf_s(str,"Vi%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Vi);
      sprintf_s(str,"Ir%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Ir);
      sprintf_s(str,"Ii%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Ii);
      sprintf_s(str,"Vq%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Vq);
      sprintf_s(str,"Vd%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Vd);
      sprintf_s(str,"Iq%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Iq);
      sprintf_s(str,"Id%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Id);
      if (strpbrk(ACptr->cont,"I")) {
	       sprintf_s(str,"Qg%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Qg);
      } else {
	       sprintf_s(str,"Ia%-d",ACptr->Num); fprintf(OutFile,"%4d %8s %-11.5g\n",++i,str,ACptr->Gen->Ia);
      }
    }
  }
  for(k=0,DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if(!strcmp(DCptrR->Type,"R")){
      for (k++,l=1;l<=11;l++){
	       sprintf_s(str,"Fdc%-d_%-d",k,l); i++;
       	fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      }
      for(l=i-11,j=1;j<=2;j++){
	       if(j==1) { DCptr=DCptrR; strcpy_s(type,"r"); }
	       else { DCptr=DCptrI; strcpy_s(type,"i"); }
	       if(strcmp(DCptr->Cont1,"VD")&&strcmp(DCptr->Cont2,"VD")) {
	         sprintf_s(str,"Vd%1s%-d",type,k);
       	  fprintf(OutFile,"%4d %8s %-11.5g\n",++l,str,DCptr->Vd);
       	}
       	if(strcmp(DCptr->Cont1,"AT")&&strcmp(DCptr->Cont2,"AT")) {
       	  sprintf_s(str,"t%1s%-d",type,k);
       	  fprintf(OutFile,"%4d %8s %-11.5g\n",++l,str,DCptr->Tap);
       	}
       	if(strcmp(DCptr->Cont1,"AL")&&strcmp(DCptr->Cont2,"AL")) {
       	  sprintf_s(str,"al%1s%-d",type,k);
       	  fprintf(OutFile,"%4d %8s %-11.5g\n",++l,str,DCptr->Alfa);
       	}
       	if(strcmp(DCptr->Cont1,"GA")&&strcmp(DCptr->Cont2,"GA")) {
       	  sprintf_s(str,"ga%1s%-d",type,k);
	         fprintf(OutFile,"%4d %8s %-11.5g\n",++l,str,DCptr->Gamma);
       	}
       	sprintf_s(str,"s%1s%-d",type,k);
       	fprintf(OutFile,"%4d %8s %-11.5g\n",++l,str,DCptr->MVA);
       	if(strcmp(DCptr->Cont1,"PA")&&strcmp(DCptr->Cont2,"PA")) {
       	  sprintf_s(str,"P%1s%-d",type,k);
       	  fprintf(OutFile,"%4d %8s %-11.5g\n",++l,str,DCptr->P);
       	}
       	if(strcmp(DCptr->Cont1,"QA")&&strcmp(DCptr->Cont2,"QA")) {
       	  sprintf_s(str,"Q%1s%-d",type,k);
       	  fprintf(OutFile,"%4d %8s %-11.5g\n",++l,str,DCptr->Q);
       	}
      }
      if(strcmp(DCptrR->Cont1,"ID")&&strcmp(DCptrR->Cont2,"ID")&&
       	 strcmp(DCptrI->Cont1,"ID")&&strcmp(DCptrI->Cont2,"ID")) {
        sprintf_s(str,"Id%-d",k);
        fprintf(OutFile,"%4d %8s %-11.5g\n",++l,str,DCptrR->Id);
      }
    }
  }

                               /*   FACTS   */
  for(k=0,SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
      k++; l=0;
      sprintf_s(str,"Qsvc%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,SVCptr->Qsvc);
      sprintf_s(str,"Fsvc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Bv%-d",k);i++;       fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,SVCptr->Bv);
      sprintf_s(str,"Fsvc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      if(!strcmp(SVCptr->Cont,"AL")){
        sprintf_s(str,"alpha%-d",k);i++;
        fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,SVCptr->alpha_svc);
      } else {
	       sprintf_s(str,"Vrefc%-d",k);i++;
        fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,SVCptr->Vvar);
      }
      sprintf_s(str,"Fsvc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
  }

  for(k=0,TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;TCSCptr=TCSCptr->Next){
      k++; l=0;
      sprintf_s(str,"Ptcsc%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,TCSCptr->Ptcsc);
      sprintf_s(str,"Ftcsc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Qtcsck%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,TCSCptr->Qtcsck);
      sprintf_s(str,"Ftcsc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Qtcscm%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,TCSCptr->Qtcscm);
      sprintf_s(str,"Ftcsc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Be%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,TCSCptr->Be);
      sprintf_s(str,"Ftcsc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"alpha%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,TCSCptr->alpha_tcsc);
      sprintf_s(str,"Ftcsc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Itcsc%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,TCSCptr->Itcsc);
      sprintf_s(str,"Ftcsc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"delta%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,TCSCptr->delta_t);
      sprintf_s(str,"Ftcsc%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
  }
  for(k=0,STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
      k++; l=0;
      if(!strcmp(STATCOMptr->Cont,"PW") || !strcmp(STATCOMptr->Cont,"AL")){
        sprintf_s(str,"Istat%-d",k);i++;
        fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,STATCOMptr->I);
      } else {
	       sprintf_s(str,"Vrefc%-d",k);i++;
        fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,STATCOMptr->Vvar);
      }
      sprintf_s(str,"Vcont%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"theta%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,STATCOMptr->theta);
      sprintf_s(str,"Cont%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Vdc%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,STATCOMptr->Vdc);
      sprintf_s(str,"Loss%-d_%-d",k,++l);  fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"k%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,STATCOMptr->k);
      sprintf_s(str,"ReS%-d_%-d",k,++l);   fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"alpha%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,STATCOMptr->alpha);
      sprintf_s(str,"ImS%-d_%-d",k,++l);   fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Pstat%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,STATCOMptr->P);
      sprintf_s(str,"Pstat%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      sprintf_s(str,"Qstat%-d",k);i++;
      fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,STATCOMptr->Q);
      sprintf_s(str,"Qstat%-d_%-d",k,++l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
  }
                             /* END FACTS */

  if (flagH) {
    fprintf(OutFile,"%4d %8s %-11.5g\n",Jac->n1,"l",lambda);
    fprintf(OutFilep,"%4d %8s %-11.5g\n",Jac->n1,"dH",dF[Jac->n1]);
  }
  else if (flag) {
    N=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;      /* FACTS */
    for (i=N,ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
      i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
      sprintf_s(str,"gP%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
      sprintf_s(str,"gQ%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      if(Acont && strpbrk(ACptr->Type,"A")){
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gPA%-d",ACptr->Area->N); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      }
      if (PQcont) for(ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next) {
	       Eptr=ELptr->Eptr;
        if(strpbrk(Eptr->Type,"PQNM")) {
          if (Eptr->From==ACptr) {
            I=Eptr->From->Num;
            J=Eptr->To->Num;
          } else {
            J=Eptr->From->Num;
            I=Eptr->To->Num;
          }
          if(!strcmp(Eptr->Type,"PM")){
            i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
            sprintf_s(str,"gP%-d_%-d",I,J); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
          } else {
            i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
            sprintf_s(str,"gQ%-d_%-d",I,J); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
          }
        }
      }
      if (ACptr->Gen!=NULL) {
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gPg%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gQg%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gEq%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gEd%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gVd%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gVq%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gId%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gIq%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gVr%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gVi%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gIa%-d",ACptr->Num); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      }
    }
    for(k=0,DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
      if(!strcmp(DCptrR->Type,"R")){
        for (k++,l=1;l<=11;l++){
          i++; sprintf_s(str,"w%-d",i-N); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
          sprintf_s(str,"gdc%-d_%-d",k,l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
        }
      }
    }

                               /* FACTS */
    for(k=0,SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
      for (k++,l=1;l<=3;l++){
        i++; sprintf_s(str,"wsvc%-d_%-d",k,l); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gsvc%-d_%-d",k,l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      }
    }
    for(k=0,TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;TCSCptr=TCSCptr->Next){
      for (k++,l=1;l<=7;l++){
        i++; sprintf_s(str,"wtcsc%-d_%-d",k,l); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gtcsc%-d_%-d",k,l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      }
    }
    for(k=0,STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
      for (k++,l=1;l<=7;l++){
        i++; sprintf_s(str,"wstat%-d_%-d",k,l); fprintf(OutFile,"%4d %8s %-11.5g\n",i,str,x0[i-N]);
        sprintf_s(str,"gstat%-d_%-d",k,l); fprintf(OutFilep,"%4d %8s %-11.5g\n",i,str,dF[i]);
      }
    }
				                       /* END OF FACTS */

    fprintf(OutFile,"%4d %8s %-11.5g\n",Jac->n1,"l",lambda);
    fprintf(OutFilep,"%4d %8s %-11.5g\n",Jac->n1,"gl",dF[Jac->n1]);
  }
  fprintf(OutFile,"%4d %8s %-11.5g\n",0,"0",0.);
  fprintf(OutFilep,"%4d %8s %-11.5g\n",0,"0",0.);
  fclose(OutFile);
  fclose(OutFilep);
}

/* --------------- WriteQmatrix ---------------------- */
#ifdef ANSIPROTO
void WriteQmatrix(INDEX count,VALUETYPE *vec)
#else
void WriteQmatrix(count,vec)
INDEX count;
VALUETYPE *vec;
#endif
{
  SparseMatrixElement *Jptr;
  ACbusData *ACptr;
  char Namebase[80],Name[80],MaxName[13];
  FILE *OutFile;
  INDEX i,j,I,J,m,MaxBus,MaxNum,Nvar;
  VALUETYPE val=0,valMax=-0.1;
  BOOLEAN flagTFbus=FALSE;

  strcpy_s(Namebase,NameParameter('0'));
  if(NullName(Namebase)) return;
  sprintf_s(Name,"%s%d.m",Namebase,count);
  OutFile=OpenOutput(Name);
  if (count==1) TFbus=IntegerParameter('1',0,1,9999);
  Nvar=NacVar+11*Ndc/2+3*Nsvc+NtcscVar+7*Nstatcom;     /* FACTS */
  for(i=1;i<=Nvar+1;i++) for(Jptr=Jac->RowHead[i];Jptr!=NULL;Jptr->Value=0,Jptr=Jptr->RowNext);
  ACFunJac(Jac,NULL,FALSE,TRUE,FALSE);
  DCFunJac(Jac,FALSE,TRUE);
  SVCFunJac(Jac,FALSE,TRUE);                  /*  FACTS  */
  TCSCFunJac(Jac,FALSE,TRUE);                 /*  FACTS  */
  STATCOMFunJac(Jac,FALSE,TRUE);              /*  FACTS  */
  fprintf(OutFile,"lambda(%d)=%-10.6g;\n",count,lambda+lambda_o);
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Full System Jacobian:\n","%%");
  fprintf(OutFile,"N=%d;\nJ=zeros(N);\n",Nvar);
  for(i=1;i<=Nvar+1;i++) {
    for(Jptr=Jac->RowHead[i];Jptr!=NULL;Jptr=Jptr->RowNext) {
      I=Jptr->Row;
      if (OldRow->p[I]!=0) I=OldRow->p[I];
      J=Jptr->Col;
      if (OldCol->p[J]!=0) J=OldCol->p[J];
      if (I<=Nvar && J<=Nvar) fprintf(OutFile,"J(%d,%d)=%-10.10g;\n",I,J,Jptr->Value);
    }
  }
  fprintf(OutFile,"J=sparse(J);\n");
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Voltage (load) Buses:\n","%%");
  fprintf(OutFile,"Vars=[(1:N)' zeros(N,1) zeros(N,1) ones(N,1)];\n");
  fprintf(OutFile,"clear Buses Vbuses\n");
  for (m=0,ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
    i=ACvar[ACptr->N];
    j=i+1;
    fprintf(OutFile,"Vars(%d,2)=%d; Vars(%d,2)=%d;\n",i,-ACptr->Num,j,ACptr->Num);
    if (!strpbrk(ACptr->Type,"S")) fprintf(OutFile,"Vars(%d,3)=1;\n",i);
    if ((ACptr->Cont!=NULL &&(QRcont || !strpbrk(ACptr->Type,"G"))) ||
       	(!Rcont && strpbrk(ACptr->Type,"T"))||
       	(!QRcont && strpbrk(ACptr->Type,"C"))){
      i++; m++;
      fprintf(OutFile,"Vbuses(%d,1)=%d; Vbuses(%d,2)=%d;\n",m,i,m,ACptr->Num);
      fprintf(OutFile,"Vars(%d,3)=1; Vars(%d,4)=0;\n",i,i);
      if (count==1) {
       	val=fabs(vec[i]);
       	if (TFbus==ACptr->Num && !flagTFbus) {
          flagTFbus=TRUE;
          TFnum=i;
          strcpy_s(TFname,ACptr->Name);
       	}
       	if (val>valMax) {
       	  valMax=val;
       	  MaxBus=ACptr->Num;
       	  MaxNum=i;
       	  strcpy_s(MaxName,ACptr->Name);
       	}
      }
    }
  }
  if (count==1 && !flagTFbus) {
    TFnum=MaxNum;
    TFbus=MaxBus;
    strcpy_s(TFname,MaxName);
  }
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Bus for test functions and tangent vector:\n","%%");
  fprintf(OutFile,"l=%d;\n",TFbus);
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Minimum singular value and real |e-value| for full J:\n","%%");
  fprintf(OutFile,"[val,svecJ]=inviter(J'*J); svJ(%d)=sqrt(val);\n",count);
  fprintf(OutFile,"[val,evecJ]=inviter(J); evJ(%d)=abs(val);\n",count);
  fprintf(OutFile,"  %s Critical bus numbers and ranking of bus l:\n","%%");
  fprintf(OutFile,"  [val,max_sv]=max(abs(svecJ));\n");
  fprintf(OutFile,"  crsvJ(%d,1)=Vars(max_sv,2);\n",count);
  fprintf(OutFile,"  [val,max_ev]=max(abs(evecJ));\n");
  fprintf(OutFile,"  crevJ(%d,1)=Vars(max_ev,2);\n",count);
  fprintf(OutFile,"  for i=1:N,\n");
  fprintf(OutFile,"    if (Vars(i,2)==l), lnumJ=i; break; end\n");
  fprintf(OutFile,"  end\n");
  fprintf(OutFile,"  crsvJ(%d,2)=rankbus(svecJ,lnumJ);\n",count);
  fprintf(OutFile,"  crevJ(%d,2)=rankbus(evecJ,lnumJ);\n",count);
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Compute J_PF matrix:\n","%%");
  fprintf(OutFile,"P0=zeros(N);\n");
  fprintf(OutFile,"j=0;\n");
  fprintf(OutFile,"for i=1:N,\n");
  fprintf(OutFile,"  if(Vars(i,3)>0),\n");
  fprintf(OutFile,"    j=j+1;\n");
  fprintf(OutFile,"    Buses(j)=Vars(i,2);\n");
  fprintf(OutFile,"    if (Vars(i,2)==l), lnumPF=j; end\n");
  fprintf(OutFile,"    P0(Vars(i,1),j)=1;\n");
  fprintf(OutFile,"  end\n");
  fprintf(OutFile,"end\n");
  fprintf(OutFile,"M=j;\n");
  fprintf(OutFile,"for i=1:N,\n");
  fprintf(OutFile,"  if(Vars(i,3)==0), j=j+1; P0(Vars(i,1),j)=1; end\n");
  fprintf(OutFile,"end\n");
  fprintf(OutFile,"P0=sparse(P0); Jp=P0'*J*P0;\n");
  fprintf(OutFile,"J_PF=Jp(1:M,1:M);\n");
  fprintf(OutFile,"%s Minimum singular value and real |e-value| for standard P.F. J_PF:\n","%%");
  fprintf(OutFile,"[val,svecPF]=inviter(J_PF'*J_PF); svPF(%d)=sqrt(val);\n",count);
  fprintf(OutFile,"[val,evecPF]=inviter(J_PF); evPF(%d)=abs(val);\n",count);
  fprintf(OutFile,"  %s Critical bus numbers and ranking of bus l:\n","%%");
  fprintf(OutFile,"  [val,max_sv]=max(abs(svecPF));\n");
  fprintf(OutFile,"  crsvPF(%d,1)=Buses(max_sv);\n",count);
  fprintf(OutFile,"  [val,max_ev]=max(abs(evecPF));\n");
  fprintf(OutFile,"  crevPF(%d,1)=Buses(max_ev);\n",count);
  fprintf(OutFile,"  crsvPF(%d,2)=rankbus(svecPF,lnumPF);\n",count);
  fprintf(OutFile,"  crevPF(%d,2)=rankbus(evecPF,lnumPF);\n",count);
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Compute J_QV matrix:\n","%%");
  fprintf(OutFile,"P1=zeros(N);\n");
  fprintf(OutFile,"Mp=%d; j=Mp;\n",m);
  fprintf(OutFile,"for i=1:Mp,\n");
  fprintf(OutFile,"  if (Vbuses(i,2)==l), lnumQV=i; end\n");
  fprintf(OutFile,"  P1(Vbuses(i,1),i)=1;\n");
  fprintf(OutFile,"end;\n");
  fprintf(OutFile,"for i=1:N,\n");
  fprintf(OutFile,"  if(Vars(i,4)>0), j=j+1; P1(Vars(i,1),j)=1; end\n");
  fprintf(OutFile,"end\n");
  fprintf(OutFile,"P1=sparse(P1); Jp=P1'*J*P1;\n");
  fprintf(OutFile,"A1=Jp(1:Mp,1:Mp);  B1=Jp(1:Mp,Mp+1:N); C1=Jp(Mp+1:N,1:Mp); D1=Jp(Mp+1:N,Mp+1:N);\n");
  fprintf(OutFile,"J_QV=A1-B1*(D1\\C1);\n");
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Minimum singular value and real |e-value| for J_QV:\n","%%");
  fprintf(OutFile,"[val,svecQV]=inviter(J_QV'*J_QV); svQV(%d)=sqrt(val);\n",count);
  fprintf(OutFile,"[val,evecQV]=inviter(J_QV,-0.001); evQV(%d)=abs(val);\n",count);
  fprintf(OutFile,"  %s Critical V bus numbers and ranking of bus l:\n","%%");
  fprintf(OutFile,"  [val,max_sv]=max(abs(svecQV));\n");
  fprintf(OutFile,"  crsvQV(%d,1)=Vbuses(max_sv,2);\n",count);
  fprintf(OutFile,"  [val,max_ev]=max(abs(evecQV));\n");
  fprintf(OutFile,"  crevQV(%d,1)=Vbuses(max_ev,2);\n",count);
  fprintf(OutFile,"  crsvQV(%d,2)=rankbus(svecQV,lnumQV);\n",count);
  fprintf(OutFile,"  crevQV(%d,2)=rankbus(evecQV,lnumQV);\n",count);
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Reduced determinant detD_ll, for bus l=%d, '%12s':\n","%%",TFbus,TFname);
  fprintf(OutFile,"P2=speye(N);\n");
  fprintf(OutFile,"v=P2(1,:); P2(1,:)=P2(%d,:); P2(%d,:)=v; \n",TFnum-1,TFnum-1);
  fprintf(OutFile,"v=P2(2,:); P2(2,:)=P2(%d,:); P2(%d,:)=v; \n",TFnum,TFnum);
  fprintf(OutFile,"Jp=P2'*J*P2;\n");
  fprintf(OutFile,"A2=Jp(1:2,1:2);  B2=Jp(1:2,3:N); C2=Jp(3:N,1:2); D2=Jp(3:N,3:N);\n");
  fprintf(OutFile,"detD_ll(%d)=det(A2-B2*(D2\\C2));\n",count);
  fprintf(OutFile,"\n");
  fprintf(OutFile,"%s Test Function t_ll, for bus l=%d, '%12s':\n","%%",TFbus,TFname);
  fprintf(OutFile,"el=zeros(N,1); el(%d)=1; el=sparse(el);\n",TFnum);
  fprintf(OutFile,"J_ll=(speye(N)-el*el')*J+el*el';\n");
  fprintf(OutFile,"t_ll(%d)=el'*J*(J_ll\\el);\n",count);
  fclose(OutFile);
}
