/* Homotopy Continuation Method:  Group 2. */

#include "homot.h"

#ifdef WINDOWS
#include "Win UWPflow.h"
#include "GraphDLG.h"
#endif


/* ------- Global Variables ------ */
extern VALUETYPE *Dx,Dparam,param0,*x0,*x0p,Kh,Htol,SD0,AngTr,
                 DxiMax,VSFone,VSFinf,SF,TVI,lambda_o,
                 TotalPl,TotalQl,TotalPg,TotalQg;
extern INDEX TVIbus;
extern AClist *Vlist,*Vlistp;
extern int field;
extern BOOLEAN flagPrintTotalPl,flagPrintTotalQl,flagPrintTotalPg,flagPrintTotalQg;

#ifdef WINDOWS
extern CString dir;
extern GraphDLG* GraphDlg;
#endif

/* --------------------------- InList --------------------------------- */
#ifdef ANSIPROTO
BOOLEAN InList(ACbusData *ACptr,AClist *Vptr)
#else
BOOLEAN InList(ACptr,Vptr)
ACbusData *ACptr;
AClist *Vptr;
#endif
/* Check whether bus in list for V profiles */
{
  AClist *Lptr;

  for(Lptr=Vptr;Lptr!=NULL;Lptr=Lptr->Next) {if(ACptr==Lptr->AC) return(TRUE);}
  return(FALSE);

}


/* --------------------------- MakeVlist --------------------------------- */
#ifdef ANSIPROTO
void MakeVlist(FILE *Out)
#else
void MakeVlist(Out)
FILE *Out;
#endif
/* Prepare list of buses/areas for voltage profiles. */
{
  ACbusData *ACptr,*ACptrM;
  AreaData *Aptr;
  AClist *Lptr,*Lptrp;
  DCbusData *DCptr;
  char Line[BUFLEN],BusName[BUFLEN],Type[BUFLEN],*ptr,*Name;
  FILE *Input;
  INDEX N,i,count,countp;
  VALUETYPE MaxV=0;

  Vlist=Vlistp=NULL;
  ACptrM=NULL;
  Name=NameParameter('i');
  flagPrintTotalPl=flagPrintTotalQl=flagPrintTotalPg=flagPrintTotalQg=FALSE;

  if (!NullName(Name) && (Input=OpenInput(Name))!=NULL) {
    for (;;) {
      strcpy(Type,"");
      if (fgets(Line,BUFLEN,Input)==NULL) break;
      if ((count=sscanf(Line,"%d %s %s",&N,BusName,Type))>3 && strncmp(Line,"C",1)) {
        fprintf(stderr,"***Warning: Line-> %s",Line);
        fprintf(stderr,"            will be ignored in file %s.\n",Name);
      }
      else if (count>=2) {
        if (BusName[0]=='\"'||BusName[0]=='\'') {
          count=2;
          for(ptr=Line;*ptr!='\"'&&*ptr!='\'';ptr++); ptr++;
          for(i=0;*ptr!='\"'&&*ptr!='\''&&*ptr!='\n';BusName[i]= *ptr,i++,ptr++);
          BusName[i]='\0';
            ptr++;
          count+=sscanf(ptr,"%s",Type);
        }
        if (N==0 && !strcmp(BusName,"0") && !strcmp(Type,"PL"))      flagPrintTotalPl=TRUE;
        else if (N==0 && !strcmp(BusName,"0") && !strcmp(Type,"QL")) flagPrintTotalQl=TRUE;
        else if (N==0 && !strcmp(BusName,"0") && !strcmp(Type,"PG")) flagPrintTotalPg=TRUE;
        else if (N==0 && !strcmp(BusName,"0") && !strcmp(Type,"QG")) flagPrintTotalQg=TRUE;
        else {
          if (!strcmp(Type,"PA")){
            ACptr=NULL;
            for (Aptr=dataPtr->Area;Aptr!=NULL;Aptr=Aptr->Next)
              if(N==Aptr->N ||!strncmp(Aptr->Name,BusName,strlen(BusName))) break;
          }
          else {
            Aptr=NULL;
            for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next)
              if(N==ACptr->Num||!strncmp(ACptr->Name,BusName,strlen(BusName))) break;
          }
          if(ACptr!=NULL || Aptr!=NULL){
            Lptr=Vlist;
#ifdef WINDOWS
            Vlist= new AClist;
#else
            Vlist=(AClist *) malloc(sizeof(AClist));
            if (Vlist==NULL) {
              fclose(Out);
              fclose(Input);
              ErrorHalt("Insufficient memory to allocate profile List.");
              exit(ERROREXIT);
            }
#endif
            Vlist->AC=ACptr;
            Vlist->Area=Aptr;
            if (ACptr!=NULL) Vlist->N=ACptr->Num;
            else             Vlist->N=Aptr->N;
            if (!strcmp(Type,"V") || !strcmp(Type,"D") || !strcmp(Type,"PL") || !strcmp(Type,"QL") ||
                !strcmp(Type,"PG") ||!strcmp(Type,"QG") ||!strcmp(Type,"PA")) strcpy(Vlist->Type,Type);
            else strcpy(Vlist->Type,"V");
            Vlist->Next=Lptr;
            Vlist->Prev=NULL;
            if(Lptr!=NULL) Lptr->Prev=Vlist;
          }
          else if (strncmp(Line,"C",1)) {
            fprintf(stderr,"***Warning: Line-> %s",Line);
            fprintf(stderr,"            will be ignored in file %s.\n",Name);
          }
        }
      }
    }
    Lptr=Vlist;
    while(Lptr!=NULL){
      Lptrp=Lptr->Next;
      Lptr->Next=Lptr->Prev;
      Lptr->Prev=Lptrp;
      if(Lptrp==NULL) Vlist=Lptr;
      Lptr=Lptrp;
    }
    fclose(Input);
  }
  if (Vlistp==NULL) {
    countp=1;
    while (countp<=8 && countp<=Nac) {
      MaxV=0;
      for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
        if(Nac<=8 && !InList(ACptr,Vlistp)) {
          ACptrM=ACptr;
          break;
        }
        else if((ACptr->Cont!=NULL &&(QRcont || !strpbrk(ACptr->Type,"G"))) ||
                (!Rcont && strpbrk(ACptr->Type,"T"))||
                (!QRcont && strpbrk(ACptr->Type,"C"))) {
          i=ACvar[ACptr->N]+1;
          if (fabs(Dx[i])>MaxV && !InList(ACptr,Vlistp)) {MaxV=fabs(Dx[i]); ACptrM=ACptr;}
        }
      }
      if(ACptrM!=NULL){
        Lptr=Vlistp;
#ifdef WINDOWS
        Vlistp= new AClist;
#else
        Vlistp=(AClist *) malloc(sizeof(AClist));
        if (Vlistp==NULL) {
          fclose(Out);
          ErrorHalt("Insufficient memory to allocate Voltage Profile Bus List.");
          exit(ERROREXIT);
        }
#endif
        Vlistp->AC=ACptrM;
        strcpy(Vlistp->Type,"V");
        Vlistp->N=ACptrM->Num;
        Vlistp->Area=NULL;
        Vlistp->Next=Lptr;
        Vlistp->Prev=NULL;
        if(Lptr!=NULL) Lptr->Prev=Vlistp;
        ACptrM=NULL;
        countp++;
      } else break;
    }
  }
  if (Vlist==NULL) Vlist=Vlistp;
  if (ExistParameter('m')){
    if (ExistParameter('O')) {
      for (DCptr=dataPtr->DCbus;DCptr!=NULL;DCptr=DCptr->Next) {
        for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) if(Lptr->AC==DCptr->AC && !strcmp(Lptr->Type,"V")) break;
        if (Lptr==NULL) {
          Lptr=Vlist;
#ifdef WINDOWS
          Vlist= new AClist;
#else
          Vlist=(AClist *) malloc(sizeof(AClist));
          if (Vlist==NULL) {
            fclose(Out);
            ErrorHalt("Insufficient memory to allocate Voltage Profile Bus List.");
            exit(ERROREXIT);
          }
#endif
          Vlist->AC=DCptr->AC;
          Vlist->Area=NULL;
          Vlist->N=DCptr->AC->Num;
          strcpy(Vlist->Type,"V");
          Vlist->Next=Lptr;
          Vlist->Prev=NULL;
          if(Lptr!=NULL) Lptr->Prev=Vlist;
        }
        for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) if(Lptr->AC==DCptr->AC && !strcmp(Lptr->Type,"D")) break;
        if (Lptr==NULL) {
          Lptr=Vlist;
#ifdef WINDOWS
          Vlist= new AClist;
#else
          Vlist=(AClist *) malloc(sizeof(AClist));
          if (Vlist==NULL) {
            fclose(Out);
            ErrorHalt("Insufficient memory to allocate Voltage Profile Bus List.");
            exit(ERROREXIT);
          }
#endif
          Vlist->AC=DCptr->AC;
          Vlist->Area=NULL;
          Vlist->N=DCptr->AC->Num;
          strcpy(Vlist->Type,"D");
          Vlist->Next=Lptr;
          Vlist->Prev=NULL;
          if(Lptr!=NULL) Lptr->Prev=Vlist;
        }
      }
    }
	//extra % because when it's read by fprintf, one % will be ignored
    fprintf(Out,"%s ", "%%");
  }
  fprintf(Out,"L.F.    ");

#ifdef WINDOWS
  GraphDlg->totalElements = 0;
  GraphDlg->show = true;
  char tmpCaption[20];
  char tmpN[20];
  //clear last graph
  GraphDlg->Reset();
#endif

  for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) 
  {
	  fprintf(Out,"%s%-5d    ",Lptr->Type,Lptr->N);

#ifdef WINDOWS
		//initialize elements of the graph here where headings are printed
		//so this is the first point where the program can know how many elements there are
	  
		//create element
		if (GraphDlg->totalElements!=0)
			GraphDlg->m_GraphCtrl.AddElement();

		long color = RGB(rand(), rand(), rand());
		GraphDlg->m_GraphCtrl.SetElementLineColor(color);

		//set annotation
		strcpy(tmpCaption,Lptr->Type);
		sprintf(tmpN, "%-5d", Lptr->N);
		strcat(tmpCaption, tmpN);
		GraphDlg->m_GraphCtrl.AddAnnotation();
		GraphDlg->m_GraphCtrl.SetAnnoLabelCaption(tmpCaption);
		GraphDlg->m_GraphCtrl.SetAnnoLabelColor(color);
		
		//increase element count
		GraphDlg->totalElements++;
	  
#endif
  }
  if (ExistParameter('e')) {
    for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) if (Lptr->AC!=NULL) {
      if (Lptr->AC->Gen!=NULL) {
        fprintf(Out,"Ia%-5d    ",Lptr->N);
        fprintf(Out,"Eq%-5d    ",Lptr->N);
        fprintf(Out,"dg%-5d    ",Lptr->N);
      }
    }
  }
  if (ExistParameter('O')) {
    fprintf(Out,"Vac    ");
    if (field>7) fprintf(Out,"    ");
  }
  if (ExistParameter('O') || ExistParameter('e')) {
    for (i=0,DCptr=dataPtr->DCbus;DCptr!=NULL;DCptr=DCptr->Next) if(!strcmp(DCptr->Type,"R")) {
      i+=2;
      if (ExistParameter('O')) {
        fprintf(Out,"a%1d    ",i-1); if (field>7) fprintf(Out,"    ");
        fprintf(Out,"b%1d    ",i-1); if (field>7) fprintf(Out,"    ");
        fprintf(Out,"d%1d    ",i-1); if (field>7) fprintf(Out,"    ");
        fprintf(Out,"a%1d    ",i);   if (field>7) fprintf(Out,"    ");
        fprintf(Out,"b%1d    ",i);   if (field>7) fprintf(Out,"    ");
        fprintf(Out,"d%1d    ",i);   if (field>7) fprintf(Out,"    ");
      } else if (ExistParameter('e')) {
        fprintf(Out,"alR_%1d    ",i);
        fprintf(Out,"gaI_%1d    ",i);
        fprintf(Out,"Id_%1d    ",i);
      }
    }
  }
  if (flagPrintTotalPl) fprintf(Out,"PL        ");
  if (flagPrintTotalQl) fprintf(Out,"QL        ");
  if (flagPrintTotalPg) fprintf(Out,"PG        ");
  if (flagPrintTotalQg) fprintf(Out,"QG        ");
  if (ExistParameter('f')) {
    fprintf(Out,"VSFone    ");
    fprintf(Out,"VSFbus    ");
    fprintf(Out,"VSFinf    ");
    fprintf(Out,"SF    ");
    if (TVI!=0) fprintf(Out,"TVI_%d    %d_Rank",TVIbus,TVIbus);

  }
  if (ExistParameter('m')) fprintf(Out,"\nx=[\n");
  else fprintf(Out,"\n");
  if (ExistParameter('d') && (Out!=stdout || !NullName(NameParameter('l')))) {
    fprintf(stderr,"L.F.    ");
    for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) fprintf(stderr,"%s%-5d    ",Lptr->Type,Lptr->N);
    for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) if (Lptr->AC!=NULL) {
      if (Lptr->AC->Gen!=NULL) {
        fprintf(stderr,"Ia%-5d    ",Lptr->N);
        fprintf(stderr,"Eq%-5d    ",Lptr->N);
        fprintf(stderr,"dg%-5d    ",Lptr->N);
      }
    }
    for (i=0,DCptr=dataPtr->DCbus;DCptr!=NULL;DCptr=DCptr->Next) if(!strcmp(DCptr->Type,"R")) {
      i++;
      fprintf(stderr,"alR_%-2d    ",i);
      fprintf(stderr,"muR_%-2d    ",i);
      fprintf(stderr,"tpR_%-2d    ",i);
      fprintf(stderr,"gaI_%-2d    ",i);
      fprintf(stderr,"muI_%-2d    ",i);
      fprintf(stderr,"tpI_%-2d    ",i);
    }
    fprintf(stderr,"\n");
  }
}


/* --------------------------- VoltProf --------------------------------- */
#ifdef ANSIPROTO
void VoltProf(BOOLEAN flag,FILE *Out)
#else
void VoltProf(flag,Out)
BOOLEAN flag;
FILE *Out;
#endif
/* Write voltage profiles. */
{
  ACbusData *ACptr;
  AreaData *Aptr;
  AClist *Lptr;
  DClist *DCLptr;
  DCbusData *DCptrI,*DCptrR;
  ElementData *Eptr;
  ElementList *ELptr;
  VALUETYPE Vi,Vj,di,dj,Gi,Gj,G,B,Gp,Bp,P;

#ifdef WINDOWS
  int currElement=0;
#endif

  Print(Out,0,6,4,lambda_o+lambda); fprintf(Out,"    ");
  for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) {
    ACptr=Lptr->AC;
    Aptr=Lptr->Area;
    if (ACptr!=NULL) {
      if (!strcmp(Lptr->Type,"V")) 
	  {
		  Print(Out,0,6,4,ACptr->V); 
		  fprintf(Out,"    "); 
#ifdef WINDOWS
		  if (GraphDlg!=NULL)
		  {
			  GraphDlg->m_GraphCtrl.PlotXY(lambda_o+lambda, ACptr->V, currElement);
			  GraphDlg->CheckRange(lambda_o+lambda, ACptr->V);
			  if (++currElement > 7) 
				  currElement = 0;

			  GraphDlg->m_GraphCtrl.SetYLabel("V");
		  }
#endif
	  }
      else if (!strcmp(Lptr->Type,"D")) 
	  {
		  Print(Out,0,6,2,ACptr->Ang*AngTr); 
		  fprintf(Out,"    ");
#ifdef WINDOWS
		  if (GraphDlg!=NULL)
		  {
			  GraphDlg->m_GraphCtrl.PlotXY(lambda_o+lambda, ACptr->Ang*AngTr, currElement);
			  GraphDlg->CheckRange(lambda_o+lambda, ACptr->Ang*AngTr);
			  if (++currElement > 7) 
				  currElement = 0;

			  GraphDlg->m_GraphCtrl.SetYLabel("D");
		  }
#endif
	  }
      else if (!strcmp(Lptr->Type,"PG")) 
	  {
		  Print(Out,0,6,2,ACptr->PG*Sn); 
		  fprintf(Out,"    ");
#ifdef WINDOWS
		  if (GraphDlg!=NULL)
		  {
			  GraphDlg->m_GraphCtrl.PlotXY(lambda_o+lambda, ACptr->PG*Sn, currElement);
			  GraphDlg->CheckRange(lambda_o+lambda, ACptr->PG*Sn);
			  if (++currElement > 7) 
				  currElement = 0;

			  GraphDlg->m_GraphCtrl.SetYLabel("PG");
		  }
#endif
	  }
      else if (!strcmp(Lptr->Type,"QG")) 
	  {
		  Print(Out,0,6,2,ACptr->Qg*Sn); 
		  fprintf(Out,"    ");
#ifdef WINDOWS
		  if (GraphDlg!=NULL)
		  {
			  GraphDlg->m_GraphCtrl.PlotXY(lambda_o+lambda, ACptr->Qg*Sn, currElement);
			  GraphDlg->CheckRange(lambda_o+lambda, ACptr->Qg*Sn);
			  if (++currElement > 7) 
				  currElement = 0;

			  GraphDlg->m_GraphCtrl.SetYLabel("QG");
		  }
#endif	  
	  }
      else if (!strcmp(Lptr->Type,"PL")) 
	  {
		  Print(Out,0,6,2,ACptr->PL*Sn); 
		  fprintf(Out,"    ");
#ifdef WINDOWS
		  if (GraphDlg!=NULL)
		  {
			  GraphDlg->m_GraphCtrl.PlotXY(lambda_o+lambda, ACptr->PL*Sn, currElement);
			  GraphDlg->CheckRange(lambda_o+lambda, ACptr->PL*Sn);
			  if (++currElement > 7) 
				  currElement = 0;

			  GraphDlg->m_GraphCtrl.SetYLabel("PL");
		  }
#endif	  
	  }
      else if (!strcmp(Lptr->Type,"QL")) 
	  {
		  Print(Out,0,6,2,ACptr->QL*Sn); 
		  fprintf(Out,"    ");
#ifdef WINDOWS
		  if (GraphDlg!=NULL)
		  {
			  GraphDlg->m_GraphCtrl.PlotXY(lambda_o+lambda, ACptr->QL*Sn, currElement);
			  GraphDlg->CheckRange(lambda_o+lambda, ACptr->QL*Sn);
			  if (++currElement > 7) 
				  currElement = 0;

			  GraphDlg->m_GraphCtrl.SetYLabel("QL");
		  }
#endif	  
	  }
    }
    else if (Aptr!=NULL) {
      Aptr->SPg=0;
      for (DCLptr=Aptr->DC;DCLptr!=NULL;DCLptr=DCLptr->Next) {
        DCptrR=DCLptr->DC->Meter;
        P= -DCptrR->P;
        if(DCptrR->Area!=Aptr) P= -P;
        Aptr->SPg=Aptr->SPg+P;
      }
      for(ELptr=Aptr->Elem;ELptr!=NULL;ELptr=ELptr->Next) {
        Eptr=ELptr->Eptr;
        Vi=Eptr->From->V;  di=Eptr->From->Ang;
        Vj=Eptr->To->V;    dj=Eptr->To->Ang;
        G=(Eptr->G*cos(Eptr->Ang)-Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
        B=(Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
        Gi=(Eptr->G1+Eptr->G)*pow(Eptr->Tap,2.0)-G;
        Gp=(Eptr->G*cos(Eptr->Ang)+Eptr->B*sin(Eptr->Ang))*Eptr->Tap;
        Bp=(-Eptr->G*sin(Eptr->Ang)+Eptr->B*cos(Eptr->Ang))*Eptr->Tap;
        Gj=Eptr->G+Eptr->G2-Gp;
        if (Eptr->From==Eptr->Meter) {
          P=Vi*Vi*(Gi+G)-Vi*Vj*(G*cos(di-dj)+B*sin(di-dj));
        } else {
          P=Vj*Vj*(Gj+Gp)-Vi*Vj*(Gp*cos(dj-di)+Bp*sin(dj-di));
        }
        if(Eptr->Meter->Area!=Aptr) P= -P;
        Aptr->SPg=Aptr->SPg+P;
      }
      Print(Out,0,6,2,Aptr->SPg*Sn); fprintf(Out,"    ");
    }
  }
  if (ExistParameter('e')) for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) {
    ACptr=Lptr->AC;
    if (ACptr!=NULL && ACptr->Gen!=NULL) {
      Print(Out,0,6,4,ACptr->Gen->Ia); fprintf(Out,"    ");
      Print(Out,0,6,4,ACptr->Gen->Eq); fprintf(Out,"    ");
      Print(Out,0,6,2,ACptr->Gen->dg*AngTr); fprintf(Out,"    ");
   }
  }
  if (ExistParameter('O')) {
    TEFac(flag);
    Print(Out,0,field,4,Vac);  fprintf(Out,"    ");
    TEFdc(Out);
  }
  else if (ExistParameter('e')) {
    for (DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next) if(!strcmp(DCptrR->Type,"R")) {
      DCptrI=DCptrR->To;
      Print(Out,0,6,2,DCptrR->Alfa*AngTr); fprintf(Out,"    ");
      Print(Out,0,6,2,DCptrI->Gamma*AngTr); fprintf(Out,"    ");
      Print(Out,0,6,2,DCptrI->Id*1000.*Sn/DCptrR->Vn); fprintf(Out,"    ");
    }
  }
  if (flagPrintTotalPl) { Print(Out,0,8,2,TotalPl*Sn); fprintf(Out,"    "); }
  if (flagPrintTotalQl) { Print(Out,0,8,2,TotalQl*Sn); fprintf(Out,"    "); }
  if (flagPrintTotalPg) { Print(Out,0,8,2,TotalPg*Sn); fprintf(Out,"    "); }
  if (flagPrintTotalQg) { Print(Out,0,8,2,TotalQg*Sn); fprintf(Out,"    "); }
  if (ExistParameter('d') && (Out!=stdout || !NullName(NameParameter('l')))) {
    Print(stderr,0,6,4,lambda_o+lambda); fprintf(stderr,"    ");
    for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) {
      ACptr=Lptr->AC;
      Aptr=Lptr->Area;
      if (ACptr!=NULL) {
        if (!strcmp(Lptr->Type,"V")) {Print(stderr,0,6,4,ACptr->V); fprintf(stderr,"    "); }
        else if (!strcmp(Lptr->Type,"D")) {Print(stderr,0,6,2,ACptr->Ang*AngTr); fprintf(stderr,"    ");}
        else if (!strcmp(Lptr->Type,"PG")) {Print(stderr,0,6,2,ACptr->PG*Sn); fprintf(stderr,"    ");}
        else if (!strcmp(Lptr->Type,"QG")) {Print(stderr,0,6,2,ACptr->Qg*Sn); fprintf(stderr,"    ");}
        else if (!strcmp(Lptr->Type,"PL")) {Print(stderr,0,6,2,ACptr->PL*Sn); fprintf(stderr,"    ");}
        else if (!strcmp(Lptr->Type,"QL")) {Print(stderr,0,6,2,ACptr->QL*Sn); fprintf(stderr,"    ");}
      }
      else if (Aptr!=NULL) {Print(stderr,0,6,2,Aptr->SPg*Sn); fprintf(stderr,"    ");}
    }
    for(Lptr=Vlist;Lptr!=NULL;Lptr=Lptr->Next) {
      ACptr=Lptr->AC;
      if (ACptr!=NULL && ACptr->Gen!=NULL) {
        Print(stderr,0,6,4,ACptr->Gen->Ia); fprintf(stderr,"    ");
        Print(stderr,0,6,4,ACptr->Gen->Eq); fprintf(stderr,"    ");
        Print(stderr,0,6,2,ACptr->Gen->dg*AngTr); fprintf(stderr,"    ");
     }
    }
    for (DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next) if(!strcmp(DCptrR->Type,"R")) {
      DCptrI=DCptrR->To;
      Print(stderr,0,6,2,DCptrR->Alfa*AngTr); fprintf(stderr,"    ");
      Print(stderr,0,6,2,(PI-DCptrR->Alfa-DCptrR->Gamma)*AngTr); fprintf(stderr,"    ");
      Print(stderr,0,6,4,DCptrR->Tap); fprintf(stderr,"    ");
      Print(stderr,0,6,2,DCptrI->Gamma*AngTr); fprintf(stderr,"    ");
      Print(stderr,0,6,2,(PI-DCptrI->Alfa-DCptrI->Gamma)*AngTr); fprintf(stderr,"    ");
      Print(stderr,0,6,4,DCptrI->Tap); fprintf(stderr,"    ");
    }
    fprintf(stderr,"\n");
  }

}


/* --------------------------- PrintDirection --------------------------------- */
#ifdef ANSIPROTO
void PrintDirection(char Option,VALUETYPE *vector,VALUETYPE Max)
#else
void PrintDirection(Option,vector,Max)
char Option;
VALUETYPE *vector,Max;
#endif
/* Print direction vector  */
{
  INDEX i,j,k,N,I,J;
  ACbusData *ACptr;
  DCbusData *DCptr,*DCptrR,*DCptrI;
  SVCbusData *SVCptr;                  /* FACTS */
  TCSCbusData *TCSCptr;                /* FACTS */
  STATCOMbusData *STATCOMptr;          /* FACTS */
  ElementData *Eptr;
  ElementList *ELptr;
  char str[80],type[2];
  FILE *Out;

  Out=OpenOutput(NameParameter(Option));
  N=Jac->n1;
  fprintf(Out,"%d 1\n",N);
  for (i=0,ACptr=dataPtr->ACbus; ACptr!=NULL; ACptr=ACptr->Next){
    if (ACptr->Cont!=NULL){
      if (strpbrk(ACptr->Type,"S")) sprintf(str,"kg%-d",ACptr->Num);
      else                          sprintf(str,"d%-d",ACptr->Num);
      i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"V%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
    else if(strpbrk(ACptr->Type,"L")){
      sprintf(str,"d%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,"l",vector[i]/Max);
    }
    else if(QRcont && strpbrk(ACptr->Type,"C")){
      sprintf(str,"d%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Q%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
    else if(Rcont && strpbrk(ACptr->Type,"T")){
      sprintf(str,"d%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      for(ELptr=ACptr->Reg;ELptr!=NULL;ELptr=ELptr->Next){
         Eptr=ELptr->Eptr;
         I=Eptr->From->Num;
         J=Eptr->To->Num;
         if(!strcmp(Eptr->Type,"R")) break;
      }
      sprintf(str,"1/t%-d_%-d",I,J); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
    else if(strpbrk(ACptr->Type,"Q") || strpbrk(ACptr->Type,"V") || (!QRcont && strpbrk(ACptr->Type,"G"))){
      if (strpbrk(ACptr->Type,"S")) sprintf(str,"kg%-d",ACptr->Num);
      else                          sprintf(str,"d%-d",ACptr->Num);
      i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Qg%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
    else if(strpbrk(ACptr->Type,"Z")) {
      sprintf(str,"d%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Qz%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
    else if(strpbrk(ACptr->Type,"S")){
      sprintf(str,"kg%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Qg%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
    if(Acont && strpbrk(ACptr->Type,"A")){
      sprintf(str,"kg%-d",ACptr->Num); i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
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
           sprintf(str,"a%-d_%-d",I,J); i++;
           fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
         } else if(strpbrk(Eptr->Type,"PM")){
           sprintf(str,"P%-d_%-d",I,J); i++;
           fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
         } else if(!strcmp(Eptr->Type,"RQ")){
           sprintf(str,"1/t%-d_%-d",I,J); i++;
           fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
         } else {
           sprintf(str,"Q%-d_%-d",I,J); i++;
           fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
         }
      }
    }
    if (ACptr->Gen!=NULL) {
      i=ACptr->Gen->Nvar;
      if (!strpbrk(ACptr->cont,"E")) sprintf(str,"Eq%-d",ACptr->Num);
      else                           sprintf(str,"Qg%-d",ACptr->Num);
      i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"dg%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Vr%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Vi%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Ir%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Ii%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Vq%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Vd%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Iq%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      sprintf(str,"Id%-d",ACptr->Num); i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      if (!strpbrk(ACptr->cont,"I")) sprintf(str,"Ia%-d",ACptr->Num);
      else                           sprintf(str,"Qg%-d",ACptr->Num);
      i++; fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
  }
  for(k=0,DCptrR=dataPtr->DCbus;DCptrR!=NULL;DCptrR=DCptrR->Next){
    DCptrI=DCptrR->To;
    if(!strcmp(DCptrR->Type,"R")){
      k++;
      for(j=1;j<=2;j++){
        if(j==1) { DCptr=DCptrR; strcpy(type,"r"); }
        else { DCptr=DCptrI; strcpy(type,"i"); }
        if(strcmp(DCptr->Cont1,"VD")&&strcmp(DCptr->Cont2,"VD")) {
          sprintf(str,"Vd%1s%-d",type,k); i++;
          fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
        }
        if(strcmp(DCptr->Cont1,"AT")&&strcmp(DCptr->Cont2,"AT")) {
          sprintf(str,"t%1s%-d",type,k); i++;
          fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
        }
        if(strcmp(DCptr->Cont1,"AL")&&strcmp(DCptr->Cont2,"AL")) {
          sprintf(str,"al%1s%-d",type,k); i++;
          fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
        }
        if(strcmp(DCptr->Cont1,"GA")&&strcmp(DCptr->Cont2,"GA")) {
          sprintf(str,"ga%1s%-d",type,k); i++;
          fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
        }
        sprintf(str,"S%1s%-d",type,k); i++;
        fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
        if(strcmp(DCptr->Cont1,"PA")&&strcmp(DCptr->Cont2,"PA")) {
          sprintf(str,"P%1s%-d",type,k); i++;
          fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
        }
        if(strcmp(DCptr->Cont1,"QA")&&strcmp(DCptr->Cont2,"QA")) {
          sprintf(str,"Q%1s%-d",type,k); i++;
          fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
        }
      }
      if(strcmp(DCptrR->Cont1,"ID")&&strcmp(DCptrR->Cont2,"ID")&&
         strcmp(DCptrI->Cont1,"ID")&&strcmp(DCptrI->Cont2,"ID")) {
        sprintf(str,"Id%-d",k); i++;
        fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
      }
    }
  }
                             /* FACTS */
  for(k=0,SVCptr=dataPtr->SVCbus;SVCptr!=NULL;SVCptr=SVCptr->Next){
    k++;
    sprintf(str,"Qsvc%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"Bv%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    if(!strcmp(SVCptr->Cont,"AL")){
      sprintf(str,"alpha%-d",k);i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
    else {
      sprintf(str,"Vrefc%-d",k);i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
  }

  for(k=0,TCSCptr=dataPtr->TCSCbus;TCSCptr!=NULL;TCSCptr=TCSCptr->Next){
    k++;
    sprintf(str,"Ptcsc%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"Qtcsck%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"Qtcscm%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"Be%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"alpha%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"Itcsc%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"delta%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
  }

  for(k=0,STATCOMptr=dataPtr->STATCOMbus;STATCOMptr!=NULL;STATCOMptr=STATCOMptr->Next){
    k++;
    if(!strcmp(STATCOMptr->Cont,"PW") || !strcmp(STATCOMptr->Cont,"AL")){
      sprintf(str,"Istat%-d",k);i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    } else {
      sprintf(str,"Vrefc%-d",k);i++;
      fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    }
    sprintf(str,"theta%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"Vdc%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"k%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"alpha%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"Pstat%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
    sprintf(str,"Qstat%-d",k);i++;
    fprintf(Out,"%4d %8s %-11.5g\n",i,str,vector[i]/Max);
  }
                             /* END FACTS */

  if (Option!='Y') {
    if (!Bl) fprintf(Out,"%4d %8s %-11.5g\n",++i,"l",Dparam);
    else {
      sprintf(str,"V%-d",BlPtr->Num);
      fprintf(Out,"%4d %8s %-11.5g\n",++i,str,Dparam);
    }
  }
  fprintf(Out,"0 0 0.0\n");
  fclose(Out);
}


/* --------------------------- IndicesMatlab --------------------------------- */
#ifdef ANSIPROTO
void IndicesMatlab(INDEX count)
#else
void IndicesMatlab(count)
INDEX count;
#endif
/* Print plotting and other Matlab commands needed by the -0 option (VS indices) */
{
  char Namebase[80],Name[80];
  FILE *OutFile;
  INDEX i;
  

  strcpy(Namebase,NameParameter('0'));
  if(NullName(Namebase)) return;
  sprintf(Name,"%s.m",Namebase);
  OutFile=OpenOutput(Name);
  fprintf(OutFile,"clear lambda evJ svJ evPF svPF evQV svQV detD_ll t_ll\n");
  fprintf(OutFile,"clear crsvJ crevJ crsvPF crevPF crsvQV crevQV\n");
  fprintf(OutFile,"warning off\n");
  for(i=1; i<=count; i++) fprintf(OutFile,"%s%d\n",Namebase,i);
  fprintf(OutFile,"figure; plot(lambda,evJ); \n");
  fprintf(OutFile,"title('Full matrix |e-val.| index');\n");
  fprintf(OutFile,"xlabel('lambda [p.u.]');\n");
  fprintf(OutFile,"figure; plot(lambda,svJ); \n");
  fprintf(OutFile,"title('Full matrix sing. val. index');\n");
  fprintf(OutFile,"xlabel('lambda [p.u.]');\n");
  fprintf(OutFile,"disp(' '); disp('Critical bus numbers and bus l rank for J indices are stored in crevJ and crsvJ')\n");
  fprintf(OutFile,"figure; plot(lambda,evPF); \n");
  fprintf(OutFile,"title('Power flow matrix |e-val.| index');\n");
  fprintf(OutFile,"xlabel('lambda [p.u.]');\n");
  fprintf(OutFile,"figure; plot(lambda,svPF); \n");
  fprintf(OutFile,"title('Power flow matrix sing. val. index');\n");
  fprintf(OutFile,"xlabel('lambda [p.u.]');\n");
  fprintf(OutFile,"disp(' '); disp('Critical bus numbers and bus l rank for J_PV indices are stored in crevPF and crsvPF')\n");
  fprintf(OutFile,"figure; plot(lambda,evQV,lambda,svQV,'-.'); \n");
  fprintf(OutFile,"title('J_{QV} matrix |e-val.| and sing. val. indices');\n");
  fprintf(OutFile,"xlabel('lambda [p.u.]');\n");
  fprintf(OutFile,"legend('e-v','s.v.');\n");
  fprintf(OutFile,"disp(' '); disp('Critical bus numbers and bus l rank the J_QV indices are stored in crevQV and crsvQV')\n");
  fprintf(OutFile,"figure; plot(lambda,detD_ll); \n");
  fprintf(OutFile,"title('Reduced det. index');\n");
  fprintf(OutFile,"xlabel('lambda [p.u.]');\n");
  fprintf(OutFile,"figure; plot(lambda,-t_ll); \n");
  fprintf(OutFile,"title('Test func. index');\n");
  fprintf(OutFile,"xlabel('lambda [p.u.]');\n");
  fprintf(OutFile,"disp(' '); disp('Bus l used for the red. det. and test func. indices is: %d %s')\n",TFbus,TFname);
  fprintf(OutFile,"warning on\n\n");
  fclose(OutFile);
  /* --------------- Create 'rankbus.m' file needed for Matlab computations --------------- */
  OutFile=OpenOutput("rankbus.m");
  fprintf(OutFile,"function num_rank=rankbus(vec,num)\n");
  fprintf(OutFile,"%s\n","%%");
  fprintf(OutFile,"%s Rank entry 'num' in abs(vec).\n","%%");
  fprintf(OutFile,"%s\n","%%");
  fprintf(OutFile,"\n");
  fprintf(OutFile,"N=length(vec);\n");
  fprintf(OutFile,"[val,I]=sort(abs(vec));\n");
  fprintf(OutFile,"I=flipud(I);\n");
  fprintf(OutFile,"for i=1:N,\n");
  fprintf(OutFile,"  if (I(i)==num), num_rank=i; break; end\n");
  fprintf(OutFile,"end\n");
  fclose(OutFile);
  /* ----------------  Create 'inviter.m' file needed for Matlab computations --------------- */
  OutFile=OpenOutput("inviter.m");
  fprintf(OutFile,"function [e_val,v]=inviter(A,e_o,tol,iter,warn)\n");
  fprintf(OutFile,"%s\n","%%");
  fprintf(OutFile,"%s Inverse iteration method to compute a real e-value.\n","%%");
  fprintf(OutFile,"%s Designed for sparse matrices, but works with full matrices.\n","%%");
  fprintf(OutFile,"%s\n","%%");
  fprintf(OutFile,"%s      [v,e_val]=inviter(A,e_o,tol,iter,warn)\n","%%");
  fprintf(OutFile,"%s \n","%%");
  fprintf(OutFile,"%s Input:  A    -> NxN matrix\n","%%");
  fprintf(OutFile,"%s         e_o  -> Optional eigenvalue guess (default 0)\n","%%");
  fprintf(OutFile,"%s          tol  -> Optional convergence tolerance (default 1e-4)\n","%%");
  fprintf(OutFile,"%s         iter -> Optional maximum number of iterations (default 30)\n","%%");
  fprintf(OutFile,"%s         warn -> Use 0 to cancel display of no convergence warning message\n","%%");
  fprintf(OutFile,"%s                 (default 1, i.e., display warning)\n","%%");
  fprintf(OutFile,"%s\n","%%");
  fprintf(OutFile,"%s OutFileput: e_val -> Eigenvalue\n","%%");
  fprintf(OutFile,"%s         v     -> Eigenvector\n","%%");
  fprintf(OutFile,"%s\n","%%");
  fprintf(OutFile,"%s\n","%%");
  fprintf(OutFile,"%s Copyright (c) Claudio Canizares, Shu Zhang, 1996, 2006.\n","%%");
  fprintf(OutFile,"%s University of Waterloo, Waterloo, Canada\n","%%");
  fprintf(OutFile,"%s\n","%%");
  fprintf(OutFile,"if nargin<2, e_o=0; sign=1;\n");
  fprintf(OutFile,"else\n");
  fprintf(OutFile,"  if e_o<0, sign=-1;  else, sign=1;  end\n");
  fprintf(OutFile,"end\n");
  fprintf(OutFile,"if nargin<3, tol=1e-4; else, tol=abs(tol); end\n");
  fprintf(OutFile,"if nargin<4, iter=30; else, iter=abs(iter); end\n");
  fprintf(OutFile,"if nargin<5, warn=1; end\n");
  fprintf(OutFile,"\n");
  fprintf(OutFile,"N=size(A);\n");
  fprintf(OutFile,"count=1; \n");
  fprintf(OutFile,"conv=0;\n");
  fprintf(OutFile,"J=sparse(A-e_o*eye(N(1)));\n");
  fprintf(OutFile,"[L,U,P]=lu(J);\n");
  fprintf(OutFile,"while (count<=2)\n");
  fprintf(OutFile,"  v=ones(N(1),1);\n");
  fprintf(OutFile,"  for I=1:iter,\n");
  fprintf(OutFile,"   e_val=sign*1/norm(v,1);\n");
  fprintf(OutFile,"   v=v*e_val;\n");
  fprintf(OutFile,"   z=L\\(P*v);\n");
  fprintf(OutFile,"   vp=U\\z;\n");
  fprintf(OutFile,"   if(norm((1/e_val)*v-vp,1)<tol) break; end,\n");
  fprintf(OutFile,"   v=vp;\n");
  fprintf(OutFile,"  end\n");
  fprintf(OutFile,"  if (I<iter)\n");
  fprintf(OutFile,"    conv=1;\n");
  fprintf(OutFile,"    break;\n");
  fprintf(OutFile,"  else\n");
  fprintf(OutFile,"    sign=-sign;\n");
  fprintf(OutFile,"  end\n");
  fprintf(OutFile,"  count=count+1;\n");
  fprintf(OutFile,"end\n");
  fprintf(OutFile,"e_val=e_o+sign*1/norm(vp,1);\n");
  fprintf(OutFile,"v=e_val*vp;\n");
  fprintf(OutFile,"\n");
  fprintf(OutFile,"if (conv==0 & warn==1)\n");
  fprintf(OutFile,"  disp(' ')\n");
  fprintf(OutFile,"  disp('Warning: Inverse iteration method failed to converge')\n");
  fprintf(OutFile,"  str=sprintf('         for tol=%s6.4e, iter=%sd.',tol,iter);\n","%%","%%");
  fprintf(OutFile,"  disp(str)\n");
  fprintf(OutFile,"  disp(' ')\n");
  fprintf(OutFile,"end\n");
  fclose(OutFile);
}


