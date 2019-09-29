#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

/*    This routine reads initial values for AC voltage phasors.
      The input format is in "i V Ang" format, i.e., bus number
      and Voltage magnitude and angle (terminate list with a 0). */

#include <limits.h>
#include "readdata.h"


/* --------------------- ReadInit ---------------------------------- */
#ifdef ANSIPROTO
BOOLEAN ReadInit(void)
#else
BOOLEAN ReadInit()
#endif
 /* Main routine. */
{
  ACbusData *ACptr;
  DCbusData *DCptr;
  AreaData *Aptr;
  INDEX N;
  VALUETYPE V,d,SPg=0,DPg,Sum,Pn,Qn,Pz,Qz,val,PgMax,Vmax,Vmin;
  VALUETYPE Ra,Xd,Xq,IaMax,EqMax,EqMin,Smax;
  VALUETYPE Pg,Qg,Eq,dg,Vr,Vi,Ir,Ii,Vq,Vd,Iq,Id,Ia,apr,api;
  char *Name,BusName[BUFLEN],Line[BUFLEN],*ptr,Vars[BUFLEN],variable[3];
  FILE *InputFile;
  int i,count,num;

  num=6;
 /* -------------------- Generation and Load Factors --------------- */
  Name=NameParameter('K');
  Sum=0;
  if (!NullName(Name) && (InputFile=OpenInput(Name))!=NULL) {
    for (;;) {
      DPg=Pn=Qn=Pz=Qz=PgMax=Vmax=Vmin=0;
      if (fgets(Line,BUFLEN,InputFile)==NULL) break;
      if (Line[0]!='C') {
        if (sscanf(Line,"%d %s",&N,BusName)!=2) count=0;
        else if (BusName[0]=='\"'||BusName[0]=='\'') {
           count=2;
           for(ptr=Line;*ptr!='\"'&&*ptr!='\'';ptr++); ptr++;
           for(i=0;*ptr!='\"'&&*ptr!='\''&&*ptr!='\n';BusName[i]= *ptr,i++,ptr++);
           BusName[i]='\0';
           ptr++;
           count+=sscanf(ptr,"%lf %lf %lf %lf %lf %lf %lf %lf",&DPg,&Pn,&Qn,&PgMax,&Smax,&Vmax,&Vmin,&Pz,&Qz);
        } else {
           count=sscanf(Line,"%d %s %lf %lf %lf %lf %lf %lf %lf %lf",&N,BusName,&DPg,&Pn,&Qn,&PgMax,&Smax,&Vmax,&Vmin,&Pz,&Qz);
        }
        if (!strcmp(BusName,"0")) BusName[0]='\n';
        if (count>=num){
           for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
              if(N==ACptr->Num||!strncmp(ACptr->Name,BusName,strlen(BusName))) break;
           }
           if(ACptr!=NULL){
              if (PgMax>ACptr->Pg && PgMax>0) ACptr->PgMax=PgMax;
              if (Smax>=ACptr->PgMax && Smax>0) ACptr->Smax=Smax;
              if (Vmax>0 && Vmin>0 && Vmax>Vmin) { ACptr->Vlmax=Vmax; ACptr->Vlmin=Vmin; }
              ACptr->DPG=ACptr->DPg=DPg;
              if (Narea>1 && Acont) ACptr->Area->Slack->Area->SPg+=DPg*DPg;
              else SPg+=DPg*DPg;
              ACptr->Pnl=Pn; ACptr->Qnl=Qn;
              ACptr->Pzl=Pz; ACptr->Qzl=Qz;
              Sum+=Pn+Qn+Pz+Qz;
           }
        } else ACptr=NULL;
        if (ACptr==NULL) {
           fCustomPrint(stderr,"***Warning: Line-> %s",Line);
           fCustomPrint(stderr,"            will be ignored in file %s.\n",Name);
        }
      }
    }
    if (Narea>1 && Acont)
      for(Aptr=dataPtr->Area;Aptr!=NULL;Aptr=Aptr->Next) {
         if(!Aptr->SPg) Aptr->Slack->DPg=Aptr->SPg=1;
      }
    else if (!ExistParameter('6') || NullName(NameParameter('6'))) {
      if (!SPg) {
         for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) if(strpbrk(ACptr->Type,"S")) break;
         ACptr->DPg=SPg=1;
      }
      for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next){
        if(Narea>1 && Acont) ACptr->DPg=ACptr->DPg/sqrt(ACptr->Area->Slack->Area->SPg);
        else ACptr->DPg=ACptr->DPg/sqrt(SPg);
      }
    }
    fclose(InputFile);
  }
  if (!flagKdirection) {
    if (ExistParameter('v') && !Sum) {
      fCustomPrint(stderr,"ERROR: The -v option will yield a singular Jacobian in voltage collapse\n");
      fCustomPrint(stderr,"       studies since Pnl, Qnl, Pzl, and Qzl are zero in all load buses.\n");
      stopExecute(ERROREXIT);
    } else if (ExistParameter('L') && !Sum) {
      fCustomPrint(stderr,"***Warning: The loading factor lambda will not yield different results\n");
      fCustomPrint(stderr,"            from the base case since Pnl, Qnl, Pzl, and Qzl are zero\n");
      fCustomPrint(stderr,"            in all load buses.\n");
    } else if ((ExistParameter('H') || ExistParameter('c')) && !Sum) {
      fCustomPrint(stderr,"ERROR: The Homotopy Continuation Method will not yield different results\n");
      fCustomPrint(stderr,"       from the base case since Pnl, Qnl, Pzl, and Qzl are zero in all\n");
      fCustomPrint(stderr,"       load buses.\n");
      stopExecute(ERROREXIT);
    } else if (ExistParameter('C') && !Sum) {
      fCustomPrint(stderr,"ERROR: The Point of Collapse Method will not yield different results\n");
      fCustomPrint(stderr,"       from the base case since Pnl, Qnl, Pzl, and Qzl are zero in\n");
      fCustomPrint(stderr,"       all load buses.\n");
      stopExecute(ERROREXIT);
    }
  }

 /* -------------------- Initialize AC/DC Variables -------------------------- */
  if (ExistParameter('V')){
    Name=NameParameter('V');
    if (!NullName(Name) && (InputFile=OpenInput(Name))!=NULL) {
      for (;;) {
        if (fgets(Line,BUFLEN,InputFile)==NULL) break;
        if (Line[0]!='C') {
          if ((count=sscanf(Line,"%d %s",&N,BusName))!=2) count=0;
          else if (BusName[0]=='\"'||BusName[0]=='\'') {
             count=2;
             for(ptr=Line;*ptr!='\"'&&*ptr!='\'';ptr++); ptr++;
             for(i=0;*ptr!='\"'&&*ptr!='\''&&*ptr!='\n';BusName[i]= *ptr,i++,ptr++);
             BusName[i]='\0';
             ptr++;
          } else {
             for (i=1,ptr=Line;i<=2;i++) {
                if (ptr!=NULL && *ptr!='\0') for(;*ptr==' '&&*ptr!='\n';ptr++);
                if (ptr!=NULL && *ptr!='\0') for(;*ptr!=' '&&*ptr!='\n';ptr++);
            }
          }
          if (count){
             for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next)
                if(N==ACptr->Num||!strncmp(ACptr->Name,BusName,strlen(BusName))) break;
             if(ACptr!=NULL){
                if (ptr!=NULL) count=sscanf(ptr,"%lf %lf",&V,&d);
                else count=0;
                if (count==2) {
                   if (V>0 && ACptr->Cont!=NULL) ACptr->V=V;
                   ACptr->Ang=d*K3;
                }
             }
             for (DCptr=dataPtr->DCbus;DCptr!=NULL;DCptr=DCptr->Next)
             if(!strncmp(DCptr->Name,BusName,strlen(BusName))) break;
             if(DCptr!=NULL){
                if (ptr!=NULL) count=sscanf(ptr,"%s",Vars);
                else count=0;
                val=0;
                if(ExistParameter('d')) CustomPrint("Read initial values for DC bus %s -> %s\n",DCptr->Name,Vars);
                if (count) for(i=1;i<=strlen(Vars) && count;i=i+2) {
                   GetStr(Vars,i,2,2,variable);
                   if (ptr!=NULL && *ptr!='\0') for(;*ptr==' '&&*ptr!='\n';ptr++);
                   if (ptr!=NULL && *ptr!='\0') for(;*ptr!=' '&&*ptr!='\n';ptr++);
                   if (ptr!=NULL) {
                      count=sscanf(ptr,"%lf",&val);
                      if(ExistParameter('d')) CustomPrint("                                           %lf\n",val);
                   } else count=0;
                   if(!strcmp(variable,"PA") && strcmp(DCptr->Cont1,"PA") && strcmp(DCptr->Cont2,"PA")) DCptr->P=val;
                   if(!strcmp(variable,"AL") && strcmp(DCptr->Cont1,"AL") && strcmp(DCptr->Cont2,"AL")) DCptr->Alfa=val;
                   if(!strcmp(variable,"GA") && strcmp(DCptr->Cont1,"GA") && strcmp(DCptr->Cont2,"GA")) DCptr->Gamma=val;
                   if(!strcmp(variable,"VD") && strcmp(DCptr->Cont1,"VD") && strcmp(DCptr->Cont2,"VD")) DCptr->Vd=val;
                   if(!strcmp(variable,"ID") && strcmp(DCptr->Cont1,"ID") && strcmp(DCptr->Cont2,"ID")) DCptr->Id=val;
                   if(!strcmp(variable,"QA") && strcmp(DCptr->Cont1,"QA") && strcmp(DCptr->Cont2,"QA")) DCptr->Q=val;
                   if(!strcmp(variable,"AT") && strcmp(DCptr->Cont1,"AT") && strcmp(DCptr->Cont2,"AT")) DCptr->Tap=val;
                }
             }
          } else { ACptr=NULL; DCptr=NULL;}
          if (ACptr==NULL && DCptr==NULL) {
             fCustomPrint(stderr,"***Warning: Line-> %s",Line);
             fCustomPrint(stderr,"            will be ignored in file %s.\n",Name);
          }
        }
      }
      fclose(InputFile);
    }
    else if (NullName(Name)) {
      for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
        if (ACptr->Cont!=NULL) ACptr->V=1;
        ACptr->Ang=0;
      }
    }
  }

 /* -------------------- Read Generator Steady-State Data -------------------- */
  Name=NameParameter('3');
  Ngen=0;
  if (!NullName(Name) && (InputFile=OpenInput(Name))!=NULL) {
   for (;;) {
      Ra=Xd=Xq=IaMax=EqMax=EqMin=0;
      if (fgets(Line,BUFLEN,InputFile)==NULL) break;
      if (Line[0]!='C') {
        if (sscanf(Line,"%d %s",&N,BusName)!=2) count=0;
        else if (BusName[0]=='\"'||BusName[0]=='\'') {
           count=2;
           for(ptr=Line;*ptr!='\"'&&*ptr!='\'';ptr++); ptr++;
           for(i=0;*ptr!='\"'&&*ptr!='\''&&*ptr!='\n';BusName[i]= *ptr,i++,ptr++);
           BusName[i]='\0';
           ptr++;
           count+=sscanf(ptr,"%lf %lf %lf %lf %lf %lf",&Ra,&Xd,&Xq,&IaMax,&EqMax,&EqMin);
        } else {
           count=sscanf(Line,"%d %s %lf %lf %lf %lf %lf %lf",&N,BusName,&Ra,&Xd,&Xq,&IaMax,&EqMax,&EqMin);
        }
        if (!strcmp(BusName,"0")) BusName[0]='\n';
        if (count==8){
           for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
              if(N==ACptr->Num||!strncmp(ACptr->Name,BusName,strlen(BusName))) break;
           }
           if(ACptr!=NULL){
              if (strpbrk(ACptr->Type,"G,Q,E,S,V")==NULL || Xd<=0) {
                 fCustomPrint(stderr,"***Warning: Line-> %s",Line);
                 fCustomPrint(stderr,"            will be ignored in file %s.\n",Name);
                 if (Xd<=0) fCustomPrint(stderr,"            The value of Xd is less than or equal to zero.\n");
                 else       fCustomPrint(stderr,"            The bus is not a generator type (BQ,BG,BE,BV,BS).\n");
              }
              else {
#ifdef WINDOWS
                 ACptr->Gen=new GenModel;
#else
                 ACptr->Gen=(GenModel *) malloc (sizeof(GenModel));
                 if (ACptr->Gen==NULL) {
                    fclose(InputFile);
                    ErrorHalt("Insufficient memory to allocate steady-state Generator data.");
                    stopExecute(ERROREXIT);
                 }
#endif
                 ACptr->Gen->Ra=fabs(Ra);
                 if (Xd<Ra) {
                    fCustomPrint(stderr,"***Warning: The generator steady-state data for bus:%s\n",ACptr->Name);
                    fCustomPrint(stderr,"            has Ra > Xd.\n");
                 }
                 ACptr->Gen->Xd=Xd;
                 if (Xq==0) Xq=Xd;
                 if (Xq>Xd) {
                    fCustomPrint(stderr,"***Warning: The generator steady-state data for bus:%s\n",ACptr->Name);
                    fCustomPrint(stderr,"            has Xq > Xd.  The program will force Xq=Xd.\n");
                    Xq=Xd;
                 }
                 ACptr->Gen->Xq=fabs(Xq);
                   /*
                 if (IaMax!=0 && strpbrk(ACptr->Type,"S")) {
                    fCustomPrint(stderr,"***Warning: The IaMax limit in slack generator:%s\n",ACptr->Name);
                    fCustomPrint(stderr,"            will be assumed large.\n");
                    IaMax=9999999999.;
                 }
                 */
                 if (IaMax<=0) {
                    fCustomPrint(stderr,"***Warning: The generator steady-state data for bus:%s\n",ACptr->Name);
                    fCustomPrint(stderr,"            has IaxMax<=0.  The program will force IaMax=9999999999.\n");
                    IaMax=9999999999.;
                 }
                 ACptr->Gen->IaMax=IaMax;
                 if (EqMax<=0) {
                    fCustomPrint(stderr,"***Warning: The generator steady-state data for bus:%s\n",ACptr->Name);
                    fCustomPrint(stderr,"            has EqMax<=0.  The program will force EqMax=9999999999.\n");
                    EqMax=9999999999.;
                 }
                 ACptr->Gen->EqMax=EqMax;
                 if (EqMin>EqMax || EqMin<0) {
                    fCustomPrint(stderr,"***Warning: The generator steady-state data for bus:%s\n",ACptr->Name);
                    fCustomPrint(stderr,"            has EqMin>EqMax or EqMin<0.  The program will force EqMin=0.\n");
                    EqMax=0.;
                 }
                 ACptr->Gen->EqMin=EqMin;
                 V=ACptr->V;
                 d=ACptr->Ang;
                 Pg=ACptr->Pg;
                 Qg=ACptr->Qg;
                 Vr=V*cos(d);
                 Vi=V*sin(d);
                 Ir=(Vr*Pg+Vi*Qg)/(V*V);
                 Ii=(Vi*Pg-Vr*Qg)/(V*V);
                 Ia=sqrt(Ir*Ir+Ii*Ii);
                 if (Ia>IaMax) Ia=IaMax;
                 else if (Ia==0) Ia=1.;
                 apr=Vr-Xq*Ii;
                 api=Vi+Xq*Ir;
                 dg=atan2(api,apr);
                 Vq=cos(dg)*Vr+sin(dg)*Vi;
                 Vd=-sin(dg)*Vr+cos(dg)*Vi;
                 Iq=cos(dg)*Ir+sin(dg)*Ii;
                 Id=-sin(dg)*Ir+cos(dg)*Ii;
                 Eq=Vq-Xd*Id;
                 if (Eq>EqMax) Eq=EqMax;
                 else if (Eq<EqMin) Eq=EqMin;
                 ACptr->Gen->Eq=Eq;
                 ACptr->Gen->dg=dg;
                 ACptr->Gen->Vr=Vr;
                 ACptr->Gen->Vi=Vi;
                 ACptr->Gen->Ir=Ir;
                 ACptr->Gen->Ii=Ii;
                 ACptr->Gen->Vd=Vd;
                 ACptr->Gen->Vq=Vq;
                 ACptr->Gen->Id=Id;
                 ACptr->Gen->Iq=Iq;
                 ACptr->Gen->Ia=Ia;
                 Ngen++;
              }
           }
        } else ACptr=NULL;
        if (ACptr==NULL) {
            fCustomPrint(stderr,"***Warning: Line-> %s",Line);
            fCustomPrint(stderr,"            will be ignored in file %s.\n",Name);
        }
      }
   }
   fclose(InputFile);
  }

 /* -------------------- Read OH Load Model --------------------- */
  if (!ReadOHload(NameParameter('D'))) {
    if (!NullName(NameParameter('D'))) {
      fCustomPrint(stderr,"Error in file %s.  The program will assume constant\n",NameParameter('D'));
      fCustomPrint(stderr,"P-Q load models.\n");
    }
    if (!flagVloads) for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
      ACptr->Pn=ACptr->Pl;
      ACptr->a=ACptr->Pz=0;
      ACptr->Qn=ACptr->Ql;
      ACptr->b=ACptr->Qz=0;
    }
  }


 /* -------------------- Voltage/Lambda Options -------------------- */
  Bl=0;
  if (ExistParameter('v')) {
    RealParameter('v',&V,0.0,10.0);
    if (V>0.001) {
      ACptr=NULL;
      N=IntegerParameter('B',0,0,9999);
      if (N) {
         for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next)
         if((!strcmp(ACptr->Type,"B")||!strcmp(ACptr->Type,"BA")) && ACptr->Num==N) break;
         if (ACptr==NULL) {
            fCustomPrint(stderr,"***Warning: The program will ignore the number in -B option (not a PQ bus).\n");
         }
      }
      if (ACptr==NULL) {
         for(ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next)
         if(!strcmp(ACptr->Type,"B")||!strcmp(ACptr->Type,"BA")) break;
      }
      if (ACptr!=NULL) {
         BlPtr=ACptr;
         BlPtr->V=V;
         return(TRUE);
      }
      else {
         fCustomPrint(stderr,"***Warning: The program will ignore the -v option (there is no PQ bus).\n");
      }
    } else {
      fCustomPrint(stderr,"***Warning: The program will ignore the -v option (mag=0).\n");
    }
  }
  return(FALSE);
}



/* --------------------- ReadOHload ---------------------------------- */
#ifdef ANSIPROTO
BOOLEAN ReadOHload(char *File)
#else
BOOLEAN ReadOHload(File)
char *File;
#endif
 /* Read Ontatio Hydro load model using SSSP (OH) format. */
{
  ACbusData *ACptr;
  FILE *Input;
  char Line[BUFLEN],Name[13];
  int LineNum=0;
  BOOLEAN flag=FALSE,flagp=FALSE;
  INDEX Num,NA,Nlow,Nhigh;
  VALUETYPE Pn,Pz,Qn,Qz,a,b,MinMW,MinMVAR;


  if (!NullName(File) && (Input=OpenInput(File))!=NULL) {
    for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
      ACptr->Pn=0;
      ACptr->Pz=1;
      ACptr->a=0;
      ACptr->Qn=0;
      ACptr->Qz=1;
      ACptr->b=0;
    }
    for(;;){ /* Reading Loop */
      if (fgets(Line,BUFLEN,Input)==NULL) {
         ErrorHalt("Missing NLBS, EDATA, or END cards on OH load data file.");
         fclose(Input);
         return(FALSE);
      }
      LineNum++;

      /* -------------- Comments ------------- */
      if (!strncmp(Line,"!",1) || !strncmp(Line,"C",1) ) continue;

      /* ------------ OH load Model ---------- */
      else if (!strncmp(Line," NLBS",5)) flag=TRUE;
      else if (flag && (!strncmp(Line,"EDATA",5) || !strncmp(Line,"  END",5))) break;
      else if (flag) {
         GetStr(Line,1,12,12,Name);
         Num=GetInt(Line,1,5);
         if (strncmp(Name,"50000",5) && strncmp(Name,"60000",5)) {

         /* -------------- Definition of each bus ------------- */
            for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
               if(!strncmp(ACptr->Name,Name,strlen(Name))){
                  flagp=TRUE;
                  ACptr->Pn=(VALUETYPE) GetInt(Line,13,5)/100.0;
                  ACptr->Pz=1-ACptr->Pn;
                  ACptr->a=GetValue(Line,23,10,5);
                  ACptr->Qn=(VALUETYPE) GetInt(Line,18,5)/100.0;
                  ACptr->Qz=1-ACptr->Qn;
                  ACptr->b=GetValue(Line,33,10,5);
                  break;
               } else if(Num==ACptr->Num) {
                  flagp=TRUE;
                  ACptr->Pn=(VALUETYPE) GetInt(Line,6,5)/100.0;
                  ACptr->Pz=1-ACptr->Pn;
                  ACptr->a=GetValue(Line,16,10,5);
                  ACptr->Qn=(VALUETYPE) GetInt(Line,11,5)/100.0;
                  ACptr->Qz=1-ACptr->Qn;
                  ACptr->b=GetValue(Line,26,10,5);
                  break;
               }
            }
            if (ACptr==NULL) {
               fCustomPrint(stderr,"***Warning: Line #%d-> %s",LineNum,Line);
               fCustomPrint(stderr,"            will be ignored in file %s.\n",File);
            }
         }

         /* -------------- Definition of buses by range ------------------- */
         else {
            Pn=(VALUETYPE) GetInt(Line,6,5)/100.0;
            Pz=1-Pn;
            a=GetValue(Line,16,10,5);
            Qn=(VALUETYPE) GetInt(Line,11,5)/100.0;
            Qz=1-Qn;
            b=GetValue(Line,26,10,5);
            LineNum++;
            if (fgets(Line,BUFLEN,Input)==NULL || !strncmp(Line,"EDATA",5) || !strncmp(Line,"  END",5)) {
               fCustomPrint(stderr,"Range card is missing in line #%d of file %s.\n",LineNum,File);
               fclose(Input);
               return(FALSE);
            }
            NA=GetInt(Line,1,5);
            if (!strncmp(Name,"50000",5) && Narea<=1 ) {
               fCustomPrint(stderr,"***Warning: Area %d doen't exist (line #%d of file %s).\n",NA,LineNum,File);
               fCustomPrint(stderr,"            The program will assume these data for all system loads.\n");
            }
            Nlow=GetInt(Line,6,5);
            if (Nlow<=0) Nlow=1;
            Nhigh=GetInt(Line,11,5);
            if (Nhigh<=0) Nhigh=INT_MAX;
            MinMW=GetValue(Line,16,10,5);
            MinMW=MinMW/Sn;
            MinMVAR=GetValue(Line,26,10,5);
            MinMVAR=MinMVAR/Sn;
            for (ACptr=dataPtr->ACbus;ACptr!=NULL;ACptr=ACptr->Next) {
               if (!strncmp(Name,"60000",5) || Narea<=1 || ACptr->Area->N==NA) {
                  flagp=TRUE;
                  if (ACptr->Num<=Nhigh && ACptr->Num>=Nlow) {
                     if (ACptr->Pl>MinMW) {
                        ACptr->Pn=Pn;
                        ACptr->Pz=Pz;
                        ACptr->a=a;
                     }
                     if (ACptr->Ql>MinMVAR) {
                        ACptr->Qn=Qn;
                        ACptr->Qz=Qz;
                        ACptr->b=b;
                     }
                  }
               }
            }
         }
      }
    }
    fclose(Input);
  }
  return(flagp);
}

