// Win UWPflowView.cpp : implementation of the CWinUWPflowView class
//
using System.Drawing;

#include "stdafx.h"
#include "Win UWPflow.h"
#include "Win UWPflowDoc.h"
#include "MainFrm.h"
#include "Win UWPflowView.h"
#include "constant.h"
#include "pfwstdio.h"
#include "process.h"
#include "GraphDLG.h"
#include "SETJMP.H"
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowView

IMPLEMENT_DYNCREATE(CWinUWPflowView, CScrollView)

BEGIN_MESSAGE_MAP(CWinUWPflowView, CScrollView)
	//{{AFX_MSG_MAP(CWinUWPflowView)
	ON_COMMAND(ID_HELP_HELP, OnHelpHelp)
	ON_COMMAND(ID_EXECUTE_BATCHSCRIPTFILE, OnExecuteBatchscriptfile)
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONDBLCLK()
	ON_WM_LBUTTONDOWN()
	ON_COMMAND(ID_VIEW_SIZE_12, OnViewSize12)
	ON_COMMAND(ID_VIEW_SIZE_13, OnViewSize13)
	ON_COMMAND(ID_VIEW_SIZE_14, OnViewSize14)
	ON_COMMAND(ID_VIEW_SIZE_15, OnViewSize15)
	ON_COMMAND(ID_VIEW_SIZE_16, OnViewSize16)
	ON_COMMAND(ID_VIEW_SIZE_17, OnViewSize17)
	ON_COMMAND(ID_VIEW_SIZE_18, OnViewSize18)
	ON_COMMAND(ID_VIEW_SIZE_19, OnViewSize19)
	ON_COMMAND(ID_VIEW_SIZE_20, OnViewSize20)
	ON_COMMAND(ID_VIEW_SIZE_11, OnViewSize11)
	ON_COMMAND(ID_FILE_EDIT, OnFileEdit)
	ON_COMMAND(ID_FILE_CHANGEDIRECTORY, OnFileChangedirectory)
	ON_WM_KEYDOWN()
	ON_COMMAND(ID_OPTIONS_SHOWGRAPH, OnOptionsShowgraph)
	ON_COMMAND(ID_VIEW_SIZE_10, OnViewSize10)
	ON_COMMAND(ID_EXECUTE_RUNUWPFLOW, OnExecuteRunuwpflow)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

//command line arguments
int Argc=0;
char **Argv;

//the previous commandline
char Options[400];

//definition of uwpflow main
int pfw_main(int argc, char **argv);

//the document to be updated and printed to screen
CWinUWPflowDoc *myDoc;

// Define exit label
jmp_buf exit_main;

//this
CWinUWPflowView *myWindow;

//the current directory
CString dir;

//graph dialog
GraphDLG* GraphDlg;	

/*--CWinUWPflowView construction/destruction--*/

CWinUWPflowView::CWinUWPflowView()
{
	myWindow = this;

	//initialize directory
	dir = __FILE__;
	dir.TrimRight("Win UWPflowView.cpp");

	//initialize graph dialog
	GraphDlg = new GraphDLG(this);
	GraphDlg->Create(IDD_GRAPH_DLG, this);

	wordHeight = 15; //default font

}



CWinUWPflowView::~CWinUWPflowView()
{

}

BOOL CWinUWPflowView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CScrollView::PreCreateWindow(cs);
}


void CWinUWPflowView::OnInitialUpdate()
{
	CScrollView::OnInitialUpdate();

	SetScrollSizes(MM_TEXT, CSize(0,0));

	//initialize myDoc
	myDoc = GetDocument();
	ASSERT_VALID(myDoc);


	//initialize font
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);

	//get the toolbar and statusbar
	TCHAR strClassName[255];
	CWnd* wndParent = GetParent();
	CWnd* wndChild = wndParent->GetWindow(GW_CHILD);
    while( wndChild != NULL )
    {
		// Get the class name of child control
		::GetClassName(wndChild->GetSafeHwnd(), strClassName, 255);

		if( _tcscmp(_T("ToolbarWindow32"), strClassName) == 0 )
		{	
			m_wndToolBar = (CMainToolBar*) wndChild;	
		}

		if( _tcscmp(_T("msctls_statusbar32"), strClassName) == 0 )
		{	
			m_wndStatusBar = (CStatusBar*) wndChild;
		}
		
		wndChild = wndChild->GetNextWindow();
	}
		//intialize streams
		fclose(stderr);
		stderr->_file=-1;
		fclose(stdout);
		stdout->_file=-1;
}

//paints the screen
//automatically called (via OnPaint())
void CWinUWPflowView::OnDraw(CDC* pDC)
{
	CString currLine="";
	myDoc->height=0;
	int index=0; //current position of currLine 

	pDC->SelectObject(MyFont);

	for (int i=0; i<myDoc->text.GetLength(); i++)
	{
		if (myDoc->text[i]=='\n')
		{
			pDC->TextOut(0,myDoc->height, currLine);

			//start drawing on the next line
			myDoc->height += wordHeight;

			//reset currLine
			currLine ="";
			index = 0;
		}
		else 
		{
			currLine.Insert(index++,myDoc->text[i]);
		}
	}
	
	if(myDoc->height>=MAXHEIGHT)
	{
		//save contents in to file
		FILE *OutStream = new FILE;
		OutStream = fopen("tmp.log", "w");
		LPCTSTR text = myDoc->allText;
		fprintf(OutStream, text); 
		fclose( OutStream );
		this->m_wndStatusBar->SetPaneText(0,"Too much text; the top half has been dumped into tmp.log", true);

		myDoc->text= myDoc->newText;
		myDoc->newText = "";
		OnDraw(pDC); //call this function again, with the updated text
	}

	pDC->TextOut(0,0, currLine);

	myDoc->unPrintedLines = 0;
}

//reset scrollbar
void CWinUWPflowView::ResetScroll()
{
	SetScrollSizes(MM_TEXT, CSize(0,myDoc->height+400));
	SetScrollPos(SB_VERT, myDoc->height+400, true);
}

/*--CWinUWPflowView diagnostics--*/

#ifdef _DEBUG
void CWinUWPflowView::AssertValid() const
{
	CScrollView::AssertValid();
}

void CWinUWPflowView::Dump(CDumpContext& dc) const
{
	CScrollView::Dump(dc);
}

CWinUWPflowDoc* CWinUWPflowView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CWinUWPflowDoc)));
	return (CWinUWPflowDoc*)m_pDocument;
}
#endif //_DEBUG

/*--Tools--*/

// Define command line arguments Argc and Argv
void CWinUWPflowView::GetArguments(char *Line,char *Program)
{
  int i,j,k;
  char string[200];
  bool FirstBlank=true;

  //clean up
  if (Argc!=0) {
    for(i=0;i<Argc;i++) 
		delete[] Argv[i];
    delete[] Argv;
    Argc=0;
  }
  i=0;

  //count the number of arguments (fill in Argc)
  while(Line[i]!='\0') {
    if (Line[i]!=' ' && FirstBlank) {
      Argc++;
      FirstBlank=FALSE;
    }
    else if (Line[i]==' ') FirstBlank=TRUE;
    i++;
  }

  //fill in Argv
  if (Argc!=0) {
    Argc++;
    Argv = new char* [Argc+1];
    Argv[0]=new char [strlen(Program)+1];
    strcpy(Argv[0],Program);
    Argv[Argc]=NULL;
    FirstBlank=TRUE;
    for(i=k=0,j=1;;i++) {
     if (Line[i]!=' ' && Line[i]!='\0') {
       string[k]=Line[i];
       k++;
       FirstBlank=FALSE;
     } else {
       if (!FirstBlank) {
         string[k]='\0';
         k++;
		 if (strcmp(string, "uwpflow")==0)
			 Argc--;	//ignore "uwpflow"
		 else{
			 Argv[j]=new char [k];
			 if (string[0]=='<' || string[0]=='>') strcpy(Argv[j],&string[1]);
			 else strcpy(Argv[j],string);
			 j++;
		 }
         FirstBlank=TRUE;
       }
       k=0;
     }
     if (Line[i]=='\0') break;
    }
  }
}

/*--Events--*/

//display help
void CWinUWPflowView::OnHelpHelp() 
{
	ShellExecute(m_hWnd, "open","Win UWPflow.chm", NULL, NULL, SW_SHOWNORMAL);
}


// Run UWPflow
void CWinUWPflowView::OnExecuteRunuwpflow() 
{
	m_wndToolBar->m_wndEdit.GetWindowText(Options, 400);

	//insert new commandline to the listbox
	m_wndToolBar->m_wndEdit.DeleteString(m_wndToolBar->m_wndEdit.FindStringExact(-1, Options));
	m_wndToolBar->m_wndEdit.InsertString (0, Options);
	m_wndToolBar->m_wndEdit.SetCurSel(0);

	GetArguments(Options,"uwpflow");

	int i;

	if(Argc!=0) {
		if (myDoc->CountRuns!=0) {
		  CustomPrint("\n");
		  for (i=1; i<BUFFER-1; i++) CustomPrint("%c",'-');
		  CustomPrint("\n");
		}
		CustomPrint("RUN #%d: ",++myDoc->CountRuns);
		for(i=0; i<Argc; i++) 
			CustomPrint("%s ",Argv[i]);

		CustomPrint("\n");

		//update screen
		Invalidate();
		UpdateWindow();


		m_wndStatusBar->SetPaneText(0,"Press [Esc] to stop processing", true);

		pfw_main(Argc, Argv);

		//reset streams
		fclose(stderr);
		stderr->_file=-1;
		fclose(stdout);
		stdout->_file=-1;

		myDoc->SetModifiedFlag(true);

		CustomPrint("\n END OF RUN #%d\n",myDoc->CountRuns);

		//update screen
		UpdateWindow();
		ResetScroll();
		Invalidate();	
	}	
}

//opens a batch/script file and execute uwpflow commands if there are any in it
void CWinUWPflowView::OnExecuteBatchscriptfile() 
{
	fclose(stdin);
	FILE *InputFile;
	char Buffer[BUFFER],*string;
	int i;
	bool end=false;
	CFileDialog OpenFileDlg (TRUE, "bat", "*.bat",
      OFN_FILEMUSTEXIST| OFN_HIDEREADONLY, "Batch/Script Files (*.bat)|*.bat|All Files (*.*)|*.*||", this);

	//set it to the correct default directory
	OpenFileDlg.m_ofn.lpstrInitialDir=dir.GetBuffer(dir.GetLength());

	if(OpenFileDlg.DoModal() == IDOK)
	{
		
		CString filename = OpenFileDlg.GetPathName();

		//reset directory
		dir="";

		InputFile=fopen(filename,"r");		
		if (InputFile==NULL) MessageBox("Unable to open file", "File Error", MB_OK | MB_ICONEXCLAMATION);
		else {
			for (i=1; i<BUFFER-1; i++) CustomPrint("%c",'-');
			myDoc->SetModifiedFlag(true);
			CustomPrint("RUN batch/script file: %s\n",filename);
			while (fgets(Buffer,BUFFER,InputFile)!=NULL && !end) if (Buffer[0]!='@'){
				if((string=strstr(Buffer,"echo"))!=NULL) {
					if ((string=strchr(string,' '))!=NULL) 
						fCustomPrint(stdout, "%s",&string[1]);
				}
				else if(strstr(Buffer,"rem")!=NULL) 
				{ 
					continue;
				}
				else if(strstr(Buffer,"pause")!=NULL) {

					//update screen
					UpdateWindow();
					this->ResetScroll();
					Invalidate();

					switch(MessageBox("Do you want to continue?", "Pause", MB_YESNO | MB_ICONQUESTION)) {
					case IDYES:
						continue;

					case IDNO:
						end = true;
					}
					CustomPrint("\n");
				}
				else if((string=strstr(Buffer,"uwpflow"))!=NULL ||
					(string=strstr(Buffer,"UWPflow"))!=NULL ||
					(string=strstr(Buffer,"UWPFLOW"))!=NULL ) 
				{
					string=strchr(string,' ');

					if (string[strlen(string)-1]==' '||string[strlen(string)-1]=='\n')
						string[strlen(string)-1]='\0';
					else
						string[strlen(string)]='\0';

					strcpy_s(Options,&string[1]);
					//set the textbox to the last commandline
					m_wndToolBar->m_wndEdit.SetWindowText(Options);
					m_wndToolBar->UpdateWindow();
					OnExecuteRunuwpflow();
				}
				else if((string=strstr(Buffer,"tomatlab"))!=NULL ||
					(string=strstr(Buffer,"Tomatlab"))!=NULL ||
					(string=strstr(Buffer,"TOMATLAB"))!=NULL ) 
				{
					CustomPrint("%s",string);
					string=strchr(string,' ');
					string[strlen(string)-1]='\0';
					strcpy_s(Options,&string[1]);
					GetArguments(Options,"tomatlab");
					if (_spawnvp(P_NOWAIT,"tomatlab.bat",Argv)==-1) {
						m_wndStatusBar->SetPaneText(0,"Problems with Tomatlab routine",true);
					}
				}
				else if((string=strstr(Buffer,"maxim"))!=NULL ||
					(string=strstr(Buffer,"Maxim"))!=NULL ||
					(string=strstr(Buffer,"MAXIM"))!=NULL ) 
				{
					CustomPrint("%s",string);
					string=strchr(string,' ');
					string[strlen(string)-1]='\0';
					strcpy_s(Options,&string[1]);
					GetArguments(Options,"maxim");
					if (_spawnvp(P_NOWAIT,"maxim.exe",Argv)==-1) {
						m_wndStatusBar->SetPaneText(0,"Problems with Maxim routine",true);
					}
				}
				else if((string=strstr(Buffer,"awk"))!=NULL ||
					(string=strstr(Buffer,"Awk"))!=NULL ||
					(string=strstr(Buffer,"AWK"))!=NULL ) 
				{
					CustomPrint("%s",string);
					string=strchr(string,' ');
					string[strlen(string)-1]='\0';
					strcpy_s(Options,&string[1]);
					GetArguments(Options,"awk");
					if (_spawnvp(P_WAIT,"awk.exe",Argv)==-1) {
						m_wndStatusBar->SetPaneText(0,"Problems with Awk routine",true);
					}
				}
				else if((string=strstr(Buffer,".bat"))!=NULL ||
					(string=strstr(Buffer,".Bat"))!=NULL ||
					(string=strstr(Buffer,".BAT"))!=NULL ) 
				{
					CustomPrint("%s",Buffer);
					string=strchr(string,' ');
					string[strlen(string)-1]='\0';
					strcpy_s(Options,&string[1]);
					GetArguments(Options,"BATCH");
					for (i=0; i<strlen(Buffer); i++) {
						if (Buffer[i]!='.') string[i]=Buffer[i];
						else { string[i]='\0'; break; }
					}
					if (_spawnvp(P_WAIT,string,Argv)==-1) {
						m_wndStatusBar->SetPaneText(0,"Problems with Batch file",true);
					}
				}
				else CustomPrint("%s",Buffer);
			}
		}
		CustomPrint("\nEND of batch/script file\n");
		UpdateWindow();
		ResetScroll();
		Invalidate();
		
		fclose(InputFile);

		m_wndStatusBar->SetPaneText(0,"Ready", true);
	}	
}

//save file when right mouse button is pressed
void CWinUWPflowView::OnRButtonDown(UINT nFlags, CPoint point) 
{
	myDoc->mouseSave();
}

//overwrite default right button click event with nothing
//(otherwise a pop up menu will appear)
void CWinUWPflowView::OnRButtonDblClk(UINT nFlags, CPoint point) 
{
	//nothing
}

void Highlight ( bool Highlight )
{

}

//when left mouse button is pressed
void CWinUWPflowView::OnLButtonDown(UINT nFlags, CPoint point) 
{
	
}

//keyboard events
void CWinUWPflowView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	//run commandline when enter is pressed
	if (nChar == VK_RETURN)
		OnExecuteRunuwpflow(); 
	CScrollView::OnKeyDown(nChar, nRepCnt, nFlags);
}

//open a file with notepad notepad
void CWinUWPflowView::OnFileEdit() 
{
	CString filename;

	CFileDialog OpenFileDlg (TRUE, "wsc", "*.wsc",
      OFN_FILEMUSTEXIST| OFN_HIDEREADONLY, "WSCC Files (*.wsc)|*.wsc|IEEE Files (*.cf)|*.cf|All Files (*.*)|*.*||", this);

	//set it to the correct default directory
	OpenFileDlg.m_ofn.lpstrInitialDir=dir.GetBuffer(dir.GetLength());


	if(OpenFileDlg.DoModal() == IDOK)
	{
		//get the file name
		filename = OpenFileDlg.GetPathName();
		
		//open that file with notepad
		filename.Insert(0,"notepad.exe ");
		if (WinExec(filename,SW_NORMAL)==(0|ERROR_BAD_FORMAT|ERROR_FILE_NOT_FOUND|ERROR_PATH_NOT_FOUND))
			m_wndStatusBar->SetPaneText(0,"Problems with opening the file", true);
	}
	
}

//Size changes

void CWinUWPflowView::OnViewSize10() 
{	
	wordHeight=10;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();
}

void CWinUWPflowView::OnViewSize12() 
{
	wordHeight=12;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();
	

}

void CWinUWPflowView::OnViewSize13() 
{
	wordHeight=13;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);	

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();
	
}

void CWinUWPflowView::OnViewSize14() 
{
	wordHeight=14;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);
	
	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();
	
}

void CWinUWPflowView::OnViewSize15() 
{
	wordHeight=15;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);	

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();
	
}

void CWinUWPflowView::OnViewSize16() 
{
	wordHeight=16;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);
	
	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();

}

void CWinUWPflowView::OnViewSize17() 
{
	wordHeight=17;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);	

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();

}

void CWinUWPflowView::OnViewSize18() 
{
	wordHeight=18;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);	

	//OnDraw(

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();

}

void CWinUWPflowView::OnViewSize19() 
{
	wordHeight=19;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);	

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();
	
}

void CWinUWPflowView::OnViewSize20() 
{
	wordHeight=20;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);	

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();
	
	
}

void CWinUWPflowView::OnViewSize11() 
{
	wordHeight=11;

	//update font
	MyFont.DeleteObject();
	MyFont.CreateFont(wordHeight,0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,0,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,
        FIXED_PITCH|FF_DONTCARE,0);	

	//update screen
	OnPaint();
	UpdateWindow();
	//ResetScroll();
	Invalidate();
	
}

//whether or not to show the graph
void CWinUWPflowView::OnOptionsShowgraph() 
{
	GraphDlg->ShowWindow(SW_SHOW);
}

//checks if escap key is pressed (indicating that the program should stop)
BOOL CWinUWPflowView::Stop()
{
  MSG msg;

  if (PeekMessage(&msg,NULL,0,0,PM_REMOVE)) {
    if (msg.wParam==VK_ESCAPE) {
      CustomPrint("\nWARNING:  Case interrupted by user.  Some files might\n");
      CustomPrint("          still be open; this could cause some problems.\n");
      Invalidate();
      UpdateWindow();
      return TRUE;
    } else return FALSE;
  }
  return FALSE;
}

//change the current directory
void CWinUWPflowView::OnFileChangedirectory() 
{

    CString strDirectory;

    LPMALLOC pMalloc;
    //set up BROWSEINFO struct for the ShBrowseForFolder dialog box
    BROWSEINFO bi;
    //Gets the Shell's default allocator (260)
    char buf[MAX_PATH];
    LPITEMIDLIST pidl;
    // Get help on BROWSEINFO struct - it's got all the bit settings.
    bi.hwndOwner = GetSafeHwnd();
    bi.pidlRoot = NULL;
    bi.pszDisplayName = NULL;
    bi.lpszTitle = _T("Change Directory");
    bi.ulFlags = BIF_RETURNFSANCESTORS | BIF_RETURNONLYFSDIRS;
    bi.lpfn = NULL;
    bi.lParam = 0;

    if (::SHGetMalloc(&pMalloc) == NOERROR)
    {
            if ((pidl = ::SHBrowseForFolder(&bi)) != NULL)
            {
                    if (::SHGetPathFromIDList(pidl, buf))
                    {
                            dir.Empty();
							dir.Insert(0,buf);
							if (dir.Right(1) != "\\")
								dir.Insert(dir.GetLength(),"\\");
                    }
                    pMalloc->Free(pidl);
            }
            pMalloc->Release();
     }
	
	
}

void CWinUWPflowView::OnActivateView(BOOL bActivate, CView* pActivateView, CView* pDeactiveView) 
{
	// TODO: Add your specialized code here and/or call the base class
	
	CScrollView::OnActivateView(bActivate, pActivateView, pDeactiveView);

	UpdateWindow();
	
	//This stops the program from putting the scroll bar to the end everytime
	//ResetScroll();
	
	Invalidate();
}
