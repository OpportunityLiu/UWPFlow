// Win UWPflowDoc.cpp : implementation of the CWinUWPflowDoc class
//

#include "stdafx.h"
#include "Win UWPflow.h"
#include "MainFrm.h"
#include "Win UWPflowDoc.h"
#include "Win UWPflowView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowDoc

IMPLEMENT_DYNCREATE(CWinUWPflowDoc, CDocument)

BEGIN_MESSAGE_MAP(CWinUWPflowDoc, CDocument)
	//{{AFX_MSG_MAP(CWinUWPflowDoc)
	ON_COMMAND(ID_FILE_SAVE, OnFileSave)
	ON_COMMAND(ID_FILE_SAVE_AS, OnFileSaveAs)
	ON_COMMAND(ID_FILE_OPEN, OnFileOpen)
	ON_COMMAND(ID_FILE_CLEAR, OnFileClear)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowDoc construction/destruction

CWinUWPflowDoc::CWinUWPflowDoc()
{
	// initialize variables
	text = "";
	newText="";
	allText ="";
	height = 0;
	CountRuns = 0;
	unPrintedLines = 0;
	IsNewFile = true;
	FileName ="";
}

CWinUWPflowDoc::~CWinUWPflowDoc()
{
}

BOOL CWinUWPflowDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowDoc serialization

void CWinUWPflowDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowDoc diagnostics

#ifdef _DEBUG
void CWinUWPflowDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CWinUWPflowDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowDoc commands

extern CWinUWPflowApp theApp;

//if there are changes, bring up a message box and ask the user to save
BOOL CWinUWPflowDoc::SaveModified() 
{
	if (IsModified())
	{
		switch(theApp.m_pMainWnd->MessageBox("Do you want to save Screen?", "Screen has changed",
						  MB_YESNOCANCEL | MB_ICONQUESTION)) {
		  case IDCANCEL:
			// Choosing Cancel means to abort the close -- return FALSE.
			return FALSE;

		  case IDYES:
			// Choosing Yes means to save the Screen.
			OnFileSave();
			if (IsModified()) return FALSE;

		  case IDNO:
			// Choosing No means clear variables
			SetModifiedFlag(false);
		}
	} 
	return TRUE;
}

//called when right mouse button is clicked, saves file
void CWinUWPflowDoc::mouseSave()
{
	OnFileSave();
}

//event for File->save
void CWinUWPflowDoc::OnFileSave() 
{
	if (IsNewFile) OnFileSaveAs();
	else           SaveFile();
}

//event for File->save as
void CWinUWPflowDoc::OnFileSaveAs() 
{
	char strFilter[] = { "Log Files (*.log)|*.log|All Files (*.*)|*.*||" };

	CFileDialog SaveFileDlg(FALSE, ".log", FileName, 0, strFilter);

	if( SaveFileDlg.DoModal() == IDOK )
	{
		FileName = SaveFileDlg.GetPathName();
		SaveFile();
	}	
}

//save file under its current filename
void CWinUWPflowDoc::SaveFile() 
{
	//crate an OutStream to write to specified file	
	FILE *OutStream = new FILE;
	OutStream = fopen(FileName, "w");
	
	//convert myDoc->text to char, and print it to OutStream
	LPCTSTR text = this->text;
	fprintf(OutStream, text); 
	
	fclose( OutStream );

	IsNewFile = false;
	SetModifiedFlag(false);
}

//opens a file, load text contents to myDoc and then display to screen
void CWinUWPflowDoc::OnFileOpen() 
{
	char line[BUFFER];

	CFileDialog OpenFileDlg (TRUE, "log", "*.log",
      OFN_FILEMUSTEXIST| OFN_HIDEREADONLY, "Log Files (*.log)|*.log|All Files (*.*)|*.*||", theApp.m_pMainWnd);

	if(OpenFileDlg.DoModal() == IDOK)
	{
		
		OnFileClear();
		FileName = OpenFileDlg.GetPathName();
		
		FILE *InStream = new FILE;
		InStream = fopen(FileName, "r");
		
		while (fgets( line, BUFFER, InStream ) != NULL){
			text.Insert(text.GetLength(), line);
		}

		fclose( InStream );

		//update screen
		theApp.m_pMainWnd->Invalidate();
		theApp.m_pMainWnd->UpdateWindow();
	}
}

//clears document/screen
void CWinUWPflowDoc::OnFileClear() 
{
	if (SaveModified())
	{
		//clear document
		text="";
		newText = "";
		height = 0;

		//reset variables
		SetModifiedFlag(false);
		IsNewFile = true;
		CountRuns = 0;
		FileName ="";

		extern CWinUWPflowView* myWindow;

		//redraw screen
		myWindow->SetScrollPos(SB_VERT, 0, true);
		myWindow->SetScrollSizes(MM_TEXT, CSize(0,0));
		myWindow->SetScrollRange(SB_VERT,0,0,true);
		myWindow->UpdateWindow();
		myWindow->Invalidate();
		
	}
	
}
