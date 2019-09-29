// Win UWPflow.cpp : Defines the class behaviors for the application.
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
// CWinUWPflowApp

BEGIN_MESSAGE_MAP(CWinUWPflowApp, CWinApp)
	//{{AFX_MSG_MAP(CWinUWPflowApp)
	ON_COMMAND(ID_APP_ABOUT, OnAppAbout)
	//}}AFX_MSG_MAP
	// Standard file based document commands
	ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowApp construction

CWinUWPflowApp::CWinUWPflowApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CWinUWPflowApp object

CWinUWPflowApp theApp;

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowApp initialization

BOOL CWinUWPflowApp::InitInstance()
{
	AfxEnableControlContainer();

	//register active controls (the graph control)
	RegisterActiveXControls();

	// Standard initialization
	// If you are not using these features and wish to reduce the size
	//  of your final executable, you should remove from the following
	//  the specific initialization routines you do not need.

#ifdef _AFXDLL
	Enable3dControls();			// Call this when using MFC in a shared DLL
#else
	Enable3dControlsStatic();	// Call this when linking to MFC statically
#endif

	// Change the registry key under which our settings are stored.
	// TODO: You should modify this string to be something appropriate
	// such as the name of your company or organization.
	SetRegistryKey(_T("Local AppWizard-Generated Applications"));

	LoadStdProfileSettings();  // Load standard INI file options (including MRU)

	// Register the application's document templates.  Document templates
	//  serve as the connection between documents, frame windows and views.

	CSingleDocTemplate* pDocTemplate;
	pDocTemplate = new CSingleDocTemplate(
		IDR_MAINFRAME,
		RUNTIME_CLASS(CWinUWPflowDoc),
		RUNTIME_CLASS(CMainFrame),       // main SDI frame window
		RUNTIME_CLASS(CWinUWPflowView));
	AddDocTemplate(pDocTemplate);

	// Parse command line for standard shell commands, DDE, file open
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);

	// Dispatch commands specified on the command line
	if (!ProcessShellCommand(cmdInfo))
		return FALSE;

	// The one and only window has been initialized, so show and update it.
	m_pMainWnd->ShowWindow(SW_SHOW);
	m_pMainWnd->UpdateWindow();
	return TRUE;
}

BOOL CWinUWPflowApp::RegisterActiveXControls()
{
	char szDrive[_MAX_DRIVE] = "";
	char szDir[_MAX_DIR] = "";
	char szFname[_MAX_FNAME] = "";
	char szExt[_MAX_EXT] = "";
	char szOcxFilePath[_MAX_PATH] = "";
	char g_szExePath[MAX_PATH];

	HRESULT hResult;
	FARPROC pfnRegServer;
	HINSTANCE hHandleToServerInstance = NULL;
	if(GetModuleFileName(GetModuleHandle("Win UWPflow.exe"), g_szExePath, _MAX_PATH))
	{
		_splitpath(g_szExePath, szDrive, szDir, szFname, szExt);
		_makepath(g_szExePath, szDrive, szDir, NULL, NULL);
	}
	strcpy_s(szOcxFilePath, g_szExePath);
	//strcpy_s(DBPath,g_szExePath);
	strcat_s(szOcxFilePath, "NTGraph.ocx");
	hHandleToServerInstance = ::LoadLibrary(szOcxFilePath);
	if(NULL == hHandleToServerInstance)
	{
		return FALSE;
	}
	if(hHandleToServerInstance)
	{
	/*(FARPROC&)*/pfnRegServer = ::GetProcAddress(hHandleToServerInstance, _T("DllRegisterServer"));
	}
	if(pfnRegServer)
	{
		hResult = (*pfnRegServer)();
	}
	::CoFreeLibrary(hHandleToServerInstance);
	return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
/*	
BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()*/

public:
	CAboutDlg();

	// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
		// No message handlers
	//}}AFX_MSG
	//DECLARE_MESSAGE_MAP()
	HTREEITEM m_hItemFirstSel;
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}


// App command to run the dialog
void CWinUWPflowApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}

extern CWinUWPflowView* myWindow;
BOOL CWinUWPflowApp::PreTranslateMessage(MSG* pMsg) 
{
	// Whenever user presses enter, run uwpflow
	if (pMsg->wParam == VK_RETURN && pMsg->message ==256)
	{
		myWindow->OnExecuteRunuwpflow();
	}
	return CWinApp::PreTranslateMessage(pMsg);
}


