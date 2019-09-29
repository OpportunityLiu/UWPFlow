// MainFrm.cpp : implementation of the CMainFrame class
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
// CMainFrame

IMPLEMENT_DYNCREATE(CMainFrame, CFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CFrameWnd)
	//{{AFX_MSG_MAP(CMainFrame)
	ON_WM_CREATE()
	ON_WM_ACTIVATE()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

static UINT indicators[] =
{
	ID_SEPARATOR,           // status line indicator
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

/////////////////////////////////////////////////////////////////////////////
// CMainFrame construction/destruction

CMainFrame::CMainFrame()
{
}

CMainFrame::~CMainFrame()
{
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CFrameWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	//create toolbar
	if (!m_wndToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC) ||
		!m_wndToolBar.LoadToolBar(IDR_MAINFRAME))
	{
		TRACE0("Failed to create toolbar\n");
		return -1;      // fail to create
	}
	
	//the following code puts a combo box in to the toolbar for commandline input
    #define SNAP_WIDTH 565 //the width of the commandline edit box

    int index;
    CRect rect;

    //A button has been placed on the toolbar, to be replaced with the edit box
    //get the index of that button's position in the toolbar
    index = 0;
    while(m_wndToolBar.GetItemID(index)!=ID_Cmd_Line) index++;

    //next convert that button to a seperator and get its position
    m_wndToolBar.SetButtonInfo(index, ID_Cmd_Line,
                               TBBS_SEPARATOR, SNAP_WIDTH);
    m_wndToolBar.GetItemRect(index, &rect);

	//reset the dimensions of the text box
	POINT BottomR;
	BottomR.x = rect.BottomRight().x;
	BottomR.y = rect.BottomRight().y +200;
	POINT TopL;
	TopL.x = rect.TopLeft().x;
	TopL.y = rect.TopLeft().y+2;

	rect.SetRect(TopL, BottomR);

    // then .Create the combo box and show it
    if (!m_wndToolBar.m_wndEdit.Create(WS_CHILD|WS_VISIBLE | CBS_AUTOHSCROLL |
										CBS_DROPDOWN,
										rect, &m_wndToolBar,
										IDC_CMD_LINE))
    {	
       TRACE0("Failed to create edit box\n");
       return FALSE;
    }
    m_wndToolBar.m_wndEdit.ShowWindow(SW_SHOW);

	//create a statusbar
	if (!m_wndStatusBar.Create(this) ||
		!m_wndStatusBar.SetIndicators(indicators,
		  sizeof(indicators)/sizeof(UINT)))
	{
		TRACE0("Failed to create status bar\n");
		return -1;      // fail to create
	}


	return 0;
}


BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	if( !CFrameWnd::PreCreateWindow(cs) )
		return FALSE;
	
	// Window Style
	cs.style = WS_OVERLAPPED | WS_CAPTION | WS_THICKFRAME | WS_SYSMENU |  WS_MINIMIZEBOX | WS_MAXIMIZEBOX | WS_MAXIMIZE;

	return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CMainFrame diagnostics

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CFrameWnd::Dump(dc);
}

#endif //_DEBUG

void CMainFrame::OnActivate(UINT nState, CWnd* pWndOther, BOOL bMinimized) 
{
	CFrameWnd::OnActivate(nState, pWndOther, bMinimized);
	
	// TODO: Add your message handler code here
	UpdateWindow();
	Invalidate();
}
