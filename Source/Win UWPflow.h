// Win UWPflow.h : main header file for the WIN UWPFLOW application
//

#if !defined(AFX_WINUWPFLOW_H__D2982A79_60ED_4BDC_850E_9E20511469F8__INCLUDED_)
#define AFX_WINUWPFLOW_H__D2982A79_60ED_4BDC_850E_9E20511469F8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// CWinUWPflowApp:
// See Win UWPflow.cpp for the implementation of this class
//

class CWinUWPflowApp : public CWinApp
{
public:
	CWinUWPflowApp();
	BOOL RegisterActiveXControls(); //called to register graph control

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWinUWPflowApp)
	public:
	virtual BOOL InitInstance();
	virtual BOOL PreTranslateMessage(MSG* pMsg);
	//}}AFX_VIRTUAL

// Implementation
	//{{AFX_MSG(CWinUWPflowApp)
	afx_msg void OnAppAbout();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_WINUWPFLOW_H__D2982A79_60ED_4BDC_850E_9E20511469F8__INCLUDED_)
