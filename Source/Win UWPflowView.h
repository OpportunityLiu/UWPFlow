// Win UWPflowView.h : interface of the CWinUWPflowView class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_WINUWPFLOWVIEW_H__D4B631F6_6E8F_4748_B405_0DCED76786F3__INCLUDED_)
#define AFX_WINUWPFLOWVIEW_H__D4B631F6_6E8F_4748_B405_0DCED76786F3__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define BUFFER 512

class CWinUWPflowView : public CScrollView
{
protected: // create from serialization only
	CWinUWPflowView();
	DECLARE_DYNCREATE(CWinUWPflowView)
	
// Attributes
public:
	CWinUWPflowDoc* GetDocument();
	CFont MyFont;
protected:
	//pointers to the toolbar and statusbar
	CMainToolBar* m_wndToolBar;
	CStatusBar* m_wndStatusBar;

	int wordHeight;	//current font


// Operations
public:
	void CWinUWPflowView::ResetScroll();
	BOOL CWinUWPflowView::Stop();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWinUWPflowView)
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	protected:
	virtual void OnInitialUpdate();
	virtual void OnActivateView(BOOL bActivate, CView* pActivateView, CView* pDeactiveView);
	//}}AFX_VIRTUAL
	


// Implementation
public:
	virtual ~CWinUWPflowView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif
	
protected:
	void SaveFile();
	void GetArguments(char *Line,char *Program);
	void CmRun();
	void ReDraw();
// Generated message map functions
protected:
	//{{AFX_MSG(CWinUWPflowView)
	afx_msg void OnHelpHelp();
	afx_msg void OnExecuteBatchscriptfile();
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDblClk(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnViewSize12();
	afx_msg void OnViewSize13();
	afx_msg void OnViewSize14();
	afx_msg void OnViewSize15();
	afx_msg void OnViewSize16();
	afx_msg void OnViewSize17();
	afx_msg void OnViewSize18();
	afx_msg void OnViewSize19();
	afx_msg void OnViewSize20();
	afx_msg void OnViewSize11();
	afx_msg void OnFileEdit();
	afx_msg void OnFileChangedirectory();
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnOptionsShowgraph();
	afx_msg void OnViewSize10();
public:
	afx_msg void OnExecuteRunuwpflow();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};



#ifndef _DEBUG  // debug version in Win UWPflowView.cpp
inline CWinUWPflowDoc* CWinUWPflowView::GetDocument()
   { return (CWinUWPflowDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_WINUWPFLOWVIEW_H__D4B631F6_6E8F_4748_B405_0DCED76786F3__INCLUDED_)
