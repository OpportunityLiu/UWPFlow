//{{AFX_INCLUDES()
#include "ntgraph.h"
//}}AFX_INCLUDES
#if !defined(AFX_GRAPHDLG_H__76F606AB_69EC_4E58_BE0E_36FDB6A2D8D6__INCLUDED_)
#define AFX_GRAPHDLG_H__76F606AB_69EC_4E58_BE0E_36FDB6A2D8D6__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// GraphDLG.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// GraphDLG dialog

class GraphDLG : public CDialog
{
// Construction
public:
	GraphDLG(CWnd* pParent = NULL);   // standard constructor
	double maxX, maxY, minX, minY;	//graph range
	void CheckRange(double x, double y);	//update the max/min x/y
	void Reset();
	int totalElements; //how many elements are in the graph
	bool show;

	//CNTGraph control taken from www.codeproject.com
	//source code to that control by Nikolai Teofilov
// Dialog Data
	//{{AFX_DATA(GraphDLG)
	enum { IDD = IDD_GRAPH_DLG };
	CNTGraph	m_GraphCtrl;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(GraphDLG)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(GraphDLG)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_GRAPHDLG_H__76F606AB_69EC_4E58_BE0E_36FDB6A2D8D6__INCLUDED_)
