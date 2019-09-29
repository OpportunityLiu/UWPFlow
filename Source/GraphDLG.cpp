// GraphDLG.cpp : implementation file
//

#include "stdafx.h"
#include "Win UWPflow.h"
#include "GraphDLG.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// GraphDLG dialog


GraphDLG::GraphDLG(CWnd* pParent /*=NULL*/)
	: CDialog(GraphDLG::IDD, pParent)
{
	//{{AFX_DATA_INIT(GraphDLG)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
	maxX = maxY = 0;
	minX = minY = 100;
	show = false;
	
}


void GraphDLG::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(GraphDLG)
	DDX_Control(pDX, IDC_GRAPHCTRL, m_GraphCtrl);
	//}}AFX_DATA_MAP
}

void GraphDLG::CheckRange(double x, double y)
{
	if (x>maxX)
		maxX=x;
	else if (x<minX)
		minX = x;

	if(y>maxY)
		maxY = y;
	else if (y<minY)
		minY = y;
}

void GraphDLG::Reset()
{
	m_GraphCtrl.ClearGraph();

	while(m_GraphCtrl.GetAnnoCount()>0)
		m_GraphCtrl.DeleteAnnotation(0);

	minY = minX = 100;
	maxY = maxX = 0;

	show = false;

}

BEGIN_MESSAGE_MAP(GraphDLG, CDialog)
	//{{AFX_MSG_MAP(GraphDLG)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

