// Win UWPflowDoc.h : interface of the CWinUWPflowDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_WINUWPFLOWDOC_H__82B714F9_E424_46DE_B039_6335736B6830__INCLUDED_)
#define AFX_WINUWPFLOWDOC_H__82B714F9_E424_46DE_B039_6335736B6830__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define BUFFER 512
#define MAXHEIGHT 30000

class CWinUWPflowDoc : public CDocument
{
protected: // create from serialization only
	CWinUWPflowDoc();
	DECLARE_DYNCREATE(CWinUWPflowDoc)
	CString text;	//the text in the document
	CString newText;	//there is only so much text that can be printed to screen
						//so after it reaches a certain height, we update this variable, and then later save text to a file and assign this to text
	CString allText; //storage for all the past text, to be printed in to a log file
	
	bool IsNewFile;	//has it been saved before?
	
	

// Attributes
public:
	CString FileName; //name of current document 
	int CountRuns;
	int height; //current height of the document in pixels
	int unPrintedLines; //indicates how much text has been placed in the document since the last screen refresh

// Operations
public:
	void mouseSave();	//to call OnFileSave when the right mouse button is pressed
	
// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWinUWPflowDoc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	virtual BOOL SaveModified();
	protected:
	
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CWinUWPflowDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	void SaveFile();

	
	

// Generated message map functions
protected:
	//{{AFX_MSG(CWinUWPflowDoc)
	afx_msg void OnFileSave();
	afx_msg void OnFileSaveAs();
	afx_msg void OnFileOpen();
	afx_msg void OnFileClear();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_WINUWPFLOWDOC_H__82B714F9_E424_46DE_B039_6335736B6830__INCLUDED_)
