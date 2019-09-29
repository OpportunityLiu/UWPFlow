#define WINVER 0x0601
#define _WIN32_WINNT_ 0x0601

//---------------------------------------------------------------------------
// Pflow for Windows:  Redefinition in C of several stdio.h routines
//                     and variables, and redefinition of exit in stdlib.h
//
// Claudio Canizares, Shu Zhang (c) 1996, 2006
// University of Waterloo
//---------------------------------------------------------------------------

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <time.h>
#include <stdio.h>
#include "pfwstdio.h"
#ifdef WINDOWS
#include <process.h>
#include "Win UWPflow.h"
#include "MainFrm.h"
#include "Win UWPflowDoc.h"
#include "Win UWPflowView.h"
#endif
#include "constant.h"

#ifdef WINDOWS
//pointers to the document and the window
extern CWinUWPflowDoc *myDoc;
extern CWinUWPflowView *myWindow;
#endif

//do fprintf if in window mode or if stream is initialized
//otherwise, print it to the document (and then to screen)
int fCustomPrint(FILE *stream, const char *format,...)
{
	char Buffer[BUFFER];
	int l=0;
	va_list ap;

	va_start(ap, format);
	l=vsprintf_s(Buffer,format,ap);
	va_end(ap);

#ifdef WINDOWS
	if(myWindow->Stop()) {
		if (stream!=NULL)
			fclose(stream);
		stopExecute(1);	
	}
#endif
#ifndef WINDOWS
	return fprintf(stream, Buffer);

#else
	if( stream!=NULL && stream->_file != -1 )
	{
		return fprintf(stream, Buffer);
	}
	else
	{
		myDoc->text.Insert(myDoc->text.GetLength(), Buffer);
		myDoc->allText.Insert(myDoc->allText.GetLength(), Buffer);

		if (myDoc->height>(MAXHEIGHT/2))
			myDoc->newText.Insert(myDoc->newText.GetLength(), Buffer);
			
		if(strstr(Buffer, "\n")!=NULL)
			myDoc->unPrintedLines++;

		//refresh screen after every 10 lines
		if(myDoc->unPrintedLines  >= 10)
		{
			myWindow->UpdateWindow();
			myWindow->ResetScroll();			
			myWindow->Invalidate();
			myDoc->unPrintedLines=0;
		}

		return 0;
	}

	

#endif	
}

//do printf if in window mode or if stdout is initialized
//otherwise, print it to the document (and then to screen)
int CustomPrint( const char *format,... )
{
	char Buffer[BUFFER];
	int l=0;
	va_list ap;

	va_start(ap, format);
	l=vsprintf_s(Buffer,format,ap);
	va_end(ap);

#ifndef WINDOWS
	return printf(Buffer);

#else
	if( stdout!=NULL && stdout->_file!=-1 )
	{
		return printf(Buffer);
	}
	else
	{
		myDoc->text.Insert(myDoc->text.GetLength(), Buffer);
		myDoc->allText.Insert(myDoc->allText.GetLength(), Buffer);

		if (myDoc->height>(MAXHEIGHT/2))
			myDoc->newText.Insert(myDoc->newText.GetLength(), Buffer);

		if(strstr(Buffer, "\n")!=NULL)
			myDoc->unPrintedLines++;

		//refresh screen after every 10 lines
		if(myDoc->unPrintedLines >= 10)
		{
			myWindow->UpdateWindow();
			myWindow->ResetScroll();
			myWindow->Invalidate();
			myDoc->unPrintedLines=0;
		}

		
		return 0;
	}
#endif
}

// Definition of Pflow clean up routine
void CleanUp();

//Stops current execution
void stopExecute(int status)
{
#ifdef WINDOWS
	//get the status bar and update it
	TCHAR strClassName[255];
	CWnd* wndParent = myWindow->GetParent();
	CWnd* wndChild = wndParent->GetWindow(GW_CHILD);
    while( wndChild != NULL )
    {
		// Get the class name of child control
		::GetClassName(wndChild->GetSafeHwnd(), strClassName, 255);

		if( _tcscmp(_T("msctls_statusbar32"), strClassName) == 0 )
		{	
			CStatusBar* m_wndStatusBar = (CStatusBar*) wndChild;
			if (status<0 || status==2) m_wndStatusBar->SetPaneText(0, "Pflow ERROR", true);
			else if (status==1)   m_wndStatusBar->SetPaneText(0,"Pflow WARNING",true);
			break;
		}
		
		wndChild = wndChild->GetNextWindow();
	}
	
	fclose(stdout);
	fclose(stderr);
	fclose(stdin);
	CleanUp();

	myWindow->Invalidate();
	myWindow->UpdateWindow();
	longjmp(exit_main,1);
#else
     exit(status);
#endif   
}

