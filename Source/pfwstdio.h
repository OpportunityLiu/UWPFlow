//---------------------------------------------------------------------------
// Pflow for Windows:  Redefinition in C of several stdio.h routines
//                     and variables
//
// Claudio Canizares, Shu Zhang (c) 1996, 2006
// University of Waterloo
//---------------------------------------------------------------------------
#ifdef WINDOWS
#include <iostream>
#include <fstream>
#include <setjmp.h>
#include "StdAfx.h"
#else
#include "stdio.h"
#endif

// Define Screen data type to work with Windows Paint and some
// other variables
#define TABSPACES 4
#define BUFFER 512
#define argc Argc
#define argv Argv

#ifdef WINDOWS
extern jmp_buf exit_main;
#endif

//stops current execution. jumps out of main
void stopExecute(int status);

//replacement for fprintf/printf
//does fprintf in unix, does fprintf or print to screen in windows
int fCustomPrint(FILE *stream, const char *format,...);

//does printf in unix, does printf or print to screen in windows
int CustomPrint( const char *format,... );


/*
// Redefinition of FILE for input and output
// NOTES: - stdio.h cannot be used in Windows
//        - ofstream operations don't work properly; file names must
//          be used instead
//        - pointers to windows don't work either
typedef struct {
  char Name[100];
  fstream ios;
} FILE;
#define _FILE_DEFINED

FILE *stdout, *stdin, *stderr;

// Redefine:
// * Output to file FileOut using Format
int fprintf(FILE *FileOut,const char *Format,...);

// * Open a file
FILE *fopen(const char *FileName,const char *mode);
// * Associate an existing stream to a file
FILE *freopen(const char *FileName,const char *mode,FILE *File);
// * Close a file
int fclose(FILE *File);

/*typedef struct ScreenType {
   char *Line;
   ScreenType *Next;
} ScreenType;



#ifndef _PFLOW_WINDOWS
// Redefinition of some std* variables and command line arguments
// and exit label

extern int Argc;
extern char **Argv;
extern jmp_buf exit_main;
#endif


/*
// * Output to stdout using Format
int printf(const char *Format,...);
// * Output to string Buffer using Format
int sprintf_s(char *Buffer,const char *Format,...);
// * Read line up to N number of characters from input file FileIn
//   and put it into string Buffer
char *fgets(char *Buffer, int n, FILE *FileIn);
// * Read from string Buffer using Format
int sscanf(const char *Buffer, const char *Format, ...);
// * Read from stdin using Format
int scanf(const char *Format, ...);

*/
