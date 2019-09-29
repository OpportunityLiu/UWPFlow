/* constant.h 030390 */

#define INDEX int
#define VALENCE unsigned int
#define LONGINT long int

#define NORMALEXIT  0
#define WARNINGEXIT 1
#define ERROREXIT   2

// comment/remove comment theses lines based on environment
#define ANSIPROTO
//#define WINDOWS
//#define WINDOWS95

#ifndef maxNp1
#define maxNp1 1001
#define maxN 1000
#define ELEMENTVALUETYPE double
#define VALUETYPE double
#define SINGULARITYZERO 1.0e-20
#define BUFLEN 1000
#define PI 3.141592654
#define UPARROW 134
#define LEFTARROW 132
#define RIGHTARROW 135
#define DOWNARROW 133
#define ESCAPE 27
#define BACKSPACE 1
#define NUMBER 2
#endif

#ifndef BOOLEAN
/* For the other compilers, use either "unsigned" or "short int" */
#define BOOLEAN unsigned
#endif
#ifndef TRUE
#define TRUE (BOOLEAN) 1
#define FALSE (BOOLEAN) 0
#endif

#ifndef RAND_MAX
#define RAND_MAX INT_MAX
#endif

