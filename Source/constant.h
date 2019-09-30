/* constant.h 030390 */
#pragma once

using INDEX = int;
using VALENCE = unsigned int;
using LONGINT = long int;

constexpr const int NORMALEXIT = 0;
constexpr const int WARNINGEXIT = 1;
constexpr const int ERROREXIT = 2;

// comment/remove comment theses lines based on environment
#define ANSIPROTO
//#define WINDOWS
//#define WINDOWS95

#define maxN 1000
using ELEMENTVALUETYPE = double;
using VALUETYPE = double;
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

using BOOLEAN = bool;

constexpr BOOLEAN TRUE = true;
constexpr BOOLEAN FALSE = false;
