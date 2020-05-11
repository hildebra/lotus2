//libload.h
//contains constants and the std libs being loaded
#pragma once
#include <stdio.h>
#include <string>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <iterator>
#include <string>
#include <map>
#include <list>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <random>
#include <assert.h>
#include <unordered_map>



#define paral//lel
#ifdef parallel
#include <future>
#endif // parallel




#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea //d
#else
#define _gzipread
#endif

#define LCAd//ebg

#define DE//BUG

#ifdef _gzipread
#include "gzstream.h"
#endif

using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;


//some static vars..
static int __default_depth = 7;
static string __unkwnTax = "";
static string __unkwnTaxWR = "?";
static string __defaultTaxSep = "\t";


