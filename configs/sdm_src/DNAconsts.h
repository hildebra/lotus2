/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef _DNAconsts_h
#define _DNAconsts_h

//compile with multi  threading support
#define _THREADE //D

//do vector based DNA matching
#define _NEWMATCH

//sum up uc file to OTU abundance matrix in seed extension step
#define matrix_sum

//match barcodes based on maps
#define _fastBCmatch

//keep a map of dereplicated file
#define _MAPDEREPLICATE

//KHASH for faster / lower mem access
#define KHAS_H

//disable win warning about fopen
#define _CRT_SECURE_NO_WARNINGS 


//read gzip'd files using zlib.h
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea//d
#else
#define _gzipread
#endif

//DEBUG mode: more output
#define DEB//UG


#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
//#include <string.h>
#include <string.h>
#include <map>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <time.h>
#include <list>
#include <math.h>
#include <unordered_map>
#include <unordered_set>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <memory>

//#include <future>
#ifdef KHASH
#include "khash.hh"
#endif

#ifdef _gzipread
#include "gzstream.h"
#endif


#ifdef _THREADED
#include <thread>
#include <mutex>
#endif
//#include <future>

#ifdef _WIN32
//#include <boost/chrono/thread_clock.hpp>
#endif // _WIN32


static const float sdm_version = 1.50f;
static const char* sdm_status = "beta";


using namespace std;

static const double SAqualP[110] = {1.000000e+00,7.943282e-01,6.309573e-01,5.011872e-01,3.981072e-01,3.162278e-01,2.511886e-01,1.995262e-01,1.584893e-01,1.258925e-01
,1.000000e-01,7.943282e-02,6.309573e-02,5.011872e-02,3.981072e-02,3.162278e-02,2.511886e-02,1.995262e-02,1.584893e-02,1.258925e-02
,1.000000e-02,7.943282e-03,6.309573e-03,5.011872e-03,3.981072e-03,3.162278e-03,2.511886e-03,1.995262e-03,1.584893e-03,1.258925e-03
,1.000000e-03,7.943282e-04,6.309573e-04,5.011872e-04,3.981072e-04,3.162278e-04,2.511886e-04,1.995262e-04,1.584893e-04,1.258925e-04
,1.000000e-04,7.943282e-05,6.309573e-05,5.011872e-05,3.981072e-05,3.162278e-05,2.511886e-05,1.995262e-05,1.584893e-05,1.258925e-05
,1.000000e-05,7.943282e-06,6.309573e-06,5.011872e-06,3.981072e-06,3.162278e-06,2.511886e-06,1.995262e-06,1.584893e-06,1.258925e-06
,1.000000e-06,7.943282e-07,6.309573e-07,5.011872e-07,3.981072e-07,3.162278e-07,2.511886e-07,1.995262e-07,1.584893e-07,1.258925e-07
,1.000000e-07, 1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 
,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 ,1.000000e-07 
, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07
, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07, 1.000000e-07 };

static const char DNA_SPACE[15] = {'A','C','G','T','N','R','Y','M','K','W','S','B','D','H','V'};
static const int DNAinMemory = 5000;
static const unsigned int maxFileStreams = 500;
static const int RDBUFFER = 4096;

typedef unsigned int uint;
typedef unsigned long ulong;
typedef int qual_score; //used for quality scores in vectors

//seeding
static const float BestLengthRatio = 0.83f;
static const float RefLengthRatio = 0.9f;
static const qual_score MinQualDiff = 5;



void ini_DNAconstants();


// first base is [ACTG], second can be IUPAC
/*code	description
A	Adenine
C	Cytosine
G	Guanine
T	Thymine
U	Uracil
R	Purine (A or G)
Y	Pyrimidine (C, T, or U)
M	C or A
K	T, U, or G
W	T, U, or A
S	C or G
B	C, T, U, or G (not A)
D	A, T, U, or G (not C)
H	A, T, U, or C (not G)
V	A, C, or G (not T, not U)
N	Any base (A, C, G, T, or U)
*/



#endif