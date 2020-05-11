#pragma once
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <iterator>
#include <string>
#include <cstring>
#include <map>
#include <list>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <random>
#include <assert.h>
#include <unordered_map>
#include <numeric>
#include <future>
#include <mutex>
#include <chrono>
#include "options.h"

//#include <tchar.h>
//#include <string.h>

//multi threading
//#include <thread>
//#include <future>


const bool verbose=0;

#define LINELENGTH 20000;

typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
               // e.g. keep one global instance (per thread)
typedef double mat_fl;
typedef float smat_fl;


using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unordered_map <uint, uint> rare_map;

ulong thr_rng(unsigned long,MyRNG&);
std::istream& safeGetline(std::istream& is, std::string& t);
template<typename T> T getMedian(vector<T>& in){
	sort(in.begin(), in.end());
	size_t size = in.size();
	if (size == 0){ return (T)0; }
	if (size == 1){ return (in[0]) ; }
	if (size == 2){ return (in[0] + in[1]) / 2; }
	T median(in[size / 2]);
	if (size % 2 == 0)	{
		median = (in[size / 2 - 1] + in[size / 2]) / 2;
	}
	return median;
}
void lineCntOut(const string inF, const string outF, const string arg4);



inline std::string stringify(double x)
 {
   std::ostringstream o;
   o << x;
   return o.str();
 }
inline std::string itos(int number) {
	std::stringstream ss;
	ss << number;
	return ss.str();
}

class DivEsts{
public:
	DivEsts():richness(0),shannon(0),
		simpson(0),invsimpson(0),chao1(0),eve(0){}
	~DivEsts(){}
	//void print2file(const string);
	//data vectors
	vector<vector<long>> richness;
	vector<vector<double>> shannon,simpson,invsimpson,chao1,eve;
	string SampleName;
	int depth;
};
void printDivMat(const string outF, vector<DivEsts*>&, bool, options*);
void printRareMat(const string outF,const vector< rare_map>& rMat, vector< string >& sampleNames, vector < string >& rowId);
string printSimpleMap(const rare_map &vec, string outF, string id, vector<string> rowNames);
void reassembleTmpMat(vector<string> inF, vector< string > rowNames,vector< string > colNames, string outF);

class smplVec{
public:
	smplVec(const string, const int);
	smplVec(const vector<mat_fl>&, const int);
	~smplVec(){
		//delete[] arr;
	}
	void rarefy(vector<double> ,string o,int rep,DivEsts*, vector<vector<rare_map>>& RareSample,
		vector<string>& retCntsSampleName, string& skippedSample, vector<vector<vector<uint>>>* ,vector<vector<vector<uint>>>* , int=0,bool=false, bool=false);
	long getRichness(rare_map& cnts);
	long getRichness(const vector<unsigned int>&);
	//int maxSiz(){return vector<unsigned short>::max_size();}
	vector < string > getRowNames(){ return(IDs); }

private:
	int binarySearch(vector<float>,const float x);
	//void shuffle();
	void shuffle_singl();

	//diversity indices
	//method: 1=shannon, 2=simpson, 3=invsimpson
	vector<double> calc_div(const vector<uint>& vec,int meth=1, float base=2.718282f);
	vector <double> calc_div(rare_map& , int meth=1, float base=2.718282f);
	double calc_chao1(const vector<uint> & vec,int corrBias=1);
	double calc_chao1(rare_map& , int corrBias=1); //corrBias: 0/1
	double calc_eveness(const vector<uint>& vec);
	double calc_eveness(rare_map& );

	void print2File(const vector<unsigned int>&,const string);
	//unsigned short * arr;
	vector<string> IDs;
	vector<unsigned int> arr;
	double totSum;
	vector<MyRNG> rng_P;
	MyRNG rng;
	int num_threads;
	long richness;
	double Shannon;
	int numFeatures;

	//vector<float> vec;
};

void computeChao2(std::vector<vector<mat_fl>>& chao2, vector<vector<vector<uint>>>& abundInRow);
// compute ace or ice, depending on input data
void computeCE(vector<vector<mat_fl>>& CE, vector<vector<vector<uint>>>& abundInRow);


void writeGlobalDiv(options* opts, vector<vector<mat_fl>>& ICE, vector<vector<mat_fl>>& ACE, vector<vector<mat_fl>>& chao2, string outF);
//void writeGlobalDiv(vector<mat_fl>& ICE, vector<mat_fl>& ACE, vector<mat_fl>& chao2, string outF);
