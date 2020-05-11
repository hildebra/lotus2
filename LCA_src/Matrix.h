#pragma once
#include "libload.h"
#include "RefTax.h"
typedef double mat_fl;
static string __taxSepMat = ";";
static string __MatSep = "\t";

class Matrix
{
public:
	Matrix(int depth, vector<string> taxs, bool reportRead);
	~Matrix();
	void add(TaxObj*);
	void writeAllLevels(const string&);
//vars
	vector< vector< mat_fl > > mat;
	vector< string > colIDs; //rows = features, cols = tax levels
	vector<unordered_map<string, size_t> > rowIDs;
	bool readReport;
};

