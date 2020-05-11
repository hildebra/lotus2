#pragma once
#include "Matrix.h"

const string path2abundance = "/assemblies/metag/ContigStats/Coverage.pergene";
const string path2counts = "/assemblies/metag/ContigStats/Coverage.count_pergene";
const string pseudoAssMarker = "/assemblies/metag/longReads.fasta.filt.sto";

//typedef std::unordered_map<uint, uint>::iterator SmplAbunIT;
//typedef std::unordered_map<uint, uint> SmplAbun;


class ContigCrossHit
{
public:
	ContigCrossHit(int sn) :smplN(sn), CtgsPerSmpl(smplN,0),
		SmplNms(0,""){}
	void setSmplNms(vector<string> &sns) { SmplNms = sns; }
	void addHit(int Smpl, int Ctg);
	//~ContigCrossHit() {}
private:
	size_t smplN;
	vector<uint> CtgsPerSmpl;
	vector<string> SmplNms;
	
};

class GeneAbundance
{
public:
	GeneAbundance(const string, const string);
	inline smat_fl getAbundance(const string);
private:
	bool isPsAss;
	SmplAbun GeneAbu;
	//vector<int> ContigID;
	
};

class ClStr2Mat
{
	//class for gene catalog creation with cd-hit
public:
	ClStr2Mat(const string inF, const string outF, const string mapF, const string baseP, bool covCalc);
	virtual ~ClStr2Mat();
private:
	void read_map(const string,bool);
	void printVec(ofstream&, vector<smat_fl>&,const string&);

	vector<GeneAbundance*> GAs;
	ContigCrossHit* CCH;
	SmplOccurMult smpls;
	vector<string> smplLoc;
	vector<string> baseP;
	size_t smplN;
	size_t curr ;
};

