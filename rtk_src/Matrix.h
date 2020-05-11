#pragma once
#include "IO.h"
#include "options.h"

typedef std::map<std::string, int> GeneIDidx;
//contains gene ID, taxa
typedef std::map<std::string, vector<string>> LvlUp;

typedef std::unordered_map<string, smat_fl>::iterator SmplAbunIT;
typedef std::unordered_map<string, smat_fl> SmplAbun;
typedef std::unordered_map<string, vector<int> >::iterator SmplOccurITmult;
typedef std::unordered_map<string, vector<int> > SmplOccurMult;
typedef std::unordered_map<string, int >::iterator SmplOccurIT;
typedef std::unordered_map<string, int> SmplOccur;
typedef std::unordered_map<string, int> ModOccur;


mat_fl median(std::vector<mat_fl> vec,bool ignoreZeros=false);
void vecPurge(vector<vector<mat_fl>>& vec, mat_fl val);//removes val from each entry in vec<vec<>>
string join(const vector<string>& in, const string &X);


class column{
	public:
		double colsum;
		string id;
};

class HMat
{
public:
	HMat(string L, vector<string> Samples, vector<string> Features);
	~HMat(){}
	//get
	//unsigned long operator [](int i) const    { return registers[i]; }
	//set
	void set(string kk, int j, mat_fl v);

	void print(ofstream&);

private:
	GeneIDidx Feat2mat;
	string LvlName;

	vector<string> FeatureNs, SampleNs;
	vector< mat_fl > empty;
	vector< vector< mat_fl > > mat;
};

class SparseMatrix
{
public:
	SparseMatrix();
	~SparseMatrix();
	void addCount(string smpl, int row, smat_fl abund);
	//void newRow(string x) { mat.push_back(vector<SmplAbun>(0)); rowIDs.push_back(x); }
	void newRow(void) { mat.push_back(SmplAbun(0)); }
private:
	vector<SmplAbun> mat;
	SmplOccur colNames;
	vector<string> rowIDs;
};

class VecFiles
{
public:
	VecFiles(const string inF, const string outF, const string xtra);
	~VecFiles(){}
private:
	void readVecFile(const string inF);
	int getIdx(const string&);

	//VARIABLES
	vector<string> infiles;
};



struct ModStep
{
public:
	ModStep() :alternates(0), redundancy(0){}
	ModStep(const string&, bool&, vector<string>&);
	void getAllKOs(list<string>&);
	void setRedund(ModOccur& m);
	void abundParts(const vector<mat_fl>& v, const unordered_map<string, int>& IDX,
		vector<mat_fl>&, vector<bool>&, vector<string>&,
		float hitComplRatio =0.8, int redund=0);



	//e.g. alt[0][0] = KO001 and requires alt[0][1] =KO002 or alternatively only alt[1][0]=K0003
	vector< vector< string > > alternates;
	//how redundant is each KO and the different steps (basically occurence of KOs across DB)
	vector< vector <int> > redundancy;
};

class Module
{
public:
	Module(vector<string>& n);
	void getAllKOs(list<string>& r) { for (size_t i = 0; i < steps.size(); i++) { steps[i].getAllKOs(r); } }
	void setReddundancy(ModOccur& m) { for (size_t i = 0; i < steps.size(); i++) { steps[i].setRedund(m); } }
	mat_fl pathAbundance(const vector<mat_fl>& v, const unordered_map<string, int>& IDX,
		const int redun, const float PathwCompl, const float enzymCompl, string & savedNmsKO,float & modScoreOut);
	string name;
	string description;
	vector<ModStep> steps;
//handles module recurrency
	vector<string> submods;
	bool containsMods,usedInOtherMods;
};



class Matrix
{
// convention: mat[smpl][feature]
public:
	//Matrix(const string inF);
	Matrix(const string inF, const string, const string xtra, vector<string>& outFName, bool highLvl = false, bool NumericRowId = false, bool writeTmpFiles = true);

	Matrix(const string inF, const string xtra, bool highLvl = false); // this reads to mem
	Matrix(const vector<string>& rnms, const vector<string>& cnms);//module abundance matrix
	Matrix(void);
	~Matrix(void);
	void addTtlSmpl(vector<mat_fl> x, int idx) { mat[idx] = x; }
	void splitOnHDD(string out_seed);
	void writeSums(string);
	void normalize();
	void transpose();
	void writeMatrix(const string ofile,bool onlyFilled=false);
	size_t smplNum(){ return colIDs.size(); }
	int rowNum(){ return rowIDs.size(); }

	smplVec* getSampleVec(uint which){ return new smplVec(mat[which],1); }
	string getSampleName(uint which){ return colIDs[which]; }

	int SmplNum() { return (int)mat.size(); }
	int FtNum() {
		if (mat.size() >= 1) { return (int)mat[0].size(); }
		else { return 0; }
	}
	void estimateModuleAbund(char ** args, int argc);
	void estimateModuleAbund(options*);
	void resizeMatRows(uint x,mat_fl def=(mat_fl)0);
	//for the R module, all used for rarefactions only
	void addRow(vector<mat_fl>);//idea is that a single row is added on to the matrix
	void setSampleNames(vector<string> in) { colIDs = in; }
	void setRowNames(vector<string> in) { rowIDs = in; }
	vector < string > getSampleNames(){ return(colIDs); }
	vector < string > getRowNames(){ return(rowIDs); }
	//void addCount(string, int, mat_fl);

	double getMinColSum();
	column getMinColumn(uint offset = 0);
	vector< pair <double, string>> getColSums(bool sorted = false);
	void writeColSums(string outF);
protected:
	//subroutines
	void read_subset_genes(const string);
	void read_hierachy(const string );
	void addColumn(string);
	void readModuleFile(const string&);
	vector<mat_fl> getRowSums();
	void ini_mat();

	//storage
	vector< vector< mat_fl > > mat;
	vector< string > rowIDs,colIDs;
	unordered_map<string, int> colID_hash, rowID_hash;
	int maxCols;//number samples
	vector<HMat*> HI;
	LvlUp LUp;
	int maxLvl;
	vector<string> LvlNms;
	string sampleNameSep;
	GeneIDidx subset;
	bool doSubsets, doHigh;
	vector<double> colSum;

	vector< pair <double, string>> colsums;
};




class SigMatrix :public Matrix {
public:
	SigMatrix(const string& inf) { Matrix(inf, ""); }
private:
	void estimateBinModel() { 
#ifdef notRpackage
cerr << "todo estimateBinModel"; exit(33); 
#endif
}
	//functions to determine parameters
	//stored model parameters

};
class Modules : public Matrix
{
public:
	Modules(const string&, vector<string>);
	~Modules() {}

	void addDescription(const string&);
	void addHierachy(const string&);
	void writeModDescr(const string&,bool onlyUsed);


	void setRedund(int x) { redund = x; }
	void setPathwCompl(float x) { PathwCompl = x; }
	void setEnzymCompl(float x) { enzymCompl = x; }

	void calcModAbund(vector<mat_fl>&, const int pos, const unordered_map<string, int>&,
		vector<string>&, vector<float>&);


	vector<string> & modNms() { return rowIDs; }
	vector<string> modNms_numbered();//puts a number behind double modules
	vector<string> & modDescr() { return moduleDescriptions; }
	vector<string> & getRedundantUsedMods() { return redundantUsedMods; }
//special implementation to collapse rows..
	void writeMatrix(const string ofile, bool onlyFilled = false, bool collapseDblFeats = false);


private:
	void calc_redund();
	//contains the modules in the DB, each entry being one module
	vector<Module> mods;
	vector<string>  moduleDescriptions, redundantUsedMods; 
	vector<vector<string>> hierachy;
	//moduleNames = rowIDs
	//list of KOs used in DB, and how often they occur
	ModOccur MO;
	//in case of double entries, track these
	unordered_map<string, vector<int>> ModPos;
	vector<int> recurrentMods;
	vector<bool> ModUsed;

	//list of options
	int redund; // max redundancy of KOs used
	float PathwCompl; //corresponds to -c
	float enzymCompl; //single enzymes complexes - how much needs to be present to trigger present

};
