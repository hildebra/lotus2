#pragma once
#include "libload.h"
#include "options.h"
void trim(string& str,const std::string& whitespace = " \t");
bool isGZfile(const string fi);


struct TaxObj
{
	TaxObj(const string&, int, bool nativeSLV, bool doNotCheckTax);
	TaxObj(TaxObj*t);
	TaxObj(int d) :SavedTaxs(d, __unkwnTax), Subj(""), perID(0.f),depth(d) {}
	//functions
	string getWriteString();
	void copy_vals(TaxObj*t) { SavedTaxs = t->SavedTaxs; depth = t->depth; }
	//get tax at depth x
	string& get(int x) { if (x > depth) {  return __unkwnTax; } return SavedTaxs[x]; }
	void set(int x, string v) { 
		if (x < (int)SavedTaxs.size()) { SavedTaxs[x] = v; } 
		if (x > depth) { depth = x; }
	}
	//check if other tax is better and copies if so these vals over itself
	bool evalAcpyTax(TaxObj* oth);
	inline void copyOver(TaxObj* oth);
	void addHitDB(string x) { SavedTaxs.push_back(x); depth = SavedTaxs.size(); }

	//int dept() { return depth; }
	//variables
	vector<string> SavedTaxs;
	string Subj;
	float perID;
	int depth;//saves explicitly the depth, taking ? etc into account

};


class RefTax
{
public:
	RefTax(const string&,int tdep,bool,bool);
	~RefTax();
	void stats();
	int depth() { return (int) tlevels.size(); }
	unordered_map <string, TaxObj*>::const_iterator find(const string& s) {
		return Tlink.find(s);
	}
	unordered_map <string, TaxObj*>::const_iterator end() {
		return Tlink.end();
	}
	void setTaxLvls(vector<string> x) { tlevels = x; }
	//const string getLevels();
	//vector<string> taxLevls() { return tlevels; }
private:
	string TaxFile;
	unordered_map <string, TaxObj*> Tlink;
	vector<string> tlevels;
};

struct BlastRes
{
	BlastRes(const string&,int);
	bool isSameQuery(const string &q) {	if (q == Query) { return true; } return false;	}
	
	string Query, Sbj;
	int alLen;
	double perID, eval,score;
	bool fail;
};

class BlastReader
{
public:
	BlastReader(const string&, const string&);
	~BlastReader() {  }
	list<BlastRes*> getResBatch();
private:

	bool openedGZ,processedBatch;
	BlastRes* lastBlast;
	istream* blast;
	bool allRead;
	int inptFmt;
};
