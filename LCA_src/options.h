#pragma once
#include "libload.h"


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);


struct options
{
public:
	options(int argc, char **argv, int);
	~options();
	const string TaxLvl2string();


	//vars
	string RefTaxFile, blastres, outF;
	string input_format;//bl8 , uc
	bool BLfilter;     //do my own blast filter before LCA
	bool calcHighMats; //calculate phylum etc level sum
	bool hitRD; //show the database entry that was hit, in hiera file
	bool isReads;
	bool annotateAll; //give out OTU / READ even if not assigned??
	bool nativeSlVdb;
	bool reportBestHit; //reports best hit, if higher than required %id
	bool checkTaxoUnkw; //check in the tax DB, if unkownn, ? etc levels are there and replaces with ?
	int numThr; // number of threads
	int taxDepth; //how deep does the taxonomy go?
	double LCAfract;
	vector<double> idThr;
	vector<string> blFiles, refDBs;
	vector<string> Taxlvls;
};

