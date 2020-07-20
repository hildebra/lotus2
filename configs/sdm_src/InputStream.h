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


#ifndef _InputStr_h
#define _InputStr_h


#include "DNAconsts.h"
#include <functional> 
#include <cctype>
#include <locale>

extern char DNA_trans[256];
extern short DNA_amb[256];
extern short NT_POS[256];
extern short DNA_IUPAC[256 * 256];
typedef float matrixUnit;


string spaceX(uint k);
int digitsInt(int x);
int digitsFlt(float x);
string intwithcommas(int value);
std::string itos(int number);
std::string ftos(float number);
bool isGZfile(const string fileS);//test if file is gzipped input

static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}


//MOCAT header fix
std::vector<std::string> header_string_split(const std::string str, const std::string sep);
void remove_paired_info(string&, short = -1);
//MOCAT header fix
std::string header_stem(string& header);
std::istream& safeGetline(std::istream& is, std::string& t);
string reverseTS2(const std::string & Seq);
void reverseTS(std::string & Seq);


bool any_lowered(const string& is);
//this function changes input string (file location) to have consistent file names
string applyFileIT(string x, int it, const string xtr = "");
bool fileExists(const std::string& name, int i=-1,bool extiffail=true);
//vector<int> orderOfVec(vector<int>&);




class ofbufstream {
public:
	ofbufstream(const string IF, int mif) :file(IF), modeIO(mif), used(0) {
		if (modeIO == ios::out) {
			remove(file.c_str());
		}
		keeper = new char[bufS];
	}
	~ofbufstream() {
		writeStream();
		delete[] keeper;
	}
	void operator<< (const string& X) {
		size_t lX(X.length());
		if (lX + used > bufS) {
			writeStream();
		}
		memcpy(keeper + used, X.c_str(), lX);
		used += lX;
	}
private:
	void writeStream() {
		if (used == 0) { return; }
		ofstream of(file.c_str(), ios::app);
		of.write(keeper, used);
		of.close();
		used = 0;
	}
	string file;
	char *keeper;
	int modeIO;
	size_t used;
	static const size_t bufS = 500000;
};



inline vector<string> splitByComma(const string& fileS,bool requireTwo, char SrchStr=','){
	string::size_type pos = fileS.find(SrchStr);
	if (pos == string::npos){
		if (requireTwo){
			cerr<<fileS<<endl;
			cerr << "Could not find \"" << SrchStr<<"\" in input file (required for paired end sequences, as these are two files separated by \",\")\n";
			exit(14);
		}else {
			return vector<string> (1,fileS);
		}
	}
	vector<string> tfas(2,"");
	tfas[0] = fileS.substr(0,pos);
	tfas[1] = fileS.substr(pos+1);
	return tfas;
} 

inline vector<string> splitByCommas(const string& fileS, char SrchStr = ',') {
	if (fileS.find(SrchStr) == string::npos) { return vector<string>(1, fileS); }
	vector<string> res = splitByComma(fileS, true, SrchStr);
	vector<string> ret(0); ret.push_back(res[0]);
	while (res[1].find(SrchStr) != string::npos) {
		res = splitByComma(res[1], true, SrchStr);
		ret.push_back(res[0]);
	}
	ret.push_back(res[1]);
	return ret;
}


//requires sorted vector with the entries being actual datapoints
template<class TYPE>
TYPE calc_median2(vector<TYPE>& in, float perc){
	size_t sum = in.size();
	size_t tar = (size_t)(((float)sum) * perc);
	return in[tar];
}


//returns "i_fna" or "i_fastq"
string detectSeqFmt(const string);


class Filters;

class DNA{
public:
	DNA(string seq, string names) :Seq(seq), SeqLength(Seq.length()),
		ID(names), NewID(names),
		Qual(0),QualTraf(""),Sample(-1),avgQual(-1.f),
		Qsum(0),AccumError(0.),goodQual(false),midQual(false),
		Read_position(-1), 
		FtsDetected(),
		IDfixed(false), tempFloat(0.f){}
	DNA():Seq(""),SeqLength(0),ID(""),NewID(""),Qual(0),QualTraf(""),
		Sample(-1), avgQual(-1.f),
		Qsum(0), AccumError(0.), goodQual(false), midQual(false),
		Read_position(-1),
		FtsDetected(),
		IDfixed(false), tempFloat(0.f) {
	}
	bool operator==(DNA i) {
		if (i.getSeqPseudo() == this->getSeqPseudo()) {
			return true;
		}	else {
			return false;
		}
	}
	bool operator==(shared_ptr<DNA> i) {
		if (i->getSeqPseudo() == this->getSeqPseudo()) {
			return true;
		} else {
			return false;
		}
	}
	
	//~DNA(){}
	void append(const string &s) { Seq += s; SeqLength = Seq.length(); }
	void Qappend(const vector<qual_score> &q);
	void setSeq(string & s) { Seq = s; SeqLength = Seq.length();/*DN = Seq.c_str();*/ }
	string &getSeq() { return Seq; }
	string getSeq_c() { return Seq; }
	string getSeqPseudo() { return Seq.substr(0, SeqLength); }
	void setQual(vector<qual_score>& Q) { Qual = Q; avgQual = -1.f; }
	const string& getID() { if (IDfixed) { return NewID; }return ID; }
	string getID_copy() { string x = getID(); string y = x; return y; }
	string getIDPosFree(); // remove /1 /2 #1:0 etc
	const string& getOldID() { return ID; }
	string getIDshort() { return ID.substr(0, getShorterHeadPos(ID)); }
	string getNewIDshort() { return NewID.substr(0, getShorterHeadPos(NewID)); }
	bool seal();
	bool isEmpty() { if (ID.length() == 0 && Seq.length() == 0) { this->setPassed(false); return true; }  return false; }

	int numACGT();
	float getAvgQual();
	unsigned int getQsum(){return Qsum;}
	float qualWinfloat(unsigned int,float,int&);
	
	
	float binomialFilter(int, float);
	//float qualWinfloat_hybr(int,float,int,float,int&);
	bool qualWinPos(unsigned int,float);
	bool qualAccumTrim(double d);
	int qualAccumulate(double d);
	double getAccumError(){
		if (AccumError == 0.f) { for (uint i = 0; i < Qual.size(); i++) { if (Qual[i] >= 0) { AccumError += SAqualP[Qual[i]]; } } }
		if (std::isinf((double)AccumError)) {
			AccumError = 5.f;
		}
		return AccumError;}
	int minQual(){int mq=50; for (uint i=0;i<Qual.size();i++){if (Qual[i]<mq){mq=Qual[i];}}return mq;}
	void NTspecQualScores(vector<long>&, vector<long>&);


	inline uint length() { return (uint)SeqLength; }
	inline uint mem_length() { return (uint) Seq.length(); }
	bool cutSeq(int start, int stop=-1, bool = false);
	bool HomoNTRuns(int);
	int matchSeq(string, int, int, int);
	void reverse_transcribe();
	int matchSeqRev(string, int, int, int=0);
	int matchSeq_tot(string, int, int, int&);
	void writeSeq(ostream&,bool singleLine=false);
	void writeQual(ostream&, bool singleLine = false);
	void writeFastQ(ostream&,bool=true);
	void writeFastQ(ofbufstream&, bool = true);
	void writeFastQEmpty(ostream&);
	void setNewID(string x) { NewID = x; }
	void newHEad(string x){NewID=x;ID=x;}
	void changeHeadPver(int ver);
	void setTA_cut(bool x) { FtsDetected.TA_cut = x; }
	bool getTA_cut() { return FtsDetected.TA_cut; }
	void setBarcodeCut() { FtsDetected.Barcode_cut = true; FtsDetected.Barcode_detected = true; }
	bool getBarcodeCut() { return FtsDetected.Barcode_cut; }
	void setBarcodeDetected(bool x){ FtsDetected.Barcode_detected = x; }
	bool getBarcodeDetected() { return FtsDetected.Barcode_detected; }
	bool isMIDseq() { if (Read_position == 3) { return true; } return false; }
	void setMIDseq(bool b){ if (b){ Read_position = 3; } }
	void setpairFWD(){ Read_position = 0; }
	void setpairREV(){ Read_position = 1; }
	int getReadMatePos() { return (int) Read_position; }
	bool sameHead(shared_ptr<DNA>);
	bool sameHead(const string&);
	//inline void reverseTranscribe();
	void setTempFloat(float i){tempFloat = i;}
	float getTempFloat(){return tempFloat;}
	void adaptHead(shared_ptr<DNA>,const int,const int);
	void failed(){goodQual=false;midQual=false;}
	bool control(){ if (Qual.size()==0){return false;}return true;}
	void setBCnumber(int i, int BCoff) { if (i < 0) { Sample = i ; FtsDetected.Barcode_detected = false; } else { Sample = i + BCoff; FtsDetected.Barcode_detected = true; } }
	int getBCnumber();//always return BC tag IDX global (no local filter idx accounted for, use getBCoffset() to correct)

	void prepareWrite(int fastQver);
	void reset();
	void resetTruncation() { SeqLength = Seq.length(); }
	void setPassed(bool b);
	void setMidQual(bool b) { midQual = b; }
	bool isPassed(void){return goodQual;}
	bool isMidQual(void){return midQual;}
	string getSubSeq(int sta, int sto){return Seq.substr(sta,sto);}
	void resetQualOffset(int off, bool solexaFmt);
	
	//control & check what happened to any primers (if)
	bool has2PrimersDetected() { return (FtsDetected.reverse && FtsDetected.forward); }
	bool getRevPrimCut() { return FtsDetected.reverse; }
	bool getFwdPrimCut() { return FtsDetected.forward; }
	void setRevPrimCut() { FtsDetected.reverse = true; }
	void setFwdPrimCut() { FtsDetected.forward = true; }
	//only used in pre best seed step
	//float getSeedScore() { return tempFloat; }
	//void setSeedScore(float i) { tempFloat = (float)i; }

	struct QualStats {
		bool maxL; bool PrimerFail; bool AvgQual; //sAvgQual
			bool HomoNT; bool PrimerRevFail; bool minL; 
			bool minLqualTrim; //<-sMinQTrim trimmed due to quality
			bool TagFail; bool MaxAmb; bool QualWin;//sQualWin 
			bool AccErrTrimmed; bool QWinTrimmed;  // either of these makes bool Trimmed; 
			bool fail_correct_BC; bool suc_correct_BC; bool
			failedDNAread; 
			//bool adapterRem; -> setTA_cut
			bool RevPrimFound; 
			bool BinomialErr; bool dblTagFail;
		QualStats() :
			maxL(false), PrimerFail(false), AvgQual(false), HomoNT(false),
			PrimerRevFail(true), minL(false), minLqualTrim(false),
			TagFail(false), MaxAmb(false), QualWin(false),
			AccErrTrimmed(false), QWinTrimmed(false),
			fail_correct_BC(false), suc_correct_BC(false),
			failedDNAread(false), RevPrimFound(false),
			BinomialErr(false),
			dblTagFail(false)
		{}
	} QualCtrl;

protected:
	size_t getShorterHeadPos(const string & x, int fastQheadVer=-1) {

		size_t pos(string::npos);
		if (fastQheadVer != 0) {
			if (Read_position == 1) {
				pos = x.find("/2");
				//			if (pos == string::npos) {	pos = x.find_first_of(" 1:");}
					//		if (pos == string::npos) { pos = x.find_first_of("/1"); }
			}
			else if (Read_position == 0) {
				pos = x.find("/1");
			}
			else {
				pos = x.find("/1");
				if (pos == string::npos) { pos = x.find("/2"); }
			}
		}
		//if (pos == string::npos){pos=x.length()-min((size_t)5,x.length());}}
		//if(pos<0){pos=0;}
		if (pos == string::npos) { pos = min(x.find(' '), x.find('\t')); }
		
		if (pos == string::npos) { pos = x.length(); }
		return pos;
	}
	//mainly used to mark if rev/Fwd primer was detected
	string xtraHdStr();
	size_t getSpaceHeadPos(const string & x) {
		size_t pos = x.find(' ');
		if (pos == string::npos) { pos = x.length(); }
		return pos;
	}
	//binomial accumulated error calc
	inline float interpolate(int errors1, float prob1, int errors2, float prob2, float alpha);
	float sum_of_binomials(const float j, int k, float n, int qual_length, const vector<float>& error_probs, const vector< vector<float>> & per_position_accum_probs);
	inline float prob_j_errors(float p, float j, float n);
	
	inline bool matchDNA(char,char);
	
	string Seq;
	size_t SeqLength;
	string ID,NewID; //original and newly constructed ID
	vector<qual_score> Qual;
	string QualTraf;
	int Sample;

	//const char* DN;
	float avgQual;
	unsigned int Qsum;
	double AccumError;
	bool goodQual,midQual;
	//bool TA_cut, Barcode_cut; //technical adapter, barcode (tag)	
	short Read_position;//-1=unkown; 0=pair1 (fwd primer); 1=pair2 (rev primer); 3=MID seq ;

	struct ElementsDetection{
		bool forward; bool reverse;//primers detected
		bool TA_cut; bool Barcode_detected;  bool Barcode_cut;
		ElementsDetection() :forward(false), reverse(false), TA_cut(false), Barcode_detected(false), Barcode_cut(false) {}
		void reset() { forward = false; reverse = false; TA_cut = false; Barcode_detected = false; Barcode_cut = false; }
	} FtsDetected;



	bool IDfixed;
	float tempFloat;
};

struct DNAHasher
{
	size_t operator()(shared_ptr<DNA> k) const
	{
		// Compute individual hash values for two data members and combine them using XOR and bit shifting
		return ((hash<string>()(k->getSeqPseudo())) >> 1);
	}
};





class DNAunique : public DNA{//used for dereplication
public:
	DNAunique() : DNA(), Count(0), pair(0){}//chimeraCnt((matrixUnit) 1), 
	DNAunique(string s, string x) :DNA(s, x), Count(1) {}
	DNAunique(shared_ptr<DNA>d, int BC) : DNA(*d), Count(0), BestSeedLength( (uint)Seq.size()),pair(0){ addSmpl(BC);  }
	~DNAunique() { ; }// if (pair != NULL) { delete pair; }
	//string Seq;	string ID;
	void Count2Head(bool);
	void addSmpl(int k);
	void writeMap(ofstream & o, const string&, vector<int>&, const vector<int>&);
	inline int getCount() { return Count; }
	uint getBestSeedLength() { return BestSeedLength; }
	void setBestSeedLength(uint i) { BestSeedLength = i; }
	void setOccurence(int smpl, int N);
	void transferOccurence(shared_ptr<DNAunique>);
	const unordered_map<int, int> & getDerepMap() { return occurence; }
	vector<int> getDerepMapSort(size_t);
	//vector<pair<int, int>> getDerepMapSort2(size_t wh);
	void getDerepMapSort(vector<int>&, vector<int>&);
	void saveMem() { QualTraf = ""; NewID = ID.substr(0, getSpaceHeadPos(ID)); ID = ""; }
	void attachPair(shared_ptr<DNAunique> d) { pair = d; pair->saveMem(); }
	shared_ptr<DNAunique> getPair(void) { return pair; }
	//estimates if one sample occurence covers the unique counts required for sample specific derep min counts
	bool pass_deprep_smplSpc( const vector<int>&);

	//matrixUnit chimeraSplitNum() { return chimeraCnt; }
	//void setChimSplitNum(matrixUnit x) { chimeraCnt = x; }
	//sort
	//bool operator < (const DNAunique& str) const {		return (Count < str.Count);	} 
private:
	int Count;
	//matrixUnit chimeraCnt;
	int BestSeedLength;
	unordered_map<int, int> occurence;
	shared_ptr<DNAunique> pair;

};

class InputStreamer{
public:
	InputStreamer(bool fnRd, int fq, string ignoreInptEr="0") :
		_fileLength(10), _max(60), _last(0),
		fna_u(3, NULL), qual_u(3, NULL), fastq_u(3,NULL),
		inFiles_fna(3, ""), inFiles_qual(3, ""), inFiles_fq(3, ""),
		//fna(3,NULL), qual(3,NULL), fastq(3,NULL),
		
#ifdef _gzipread	
		//gzfna(3,NULL), gzqual(3,NULL), gzfastq(3,NULL),
#endif
		tdn1(3, NULL), tdn2(3, NULL),
		fnaRead(fnRd), hasMIDs(false),
		lnCnt(3, 0), fastQver(fq),
		minQScore(1000), maxQScore(-1), 
		QverSet(true), numPairs(1),
		pairs_read(3, 0), opos(3,0), 
		currentFile(0), totalFiles(0), BCnumber(0),
		qualAbsent(false),
		fqReadSafe(true), fqPassedFQsdt(true),
		fqSolexaFmt(false), openedGZ(false),
		ErrorLog(0), DieOnError(true)
	{
		opos[0] = 1; if (fastQver == 0) { QverSet = false; }
		if (ignoreInptEr=="1") { DieOnError = false; }
	}
	~InputStreamer();
	//path, fasta, qual, pairNum
	string setupInput(string path, int i, int tarID,
		const vector<string>& uF, const vector<string>& FQ, const vector<string>& Fas, const vector<string>& Qual,
		const vector<string>& midf, int &paired, string onlyPair, 
		string& shortMainFile, bool simu = false);
	bool setupFastaQual(string,string, string, int&, string,bool=false);
	void setupFna(string);
	//path, fastq, fastqVer, pairNum
	bool setupFastq(string,string, int&,string,bool = false);
	//0=pair 1; 1=pair 2; 2=midSeq; sync=synchronize read pairs (ie only first pair read so far, jump to same DNA reads with second pair)
	shared_ptr<DNA> getDNA(bool&,int,bool& sync);
	void jumpToNextDNA(bool&, int);
	//shared_ptr<DNA> getDNA2(bool&);
	//shared_ptr<DNA> getDNA_MID(bool&);
	bool hasMIDseqs(){return hasMIDs;}
	void allStreamClose();
	void allStreamReset();
	void openMIDseqs(string,string);
	int pairNum() { return numPairs; }
	bool qualityPresent() { return !qualAbsent; }
	bool checkInFileStatus();
	void atFileYofX(uint cF, uint tF, uint BCn) { currentFile = cF; totalFiles = tF; BCnumber = BCn; }
	uint getCurFileN() { return currentFile; }

private:
	inline qual_score minmaxQscore(qual_score t);// , int lnCnt);
	int parseInt(const char** p1);// , int &pos);// , const char ** &curPos);
	bool setupFastq_2(string, string, string);
	bool setupFastaQual2(string, string, string = "fasta file");
	shared_ptr<DNA> read_fastq_entry(istream & fna, int &minQScore,
		int&,bool&,bool);
	shared_ptr<DNA> read_fastq_entry_fast(istream & fna, int&,bool&);
	void jmp_fastq(istream &, int&);
	bool read_fasta_entry(istream&fna,istream&qual,shared_ptr<DNA> in,shared_ptr<DNA>,int&);
	bool getFastaQualLine(istream&fna,  string&);
	void maxminQualWarns_fq();
	int auto_fq_version();
	int auto_fq_version(int minQScore, int maxQScore=0);
	void resetLineCounts(){ lnCnt[0] = 0;	lnCnt[1] = 0;	lnCnt[2] = 0; }
	bool desync(int pos) { if ( abs(pairs_read[pos] - pairs_read[opos[pos]]) > 1 ) {return true; } return false; }
	void IO_Error(string x);
	//bar on file read progress
	void _measure(istream &);
	inline bool _drawbar(istream &);
	inline void _print(int cur, float prog);
	int _fileLength, _max, _last;



	//abstraction to real file type
	vector<istream*> fna_u, qual_u, fastq_u;
	vector<string> inFiles_fna, inFiles_qual, inFiles_fq;
	//0,1,2 refers to pairs / MID fasta files
	//vector<ifstream> fna, qual, fastq;
	//ifstream qual, fastq,
		//second pair
		//ifstream fna2, qual2, fastq2,
		//usually used for MID
		//fna3, qual3, fastq3;

	//required for Fasta in term storage
	vector<shared_ptr<DNA>> tdn1; vector<shared_ptr<DNA>> tdn2;
	//shared_ptr<DNA> tdn21; shared_ptr<DNA> tdn22;
	//shared_ptr<DNA> tdn31; shared_ptr<DNA> tdn32;
	bool fnaRead, hasMIDs;
	vector<int> lnCnt;// , lnCnt2, lnCnt3;//line count
	int fastQver,minQScore,maxQScore;//which version of Fastq? minima encountered Qscore
	bool QverSet;
	//1 or 2?
	int numPairs;
	//keep track of sequences read for each pair; other position (1=2,2=1)
	vector<int> pairs_read, opos;
	//some stats to print, nothing really relevant
	uint currentFile, totalFiles, BCnumber;
	//is quality information even available?
	bool qualAbsent;
	//fq format not checked for completeness
	bool fqReadSafe, fqPassedFQsdt, fqSolexaFmt;
	bool openedGZ;

	//collects errors, handles errors
	vector<string> ErrorLog;
	bool DieOnError;
};






#ifdef _gzipread2
std::vector< char > readline(gzFile f);
#endif 

#endif