/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand
email: Falk.Hildebrand@gmail.com

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

#include "InputStream.h"
string spaceX(uint k){
	string ret = "";
	for (uint i = 0; i < k; i++){
		ret += " ";
	}
	return ret;
}

int digitsInt(int x){
	int length = 1;
	while (x /= 10)
		length++;
	return length;
}
int digitsFlt(float x){ 
	std::stringstream s;
	s << x;
	return (int)s.str().length(); 
}
string intwithcommas(int value) {
	string numWithCommas = std::to_string((long long)value);
	int insertPosition = (int)numWithCommas.length() - 3;
	while (insertPosition > 0) {
		numWithCommas.insert(insertPosition, ",");
		insertPosition -= 3;
	}
	return (numWithCommas);
}

std::string itos(int number) {
	std::stringstream ss;
	ss << number;
	return ss.str();
}
std::string ftos(float number) {
	std::stringstream ss;
	ss << number;
	return ss.str();
}
bool isGZfile(const string fi) {
	string subst = fi.substr(fi.length() - 3);
	if (subst == ".gz") {
		return true;
	}
	return false;
}

std::istream& safeGetline(std::istream& is, std::string& t) {
	t.clear();
	//from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();

	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:
			// Also handle the case when the last line has no line ending
			if (t.empty())
				is.setstate(std::ios::eofbit);
			return is;
		default:
			t += (char)c;
		}
	}
}
//MOCAT
std::vector<std::string> header_string_split(const std::string str, const std::string sep) {
	std::vector<std::string> tokens;
	tokens.reserve(13);
	size_t start = 0;
	size_t pos = 0;
	while ((pos = str.find_first_of(sep, start)) != std::string::npos) {
		tokens.push_back(str.substr(start, pos - start));
		start = pos + 1;
	}
	if (start < str.length()) {
		tokens.push_back(str.substr(start));
	} else if (start == str.length()) {
		tokens.push_back("0");
	}
	return tokens;
}
void remove_paired_info(string &s, short RP) {
	size_t f1 = s.find(" ");
	if (f1 != string::npos) {
		s = s.substr(0, f1);
		f1 = string::npos;
	}
	if (RP == 0) {
		f1 = s.find("/1");
	}
	else if (RP == 1) {
		f1 = s.find("/2");
	}
	else {
		f1 = s.find("/1");
		if (f1 == string::npos) { f1 = s.find("/2"); }
	}
	if (f1 != string::npos && f1 == s.size() - 2) {
		s = s.substr(0, f1);
	}
}
std::string header_stem(string& header) {
	const size_t slash = header.find('/');
	if (slash != std::string::npos) {
		return header.substr(0, slash);
	}
	
	std::vector<std::string> tokens = header_string_split(header, ": ");

	if (tokens.size() == 11) {
		//if (tokens[8] == "Y") seq.clear();
		return tokens[0] + ":" + tokens[3] + ":" + tokens[4]
			+ ":" + tokens[5] + ":" + tokens[6];// +"#" + tokens[10];
	} else if (tokens.size() == 6 || tokens.size() == 7) {
		return tokens[0] + ":" + tokens[2] + ":" + tokens[3]
			+ ":" + tokens[4] + ":" + tokens[5];// + "#0";
	} else {
		cerr << "fastq header format error\n"; exit(98);
	}
		//throw std::runtime_error("fastq header format error\n" + header); }
}
void reverseTS(string & Seq) {
	int qs = (int)Seq.length() - 1;
	string S2 = Seq.c_str();
	for (int i = qs; i >= 0; i--) {
		Seq[i] = DNA_trans[(int)S2[qs - i]];
	}
}
string reverseTS2(const string & Seq) {
	int qs = (int)Seq.length() - 1;
	string S2 = Seq.c_str();
	for (int i = qs; i >= 0; i--) {
		S2[i] = DNA_trans[(int)Seq[qs - i]];
	}
	return S2;
}
bool any_lowered(const string& is) {
	for (uint i = 0; i < is.length(); i++) {
		char c = is[i];
		if (islower(c)) { return true; }
	}
	return false;
}
string applyFileIT(string x, int it,const string xtr){

	size_t pos = x.find_last_of(".");
	if (pos != string::npos && isGZfile(x)) {
		pos = x.find_last_of(".", pos-1);
	}
	if (it == 0) {
		if (pos == string::npos) {
			return x + xtr;
		}
		return x.substr(0, pos) + xtr + x.substr(pos);
	} 
	ostringstream ss;
	ss << it;
	if (pos == string::npos) {
		return x + "." + ss.str() + xtr ;
	}
	return x.substr(0,pos) + "."+ss.str() + xtr + x.substr(pos);
	
	
}
bool fileExists(const std::string& name, int i, bool extiffail) {
	if (name == "") { return true; }
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);		return true;
	} else {
		if (extiffail) {
			cerr << "ERROR: Could not find file " << name << endl;
			if (i >= 0) {
				cerr << "on mapping file line " << i << endl;
			}
			exit(92);
		}
		return false; 
	}
}

string detectSeqFmt(const string inF) {

	vector<string> tfasP = splitByCommas(inF, ';');
	for (size_t i = 0; i < tfasP.size(); i++) {
		vector<string> tfas = splitByCommas(tfasP[i]);
		string fileS = tfas[0];
		istream* fnax(NULL);
		string file_type = "test file";
		if (fileS != "") {
			if (isGZfile(fileS)) {
#ifdef _gzipread
				file_type = "gzipped fasta file";
				fnax = new igzstream(fileS.c_str(), ios::in);
#else
				cerr << "gzip not supported in your sdm build\n" << fileS; exit(50);
#endif
			}
			else {
				fnax = new ifstream(fileS.c_str(), ios::in);

			}
			if (!*(fnax)) { cerr << "\nCouldn't find " << file_type << " file \"" << fileS << "\"!\n Aborting..\n"; exit(4); }
			//char Buffer[RDBUFFER];
			//fna_u[0]->rdbuf()->pubsetbuf(Buffer, RDBUFFER);
		}
		string tmp("");
		string ret = "";
		while (safeGetline(*fnax, tmp)) {
			if (tmp[0] == '>') {
				ret = "-i_fna"; break;
			}
			else if (tmp[0] == '@') {
				ret = "-i_fastq"; break;
			}
			else if (tmp.length() == 0) {//do nothing
				;
			}
			else {
				cerr << " Could not auto detect input format. First non-empty line of your file looked like:\n" << tmp << endl;
				exit(888);
			}
		}
		if (ret == "") {
			cerr << "Empty input file detected:\n" << fileS << endl;
		} else {
			return ret;
		}

		delete fnax;
	}

	return "empty";
}


/*vector<int> orderOfVec(vector<int>& vin) {
	struct MyStruct
	{
		int key;
		int Value;
		MyStruct() :key(0), Value(0) {}
		MyStruct(int k, const int s) : key(k), Value(s) {}

		bool operator < (const MyStruct& str) const {
			return (key > str.key);
		}
	};
	std::vector < MyStruct > vec(vin.size());
	//fill vector
	for (int i = 0; i < (int)vin.size(); i++) {
		vec[i] = MyStruct(vin[i], i);
	}

	sort(vec.begin(), vec.end());
	//extract from sorted vector
	vector<int> ret(vin.size(), 0);
	for (size_t i = 0; i < vin.size(); i++) {
		ret[i] = vec[i].Value;
	}
	return ret;
	}
	*/
bool DNA::seal() {//DN = Seq.c_str();
	size_t QSi = Qual.size();
	if (QSi == 0 && Seq == "" ) {

		this->setPassed(false);
		return true;//nothing to be done, just empty DNA
	}
	if (QSi != Seq.length()) {
		cerr << "Unequal length of seq and quality for name " << this->getID() << "\n";
		this->setPassed(false);
		return false;
	}
	//uppercase DNA
	std::transform(Seq.begin(), Seq.end(), Seq.begin(), ::toupper);
	SeqLength = Seq.length();
	return true;
}

string DNA::getIDPosFree() { // remove /1 /2 #1:0 etc
	string s = this->getIDshort();
	remove_paired_info(s,Read_position);
	return s;
} 


int DNA::numACGT(){
	int DNAch = 0;
	for (unsigned int i = 0; i < length(); i++){

		DNAch += DNA_amb[(int)Seq[i]];
//		if (tmp == 'A' || tmp == 'C' || tmp == 'G' || tmp == 'T'){			DNAch++;		}
	}
	return static_cast<int> (length()) - DNAch;
}
void DNA::Qappend(const vector<int> &q)
{
	Qual.insert(Qual.end(), q.begin(), q.end());
	/*unsigned int Qsiz = (unsigned int) Qual.size();
	Qual.resize(Qsiz+q.size(),0);
	for (register unsigned int i=0; i<q.size();i++){
	Qual[i+Qsiz] = q[i];
	}
	*/
	avgQual = -1.f;
}

float DNA::getAvgQual(){
	if (avgQual == -1.f){
		if (Qsum == 0){
			for (unsigned int i = 0; i < this->length(); i++){
				Qsum += Qual[i];
			}
		}
		avgQual = static_cast<float> (Qsum) / static_cast<float> (this->length());

	}
	return avgQual;
}

/*float DNA::qualWinfloat_hybr(int W, float T, int W2, float T2, int& reason){//not used
	//if (T==0.f){return true;}
	int AQS=0, AQL;
	int TotQ = 0;
	int upTs = static_cast<int>(T * W);
	int upTl = static_cast<int>(T2 * W2);
	int QS = int (Qual.size());
	if (W>=QS){W = QS;} // too short
	int smallerW = W, largerW = W2;

	bool W1IsSmall = true;
	if (smallerW > W2){ largerW=W; smallerW = W2; W1IsSmall=false;
	std::swap(upTl,upTs); }

	for (unsigned int i=0; i<(unsigned int) smallerW; i++){
	AQS += Qual[i];
	}
	AQL = AQS;

	//hybrid schleife
	for (unsigned int i=smallerW; i<(unsigned int) largerW; i++){
	AQL += Qual[i];
	AQS += Qual[i]; AQS -=  Qual[i-smallerW];
	if (AQS < upTs){
	if (W1IsSmall){		reason=0;	return 0.f;
	} else {reason=1; this->cutSeq(i/2,this->length()); return float(AQL)/float(QS);}
	}
	}

	TotQ = AQL;
	for (unsigned int i=W; i<(unsigned int) QS; i++){
	AQS += Qual[i]; AQL += Qual[i];
	TotQ += Qual[i];
	AQS -= Qual[i-W];AQL -= Qual[i-W];
	if (AQS < upTs || AQL < upTl){
	if (W1IsSmall){		reason=0;	return 0.f;
	} else {reason=1; this->cutSeq(i/2,this->length()); return float(AQL)/float(QS);}
	}
	}
	//if (averageQ > static_cast<float> (TotQ) /static_cast<float> ( Qual.size())){return false;}
	return static_cast<float> (TotQ) /static_cast<float> ( QS);
	}
	*/

// modified from https://github.com/fpusan/moira/blob/master/moira/bernoullimodule.c
float DNA::interpolate(int errors1, float prob1, int errors2, float prob2, float alpha)
{
	float result = errors1 + ((errors2 - errors1) * ((1 - alpha) - prob1) / (prob2 - prob1));
	if (result < 0) //Happens only for very-short high qual sequences in which the probability of having 0 errors is higher than 1 - alpha.
	{
		result = 0;
	}
	return result;
}
float DNA::prob_j_errors(float p, float j, float n) //Where p is the error probability, j is the number of errors and n the number of observations.
{
	float per_position_accum_probs;
	if (j > n)	{
		return 0.0f; //The formula below would also return 0.
	}
	per_position_accum_probs = pow((1 - p), n);	//For j == 0.
	float i(1);
	for (; i <= j; i += 1.f)	{//For j > 0.
		per_position_accum_probs = ((n - i + 1.f) / (1.0f*i)) * (p / (1.f - p)) * per_position_accum_probs;
	}
	return per_position_accum_probs;
	
}
float DNA::sum_of_binomials(const float j, int k, float n, int qual_length, const vector<float>& error_probs, const vector< vector<float>> & per_position_accum_probs)
//#Where j is the number of errors and k is the position in the sequence.
{
	float probability = 0;
	float i(0); int k1 = (int) k - 1;

	for (; i <= j; i+=1.f)
	{
		probability += DNA::prob_j_errors(error_probs[k], i, n) *  per_position_accum_probs[int(j - i)][k1];
		//Where error_probs[k] is the error probability of the k-th position.
		//Where per_position_accum_probs[j-i][k-1] is the probability that all the bases from position 0 to k-1 had a total of j-i errors.
	}

	return probability;
}

float DNA::binomialFilter(int maxErr, float alpha){

	if (alpha == -1.f|| this->length()<3){ return 0; }//deactivated

	///Initialize some variables.
	
	int SeqLengthI = (int)SeqLength;
	int n = 1; //Since we have a Bernoulli distribution.
	float n_f = (float)n;
	float alpha1 = 1.f - alpha;

	///Translate quality scores into error probabilities.
	vector<float> error_probs (SeqLength,0.f);
	for (size_t i = 0; i < SeqLength; i++){
		error_probs[i] = (float)SAqualP[Qual[i]];//pow(10, (contig_quals[i] / -10.0)); //Since we want a continuous list of non-N error_probs.
	}

	///Actually get the job done.
	int max_expected_errors = maxErr + 3;// (int)SeqLength + 1;
	int expected_errors = 0;
	float probability;
	vector<float> accumulated_probs (max_expected_errors,0.f);
	//int j;
	int k;
	vector<float> empty(SeqLength, 0.f);
	vector< vector<float>> per_position_accum_probs(max_expected_errors, empty);

	while (1)
	{
		float expected_errors_f = (float)expected_errors; 
		//vector<float > per_position_accum_probs(SeqLength, 0.f);
		for (k = 0; k < (int)SeqLength; k++) {
			if (k == 0)	{
				per_position_accum_probs[expected_errors][k] = DNA::prob_j_errors(error_probs[k], expected_errors_f, n_f);
			} else {
				
				per_position_accum_probs[expected_errors][k] = DNA::sum_of_binomials((float)expected_errors, k, n_f, SeqLengthI, error_probs, per_position_accum_probs);
			}
		}
		probability = per_position_accum_probs[expected_errors][SeqLengthI - 1];

		if (expected_errors == 0){
			accumulated_probs[expected_errors] = probability;
		}else{
			accumulated_probs[expected_errors] = accumulated_probs[expected_errors - 1] + probability;
		}

		if (accumulated_probs[expected_errors] > (alpha1) || expected_errors >= (max_expected_errors-1)){
			break;
		}else{
			expected_errors++;
		}
	}
	if (expected_errors == 0){
		return 0;
	}
	float EXE = interpolate(expected_errors - 1, accumulated_probs[expected_errors - 1], expected_errors, accumulated_probs[expected_errors], alpha);
	return EXE;
}

float DNA::qualWinfloat(unsigned int W, float T, int& reason){
	//if (T==0.f){return true;}
	int AQS=0;
	int TotQ = 0;
	int upTs = static_cast<int>(T * W);
	unsigned int QS = this->length();//static_cast<unsigned int> (Qual.size());
	if (W>=QS){W = QS;} // too short

//1st loop to ini window
	for (unsigned int i=0; i<(unsigned int) W; i++){
		AQS += Qual[i];
	}
	TotQ = AQS;

	for (unsigned int i=W; i<(unsigned int) QS; i++){
		AQS += Qual[i] - Qual[i - W];
		TotQ += Qual[i];
		if (AQS < upTs ){
			reason=1;	return 0.f;
		}
	}
	//if (averageQ > static_cast<float> (TotQ) /static_cast<float> ( Qual.size())){return false;}
	return static_cast<float> (TotQ) /static_cast<float> ( QS);
}

int DNA::qualAccumulate(double d){
	unsigned int i(0);double accErr(0.0);
	for (; i<this->length(); i++) {
		accErr+=SAqualP[Qual[i]];
		if (accErr >= d){break;}
	}

	this->AccumError = accErr;

	return i;
}
void DNA::NTspecQualScores(vector<long>& qsc, vector<long>& ntcnt) {
	size_t sql = Seq.length();
	if (qsc.size() < sql) {
	//	cerr << "qsc";
		qsc.resize(sql,0);
	}
	if (ntcnt.size() < sql) {
	//	cerr << "ntcnt";
		ntcnt.resize(sql,0);
	}
	for ( uint i = 0; i < sql; i++ ) {
		short p = NT_POS[(int)Seq[i]];
		qsc[p] += Qual[i];
		ntcnt[p]++;
	}
}

bool DNA::qualAccumTrim(double d){
	if (d == -1.) {
		return true;
	}
	unsigned int i(qualAccumulate(d));
	if (i != this->length()){
		//cut 3' end
		this->cutSeq(i,this->length());
		this->QualCtrl.AccErrTrimmed = true;
		return false;
	}
	//did not cut this sequence:
	return true;
}
bool DNA::qualWinPos(unsigned int W, float T){
	if (T==0.f){return true;}
	int AQ=0;
	int unT = static_cast<int>((float)W*T);
	unsigned int QS = this->length();// (unsigned int)Qual.size();
	unsigned int QSh = QS >> 2;
	QSh = max(QSh,W);
	if (W>=QS){return true;} // too short

	vector<float> WQ((int) W);
	//TODO: check that the right num of nucs is taken..
	int cnt=0;
	if (QS < W) {return true;}
	for (unsigned int i=QS-1; i> QS-(unsigned int) W-1; i--){
		AQ += Qual[i];
		cnt++;
	}
	if (AQ > unT){return true;}
	int curW = QS-(unsigned int) W;
	for (uint i=QS-(unsigned int) W-1; i > QSh; i--){
		AQ += Qual[i];
		AQ -= Qual[i+W];

		if (AQ < unT){ //min Window qual was broken.. kick seq
			curW = i;
		} else {
			break;
		}
	}
	
	//partial seq  removal
	int pos = curW - (W>>1);
	this->cutSeq(pos);
	this->QualCtrl.QWinTrimmed = true;
	return false;
}

//removes part of seq and qual indexes between start & stop
bool DNA::cutSeq(int start, int stop,bool Pseudo){

	if (stop == -1) {
		if (start >= (int) SeqLength) { return false; }
	} else if (start >= stop || stop > (int)Qual.size() || start >= (int)Qual.size()) {
		return false;
	}
	
	//pseudo deactivates cutting of 3'
	if (Pseudo) {
		if (stop == -1 || stop <= (int) SeqLength) {
			SeqLength = start;
			return true;
		}
	}

	string se = Seq;
	if (stop == -1) {
		stop = (int)Seq.length();
	}
	if (start == 0) {
		Seq = se.substr(stop);
	} else {
		Seq = se.substr(0, start) + se.substr(stop);
	}

	//DN = Seq.c_str();
	//Quali
	Qual.erase(Qual.begin()+start, Qual.begin()+stop); 
	
	SeqLength = Seq.length();

	return true;
}

int DNA::matchSeq(std::string PrSt,int Err,int tolerance, int startPos){
	//const char* DN = Seq.c_str();
	//const char* Pr = PrSt.c_str();
	int PrL = (int) PrSt.length();
	int mthL = this->length() - PrL;
	//int wantSc = PrL - Err;
	int endPos(-1),pos(startPos), Prp(0), c_err(0),Prp2(0);
	//bool res(false);
	for (; pos< tolerance; pos++){
		if (pos > mthL) {	break;	}
		c_err=0;Prp=0; Prp2=pos;
		do {
			
#ifdef _NEWMATCH			
			//new vector based matching
			c_err += DNA_IUPAC[Seq[Prp2]+256*PrSt[Prp]]; if (c_err > Err){break;}
#else
			//old, direct match
			if (!matchDNA(Seq[Prp2],PrSt[Prp])){c_err++;if (c_err > Err){break;}}
#endif
			Prp++; Prp2++;
		} while ( Prp < PrL);
		if (c_err<=Err ){endPos=pos;break;}
	}
	//if(!suc){pos=-1;}
	return endPos;
}
void DNA::reset() { 
	AccumError = 0.; goodQual = false; midQual = false;
	//FtsDetected.reset();
	avgQual = -1.f; Qsum = 0; tempFloat = 0.f;
	QualTraf = "";

	this->resetTruncation(); 
}

void DNA::reverse_transcribe() {
	reverseTS(Seq);
	std::reverse(Qual.begin(), Qual.end());
	AccumError = 0.; goodQual = false; midQual = false;
	avgQual = -1.f; Qsum = 0; tempFloat = 0.f;
	QualTraf = "";
	
}

//match from end of Seq to find rev primer
int DNA::matchSeqRev(std::string PrSt,int Err, int check_l,  
				  int coverage){
	//fail::ret -1
	int PrL = (int) PrSt.length();
	if (coverage==0){coverage=5;} //default seed set to 5
	int SeL = (int) Seq.size();
	//int wantSc = PrL - Err;
	int pos(SeL-coverage), Prp(0), c_err(0),endPos(-1);
	for (; pos> check_l; pos--){
		c_err=0;Prp=0;
		int PrL2 = min(PrL,SeL-pos);
		do {
#ifdef _NEWMATCH
			c_err += DNA_IUPAC[Seq[pos+Prp]+256*PrSt[Prp]]; if (c_err > Err){break;}
#else
			if(!matchDNA(Seq[pos+Prp],PrSt[Prp])){c_err++;if (c_err > Err){break;}	}
#endif
			Prp++; 
		} while (Prp < PrL2);
		if (c_err<=Err ){
			endPos=pos;break;}
	}
	//secondary check for last few NT's
	if (endPos==-2){
		pos = (SeL-1);
		for (; pos> (SeL-coverage); pos--){
			c_err=0;Prp=0;
			int PrL2 = min(PrL,SeL-pos);
			do {
#ifdef _NEWMATCH
				c_err += DNA_IUPAC[Seq[pos+Prp]+256*PrSt[Prp]]; if (c_err > Err){break;}
#else
				if(!matchDNA(Seq[pos+Prp],PrSt[Prp])){c_err++;	if (c_err > Err){break;}}
#endif
				Prp++; 
			} while (Prp < PrL2);
			if (c_err<=Err ){
				endPos=pos;break;}
		}
	}
	return endPos;
}
// looks through total DNA seq
int DNA::matchSeq_tot(std::string Pr,int Err,int MaxPos, int& c_err){
	//const char* Pr = PrSt.c_str();
	int PrL = (int) Pr.length();
	int pos(0), Prp(0), Prp2(0);
	bool suc(false);
	for (pos=0; pos< MaxPos; pos++){
		c_err=0;Prp=0;Prp2=pos;
		do {
			if(Seq[Prp2]!=Pr[Prp]){
				c_err++;
			}
			Prp++; Prp2++;
		} while (c_err <= Err && Prp< PrL);
		if (c_err<=Err){suc=true;break;}
	}
	if(!suc){pos=-1;}
	return pos;
}


bool DNA::matchDNA(char t1,char t2){
	if (t1==t2){
		return true;
	}
	switch (t2){
		case 'N': return true;
		case 'R': if (t1=='A' || t1=='G' ) {return true;}break;
		case 'Y': if (t1=='T' || t1=='C' ) {return true;}break;
		case 'M': if (t1=='C' || t1=='A' ) {return true;}break;
		case 'K': if (t1=='T' || t1=='G' ) {return true;}break;
		case 'W': if (t1=='T' || t1=='A' ) {return true;}break;
		case 'S': if (t1=='C' || t1=='G' ) {return true;}break;
		case 'B': if (t1!='A') {return true;}break;
		case 'D': if (t1!='C' ) {return true;}break;
		case 'H': if (t1!='G' ) {return true;}break;
		case 'V': if (t1!='T' ) {return true;}break;
	}
	return false;
}
bool DNA::HomoNTRuns(int max){
	char lastC = Seq[0];
	int rowC=1;
	for (unsigned int i=1;i<length();i++){
		if (Seq[i]==lastC){
			rowC++;
			if (rowC>= max){
				return false;
			}
		} else {
			rowC = 1;
			lastC = Seq[i];
		}
	}
	return true;
}

/*
void DNA::writeSeq(ofstream& wr){
	int cnt=0;
	if (Seq.size()==0){return;}
	wr<<">"<<NewID<<std::endl;
	for (unsigned int i=0;i<Seq.size();i++){
		cnt++;
		if (cnt<80){wr<<Seq[i];
		} else {wr<<Seq[i] << endl;	cnt=0;
		}
	}
	if (cnt>0){
		wr<<endl;
	}
}
*/

string DNA::xtraHdStr(){
	string xtr = " ";
	if (FtsDetected.forward){ xtr += "F"; }
	if (FtsDetected.reverse){ xtr += "R"; }
	return xtr;
}
void DNA::writeSeq(ostream& wr, bool singleLine) {
	if (Seq.size()==0){return;}
	wr<<">"<<NewID<<std::endl;
	if (singleLine) {
		wr << Seq.substr(0,length()) << std::endl;
		return;
	}
	unsigned int SeqS = length();
	unsigned int leftOver = SeqS%80;
	int last = 0;
	for (unsigned int i=0;i<(SeqS-leftOver)/80;i++){
		wr<<Seq.substr(last,80)<<endl;last+=80;
	}
	if (leftOver>0){
		wr<<Seq.substr(last)<<endl;
	}
}
/**/
void DNA::writeQual(ostream& wr, bool singleLine) {
	int cnt=0;
	if (Qual.size()==0){return;}
	wr<<">"<<NewID<<endl;
	//wr<< QualTraf<<endl;
	string endlChr("\n");
	if (singleLine) { endlChr = " "; }
	for (unsigned int i=0;i<length();i++){
		cnt++;
		if (cnt<80){
			wr<<Qual[i] << " ";
		} else {
			wr << Qual[i] << endlChr;
			cnt=0;
		}
	}
	if (cnt>0){
		wr<<endl;
	}

}

void DNA::writeFastQ(ostream& wr,bool newHD){//, int fastQver){

	if (Qual.size()==0 || Seq.length()==0 || length()==0){return;}
	string xtr = xtraHdStr(); 
	if (FtsDetected.forward){ xtr += "F"; }
	if (FtsDetected.reverse){ xtr += "R"; }
	if (newHD) {
		wr << "@" << NewID << endl;
	} else {
		wr << "@" << ID << endl;
	}
	wr<<Seq.substr(0,length())<<endl;
	wr <<"+"<<endl;//NewID<<endl;
	//char* QualTraf = new char[Qual.size()+1];
	wr<< QualTraf<<endl;
	//delete [] QualTraf;
}
void DNA::writeFastQ(ofbufstream& wr, bool newHD) {//, int fastQver){

	if (Qual.size() == 0 || Seq.length() == 0 || length() == 0) { return; }
	string xtr = xtraHdStr();
	if (FtsDetected.forward) { xtr += "F"; }
	if (FtsDetected.reverse) { xtr += "R"; }
	if (newHD) {
		wr << "@" + NewID + "\n";
	}else {
		wr << "@" + ID + "\n";
	}
	wr << Seq.substr(0, length()) + "\n";
	wr << "+\n";//NewID<<endl;
	wr << QualTraf + "\n";
}
void DNA::writeFastQEmpty(ostream& wr) {

	wr << "@" << NewID << endl;
	wr << "" << endl;
	wr << "+" << endl;//NewID<<endl;
	//char* QualTraf = new char[Qual.size()+1];
	wr << "" << endl;
	//delete [] QualTraf;
}

void DNA::changeHeadPver(int ver){
	string& oID(ID);
	if (IDfixed){
		oID = NewID;
	} 
	
	if (ver==1){ //change from XX 1: to XX/1
		size_t pos = oID.find_first_of(" ");
		//unsigned int pos2 = ID.find(" ",pos);
		if (pos != string::npos){			
			NewID = oID.substr(0,pos) + "/" + oID.substr(pos+1,1);
			//if (pos2 != string::npos){NewID += ID.substr(pos2);}
		}
	} else if (ver==2){
		cerr<<"Head change Not implemented\n";exit(23);
	}
	IDfixed=true;
}

/*void DNA::reverseTranscribe(){
	int qs = (int)Qual.size()-1;
	vector<int> Q2(Qual);
	for (int i=qs;i>=0;i--){
		Qual[i] = Q2[qs-i];
	}
	reverseTS(Seq);
}
*/
bool DNA::sameHead(shared_ptr<DNA> d){
	if (d == NULL){return false;}
	return sameHead(d->getIDshort());
}
bool DNA::sameHead(const string& oID) {
	size_t pos = getShorterHeadPos(ID);
	if (oID.size() < pos) { return false; }
	if (ID.substr(0, pos) == oID.substr(0, pos)) {
		return true;
	}
	return false;
}


void DNA::setPassed(bool b){
	goodQual=b;
	if (goodQual && midQual) {
		midQual = false;
	}
}
int DNA::getBCnumber() {
	//if (Sample==-1 )
	return Sample; 
}

void DNA::prepareWrite(int ofastQver) {
	uint len = length();
	if (QualTraf.size() == len) {
		return;
	}
	QualTraf.resize(len);
	unsigned int i = 0;
	for (; i < len; i++) {
		QualTraf[i] = char(Qual[i] + ofastQver);
	}
	QualTraf[i] = '\0';
}
void DNA::resetQualOffset(int x, bool fqSol) {
	for (size_t i = 0; i < Qual.size(); i++) { Qual[i] += x; }
	if (fqSol) {//quick hack, since q<13 scores are the only deviants and uninteresting in most cases..
		for (size_t i = 0; i < Qual.size(); i++) {
			if (Qual[i] < 0) {
				Qual[i] = 0;
			}
		}
	}
}

///////////////////////////////////////////////////////////////
//INPUT STREAMER



void DNAunique::addSmpl(int k) {
	if (k < 0) {
		return;
	}
	Count++;
#ifdef _MAPDEREPLICATE
	unordered_map<int, int>::iterator smID = occurence.find(k);
	if (smID == occurence.end()) {
		occurence[k] = 1;
	} else {
		smID->second++;
	}
#endif
}
void DNAunique::setOccurence(int smpl, int N) {
	unordered_map<int, int>::iterator smID = occurence.find(smpl);
	if (smID == occurence.end()) {
		occurence[smpl] = N;
	} else {
		smID->second += N;
	}
}
/*vector<pair<int, int>> DNAunique::getDerepMapSort2(size_t wh ){
	typedef std::pair<int, int> mypair;
	size_t siz = occurence.size();
	if (wh > siz) { wh = siz; }

	struct IntCmp {
		bool operator()(const mypair &lhs, const mypair &rhs) {
			return lhs.second > rhs.second;
		}
	};


	vector<mypair> myvec(occurence.begin(), occurence.end());
	std::partial_sort(myvec.begin(), myvec.begin() + wh, myvec.end(), IntCmp());

	return myvec;
}*/

bool sortDescending(int i, int j) { return (i>j); }//descending sort
vector<int>  DNAunique::getDerepMapSort(size_t wh) {
	vector<int> vals;
	size_t siz = occurence.size();
	if (wh > siz) {	wh = siz;}
	vals.reserve(siz);
	for (auto kv = occurence.begin(); kv != occurence.end(); kv++) {
		vals.push_back(kv->second);
	}
	//partial sort doesn't make sense, as I want to break border asap
	//partial_sort(vals.begin(), vals.begin() + wh, vals.end(), sortDescending);
	sort(vals.begin(), vals.end(), sortDescending);
	return vals;
}

bool DNAunique::pass_deprep_smplSpc(const vector<int>& cv) {
	unordered_map<int, int> occ;
	//combined samples will not be considered
	//occ = occurence;
	for (std::unordered_map<int, int>::iterator iter = occurence.begin(); iter != occurence.end(); ++iter) {
		//int cnts = iter->second;
		int ref = cv[iter->first];
		if (ref != -1 && iter->second >= ref ) {
			return true;
		}
	}
	return false;

}


void DNAunique::transferOccurence(shared_ptr<DNAunique> odu) {
	if (occurence.size() == 0) {
		occurence = odu->occurence;
		Count = odu->Count;
	} else {
		//which sample contains this dna?
		unordered_map<int, int> oldMap = odu->getDerepMap();
		unordered_map<int, int>::iterator smID;
		for (std::unordered_map<int, int>::iterator oID = oldMap.begin(); oID != oldMap.end(); ++oID) {
			smID = occurence.find(oID->first);
			if (smID == occurence.end()) {
				occurence[oID->first] = oID->second;
			} else {
				smID->second += oID->second;
			}
		}
		//size track
		Count += odu->Count;
	}
}
void DNAunique::writeMap(ofstream & o, const string & hd, 
	vector<int> & cntspersmpl, const vector<int>& combiID) {
	if (occurence.size() == 0) { return; }
	int totCnt(0);
	unordered_map<int, int> occ;
	if (combiID.size() > 0){//combine all counts on combined categories
		std::unordered_map<int, int>::iterator fnd;
		for (std::unordered_map<int, int>::iterator iter = occurence.begin(); iter != occurence.end(); ++iter) {
			//aim: occ[combiID[iter->first]] += iter->second;
			fnd = occ.find(combiID[iter->first]);
			if (fnd == occ.end()){
				occ[combiID[iter->first]] = iter->second;
			}else{
				fnd->second += iter->second;
			}
		}
	}
	else {
		occ = occurence;
	}

	//prints combined sample counts
	o << hd;
	for (std::unordered_map<int, int>::iterator iter = occ.begin(); iter != occ.end(); ++iter) {
		int cnts = iter->second;
		totCnt += cnts;
		o << "\t";
		o << iter->first << ":" << cnts;
	}
	//counts non-combined sample counts
	for (std::unordered_map<int, int>::iterator iter = occurence.begin(); iter != occurence.end(); ++iter) {
		//int cnts = iter->second;
		cntspersmpl[iter->first] += iter->second;
	}
	o << endl;
	
	if (totCnt != Count) {
		cerr << "Unequal Counts in Map("<<totCnt<<") and HeadCnt(" << Count << "):"<<endl<< this->getID()<<endl;
		exit(82);
	}
}
void DNAunique::Count2Head(bool usFmt) {
	//string tmp = this->getIDshort();
	if (!usFmt) {
		NewID += "_" + itos(Count);
	} else {
		NewID += ";size=" + itos(Count) + ";";
	}
	IDfixed = true;
}


///////////////////////////////////////////////////////////////
//INPUT STREAMER

InputStreamer::~InputStreamer(){
	allStreamClose();
//	for (uint i = 0; i < tdn1.size(); i++) { if (tdn1[i] != NULL) { delete tdn1[i]; } }
//	for (uint i = 0; i < tdn2.size(); i++) { if (tdn2[i] != NULL) { delete tdn2[i]; } }
}

bool InputStreamer::getFastaQualLine(istream&fna,  string&line) {
	
	if (!safeGetline(fna, line)) { return false; }
	while (line[0] == '$') { //$ marks comment
		safeGetline(fna, line);
	}
	return true;
}
bool InputStreamer::read_fasta_entry(istream&fna,istream&qual,shared_ptr<DNA> tdn1, shared_ptr<DNA>tdn2,int &cnt){
	
	if(fna.eof()){return false;}
	
	//int in_int; //char in_char;
	//int cnt=0;
	string tqual(""),tseq("");
	string line(""), lineQ(""); 
	if (!getFastaQualLine(fna,  line)) { return false; }
	if (line == "" && fna.eof()){return false;}
	if (!qualAbsent) {		getFastaQualLine( qual,  lineQ);	}
	cnt++;

	if (cnt == 1) { //fasta description
		if (line[0] != '>' && fna) { cerr << "ERROR: Line 1 in fasta file does not start with \">\" \n"; exit(23); }
		//new DNA, set up in tdn1
			tdn1->newHEad(line.substr(1));
			if (!getFastaQualLine(fna, line)) { return false; }
			if (!qualAbsent) { getFastaQualLine(qual, lineQ); }
			cnt++;
	}
	//continous read in until ">" is hit
	while (line[0] != '>') {
		tseq += line;
		if (!getFastaQualLine(fna, line)) { break; }
	}
	if (!qualAbsent) { 
		while (lineQ[0] != '>') {
			tqual += " " + lineQ;
			if (!getFastaQualLine(qual, lineQ)) { break; }
		}
	}

	//fna
	tdn1->setSeq(tseq);
	size_t lsize = tseq.size();
	vector<int> Iqual(lsize, 0);
	//qual
	if (!qualAbsent) {
		const char* lQ = rtrim(tqual).c_str();
		uint ii = 0; int nn(0);
		for (; ii < lsize; ii++) {
			nn = parseInt(&lQ);// , posStr);
			Iqual[ii] = nn;
			if (*lQ == '\0') {
				break;
			}
			//issQ >>  Iqual[ii];
		}
		if (ii != (lsize - 1)) {
			cerr << "ERROR: quality counts (" << ii << ") not the same length as DNA base counts in Sequence (" << (lsize - 1) << ")\n" << tdn1->getID() << "\n";
			exit(54);
		}
	}
	tdn1->Qappend(Iqual);


	//since already read in 1 more line, this line needs to be used on new object
	if (fna ) {
		if (!qualAbsent && (line[0] != '>' || lineQ[0] != '>')) { cerr << "ERROR: Desynced fasta reader\n"; exit(23); }
		tdn2->newHEad(line.substr(1));
	} else {
		return false;
	}
	//a new DNA obj was set up, return to process tdni, tdno will be completed next call
	return true;
}
inline int InputStreamer::parseInt(const char** p1){//,int& pos){//, const char ** curPos) {
	/*from http://stackoverflow.com/questions/5830868/c-stringstream-is-too-slow-how-to-speed-up*/
	//size_t nxtPos = input.find_first_of(' ', curPos);
	//const char *p = input.substr(curPos,nxtPos).c_str();
	//if (!*p || *p == '?')		return 0;
	//int s = 1;
	//p = (const char*)curPos;
	const char* p = *p1;
	while (*p == ' ') p++;

	int acc = 0;
	while (*p >= '0' && *p <= '9')
		acc = acc * 10 + *p++ - '0';


	*p1 = p;
	//curPos = (size_t) p;
	return acc;
}
void InputStreamer::jmp_fastq(istream & fna, int &lnCnt) {
	string line;
	string tname = "", tseq = ""; //temporary storage
	size_t cnt = 0, qcnt = 0, DNAlength = 0;
	bool mode = true; //mode=T:fna,F:qual
	bool needsAT = true, needsPlus = false; // checks if quality was completly read in and a '@' char is now expected in the next line
	//	getline(fna2,line2,'\n');

	while (getline(fna, line, '\n')) {
		lnCnt++;
		if (line == "") { return ; }
		if (line[0] == '@' && needsAT) { //fasta description
			if (cnt != 0) {
				cerr << "Line " << lnCnt << ": Could not find \'@\' when expected: on line \n" << line << endl;
			}
			//tname = line;
			needsAT = false;
			needsPlus = true;
			continue;
		} else if (line[0] == '+' && needsPlus) {
			needsPlus = false;
			mode = false;
			continue;
		} else if (needsAT) {
			cerr << "Line " << lnCnt << " :Could not find \'@\' symbol where expected on line \n " << line << endl;
			exit(6);
		}

		istringstream iss(line);
		iss >> tseq;
		if (mode) {
			DNAlength += tseq.length();
		} else {
			qcnt += tseq.length();
			if (DNAlength == qcnt) { return; }
		}

		cnt++;
	}
}
//reads out single seq + quality entry from fastq file
shared_ptr<DNA> InputStreamer::read_fastq_entry(istream & fna,  int &minQScore, int& lnCnt,
			bool&corrupt,bool repairInStream) {
	string line;
	//string tseq = ""; //temporary storage
	uint qcnt = 0;
	bool mode = true; //mode=T:fna,F:qual
	shared_ptr<DNA> tdn = make_shared<DNA>("", ""); 
	vector<qual_score> Iqual(0);
	bool needsAT = true, needsPlus = false; // checks if quality was completly read in and a '@' char is now expected in the next line

	if (!fna) { return NULL;  cerr << "Read Fastq: Input stream does not exist" << lnCnt << endl; exit(53); }

	while (safeGetline(fna, line)) {
		lnCnt++;
		while (repairInStream) {
			if (line[0] == '@') {repairInStream = false;}
			else {safeGetline(fna, line);}
		}
		//if (line.length() == 0) { return tdn; }
		if (needsAT) { //fasta description
			if (line[0] == '@') {
				if ((lnCnt-1) % 4 != 0) {
					fqPassedFQsdt = false;
				}
				tdn->newHEad(line.substr(1));
				mode = true; qcnt = 0;
				needsAT = false;
				needsPlus = true;
				continue;
			} else {
				corrupt = true;//try again,could be last empty line in file..
				if (line.length() != 0) {
					IO_Error("Line " + itos(lnCnt) + " :Could not find \'@\' symbol where expected on line \n " +line);// << endl;
				}				
				return tdn;
			}
		} else if (needsPlus && line[0] == '+') {
			if ((lnCnt-3) % 4 != 0) {
				fqPassedFQsdt = false;
			}

			
			Iqual.resize(tdn->length(), 0);
			needsPlus = false;
			mode = false;
			continue;
		}

		//istringstream iss(line);
		//iss >> tseq;
		if (mode) {
			//fna
			if ((lnCnt-2) % 4 != 0) {
				fqPassedFQsdt = false;
			}

			//if (std::any_of(line.begin(), line.end(), [](char c) {return (islower(c)); })) {
			if (any_lowered(line)){
				fqPassedFQsdt = false;
				std::transform(line.begin(), line.end(), line.begin(), ::toupper);
			}
			tdn->append(line);
		} else if (!mode) {
			//qual	
			if ((lnCnt ) % 4 != 0) {
				fqPassedFQsdt = false;
			}

			for (size_t i = 0; i < line.length(); i++) {
				//really 33?
				Iqual[qcnt] = minmaxQscore((qual_score)line[i] - fastQver);
				qcnt++;
			}
			if (qcnt == tdn->length()) {
				//needsAT=true;
				tdn->Qappend(Iqual);
				
			} else if (line.length() + qcnt > tdn->length()) {
				//check that quality gets not more length than DNA
				IO_Error("ERROR: More quality positions than nucleotides detected for sequence\n " + tdn->getID());
				//tdn->setPassed(false);
				corrupt = true;
				return tdn;
			}
			break;

		}
	}
	//if (tdn!= NULL && !tdn->control()){delete tdn; tdn=NULL;}
	corrupt = false;
	return tdn;

	//to check for fast fastq reader: 1. DNA in uppercase? 2. DNA/QUAL in single line?
}
void InputStreamer::IO_Error(string x) {
	cerr << x << endl;
	if (DieOnError) {
		exit(632);
	}
	ErrorLog.push_back(x);
}
//reads out single seq + quality entry from fastq file
shared_ptr<DNA> InputStreamer::read_fastq_entry_fast(istream & fna, int& lnCnt, bool& corrupt) {
	string line;
	if (!fna) { return NULL; cerr << "Read Fastq_f: Input stream does not exist " << lnCnt << ".\n"; exit(53); }

	if (!safeGetline(fna, line)) { return NULL; }
	shared_ptr<DNA> tdn = make_shared<DNA>("", "");
	if (line.length() == 0) { return tdn; }
	while (line[0] != '@') {
		IO_Error("ERROR on line " + itos(lnCnt) + ": Could not find \'@\' when expected (file likely corrupt, trying to recover):\n" + line);// << endl;
		//exit(55);
		//recover instead and go to next entry..
		corrupt = true;
		if (!safeGetline(fna, line)) { corrupt = true;  return NULL; }//delete tdn;
	}
	tdn->newHEad(line.substr(1));
	//cerr << line.substr(1);
	if (!safeGetline(fna, tdn->getSeq())) { corrupt = true;  return NULL; }//delete tdn;
	//if (line.length() == 0) { return NULL; }
	//std::transform(line.begin(), line.end(), line.begin(), ::toupper);
	//tdn->append(line);
	//"+"
	if (!safeGetline(fna, line)) { corrupt = true; return NULL; }//delete tdn;
	while (line[0] != '+') {
		//recovery is hard, just give up this read
		IO_Error("Error input line " + itos(lnCnt + 2) + ": Could not find \'+\' when expected (file likely corrupt, aborting):\n" + line);// << endl;
		corrupt = true;
		return tdn;
		//if (!safeGetline(fna, line)) { delete tdn; return NULL; }
	}

	//qual score
	vector<qual_score> Iqual(tdn->mem_length(), 0);
	if (!safeGetline(fna, line)) { corrupt = true;  return NULL; }//delete tdn;
	uint qcnt(0); uint lline = (uint)line.length();
	for (; qcnt < lline; qcnt++) {
		Iqual[qcnt] = minmaxQscore((qual_score)line[qcnt] - fastQver);

	}

	if (qcnt == tdn->mem_length()) {
		//needsAT=true;
		tdn->setQual(Iqual);
	} else if (line.length() + qcnt != tdn->length()) {
		//check that quality gets not more length than DNA
		corrupt = true;
		IO_Error("Error input line " + itos(lnCnt + 3) + ": More quality positions than nucleotides detected for sequence\n " +tdn->getID());// << endl;
//		exit(7);
	}
	lnCnt+=4;
	//if (tdn!= NULL && !tdn->control()){delete tdn; tdn=NULL;}
	corrupt = false;
	return tdn;
}

int InputStreamer::minmaxQscore(qual_score t) {
		if (t < 0) { ////quick hack, since q<13 scores are the only deviants and uninteresting in most cases..
			if (fqSolexaFmt){
				if (t < -5) {
					cerr << "Unusually low sloexa quality score (" << t << "); setting to 0.\n";
				}
			} else {
				if (t >= -5) {
					cerr << "Resetting auto format to Solexa (illumina 1.0-1.3) format.\n";
					fqSolexaFmt = true;
				} else {
					cerr << "Unusually low quality score (" << t << "); setting to 0.\n";
				}
			}
			t = 0; 
		}
	if (minQScore > t) { 
		minQScore = t; 
		if (minQScore < 0) {
		}
	} else if (maxQScore < t) { 
		maxQScore = t; 
	}
	return t;
}
bool InputStreamer::checkInFileStatus() {
	for (uint i = 0; i < 3; i++) {
		if (fnaRead) {
			if (fna_u[i] != NULL && *fna_u[i]) {
				return true;
			}
		} else {
			if (fastq_u[i] != NULL && *fastq_u[i]) {
				return true;
			}
		}
	}
	return false;
}
void InputStreamer::allStreamReset() {
	resetLineCounts();
#ifdef DEBUG
	cerr << "Resetting input streams" << endl;
#endif
	//reopen streams in gz case // sdm 1.01: make default
	if (true || openedGZ) { 
		allStreamClose();
		setupFastq_2(inFiles_fq[0], inFiles_fq[1], inFiles_fq[2]); 
		setupFastaQual2(inFiles_fna[0], inFiles_qual[0]);
	} else {
		for (uint i = 0; i < 3; i++) {
			if (fna_u[i] != NULL && * (fna_u[i])) { fna_u[i]->clear(); fna_u[i]->seekg(0, ios::beg); }
			if (qual_u[i] != NULL && *(qual_u[i])) { qual_u[i]->clear(); qual_u[i]->seekg(0, ios::beg); }
			if (fastq_u[i] != NULL && *(fastq_u[i])) { fastq_u[i]->clear(); fastq_u[i]->seekg(0, ios::beg); }
		}
	}
	//checkInFileStatus();
}
void InputStreamer::allStreamClose(){
	for (uint i = 0; i < 3; i++){
/*		if (*(fna[i])){ fna[i]->close(); }
		if (*(qual[i])){ qual[i]->close(); }
		if (*(fastq[i])){ fastq[i]->close(); }
#ifdef _gzipread	
		if (gzfna[i]){ gzfna[i].close(); }
		if (gzqual[i]){ gzqual[i].close(); }
		if (gzfastq[i]){ gzfastq[i].close(); }
#endif
		*/
		if (fna_u[i] != NULL) { delete fna_u[i]; }  fna_u[i] = NULL;
		if (qual_u[i] != NULL) { delete qual_u[i]; } qual_u[i] = NULL;
		if (fastq_u[i] != NULL) { delete fastq_u[i];   }fastq_u[i] = NULL;
	}

	if (!fnaRead && minQScore < 1000) {
			maxminQualWarns_fq( );
	}
}
void InputStreamer::jumpToNextDNA(bool&stillMore, int pos) {
	if (fnaRead) {//get DNA from fasta + qual files
		//shared_ptr<DNA> ret;
		stillMore = read_fasta_entry(*(fna_u[pos]), *(qual_u[pos]), tdn1[pos], tdn2[pos], lnCnt[pos]);
		//tdn1 will be completed, tdn2 will have a header set up
		//ret = tdn1[pos];
		tdn1[pos] = tdn2[pos]; tdn2[pos].reset(new DNA("", ""));
	} else {
		jmp_fastq(*fastq_u[pos], lnCnt[pos]);
		if (!*(fastq_u[pos])) {
			stillMore = false;
		}
	}
}

shared_ptr<DNA> InputStreamer::getDNA(bool& stillMore, int pos, bool& sync){
	//if (sync) {
	//	while (desync(pos)) {
	//		jumpToNextDNA(stillMore, pos);
	//	}
	//}
	if (pos == 1 && numPairs <= 1) {
		return NULL;
	}
	shared_ptr<DNA> ret;
	bool corrupt(true); //corrupt state isn't implemented for fnaread
	
	bool repairInStream(false);
	while (corrupt) {
		if (fnaRead) {//get DNA from fasta + qual files
			stillMore = read_fasta_entry(*(fna_u[pos]), *(qual_u[pos]), tdn1[pos], tdn2[pos], lnCnt[pos]);
			corrupt = false;
			//tdn1 will be completed, tdn2 will have a header set up
			ret = tdn1[pos];
			tdn1[pos] = tdn2[pos]; tdn2[pos].reset( new DNA("", ""));
			if (!ret->seal() || ret->isEmpty()) { ret = NULL; }
			if (pos == 0 && pairs_read[pos] % 100 == 0) {
				_drawbar(*(fna_u[pos]));
			}
			else 	if (!stillMore) { _drawbar(*(fna_u[pos])); }
		}
		else { //fqRead
			if (fqReadSafe) {
				ret = read_fastq_entry(*(fastq_u[pos]), minQScore, lnCnt[pos], corrupt, repairInStream);
				if (fastQver == 0 && ret->length() > 5 && !corrupt) {//autodetect
					ret->resetQualOffset(auto_fq_version(), fqSolexaFmt);
					//reset streams
				}
				if (lnCnt[pos] > 100) {//tmp set back to 500
					if (fqPassedFQsdt) {
						fqReadSafe = false;
#ifdef DEBUG
						cerr << "Switching to fast fastq reader..\n ";
#endif
					}
				}
			}
			else {
				ret = read_fastq_entry_fast(*(fastq_u[pos]), lnCnt[pos], corrupt);
			}
			if (!stillMore || fastq_u[pos]->eof() || (!*(fastq_u[pos])) ) {
				if (ret != NULL) { if (!ret->seal() || ret->isEmpty()) { ret = NULL; } } //delete ret;
				stillMore = false; break;
			} else if (ret == NULL || !ret->seal() || ret->isEmpty()) {
				corrupt = true;
			}
			if (pos == 0 && pairs_read[pos] % 100 == 0) {
				if (_drawbar(*(fastq_u[pos]))) { stillMore = false; break; }
			}
			else 	if (!stillMore) { _drawbar(*(fastq_u[pos])); }

			
			if (corrupt) {
				//delete ret; 
				ret = NULL;
				sync = true;
				repairInStream = true;
			}
		}
		pairs_read[pos]++;
		//last check
	}
	//
	return ret;
}

void InputStreamer::openMIDseqs(string p,string in){
	if (in==""){return;}
	
	if (fastq_u[2] != NULL){
		cerr << "MID file was already initialized" << endl;
	}
#ifdef DEBUG
	cerr << "Open Mid Seq file" << endl;
#endif

	string file_type = "MID specific fastq";
	string tmp = (p + in);
	inFiles_fq[2] = tmp;
	if (isGZfile(tmp)){
#ifdef _gzipread
		fastq_u[2] = new igzstream(tmp.c_str(), ios::in);
		file_type = "MID specific gzipped fastq";
#else
		cerr << "gzip not supported in your sdm build\n" << tmp; exit(50);
#endif
	}
	else {
		fastq_u[2] = new ifstream(tmp.c_str(), ios::in);
	}
	if (!*(fastq_u[2])){ cerr << "\nCouldn't find " << file_type<<": " << in << " !\n Aborting..\n";		exit(4); }

	hasMIDs=true;
}

bool InputStreamer::setupFastq(string path, string fileS, int& pairs, string subsPairs,bool simu) {
	allStreamClose(); openedGZ = false;
	minQScore = 1000;	maxQScore = -1;	fqSolexaFmt = false;
	resetLineCounts();
	vector<string> tfas = splitByCommas(fileS);
	if ( pairs == -1 ) {
		pairs = (int) tfas.size();
		if ( pairs > 1 ) { cerr << "Paired input (\",\" separated) detected\n"; }
	}
	numPairs = pairs;
	string p1(""), p2(""), midp("");
	string xtraMsg = "";
	if (getCurFileN() > 0) {
		xtraMsg = " " + itos(getCurFileN()) + " of " + itos(totalFiles);
		if (BCnumber > 1) {
			xtraMsg = ", looking for " + itos(BCnumber) + "BCs.\n";
		} else {
			xtraMsg = ".\n";
		}
	}

	if (tfas.size() != (uint)pairs && subsPairs == "") {
		cerr << "Unequal number of files (" << tfas.size() << ") and option-set paired files (" << pairs << ").\n Aborting...\n"; exit(52);
	}
	if (tfas.size() == 3) {
		if (tfas.size() != 3) { cerr << "Could not detect 3 input files in string\n" << fileS << "\n Aborting.." << endl; exit(76); }
		midp = path + tfas[1];
		p1 = path + tfas[0];
		p2 = path + tfas[2];
//		cerr << p1 << " + " << p2 << " and " << midp << endl;
	} else if (tfas.size() == 2) {
		p1 = path + tfas[0];
		p2 = path + tfas[1];
	} else if (tfas.size() == 1) {
		p1 = path + fileS;
		//cerr << "Reading fastq " << p1 << endl;
	}
	if (subsPairs == "1") {
		p2 = "";
	} else 	if (subsPairs == "2") {
		p1 = p2; p2 = "";
	}

	if (simu) {
		return fileExists(p1) && fileExists(p2) && fileExists(midp);
	}


	if (p1 != ""&&p2 != "") {
		if (midp != "") {
			cerr << "Reading paired fastq + MID file" << xtraMsg<<"."<<endl;
		} else {
			cerr << "Reading paired fastq" << xtraMsg << "." << endl;
			cerr << p1 << " + " << p2 << endl;
		}
	} else  {
		cerr << "Reading fastq" << xtraMsg << "." << endl;
		cerr << p1 << endl;
	}

	bool suc = setupFastq_2(p1, p2, midp);
	//read progress  bar setup
	_measure(*fastq_u[0]);
	//allStreamClose();
	//setupFastq_2(p1, p2, midp);
	return suc;
}

// Prints out the progress bar
inline void InputStreamer::_print(int cur, float prog) {

	std::cerr << std::fixed << std::setprecision(2)
		<< "\r   [" << std::string(cur, '#')
		<< std::string(_max + 1 - cur, ' ') << "] " << 100 * prog << "%";

	if (prog == 1) std::cerr << std::endl;
	else std::cerr.flush();

}
inline bool InputStreamer::_drawbar(istream & tar) {
	if (_fileLength <= 0) { return false; }
	int pos((int)tar.tellg());
	float prog(pos / float(_fileLength)); // percentage of infile already read
	if (pos == -1 || prog > 1.f) { _print(_max + 1, 1); _fileLength = 0;  return false; }

	// Number of #'s as function of current progress "prog"
	int cur((int) ceil(prog * (float) _max));
	if (_last != cur) _last = cur, _print(cur, prog);
	if (prog == 1.f) {
		return true;
	}
	return false;

}
inline void InputStreamer::_measure(istream& tar) {
	tar.seekg(0, ios_base::end);
	_fileLength = (int) tar.tellg();
	tar.seekg(0, ios_base::beg);
	tar.clear();
}

bool InputStreamer::setupFastq_2(string p1, string p2, string midp) {
	//setupFastq_2(inFiles_fq[0],inFiles_fq[1],inFiles_fq[2]);
#ifdef DEBUG
	cerr << "setupFastq2 " << p1<<endl<<midp << endl;
	//cerr << inFiles_fq.size() << endl;
#endif
	string file_type = "";

	//        INPUT   files
	if (p1 != ""){ // first file exists
		inFiles_fq[0] = p1;
		file_type = "fastq file 1";
		if (isGZfile(p1)){
			openedGZ=true;
#ifdef _gzipread
			fastq_u[0] = new igzstream(p1.c_str(), ios::in);
			file_type = "gzipped fastq file 1";
#else
			cerr << "gzip not supported in your sdm build\n" << p1; exit(50);
#endif
		}
		else { fastq_u[0] = new ifstream(p1.c_str(), ios::in); }
		if (!*(fastq_u[0])){ cerr << "\nCouldn't find "<<file_type<<" " << p1 << "!\n Aborting..\n";		exit(4); }
	}
	//second pair
	if (p2 != ""){
		inFiles_fq[1] = p2;
		file_type = "fastq file 2";
		if (isGZfile(p2)){
			openedGZ=true;
#ifdef _gzipread
			fastq_u[1] = new igzstream(p2.c_str(), ios::in);
			file_type = "gzipped fastq file 2";
#else
			cerr << "gzip not supported in your sdm build\n" << p2; exit(50);
#endif
		}
		else { fastq_u[1] = new ifstream(p2.c_str(), ios::in); }
		if (!*(fastq_u[1])){ cerr << "\nCouldn't find "<<file_type<<" " << p2 << "!\n Aborting..\n";		exit(4); }
	}
	//MID file
	if (midp != ""){
		this->openMIDseqs("", midp);
	} 
	return true;
}
//mainFile = IS->setupInput(path, uniqueFas[i], FastqF[tarID], FastaF[tarID], QualF[tarID], MIDfq[tarID], fil->isPaired(), cmdArgs["-onlyPair"]);
string InputStreamer::setupInput(string path, int i, int t, const vector<string>& uF, const vector<string>& FQ, const vector<string>& Fas,
	const vector<string>& Qual, const vector<string>& midf, int &paired, string onlyPair, 
	string& shMain, bool simu) {
	string mainFile("");
	if (fnaRead) {
		if (Fas[t] != uF[i]) {
			cerr << "Error in matching FASTA target filenames.\n";
			exit(11);
		}
		this->setupFastaQual(path, Fas[t], Qual[t], paired, onlyPair,simu);
		mainFile = path + Fas[t];
		shMain = Fas[t];
	} else {
		if (FQ[t] != uF[i]) {
			cerr << "Error in matching target filenames.\n";
			exit(11);
		}
		this->setupFastq(path, FQ[t], paired, onlyPair,simu);
		mainFile = path + FQ[t];
		shMain = FQ[t];
	}
	this->openMIDseqs(path, midf[t]);
	return mainFile;
}


bool InputStreamer::setupFastaQual(string path, string Sfil, string Qfil, int& pairs, string subsPairs, bool simu) {
	resetLineCounts();
	allStreamClose();
	vector<string> tfas = splitByCommas(Sfil);
	if ( pairs == -1 ) {
		pairs = (int) tfas.size();
		if ( pairs > 1 ) { cerr << "Paired input (\",\" separated) detected\n"; }
	}
	if (pairs > 1) { cerr << "Paired fasta+qual is currently not implemented\n"; exit(72); }
	if (subsPairs == "1") {
		//
	} else 	if (subsPairs == "2") {
		//
	}

	numPairs = pairs;
	tdn1[0].reset( new DNA("", "")); tdn2[0].reset( new DNA("", ""));

	inFiles_fna[0] = path + Sfil;
	inFiles_qual[0] = path + Qfil;
	if (simu) {
		if (!fileExists(inFiles_qual[0],-1,false)) {
			cerr << "Warning: corresponding quality file is missing: " << inFiles_qual[0]<<endl;
		}
		return fileExists(inFiles_fna[0]) ;
	}

	bool suc = setupFastaQual2(inFiles_fna[0], inFiles_qual[0]);
	//activate length measure
	_measure(*fna_u[0]);
	//xtra safety for cases where pointer isn't reset properly
	//allStreamClose();
	//setupFastaQual2(inFiles_fna[0], inFiles_qual[0]);
	return suc;
}
bool InputStreamer::setupFastaQual2(string fileS, string fileQ, string file_type) {
	string file_typeq = "quality file";

	//        INPUT   file
	if (fileS != "") {
		if (isGZfile(fileS)) {
#ifdef _gzipread
			file_type = "gzipped fasta file";
			fna_u[0] = new igzstream(fileS.c_str(), ios::in);
#else
			cerr << "gzip not supported in your sdm build\n" << fileS; exit(50);
#endif
		} else {
			//fna[0].open(fileS.c_str(), ios::in);
			fna_u[0] = new ifstream(fileS.c_str(), ios::in);
			
		}
		if (!*(fna_u[0])) { cerr << "\nCouldn't find " << file_type << " file \"" << fileS << "\"!\n Aborting..\n"; exit(4); }
		//char Buffer[RDBUFFER];
		//fna_u[0]->rdbuf()->pubsetbuf(Buffer, RDBUFFER);
	}
	//quality file
	if (fileQ != "") {
		if (isGZfile(fileQ)) {
#ifdef _gzipread
			qual_u[0] = new igzstream(fileQ.c_str(), ios::in);
			file_typeq = "gzipped quality file";
#else
			cerr << "gzip not supported in your sdm build \n" << fileQ; exit(50);
#endif
		} else if (fileQ == "") {
			//setup empty stream
			qual_u[0] = new ifstream();
			qualAbsent = true;
		} else {
			qual_u[0] = new ifstream(fileQ.c_str(), ios::in);
		}
		if (!*(qual_u[0])) {
			cerr << "\nCouldn't find " << file_typeq << " file \"" << fileQ << "\"!\n Running in no qual filter mode\n";
			qualAbsent = true;
			//exit(55);
		} 
		if (getCurFileN() > 0) {
			cerr << "Reading Fasta + Quality file " << getCurFileN() << " of " << totalFiles << ".\n";
		} else {
			cerr << "Reading Fasta + Quality file.\n";
		}
		cerr<< file_type << " : " << fileS << endl << file_typeq << " : " << fileQ << endl;
	} else if (fileS != "") {
		if (getCurFileN() > 0) {
			cerr << "Reading Fasta file " << getCurFileN() << " of " << totalFiles << ".\n";
		} else {
			cerr << "Reading Fasta.\n";
		}
		qual_u[0] = new ifstream();
		qualAbsent = true;
	}
	return true;
}
void InputStreamer::setupFna(string fileS){
	allStreamClose();
	resetLineCounts();
	numPairs = 1;
	tdn1[0].reset(new DNA("", "")); tdn2[0].reset(new DNA("", ""));
	setupFastaQual2(fileS, "","seq");
	/*cerr << "Reading Fasta file.\n" << fileS << endl;
	string file_type = "seq";
	//        INPUT   file
	if (isGZfile(fileS)){
#ifdef _gzipread
		fna_u[0] = new igzstream(fileS.c_str(), ios::in);
		file_type = "gzipped seq";
#else
		cerr << "gzip not supported in your sdm build"; exit(50);
#endif
	}
	else {
		//fna[0].open(fileS.c_str(), ios::in);
		fna_u[0] = new ifstream(fileS.c_str(), ios::in);
	}
	if (!fna_u[0]){ cerr << "\nCouldn't find "<<file_type<<" file \"" << fileS << "\"!\n Aborting..\n"; exit(4); }
	//set quality to empty read buffer
	qual_u[0] = new ifstream();
	qualAbsent = true;
	*/
}

int InputStreamer::auto_fq_version() {
	int fqDiff(0);
	if (minQScore >= 100 || maxQScore < 2) {
		return fqDiff;
	}
	fqSolexaFmt = false;
	if (minQScore >= 59 && maxQScore > 74){
		fqDiff = (fastQver - 64); fastQver = 64;
		if (minQScore < 64) { //set to illumina1.0 (solexa)
			fqSolexaFmt = true;
			cerr << "\nSetting to illumina 1.0-1.3 (solexa) fastq version (q offset = 64, min Q=-5).\n\n";
		} else {
			cerr << "\nSetting to illumina 1.3-1.8 fastq version (q offset = 64).\n\n";
		}
	} else if (minQScore >= 33 && maxQScore <= 74) {
		fqDiff = (fastQver - 33); fastQver = 33;
		cerr << "\nSetting to Sanger fastq version (q offset = 33).\n\n";
	} else {
		cerr << "\nUndecided fastq version..\n";
		fqDiff = (fastQver - 33); fastQver = 0;
		//exit(53);
	}
	QverSet = true;
	return fqDiff;
}

void InputStreamer::maxminQualWarns_fq(){
	if (minQScore >= 59-33 && fastQver==33 ){ //set to sanger version, but no low qual over whole dataset -> probably illumina version
		cerr << "     WARNING :: \nQuality scores in your dataset are unusually high (min Q=" << minQScore<<"). Please check that you have a fastQ file in NCBI SRA, Sanger or Illumina 1.8+ version.\nIf not, set fastqVersion in option file to \"3\" (Illumina 1.0) or \"2\" (Illumina 1.3 < 1.8) .\n\n";
	} 
	if (minQScore < 0 ){ //set to sanger version, but no low qual over whole dataset -> probably illumina version
		cerr << "     WARNING :: \nQuality scores in your dataset are unusually low (min Q=" << minQScore << "). Please check that you have a fastQ file in Illumina 1.0 or Illumina 1.3 < 1.8 version.\nIf not, set fastqVersion in option file to \"1\" (NCBI SRA, Sanger or Illumina 1.8+ version).\n\n";
	} 
}

#ifdef _gzipread2

std::vector< char > readline(gzFile f) {
	//  gzFile fp =gzopen(fname,"r");

	std::vector< char > v(2056);
	int pos = 0;
	for (;;) {
		if (gzgets(f, &v[pos], (int)v.size() - pos) == 0) {
			// end-of-file or error
			int err;
			//const char *msg = gzerror(f, &err);
			if (err != Z_OK) {
				// handle error
			}
			break;
		}
		unsigned read = strlen(&v[pos]);
		if (v[pos + read - 1] == '\n') {
			if (v[pos + read - 2] == '\r') {
				pos = pos + read - 2;
			}
			else {
				pos = pos + read - 1;
			}
			break;
		}
		if (read == 0 || pos + read < v.size() - 1) {
			pos = read + pos;
			break;
		}
		pos = v.size() - 1;
		v.resize(v.size() * 2);
	}
	v.resize(pos);
	return v;
}
#endif