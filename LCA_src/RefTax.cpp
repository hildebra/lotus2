#include "RefTax.h"

//generic functions
void trim(string& str,
	const std::string& whitespace)
{
	auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		strBegin = 0;

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	str = str.substr(strBegin, strRange);
}

bool isGZfile(const string fi) {
	string subst = fi.substr(fi.length() - 3);
	if (subst == ".gz") {
		return true;
	}
	return false;
}



TaxObj::TaxObj(TaxObj* t):SavedTaxs(t->SavedTaxs), Subj(t->Subj), perID(t->perID),depth(t->depth) {
	;
}

string TaxObj::getWriteString() {
	string ret(SavedTaxs[0]);
	for (int i = 1; i < (int) SavedTaxs.size(); i++) {
		if (i >= depth || SavedTaxs[i] == __unkwnTax) {
			ret += __defaultTaxSep + __unkwnTaxWR;
		} else {
			ret += __defaultTaxSep + SavedTaxs[i];
		}
	}
	return ret;
}



TaxObj::TaxObj(const string& X,int d, bool nativeSLV, bool doNotCheckTax):SavedTaxs(d, __unkwnTax), Subj(""), perID(0.f),depth(0){
	size_t fnd = 0;
	//size_t offset = 1;
	if (!nativeSLV) {
		fnd = X.find("__", 0) + 2;
	}
	
	
	int cnt(0);
	size_t strL = X.length();
	while (fnd != string::npos) {
		size_t f2 = X.find(";", fnd );
		string substr = X.substr(fnd , (f2-fnd));
		if (!nativeSLV) { fnd = X.find("__", fnd + 2)+2; } else {
			fnd = f2 + 1;
			if (fnd >= strL) { fnd = string::npos; }
		}
		trim(substr);
		string lcsubstr(substr); transform(lcsubstr.begin(), lcsubstr.end(), lcsubstr.begin(), ::tolower);
		bool taxKnown = substr.length() > 0 && lcsubstr != "unclassified" && lcsubstr != "uncultured bacterium"
			&& lcsubstr != "uncultured" && substr != "?";
		if ( doNotCheckTax || taxKnown) {
			depth = cnt+1;
			SavedTaxs[cnt] = substr;
		} else {
			//is already "?"
			;
		}
		cnt++;
		if (cnt >= d) {
			break;
		}
	}
}

bool TaxObj::evalAcpyTax(TaxObj* oth) {
	if (perID != 0.f && oth->perID != 0.f) {
		if (perID > (oth->perID *0.99) ) {
			return false;
		}
		else if ( (oth->perID * 0.985) > (perID ) ) { // real advantage
			if (oth->depth >= depth) { // also make sure that this is not a hit to "?"
				copyOver(oth);
				return true;
			}
		}
	}
	if (oth->depth > depth) {//replace tax
		copyOver(oth);
		return true;
	}
	return false;
}

void TaxObj::copyOver(TaxObj* oth) {
	SavedTaxs = oth->SavedTaxs;
	depth = oth->depth;
	perID = oth->perID;
}



RefTax::RefTax(const string& inF, int tdep, bool nativeSLV,bool checktaxStr):TaxFile(inF),
tlevels(tdep,"")
{
	//ini constants
	//#my @taxLvls = ("domain", "phylum", "class", "order", "family", "genus");
	//cout << "DEBUG\n\n"; return;
	cout << "Loading tax DB.." << endl;

	string line;
	ifstream in(inF.c_str());
	if (!in) { cerr << "Cant open file " << inF << endl; std::exit(11); }
	int TaxDbl(0), TaxSingl(0);
	while (getline(in, line, '\n')) {
		size_t dlmt = line.find("\t");
		string ID = line.substr(0, dlmt);
		//string ttax = line.substr(dlmt+1);
		TaxObj* t = new TaxObj(line.substr(dlmt + 1),7, nativeSLV, !checktaxStr);
		auto fnd = Tlink.find(ID);
		if (fnd == Tlink.end()){//all good
			Tlink[ID] = t; TaxSingl++;
		} else {//not good: tax is double annotated
			//cerr << "Tax ID: " << ID << " is double used!\n";
			//exit(930);
			TaxDbl++;
		}
	}
	//cerr << "C1\n";
	cout << TaxDbl << " of " << TaxDbl + TaxSingl << " are duplicate entries\n";
	//cerr << "C2\n";
	this->stats();
	//cerr << "C3\n";
}

RefTax::~RefTax()
{
	for (auto it = Tlink.begin(); it != Tlink.end(); ++it){
		//std::cout << " " << it->first << ":" << it->second;
		delete it->second;
	}
}
void RefTax::stats() {
	int cnt = 0; vector<int> hist(10, 0);
	int maxD = 0;
	for (auto it = Tlink.begin(); it != Tlink.end(); ++it) {
		int dep = it->second->depth -1; if (maxD < dep) { maxD = dep; }
		if (dep > 10) { cerr << "Tax depth " << dep << " of object " << it->first << " is too great!\n"; }
		hist[dep]++;
		cnt++;
	}
	cout << "TaxDB " << TaxFile << " contained " << cnt << " entries, depth distribution is:\n";
	for (int i = 0; i < maxD+1; i++) {
		cout << i << ":" << hist[i] << " ";
	}
	cout << endl;
}


//*******************************************************
//        BlastRes
//*******************************************************

BlastRes::BlastRes(const string& line, int inptFmt):
		Query(""), Sbj(""), alLen(0), perID(0.f), eval(100.f), score(0.f),
		fail(false){

	if (line.length()==0){ fail = true; return; }
	size_t f1 = line.find("\t");
	if (f1 == string::npos) {fail = true; return;}
	Query = line.substr(0, f1);
	//subject
	size_t fp = f1 + 1; f1 = line.find("\t", f1 + 1);
	if (f1 == string::npos) { fail = true; return; }
	Sbj = line.substr(fp, f1-fp);
	
	//id
	fp = f1 + 1; f1 = line.find("\t", f1 + 1);
	if (f1 == string::npos) { fail = true; return; }
	perID = atof(line.substr(fp, f1-fp).c_str());

	//alignmentLength
	fp = f1 + 1; f1 = line.find("\t", f1 + 1);
	if (f1 == string::npos) { fail = true; return; }
	alLen = atoi(line.substr(fp, f1-fp).c_str());
	//inserts
	f1 = line.find("\t", f1 + 1);
	//qstart
	f1 = line.find("\t", f1 + 1);
	//qstop
	f1 = line.find("\t", f1 + 1);
	//sstart
	f1 = line.find("\t", f1 + 1);
	//sstop
	f1 = line.find("\t", f1 + 1);
	f1 = line.find("\t", f1 + 1);
	//eval
	fp = f1 + 1; f1 = line.find("\t", f1 + 1);
	if (f1 == string::npos) { fail = true; return; }
	eval = atof(line.substr(fp, f1-fp).c_str());
	//score
	score = atof(line.substr(f1+1).c_str());

}


//*******************************************************
//        BlastReader
//*******************************************************


BlastReader::BlastReader(const string& inf, const string& inFmt):openedGZ(false), processedBatch(false),
lastBlast(NULL), allRead(false), inptFmt(-1){
#ifdef DEBUG
	cerr << "ini blast file\n";
#endif // DEBUG

	if (isGZfile(inf)) {
		openedGZ = true;
#ifdef _gzipread
		blast = new igzstream(inf.c_str(), ios::in);
#else
		cerr << "gzip not supported in your LCA build\n"; exit(50);
#endif
	} else { blast = new ifstream(inf.c_str(), ios::in); }
	if (!*blast) {
		cerr << "Blast input file " << inf << " could not be opened. exiting..\n";
		exit(23);
	}

	if (inFmt == "bl8") {
		inptFmt = 0;
	} else if (inFmt == "uc") {
		inptFmt = 1;
	}

}

list<BlastRes*> BlastReader::getResBatch() {
	string line("");
	if (!processedBatch) {
		getline(*blast, line, '\n');
		lastBlast = new BlastRes(line,inptFmt);
		processedBatch = true;
	}
	list <BlastRes*> ret(0);
	if (lastBlast == NULL || allRead || lastBlast->fail ) { return ret; }
	ret.push_back(lastBlast);
	string cmpQu = lastBlast->Query;
	while (1) {
		if (!getline(*blast, line, '\n')) {
			lastBlast = NULL; allRead = true;  break;
		}
		if (line.length() == 0) { continue; }
		BlastRes* blr = new BlastRes(line, inptFmt);
		if (!blr->isSameQuery(cmpQu)) {
			//abort loop
			lastBlast = blr;
			break;
		}
		ret.push_back(blr);
	}

	return ret;
}
