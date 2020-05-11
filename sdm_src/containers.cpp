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


#include "containers.h"
using namespace std;


void trim(string& str){
	// trim trailing spaces
	size_t endpos = str.find_last_not_of(" \t");
	if( string::npos != endpos )
	{
	    str = str.substr( 0, endpos+1 );
	}

	// trim leading spaces
	size_t startpos = str.find_first_not_of(" \t");
	if( string::npos != startpos )
	{
	    str = str.substr( startpos );
	}
}
//from http://stackoverflow.com/questions/8888748/how-to-check-if-given-c-string-or-char-contains-only-digits
bool is_digits(const std::string &str)
{
	return std::all_of(str.begin(), str.end(), ::isdigit); // C++11
}



bool betterPreSeed(shared_ptr<DNA> d1, shared_ptr<DNA> d2, shared_ptr<DNAunique> ref) {
	//0.2% difference is still ok, but within 0.5% of the best found seed (prevent detoriating sequence match)
	//float blen = (float)ref->length() + (float)d1->length();
	shared_ptr<DNAunique> ref2 = ref->getPair();
	uint curL = d1->mem_length();
	if (d2 != NULL) { curL += d2->mem_length(); }
	else {
		if (d1->has2PrimersDetected() && !ref->has2PrimersDetected()) { return true; }
		if (!d1->has2PrimersDetected() && ref->has2PrimersDetected()) { return false; }
	}
	uint bestL = ref->getBestSeedLength();
	if (d1->getFwdPrimCut() && !ref->getFwdPrimCut()){ return true; }//hard reason
	 if (!d1->getFwdPrimCut() && ref->getFwdPrimCut()) { return false; }

	if (float(curL) / float(bestL) < BestLengthRatio) { return false; }

	//at least 90% length of "good" hit
	if (d1->mem_length() / ref->mem_length() < RefLengthRatio) { return false; }

	//checks if the new DNA has a better overall quality
	//1 added to qual, in case no Qual DNA is used
	float thScore = (1 + d1->getAvgQual())* log((float)d1->mem_length());
	float rScore = (1 + ref->getAvgQual())* log((float)ref->mem_length());
	if (thScore > rScore) {
		//also check for stable lowest score
		if (d1->minQual() > ref->minQual() - MinQualDiff && (d2 == NULL || ref2 == NULL)) { 
			if (curL > bestL) { ref->setBestSeedLength(curL); }
			return true; 
		}
	}
	if (d2 == NULL || ref2 == NULL) {
		return false;
	}
	if (d2->getRevPrimCut() && !ref2->getRevPrimCut()) { return true; }//hard reason
	if (!d2->getRevPrimCut() && ref2->getRevPrimCut()) { return false; }//hard reason

	//at least 90% length of "good" hit
	if (d2->mem_length() / ref2->mem_length() < RefLengthRatio) { return false; }

	//checks if the new DNA has a better overall quality
	//weigh with average id to OTU seed
	thScore += (1 + d2->getAvgQual()) * log((float)d2->mem_length()) * 97;
	rScore += (1 + ref2->getAvgQual()) * log((float)ref2->mem_length()) * 97;
	if (thScore > rScore) {
		//update best seed length score
		if (curL > bestL) { ref->setBestSeedLength(curL); }
		return true;
	}

	return false;
}

dualPrimerDistrStats::dualPrimerDistrStats(const vector<string>&, const vector<string>&){

}

ReadSubset::ReadSubset(const string inf, const string default_outfile):
RemainderStrPos(-1), newHD(0), outFiles(0), outFilesIdx(0) {
	string line;
	ifstream in(inf.c_str());
	if (!in){
		cerr << "Could not find " << inf << " read subset file. Exiting.\n"; exit(90);
	}
	int ini_ColPerRow(0), cnt(0), skips(0);

	//check read subset format
	while (!safeGetline(in, line).eof()) {
		if (line.substr(0, 1) == "#"){ skips++; continue; }
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss, segments, '\t')) {
			ColsPerRow++;
		}

		if (cnt == 0){
			ini_ColPerRow = ColsPerRow;
		}
		else {
			if (ColsPerRow != ini_ColPerRow){
				cerr << "Number of columns in read subset file on line " << cnt + skips << " is " << ColsPerRow << ". Expected " << ini_ColPerRow << " columns.\n";
				exit(91);
			}
			if (cnt > 1000){
				break;
			}
		}
		cnt++;
	}
	if (ini_ColPerRow == 0){
		cerr << "Read Subset File exists, but appears to be badly formated (0 columns detected). Exiting\n"; exit(92);
	}
	if (cnt == 0){
		cerr << "Read Subset File exists, but appears to be badly formated (0 lines detected). Exiting\n"; exit(92);
	}
	in.clear();
	in.seekg(0, ios::beg);
	//extract read subset content
	cnt = 0;
	map<string, int> tmpFiles;
	map<string, int>::iterator tmpFilIT;
	unordered_map <string, int>::iterator TarIT;
	//parameters were set for in matrix, now read line by line
	while (!safeGetline(in, line).eof()) {
		//	while(getline(in,line,'\n')) {
		if (cnt != 0 && line.substr(0, 1) == "#"){ continue; }
		if (line.length() < 5){ continue; }
		stringstream ss; string segments; int cnt2 = 0;
		ss << line;
		while (getline(ss, segments, '\t')) {
			if (cnt2 == 0){
				//free target RD from paired end info
				remove_paired_info(segments);
				TarIT = Targets.find(segments);
				if (TarIT == Targets.end()){
					Targets[segments] = cnt;
				} else {
					break;
				}
			}
			else if (cnt2 == 1){
				newHD.push_back(segments);
			}
			else if (cnt2 == 2){
				uint Idx(0);
				//1: test if outfile already exists
				tmpFilIT = tmpFiles.find(segments);
				if (tmpFilIT != tmpFiles.end()){
					Idx = (*tmpFilIT).second;
				}
				else {//create new key
					Idx = (uint) outFiles.size();
					outFiles.push_back(segments);
					tmpFiles[segments] = Idx;
				}
				//2: add the index to outfile to this position
				outFilesIdx.push_back(Idx);
			}
			cnt2++;
		}
		if (cnt2 == 1) {//no specific header
			newHD.push_back("");
			tmpFilIT = tmpFiles.find("Default");
			uint Idx(0);
			if (tmpFilIT != tmpFiles.end()) {
				Idx = (*tmpFilIT).second;
			} else {//create new key
				Idx = (uint)outFiles.size();
				outFiles.push_back("Default");
				tmpFiles["Default"] = Idx;
			}

			outFilesIdx.push_back(Idx);
		}
		if (cnt2>0) { cnt++; }
	}
	//is extra file info in read subset file?
	if (ini_ColPerRow < 3){
		outFilesIdx.resize(Targets.size(), 0);
		outFiles.resize(1, default_outfile);
	}

}

void ReadSubset::findMatches(shared_ptr<InputStreamer> IS, shared_ptr<MultiDNA> MD, bool mocatFix) {
	shared_ptr<DNA> match(NULL); shared_ptr<DNA> match2(NULL);
	bool cont(true), cont2(true);
	int idx(0);
	int pairs = IS->pairNum();
	unordered_map <string, int>::iterator SEEK;
	bool b_doHD = newHD.size() > 0;
	bool sync(false);//meaningless placeholder
	while (cont) {
		match = IS->getDNA(cont, 0, sync);
		if (match == NULL) {break;}
		if (pairs > 1) {
			match2 = IS->getDNA(cont2, 1, sync);
		}
		string curID = match->getIDPosFree();
		if (mocatFix && curID.length()>5 ) {
			curID = header_stem(curID);
			//cerr << curID << endl;
			//exit(0);
		} 
		//cerr << curID << endl;
		SEEK = Targets.find(curID);
		if (SEEK == Targets.end()) {//no hit
			if (RemainderStrPos != -1) {
				MD->writeSelectiveStream(match, 0, RemainderStrPos);
				if (pairs > 1) {
					MD->writeSelectiveStream(match2, 1, RemainderStrPos);
				}
			} else {// nothing written, but still need to delete
//				delete match; if (pairs > 1) {delete match2;}
			}
			continue;
		} 
		//serious work 
		//cerr << "H";
		idx = SEEK->second;
		if (b_doHD && newHD[idx] != "" && pairs>1) {
			match->newHEad(newHD[idx] + "#1:0");
			match2->newHEad(newHD[idx] + "#2:0");
		}
		MD->writeSelectiveStream(match, 0, outFilesIdx[idx]);
		if (pairs > 1) {
			MD->writeSelectiveStream(match2, 1, outFilesIdx[idx]);
		}

		//finished, clean up to reduce search space (faster??)
		Targets.erase(SEEK);
		if (Targets.size() == 0 && RemainderStrPos==-1) {
			return;
		}
	}
	cerr << Targets.size() << " seqs remaining (not found in current file)\n";
}





MultiDNA::MultiDNA(shared_ptr<Filters> fil,OptContainer& cmdArgs, 
	std::ios_base::openmode writeStatus, shared_ptr<ReadSubset> RDSset, 
	string fileExt,int forceFmt) :
		MFil(fil),subFilter(0),DNAsP1(0),DNAsP2(0),DNAsS1(0),DNAsS2(0),
		DNAsNoHead(0),DNAsP1_alt(0),DNAsP2_alt(0),DNAsS1_alt(0),DNAsS2_alt(0),
		suppressOutWrite(0),write2File(true), mem_used(false),
		DNAinMem(0),writeThreadStatus(0),
		fastQver(fil->getuserReqFastqVer()),
		fastQoutVer(fil->getuserReqFastqOutVer()),BWriteQual(false),
		BWriteFastQ(false), b_multiOutStream(false), b_changeFQheadVer(false),
		b_oneLinerFasta(false), b_doDereplicate(false), b_writePassed(true), b_writeMidPass(true),
		maxReadsPerOFile(fil->maxReadsOutput()),ReadsWritten(fil->writtenReads()),
		maxRdsOut(-1), stopAll(false),
		leadingOutf(""), locCmdArgs(cmdArgs), Derepl(NULL), cntDerep(0), wrMode(ios::out),
		sFile(0), qFile(0), fqFile(0),
		sFileStr(0), qFileStr(0), fqFileStr(0), fqNoBCFile(0), totalFileStrms(0)
/*	qFilePos(0), sFilePos(0), fqFilePos(0),
	qFile2Pos(0), sFile2Pos(0), fqFile2Pos(0),//second pair
	qFileSPos(0), sFileSPos(0), fqFileSPos(0),//singleton
	qFileS2Pos(0), sFileS2Pos(0), fqFileS2Pos(0)//singleton*/
{
	
 	pairedSeq = MFil->isPaired();
	if (cmdArgs.find("-suppressOutput") != cmdArgs.end()) {
		suppressOutWrite = atoi( cmdArgs["-suppressOutput"].c_str() );
	}
	if (suppressOutWrite == 3 || suppressOutWrite == 1){
		b_writePassed = false;
	}
	if (suppressOutWrite == 3 || suppressOutWrite == 2){
		b_writeMidPass = false;
	}
	if (fil->doSubselReads() && RDSset != nullptr ) {
		this->openSeveralOutstreams(locCmdArgs,RDSset, writeStatus);
	} else {//standard one file output stream
		this->openOutStreams(locCmdArgs, MFil->getFileIncrementor(), writeStatus, fileExt, forceFmt);
	}
	if (cmdArgs["-o_fastq_noBC"] != "") {
		openNoBCoutstrean(cmdArgs["-o_fastq_noBC"]);
	}

	//threads = futures(num_threads);
	//if (pairedSeq){cerr<<"paired MFil\n";}
}
MultiDNA::~MultiDNA(){
		//delete MFil;
#ifdef DEBUG
	cerr << "Destr MultiDNA" ;
#endif
	delAllDNAvectors();
#ifdef DEBUG
	cerr << ".. done" << endl;
#endif

/*	for (uint i=0; i<subFilter.size();i++){
		delete subFilter[i];
	}*/
	this->closeOutStreams(true);
#ifdef DEBUG
	cerr << "Subfilters deleted, streams closed" << endl;
#endif
	//delete optim;
}
void MultiDNA::delAllDNAvectors(){
#ifdef DEBUG
	cerr << "cleaning MD..";
#endif
	/*for (unsigned int i=0; i < DNAsP1.size(); i++){ delete DNAsP1[i]; }
	for (unsigned int i=0; i < DNAsP2.size(); i++){ delete DNAsP2[i]; }
	for (unsigned int i=0; i < DNAsS1.size(); i++){ delete DNAsS1[i]; }
	for (unsigned int i=0; i < DNAsS2.size(); i++){ delete DNAsS2[i]; }
	for (unsigned int i=0; i < DNAsNoHead.size(); i++){ delete DNAsNoHead[i]; }
	for (unsigned int i=0; i < DNAsP1_alt.size(); i++){ delete DNAsP1_alt[i]; }
	for (unsigned int i=0; i < DNAsP2_alt.size(); i++){ delete DNAsP2_alt[i]; }
	for (unsigned int i=0; i < DNAsS1_alt.size(); i++){ delete DNAsS1_alt[i]; }
	for (unsigned int i=0; i < DNAsS2_alt.size(); i++){ delete DNAsS2_alt[i]; }
	*/
	DNAsP1.resize(0);DNAsP2.resize(0);DNAsS1.resize(0);
	DNAsS2.resize(0);DNAsNoHead.resize(0);
	DNAsP1_alt.resize(0);DNAsP2_alt.resize(0);
	DNAsS1_alt.resize(0);DNAsS2_alt.resize(0);
#ifdef DEBUG
	cerr << " finished\n";
#endif

}


void MultiDNA::analyzeDNA(shared_ptr<DNA> d, int FilterUse, int pair,int& idx) {
	if (d==NULL){
		return ;
	}
	//collect some info on general run parameters
	MFil->preFilterSeqStat(d, pair);


	if ( !MFil->doFilterAtAll() ) {
		bool isP1 = max(0, pair) == 0;
		if (idx < 0 && !isP1 && !MFil->doubleBarcodes()) {
			;
		} else if (idx < 0) {
			idx = MFil->cutTag(d, isP1); //still need to check for BC
		}
		if (idx >= 0) { //prevent second read pair from being flagged as true
			d->setPassed(true);
		}
		return;
	}

	if (MFil->secondaryOutput() ){
		MFil->checkXtra(d, pair, idx);
	} else if (FilterUse==-1){
		MFil->check(d, false, pair, idx);
	} else {
		cerr << "Invalid control path in analyzeDNA\n"; exit(55);
		subFilter[FilterUse]->check(d, false, pair, idx);
	}
	//count this as failure if BC was present
	//d->prepareWrite(fastQoutVer);
}
vector<bool> MultiDNA::analyzeDNA(shared_ptr<DNA> p1, shared_ptr<DNA> p2, shared_ptr<DNA> mid,
	bool changePHead, int FilterUse){
	cerr << "deprecated analuze DNA"; exit(2323);
	vector<bool> ret(2, true);
	//1st: check if DNA pointer valid
	if (p1 == NULL){
		ret[0] = false;
	} else {
		MFil->preFilterSeqStat(p1, 0);
	}
	if (p2 == NULL){
		ret[1] = false;
	} else {
		MFil->preFilterSeqStat(p2, 1);
	}
	if (mid == NULL){//no MID? kill
		ret[1] = false; ret[0] = false; return ret;
	}
	else {
		mid->setMIDseq(true);
	}

	//2nd: check if basic DNA signatures are valid

	//3rd: check if all quality criteria are met, BC is true etc.
	//if (FilterUse == -1){
	//ret = MFil->check_pairs(p1, p2, mid, ret, changePHead);
	//}	else { //mutithread, to avoid race conditions etc use separate filter
	//	ret = subFilter[FilterUse]->check_pairs(p1, p2, mid, ret, changePHead);
	//}
	//if (ret[0]){ p1->prepareWrite(fastQoutVer); }
	//if (ret[1]){	p2->prepareWrite(fastQoutVer);}

	return ret;
}

bool MultiDNA::checkFastqHeadVersion(shared_ptr<DNA> d,bool disable){
	b_changeFQheadVer = false;
	if (pairedSeq==1){return false;}
	int fastQheadVer = 0;
	string head = d->getID();
	int shouldVer = MFil->FQheadV();
	if (head.find("/1") != string::npos || head.find("/2") != string::npos){
		fastQheadVer = 1;
	}
	if (head.find(" 1:") != string::npos|| head.find(" 2:") != string::npos){
		fastQheadVer = 2;
	}
	if (shouldVer!=0 && shouldVer != fastQheadVer){
		b_changeFQheadVer=true;
	}
	bool ret = b_changeFQheadVer;
	if (disable){
		b_changeFQheadVer = false;
	}
	return ret;
}



void MultiDNA::writeAllStoredDNA(){
#ifdef DEBUG
	printStorage();
	cerr << "Writting stored DNA" << DNAsP1.size() <<" " <<DNAsP2.size()<<endl;
#endif

#ifdef _THREADED
	if (writeThreadStatus>0){
		if (writeThreadStatus>1){wrThread.join();}
		writeThreadStatus++;
		wrThread = thread(&MultiDNA::writeAllStoredDNA2t,this);
	} else {
		writeAllStoredDNA2();
	}
#else
		writeAllStoredDNA2();
#endif
}
void MultiDNA::writeAllStoredDNA2(){
	mem_used=false;
	//
	if (b_writePassed) {

		//
		if (DNAsP1.size() != DNAsP2.size()) {
			for (unsigned int i = 0; i < DNAsP1.size(); i++) {
				writeAndDel(DNAsP1[i], 0);
				ReadsWritten++;
			}
			for (unsigned int i = 0; i < DNAsP2.size(); i++) {
				writeAndDel(DNAsP2[i], 1);
			}
		} else {
			for (unsigned int i = 0; i < DNAsP1.size(); i++) {
				writeAndDel(DNAsP1[i], 0);
				writeAndDel(DNAsP2[i], 1);
				ReadsWritten++;
			}
		}
		for (unsigned int i = 0; i < DNAsS1.size(); i++) {
			writeAndDel(DNAsS1[i], 2);
		}

		for (unsigned int i = 0; i < DNAsS2.size(); i++) {
			writeAndDel(DNAsS2[i], 3);
		}
		DNAsP1.resize(0); DNAsP2.resize(0);
		DNAsS1.resize(0); DNAsS2.resize(0);
	} 
	if (b_writeMidPass) {

		//Xtra file
		if (DNAsP1_alt.size() != DNAsP2_alt.size()) {
			for (unsigned int i = 0; i < DNAsP1_alt.size(); i++) {
				writeAndDel(DNAsP1_alt[i], 0);
				ReadsWritten++;
			}
			for (unsigned int i = 0; i < DNAsP2_alt.size(); i++) {
				writeAndDel(DNAsP2_alt[i], 1);
			}
		} else {
			for (unsigned int i = 0; i < DNAsP1_alt.size(); i++) {
				writeAndDel(DNAsP1_alt[i], 0);
				writeAndDel(DNAsP2_alt[i], 1);
				ReadsWritten++;
			}
		}
		for (unsigned int i = 0; i < DNAsS1_alt.size(); i++) {
			writeAndDel(DNAsS1_alt[i], 2);
		}

		for (unsigned int i = 0; i < DNAsS2_alt.size(); i++) {
			writeAndDel(DNAsS2_alt[i], 3);
		}
		DNAsP1_alt.resize(0); DNAsP2_alt.resize(0);
		DNAsS1_alt.resize(0); DNAsS2_alt.resize(0);
	} 
	//just to be on the safe side..
	delAllDNAvectors();
#ifdef DEBUG
	cerr <<  " .. Finished" << endl;
#endif
}
#ifdef _THREADED
void MultiDNA::writeAllStoredDNA2t(){
    std::lock_guard<std::mutex> guard(mutex);
	write2File=true;mem_used=false;
	if (DNAsP1.size() != DNAsP2.size()){
		for (unsigned int i=0; i<DNAsP1.size();i++){
			writeAndDel(DNAsP1[i],0);
		}
		for (unsigned int i=0; i<DNAsP2.size();i++){
			writeAndDel(DNAsP2[i],1);
		}
	}else {
		for (unsigned int i=0; i<DNAsP1.size();i++){
			writeAndDel(DNAsP1[i],0);
			writeAndDel(DNAsP2[i],1);
		}
	}
	for (unsigned int i=0; i<DNAsS1.size();i++){
		writeAndDel(DNAsS1[i],2);
	}

	for (unsigned int i=0; i<DNAsS2.size();i++){
		writeAndDel(DNAsS2[i],3);
	}
	DNAsP1.resize(0);DNAsP2.resize(0);
	DNAsS1.resize(0);DNAsS2.resize(0);
}
#endif
void MultiDNA::incrementOutputFile(){
	MFil->incrementFileIncrementor();
	this->closeOutStreams(true);
	//this is definetely a new file
	
	this->openOutStreams(locCmdArgs,MFil->getFileIncrementor(),ios_base::out);
	ReadsWritten = 0;
}
void Filters::addDNAtoCStats(shared_ptr<DNA> d,int Pair) {
	//here should be the only place to count Barcodes!
	int easyPair = Pair < 3 ? Pair - 1 : Pair - 3;
	colStats[easyPair].total2++;
	if (d->isPassed() || d->isMidQual()) {
		this->DNAstatLQ(d, easyPair, d->isMidQual());
		colStats[easyPair].totalSuccess++;
	} else {
		colStats[easyPair].totalRejected++;

	}

	//some general stats that always apply:
	if (d->QualCtrl.PrimerFail) {
		colStats[easyPair].PrimerFail++;
	}
	if (d->QualCtrl.PrimerRevFail) {
		colStats[easyPair].PrimerRevFail++;
	}
	if (d->QualCtrl.minLqualTrim) {
		colStats[easyPair].minLqualTrim++;
	}
	if (d->QualCtrl.TagFail) {
		colStats[easyPair].TagFail++;
	}
	if (d->QualCtrl.fail_correct_BC) {
		colStats[easyPair].fail_correct_BC++;
	}
	if (d->QualCtrl.suc_correct_BC) {
		colStats[easyPair].suc_correct_BC++;
	}
	if (d->QualCtrl.RevPrimFound) {
		colStats[easyPair].RevPrimFound++;
	}
	if (d->QualCtrl.QWinTrimmed || d->QualCtrl.AccErrTrimmed) {
		colStats[easyPair].Trimmed++;
	}
	if (d->getTA_cut()) {
		colStats[easyPair].adapterRem++;
	}


	if (d->isPassed() || d->isMidQual()) {
		countBCdetected(d->getBCnumber(), easyPair, false);
		//and register as success
	} else {
		if (d->getBarcodeDetected()) {
			//DNA is no longer useful
			failedStats2(d, easyPair);
		}
		//delete d; 
		if (d->QualCtrl.AvgQual) {
			colStats[easyPair].AvgQual++;
		}
		if (d->QualCtrl.minL) {
			colStats[easyPair].minL++;
		}
		if (d->QualCtrl.maxL) {
			colStats[easyPair].maxL++;
		}
		if (d->QualCtrl.HomoNT) {
			colStats[easyPair].HomoNT++;
		}
		if (d->QualCtrl.MaxAmb) {
			colStats[easyPair].MaxAmb++;
		}
		if (d->QualCtrl.BinomialErr) {
			colStats[easyPair].BinomialErr++;
		}
		if (d->QualCtrl.QualWin) {
			colStats[easyPair].QualWin++;
		}
	}
	

}
bool  MultiDNA::saveForWrite(shared_ptr<DNA> d,int Pair){
	//second most important part: collect stats on DNA passing through here (should be all read)
	//most important part: save DNA to be written later (or discard)
	if (d == NULL || stopAll) {
		return !stopAll;
	}
	//int easyPair = Pair < 3 ? Pair - 1 : Pair - 3;
	MFil->addDNAtoCStats(d, Pair);

	if( d->isPassed()){

		if (b_writePassed){
			d->prepareWrite(fastQoutVer);
			//lock MultiDNA
#ifdef _THREADED
			std::lock_guard<std::mutex> guard(mutex);
#endif
			//dereplicate & create copy of DNA?

			mem_used = true;
			if (Pair == 1){
				DNAsP1.push_back(d);
			}
			else if (Pair == 2){
				DNAsP2.push_back(d);
			}
			else if (Pair == 3){
				DNAsS1.push_back(d);
				MFil->colStats[0].singleton++;
			}
			else if (Pair == 4){
				DNAsS2.push_back(d);
				MFil->colStats[1].singleton++;
			}
			DNAinMem++;
		}
		
	} else if (d->isMidQual()){

		if (b_writeMidPass){
			d->prepareWrite(fastQoutVer);
			mem_used = true;
			if (Pair == 1){
				DNAsP1_alt.push_back(d);
			}
			else if (Pair == 2){
				DNAsP2_alt.push_back(d);
			}
			else if (Pair == 3){
				MFil->statAddition.singleton++;
				DNAsS1_alt.push_back(d);
			}
			else if (Pair == 4){
				MFil->statAddition.singleton++;
				DNAsS2_alt.push_back(d);
			}
			DNAinMem++;
		}

	} 
	//automatic mechanism to write to File, once enough DNA is in memory
	if (write2File && DNAinMem>DNAinMemory){
		writeAllStoredDNA();
		DNAinMem=0;
	}
	if (maxReadsPerOFile>0 && ReadsWritten+DNAinMem >= maxReadsPerOFile){
		//cerr << "ReadsWritten " << ReadsWritten << " DNAinMem " << DNAinMem << endl;
		DNAinMem=0;
		incrementOutputFile();
	}
	if (maxRdsOut > 0 && ReadsWritten + DNAinMem >= maxRdsOut) {
		writeAllStoredDNA();
		stopAll = true;
		
	}
	return !stopAll;
}
void MultiDNA::writeAndDel(shared_ptr<DNA> d,int Pair){
	//ofstream tmpS, tmpQ, tmpFQ;
//	int PairC = Pair;
//	if (Pair > 1) {
//		PairC = Pair - 2;
//	}

	
	int ClsVec = -1;
	if (d != NULL) {
		if (MFil->doFilterAtAll()) {
			if (d->isPassed()) {
				ClsVec = 0;
			} else if (d->isMidQual()) {
				ClsVec = 1;
			}
		}
		else {
			ClsVec = 0;
		}
	}
	if (ClsVec >= 0) {
		if (BWriteFastQ && b_changeFQheadVer) {//check if header PE naming needs to be changed
			d->changeHeadPver(MFil->FQheadV());
		}
		if (BWriteFastQ) {
			if (fqFile[ClsVec][Pair ] == NULL ) {
				//cerr << "_";
				openOFstreamFQ(fqFileStr[ClsVec][Pair], wrMode, ClsVec, Pair , "Appending");
			}
			//cerr << "X";
			d->writeFastQ(*(fqFile[ClsVec][Pair]));
		} else {
			if (sFile[ClsVec][Pair]==NULL ) {
				openOFstreamFNA(sFileStr[ClsVec][Pair], wrMode, ClsVec, Pair, "Appending");
			}

			d->writeSeq(*sFile[ClsVec][Pair ], b_oneLinerFasta);
			if (BWriteQual) { 
				if (qFile[ClsVec][Pair]==NULL) {
					openOFstreamQL(qFileStr[ClsVec][Pair], wrMode, ClsVec, Pair, "Appending");
				}
				d->writeQual(*qFile[ClsVec][Pair], b_oneLinerFasta);
			}
		}
	}

//	delete d;

}
void MultiDNA::writeSelectiveStream(shared_ptr<DNA> d, int Pair,int FS) {
	//ofstream tmpS, tmpQ, tmpFQ;
	if (d == 0) {
		return;
	}
	d->prepareWrite(fastQoutVer);
	int PairC = Pair;
	if (Pair > 1) {
		PairC = Pair - 2;
	}
	if  (d->isPassed()) {
		MFil->DNAstatLQ(d, PairC, false);
	} else if ( d->isMidQual()) {
		MFil->DNAstatLQ(d, PairC, true);
	}
	if (BWriteFastQ && b_changeFQheadVer) {//check if header PE naming needs to be changed
		d->changeHeadPver(MFil->FQheadV());
	}
	if (BWriteFastQ) {
		if (fqFile[FS][Pair]==NULL||!*fqFile[FS][Pair]) {
			openOFstreamFQ(fqFileStr[FS][Pair], ios_base::app, FS, Pair, "Appending");
		}
		d->writeFastQ(*(fqFile[FS][Pair]));
		delete fqFile[FS][Pair]; fqFile[FS][Pair] = NULL;
	} else {
		if (sFile[FS][Pair]==NULL||!*sFile[FS][Pair]) {
			openOFstreamFNA(sFileStr[FS][Pair], ios_base::app, FS, Pair, "Appending");
		}
		d->writeSeq(*sFile[FS][Pair],b_oneLinerFasta);
		if (BWriteQual) { 
			if (qFile[FS][Pair]==NULL||!*qFile[FS][Pair]) {
				openOFstreamQL(sFileStr[FS][Pair], ios_base::app, FS, Pair, "Appending");
			}
			d->writeQual(*qFile[FS][Pair], b_oneLinerFasta);
		}
	}


//	delete d;

}
void MultiDNA::writeNonBCReads(shared_ptr<DNA> d, shared_ptr<DNA> d2) {
	if (fqNoBCFile.size() == 2) {
		if (!d->getBarcodeDetected() && !d2->getBarcodeDetected()) {
			if ((!d->getBarcodeDetected() && d2->getBarcodeDetected()) ||
				(d->getBarcodeDetected() && !d2->getBarcodeDetected())) {
				cerr << "Barcode only set in 1 reads.. something wrong!\n";
			}
			d->writeFastQ(*(fqNoBCFile[0]));
			d2->writeFastQ(*(fqNoBCFile[1]));
		}
	}
}

void MultiDNA::setSubfilters(int num){
	if (num<2){return;}
	subFilter.resize(num);
	for (uint i=0;i<subFilter.size();i++){
		subFilter[i].reset(new Filters(MFil,true));
	}
}
void MultiDNA::mergeSubFilters(){
	vector<int> idx (MFil->Barcode.size(),0);
	for (uint i=0;i<idx.size();i++){idx[i]=i;}
	for (uint i=0; i<subFilter.size();i++){
		MFil->addStats(subFilter[i],idx);
	}
}
void MultiDNA::attachDereplicator(shared_ptr<Dereplicate> de) {
	if (de != NULL) { 
		Derepl = de; b_doDereplicate = true; 
		Derepl->setPaired(pairedSeq>1); 
		//insert code here to fix BC offset in filter & add to derep info on sample names from the current filter
		int curBCOffset = Derepl->getHighestBCoffset();
#ifdef DEBUG
		cerr << "BARCODE INCREMENT: " << curBCOffset << endl;
#endif
		
		MFil->setBCoffset(curBCOffset);
		Derepl->BCnamesAdding(MFil);
	}
}

/*void MultiDNA::depPrep(shared_ptr<DNA> tdn) {

	if (b_doDereplicate && tdn->isPassed()) {
		cntDerep++;
		Derepl->addDNA(tdn);
	}
}*/
void MultiDNA::depPrep(shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2) {
	bool added = false;
	if (!b_doDereplicate) {
		return;
	}
	Derepl->addDNA(tdn, tdn2, added);
	if (added){ 
		cntDerep++; 
		if (tdn->getBarcodeDetected() && !tdn->isPassed() && !tdn->isMidQual() ){
			MFil->statAddDerepBadSeq(tdn->getBCnumber());
					
		}
	}
	
}

void MultiDNA::resetOutFilesAndFilter(){
#ifdef DEBUG
	cerr << "resetOutFilesAndFilter" << endl;
#endif
	//reset count stats
	MFil->resetStats();
	//reset all stored DNA
	delAllDNAvectors();
	Derepl->reset();
	//close streams -- why? no reason, since nothing has been writen
	//this->closeOutStreams();
	totalFileStrms = 0;
}

void MultiDNA::closeOutStreams(bool wr){


	write2File = true;
	if (wr){
		this->writeAllStoredDNA();
	}

#ifdef DEBUG
	cerr << "closing output streams";
#endif

	for (int i = 0; i < (int)sFile.size(); i++) {
		for (int j = 0; j < (int) sFile[i].size(); j++) {
			if (sFileStr[i][j] != "T") {
				delete sFile[i][j]; sFile[i][j] = NULL;
			}
		}
	}
	for (int i = 0; i < (int)qFile.size(); i++) {
		for (int j = 0; j < (int)qFile[i].size(); j++) {
			if (qFileStr[i][j] != "T") {
				delete qFile[i][j]; qFile[i][j] = NULL;
			}
		}
	}
	for (int i = 0; i < (int)fqFile.size(); i++) {
		for (int j = 0; j < (int)fqFile[i].size(); j++) {
			if (fqFileStr[i][j] != "T") {
				delete fqFile[i][j]; fqFile[i][j] = NULL;
			}
		}
	}
	//if(qFile){qFile.close();}if(sFile){sFile.close();}if(fqFile){fqFile.close(); }
	//other housekeeping tasks
	this->mergeSubFilters();
	MFil->setWrittenReads(ReadsWritten);
#ifdef DEBUG
	cerr << ".. closed\n";
#endif

}
/*void MultiDNA::resetOutStreams(){
	if(qFile){qFile.seekp(qFilePos);}if(sFile){sFile.seekp(sFilePos);}if(fqFile){fqFile.seekp(fqFilePos); }
	if(qFile2){qFile2.seekp(qFile2Pos);}if(sFile2){sFile2.seekp(sFile2Pos);}if(fqFile2){fqFile2.seekp(fqFile2Pos); }
	if(qFileS){qFileS.seekp(qFileSPos);}if(sFileS){sFileS.seekp(sFileSPos);}if(fqFileS){fqFileS.seekp(fqFileSPos); }
	if(qFileS2){qFileS2.seekp(qFileS2Pos);}if(sFileS2){sFileS2.seekp(sFileS2Pos);}if(fqFileS2){fqFileS2.seekp(fqFileS2Pos); }
}*/
void MultiDNA::openOFstream(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep, int wh) {
	switch (wh) {
	case 1:
		openOFstreamFNA(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	case 0:
		openOFstreamFQ(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	case 2:
		openOFstreamQL(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	default:
		cerr << "Wrong wh specified"; exit(1002);
	}
}
void MultiDNA::openNoBCoutstrean(const string inS) {
	vector<string> tfnaout = splitByCommas(inS);
	fqNoBCFile.resize(2, NULL);
	fqNoBCFile[0] = new ofstream(tfnaout[0], wrMode);
	fqNoBCFile[1] = new ofstream(tfnaout[1], wrMode);

}

void MultiDNA::openOFstreamFQ(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep) {
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1+1 >= (int)fqFile.size()) {
		vector<ostream*> nullVec(4, NULL);
		fqFile.resize(p1+1, nullVec);
	}
	if ((int)fqFileStr.size() - 1 <= p1) {
		fqFileStr.resize(p1 + 1, vector<string>(4, ""));
	}
	fqFileStr[p1][p2] = opOF;
	if (onlyPrep) { return; }
	if (p1 == 1 && !b_writeMidPass ){ return; }//p1==1: mid passed suppressOutWrite >= 2
	if (p1 == 0 && !b_writePassed){ return; }//suppressOutWrite == 1
	if (opOF == "T") {
		fqFile[p1][p2] = &std::cout;
	} else if (isGZfile(opOF)) {
#ifdef DEBUG
		cerr << "open fq.gz with wrmode" << wrMode << endl;
#endif
#ifdef _gzipread
		fqFile[p1][p2] = new ogzstream(opOF.c_str(), wrMode);// wrMode);
#else
		cerr << "gzip outpout not supported in your sdm build\n" << opOF; exit(50);
#endif

	}else{
		fqFile[p1][p2] = new ofstream(opOF, wrMode);
	}

	if (!*fqFile[p1][p2]) {
		cerr << "Could not open " << errMsg << " fastq output file " << opOF << endl << p1 << " " << p2 << " " << totalFileStrms << endl;
		exit(4);
	}
}
void MultiDNA::openOFstreamFNA(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep) {
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1+1 >= (int)sFile.size()) {
		vector<ostream*> nullVec(4, NULL);
		sFile.resize(p1+1, nullVec);
	}
	if ((int)sFileStr.size() - 1 <= p1) {
		sFileStr.resize(p1 + 1, vector<string>(4, ""));
	}
	sFileStr[p1][p2] = opOF;
	if (onlyPrep ) { return; }
	if (p1 == 1 && !b_writeMidPass){ return; }//p1==1: mid passed suppressOutWrite >= 2
	if (p1 == 0 && !b_writePassed){ return; }//suppressOutWrite == 1
	if (opOF == "T") {
		sFile[p1][p2] = &std::cout;
	} else if (isGZfile(opOF)) {
#ifdef _gzipread
		sFile[p1][p2] = new ogzstream(opOF.c_str(), wrMode);
#else
		cerr << "gzip outpout not supported in your sdm build\n" << opOF; exit(50);
#endif
	}
	else {
		sFile[p1][p2] = new ofstream(opOF, wrMode);
	}
	if (!*sFile[p1][p2]) {
		cerr << "Could not open " << errMsg << " fasta output file " << opOF << endl;
		exit(4);
	}
}
void MultiDNA::openOFstreamQL(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep) {
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1+1 >= (int)qFile.size()) {
		vector<ostream*> nullVec(4, NULL);
		qFile.resize(p1+1, nullVec);
	}
	if ((int)qFileStr.size() - 1 <= p1) {
		qFileStr.resize(p1 + 1, vector<string>(4, ""));
	}
	qFileStr[p1][p2] = opOF;
	if (onlyPrep) { return; }
	if (p1 == 1 && !b_writeMidPass){ return; }//p1==1: mid passed suppressOutWrite >= 2
	if (p1 == 0 && !b_writePassed){ return; }//suppressOutWrite == 1
	if (opOF == "T") {
		qFile[p1][p2] = &std::cout;
	}
	else if (isGZfile(opOF)) {
#ifdef _gzipread
		qFile[p1][p2] = new ogzstream(opOF.c_str(), wrMode);
#else
		cerr << "gzip outpout not supported in your sdm build\n" << opOF; exit(50);
#endif
	} else {
		qFile[p1][p2] = new ofstream(opOF, wrMode);
	}
	if (!*qFile[p1][p2]) {
		cerr << "Could not open " << errMsg << " quality output file " << opOF << endl;
		exit(4);
	}
}

void MultiDNA::openSeveralOutstreams(OptContainer& cmdArgs, shared_ptr<ReadSubset> RDS, std::ios_base::openmode wrMode) {
#ifdef DEBUG
	cerr << " openining multiple out streams" << endl;
#endif
	string path = "", fileEnd(".fna");
	vector<string> outFile = RDS->getOFiles();
	bool openStrms = true; int omode(1);
	vector<vector<string>>& tmp = sFileStr;
	if (cmdArgs.find("-o_fastq") != cmdArgs.end() && cmdArgs["-o_fastq"] != "" && cmdArgs["-o_fna"] == "") { //write fastq
		BWriteFastQ = true;
		path = cmdArgs["-o_fastq"];
		tmp = fqFileStr; omode = 0;
		fileEnd = ".fastq";
	} else  if (cmdArgs["-o_fna"] != "") {
		BWriteFastQ = false;
		path = cmdArgs["-o_fna"];
	} else {
		cerr << "Could not find valid sequence output file path. Exiting\n";
		exit(55);
	}
	if (RDS->multiFile()) {
		b_multiOutStream = true;
	}
	int i(0);
	for (i = 0; i < (int)outFile.size(); i++) {
		if (totalFileStrms >= maxFileStreams && openStrms) {
			cerr << "Too many output file streams\nSwitching to dynamical file appending\n";
			openStrms = false;
			//this->closeOutStreams();
		}
		string baseFile = path + removeFileEnding(outFile[i]);
		if (pairedSeq > 1) {
			if (i < 4) {
				cerr << "Outfile " << i << ": "<<baseFile + ".1" << fileEnd << "\n";
			}
			openOFstream(baseFile + ".1" + fileEnd, wrMode, i, 0, "paired 1st " + itos(i), !openStrms, omode);
			//2nd pair
			openOFstream(baseFile + ".2" + fileEnd, wrMode, i, 1, "paired 2nd " + itos(i), !openStrms, omode);
			//1st singleton
			//if (openStrms) { openOFstream(baseFile + ".1" + fileEnd + SingletonFileDescr, wrMode, i, 2, "Singleton 1 " + itos(i), omode); }
			//2nd singleton
			//if (openStrms) { openOFstream(baseFile + ".2" + fileEnd + SingletonFileDescr, wrMode, i, 3, "Singleton 2 " + itos(i), omode); }
			totalFileStrms += 2;
		} else {
			openOFstream(baseFile + fileEnd, wrMode, i, 0, "main " + itos(i), !openStrms, omode);
			totalFileStrms++;
		}

	}
	//pipe remaining reads into this file
	if (cmdArgs["-excludeFile"] != "") {
		//set a flag
		RDS->setRemainingFilepipe(i);
		if (pairedSeq == 1) {
			openOFstream(cmdArgs["-excludeFile"], wrMode, i, 0, "excludeFile file ", !openStrms, omode);
		} else {
			string baseFile = path + removeFileEnding(cmdArgs["-excludeFile"]);
			if (i < 4) {
				cerr << "out file " << i << baseFile + ".1" << fileEnd << "\n";
			}
			openOFstream(baseFile + ".1" + fileEnd, wrMode, i, 0, "paired 1st " + itos(i), !openStrms, omode);
			//2nd pair
			openOFstream(baseFile + ".2" + fileEnd, wrMode, i, 1, "paired 2nd " + itos(i), !openStrms, omode);
			//1st singleton
			//if (openStrms) { openOFstream(baseFile + ".1" + fileEnd + SingletonFileDescr, wrMode, i, 2, "Singleton 1 " + itos(i), omode); }
			//2nd singleton
			//if (openStrms) { openOFstream(baseFile + ".2" + fileEnd + SingletonFileDescr, wrMode, i, 3, "Singleton 2 " + itos(i), omode); }
			totalFileStrms += 2;
		}
	}
}

void MultiDNA::openOutStreams(OptContainer& cmdArgs,int fileIt,std::ios_base::openmode wrMode_i, 
	string fileExt, int forceFmt){
	this->setwriteMode(wrMode_i);
	if ( suppressOutWrite == 3 || (cmdArgs["-o_fastq"] == "" && cmdArgs["-o_fna"] == "" && cmdArgs["-o_qual"] == "") ){
		suppressOutWrite = 3; b_writePassed = false; b_writeMidPass = false; return;
	}
	if (forceFmt != -1){
		if (forceFmt == 1 && (cmdArgs.find("-o_fna") == cmdArgs.end() || cmdArgs["-o_fna"] == "") ){//force fna, no qual, required for seed ref fastas
			cmdArgs["-o_fna"] = cmdArgs["-o_fastq"];
			cmdArgs["-o_fastq"] = "";
		}
	}
	

	if (cmdArgs.find("-o_fastq")  != cmdArgs.end() && cmdArgs["-o_fastq"] != "" && cmdArgs["-o_fna"] == ""){ //write fastq
		this->setFastQWrite(true);
		if (pairedSeq>1){ //open second pair + singleton
#ifdef DEBUG
			cerr << " paired fastq out " << endl;
#endif
			vector<string> tfnaout = splitByCommas(cmdArgs["-o_fastq"]);
			if (tfnaout.size()!=2){
				cerr<<"Paired sequences given as input, requires paired output file (2 files separated by \",\"). Given output file = "<<cmdArgs["-o_fastq"]<<endl; exit(57);
			}
			leadingOutf = applyFileIT(tfnaout[0] + fileExt, fileIt);
			openOFstreamFQ(leadingOutf.c_str(), wrMode, 0, 0, "paired 1st");
			//2nd pair
			openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
			//1st singleton
			openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 0, 2, "Singleton 1", true);
			//2nd singleton
			openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 0, 3, "Singleton 2", true);
			//additional file
			if (cmdArgs.find("-o_fastq2")  != cmdArgs.end() && cmdArgs["-o_fastq2"].length()>1){ //write fastq
				tfnaout = splitByCommas(cmdArgs["-o_fastq2"]);
				openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st");
				openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd");
				openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt , fileIt, SingletonFileDescr).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
				openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt , fileIt, SingletonFileDescr).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
			}
		} else {
#ifdef DEBUG
			cerr << " single fastq out " << endl;
#endif
			openOFstreamFQ(applyFileIT(cmdArgs["-o_fastq"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "the main");
			leadingOutf = applyFileIT(cmdArgs["-o_fastq"] + fileExt, fileIt);
			//additional file
			if (cmdArgs.find("-o_fastq2") != cmdArgs.end() && cmdArgs["-o_fastq2"].length()>1) { //write fastq
				openOFstreamFQ(applyFileIT(cmdArgs["-o_fastq2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "the additional");
			}
		}

		return;
	}
	BWriteFastQ=false;
	if (pairedSeq==1){
#ifdef DEBUG
		cerr << " fasta singleton output " << endl;
#endif
		openOFstreamFNA(applyFileIT(cmdArgs["-o_fna"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "main");
		leadingOutf = applyFileIT(cmdArgs["-o_fna"] + fileExt, fileIt);
		if (cmdArgs["-o_qual"] != ""){
			openOFstreamQL(applyFileIT(cmdArgs["-o_qual"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "main");
			this->setQualWrite(true);
		} else {
			this->setQualWrite(false);
		}
		//additional file (secondary filter)
		if (cmdArgs.find("-o_fna2") != cmdArgs.end() && cmdArgs["-o_fna2"].length()>1) { //add file
			openOFstreamFNA(applyFileIT(cmdArgs["-o_fna2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional");
		}
		if (cmdArgs.find("-o_qual2") != cmdArgs.end() && cmdArgs["-o_qual2"].length()>1) { //add file
			openOFstreamQL(applyFileIT(cmdArgs["-o_qual2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional");
			this->setQualWrite(true);
		}

	} else {
#ifdef DEBUG
		cerr << " fasta paired output " << endl;
#endif

		vector<string> tfnaout = splitByCommas(cmdArgs["-o_fna"]);
		if (tfnaout.size()!=2){
			cerr<<"Paired sequences given as input, requires paired output file (2 files separated by \",\"). Given output file = "<<cmdArgs["-o_fna"]<<endl; exit(57);
		}
		
		openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 0, 0, "paired 1st");
		leadingOutf = applyFileIT(tfnaout[0] + fileExt, fileIt);

		openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
		openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 2, "Singleton 1", true);
		openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 3, "Singleton 2", true);
		/*setFilePos(sFile,sFilePos);
		setFilePos(sFile2,sFile2Pos);
		setFilePos(sFileS,sFileSPos);
		setFilePos(sFileS2,sFileS2Pos);*/
		if (cmdArgs.find("-o_fna2")  != cmdArgs.end() && cmdArgs["-o_fna2"].length()>1){ //write fastq
			tfnaout = splitByCommas(cmdArgs["-o_fna2"]);

			openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st", true);
			openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd", true);
			openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
			openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
		}
		if (cmdArgs["-o_qual"] != ""){
			this->setQualWrite(true);
			vector<string> tqout = splitByComma(cmdArgs["-o_qual"],true);
			openOFstreamQL(applyFileIT(tqout[0] + fileExt, fileIt).c_str(), wrMode, 0, 0, "paired 1st");
			openOFstreamQL(applyFileIT(tqout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
			openOFstreamQL(applyFileIT(tqout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 2, "Singleton 1", true);
			openOFstreamQL(applyFileIT(tqout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 3, "Singleton 2", true);
			if (cmdArgs["-o_qual2"] != "") {
				vector<string> tqout = splitByComma(cmdArgs["-o_qual2"], true);
				openOFstreamQL(applyFileIT(tqout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st", true);
				openOFstreamQL(applyFileIT(tqout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd", true);
				openOFstreamQL(applyFileIT(tqout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
				openOFstreamQL(applyFileIT(tqout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
			} 
		} else {
			this->setQualWrite(false);
		}

	}


}
//*******************************************
//*        DEREPLICATE OBJECT
//*******************************************
Dereplicate::Dereplicate(OptContainer& cmdArgs):
BCN2SmplID(0), b_usearch_fmt(true), b_singleLine(true), b_pairedInput(false), 
minCopies(1,0), minCopiesStr("0"), //default minCopies accepts every derep 
totSize(0), tmpCnt(0), curBCoffset(0){
	outfile = cmdArgs["-o_dereplicate"];
	if (cmdArgs.find("-dere_size_fmt") != cmdArgs.end() && cmdArgs["-dere_size_fmt"] == "1") {
		b_usearch_fmt = false;
	}
	if (cmdArgs.find("-min_derep_copies") != cmdArgs.end()) {
		minCopies[0] = -1;//in this case reset to -1 the first entry..
		minCopiesStr = cmdArgs["-min_derep_copies"];
		string x = cmdArgs["-min_derep_copies"];

		vector<string> xs = splitByCommas(x, ',');
		for (size_t i = 0; i < xs.size(); i++) {
			vector<string> tmp = splitByComma(xs[i], false, ':');
			if (tmp.size() > 2) {
				cerr << "Derep string " << xs[i] << " has to be two integers seqparated by \":\"\n"; exit(623);
			}
			int pos = 1; int cnt = -1;
			if (tmp.size() == 1) {//this is the 1 sample position, if not set
				cnt = atoi(tmp[0].c_str());
			}else{
				pos = atoi(tmp[1].c_str());//number of samples
				cnt = atoi(tmp[0].c_str());//16s copy numbers
			}
			if (pos <= 0) { cerr << "wrong derplicate position (<=0):" << xs[i] << endl; exit(312); }
			if (pos > 3000) { cerr << "too large derplicate position (>3000):" << xs[i] << endl; exit(313); }
			if (minCopies.size() < (size_t)pos) {
				minCopies.resize(pos, -1);
			}
			minCopies[pos-1] = cnt;
		}
		minCopiesSiz = minCopies.size();
		//minCopies = atoi(cmdArgs["-min_derep_copies"].c_str());

	}
}
/*bool Dereplicate::addDNA(shared_ptr<DNA> d) {
	//1st build hash of DNA
	string seq = d->getSeqPseudo();
	int BCN = d->getBCnumber();
	if (BCN >= 0) { tmpCnt++; }
	HashDNAIT spotted = Tracker.find(seq);
	if (spotted != Tracker.end()) {// found something
		int hpos ( spotted->second );
		Dnas[hpos]->addSmpl(BCN);
		if (betterPreSeed(d, NULL, Dnas[hpos])) {
			//replace old DNA
			DNAunique *du = new DNAunique(d, -1);
			du->saveMem();
			du->setBestSeedLength(Dnas[hpos]->getBestSeedLength());
			du->transferOccurence(Dnas[hpos]);
			delete Dnas[hpos];
			Dnas[hpos] = du;
		}
		return false;
	} else {
		//create entry
		Tracker[seq] = (int)Dnas.size();
		DNAunique *tmp = new DNAunique(d, BCN);
		tmp->saveMem();
		Dnas.push_back(tmp);
		return true;
	}
	return true;
}*/

bool Dereplicate::addDNA( shared_ptr<DNA> d, shared_ptr<DNA> d2, bool& added) {
	//1st build hash of DNA
	if (!d->getBarcodeDetected()) {
		return false;
	}
	string seq = d->getSeqPseudo();
	int BCN = d->getBCnumber();
	bool pass = d->isPassed();
	shared_ptr<DNAunique> tmp = make_shared<DNAunique>(d, BCN);

	//if (!) { return false; }
	auto spotted = Tracker.find(tmp);
	if (spotted != Tracker.end()) {// found something
		added = true; //d->setMidQual(false); tmp->setMidQual(false);
		//int idx (spotted->second);
		(*spotted)->addSmpl(BCN);
		if (pass && betterPreSeed(d, d2, (*spotted))) {
			//replace old DNA
			//shared_ptr<DNAunique> du = make_shared<DNAunique>(d, -1);
			tmp->saveMem();
			tmp->setBestSeedLength((*spotted)->getBestSeedLength());
			tmp->transferOccurence((*spotted));
			if (d2 != NULL) {
				tmp->attachPair( make_shared<DNAunique>(d2, -1));
			} else {//do nothing
				//du->attachPair(new DNAunique());
			}
			//delete Dnas[idx];
			Tracker.erase(*spotted);
			Tracker.insert(tmp);
		}
		return false;//added to seeds
	} else if (pass){
		added = true; d->setMidQual(false);
		//create entry
		//Tracker[seq] = (int)Dnas.size();
		//shared_ptr<DNAunique> tmp = make_shared<DNAunique>(d, BCN);
		tmp->saveMem();
		if (d2 != NULL) {
			tmp->attachPair( make_shared<DNAunique>(d2, BCN));
		} else {
			//tmp->attachPair(new DNAunique());
		}
		Tracker.insert(tmp);
		//Dnas.push_back(tmp);
		//size_t idx = Dnas.size();
		return true;//new seed
	}
	
	return true;//didn't do a thing..
}
void Dereplicate::BCnamesAdding(shared_ptr<Filters> fil) {
	vector<string> refSID = fil->SampleID;
	for (size_t i = 0; i < refSID.size(); i++) {
		BCN2SmplID.push_back(refSID[i]);
	}
	//this will in the end be used to synchronize the different barcodes between samples
	curBCoffset = (int)BCN2SmplID.size();
	//cerr << curBCoffset << " Added BC list\n";
}
void Dereplicate::reset() {
//	for (size_t i = 0; i < Dnas.size(); i++) { delete Dnas[i]; }
//	Dnas.resize(0);
	totSize = 0; tmpCnt = 0;
	//BCN2SmplID.resize(0);
	Tracker.clear();
}
bool DNAuPointerCompare(shared_ptr<DNAunique> l, shared_ptr<DNAunique> r) { return l->getCount() > r->getCount(); }

string Dereplicate::writeDereplDNA(shared_ptr<Filters> mf) {
	ofstream of, omaps, of2, ofRest, of2p2;
	cerr << "Evaluating and writing dereplicated reads..\n";
	int fastqVer = mf->getuserReqFastqOutVer();
	string mapF = outfile.substr(0, outfile.find_last_of('.')) + ".map";
	string outHQf = outfile.substr(0, outfile.find_last_of('.')) + ".hq.fq";
	string outHQf_p2 = outfile.substr(0, outfile.find_last_of('.')) + ".2.hq.fq";
	string outRest = outfile + ".rest";
	//setup map header
	omaps.open(mapF.c_str(), ios::out);//| ios::binary
	of2.open(outHQf.c_str(), ios::out);
	of.open(outfile.c_str(), ios::out);
	ofRest.open(outRest.c_str(), ios::out);
	if (b_pairedInput) {
		of2p2.open(outHQf_p2, ios::out);
	}

	//sample specific derep filter
	bool bDerepSmpSpcfc(false);
	vector<int> derepSmpSpcfc = mf->getDrerepSampleSpecifity();
	if (derepSmpSpcfc.size() > 0) { bDerepSmpSpcfc = true; }

	//combiner relevant vars
	bool bCombiSmpl = mf->combineSamples();
	vector<int> smplId2comb(0,0);
	unordered_map<string, int> & combiMapCollectGrp = mf->combiMapCollectGrp;

	//vector<string> refSID = mf->SampleID;
	omaps << "#SMPLS";
	if (!bCombiSmpl){
		for (size_t i = 0; i < BCN2SmplID.size(); i++) {
			omaps << "\t" + itos((int)i) + ":" + BCN2SmplID[i];
		}
	} else {
		smplId2comb = mf->combiSmplConvergeVec(BCN2SmplID);
		if (BCN2SmplID.size() != smplId2comb.size()){
			cerr << "FATAL: BCN2SmplID != smplId2comb\n"; exit(234);
		}
		for(unordered_map<string, int>::iterator IT = combiMapCollectGrp.begin(); IT != combiMapCollectGrp.end(); IT++){
			omaps << "\t" << itos(IT->second) + ":" + IT->first;
		}
	}
	omaps << "\n";
	
	//convert to vector, that can than be written out
	vector<shared_ptr<DNAunique>> Dnas (Tracker.size());
	int cnt = 0;
	for (const auto& dd : Tracker) {
		Dnas[cnt] = dd;
		cnt++;
	}
//	bool DNAuPointerCompare(shared_ptr<DNAunique> l, shared_ptr<DNAunique> r) {	return l->getCount() < r->getCount();}
	sort(Dnas.begin(), Dnas.end() , DNAuPointerCompare);
	totSize = 0; size_t passed_hits(0);
	//bool thrHit = false;

	//sanity check
	vector<int> cntspersmpl(BCN2SmplID.size(), 0);
	//print unique DNAs
	for (size_t i = 0; i < Dnas.size(); i++) {
	//for (const auto& dd : Tracker){
		//reformat header
		shared_ptr<DNAunique> dd = Dnas[i];
		dd->Count2Head(b_usearch_fmt);

		if (( bDerepSmpSpcfc && dd->pass_deprep_smplSpc(derepSmpSpcfc)) ||
			pass_deprep_conditions(dd) ) {
			
			passed_hits++;
			totSize += dd->getCount();
			dd->writeSeq(of, b_singleLine);
			
		} else {
			dd->writeSeq(ofRest, b_singleLine);
		}
		dd->writeMap(omaps, dd->getID(), cntspersmpl, smplId2comb);

		//write out full seq + fastq
		dd->resetTruncation();
		dd->prepareWrite(fastqVer);
		dd->writeFastQ(of2);
		if (b_pairedInput) {
			shared_ptr<DNAunique> oD = dd->getPair();
			if (oD != NULL) {
				oD->resetTruncation();
				oD->prepareWrite(fastqVer);
				oD->writeFastQ(of2p2);
			} else {
				shared_ptr<DNA> tmp = make_shared< DNA>("", dd->getID());
				tmp->writeFastQEmpty(of2p2);

			}
		}
	}
//	if (tmpCnt != totSize) {
//		cerr << "Counting failed\n" << tmpCnt << " " << totSize<<endl;
//	}
	of.close(); omaps.close();  ofRest.close();
	float avgSize = (float)totSize / (float)(passed_hits);
	string report = "";
	report += "Dereplication:\nAccepted " + intwithcommas((int)passed_hits) +" unique sequences ( " + minCopiesStr;
	if (bDerepSmpSpcfc) { report += " & sample specific restrictions"; }
	report += " )";
	if (passed_hits > 0) {
		report += "; average size in this set is " + ftos(avgSize) + ".\nUniques with insufficient abundance : " + intwithcommas(int(Tracker.size() - passed_hits)) + " not passing derep conditions\n";
	}
	//cerr << tmpCnt << endl;
	cerr << "\n" << report << endl << endl;
	int uneqCnts(0);
	//check cntspersmpl vector
	vector<string> SampleID = mf->SampleID;
	for (unsigned int i = 0; i<SampleID.size(); i++) {
		size_t j(0); bool detected = false;
		for (; j < BCN2SmplID.size(); j++) {
			if (SampleID[i] == BCN2SmplID[j]) { detected = true;  break; }
		}
		if (detected && mf->colStats[0].BarcodeDetected[i] != cntspersmpl[j] && mf->colStats[0].BarcodeDetected[i] != -1) {
			//cerr << "ERROR: Unequal counts for " << SampleID[i] << ": " << mf->colStats.BarcodeDetected[i] << " vs " << cntspersmpl[j] << endl;
			uneqCnts++;
			
		} else if (!detected) {
			cerr << "Could not detect Sample " << BCN2SmplID[j] << " in ref set (check that mapping file is correctly formatted?).\n"; exit(87);
		}
	}
#ifdef DEBUG
	cerr << "Derep Fin" << endl;
#endif
	if (b_pairedInput) { of2.close(); }
	if (uneqCnts>0) {
		cerr << "Unequal counts in " << uneqCnts << " cases. \n";
		//exit(66);
	}
	return report;
}


bool Dereplicate::pass_deprep_conditions(shared_ptr<DNAunique> d) {
	vector<int> x = d->getDerepMapSort(minCopiesSiz);
	int cumSum(0);
	for (size_t i = 0; i < x.size(); i++) {
		cumSum += x[i];
		//if (minCopies[i] == -1) { continue; }
		size_t yy(minCopiesSiz);
		if (yy > i) { yy = i+1; }//if checking for 1 smpl copy, don't need to check for 2 samlpe allowed copy number..
		for (size_t j = 0; j < yy; j++) {
			if (minCopies[j] == -1) { continue; }

			//if (i > minCopiesSiz) { break; }
			if (cumSum >= minCopies[j]) { return true; }
		}
	}
	return false;
}

void writeLog(string& logf) {
	ofstream of;
	of.open(logf.c_str());
	of.close();
}

string additionalFileName(const string& in){
	if (in.find(",") == std::string::npos){
		return additionalFileName2(in);
	}
	else {
		vector<string> vi = splitByCommas(in);
		string out = additionalFileName2(vi[0]);
		for (uint i = 1; i < vi.size(); i++){
			out += "," + additionalFileName2(vi[i]);
		}
		return out;
	}
}
string additionalFileName2(const string& in){
	size_t point = in.find_last_of('.');
	if (point == (size_t)-1){ return in + ".add"; }
	return in.substr(0, point) + ".add." + in.substr(point + 1);
}


//*******************************************
//*        FILTERS OBJECT
//*******************************************


Filters::Filters(OptContainer& cmdArgs) :
		PrimerL(0), PrimerR(0), PrimerL_RC(0), PrimerR_RC(0), PrimerIdx(0),
		Barcode(0), revBarcode(0), Barcode2(0), revBarcode2(0),
		HeadSmplID(0),
		hetPrimer(2,vector<string>(0)),
		demultiSinglFiles(0), demultiSinglFilesF(0),
		colStats(2),
		FastaF(0), QualF(0), FastqF(0), MIDfqF(0),
		derepMinNum(0),
		lMD(NULL), RepStat(2, NULL), RepStatAddition(2, NULL),
		tAdapter(""), tAdapterLength(0),
		bDoAdapter(false), bDoMultiplexing(true), bDoBarcode(true),
		bDoBarcode2(false),
		bDoHeadSmplID(false), bBarcodeSameSize(false),
		bOneFileSample(false), curBCnumber(-1), BCoffset(0),
		bAdditionalOutput(false), b2ndRDBcPrimCk(false),
		bRevRdCk(false), bChkRdPrs(true),
		min_l(0), alt_min_l(0), min_l_p(-1.f), alt_min_l_p(-1.f),
		maxReadLength(0), norm2fiveNTs(false),
		max_l(10000), min_q(0.f), alt_min_q(0.f),
		BcutPrimer(true), alt_BcutPrimer(true), bPrimerR(false),
		bRequireRevPrim(false), alt_bRequireRevPrim(false),
		bRequireFwdPrim(false), alt_bRequireFwdPrim(false), BcutTag(true),
		bCompletePairs(false), bShortAmplicons(false),
		MinTagLen(0), MinTagLen2(0),MaxTagLen(0), MaxTagLen2(0), MinPrimLen(0), maxHomonucleotide(0),
		PrimerErrs(0), alt_PrimerErrs(0), TagErrs(0),
		MaxAmb(-1), alt_MaxAmb(-1), 
		FQWwidth(0), EWwidth(0),
		RevPrimSeedL(5), 
		b_BinFilBothPairs(false),
		BinFilErr(2.5), BinFilP(-1.f),
		alt_FQWthr(0), alt_EWthr(0),
		PEheaderVerWr(0), TrimStartNTs(0), TruncSeq(-1),
		userReqFastqVer(0), userReqFastqOutVer(33), maxAccumQP(-1),
		alt_maxAccumQP(-1),
		pairedSeq(-1),
		revConstellationN(0),
		BCdFWDREV(2),
		Xreads(-1),
		restartSet(false),b_optiClusterSeq(false),
		b_subselectionReads(false), b_doQualFilter(true),
		b_doFilter(true),
		bDoDereplicate(false), bDoCombiSamples(false),
		bDoDemultiplexIntoFiles(false), 
		maxReadsPerOFile(0),ReadsWritten(0),OFileIncre(0),
		Barcode_len(0), Barcode2_len(0)
//		dPDS(make_unique()), dHDS(make_unique())
		{
			
	colStats[0] = collectstats(); colStats[1] = collectstats();
	if (bAdditionalOutput){
		statAddition = collectstats();
	}
	bool alt_bRequireRevPrimSet=false;


	string optF ("");
	if (cmdArgs.find("-options") != cmdArgs.end()) {
		optF = cmdArgs["-options"];
	}

	iniSpacer = cmdArgs["-sample_sep"];
	//***************************************
	//default options
	int maxAmb(0),PrimerErrs(1),TagErrs(0);
	float minQual(25);
	float minL(250.f);
	int maxL(1000);
	int QualWinWidth = 50;
	float QualWinThr = 0;
	int EndWinWidth = 15;
	float EndWinThr = 20;
	int maxHomoNT(12);
	bool keepTag(false),keepPrimer(false);
	bool addModConf = false;

	//set up some basic objects
	if ( cmdArgs.find("-paired") != cmdArgs.end() ) {
		pairedSeq = atoi(cmdArgs["-paired"].c_str()); //fakeEssentials();
		if ( pairedSeq<1 || pairedSeq>3 ) { cerr << "Argument \"-paired\" supplied with unknown parameter. Aborting.\n"; exit(28); }
		if ( cmdArgs["-onlyPair"] == "1" || cmdArgs["-onlyPair"] == "2" ) {
			pairedSeq = 1;
		}
	}
	if (cmdArgs.find("-normRdsToFiveNTs") != cmdArgs.end()) {
		norm2fiveNTs = true;
		cerr << "Warning: normRdsToFiveNTs is not implemented!\n";
	}
	//delimit output file size to X reads
	if ( cmdArgs.find("-maxReadsPerOutput") != cmdArgs.end() ) {
		maxReadsPerOFile = atoi(cmdArgs["-maxReadsPerOutput"].c_str());
	}
	//important for fastq format
	if ( cmdArgs.find("-i_qual_offset") != cmdArgs.end() ) {
		if ( cmdArgs["-i_qual_offset"] == "auto" ) {
			userReqFastqVer = 0;
		} else {
			userReqFastqVer = atoi(cmdArgs["-i_qual_offset"].c_str());
		}
	}
	//cerr<<cmdArgs["-o_qual_offset"]<<endl;
	userReqFastqOutVer = atoi(cmdArgs["-o_qual_offset"].c_str());
	//statistic tracker
	for ( size_t i = 0; i < 2; i++ ) {
		RepStat[i] = make_shared<ReportStats>(true);
		RepStatAddition[i] = make_shared<ReportStats>(true);
	}
	PreFiltP1 = make_shared<ReportStats>(true);
	PreFiltP2 = make_shared<ReportStats>(true);
	//do new SEED sequence selection?
	if ( cmdArgs.find("-optimalRead2Cluster") != cmdArgs.end() ) {
		b_optiClusterSeq = true;
	}
	//do selection of specific reads?
	if ( cmdArgs["-specificReads"] != "" ) {
		b_subselectionReads = true;
	}
	if (cmdArgs.find("-binomialFilterBothPairs") != cmdArgs.end() && cmdArgs["-binomialFilterBothPairs"] == "1") {
		b_BinFilBothPairs = true;
	}
	//***************************************
	//read options
	ifstream opt;
	opt.open(optF.c_str(),ios::in);
	if ( !opt || optF=="") {
		cerr << "NO filtering will be done on your reads (just rewriting / log files created)." << endl;
		b_doFilter = false;
		return;
	}
	string line;
	while (getline(opt,line,'\n')){
		
		if (line.length()<=1 || line.substr(0,1)=="#"){
			continue;
		}
		
		bool addMod = false;
		if (line.substr(0,1)=="*"){
			addMod= true;
			line = line.substr(1);
		}
		string segs;
		string segs2;
		stringstream ss;
		ss << line;
		getline(ss,segs,'\t');
		getline(ss,segs2,'\t');

		if (cmdArgs["-XfirstReads"] != "") {
			Xreads = atoi(cmdArgs["-XfirstReads"].c_str());
		}


		if (strcmp(segs.c_str(),"minSeqLength") == 0){
			if (addMod){ 
				 float tmp = (float)atof(segs2.c_str());
				 if (tmp>1.f) {
					 alt_min_l = (int)tmp;
				 } else {
					 alt_min_l = -1;
					 alt_min_l_p = tmp;
				 }
				if (alt_min_l != minL){	addModConf = true;	}
			} else {
				minL = (float)atof(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"maxSeqLength") == 0){
			maxL = atoi(segs2.c_str());
		} else if (strcmp(segs.c_str(),"minAvgQuality") == 0){
			if (addMod){ 
				alt_min_q = (float) atof(segs2.c_str());
				if (alt_min_q!= minQual){addModConf = true;}
			} else {
				minQual = (float) atof(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"maxAmbiguousNT") == 0){
			if (addMod){ 
				alt_MaxAmb = atoi(segs2.c_str());
				if (MaxAmb!= alt_MaxAmb){addModConf = true;}
			} else {
				maxAmb = atoi(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"QualWindowThreshhold") == 0){
			if (addMod){ 
				alt_FQWthr = (float) atof(segs2.c_str());
				if (alt_FQWthr!= QualWinThr){addModConf = true;}
			} else {
				QualWinThr = (float) atof(segs2.c_str());
			}
		}
		else if (strcmp(segs.c_str(), "QualWindowWidth") == 0){
			QualWinWidth = atoi(segs2.c_str());
		}
		else if (strcmp(segs.c_str(), "BinErrorModelMaxExpError") == 0){
			BinFilErr = (float) atof(segs2.c_str());
			if (BinFilErr < 0){
				cerr << "BinErrorModelMaxExpError was set to <0. Set to 0 instead.\n";
				BinFilErr = 0;
			}
		}
		else if (strcmp(segs.c_str(), "BinErrorModelAlpha") == 0){
			BinFilP = (float)atof(segs2.c_str());
			if (BinFilP != -1.f && (BinFilP<0.f || BinFilP>1.f)){
				cerr << "BinErrorModelAlpha has to be between 0 and 1 (or -1 to deactivate).\nAborting..\n";
				exit(542);
			}

		} else if (strcmp(segs.c_str(),"TrimWindowWidth") == 0){
			EndWinWidth = atoi(segs2.c_str());
		} else if (strcmp(segs.c_str(),"TrimWindowThreshhold") == 0){
			if (addMod){ 
				alt_EWthr = (float) atof(segs2.c_str());
				if (alt_EWthr != EndWinThr){addModConf = true;}
			} else {
				EndWinThr = (float) atof(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"maxBarcodeErrs") == 0){
			TagErrs = atoi(segs2.c_str());
		} else if (strcmp(segs.c_str(),"maxPrimerErrs") == 0){
			if (addMod){ 
				alt_PrimerErrs = atoi(segs2.c_str());
				if (alt_PrimerErrs!=PrimerErrs){ addModConf = true;}
			} else {
				PrimerErrs = atoi(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"keepBarcodeSeq") == 0){
			atoi(segs2.c_str())==0 ? keepTag=false : keepTag=true;
		} else if (strcmp(segs.c_str(),"keepPrimerSeq") == 0){
			if (addMod){ 
				atoi(segs2.c_str())==0 ? alt_BcutPrimer=false : alt_BcutPrimer=true;
				if(alt_BcutPrimer!=keepPrimer) {addModConf = true;}
			} else {
				atoi(segs2.c_str())==0 ? keepPrimer=false : keepPrimer=true;
			}
		} else if (strcmp(segs.c_str(),"maxHomonucleotide") == 0){
			maxHomoNT = atoi(segs2.c_str());
		} else if (strcmp(segs.c_str(),"maxAccumulatedError") == 0){
			if (addMod){
				alt_maxAccumQP = double(atof(segs2.c_str()));
				if(alt_maxAccumQP!=maxAccumQP) {addModConf = true;}
			} else{
				maxAccumQP = double(atof(segs2.c_str()));
			}
		} else if (strcmp(segs.c_str(),"TechnicalAdapter") == 0){
			tAdapter = segs2.c_str();
			transform(tAdapter.begin(), tAdapter.end(),tAdapter.begin(), ::toupper);
			tAdapterLength = (int) tAdapter.length();
			bDoAdapter = true;
		}else if (segs == "PEheaderPairFmt"){
			PEheaderVerWr = atoi(segs2.c_str());
		} else if (segs == "TrimStartNTs"){
			TrimStartNTs = atoi(segs2.c_str());
		} else if (segs == "fastqVersion"){
			if (segs2 == "auto") {
				userReqFastqVer = 0;
			} else {
				userReqFastqVer = FastqVerMod(atoi(segs2.c_str()));
			}
		} else if (segs == "RejectSeqWithoutRevPrim"){
			if (addMod){ 
				alt_bRequireRevPrimSet=true;
				if (segs2=="T"){alt_bRequireRevPrim=true;
				} else {alt_bRequireRevPrim=false;}				
				if(alt_bRequireRevPrim!=bRequireRevPrim){addModConf = true;}
			} else {
				if (segs2=="T"){bRequireRevPrim=true;
				} else {bRequireRevPrim=false;}
			}
		} else if (segs == "RejectSeqWithoutFwdPrim"){
			if (addMod){ 
				alt_bRequireFwdPrim=true;
				if (segs2=="T"){alt_bRequireFwdPrim=true;
				} else {alt_bRequireFwdPrim=false;}				
				if(alt_bRequireFwdPrim!=bRequireFwdPrim){addModConf = true;}
			} else {
				if (segs2=="T"){bRequireFwdPrim=true;
				} else {bRequireFwdPrim=false;}
			}
		} else if (segs == "TruncateSequenceLength") {
			TruncSeq = atoi(segs2.c_str());
			if (TruncSeq != -1 && TruncSeq < (int)minL) { minL = (float)TruncSeq; }
		} else if (segs == "AmpliconShortPE") {
			if (segs2 == "T") {
				bShortAmplicons = true;
			} else { bShortAmplicons = false; }
		} else if (segs == "CheckForMixedPairs") {
			if (segs2 == "T") {
				b2ndRDBcPrimCk = true;
			} else { b2ndRDBcPrimCk = false; }
		} else if ( segs == "CheckForReversedSeqs" ) {
			if ( segs2 == "T" ) {
				bRevRdCk = true;
			} else { bRevRdCk = false; }
		}
		else if (segs == "SyncReadPairs") {
			if (segs2 == "T") {
				bChkRdPrs = true;
			}
			else { bChkRdPrs = false; }
		}
		

	}
	
	//report some non-std options
	if (bShortAmplicons){
		cerr << "Checking for reverse primers on 1st read.\n";
	}
	if (b2ndRDBcPrimCk){
		cerr << "Checking for switched pairs.\n";
	}

	opt.close();
	//set in filter object
	this->setSeqLength(minL,maxL);
	this->setPrimerErrs(PrimerErrs);
	this->setTagErrs(TagErrs);
	this->removePrimer(!keepPrimer);
	this->removeTag(!keepTag);
	this->setMaxAmb(maxAmb);
	this->setAvgMinQual(minQual);
	this->setFloatingQWin(QualWinWidth,QualWinThr);
	this->setFloatingEWin(EndWinWidth,EndWinThr);
	this->setMaxHomo(maxHomoNT);
	//alternative options (mid qual filtering)
	if (addModConf){
		if (!alt_bRequireRevPrimSet){alt_bRequireRevPrim = bRequireRevPrim;}
		if (cmdArgs.find("-o_fna")  != cmdArgs.end() && cmdArgs["-o_fna"].length()>1){
			if (cmdArgs.find("-o_fna2")  == cmdArgs.end()){
				cmdArgs["-o_fna2"] = additionalFileName(cmdArgs["-o_fna"]);
				//cmdArgs["-o_fna2"] = cmdArgs["-o_fna"].substr(0,cmdArgs["-o_fna"].length()-4)+".add.fna";
			}
		} else if (cmdArgs.find("-o_fastq")  != cmdArgs.end() && cmdArgs["-o_fastq"].length()>1){
			if (cmdArgs.find("-o_fastq2")  == cmdArgs.end()){
				cmdArgs["-o_fastq2"] = additionalFileName(cmdArgs["-o_fastq"]);
			}
		}
		bAdditionalOutput=true;
	}
	

	

}

//copy of filter, but leave stat vals empty
Filters::Filters(shared_ptr<Filters> of, int BCnumber, bool takeAll) :
		PrimerL(0,""), PrimerR(0,""),
		PrimerL_RC(0,""), PrimerR_RC(0,""),
		PrimerIdx(of->PrimerIdx), 
		Barcode(0), revBarcode(0), Barcode2(0), revBarcode2(0),
		HeadSmplID(0), hetPrimer(2, vector<string>(0)), colStats(2),
		FastaF(of->FastaF), QualF(of->QualF), FastqF(of->FastqF),
		MIDfqF(of->MIDfqF), derepMinNum(of->derepMinNum),
		lMD(NULL), RepStat(2,NULL), RepStatAddition(2,NULL),
		tAdapter(of->tAdapter),tAdapterLength(of->tAdapterLength),
		bDoAdapter(of->bDoAdapter),bDoMultiplexing(of->bDoMultiplexing),
		bDoBarcode(of->bDoBarcode), bDoBarcode2(of->bDoBarcode2),
		bDoHeadSmplID(of->bDoHeadSmplID),
		bBarcodeSameSize(of->bBarcodeSameSize),
		bOneFileSample(of->bOneFileSample), curBCnumber(BCnumber), BCoffset(0),
		bAdditionalOutput(of->bAdditionalOutput), b2ndRDBcPrimCk(of->b2ndRDBcPrimCk),
		bRevRdCk(of->bRevRdCk), bChkRdPrs(of->bChkRdPrs),
		min_l(of->min_l), alt_min_l(of->alt_min_l), min_l_p(of->min_l_p), alt_min_l_p(of->alt_min_l_p),
		maxReadLength(0), norm2fiveNTs(of->norm2fiveNTs),
		max_l(of->max_l), min_q(of->min_q), alt_min_q(of->alt_min_q),
		BcutPrimer(of->BcutPrimer),alt_BcutPrimer(of->alt_BcutPrimer),
		bPrimerR(of->bPrimerR),
		bRequireRevPrim(of->bRequireRevPrim),alt_bRequireRevPrim(of->alt_bRequireRevPrim),
		bRequireFwdPrim(of->bRequireFwdPrim),alt_bRequireFwdPrim(of->alt_bRequireFwdPrim),
		BcutTag(of->BcutTag),
		bCompletePairs(of->bCompletePairs), bShortAmplicons(of->bShortAmplicons),
		MinTagLen(of->MinTagLen), MinTagLen2(of->MinTagLen2), MaxTagLen(of->MaxTagLen), MaxTagLen2(of->MaxTagLen2), MinPrimLen(of->MinPrimLen), maxHomonucleotide(of->maxHomonucleotide),
		PrimerErrs(of->PrimerErrs),alt_PrimerErrs(of->alt_PrimerErrs),TagErrs(of->TagErrs),
		MaxAmb(of->MaxAmb),	alt_MaxAmb(of->alt_MaxAmb),
		FQWwidth(of->FQWwidth),EWwidth(of->EWwidth),
		RevPrimSeedL(of->RevPrimSeedL),
		b_BinFilBothPairs(of->b_BinFilBothPairs),
		BinFilErr(of->BinFilErr), BinFilP(of->BinFilP),
		FQWthr(of->FQWthr),EWthr(of->EWthr),
		alt_FQWthr(of->alt_FQWthr),alt_EWthr(of->alt_EWthr),
		PEheaderVerWr(of->PEheaderVerWr),TrimStartNTs(of->TrimStartNTs),
		TruncSeq(of->TruncSeq),
		iniSpacer(of->iniSpacer),userReqFastqVer(of->userReqFastqVer),
		userReqFastqOutVer(of->userReqFastqOutVer),maxAccumQP(of->maxAccumQP),
		alt_maxAccumQP(of->alt_maxAccumQP),
		//BChit, BCrevhit initialize to 0 - new set, new luck
		pairedSeq(of->pairedSeq),
		revConstellationN(0),
		BCdFWDREV(of->BCdFWDREV),
		restartSet(false),
		b_optiClusterSeq(of->b_optiClusterSeq), b_subselectionReads(of->b_subselectionReads),
		b_doQualFilter(of->b_doQualFilter), 
		b_doFilter(of->b_doFilter),
		bDoDereplicate(of->bDoDereplicate),
		bDoCombiSamples(of->bDoCombiSamples),
		bDoDemultiplexIntoFiles(of->bDoDemultiplexIntoFiles),
		maxReadsPerOFile(of->maxReadsPerOFile),
		ReadsWritten(of->ReadsWritten), OFileIncre(of->OFileIncre),
		Barcode_len(0), Barcode2_len(0)
//		dPDS(NULL), dHDS(NULL)
		{
	RepStat[0] = make_shared<ReportStats>(of->RepStat[0]->bMedianCalcs);
	RepStat[1] = make_shared<ReportStats>(of->RepStat[1]->bMedianCalcs);
	BCdFWDREV[0].fix(); BCdFWDREV[1].fix();
	RepStatAddition[0] = make_shared<ReportStats>(of->RepStatAddition[0]->bMedianCalcs);
	RepStatAddition[1] = make_shared<ReportStats>(of->RepStatAddition[1]->bMedianCalcs);
	PreFiltP1 = make_shared<ReportStats>(of->PreFiltP2->bMedianCalcs);
	PreFiltP2 = make_shared<ReportStats>(of->PreFiltP2->bMedianCalcs);
	colStats[0] = collectstats(); colStats[1] = collectstats();
	if (bAdditionalOutput) {
		statAddition = collectstats();
	}
	if (takeAll) {
		this->allResize((uint)of->PrimerIdx.size());
		PrimerIdxRev = of->PrimerIdxRev;
		PrimerIdx = of->PrimerIdx;
		Barcode = of->Barcode;
		Barcode2 = of->Barcode2;
		SampleID = of->SampleID;
		SampleID_Combi = of->SampleID_Combi;
		HeadSmplID = of->HeadSmplID;
		PrimerL = of->PrimerL;
		PrimerR = of->PrimerR;
		PrimerL_RC = of->PrimerL_RC;
		PrimerR_RC = of->PrimerR_RC;
		hetPrimer = of->hetPrimer;
		lMD = of->lMD;
		Barcode_len = of->Barcode_len;
		Barcode2_len = of->Barcode2_len;
		demultiSinglFiles = of->demultiSinglFiles;
		demultiSinglFilesF = of->demultiSinglFilesF;
		BarcodePreStats();

	}

}

Filters::~Filters() {
#ifdef DEBUG
	cerr << "Deleting filter .. " ;
#endif

//	for (size_t i = 0; i < 2; i++) { delete RepStat[i]; delete RepStatAddition[i]; }
//	delete PreFiltP1; delete PreFiltP2;
//	delete dPDS; delete dHDS;
#ifdef DEBUG
	cerr << "Done" << endl;
#endif
}

//simulates that in mapping file links to sequence file was given.
bool Filters::setcmdArgsFiles(OptContainer& cmdArgs){

	if (FastqF.size()==0 && QualF.size()==0 && FastaF.size() > 0){
		//fasta entry but no qual entries

		string path="";
		if (cmdArgs.find("-i_path")  != cmdArgs.end() && cmdArgs["-i_path"].length() > 2){
			path=cmdArgs["-i_path"] + string("/");
		}

		QualF.resize(FastaF.size());
		for (unsigned int i=0; i< FastaF.size(); i++){
			string newQ = FastaF[i];
			int pos = (int) newQ.find_last_of(".");
			newQ = newQ.substr(0,pos);
			newQ += string(".qual");
			fstream fin;
			string fullQ = path + newQ;
			fin.open(fullQ.c_str(),ios::in);
			if( fin.is_open() )	{
				cerr<<"Using quality file: "<<fullQ <<endl;
			} else if (cmdArgs.find("-number")!=  cmdArgs.end() && cmdArgs["-number"] =="T"){
				;
			}else {
				cerr<<"You did not supply a quality file for"<<path+ FastaF[i]<<". \nPlease give the path to your quality file as command line argument:\n  -i_qual <PathToQualityFile>\n";
				newQ = "";
				//fin.close();return false;
			}
			fin.close();
			QualF[i] = newQ;
		}
	}

	int fileSiz = (int)Barcode.size();
	//instead of max:
	if (!bDoMultiplexing){fileSiz=1;}

	if (FastaF.size()==0 && FastqF.size()==0){
		//set up fasta/fastq vector specific to corresponding BC (that should be in this file)
		if (cmdArgs.find("-i_fastq")  == cmdArgs.end()){
			FastaF.resize(fileSiz);
			QualF.resize(fileSiz);
			for (unsigned int i=0; i< FastaF.size(); i++){
				FastaF[i] = cmdArgs["-i_fna"];
				QualF[i] = cmdArgs["-i_qual"];
			}
		} else {//fastq input
			vector<string> fqTmp (1,cmdArgs["-i_fastq"]);
			if (cmdArgs["-i_fastq"].find(";") != string::npos) {//";" denotes several files
				if (fileSiz == 1) {//no BC, 
					fqTmp = splitByCommas(cmdArgs["-i_fastq"], ';');
					this->allResize((uint) fqTmp.size());
					fileSiz = (int) fqTmp.size();
					cerr << "Detected " << fileSiz << " input files (pairs)." << endl;
					FastqF = fqTmp;
				} else {
					cerr << "Fastq string contains symbol \";\". Not allowed in input string"; exit(32);
				}
			} else {
				FastqF.resize(fileSiz, cmdArgs["-i_fastq"]);
			}
		}
	}


	if (MIDfqF.size() == 0)
		if (cmdArgs.find("-i_MID_fastq") != cmdArgs.end()) {
		MIDfqF.resize(fileSiz, "");
		for (unsigned int i = 0; i < MIDfqF.size(); i++) {
			MIDfqF[i] = cmdArgs["-i_MID_fastq"];
		}
	}
	if (cmdArgs.find("-o_demultiplex") != cmdArgs.end()) {	
		generateDemultiOutFiles(cmdArgs["-o_demultiplex"]);
	}
		

	if (cmdArgs["-o_dereplicate"] != "") {
		//check if file could exist
		ofstream temp;
		temp.open(cmdArgs["-o_dereplicate"].c_str(), ios::out);
		if (!temp) { cerr << "Could not open outstream to dereplicated sequences:\n" << cmdArgs["- o_dereplicate"] << endl; exit(78); }
		temp.close();
		bDoDereplicate = true;
	}

	return true;
}


void Filters::generateDemultiOutFiles( string path) {
	//fill in demultiSinglFiles vector
	vector<ofbufstream*> empVec(2, NULL);
	vector<string> empVec2(2, "");
	
	struct stat info;
	if (stat(path.c_str(), &info) != 0) {
		cerr << "Output path for demultiplexed files does not exist, please create this directory:\n" << path << endl;
		exit(833);
	}else if (info.st_mode & S_IFDIR) { // S_ISDIR() doesn't exist on my windows 
		cerr<<"Writing demultiplexed files to: "<<path<<endl;// printf("%s is a directory\n", path);
	}else {
		cerr << path << " is no directory\n"; exit(834);
	}
	path += "/";

	bDoDemultiplexIntoFiles = true;
	bool openOstreams = true; uint ostrCnt(0);
	demultiSinglFiles.resize(SampleID.size(), empVec);
	demultiSinglFilesF.resize(SampleID.size(), empVec2);
	for (size_t i = 0; i < SampleID.size(); i++) {
		//actually needs to know if paired files..
		//if (ostrCnt > maxFileStreams) {openOstreams = false;}
		if (pairedSeq == 1|| pairedSeq == -1) {
			string nfile = path + SampleID[i] + ".fq";
			if (openOstreams) { demultiSinglFiles[i][0] = new ofbufstream(nfile.c_str(), ios::out); }
			demultiSinglFilesF[i][0] = nfile;
			ostrCnt++;
		}
		else {
			string nfile = path + SampleID[i] + ".1.fq";
			if (openOstreams) { demultiSinglFiles[i][0] = new ofbufstream(nfile.c_str(), ios::out); }
			demultiSinglFilesF[i][0] = nfile;
			nfile = path + SampleID[i] + ".2.fq";
			if (openOstreams) { demultiSinglFiles[i][1] = new ofbufstream(nfile.c_str(), ios::out); }
			demultiSinglFilesF[i][1] = nfile;
			ostrCnt += 2;
		}
	}

}
//only does BC 1
void Filters::reverseTS_all_BC(){
//	for (int i=0; i<Barcode.size();i++){
//		reverseTS(Barcode[i]);
//	}
	Barcode = revBarcode;
	revBarcode.resize(0);
	BCList.clear();
	for (uint i = 0; i < Barcode.size(); i++){
		BCList[Barcode[i]] = i;
		Barcode_len[i] = (int) Barcode[i].length();
	}

}
void Filters::reverseTS_all_BC2() {
	//	for (int i=0; i<Barcode.size();i++){
	//		reverseTS(Barcode[i]);
	//	}
	Barcode2 = revBarcode2;
	revBarcode2.resize(0);
	BCList2.clear();
	for ( uint i = 0; i < Barcode2.size(); i++ ) {
		BCList2[Barcode2[i]] = i;
		Barcode2_len[i] = (int)Barcode2[i].length();
	}

}


void Filters::preFilterSeqStat(shared_ptr<DNA> d, int pair) {
	if (pair <= 0) {
		PreFiltP1->addDNAStats(d);
	} else if (pair == 1) {
		PreFiltP2->addDNAStats(d);
	}
	updateMaxSeqL(d->length());
}
void Filters::updateMaxSeqL(int x) {
	if (x < maxReadLength) { return; }
	maxReadLength = x;
	if (min_l_p != -1.f) {
		min_l = (int)((float)maxReadLength * min_l_p);
	}
}
void Filters::setSeqLength(float minL, int maxL) {
	if (minL>1.f) {
		min_l = (int)minL;
		min_l_p = -1.f;
	} else {
		min_l = -1;
		min_l_p = minL;
	}
	max_l = maxL;
}


//ever_best is the best %id that was ever observed for this cluster match
bool Filters::betterSeed(shared_ptr<DNA> d1, shared_ptr<DNA> d2, shared_ptr<DNA> ref, shared_ptr<DNA> ref2, float ever_best, 
	uint bestL, int usePair, bool checkBC) {
	float d1pid(d1->getTempFloat()), refpid(ref->getTempFloat());
	int TagIdx(0);
	if (checkBC) {
		TagIdx = -2;
	}
	//0.2% difference is still ok, but within 0.5% of the best found seed (prevent detoriating sequence match)
	//float blen = (float)ref->length() + (float)d1->length();
	uint curL = d1->length();
	if (d2 != NULL) {		curL += d2->length();	}
	if (float(curL) / float(bestL) < BestLengthRatio) { return false; }
	if (d1pid<refpid - 0.4f || d1pid < ever_best - 1){ return false; }

	//*** DNA1
	//needs to quality filter first
	if (!check(d1, true, usePair, TagIdx)) {
		return false;
	}
	//at least 90% length of "good" hit
	if (d1->length() / ref->length() < RefLengthRatio) { return false; }

	//checks if the new DNA has a better overall quality
	//1 added to qual, in case no Qual DNA is used
	float thScore = (1+d1->getAvgQual())*(d1pid ) * log((float)d1->length() );
	float rScore = (1+ref->getAvgQual())*(refpid ) * log((float)ref->length() );
	if (thScore > rScore){
		//also check for stable lowest score
		if (d1->minQual() > ref->minQual() - MinQualDiff && (d2 == NULL || ref2 == NULL)) { return true; }
	}
	if (d2 == NULL || ref2 == NULL) {
		return false;
	}
	//*** DNA2
	//second pair likely to be of worse qual, but only direct comparison relevant here

	check(d2, true, 1, TagIdx);
	//at least 90% length of "good" hit
	if (d2->length() / ref2->length() < RefLengthRatio) { return false; }

	//checks if the new DNA has a better overall quality
	//weigh with average id to OTU seed
	thScore +=(1+ d2->getAvgQual()) * log((float)d2->length()) * 97;
	rScore += (1+ref2->getAvgQual()) * log((float)ref2->length()) * 97;
	if (thScore > rScore) {
		return true;
	}

	return false;
}
bool Filters::check(shared_ptr<DNA> d, bool doSeeding, int pairPre,
	int &tagIdx) {
	//bool qualWinTrim(false);// , TecAdap(false); //RevPrimFound(false), 	bool AccErrTrim(false);
	int pair = max(0, pairPre);//corrects for -1 (undefined pair) to set to 0
	unsigned int hindrance = 0;

	//sTotalPlus(pair);
	if (check_length(d->length())) {
		d->QualCtrl.minL = true; //sMinLength(pair);	
		return false;
	}
	//remove technical adapter
	if (pair == -1) {
		if (bDoAdapter){
			remove_adapter(d);
		}
	}

	//BC already detected (e.g. MID)?
	if (tagIdx == -2){
		tagIdx = cutTag(d, pair == 0); //barcode 2nd part
	}
	if ((bDoBarcode2||bDoBarcode)&& tagIdx < 0) {
		 d->QualCtrl.TagFail = true; d->failed();	return false;//sTagFail(pair);
	}
	

	//cerr<<tagIdx<<" "<<PrimerIdx[tagIdx]<<endl;	cerr<<"next\n";
	//fwd primer only on pair 0
	if (BcutPrimer){
		if (pair != 1) {//0 or -1
			if (!cutPrimer(d, PrimerIdx[tagIdx], false,pair) && bRequireFwdPrim) {
				d->failed(); return false;
			}
		}
		else if (bShortAmplicons) {//pair == 1, check for fwd primer in pair 2 (rev-compl)
			cutPrimer(d, PrimerIdx[tagIdx], true,pair);
		}
	}


	if (check_length(d->length())){
		d->QualCtrl.minL = true; //sMinLength(pair);	
		return false;
	}
	if (max_l!=0 && d->length()-hindrance > max_l){
		d->QualCtrl.maxL = true; //sMaxLength(pair);	
		return false;
	}


	
	//rev primer is the first that needs to be looked for
	//makes it slower, as higher chance for low qual and this routine is costly.. however more important to get good lock on rev primer
	if ((pair != 0 || bShortAmplicons) && bPrimerR) {
		//removal of reverse primer
		bool revCheck = pair == -1 || pair == 0;//1:false for RC, else always a reverse check
		cutPrimerRev(d, PrimerIdxRev[tagIdx], revCheck);
		if (d->getRevPrimCut()) {
			//RevPrimFound = true;
			d->QualCtrl.PrimerRevFail = false;
			//check length
			if (check_lengthXtra(d, hindrance)) {
				d->QualCtrl.minL = true; //sMinLength(pair);	
				d->failed();	return false;
			}
		} else  {//stats, but only for 2nd pair
			d->QualCtrl.PrimerRevFail = true; //sRevPrimerFail(pair);
			if (pair == 1 && bRequireRevPrim) {//failed to find reverse primer
				return false;
			}
		}
	}


	if (doSeeding){
		//cut off low qual, hard limits
		//int TWwidth = 10; float TWthr = 20;
		d->qualWinPos(EWwidth, EWthr);//) { d->QualCtrl.Trimmed = true; }// sTrimmed(pair); }
		//SEED is needed for building trees, accumulated error not interesting
		//double tmpAccumQP = 2.;
		//if (!d->qualAccumTrim(tmpAccumQP)){ colStats.AccErrTrimmed++; }
//		if (Trim){ colStats.Trimmed++; }
//		if (AccErrTrim){ colStats.AccErrTrimmed++; }

		return true;
	}
	
	
	//if seq needs to be cut, than here
	if (TruncSeq>0){d->cutSeq(TruncSeq,-1,true);
		if ( check_length(d->length()) ){
			d->QualCtrl.minL = true; //sMinLength(pair);	
			return false;
		}
	}
	if (b_doQualFilter) {
		//second cut off low qual
		d->qualWinPos(EWwidth, EWthr);
		//cut off accumulation error larger than maxAccumQP
		d->qualAccumTrim(maxAccumQP);
		//if (check_length(d->length(),hindrance) ){sMinLength();	return false;}
		if (check_length(d->length())) {
			d->QualCtrl.minLqualTrim = true; return false;//sMinQTrim(pair);	
		}
		int rea = 2;
		if ((min_q > 0 || FQWthr > 0) && d->qualWinfloat(FQWwidth, FQWthr, rea) < min_q) {
			d->QualCtrl.AvgQual = true;//sAvgQual(pair);
			return false;
		}
		if (rea == 1) {
			d->QualCtrl.QualWin = true; //sQualWin(pair);
			return false;
		}
		//binomial filter here
		if (b_BinFilBothPairs || pair != 1 ){
			float ExpErr = d->binomialFilter((int)BinFilErr, BinFilP);
			if (ExpErr > BinFilErr){
				d->QualCtrl.BinomialErr = true;
				//sBinomError(pair, ExpErr);
				d->failed(); return false;
			}
		}
	}

	if (MaxAmb!=-1 && d->numACGT() > MaxAmb){
		d->QualCtrl.MaxAmb = true;//sMaxAmbig(pair);
		return false;
	}
	if (maxHomonucleotide!=0 && !d->HomoNTRuns(maxHomonucleotide)){
		d->QualCtrl.HomoNT = true;// sHomoNT(pair);
		return false;
	}

	//cerr<<d->getID()<<endl;
	//if (RevPrimFound) {//this sequence has been filtered for the reverse primer
	//	sReversePrimerFnd(pair);
	//}

	//adapter removed, quality filtering done. If no map is provided, that is all that is needed
	if (!bDoMultiplexing){
		if (TrimStartNTs>0){
			if (d->length()-TrimStartNTs > max_l){//length check
				d->QualCtrl.maxL = true; //sMaxLength(pair);	
				return false;
			}
			//remove start NTs
			d->cutSeq(0,TrimStartNTs);

		}
		//MUTEX
		//colStats[pair].totalRejected--;
		//if (qualWinTrim || AccErrTrim) {d->QualCtrl.Trimmed = true;}//sTrimmed(pair); }
		d->setPassed(true);
		return true;
	}
	d->setPassed(true);

	/*if (qualWinTrim || AccErrTrim) {
		d->QualCtrl.Trimmed = true;//sTrimmed(pair);
		if (AccErrTrim) {d->QualCtrl.AccErrTrimmed = true;//colStats[pair].AccErrTrimmed++; }
		if (qualWinTrim) { d->QualCtrl.QWinTrimmed = true; }//colStats[pair].QWinTrimmed++;
		}
	}*/
	//if (TecAdap) { colStats[pair].adapterRem++; }
	//colStats[pair].totalRejected--;

	//keep control over passed / not as close as possible to source
	return true;
}
//DNA qual check, and some extra parameters
bool Filters::checkXtra(shared_ptr<DNA> d, int pairPre, int &tagIdx) {

	//bool qualWinTrim(false), AccErrTrim(false); // RevPrimFound(false),TecAdap(false),  
	unsigned int hindrance = 0;
	int pair = max(0, pairPre);//corrects for -1 (undefined pair) to set to 0

//	sTotalPlus(pair);
//	sTotalPlusXtra(pair);
	//remove technical adapter
	if (pairPre == -1) {
		if (bDoAdapter){
			remove_adapter(d);
		}
	}
	//set in outer routines that check mid (-1) or needs to be checked here(-2)
	if (tagIdx == -2){
		if ( bDoBarcode2 && pair == 1 ) {
			tagIdx = cutTag(d, false); //barcode 2nd part
		} else if(bDoBarcode && pair == 0) {
			tagIdx = cutTag(d,true); //barcode
		}
		else {
			tagIdx = 0;
		}
	}
	if ((bDoBarcode || bDoBarcode2) && tagIdx < 0) { d->QualCtrl.TagFail = true;return false; }// sTagFail(pair);	

	
	

	//cerr<<tagIdx<<" "<<PrimerIdx[tagIdx]<<endl;	cerr<<"next\n";
	if (BcutPrimer){
		if (pair != 1) {//0 or -1
			if ( !cutPrimer(d, PrimerIdx[tagIdx], false,pair) && bRequireFwdPrim) {
				d->failed(); return false;
			}
		}
		else if (bShortAmplicons) {//pair == 1, check for fwd primer in pair 2 (rev-compl)
			cutPrimer(d, PrimerIdx[tagIdx], true,pair);
		}
	}
	if (check_lengthXtra(d)){
		d->QualCtrl.minL = true;//sMinLength(pair);	
		d->failed(); return false;
	}
	if (max_l!=0 && d->length()-hindrance > max_l){
		d->QualCtrl.maxL = true; //sMaxLength(pair);	
		d->failed(); return false;
	}

	
	//rev primer is the first that needs to be looked for
	//makes it slower, as higher chance for low qual and this routine is costly.. however more important to get good lock on rev primer
	if ((pair != 0 || bShortAmplicons) && bPrimerR) {
		//removal of reverse primer
		//bool revPrm(false);
		bool revCheck = pair == -1 || pair == 0;//1:false for RC, else always a reverse check
		cutPrimerRev(d, PrimerIdxRev[tagIdx], revCheck);
		if (d->getRevPrimCut()) {
			//RevPrimFound = true;
			d->QualCtrl.PrimerRevFail = false;
			//check length
			if (check_lengthXtra(d, hindrance)) {
				d->QualCtrl.minL = true; //sMinLength(pair);	
				d->failed(); return false;
			}
		} else  {//stats, but only for 2nd pair
			//sRevPrimerFail(pair);
			d->QualCtrl.PrimerRevFail = true; // maybe replace later with routine that checks for FtsDetected
			if (pair != 0 && bRequireRevPrim) {//failed to find reverse primer
				if (alt_bRequireRevPrim) {
					d->failed(); return false;
				} else {
					//RevPrimFound = false;
					d->QualCtrl.PrimerRevFail = false;
					d->setMidQual(true);
				}
			}
		}
	}

	//if seq needs to be cut, than here
	if (TruncSeq>0){d->cutSeq(TruncSeq,-1,true);
		if ( check_lengthXtra(d) ){
			//d->QualCtrl.minL = true; //sMinLength(pair);
			d->failed(); return false;
		}
	}

	if (b_doQualFilter) {
		//second cut off low qual
		d->qualWinPos(EWwidth, EWthr);	// { qualWinTrim = true; }
		//cut off accumulation error larger than maxAccumQP
		if (maxAccumQP != -1.0) {
			int cP = d->qualAccumulate(maxAccumQP);
			if (check_lengthXtra(d, 0, cP)) {
				d->isMidQual();  d->QualCtrl.minLqualTrim = true;//sMinQTrim(pair);
				cP = d->qualAccumulate(alt_maxAccumQP);
				if (check_lengthXtra(d, 0, cP)) {//check if passes alt
					d->failed(); return false;
				}
			} else {
				d->qualAccumTrim(maxAccumQP);//) {  AccErrTrim = true; }
			}
		}


		int rea = 2, rea2 = 2;
		if ((min_q > 0 || FQWthr > 0) && d->qualWinfloat(FQWwidth, FQWthr, rea) < min_q) {
			d->QualCtrl.AvgQual = true; //sAvgQual(pair);
			if ((alt_min_q > 0 || alt_FQWthr > 0) && d->qualWinfloat(FQWwidth, alt_FQWthr, rea2) < alt_min_q) {
				d->QualCtrl.AvgQual= true; //statAddition.AvgQual++;
				d->failed(); return false;
			} else {
				d->QualCtrl.AvgQual = false;
				d->setMidQual(true);
			}
		}
		if (rea == 1) {
			d->QualCtrl.QualWin = true; //sQualWin(pair);
			if (rea2 == 1) {
				//statAddition.QualWin++;
				//d->QualCtrl.QualWin = true;
				d->failed(); return false;
			}
		}
		if (b_BinFilBothPairs || pair == 0){
			float ExpErr = d->binomialFilter((int)BinFilErr, BinFilP);
			if (ExpErr > BinFilErr){
				d->QualCtrl.BinomialErr = true;
				//sBinomError(pair, ExpErr);
				d->failed(); return false;
			}
		}
	}
	int ambNTs = d->numACGT();
	if (MaxAmb!=-1 && ambNTs > MaxAmb){
//		sMaxAmbig(pair);
		d->QualCtrl.MaxAmb = true;
		if (alt_MaxAmb!=-1 && ambNTs>= alt_MaxAmb){
			d->QualCtrl.MaxAmb = true; //statAddition.MaxAmb++;
			d->failed(); return false;
		}
	}
	if (maxHomonucleotide!=0 && !d->HomoNTRuns(maxHomonucleotide)){
		d->QualCtrl.HomoNT = true;//sHomoNT(pair);
		d->failed(); return false;
	}

	//cerr<<d->getID()<<endl;
	//if (RevPrimFound) {//this sequence has been filtered for the reverse barcode
	//	sReversePrimerFnd(pair);
	//}


	//adapter removed, quality filtering done. If no map is provided, that is all that is needed
	if (!bDoMultiplexing){
		if (TrimStartNTs>0){
			if (d->length()-TrimStartNTs > max_l){//length check
				d->QualCtrl.maxL = true; //sMaxLength(pair);	
				d->failed();
				return false;
			}
			//remove start NTs
			d->cutSeq(0,TrimStartNTs);

		}
		//colStats[ pair].totalRejected--;
		//if (qualWinTrim || AccErrTrim) { colStats[pair].Trimmed++; }
		d->setPassed(true);
		return true;
	}

	if (!d->isMidQual()) {
		d->setPassed(true);
	}

	/*	if (pairPre <= 0) {
			statAddition.totalRejected--;
			if (qualWinTrim || AccErrTrim) {
				statAddition.Trimmed++;
				if (AccErrTrim) { statAddition.AccErrTrimmed++; }
				if (qualWinTrim) { statAddition.QWinTrimmed++; }
			} 
			//if (TecAdap) { statAddition.adapterRem++; }
		}
	} else {
		if (qualWinTrim || AccErrTrim) {
			colStats[pair].Trimmed++;
			if (AccErrTrim) { colStats[pair].AccErrTrimmed++; }
			if (qualWinTrim) { colStats[pair].QWinTrimmed++; }
		}
		//if (TecAdap) { colStats[pair].adapterRem++; }
		//keep control over passed / not as close as possible to source
		colStats[pair].totalRejected--;
		
	}*/

	return true;
}
void Filters::noMapMode(OptContainer& cmdArgs){
	string noMapTxt = "sdm run in No Map Mode.";
	if (cmdArgs.find("-paired")  != cmdArgs.end() && (cmdArgs["-paired"]=="2" || cmdArgs["-paired"]=="2")){
		pairedSeq = 2; //fakeEssentials();
		noMapTxt += " Using paired end sequencing files.";
	}

	BcutPrimer = false; bDoBarcode = false; bDoBarcode2 = false;
	bDoAdapter=false;bDoMultiplexing=false;
	bDoHeadSmplID=false;
	fakeEssentials();
	MinTagLen = 0; MinTagLen2 = 0; MaxTagLen = 0; MaxTagLen2 = 0; MinPrimLen = 0;
	cerr<<noMapTxt<<endl;
}
void Filters::fakeEssentials(){
	//create fake entries
	PrimerIdx.push_back(0);Barcode.push_back("NA");
	Barcode_len.push_back(0);
	Barcode2_len.push_back(0);
	PrimerL.push_back(""); PrimerL_RC.push_back(""); SampleID.push_back("NA"); SampleID_Combi.push_back("NA");
	HeadSmplID.push_back("");bDoHeadSmplID=false;
	colStats[0].BarcodeDetected.push_back(-1);
	colStats[1].BarcodeDetected.push_back(-1);
	colStats[0].BarcodeDetectedFail.push_back(-1);
	colStats[1].BarcodeDetectedFail.push_back(-1);
	
}
void Filters::allResize(unsigned int x){
	//cerr<<"resize "<<x<<endl;
	PrimerIdx.resize(x,0);
	PrimerIdxRev.resize(x,0);
	Barcode.resize(x, "");
	Barcode2.resize(x, "");
	Barcode_len.resize(x, 0);
	Barcode2_len.resize(x, 0);
	SampleID.resize(x, "");
	SampleID_Combi.resize(x, "");
	
	colStats[0].BarcodeDetected.resize(x, 0);
	colStats[1].BarcodeDetected.resize(x, 0);
	colStats[0].BarcodeDetectedFail.resize(x, 0);
	colStats[1].BarcodeDetectedFail.resize(x, 0);

	statAddition.BarcodeDetected.resize(x, 0);
	statAddition.BarcodeDetectedFail.resize(x, 0);
	HeadSmplID.resize(x, "");
	vector<ofbufstream*> emptVec(2, NULL);
	vector<string> emptVec2(2, "");
	demultiSinglFiles.resize(x, emptVec);
	demultiSinglFilesF.resize(x, emptVec2);
}

bool Filters::remove_adapter(shared_ptr<DNA> d){ //technical adapter
	//allows for 0 errors, no shifts
	const string& se = d->getSeq();
	for (unsigned int i=0;i<tAdapterLength;i++){
		if (se[i]!=tAdapter[i] ){
			return false;
		}
	}
	d->cutSeq(0,tAdapterLength);
	d->setTA_cut(true);
	return true;
}
//only identifies based on dual BCding
void Filters::dblBCeval(int& tagIdx, int& tagIdx2, string presentBC, shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2) {
	//bool BCfail = false;// , BCfail2 = false;
	if ( tagIdx < 0 || tagIdx2 < 0 || !tdn->getBarcodeDetected() || !tdn2->getBarcodeDetected()) {
		tagIdx = -1; tagIdx2 = -1;
		if (tdn != NULL) { 
			tdn->setPassed(false); /*BCfail = true; */
			tdn->setBCnumber(tagIdx, BCoffset); tdn->setMidQual(false);
		} 
		if (tdn2 != NULL) { tdn2->setPassed(false); tdn2->setMidQual(false); tdn2->setBCnumber(tagIdx2, BCoffset);}
		
		colStats[0].dblTagFail++;
		return;
	}
	string BC1 = Barcode[tagIdx];
	string BC2 = Barcode2[tagIdx2];
	bool hit(false);
	//this routine finds two matching barcodes (as several combinations are possible)
	for ( uint i = 0; i < Barcode.size(); i++ ) {
		if ( Barcode[i] == BC1 && Barcode2[i] == BC2 ) {
			tagIdx = i; tagIdx2 = i; hit = true; break;
		}
	}

	if ( !hit ) {
		//no BC, useless
		tagIdx = -1; tagIdx2 = -1;
		if (tdn != NULL) { tdn->setPassed(false); tdn->setMidQual(false); tdn->setBCnumber(tagIdx, BCoffset);	}
		if (tdn2 != NULL) { tdn2->setPassed(false); tdn2->setMidQual(false); tdn2->setBCnumber(tagIdx2, BCoffset);	}
		return;
	}
	presentBC = BC1 + "|" + BC2;
	//add new BC info to DNA
	//also reset BC in DNA
	if ( tdn != NULL ) {
		BCintoHead(tagIdx, tdn, presentBC, -1, false, true);
		//already done in BCintoHead
	}
	if ( tdn2 != NULL ) {
		BCintoHead(tagIdx2, tdn2, presentBC, -1, true, true);
	}
}

//cuts & identifies - version is just for mid sequences
int Filters::cutTag(shared_ptr<DNA> d, string&presentBC, int& c_err, bool isPair1) {
	if (bDoHeadSmplID) {
		for (unsigned int i = 0; i<HeadSmplID.size(); i++) {
			size_t pos = d->getOldID().find(HeadSmplID[i]);
			if (pos != string::npos) {
				SampleIntoHead(i, d, pos);
				return i;
			}
		}
		return -1;
	}
	/*BCdecide & locBCD(BCdFWD);
	if ( !isPair1 ) {
	locBCD = BCdREV;
	}*/
	int start(-1), stop(-1);
	int idx(-1);
	int scanRegion = 4; //dna region to scan for Tag Seq
	if (!d->getTA_cut() && isPair1) {//no technical adapter found / given by user: scan wider region for barcode
		scanRegion = 14; //arbitary value
	}
	if (d->isMIDseq()) {
		if (d->length()<MinTagLen) { return -1; }
		scanRegion = d->length() - MinTagLen + 1;
	}

	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
	if (!BCdFWDREV[!isPair1].b_BCdirFix) {
		if (start == -1) {//check reverse transcription
						  //d->reverseTranscribe();
			scanBC_rev(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
			if (start != -1) {
				BCdFWDREV[!isPair1].BCrevhit++;
			}
		}
		else {
			BCdFWDREV[!isPair1].BChit++;
		}
		//check if BC direction can be fixed
		if (BCdFWDREV[!isPair1].BCrevhit + BCdFWDREV[!isPair1].BChit > DNAinMemory) {
			if (!eval_reversingBC(isPair1)) { return -1; }
		}
	}
	if (start != -1) {
		if (BcutTag && !d->isMIDseq()) {
			//remove tag from DNA
			d->cutSeq(start, stop);
			d->setBarcodeCut();
		}
	}
	else {
		idx = -1;
	}
	return idx;
}
int Filters::findTag(shared_ptr<DNA> d, string&presentBC, int& c_err, bool isPair1) {
	if (bDoHeadSmplID) {
		for (unsigned int i = 0; i<HeadSmplID.size(); i++) {
			size_t pos = d->getOldID().find(HeadSmplID[i]);
			if (pos != string::npos) {
				SampleIntoHead(i, d, pos);
				return i;
			}
		}
		return -1;
	}
	/*BCdecide & locBCD(BCdFWD);
	if ( !isPair1 ) {
	locBCD = BCdREV;
	}*/
	int start(-1), stop(-1);
	int idx(-1);
	int scanRegion = 14; //dna region to scan for Tag Seq

	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
	if (!BCdFWDREV[!isPair1].b_BCdirFix) {
		if (start == -1) {//check reverse transcription
						  //d->reverseTranscribe();
			scanBC_rev(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
			if (start != -1) {
				BCdFWDREV[!isPair1].BCrevhit++;
			}
		}
		else {
			BCdFWDREV[!isPair1].BChit++;
		}
		//check if BC direction can be fixed
		if (BCdFWDREV[!isPair1].BCrevhit + BCdFWDREV[!isPair1].BChit > DNAinMemory) {
			if (!eval_reversingBC(isPair1)) { return -1; }
		}
	}
	if (start == -1) {
		idx = -1;
	}
	return idx;
}
int Filters::cutTag(shared_ptr<DNA> d, bool isPair1) {

	if (d->length() < MinTagLen ) { return -1; }
	if ((isPair1 && !bDoBarcode)||
		(!isPair1 && !bDoBarcode2)){
		d->setBCnumber(0, BCoffset);
		return BCoffset; //not failed, just not requested
	}
	/*BCdecide & locBCD(BCdFWD);
	if ( !isPair1 ) {
		locBCD = BCdREV;
	}*/

	int idx(-1);
	if (bDoHeadSmplID){
		for (unsigned int i=0; i<HeadSmplID.size(); i++){
			size_t pos = d->getOldID().find(HeadSmplID[i]);
			if (pos != string::npos){			
				if ( !bDoBarcode2 ) { SampleIntoHead(i, d, pos); }//this has to be done AFTER two BCs are read (on a higher lvl)
				else {d->setBCnumber(i, BCoffset);}
				return i;
			}
		}
		return idx;
	}
	else if (bOneFileSample){
		d->setBCnumber(0,this->currentBCnumber());
		BCintoHead(0, d, "FileName", isPair1, false);
		return 0;
	}
	int start(-1),stop(-1);
	string presentBC(""); int c_err(0);
	int scanRegion=4; //dna region to scan for Tag Seq
	if (!d->getTA_cut() && isPair1){//no technical adapter found / given by user: scan wider region for barcode
		scanRegion = 14; //arbitary value
	}
	if (d->isMIDseq()){
		scanRegion = d->length() - MinTagLen+1;
	}
	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
	if ( !BCdFWDREV[!isPair1].b_BCdirFix ) {
		if (start == -1){//check reverse transcription
			//d->reverseTranscribe();
			scanBC_rev(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
			if (start!=-1){
				BCdFWDREV[!isPair1].BCrevhit++;
			}
		} else {
			BCdFWDREV[!isPair1].BChit++;
		}
		//check if BC direction can be fixed
		if ( BCdFWDREV[!isPair1].BCrevhit + BCdFWDREV[!isPair1].BChit > 5000 ) {
			eval_reversingBC(isPair1);//){return -1;}
		}
	}
	if (start < 0) {
		return (-1);
	}
	d->setBCnumber(idx, BCoffset);

	if (BcutTag && !d->isMIDseq()) {
		//remove tag from DNA
		d->cutSeq(start, stop);
		d->setBarcodeCut();
		BCintoHead(idx, d, presentBC, c_err, isPair1);
	}
	return idx;
}

void Filters::BCintoHead(int idx, shared_ptr<DNA> d,const string presentBC,
						 const int c_err, bool pair1, bool atEnd){
	vector<string> * locBC = (!pair1) ? &Barcode2 : &Barcode;
	string on (d->getID());
	string spacee (" ");
	//keep only until first space
	//remove">"
	if (!atEnd){
		on = on.substr(0,on.find_first_of(' ',0));
	}

	string nID ( SampleID[idx]);
	if (bDoCombiSamples){
		nID = SampleID_Combi[idx];
	}
	nID += iniSpacer +
		on + spacee + string("orig_bc=") + (*locBC)[idx];
	if (c_err > 0 || atEnd) {//atEnd: dbl barcode
		//convert c_err to s_c_err;
		string s_c_err; stringstream conve;
		conve << c_err;
		s_c_err = conve.str();
		nID += spacee + string("new_bc=") + presentBC +
			spacee + string("bc_diffs=") + s_c_err;
	}
	d->setNewID(nID);
	d->setBCnumber(idx, BCoffset);
}

void Filters::SampleIntoHead(const int idx, shared_ptr<DNA> d, const size_t pos){
	string on = d->getID(),spacee=" ";
	size_t pos2 = on.find_first_of(" ",pos);
	string on2 = on.substr(0,pos)+on.substr(pos2+1);
	string on3 = on.substr(pos,pos2);

	string nID(SampleID[idx]);
	if (bDoCombiSamples){
		nID = SampleID_Combi[idx];
	}

	nID +=  iniSpacer +
		on2 + spacee + string("orig_hdPart=")+on3;
	d->setNewID(nID);
	d->setBCnumber(idx, BCoffset);
}
void Filters::countBCdetected(int BC, int Pair, bool MidQ) {
	if (!bDoMultiplexing) { return; }
	if (BC - BCoffset < 0) {
		cerr << "BC below 0("<< BC - BCoffset <<") pair "<< Pair<<" in countBCdetected\n"; exit(132);
	}
	if (!MidQ) {
		if (Pair < 0) { Pair = 0; }
		colStats[Pair].BarcodeDetected[BC - BCoffset]++;
	}
	else {
		if (Pair <= 0) {
			statAddition.BarcodeDetected[BC - BCoffset]++;
		}
	}
}
bool Filters::eval_reversingBC(bool fwd){
	if ( !fwd && !bDoBarcode2 ) {
		return true;
	}
	/*BCdecide lbcd(BCdFWD);
	if ( !fwd ) {
		lbcd = BCdREV;
	}*/
	if ( BCdFWDREV[!fwd].b_BCdirFix ) { return true; }
	BCdFWDREV[!fwd].b_BCdirFix = true; lMD->setBCfixed(true, fwd);
	if ( BCdFWDREV[!fwd].BCrevhit> BCdFWDREV[!fwd].BChit * 8 ) {//use reversed BC ..
		BCdFWDREV[!fwd].reversedBCs = true;
		if ( fwd ) {
			reverseTS_all_BC();
		} else {
			reverseTS_all_BC2();
		}
		if ( BCdFWDREV[!fwd].BChit > 0 ) {
			restartSet=true;
			return false;
		}
	} else if ( BCdFWDREV[!fwd].BCrevhit>0 ) {
		restartSet=true;
		return false;
	}
	return true;
}
void Filters::scanBC_rev(shared_ptr<DNA> d,int& start,int& stop,int& idx,int c_err, 
					 int scanRegion,string & presentBC,
					 bool fwdStrand) {
	vector<string> emptyV(0), emptyV2(0);
	vector<string>& locRevBC(emptyV2);
	vector<string>&  locBC(emptyV);
	if ( !fwdStrand ) {
		locRevBC = revBarcode2;
		locBC = Barcode2;
	} else {
		locRevBC = revBarcode;
		locBC = Barcode;
	}
	int BCs = (int) locRevBC.size();
	if (BCs==0){
		locRevBC = locBC;
		BCs = (int)locRevBC.size();
		for (int i=0; i< BCs; i++){
			reverseTS(locRevBC[i]);
		}
		if ( !fwdStrand ) {//copy over BC
			revBarcode2 = locRevBC;
		} else {
			revBarcode = locRevBC;
		}
	}
	//check each possible BC for a match
	if (TagErrs == 0){
		for (idx=0; idx< BCs; idx++){
			start = d->matchSeq_tot(locRevBC[idx],0,scanRegion,c_err);
			if (start!=-1){
				presentBC = locRevBC[idx];
				stop = start+ (int)locRevBC[idx].length();
				break;
			}
		}
	} else {
		vector<int> stars(0),idxses(0);
		bool zeroErr = false;
		//this version tries all BC's and if there are more than one possible match, will reject all matches
		for (idx=0; idx< BCs; idx++){
			start = d->matchSeq_tot(locRevBC[idx],TagErrs,scanRegion,c_err);
			if (start!=-1){
				if (c_err==0){
					stop = start+ (int) locRevBC[idx].length();
					presentBC = locRevBC[idx];
					zeroErr = true;
					break;
				}
				stars.push_back(start);
				idxses.push_back(idx);
			}
		}
		if (!zeroErr && stars.size()>0){
			//int pair = (int)!fwdStrand;//d->getReadMatePos();
			if (stars.size() > 1){//too many matches, thus true seq can't be found
				//currently only have only one BC, could be changed in future
				//sTagNotCorrected(pair);
				d->QualCtrl.fail_correct_BC = true;
				idx=-1; start = -1;
				return;
			}
			d->QualCtrl.suc_correct_BC = true;
			//sTagCorrected(pair);// colStats.suc_correct_BC++;
			start = stars[0];
			idx = idxses[0];
			stop = start+(int)locRevBC[idx].length();
			presentBC = d->getSubSeq(start,stop);
		}

	}
}
void Filters::scanBC(shared_ptr<DNA> d, int& start, int& stop, int& idx, int c_err,
	int scanRegion, string & presentBC, bool fwdStrand) {
	//check each possible BC for a match
	//TODO: check for BC using suffix tree
	vector < string > emptyV(0);
	//vector<string> * locBC = (!fwdStrand) ? &Barcode2 : &Barcode;
	vector<int> * locBClen = (!fwdStrand) ? &Barcode2_len : &Barcode_len;
	//int pair = d->getReadMatePos();
	//int pair = int( !fwdStrand);

	BarcodeList & locBCL(emptyBCs);
	uint MTL (MaxTagLen);
	uint ITL (MinTagLen);
	if ( !fwdStrand ) {
		locBCL = BCList2;
		MTL = MaxTagLen2;
		ITL = MinTagLen2;
	} else {
		locBCL = BCList;
	}

	for (start = 0; start < scanRegion; start++){
		string test = d->getSubSeq(start, MTL);
		std::unordered_map<string, int>::iterator itFnd = locBCL.find(test);
		if ( itFnd != locBCL.end() ) {
			idx = (*itFnd).second;
			stop = start + (int)(*locBClen)[idx];
			presentBC = test; // (*locBC)[idx];
			return;
		}
	}
	start = -1;
	if ( TagErrs != 0 || ITL != MTL ) {
		vector<int> stars(0), idxses(0);
		bool zeroErr = false;
		
		//this version tries all BC's and if there are more than one possible match, will reject all matches
		for (auto jx = locBCL.begin(); jx != locBCL.end();jx++) {
			start = d->matchSeq_tot((*jx).first, TagErrs, scanRegion, c_err);
			if (start!=-1){
				idx = (*jx).second;
				if (c_err==0){
					stop = start + (int)(*locBClen)[idx];
					presentBC = d->getSubSeq(start, MTL); // (*locBC)[idx];
					zeroErr = true;
					break;
				}
				stars.push_back(start);
				idxses.push_back(idx);
			}
		}
		if (!zeroErr && stars.size()>0){
			if (stars.size() > 1){//too many matches, thus true seq can't be found
				//sTagNotCorrected(pair);
				d->QualCtrl.fail_correct_BC = true;
				idx = -1; start = -1;
				return;
			}
			d->QualCtrl.suc_correct_BC = true;
			//sTagCorrected(pair);// colStats.suc_correct_BC++;

			start = stars[0];
			idx = idxses[0];
			stop = start + (int)(*locBClen)[idx];
			presentBC = d->getSubSeq(start,stop);
		}

	}
}
//cuts primers, tags
bool Filters::cutPrimer(shared_ptr<DNA> d,int primerID,bool RC,int pair){
	//only adapted to singular BC
	if (PrimerL[0].length()==0){return true;}
	int start(-1) ,stop(-1);
	int tolerance(30), startSearch(0);
	if (!d->getBarcodeCut() && MaxTagLen > 0) { tolerance = MaxTagLen + 4; 
	} else { tolerance = 22; }//in this case nothing is known about 5' end

	if (!BcutTag){
		//Tag was not cut out of Seq, take this into account
		startSearch = MinTagLen-2;
		tolerance += (MaxTagLen-MinTagLen)+4;
	}
	if (!RC) {
		start = d->matchSeq(PrimerL[primerID], PrimerErrs, tolerance, startSearch);
		stop = start + (int)PrimerL[primerID].length();
	} else {
		int QS = d->length();int limit = max(QS >> 1, QS - 150); stop = QS;
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit, startSearch);
	}
	if (start == -1){//failed to match primer
		d->QualCtrl.PrimerFail = true;
		//sPrimerFail(pair);// max(0, (int)d->getReadMatePos()));
		if (alt_PrimerErrs!= 0 && PrimerErrs < alt_PrimerErrs){
			if (!RC) {
				start = d->matchSeq(PrimerL[primerID], alt_PrimerErrs, tolerance, startSearch);
				stop = start + (int)PrimerL[primerID].length();
			} else {
				int QS = d->length(); int limit = max(QS >> 1, QS - 150); stop = QS;
				start = d->matchSeqRev(PrimerL_RC[primerID], alt_PrimerErrs, limit, startSearch);
			}
		}
		if (start== -1){
			//statAddition.PrimerFail++;
			return false;
		}else if (pair!=1){ //2nd read shouldnt be affected by fwd primer (but still checked in short read mode)
			d->setMidQual(true);
		}
	}
	if (!BcutPrimer){
		if ( !RC ) { d->cutSeq(0, start); }
		else { d->cutSeq(stop,-1); }
		return true;
	}
	
	//remove the primer, if confimed before
	if ( !RC ) { d->cutSeq(0, stop); }
	else { d->cutSeq(start, -1); }
	d->setFwdPrimCut();
	return true;
}
bool Filters::findPrimer(shared_ptr<DNA> d, int primerID, bool RC, int pair) {
	//only adapted to singular BC
	if (PrimerL[0].length() == 0) { return true; }
	int start(-1);// , 
	//int stop(-1);
	int tolerance(22), startSearch(0);
	if (!d->getBarcodeCut() && MaxTagLen > 0) {
		tolerance = MaxTagLen + 4;
		startSearch = MinTagLen-4;
	}
	else { tolerance = 16; }//in this case nothing is known about 5' end
	if (!RC) {
		start = d->matchSeq(PrimerL[primerID], PrimerErrs, tolerance, startSearch);
		//stop = start + (int)PrimerL[primerID].length();
	}
	else {
		int QS = d->length(); int limit = max(QS >> 1, QS - 150); //stop = QS;
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit, startSearch);
	}
	if (start == -1) {//failed to match primer
		return false;
	}
	return true;
}
bool Filters::cutPrimerRev(shared_ptr<DNA> d,int primerID,bool RC){
	//const string& se = d->getSeq();
	int start(-1) ,stop(d->length());
	int QS = d->length(); 
	int limit=max(QS>>1,QS-150);
	//int limit = QS>>1;

	if (!RC) {
		start = d->matchSeq(PrimerR[primerID] , PrimerErrs, 15,0);
		stop = start + (int)PrimerR[primerID].length();
	} else {
		start = d->matchSeqRev(PrimerR_RC[primerID] , PrimerErrs, limit);
	}
	


	if (start == -1){//failed to match primer
		return false;
	} 

	if ( !BcutPrimer ) { //found it, but no cut
		if ( !RC ) { d->cutSeq(0, start); }
		else { d->cutSeq(stop,-1); }
		return true; 
	}

	//remove the primer, if confimed before
	if ( !RC ) {
		d->cutSeq(0, stop);//start  everything in front has to be removed
	} else {
		d->cutSeq(start, -1); // everything in the end has to be removed
	}
	//string neSe = se.substr(0,start) + se.substr(stop);
	d->setRevPrimCut();

	return true;
}
bool Filters::readMap(OptContainer& cmdArgs){

	if (cmdArgs.find("-map")  == cmdArgs.end()){
		this->noMapMode(cmdArgs);
		return true;
	}

	string MapF = cmdArgs["-map"];

	string path = ""; bool pathMode = false;
	if (cmdArgs.find("-i_path")  != cmdArgs.end() && cmdArgs["-i_path"].length() > 2){
		path=cmdArgs["-i_path"] + string("/");
		pathMode = true;//check later if mapping file contains fasta/fastq
	}

	MinTagLen = 100000; MinTagLen2 = 1000000; MaxTagLen = 0; MaxTagLen2 = 0; MinPrimLen = 100000;
	string line;
	ifstream in(MapF.c_str());
	if (!in){
		cerr<<"Could not find "<<MapF<<" mapping file. Exiting.\n";exit(2);
	}
	int ini_ColPerRow(0),cnt(0),skips(0);

	//check MAP format
	//while(getline(in,line,'\n')) {
	while(!safeGetline(in,line).eof()) {
		if(line.substr(0,1) == "#"){skips++;continue;}
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss,segments,'\t')) {
			ColsPerRow++;
		}
		if (segments == "") { ColsPerRow++; }

		if (cnt==0){
			ini_ColPerRow = ColsPerRow;
		} else {
			if (ColsPerRow != ini_ColPerRow){
				cerr<<"Number of columns on line "<<cnt+skips<<" is "<<ColsPerRow<<". Expected "<<ini_ColPerRow<<" columns.\n";
				return false;
			}
		}
		cnt++;
	}
	if (ini_ColPerRow==0){
		cerr<<"Mapping File exists, but appears to be badly formated (0 columns detected). Exiting\n";exit(2);
	}
	if (cnt==0){
		cerr<<"Mapping File exists, but appears to be badly formated (0 lines detected). Exiting\n";exit(2);
	}
	

	PrimerIdx.resize(cnt,-1); PrimerIdxRev.resize(cnt,-1);
	Barcode.resize(cnt, "");
	Barcode2.resize(cnt, "");
	SampleID.resize(cnt, "");
	SampleID_Combi.resize(cnt, "");
	Barcode_len.resize(cnt, 0);
	Barcode2_len.resize(cnt, 0);

	colStats[0].BarcodeDetected.resize(cnt, 0);
	colStats[1].BarcodeDetected.resize(cnt, 0);
	colStats[0].BarcodeDetectedFail.resize(cnt, 0);
	colStats[1].BarcodeDetectedFail.resize(cnt, 0);
	HeadSmplID.resize(cnt, "");
	statAddition.BarcodeDetected.resize(cnt, 0);
	statAddition.BarcodeDetectedFail.resize(cnt, 0);
	hetPrimer[0].resize(cnt, ""); hetPrimer[1].resize(cnt, "");
	in.clear();
	in.seekg(0, ios::beg);   

	//extract MAP content
	cnt=0;
	bool hasQualityColumn=false;
	vector<string> terms(15);terms[0]="SampleID";
	terms[1] = "BarcodeSequence"; terms[2]="LinkerPrimerSequence";
	terms[3] = "ReversePrimer"; terms[4] = "fastqFile";
	terms[5] = "fnaFile"; terms[6] = "qualFile";
	terms[7] = "SampleIDinHead"; terms[8] = "MIDfqFile";
	terms[9] = "CombineSamples"; terms[10] = "ForwardPrimer";
	terms[11] = "Barcode2ndPair";
	terms[12] = "HetSpacerFwd";	terms[13] = "HetSpacerRev";
	terms[14] = "derepMin";
	bool hetOneSide = false;

	vector<int> termIdx(terms.size(),-1);

	while(!safeGetline(in,line).eof()) {
//	while(getline(in,line,'\n')) {
		if(cnt!=0 && line.substr(0,1) == "#"){continue;}
		if (line.length()<10){continue;}
		if (cnt==0){
			line = line.substr(1);
		}

		string segments;
		stringstream ss;
		ss << line;
		int tbcnt=0;

		//cmdArgs["-i_MID_fastq"]
		while (getline(ss,segments,'\t')) {
			trim(segments);
			if (cnt==0){ //search for header
				//Primer, BarcodeSequence, LinkerPrimerSequence
				//PrLCol(-1), PrRCol(-1), BCCol(-1), SIDCol(-1);
				for (unsigned int i=0; i<terms.size();i++){
					if (segments == terms[i]){
						if (i==6){hasQualityColumn=true;}
						if (i == 12){ if (hetOneSide){ doHetPrimerExplicit = true; } else{ hetOneSide = true; } }
						if (i == 13){ if (hetOneSide){ doHetPrimerExplicit = true; } else{ hetOneSide = true; } }
						if (termIdx[i] != -1){
							cerr<<"Header contains ambiguous entries: "<<segments<<"\n";
							exit(9);
						}
						termIdx[i] = tbcnt;
					}
				}
			} else { //read data into entries
				for (uint k=0; k<terms.size();k++){
					if (termIdx[k]==tbcnt){
						extractMap(k, cnt - 1, tbcnt, segments, hasQualityColumn);
					}
				}
			}
			tbcnt++;
		}
		cnt++;
	}
	//check some prerequisites 
	if (HeadSmplID.size() ==0 && Barcode.size() ==0){
		cerr<<"Could not find in mapping file either (1) valid Barcodes or (2) valid SampleID's in Sequence header. Exiting\n"; 
		exit(2);
	}
	if (pathMode) {
		if (termIdx[5] == -1 && termIdx[4] == -1) {
			cerr << "The defined input through a directory (-i_path) requries either \"fnaFile\" or \"fastqFile\" columns in the mapping file. \nAborting..\n";
			exit(55);
		}
	}

	decideHeadBC();
	
	//check for duplicate barcodes, but only if no list provided
	if(termIdx[4]!=-1 &&termIdx[5]!=-1 &&termIdx[6]!=-1){
		checkDoubleBarcode();
	}
	if (termIdx[9] != -1){ bDoCombiSamples = true; }
	checDoubleSampleID();

	//eleminate required checks for primers, if there is simply no primer given
	if (PrimerL.size() == 0 || PrimerL[0].length() == 0 || PrimerL[0].substr(0, 1) == " ") {
		BcutPrimer = false;
		bRequireFwdPrim = false;
		alt_bRequireFwdPrim = false;
	}
	if (PrimerR.size() == 0 || PrimerR[0].length() == 0 || PrimerR[0].substr(0, 1) == " ") {
		alt_bRequireRevPrim = false;
		bRequireRevPrim = false;
	}

	this->BarcodePreStats();

	return true;
}
void Filters::decideHeadBC(){
	bDoMultiplexing=true;
	if (HeadSmplID[0].length()>0 && Barcode[0].length()==0){
		bDoBarcode = true; bDoHeadSmplID = true; MinPrimLen = 0; return;
	} else if ( HeadSmplID[0].length() == 0 && Barcode[0].length() > 0 && Barcode2[0].length() > 0 ) {
		bDoBarcode = true; bDoHeadSmplID = false; bDoBarcode2 = true;  return;
	} else if ( HeadSmplID[0].length() == 0 && Barcode[0].length() > 0 ) {
		bDoBarcode = true; bDoHeadSmplID = false; bDoBarcode2 = false;  return;
	} else if ( HeadSmplID[0].length() == 0 && Barcode[0].length() == 0 ) {
		//simply check if each filename is different
		if (FastqF.size() > 0 ){
			for (uint i = 0; i<FastqF.size(); i++){
				for (uint j = i+1; j<FastqF.size(); j++){
					if (FastqF[i] == FastqF[j]){
						cerr << "File names " << i << " and " << j << " are equal - no identification by filename possible.\n   Aborting..\n"; exit(55);
					}
				}
			}
			bOneFileSample = true; bDoBarcode = true; bDoHeadSmplID = false;
			return;
		}
		if (FastaF.size() > 0){
			for (uint i = 0; i < FastaF.size(); i++){
				for (uint j = i + 1; j < FastaF.size(); j++){
					if (FastaF[i] == FastaF[j]){
						cerr << "File names " << i << " and " << j << " are equal - no identification by filename possible.\n   Aborting..\n"; exit(55);
					}
				}
			}

			bOneFileSample = true; bDoBarcode = true; bDoHeadSmplID = false;
			return;
		}
	}

	bDoMultiplexing = false;
	cerr << "No Barcode and no ID in header defined.. aborting\n";
	exit(53);

}


void Filters::checkDoubleBarcode(){
	if (!bDoBarcode || bDoHeadSmplID || Barcode.size()==0){ return; }
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	if ( bDoBarcode2 ) {
		if ( Barcode.size() != Barcode2.size() ) {
			cerr << "Unequal Barcode vector sizes in dual barcoding controls. Exiting.."; exit(45);
		}
		for ( unsigned int i = 0; i < Barcode.size(); i++ ) {
			for ( unsigned int j = i + 1; j < Barcode.size(); j++ ) {
				if ( strcmp(Barcode[i].c_str(), Barcode[j].c_str()) == 0 && strcmp(Barcode2[i].c_str(), Barcode2[j].c_str()) == 0 ) {
					empty[0] = i; empty[1] = j ; doubles.push_back(empty);
				}
			}
		}
		if (doubles.size() > 0){
			for (uint x = 0; x < doubles.size(); x++){
				int i = doubles[x][0]; int j = doubles[x][1];
				cerr << "Duplicate dual Barcode detected: Barcode1 " << i + 1 << " (" << Barcode[i] << ") and " << j + 1 << " (" << Barcode[j] << ")  as well as Barcode1 " << i + 1 << " (" << Barcode2[i] << ") and " << j + 1 << " (" << Barcode2[j] << ") are equal.\n";
			}
			exit(8);
		}
	}
	else {
		for ( unsigned int i = 0; i < Barcode.size(); i++ ) {
			for ( unsigned int j = i + 1; j < Barcode.size(); j++ ) {
				if ( strcmp(Barcode[i].c_str(), Barcode[j].c_str()) == 0 ) {
					empty[0] = i; empty[1] = j; doubles.push_back(empty);
				}
			}
		}
		if (doubles.size() > 0){
			for (uint x = 0; x < doubles.size(); x++){
				int i = doubles[x][0]; int j = doubles[x][1];
				cerr << "Duplicate Barcode detected: Barcode " << i + 1 << " (" << Barcode[i] << ") and " << j + 1 << " (" << Barcode[j] << ") are equal.\n";
			}
			exit(8);
		}

	}
}
void Filters::checkDoubleSampleIDHead(){
	if (!bDoHeadSmplID){return;}
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	for (unsigned int i=0; i<HeadSmplID.size();i++){
		for (unsigned int j=i+1;j<HeadSmplID.size();j++){
			if(strcmp(HeadSmplID[i].c_str(),HeadSmplID[j].c_str()) == 0){
				empty[0] = i; empty[1] = j; doubles.push_back(empty);
			}
		}
	}
	if (doubles.size() > 0){
		for (uint x = 0; x < doubles.size(); x++){
			int i = doubles[x][0]; int j = doubles[x][1];
			cerr << "Duplicate Header2split detected: pattern " << i + 1 << " and " << j + 1 << " are equal.\n";
		}
		exit(8);
	}

}

void Filters::checDoubleSampleID(){
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	for (unsigned int i = 0; i<SampleID.size(); i++){
		for (unsigned int j=i+1;j<SampleID.size();j++){
			if(strcmp(SampleID[i].c_str(),SampleID[j].c_str()) == 0){
				empty[0] = i; empty[1] = j; doubles.push_back(empty);
			}
		}
	}
	if (doubles.size() > 0){
		for (uint x = 0; x < doubles.size(); x++){
			int i = doubles[x][0]; int j = doubles[x][1];
			cerr << "Duplicate SampleID detected: SampleID " << i + 1 << " and " << j + 1 << " are equal.\n";
		}
		exit(8);
	}
	if (!bDoCombiSamples){
		return;
	}
	bDoCombiSamples = false;
	if (SampleID_Combi.size() <= 1){
		return;
	}
	
	string prevCSID = SampleID_Combi[0];
	for (unsigned int i = 1; i < SampleID_Combi.size(); i++){
		if (SampleID_Combi[i] != prevCSID){
			bDoCombiSamples = true; 
		}
		if (SampleID_Combi[i] == ""){
			SampleID_Combi[i] = SampleID[i];
		}
	}
}

void Filters::BarcodePreStats(){
	MinTagLen=100000;MaxTagLen=0;
	for (unsigned int i=0; i<Barcode.size();i++ ){
		if ( Barcode[i].length()<MinTagLen) MinTagLen= (unsigned int) Barcode[i].length();
		if ( Barcode[i].length()>MaxTagLen) MaxTagLen= (unsigned int) Barcode[i].length();
		//initialize Barcodes
		BCList[Barcode[i]] = i;
		Barcode_len[i] = (int) Barcode[i].length();
	}
	if (MinTagLen == MaxTagLen){
		bBarcodeSameSize = true;
	}
	MinTagLen2 = 100000; MaxTagLen2 = 0;
	for ( unsigned int i = 0; i < Barcode2.size(); i++ ) {
		//create index
		BCList2[Barcode2[i]] = i;
		if ( Barcode2[i].length()<MinTagLen2 ) MinTagLen2 = (unsigned int)Barcode2[i].length();
		if ( Barcode2[i].length()>MaxTagLen2 ) MaxTagLen2 = (unsigned int)Barcode2[i].length();
		Barcode2_len[i] = (int)Barcode2[i].length();
	}
	//unique_ptr<dualPrimerDistrStats> dPDS = make_unique<dualPrimerDistrStats>(Barcode, Barcode2)
	dPDS = make_shared<dualPrimerDistrStats>(Barcode, Barcode2);
	dHDS = make_shared<dualPrimerDistrStats>(hetPrimer[0], hetPrimer[1]);
	//fix empty last column specifically for derepMinNum
	if (derepMinNum.size() > 0) {
		derepMinNum.resize(Barcode.size(), -1);
	}

}
void Filters::resetStats(){
	statAddition.reset();
	PreFiltP1->reset(); PreFiltP2->reset();
	dPDS->reset(); dHDS->reset();
	for (size_t i = 0; i < 2; i++) {
		RepStat[i]->reset(); RepStatAddition[i]->reset(); colStats[i].reset();
	}

}

void Filters::failedStats2(shared_ptr<DNA> d,int pair){
	int pa = max(pair, 0);
	if (bDoMultiplexing){
		int idx = d->getBCnumber() - BCoffset;
		if ( bOneFileSample ) {
			colStats[pa].BarcodeDetectedFail[0]++;
		} else if (idx >= 0) {

			
#ifdef DEBUG
			if (pa < 0 || pa>1) { cerr << "Pair in failedStats2 set to:" << pa << endl; }
			if (idx >= (int)colStats[pa].BarcodeDetectedFail.size()) {
				cerr << "idx in failedStats2 too big:" << idx << endl;
			}
#endif // DEBUG
				colStats[pa].BarcodeDetectedFail[idx]++;
		}
	}

}
void Filters::prepStats() {
	float remSeqs = float(colStats[0].total - colStats[0].totalRejected);
	RepStat[0]->calcSummaryStats(remSeqs, min_l, min_q);
	RepStat[1]->calcSummaryStats(remSeqs, min_l, min_q);
	if (bAdditionalOutput) {
		remSeqs = float(statAddition.total - statAddition.totalRejected);
		RepStatAddition[0]->calcSummaryStats(remSeqs, min_l, min_q);
		RepStatAddition[1]->calcSummaryStats(remSeqs, min_l, min_q);
	}
	PreFiltP1->calcSummaryStats(1, min_l, min_q);
	PreFiltP2->calcSummaryStats(1, min_l, min_q);
}


void Filters::addPrimerL(string segments, int cnt){
	int used = -1;
	for (unsigned int i=0; i<PrimerL.size();i++){
		if (segments==PrimerL[i]){used=i;}
	}
	if (used == -1){
		PrimerL.push_back(segments);
		PrimerL_RC.push_back(reverseTS2(segments));
		PrimerIdx[cnt] = (int)PrimerL.size() - 1;
		if ( segments.length()<MinPrimLen) MinPrimLen= (unsigned int) segments.length();
	} else {
		PrimerIdx[cnt] = used;
	}
}
void Filters::addPrimerR(string segments, int cnt){
	bPrimerR=true;
	int used = -1;
	for (unsigned int i=0; i<PrimerR.size();i++){
		if (segments==PrimerR[i]){used=i;}
	}
	if (used == -1){
		PrimerR.push_back(segments);
		PrimerR_RC.push_back(reverseTS2(segments));
		PrimerIdxRev[cnt] = (int)PrimerR.size() - 1;
	} else {
		PrimerIdxRev[cnt] = used;
	}
}


void Filters::extractMap(int k, int cnt, int tbcnt, string & segments,
	bool hasQualityColumn){
	//terms[4] = "fastqFile";terms[5] = "fnaFile"; terms[6] = "qualFile";
	switch(k)
	{
	case 4: // fastq file
		trim(segments);
		FastqF.push_back(segments);
		break;
	case 5: // fna file
		trim(segments);
		FastaF.push_back(segments);
		if (!hasQualityColumn  ){//code to create artificial quality file name
				string newQ = segments;
				size_t pos = newQ.find_last_of(".");
				newQ = newQ.substr(0,pos);
				newQ += string(".qual");
				QualF.push_back(newQ);
		}
		break;
	case 6: // qual file
		trim(segments);
		QualF.push_back(segments);
		break;
	case 2: //left primer 		
	case 10:
		trim(segments);
		transform(segments.begin(), segments.end(),segments.begin(), ::toupper);
		this->addPrimerL(segments,cnt);
		break;
	case 3: // right primer
		trim(segments);
		transform(segments.begin(), segments.end(),segments.begin(), ::toupper);
		this->addPrimerR(segments,cnt);
		break;
	
	case 0: //ID
		trim(segments);
		SampleID[cnt] = segments;
		break;

	case 1: //Barcode
		trim(segments);
		transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
		Barcode[cnt] = segments;
		break;
	case 11: //Barcode rev
		trim(segments);
		transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
		Barcode2[cnt] = segments;
		break;
	case 7: //sample id in head
		trim(segments);
		HeadSmplID[cnt] = segments;
		break;
	case 8://mid xtra fq
		trim(segments);
		MIDfqF.push_back(segments);
		break;
	case 9://combine samples
		trim(segments);
		SampleID_Combi[cnt] = segments;
		break;
	case 14://demultiplex num
		if (segments.length() == 0) {
			derepMinNum.push_back(-1);
		}else if (!is_digits(segments)){
			cerr << "Wrong map entry \"" << segments << "\". For header derepMin only number can be used.\n"; exit(313);
		} else {
			int nint = atoi(segments.c_str());
			derepMinNum.push_back( nint);
		}
		break;
	case 12://het primer fw
		if (!doHetPrimerExplicit){break;}
		hetPrimer[0][cnt] = segments;
	case 13://het primer rv
		if (!doHetPrimerExplicit){ break; }
		hetPrimer[1][cnt] = segments;
	}

	if (k==6 || k==5){//a qual pushback was "" (was empty); replace
		if (QualF.size() == FastaF.size() && QualF.back()==""){
				string newQ = FastaF.back();
				size_t pos =  newQ.find_last_of(".");
				newQ = newQ.substr(0,pos);
				newQ += string(".qual");
				QualF.back() = newQ;
		}
	}

}

void Filters::printHisto(ostream& give,int which, int set){
	bool p2stat = pairedSeq > 1 ;

	if (set == 0) {
		vector<uint> colStats(RepStat[0]->get_rstat_Vmed(which));
		vector<size_t> ra( RepStat[0]->getVrange(which) );

		if (which == 1) {
			give << "Qual\tFilterObs" << endl;
		} else {
			give << "Length\tFilterObs" << endl;
		}
		for (size_t i = ra[0]; i < ra[1] ; i++) {
				give << i << "\t" << colStats[i] << endl;
		}
	} else if (set == 1) {
		vector<size_t> ra(2,0), tra;
		vector<uint> stat;
		vector<bool> skips(6, false);
		vector<vector<uint>> matHist;
		if (which == 1) {	give << "#Qual\t";
		} else { give << "#Length\t"; }
		if (p2stat && b_doFilter) { give << "FilteredP1\tFilteredP2\t"; } 
		else if ( b_doFilter){ give << "Filtered\t"; skips[1] = true; }
		else { skips[1] = true; skips[0] = true; }
		
		if ( bAdditionalOutput && b_doFilter ) {
			if (p2stat) { give << "AddFilterP1\tAddFilterP2\t"; }
			else { give << "AddFilter\t"; skips[3] = true; }
		} else {
			skips[2] = true; skips[3] = true;
		}
		if (p2stat) { give << "RawReadsP1\tRawReadsP2"; }
		else { give << "RawReads"; skips[5] = true; }
		give << endl;
		ra = RepStat[0]->getVrange(which);
		if (p2stat) { tra = RepStat[1]->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		if (!skips[2]) {
			tra = RepStatAddition[0]->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]);
			if (p2stat) { tra = RepStatAddition[1]->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		} 
		tra = PreFiltP1->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]);
		if (p2stat) { tra = PreFiltP2->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		vector<uint> empt(ra[1], 0);
		matHist = vector<vector<uint>>(6, empt);
		for (size_t kk = 0; kk < 6; kk++) {
			if (skips[kk]) { continue; }
			switch (kk) {
			case 0:	stat = RepStat[0]->get_rstat_Vmed(which); break;
			case 1:	stat = RepStat[1]->get_rstat_Vmed(which); break;
			case 2:	stat = RepStatAddition[0]->get_rstat_Vmed(which); break;
			case 3:	stat = RepStatAddition[1]->get_rstat_Vmed(which); break;
			case 4:	stat = PreFiltP1->get_rstat_Vmed(which); break;
			case 5:	stat = PreFiltP2->get_rstat_Vmed(which); break;
			}
			for (size_t i = 0; i < stat.size(); i++) {
				if (i>=ra[1]) {break;}
				matHist[kk][i] = stat[i];
			}
		}
		for (size_t i = ra[0]; i < ra[1] ; i++) {
			give << i ;
			for (size_t kk = 0; kk < matHist.size(); kk++) {
				if (skips[kk]) { continue; }
				give << "\t" << matHist[kk][i];
			}
			give << endl;
		}

	}
}

vector<int> Filters::combiSmplConvergeVec(const vector<string>& inNames){
	vector<int> retV(inNames.size(), -1);
	unordered_map<string, int> smpl2combi;
	unordered_map<string, int>::iterator s2cIT;
	int cntGrps(-1);
	for (size_t i = 0; i < SampleID_Combi.size(); i++){
		s2cIT = combiMapCollectGrp.find(SampleID_Combi[i]);
		if (s2cIT == combiMapCollectGrp.end()){
			cntGrps++;
			combiMapCollectGrp[SampleID_Combi[i]] = cntGrps;
		}
		smpl2combi[SampleID[i]] = combiMapCollectGrp[SampleID_Combi[i]];
	}
	for (size_t i = 0; i < inNames.size(); i++){
		s2cIT = smpl2combi.find(inNames[i]);
		if (s2cIT == smpl2combi.end()){
			cerr << "Can't find SampleID " << inNames[i] << " in reference Sample Names\n"; exit(113);
		}
		retV[i] = s2cIT->second;
	}
	return retV;
}


string Filters::shortStats( const string & file) {
	collectstats& cst = colStats[0];
	string ret("");
	if (file != ""){
		ret+= file + "\n";
	}
	if (pairedSeq > 1) {
		ret+= "Pair 1: ";
	}
	char buffer[50];
	float tmp = (100.f*float(cst.total - cst.totalRejected) / (float)cst.total);
	sprintf(buffer, "%.3f%% of %d", tmp, cst.total); ret += buffer;
	sprintf(buffer," reads accepted (%.3f%% end-trimmed)\n", (100.f* float(cst.total - cst.Trimmed) / (float)cst.total)); ret += buffer;
	
	if (pairedSeq > 1) {
		collectstats& cst = colStats[1];
		sprintf(buffer,"Pair 2: %.3f%% of %d", (100.f*float(cst.total - cst.totalRejected) / (float)cst.total), cst.total); ret += buffer;
		sprintf(buffer," reads accepted (%.3f%% end - trimmed)\n", (100.f* float(cst.total - cst.Trimmed) / (float)cst.total)); ret += buffer;
	}
	return ret;
}
void Filters::printGC(ostream& os,int Npair) {
	os << "Subset\t\tOccurence\t\tAvg.Quality\n";
	os << "\tA\tT\tG\tC\tA\tT\tG\tC\n";
	os << "R1 pre-filter";
	PreFiltP1->printGCstats(os);
	if ( Npair > 1 ) {
		os << "R2 pre-filter";
		PreFiltP2->printGCstats(os);
	}
	if ( !b_doFilter ) {return;}
	os << "R1 filtered";
	RepStat[0]->printGCstats(os);
	if ( Npair > 1 ) {
		os << "R2 filtered";
		RepStat[1]->printGCstats(os);
	}
}
void Filters::write2Demulti(shared_ptr<DNA> d, int p, int fqOvr) {
	if (!this->Demulti2Fls()) {
		return;
	}
	int idx = d->getBCnumber() - this->getBCoffset(); //correct for BC offset as well..
	if (idx < 0 || !d->isPassed()) {
		return;
	}
	d->prepareWrite(fqOvr);
	ofbufstream * tar = (demultiSinglFiles[idx][p]);
	if (tar == NULL) {
		ofstream tar2;
		tar2.open(demultiSinglFilesF[idx][p].c_str(), ios::app);
		d->writeFastQ(tar2,false);
		tar2.close();
	}
	else {
		d->writeFastQ(*(tar),false);
	}
}
void Filters::printStats(ostream& give, string file, string outf, bool main) {
	//TODO switch min_l to min_l_add
	collectstats& cst = colStats[0];
	collectstats& cst2 = colStats[1];
	if (cst.total != cst.total2) {
		cerr << "Unequal read numbers recorded " << cst.total << "," << cst.total2 << endl;
	}
	if (!main) {
		cst = statAddition;
	}
	bool p2stat = pairedSeq > 1 && main;
	give << "sdm " << sdm_version << " " << sdm_status << endl;
	if (file.length()>0){
		give<<"Input File:  "<<file<<endl;
	}
	if (outf == "-") {
		give << "Output File: stdout\n";
	}else if (outf.length() > 0) {
		give<<"Output File: "<<outf<<endl;
	}
	if ( !b_doFilter ) {
		give << "No valid Filter file provided; no filtering done on files\n";
		return;
	}
	if (!main) {
		give << "Statistics of reads that passed the mid Qual filter\n";
	} else {
		give << "Statistics of high quality reads\n";
	}
	float remSeqs = float (cst.total-cst.totalRejected);
	give << endl;
	if (!main){
		give << "Reads not High Qual: " << intwithcommas((int)cst.totalRejected);
	} else {
		give << "Reads processed: " << intwithcommas((int)cst.total);
		if (p2stat) {
			give << "; " << intwithcommas((int)cst2.total) << " (pair 1;pair 2)";
		}
	}
	give << endl;
	//int numAccept = (int)(cst.total - cst.totalRejected);
	int numAccept = (int)(cst.totalSuccess );
	if (!main){
		give << "Rejected:" << intwithcommas((int)(colStats[0].totalRejected - numAccept)) << endl << "Accepted: " << intwithcommas((int)numAccept) << " (" << intwithcommas((int)cst.Trimmed) << " were end-trimmed";
	} else {
		if (p2stat) {
			give << "Rejected: " << intwithcommas((int)cst.totalRejected) << "; " << intwithcommas((int)cst2.totalRejected) << endl << "Accepted: " << intwithcommas((int)numAccept) << "; " << intwithcommas((int)cst2.totalSuccess) << " (" << intwithcommas((int)cst.Trimmed) << "; " << intwithcommas((int)cst2.Trimmed) << " were end-trimmed";
		}
		else {
			give << "Rejected: " << intwithcommas((int)cst.totalRejected) << endl << "Accepted: " << intwithcommas((int)numAccept) << " (" << intwithcommas((int)cst.Trimmed) << " were end-trimmed";
		}
	}


	if ( false && bPrimerR ) { //confusing colStats
		give << ", with rev. primer: " << intwithcommas((int)cst.RevPrimFound); if ( p2stat ) { give << "; " << intwithcommas((int)cst2.RevPrimFound); }
	}
	give<<")"<<endl;

	if (pairedSeq>1) {
		give <<"Singletons among these: " << intwithcommas((int)cst.singleton) << "; " << intwithcommas((int)cst2.singleton) << endl;
	}
	give << "Bad Reads recovered with dereplication: " << intwithcommas((int)cst.DerepAddBadSeq) << endl;

	if ( bShortAmplicons ) {
		give << "Short amplicon mode.\n";
	}

	if ( checkBC2ndRd() ) {
		give << "Looked for switched read pairs (" << intwithcommas(revConstellationN) << " detected)" << endl;
	}
	if (main) {
		RepStat[0]->printStats2(give, remSeqs,0);
		RepStat[1]->printStats2(give, remSeqs,1);
	} else {
		RepStatAddition[0]->printStats2(give,remSeqs,0);
	}

	give << "Trimmed due to:\n";
	//EWwidth, EWthr  no stat for this so far
	float dval = (float)EWthr;	if (!main) { dval = (float)alt_EWthr; }
	int Xval = EWwidth; 
	if (EWthr > 0) {
		give << "  > " << EWthr << " avg qual in " << Xval << " bp windows : " << spaceX(10 - digitsFlt(dval)) << intwithcommas(cst.QWinTrimmed);
			if (p2stat) { give << "; " << intwithcommas((int)cst2.QWinTrimmed); } give << endl;
	}
	dval = (float)maxAccumQP;	if (!main) { dval = (float)alt_maxAccumQP; }
	if (maxAccumQP>0.0) {
		give << "  > (" << dval << ") acc. errors, trimmed seqs : " << spaceX(8 - digitsFlt(dval)) << intwithcommas((int)cst.AccErrTrimmed);
		if (p2stat) { give << "; " << intwithcommas((int)cst2.AccErrTrimmed); } give << endl;
	}

	give << "Rejected due to:\n";
	float val = (float)min_l;
	if (val == -1.f) {val = min_l_p;}
	if (!main){ val = (float)alt_min_l; }

	give << "  < min Seq length (" << val << ")  : " << spaceX(18 - digitsFlt(val)) << intwithcommas((int)cst.minL);
	if (p2stat) { give << "; " << intwithcommas((int)cst2.minL); } give << endl;
	if (cst.minLqualTrim>0){//this is failed because seq was too short after trimming
		give << "       -after Quality trimming : " << spaceX(10) << intwithcommas((int)cst.minLqualTrim);
		if (p2stat) { give << "; " << intwithcommas((int)cst2.minLqualTrim); } give << endl;
	}
	float valf = min_q;	if (!main){ valf = alt_min_q; }
	give << "  < avg Quality (" << valf << ")  : " << spaceX(21 - digitsInt((int)min_q)) << intwithcommas((int)cst.AvgQual);
	if (p2stat) { give << "; " << intwithcommas((int)cst2.AvgQual); } give << endl;
	give << "  < window (" << FQWwidth << " nt) avg. Quality (" << FQWthr << ")  : " << spaceX(5 - digitsInt(FQWwidth)) << intwithcommas((int)cst.QualWin);
	if (p2stat) { give << "; " << intwithcommas((int)cst2.QualWin); } give << endl;
	give << "  > max Seq length (" << max_l << ")  : " << spaceX(18 - digitsInt(max_l)) << intwithcommas((int)cst.maxL);
	if (p2stat) { give << "; " << intwithcommas((int)cst2.maxL); } give << endl;
	give << "  > (" << maxHomonucleotide << ") homo-nt run  : " << spaceX(21 - digitsInt(maxHomonucleotide)) << intwithcommas((int)cst.HomoNT);
	if (p2stat) { give << "; " << intwithcommas((int)cst2.HomoNT); } give << endl;
	int val2 = MaxAmb;	if (!main){ val2 = alt_MaxAmb; }
	give << "  > (" << val2 << ") amb. Bases  : " << spaceX(22 - digitsInt(val2)) << intwithcommas((int)cst.MaxAmb);
	if (p2stat) { give << "; " << intwithcommas((int)cst2.MaxAmb); } give << endl;
	if (BinFilP >= 0.f){
		give << "  > (" << BinFilErr << ") binomial est. errors : " << spaceX(13 - digitsFlt(BinFilErr)) << intwithcommas((int)cst.BinomialErr);
		if (p2stat) { give << "; " << intwithcommas((int)cst2.BinomialErr); } give << endl;
	}
	if ((bDoAdapter && tAdapter != "") || (bDoMultiplexing || cst.PrimerFail > 0) || ((!main && alt_bRequireFwdPrim) || bRequireFwdPrim)
		|| bPrimerR) {
		give << "Specific sequence searches:\n";
	}
	if (bDoAdapter && tAdapter!= "") {
		give << "  -removed adapter (" << tAdapter << ")  : " << spaceX(18 - (uint)tAdapter.length()) << intwithcommas((int)cst.adapterRem);
		if (p2stat) { give << "; " << intwithcommas((int)cst2.adapterRem); } give << endl;
	}
	if ( (bDoMultiplexing || cst.PrimerFail>0) || ((!main && alt_bRequireFwdPrim) || bRequireFwdPrim) ){
		give << "  -With fwd Primer remaining (<= " << PrimerErrs << " mismatches";
		if ((!main && alt_bRequireFwdPrim) || bRequireFwdPrim){
			give << ", required) : ";
			give << spaceX(1 - digitsInt(PrimerErrs));
		}
		else  {
			give <<") : "<< spaceX(11 - digitsInt(PrimerErrs));
		}
		give << intwithcommas((int)cst.PrimerFail);
		if ( p2stat ) { give << "; " << intwithcommas((int)cst2.PrimerFail) << endl; }
		else { give << endl; }
	} 
	if (bPrimerR){
		give<<"  -With rev Primer remaining (<= "<<PrimerErrs <<" mismatches";
		if ((!main && alt_bRequireRevPrim) || bRequireRevPrim){ give << ", required) : "; 
			give << spaceX(1 - digitsInt(PrimerErrs));
		}	else  {
			give << ") : "<<spaceX(11 - digitsInt(PrimerErrs));
		}
		give << intwithcommas((int)cst.PrimerRevFail);
		if ( p2stat ) { give << "; " << intwithcommas((int)cst2.PrimerRevFail) << endl; }
		else { give << endl; }
	}
	if (bDoMultiplexing){
		if (bDoBarcode){
			give << "  -Barcode unidentified (max " << TagErrs << " errors) : " << spaceX(19 - digitsInt(TagErrs)) << intwithcommas((int)cst.TagFail);
			if (p2stat && (cst2.TagFail > 0 || doubleBarcodes())) { give << "; " << intwithcommas((int)cst2.TagFail); give << " (" << intwithcommas((int)cst.dblTagFail) << " pairs failed)"; }
			give << endl;

			if (TagErrs>0){
				give << "    -corrected barcodes: " << spaceX(18) << intwithcommas((int)cst.suc_correct_BC);
				if (p2stat){give << "; " << intwithcommas((int)cst2.suc_correct_BC);	}
				give << endl;
				//<< ", failed to correct barcode: " << spaceX(5 - digitsInt(FQWwidth)) << intwithcommas((int)cst.fail_correct_BC) << endl;
			}
			
			if ( bDoBarcode2 ) {
				give << "    -used dual index barcodes";
				if ( BCdFWDREV[0].reversedBCs || BCdFWDREV[1].reversedBCs ) {
					give << " (reversed ";
					if ( BCdFWDREV[1].reversedBCs && BCdFWDREV[0].reversedBCs ) {
						give << " fwd & rev";
					} else	if ( BCdFWDREV[0].reversedBCs ) {
						give << " fwd";
					} else	if ( BCdFWDREV[1].reversedBCs ) {
						give << " rev";
					}
					give << " BCs)" << endl;
				}
				
			} else if ( BCdFWDREV[0].reversedBCs ) {
				give << "    -reversed all barcodes" << endl;
			}
			give << endl << "SampleID";
			if (bDoCombiSamples){
				give << "\tSampleGroup";
			}
			give << "\tBarcode";
			if ( bDoBarcode2 ) {give << "\tBarcode2";}
			give << "\tInstances\n";
			for (unsigned int i =0; i<Barcode.size();i++){
				give << SampleID[i] << "\t";
				if (bDoCombiSamples){ give << SampleID_Combi[i] << "\t"; }
				give << Barcode[i];
				if ( bDoBarcode2 ) {
					give << "\t"<<Barcode2[i];
				}
				give << "\t" << intwithcommas((int)cst.BarcodeDetected[i]) << endl;
			}
		} else if (bDoHeadSmplID){
			give << "  -Failed to assign sequences to header tag : " << intwithcommas((int)TagErrs )<< endl;
			give << endl << "SampleID\t";
			if (bDoCombiSamples){	give << "\tSampleGroup";	}
			give << "\tSampleID\tInstances\n"; 
			for (unsigned int i =0; i<Barcode.size();i++){
				give << SampleID[i] << "\t";
				if (bDoCombiSamples){ give << SampleID_Combi[i] << "\t"; }
				give << HeadSmplID[i] << "\t" << intwithcommas((int)cst.BarcodeDetected[i]) << endl;
			}
		}
	}
}

//statistics for each single sample 
void Filters::SmplSpecStats(ostream & give){
	collectstats& cst = colStats[0];
	collectstats& cst2 = colStats[1];
	bool p2stat = pairedSeq > 1 ;
	give << std::setprecision(3);
	give  << "SampleID";
	if ( bDoCombiSamples ) {
		give << "\tSampleGroup";
	}
	if ( bDoHeadSmplID ) {
		give << "\tSampleID";

	} else {
		if ( bDoBarcode2 ) { give << "\tBarcode\tBarcode2"; }
		else { give << "\tBarcode"; }
	}
	if ( p2stat ) {
		give << "\tRead1Accepted\tRead1Filtered\tRead1PassedFrac\tRead2Accepted\tRead2Filtered\tRead2PassedFrac\n";
	} else {
		give << "\tReadsAccepted\tReadsFailed\tPassed%\n";
	}
	for ( unsigned int i = 0; i<Barcode.size(); i++ ) {
		give << SampleID[i] << "\t";
		if ( bDoCombiSamples ) { give << SampleID_Combi[i] << "\t"; }
		if ( p2stat ) {
			give << Barcode[i] << "\t";
			if ( bDoBarcode2 ) { give << Barcode2[i] << "\t"; }
			give <<	cst.BarcodeDetected[i] << "\t" << cst.BarcodeDetectedFail[i] << "\t";
			float totSum = (float(cst.BarcodeDetected[i]) + float(cst.BarcodeDetectedFail[i]));
			if ( totSum > 0 ) {
				give << float(cst.BarcodeDetected[i]) / totSum << "\t";
			} else {give << "NA\t";	}

			give << cst2.BarcodeDetected[i] << "\t" << cst2.BarcodeDetectedFail[i] <<"\t" ;
			totSum = (float(cst2.BarcodeDetected[i]) + float(cst2.BarcodeDetectedFail[i]));
			if ( totSum > 0 ) {
				give << float(cst2.BarcodeDetected[i]) / totSum << "";
			} else { give << "NA"; }
			give<<endl;
		} else {
			give << Barcode[i] << "\t";
			if ( bDoBarcode2 ) { give << Barcode2[i] << "\t"; }
			give << (int)cst.BarcodeDetected[i]
				<< "\t" << cst.BarcodeDetectedFail[i] << "\t";
				
			float totSum = (float(cst.BarcodeDetected[i]) + float(cst.BarcodeDetectedFail[i]));
			if ( totSum > 0 ) {
				give << float(cst.BarcodeDetected[i]) / totSum << "";
			} else { give << "NA"; }
				give<< endl;
		}
	}



}

void ReportStats::calcSummaryStats(float remSeqs, unsigned int min_l, float min_q){
	if (remSeqs == 0){ return; }
	if (bMedianCalcs){
		rstat_Smed = (int) calc_median(rstat_VSmed,0.5f);
		rstat_Qmed = (int) calc_median(rstat_VQmed,0.5f);
		USQS=0.f;
	}
	RSQS = ( ( (float(rstat_NTs)/remSeqs)/(float)min_l ) +
			( (float(rstat_qualSum)/remSeqs) / min_q) ) / 2.f;
}

void Filters::addStats(shared_ptr<Filters> fil, vector<int>& idx){
	colStats[0].addStats(fil->colStats[0], idx); colStats[1].addStats(fil->colStats[1], idx);
	RepStat[0]->addStats(fil->RepStat[0]);	RepStat[1]->addStats(fil->RepStat[1]);
	if ( bAdditionalOutput){
		statAddition.addStats(fil->statAddition,idx);
		RepStatAddition[0]->addStats(fil->RepStatAddition[0]);
		RepStatAddition[1]->addStats(fil->RepStatAddition[1]);
	}
	PreFiltP1->addStats(fil->PreFiltP1);
	PreFiltP2->addStats(fil->PreFiltP2);
	maxReadsPerOFile = fil->maxReadsPerOFile;
	ReadsWritten = fil->ReadsWritten;//the idea here is to have a number of reads in CURRENT file, not total reads
	OFileIncre = fil->OFileIncre;
	revConstellationN += fil->revConstellationN;
}

void collectstats::addStats(collectstats& cs, vector<int>& idx){
	if (BarcodeDetected.size() > (uint)10000){
		cerr<<"Unrealistic number of barcodes (>10000) in addStats\n"; exit(79);}
	int BCS = (int)BarcodeDetected.size();
	for (unsigned int i=0;i<idx.size();i++){
		if (idx[i] >= BCS){ return; }
		//assert(idx[i] < BCS);
		BarcodeDetected[idx[i]] += cs.BarcodeDetected[i];
		BarcodeDetectedFail[idx[i]] += cs.BarcodeDetectedFail[i];
	}
	maxL += cs.maxL;	PrimerFail += cs.PrimerFail ;
	AvgQual += cs.AvgQual; HomoNT += cs.HomoNT;
	PrimerRevFail += cs.PrimerRevFail;
	minL += cs.minL ; minLqualTrim+= cs.minLqualTrim; TagFail += cs.TagFail;
	MaxAmb += cs.MaxAmb ; QualWin += cs.QualWin;
	Trimmed += cs.Trimmed ; AccErrTrimmed+= cs.AccErrTrimmed; total += cs.total;
	QWinTrimmed += cs.QWinTrimmed;
	totalRejected += cs.totalRejected;
	fail_correct_BC += cs.fail_correct_BC; suc_correct_BC += cs.suc_correct_BC ;
	failedDNAread += cs.failedDNAread; adapterRem += cs.adapterRem ;
	RevPrimFound += cs.RevPrimFound;
	singleton += cs.singleton;
	BinomialErr += cs.BinomialErr;
	dblTagFail += cs.dblTagFail;
	DerepAddBadSeq += cs.DerepAddBadSeq;
	total2 += cs.total2; totalSuccess += cs.totalSuccess;
}
void collectstats::reset(){
	singleton=0;
	size_t BCsiz = BarcodeDetected.size();
	for (unsigned int i=0; i<BCsiz;i++){
		BarcodeDetected[i]=0;
		BarcodeDetectedFail[i] = 0;
	}
	maxL=0; PrimerFail=0;AvgQual=0; HomoNT=0;
	PrimerRevFail=0;
	minL=0; TagFail=0; MaxAmb=0; QualWin=0;
	Trimmed=0; total=0; totalRejected=0;
	fail_correct_BC=0; suc_correct_BC=0;failedDNAread=0;
	adapterRem = 0; RevPrimFound = 0; DerepAddBadSeq = 0;
	total2 = 0; totalSuccess = 0;
}

inline void ReportStats::addDNAStats(shared_ptr<DNA> d){
	//pretty fast
	addMeanStats(d->length(),(int) d->getAvgQual(), (float)d->getAccumError());
	//NT specific quality scores
	d->NTspecQualScores(QperNT, NTcounts);
	//more memory intensive
	if (bMedianCalcs){
		//quali
		uint avq = (uint) (d->getAvgQual()+0.5f);
		//if (avq >= 10000){cerr << d->getID() << "high q: "<< d->getAvgQual()<<endl;}
		add_median2histo(avq,rstat_VQmed);
		//read length
		//if ((uint) d->length() >= 10000){cerr << d->getID() << "high l: "<< d->length()<<endl;}
		add_median2histo(d->length() ,rstat_VSmed);
	}
}
ReportStats::ReportStats(bool MedianDo):
	bMedianCalcs(MedianDo),rstat_totReads(0),rstat_NTs(0),rstat_qualSum(0),
	rstat_Qmed(0),rstat_Smed(0),RSQS(0.f),USQS(0.f),rstat_accumError(0.f),
	QperNT(1000,0), NTcounts(1000,0),
	rstat_VQmed(0), rstat_VSmed(0)
{}
void ReportStats::reset() {
	rstat_totReads = 0; rstat_NTs = 0; rstat_qualSum=0;
	rstat_Qmed = 0; rstat_Smed = 0; 
	RSQS = 0.f; USQS = 0.f; rstat_accumError = 0.f;
	QperNT.resize(1000,0); NTcounts.resize(1000,0);
	std::fill(QperNT.begin(), QperNT.end(), 0);
	std::fill(NTcounts.begin(), NTcounts.end(), 0);
	rstat_VQmed.resize(0); rstat_VSmed.resize(0);
}
unsigned int ReportStats::lowest(const vector<uint>& in){
	for (uint i=0;i<(uint)in.size();i++){
		if (in[i]>0){return i;}
	}
	return(0);
}
unsigned int ReportStats::highest(const vector<uint>& in){
	if (in.size()==0){return 0;}
	for (uint i=(uint)in.size()-1;i>=0;i--){
		if (in[i]>0){return i;}
	}
	return 0;
}
void ReportStats::printGCstats(ostream& give) {
	//NT_POS['A'] = 0; NT_POS['T'] = 1; NT_POS['G'] = 2; NT_POS['C'] = 3;	NT_POS['N'] = 4;
	vector<string> NTs(6, "X"); NTs[0] = "A"; NTs[1] = "T";
	NTs[2] = "G"; NTs[3] = "C"; NTs[4] = "N";
	for ( uint i = 0; i < 4; i++ ) {
		give << "\t" << NTcounts[i];
	}
	//give << endl;
	for ( uint i = 0; i < 4; i++ ) {
		give << "\t" << float(QperNT[i]) / float(NTcounts[i]);
	}
	give << endl;
}
void ReportStats::printStats2(ostream& give, float remSeqs,int pair){
	if ( pair == 1 ) {
		return;//deactivate for now
	}
	if ( pair == 0 ) {
		if ( bMedianCalcs ) {
			unsigned int minS = lowest(rstat_VSmed);
			unsigned int maxS = highest(rstat_VSmed);
			unsigned int minQ = lowest(rstat_VQmed);
			unsigned int maxQ = highest(rstat_VQmed);
			give << "Min/Avg/Max stats Pair 1";// -RSQS : "<<RSQS;
			if ( remSeqs == 0 ) {
				give << "\n     - Seq Length : " << "0/0/0"
					<< "\n     - Quality :   " << "0/0/0";
			} else {
				give << "\n     - Seq Length : " << minS << "/" << float(rstat_NTs) / (float)rstat_totReads << "/" << maxS
					<< "\n     - Quality :   " << minQ << "/" << float(rstat_qualSum) / (float)rstat_totReads << "/" << maxQ;
			}
		} else {
			give << "Average Stats - RSQS : " << RSQS;
			if ( remSeqs == 0 ) {
				give << "\n     - Seq Length : " << "0/0/0"
					<< "\n     - Quality :   " << "0/0/0";
			} else {
				give << "\n     - Seq Length : " << float(rstat_NTs) / (float)rstat_totReads
					<< "\n     - Quality :   " << float(rstat_qualSum) / (float)rstat_totReads;
			}
		}
	} else {
		give << "Pair 2 stats";
	}
	
		if (bMedianCalcs){
			//give << "Median Stats Pair 1";// -USQS : " << USQS;
			give << "\n     - Median Seq Length : " << rstat_Smed << ", Quality : " << rstat_Qmed;// << "\n";
		}
		give << "\n     - Accum. Error " << (rstat_accumError / (float)rstat_totReads) << "\n";
}

//add the stats from a different ReportStats object
void ReportStats::addStats(shared_ptr<ReportStats> fil){
	//report stats:
	rstat_NTs += fil->rstat_NTs; rstat_totReads += fil->rstat_totReads;
	rstat_qualSum += fil->rstat_qualSum;
	rstat_accumError += fil->rstat_accumError;
	for ( uint i = 0; i < 6; i++ ) {
		QperNT[i] += fil->QperNT[i];
		NTcounts[i] += fil->NTcounts[i];
	}
	if (bMedianCalcs){
		//vectors for median calcs
		if (rstat_VQmed.size() < fil->rstat_VQmed.size()){
			rstat_VQmed.resize(fil->rstat_VQmed.size(),0);
			assert(rstat_VQmed.size() < 10000);
		}
		for (unsigned int i=0; i<fil->rstat_VQmed.size(); i++){
			rstat_VQmed[i] += fil->rstat_VQmed[i];
		}
		if (rstat_VSmed.size() < fil->rstat_VSmed.size()){
			rstat_VSmed.resize(fil->rstat_VSmed.size(),0);
			//times change.. wrong assert here
			//assert(rstat_VSmed.size() < 10000);
		}
		for (unsigned int i=0; i<fil->rstat_VSmed.size(); i++){
			rstat_VSmed[i] += fil->rstat_VSmed[i];
		}
	}
}
//calculate median value from data stored as histogram-vector
// for median use perc = 0.5f
float ReportStats::calc_median(vector<unsigned int>& in, float perc){
	unsigned int sum = 0;
	for (unsigned int i=0; i<in.size(); i++){
		sum += in[i];
	}
	unsigned int tar = (unsigned int) (((float)sum) * perc);
	sum = 0;
	for (unsigned int i=0; i<in.size(); i++){
		sum += in[i];
		if (sum >= tar){
			return (float) i;
		}
	}

	return 0.f;
}
void ReportStats::add_median2histo(vector<unsigned int>& in, vector<unsigned int>& histo)
{
	unsigned int max = *max_element(in.begin(),in.end());
	if (max> histo.size()){
		if (max > 10000){cerr<<"max bigger 10000.\n"; exit(77);}
		histo.resize(max,0);
	}
	for (unsigned int i=0; i<histo.size(); i++){
		histo[ in[i] ] ++;
	}
}
vector<size_t> ReportStats::getVrange(int which) {
	if (which == 1) {
		return medVrange(rstat_VQmed);
	} else {
		return medVrange(rstat_VSmed);
	}
}
vector<size_t> ReportStats::medVrange(const vector<uint> x) {
	vector<size_t> ret(2, 0); ret[1] = 0; bool empty = true;
	for (size_t i = 0; i < x.size(); i++) {
		if (x[i]>0) {
			if (empty) { ret[0] = i; empty = false; }
			ret[1] = i;
		}
	}
	ret[1]++;
	return ret;
}

void ReportStats::add_median2histo(unsigned int in, vector<unsigned int>& histo)
{
	if (in >= histo.size()){
		//int oldS=histo.size();
		histo.resize(in+3,0);
		assert(in < 1e6);
		//explicit ini
		/*for (int i=oldS; i<histo.size(); i++){
			histo[i]=0;
		}*/
	}
	histo[ in ] +=1;
}

////////////  UCLINKS  ///////////////
UClinks::~UClinks(){
	ucf.close();
	mapdere.close();
/*
	for (uint i=0; i< bestDNA.size();i++){
		delete bestDNA[i];
	}
//	for (std::map<string, shared_ptr<DNA>>::iterator iterator = unusedID.begin(); iterator != unusedID.end(); iterator++) {
//		delete (*iterator).second;
//	}
	for (uint i = 0; i<oldDNA.size(); i++) {
		if (oldDNA[i] != NULL) { delete oldDNA[i]; }
	}
	for (uint i = 0; i<oldDNA2.size(); i++) {
		if (oldDNA2[i] != NULL) { delete oldDNA2[i]; }
	}
	*/
}
UClinks::UClinks( OptContainer& cmdArgs):
	CurSetPair(-1),maxOldDNAvec(20000),
	oldDNAid(maxOldDNAvec,""),oldDNA(maxOldDNAvec,NULL),
	oldDNA2(maxOldDNAvec, NULL), 
	DNAunusedPos(0), derepMapFile(""),
	bestDNA(0, NULL), oriKey(0), bestPID(0), bestLEN(0),
	clusCnt(0), uclines(0),
	SEP(""), 
	UCread(false), pairsMerge(false), MAPread(false),
	b_derepAvailable(false),
	UPARSE8up(false), UPARSE9up(false), UPARSE11up(false),
	UpUcFnd(false),
	otuTerm("OTU"),
	RefDBmode(false), RefDBotuStart(-1),
	SeedsAreWritten(false),
	OTUmat(0), unregistered_samples(false),
	doChimeraCnt(false), OTUnumFixed(true){
	//read in UC file and assign clusters
	SEP = cmdArgs["-sample_sep"];
	string str = cmdArgs["-optimalRead2Cluster"];
	

	ucf.open(str.c_str(),ios::in);
	if (!ucf){
		UCread = true;
		cerr<<"Could not open uc file\n"<< str<<endl; exit(46);
	}
	if (cmdArgs.find("-derep_map") != cmdArgs.end()) {
		derepMapFile = cmdArgs["-derep_map"];
	}
	if (cmdArgs.find("-count_chimeras") != cmdArgs.end() &&
		cmdArgs["-count_chimeras"] == "T") {
		doChimeraCnt = true;
	}
	if (cmdArgs["-uparseVer"] != "") {
		int upVer = atoi(cmdArgs["-uparseVer"].c_str());
		if (upVer >= 8 && upVer < 9) {
			UPARSE8up = true; UpUcFnd = true;
		}	else if (upVer >= 9 && upVer <11) {
			UPARSE8up = true; UPARSE9up = true; UpUcFnd = true;
		}	else if (upVer >= 11) {
			UPARSE8up = true; UPARSE9up = true; UPARSE11up = true; UpUcFnd = true;
			otuTerm = "otu";
		}
	}
}


void UClinks::readDerepInfo(string dereM) {
	int cnt(0);
	mapdere.open(dereM.c_str(), ios::in);
	if (!mapdere) {
		cerr << "Can't open " << dereM << ". \nAborting\n"; exit(54);
	}
	b_derepAvailable = true;
	string line("");
	//only read header
	while (!safeGetline(mapdere, line).eof()) {
		//	while(getline(in,line,'\n')) {
		if (line.length() < 3) { continue; }
		string segments;
		stringstream ss;
		ss << line;
		int tbcnt = -1;

		//cmdArgs["-i_MID_fastq"]
		while (getline(ss, segments, '\t')) {
			tbcnt++;
			if (cnt == 0) { //search for header
				if (tbcnt == 0) {
					if (segments != "#SMPLS") { cerr << "First line in dereplicate map has to start with #SMPLS, corrupted file.\n" << dereM << endl; exit(98); }
					continue;
				}
				vector<string> spl = header_string_split(segments, ":");
				SmplIDs[spl[1]] = stoi(spl[0]);				
			} else {//actual derep info per sorted line
				cerr << "wrong mapping reading\n"; exit(64);
			}
		}
		if (cnt == 0) {
			OTUmat.resize(SmplIDs.size(), vector<matrixUnit>(clusCnt + 1, 0));
		}
		break;
	}
	
}

//read in dereplicated info feom derp.map
void UClinks::oneDerepLine(shared_ptr<DNAunique> d) {
	if (MAPread){
		return;
	}
	string line("");
	while (!safeGetline(mapdere, line).eof()) {
		//if (line.length() < 3) { continue; }
		string segments;
		//stringstream ss;
		//string head("");
		//ss << line;
		int tbcnt = -1;
		int curCnt(0),curCnt2(0);
		size_t strpos(line.find_first_of('\t')), lastpos(0);
		while (strpos != string::npos){
		//while (getline(ss, segments, '\t')) {
			segments = line.substr(lastpos, strpos - lastpos);
			lastpos = strpos+1;
			strpos = line.find_first_of('\t', lastpos);
			tbcnt++;
			if (tbcnt == 0) {
				if (!d->sameHead(segments)) {
					cerr << segments << " is not " << d->getID() << endl;
				}
				size_t idx = segments.find(";size=");
				if (idx != string::npos) {
					curCnt = atoi(segments.substr(idx + 6, segments.find(";", idx + 5) - (idx + 6)).c_str());
					//curCnt
				}
				continue;
			}
			size_t spl = segments.find(":"); int occur = stoi(segments.substr(spl+1));
			d->setOccurence(stoi(segments.substr(0,spl)), occur);
			curCnt2 += occur;
			//string tmp = segments.substr(0, spl);
		}
		//last round
		segments = line.substr(lastpos);
		size_t spl = segments.find(":"); int occur = stoi(segments.substr(spl + 1));
		d->setOccurence(stoi(segments.substr(0, spl)), occur);
		curCnt2 += occur;

		if (curCnt2 != curCnt) {
			cerr << "ERROR: Mapping file unique abundance reconstruction: failed to find correct count (" << curCnt2 << " vs " << curCnt<<"):\n" << line << endl;
			exit(67);
		}
		break;
	}

	if (mapdere.eof()) {// end of file
		mapdere.close();
		MAPread = true;
	}
}

//finishes reading the map, loading the sequences as DNAunique into mem (no DNA, mind)
void UClinks::finishMAPfile(){
	if (MAPread){
		return;
	}
	string line("");
	bool hardAdd = false;

	//go through ucl file line by line
	while (!safeGetline(mapdere, line).eof()) {	
		string segments;
		int tbcnt = -1;
		int curCnt(0), curCnt2(0);
		size_t strpos(line.find_first_of('\t')), lastpos(strpos + 1);
		segments = line.substr(0, strpos);
		//1st setup empty DNA for each remaining line
		shared_ptr<DNAunique> d = make_shared<DNAunique>("", segments);
		d->seal();
		size_t idx = segments.find(";size=");
		if (idx != string::npos) {
			curCnt = atoi(segments.substr(idx + 6, segments.find(";", idx + 5) - (idx + 6)).c_str());
		}
		strpos = line.find_first_of('\t', lastpos);

		while (strpos != string::npos){
			segments = line.substr(lastpos, strpos - lastpos);
			lastpos = strpos + 1;
			strpos = line.find_first_of('\t', lastpos);
			tbcnt++;
			size_t spl = segments.find(":"); int occur = stoi(segments.substr(spl + 1));
			d->setOccurence(stoi(segments.substr(0, spl)), occur);
			curCnt2 += occur;
			//string tmp = segments.substr(0, spl);
		}
		//last round
		segments = line.substr(lastpos);
		size_t spl = segments.find(":"); int occur = stoi(segments.substr(spl + 1));
		d->setOccurence(stoi(segments.substr(0, spl)), occur);
		curCnt2 += occur;

		if (curCnt2 != curCnt) {
			cerr << "ERROR: Mapping file unique abundance reconstruction: failed to find correct count (" << curCnt2 << " vs " << curCnt << "):\n" << line << endl;
			exit(67);
		}

		//add d to oldDNA
		string curID = d->getID();
		curID = curID.substr(0, curID.find_first_of(' '));
		unusedID[curID] = DNAunusedPos;
		if (hardAdd){
			oldDNA.push_back(d); oldDNAid.push_back(curID);
			DNAunusedPos++;
		} else {
			oldDNA[DNAunusedPos] = d;
			oldDNAid[DNAunusedPos] = curID;
			DNAunusedPos++;
			if (DNAunusedPos >= maxOldDNAvec){
				hardAdd = true;
			}
		}
	}

	if (mapdere.eof()) {// end of file
		mapdere.close();
		MAPread = true;
	}
	oldDNA2.resize(oldDNA.size(), NULL);

}

//functions selects optimal read to represent OTU
void UClinks::findSeq2UCinstruction(shared_ptr<InputStreamer> IS, bool readFQ, 
	shared_ptr<Filters> fil){
	if (UCread){return;}

	if (IS->hasMIDseqs()){
		CurSetPair = 0;
	} else {
		CurSetPair = -1;
	}
#ifdef DEBUG
	cerr << "UC seed extension" << endl;
#endif

	DNAidmapsIT fndDNAinOld;
	std::list<string>::iterator listIT;
	shared_ptr<DNAunique> match(NULL); shared_ptr<DNA>  match2(NULL);
	bool cont(true),cont2(true);
	string segs("");
	string segs2;
	float perID;
	vector<int> curCLID(0,0);
	bool sync(false); // syncing of 2 read pairs; not implemented for this function yet

	while ( getUCFlineInfo(segs, segs2, perID, curCLID, !b_derepAvailable) ) {
		//int subcnt = 0;
		if ( curCLID.size() == 0 ) { continue; }
		if ( uclInOldDNA(segs, curCLID, perID, fil) ) {
			curCLID.resize(0);
			continue;
		}

		while(cont){

			shared_ptr<DNA> tmpDNA = IS->getDNA(cont, 0, sync);
			if (tmpDNA == NULL) { break; }//signal that at end of file
			match.reset( new DNAunique(tmpDNA, -1));
			//delete tmpDNA;
			oneDerepLine(match);
			match2 = IS->getDNA(cont2, 1, sync);
			string curID = match->getID();
			curID = curID.substr(0,curID.find_first_of(' '));
			//check if tdn1 a) matches ID b) is better
			if (curID != segs){
				//block to store unused DNA & find this id in this block
				unusedID[curID] = DNAunusedPos;
				if (oldDNAid[DNAunusedPos]!=""){//delete old DNA at position
					fndDNAinOld = unusedID.find(oldDNAid[DNAunusedPos]); unusedID.erase(fndDNAinOld);
					//delete oldDNA[DNAunusedPos];	delete oldDNA2[DNAunusedPos];
				}
				oldDNA[DNAunusedPos] = match;	oldDNA2[DNAunusedPos] = match2;
				oldDNAid[DNAunusedPos] = curID;
				DNAunusedPos++;
				if (DNAunusedPos>= maxOldDNAvec){	DNAunusedPos=0;		}
			} else {
				//assign % identity score to DNA object
#ifdef DEBUG
				cerr << "UC Hit";
#endif
				match->setTempFloat(perID);
				besterDNA(curCLID, match, match2, fil);
				curCLID.resize(0);
				break;
			}
			
			//if (tdn!=NULL && ch1 != tdn->isPassed()){cerr<<"isPassed is != ch1! Aborting..\n";exit(12);}
		}
	}
#ifdef DEBUG
	cerr << "UC seed initialized" << endl;
#endif

}

void UClinks::finishUCfile(shared_ptr<Filters> fil, string addUC, bool bSmplHd){
	if (UCread){
		addUCdo(addUC,bSmplHd);
		UCread = true;
		return;
	}
	string segs;
	string segs2;
	float perID(0.f);
	vector<int> curCLID(0);

	while ( getUCFlineInfo(segs, segs2, perID, curCLID, !b_derepAvailable) ) {
		if ( uclInOldDNA(segs, curCLID, perID, fil) ) {
			curCLID.resize(0);
			continue;
		}
		cerr << segs << " ";
	}
	
	addUCdo(addUC, bSmplHd);
	UCread = true;
}
void UClinks::addUCdo(string addUC,bool SmplHd) {
	if (addUC == "") {
		return;
	}
	string segs;
	string segs2;
	float perID;
	vector<int> curCLID(0);
	ucf.open(addUC.c_str(), ios::in);
	if (!ucf) {
		UCread = true;
		cerr << "Could not find additional uc file: " << addUC << endl;
	} else {
		std::cerr << "Reading " << addUC << endl;
	}
	UCread = false;
	if (SmplHd){
		while (getUCFlineInfo(segs, segs2, perID, curCLID, SmplHd)) { curCLID.resize(0); }
	}else {//complicated..
		int cnt(0);
		while (getUCFlineInfo(segs, segs2, perID, curCLID, SmplHd)) { 
			//find segs in remaining dereps
			cnt++; //if (cnt < 61403) { continue; }
			if (curCLID.size() == 0) { continue; }
			if (uclInOldDNA_simple(segs, curCLID)) {
				curCLID.resize(0);
				continue;
			}

			curCLID.resize(0); 
		}
	}

}

//used for specific situation where only empty DNA string was read with derep info attached
bool UClinks::uclInOldDNA_simple(const string& segs,const vector<int>& curCLID) {
	//check if DNA is in group of unmatched DNA's ?
	if (segs == ""){ return false; }//empty ucl line
	DNAidmapsIT unusedIT = unusedID.find(segs);
	if (unusedIT != unusedID.end()){// found something
		int mID = (*unusedIT).second;
		if (mID >= (int)oldDNA.size()) { cerr << "SEC MID too high\n"; exit(55); }
		if (oldDNA[mID] == NULL){ return false; }
		//give sequence a chance to be selected
		matrixUnit matchSiz = (matrixUnit)curCLID.size();
		for (int k = 0; k < matchSiz; k++) {
			add2OTUmat(oldDNA[mID], curCLID[k], matchSiz);
		}
		//was matched once to an OTU seed. Even if several later matches, doesn't mater - delete

		unusedID.erase(unusedIT);
//		delete oldDNA[mID]; if (oldDNA2[mID] != NULL){ delete oldDNA2[mID]; }

		oldDNA[mID] = NULL;
		oldDNA2[mID] = NULL;
		oldDNAid[mID] = "";
		return true;
	}
	return false;
}
bool UClinks::uclInOldDNA(const string& segs,const vector<int>& curCLID, float perID,
	shared_ptr<Filters> fil) {
	//check if DNA is in group of unmatched DNA's ?
	if (segs == ""){ return false; }//empty ucl line
	DNAidmapsIT unusedIT = unusedID.find(segs);
	if (unusedIT != unusedID.end()){// found something
		//shared_ptr<DNA> fndDNA = (*unusedIT).second;
		//fndDNA->setTempFloat(perID);
		//besterDNA(curCLID, fndDNA, fil);
		//unusedID.erase(unusedIT);

		int mID = (*unusedIT).second;
		//give sequence a chance to be selected
		oldDNA[mID]->setTempFloat(perID);
		besterDNA(curCLID, oldDNA[mID], oldDNA2[mID], fil);
		//remove all trace
		unusedID.erase(unusedIT);
		oldDNA[mID] = NULL;
		oldDNA2[mID] = NULL;
		oldDNAid[mID] = "";
		return true;
	}
	return false;
}

bool UClinks::getUCFlineInfo(string& segs, string& segs2,float& perID, 
	vector<int>& curCLID,  bool addFromHDstring) {
	//reads UC file line by line
	//can also be used to delineate UC's


	if (UCread){return false;}
	//close all file streams
	if (ucf.eof()){
		UCread=true;
		ucf.close();
		return false;
	}
	string line; 
	std::unordered_map<string, int>::iterator itCL;
	while (getline(ucf, line, '\n')) {
	//	cerr<<line<<endl;
		uclines++;
		if (line.length() <= 1){ continue; }
		
		if (!UpUcFnd){
			if (line.substr(1, 1) == "\t"){ UPARSE8up = false; 
			}else {	UPARSE8up = true;	}

			if (UPARSE8up){//could also be uparse9

				stringstream ss2;
				ss2 << line;   int tabCnt = 1;
				while (getline(ss2, segs, '\t')) {
					tabCnt++;
				}
				if (tabCnt >= 3) {
					UPARSE9up = true;
					cerr << "Switching to Uparse 10+ style map file.\n";
				}
			}


			UpUcFnd = true;
		}
		
		stringstream ss;
		ss << line;
		bool chimera = false;
		vector<string>tarsV; //saves hits to OTUs
		//2 ways to get to a) hit info b) query & otu
		if (!UPARSE8up){
			if ( (line.substr(0, 1) != "H")) {
				continue;
			}

			for (uint i = 0; i < 4; i++){//jump to pos X
				getline(ss, segs, '\t');
			}
			perID = (float)atof(segs.c_str());
			for (uint i = 0; i < 5; i++){//jump to pos X
				getline(ss, segs, '\t');
			}
			getline(ss, segs2, '\t');
			
			

		} else if (!UPARSE9up){ // 
			//query first entry
			string tmp;
			getline(ss, segs, '\t');//0
			getline(ss, tmp, '\t');//1
			//should be "match"
			if ( tmp == "chimera") {
				if ( !doChimeraCnt ) {continue;}
				chimera = true; }// segs = ""; continue;}
			if (tmp == "otu"){
				segs2 = segs;
				perID = 100.f;
			} else {
				getline(ss, tmp, '\t');//2
				perID = (float)atof(tmp.c_str());
				//indicator if hit
				//OTU last entry
				getline(ss, segs2, '\t');//3
				getline(ss, segs2, '\t');//4
			}
		} else { //UP9, uparse 10
				 //query first entry
			string tmp;
			getline(ss, segs, '\t');//0
			getline(ss, tmp, '\t');//1
								   //should be "match"
			if (tmp == "chimera" ) {
				continue;
			} else if ( tmp == "noisy_chimera" ||  tmp == "good_chimera") { //tmp == "perfect_chimera" ||
				if (!doChimeraCnt) { continue; }
				chimera = true;
			}else if (tmp == "perfect_chimera") {
				removeSizeStr(segs);
				perfectChims.insert(segs);
				continue;
			}// segs = ""; continue;}
			if (tmp.substr(0,3) == otuTerm ){
				segs2 = segs;
				perID = 100.f;
			} 
			else {//match or perfect_chimera case
				getline(ss, tmp, '\t');//2
				//dqt=1;top=GZV0ATA01ANJXZ;size=14;(99.6%);
				size_t p1(tmp.find("top=")+4);//4
				size_t p2(tmp.find(";(", p1)+2);
				size_t p3(tmp.find("%);", p2));
				segs2 = tmp.substr(p1, p2 - p1-1);
				//string xx = tmp.substr(p2, p3 - p2);
				perID = (float)atof(tmp.substr(p2,p3-p2).c_str());
				if (false && chimera) {//just use up9 top hit
//					p1(tmp.find(";top=") + 5,p2);
//					p2(tmp.find(";(", p1) + 2);

				}
			}
		}
		//remove spaces
		segs = segs.substr(0,segs.find_first_of(' '));
		//also remove sample identifier in string
		string smplID = "";
		removeSampleID(segs, SEP, smplID);

		if ( chimera && UPARSE8up) {
			tarsV = splitByComma(segs2, false, '+');
			for ( uint kk = 0; kk < tarsV.size(); kk++ ) {
				tarsV[kk] = tarsV[kk].substr(0,tarsV[kk].find_last_of("("));
			}
		} else if (!chimera){
			tarsV.push_back(segs2);
		}
		matrixUnit splChim = (matrixUnit)tarsV.size();
		curCLID.resize(tarsV.size());

		for ( uint kk = 0; kk < tarsV.size(); kk++ ) {
			//goes over all chimeric hits
			string oriClKey (tarsV[kk]);
			segs2 = oriClKey;
			//remove ;size= argument
			removeSizeStr(segs2);
			removeSampleID(segs2, SEP);
			//identifiers are ready, no identify first which cluster this is
			//curCLID = -1;

			//find the cluster in the list of registered clusters
			itCL = seq2CI.find(segs2);

			if ( itCL == seq2CI.end() ) {//new OTU?
				if (OTUnumFixed) {
					//cerr << "XX\n";
					if (perfectChims.find(segs2) == perfectChims.end()) {
						cerr << "Unkown OTU entry found:" << segs2 << endl;
					}
				}	else {
					//not found in known clusters.. create entry
					bestDNA.push_back(NULL);
					bestDNA2.push_back(NULL);
					oriKey.push_back(oriClKey);
					bestPID.push_back(0.f);
					bestLEN.push_back(0);
					clusCnt = (int)bestDNA.size() - 1;
					curCLID[kk] = clusCnt;
					seq2CI[segs2] = clusCnt;
				}
			} else {//cluster exists
				curCLID[kk] = ((*itCL).second);
			}
#ifdef matrix_sum
			if ( addFromHDstring ) {
				//cerr << segs << "\t" << segs2 << endl;
				add2OTUmat(smplID, (*itCL).second, splChim);
			} 
#endif
		}
		return true;
	}
	return true;
}
void UClinks::writeOTUmatrix(string outf,shared_ptr<Filters> fil) {
	cerr << "Writing OTU matrix to " << outf << endl;
	this->setOTUnms();
	ofstream MA;
	MA.open(outf);
	//first write all sample IDs
	std::unordered_map<string, int>::iterator OTUid;
	MA << "OTU";
	matrixUnit totalCnts = 0;
	/*const vector<string>& smpls_t = fil->SampleID;
	
	if (!unregistered_samples) {
		smpls = smpls_t;
	} else {*/
	vector<string> smpls(SmplIDs.size(),"");
	unordered_map<string, int>::iterator smplNit;
	for (smplNit = SmplIDs.begin(); smplNit != SmplIDs.end(); smplNit++) {
		smpls[smplNit->second] = smplNit->first;
	}
	//}
	for (size_t i = 0; i < smpls.size(); i++) {
		MA<<"\t"<<smpls[i];
	}
	//write OTU values
	for (OTUid = seq2CI.begin(); OTUid != seq2CI.end(); OTUid++) {
		MA << endl << oriKey[OTUid->second];
		int rowI = OTUid->second;
		for (size_t i = 0; i < smpls.size(); i++) {
//		for (smplNit = SmplIDs.begin(); smplNit != SmplIDs.end(); smplNit++) {
			matrixUnit val = OTUmat[ SmplIDs[smpls[i]] ][rowI];
			totalCnts += val;
			//char num[24];sprintf(num, "\t%.2f", val);
			//MA << num;
			if ( isnan(val) ) {
				MA << "\t0";
			} else {
				MA << "\t"<<round(val);
			}
		}
	}
	MA.close();
	cerr << "Recruited " << (int)totalCnts << " reads in OTU matrix\n";
}
void UClinks::setOTUnms(){
	int cnt = 0;
	std::unordered_map<string, int>::iterator OTUid;
	for (OTUid = seq2CI.begin(); OTUid != seq2CI.end(); OTUid++) {
		string newlySetID = "OTU_" + itos(cnt);
		/*if (newlySetID == "OTU_7711"){
			int x = 0;	cerr << "7711 ID " << OTUid->second; if (bestDNA[OTUid->second] == NULL){ cerr << "no DNA\n"; }
			else { cerr << "YES\n"; }
		}*/
		oriKey[OTUid->second] = newlySetID;
		cnt++;
	}

}
void UClinks::add2OTUmat(const string& smplID, int curCLID, matrixUnit spl) {

	std::unordered_map<string, int>::iterator smplNit = SmplIDs.find(smplID);
	if (smplNit == SmplIDs.end()) {//create entry & expand matrix
		SmplIDs[smplID] = (int) OTUmat.size();
		OTUmat.push_back(vector<matrixUnit>(clusCnt + 1, (matrixUnit)0));
		smplNit = SmplIDs.find(smplID);
		cerr << "New Sample ID in uc file detected, that is not present in map: " << smplID<< endl;
		unregistered_samples = true;
	}
	//easy, now add in
	if ( spl <1 ) { spl = 1; }
	OTUmat[(*smplNit).second][curCLID] += (matrixUnit)1 / spl;

}
void UClinks::add2OTUmat(shared_ptr<DNAunique> d, int curCLID, matrixUnit rep) {
	if (d == NULL) { cerr << " add2OTUmat::d is NULL\n"; exit(85); }
	unordered_map <int, int> map = d->getDerepMap();
	if ( rep <1 ) { rep = 1; }
	for (auto iterator = map.begin(); iterator != map.end(); iterator++) {
		OTUmat[iterator->first][curCLID] += (matrixUnit)iterator->second / rep;
	}
}

void UClinks::setupDefSeeds(shared_ptr<InputStreamer> FA, shared_ptr<Filters> fil) {
	bool contRead = true; bool sync(false);
	while (contRead) {
		shared_ptr<DNA> tmpDNA = FA->getDNA(contRead, 0,sync);
		if (tmpDNA == NULL) { break; }
		shared_ptr<DNAunique> tmp = make_shared<DNAunique>(tmpDNA, -1);
//		delete tmpDNA;
		//second pair
		//shared_ptr<DNA> tmp2 = FA->getDNA(contRead, 1);
		string segs2 = tmp->getID();
		string oriClKey = segs2;
		//remove ;size= argument
		size_t idx = segs2.find(";size=");
		segs2 = segs2.substr(0, idx);
		//also remove sample identifier in string
		removeSampleID(segs2, SEP);

		//cluster should not exist, test
		std::unordered_map<string, int>::iterator itCL;
		itCL = seq2CI.find(segs2);

		if (itCL == seq2CI.end()) {
			//not found in known clusters.. create entry
			bestDNA.push_back(tmp);
			bestDNA2.push_back(NULL);
			oriKey.push_back(oriClKey);
			bestPID.push_back(100.f);
			bestLEN.push_back(tmp->length());
			clusCnt = (int)bestDNA.size() - 1;
			seq2CI[segs2] = clusCnt;
		} else {//cluster exists
			cerr << "Found double ID " << segs2 << endl << "Aborting.." << endl;
			exit(74);
		}
	}
	//std::map<string, int>::iterator smplNit;
	if (derepMapFile != "") {
		readDerepInfo(derepMapFile);
	} else {
		//TODO: fix this smplnew
		const vector<string>& smpls = fil->SampleID;
		for (uint i = 0; i < smpls.size(); i++) {
			SmplIDs[smpls[i]] = (int)OTUmat.size();
			if (i != OTUmat.size()) { cerr << "Err in setupDefSeeds\n"; exit(453); }
			OTUmat.push_back(vector<matrixUnit>(clusCnt + 1, (matrixUnit)0));
		}
	}
}

void UClinks::addDefSeeds(shared_ptr<InputStreamer> FA, shared_ptr<Filters> fil) {
	bool contRead = true;
	int addCnt = 0; bool sync(false);
	while (contRead) {
		shared_ptr<DNA> tmpDNA = FA->getDNA(contRead, 0, sync);
		if (tmpDNA == NULL) { break; }
		shared_ptr<DNAunique> tmp = make_shared<DNAunique>(tmpDNA, -1);
		//delete tmpDNA;
		//second pair
		//shared_ptr<DNA> tmp2 = FA->getDNA(contRead, 1);
		string segs2 = tmp->getID();
		string oriClKey = segs2;
		//remove ;size= argument
		size_t idx = segs2.find(";size=");
		segs2 = segs2.substr(0, idx);
		//also remove sample identifier in string
		removeSampleID(segs2, SEP);

		//cluster should not exist, test
		std::unordered_map<string, int>::iterator itCL;
		itCL = seq2CI.find(segs2);

		if (itCL == seq2CI.end()) {
			//not found in known clusters.. create entry
			bestDNA.push_back(tmp);
			bestDNA2.push_back(NULL);
			oriKey.push_back(oriClKey);
			bestPID.push_back(100.f);
			bestLEN.push_back(tmp->length());
			clusCnt = (int)bestDNA.size() - 1;
			seq2CI[segs2] = clusCnt;
		}
		else {//cluster exists
			cerr << "Found double ID " << segs2 << endl << "Aborting.." << endl;
			exit(74);
		}
		addCnt++;
	}
	int siz = (int)bestDNA.size();
	//resize complete matrix to take these up
	
	for (uint i = 0; i < OTUmat.size(); i++){
		OTUmat[i].resize(siz, (matrixUnit)0);
	}
}

void UClinks::besterDNA(const vector<int> curCLIDpre, shared_ptr<DNAunique> tdn1, shared_ptr<DNA> tdn2, shared_ptr<Filters> fil) {
	bool checkBC = true; 
	int TagIdx(-2);
	if (tdn2 != NULL) {//fix for pairs assuming midSeqs
		checkBC = false; //TagIdx = 0;
	}
	matrixUnit matchSiz = (matrixUnit)curCLIDpre.size();
	if ( matchSiz>1 ) {
		for ( int k = 0; k < matchSiz; k++ ) {
			add2OTUmat(tdn1, curCLIDpre[k], matchSiz);
		}
		return;
	}
	int curCLID = (int) curCLIDpre[0];
	add2OTUmat(tdn1, curCLID,1);
	if (SeedsAreWritten){
		//delete tdn1; if (tdn2 != NULL){ delete tdn2; }
		return;
	}
	//general routine, matching tdn found
	if (bestDNA[curCLID] == NULL) { //just fill with current DNA
		if (tdn2 == NULL) {
			if ( fil->doReversePrimers() && !fil->check(tdn1, true, CurSetPair, TagIdx) ) {
				return;//delete tdn1; 
			}
			bestDNA[curCLID] = tdn1;
			bestPID[curCLID] = tdn1->getTempFloat();
			bestLEN[curCLID] = tdn1->length();
		} else {
			if (fil->doReversePrimers() && !fil->check(tdn1, true, CurSetPair, TagIdx)) {
				 return;//delete tdn2; delete tdn1; 
			}
			bestDNA[curCLID] = tdn1;
			bestPID[curCLID] = tdn1->getTempFloat();
			bestLEN[curCLID] = tdn1->length() + tdn2->length();
			fil->check(tdn2, true, 1, TagIdx);
			bestDNA2[curCLID] = tdn2;
		}
	}//already a candidate sequence? check who is better..
	else if (
		fil->betterSeed(tdn1, tdn2, bestDNA[curCLID], bestDNA2[curCLID], bestPID[curCLID], bestLEN[curCLID], CurSetPair, checkBC)
		){
//		delete bestDNA[curCLID];
		bestDNA[curCLID] = tdn1; 
		if (tdn2 != NULL) {
//			if (bestDNA2[curCLID] != NULL) {delete bestDNA2[curCLID];	}
			bestDNA2[curCLID] = tdn2;
		}
		if (bestPID[curCLID] < tdn1->getTempFloat()){
			bestPID[curCLID] = tdn1->getTempFloat();
		}
		uint curL = tdn1->length();
		if (tdn2 != NULL) { curL += tdn2->length(); }
		if (bestLEN[curCLID] < curL){
			bestLEN[curCLID] = curL;
		}
	} else {
//		delete tdn1;delete tdn2;
	}
}

void UClinks::removeSampleID(string& w, const string &SEP) {
	size_t pos = w.find(SEP);
	if (pos != std::string::npos) {
		w = w.substr(pos + SEP.length());
	}
}
void UClinks::removeSampleID(string& w, const string &SEP, string & SMplID) {
	size_t pos = w.find(SEP);
	if (pos != std::string::npos) {
		SMplID = w.substr(0,pos);
		w = w.substr(pos + SEP.length());
	}
}
void UClinks::removeSizeStr(string& w) {
	size_t idx = w.find(";size=");
	/*#ifdef matrix_sum   //not required, sanity check that is not really working out with usearch mappings
	int OTUsize = atoi( segs2.substr(idx+6).substr(0,-1).c_str() );
	#endif*/
	w = w.substr(0, idx);
}

void UClinks::writeNewSeeds(shared_ptr<MultiDNA> MD, shared_ptr<Filters> fil, bool refSeeds, bool printLnk) {
	if (!RefDBmode && refSeeds){ return; }
	//ofstream O(outf.c_str());
	int paired = fil->isPaired();
	//MD->printStorage();
	ofstream links;
	uint st (0), to ((uint) oriKey.size());
	//in case of refDB mode, need to start from point where refDBs were added in..
	if (refSeeds && RefDBmode && RefDBotuStart >= 0){
		cerr << "Writing ref DB sequences..";
		st = RefDBotuStart;
		//and also write a link between OTU_name and refSeq
		if (printLnk){
			string refLinkF = MD->leadOutFile() + ".lnks";
			links.open(refLinkF.c_str(), ios::out);
		}
	}else if (!refSeeds && RefDBmode && RefDBotuStart >= 0){
		cerr << "Writing new OTU seeds..";
		to = RefDBotuStart;
	}
	shared_ptr<DNA> d;
	for (uint i = st; i < to; i++) {
		//if (i == 6628){			int x = 0;		}
		//if check DNA was done before (remove rev primer), do it now
		if (bestDNA[i] == NULL) {
			MD->writeAllStoredDNA();
			cerr << "No seed sequence found for DNA " << oriKey[i] << ". Aborting..\n"; // " << bestDNA[i]->getOldID()<<"(
			//continue;
			exit(54);
		}
		//fil->check(bestDNA[i],true);
		string newH = oriKey[i];
		d = bestDNA[i];
		if (printLnk){
			string oriH = d->getIDshort(); removeSizeStr(oriH);
			links << newH << "\t" << oriH << endl;
		}
		//replace ID to allow for later linkup with cluster
		if (paired == 2) {
			d->setNewID(newH + ".1");
			d->setPassed(true);
			if (bestDNA2[i] != NULL) {
				MD->saveForWrite(d, 1);
				bestDNA2[i]->setPassed(true);
				bestDNA2[i]->setNewID(newH + ".2");
				MD->saveForWrite(bestDNA2[i], 2);
			} else {
				MD->saveForWrite(d, 3);
			}
		} else {
			d->setNewID(newH);
			d->setPassed(true);
			MD->saveForWrite(d, 1);
		}
	}
	MD->writeAllStoredDNA();
	SeedsAreWritten = true;
	if (printLnk){ links.close(); }
	cerr << "Done" << endl;
}
void UClinks::printStats(ostream& os){
	if (oriKey.size() == 0) {
		os << "No OTU's to re-seed\n";
		return;
	}
	uint numCls = 0;
	float avgQ(0.f), minQ(10000.f), maxQ(0.f),
		avgA(0.f), minA(10000.f), maxA(0.f),
		avgS(0.f), minS(10000.f), maxS(0.f);
	uint avgL = 0, minL(10000), maxL(0);
	
	//uint div(0);// = (uint)oriKey.size();
	vector<uint> lengths;
	vector<float> quals,accums,sims;
	//oriKey.size() - not if there are only refs
	uint to((uint)oriKey.size());
	if (RefDBmode && RefDBotuStart >= 0){
		to = RefDBotuStart;
	}

	for (uint i = 0; i < to; i++){
		if (bestDNA[i] == NULL){ continue; }
		float curQ = bestDNA[i]->getAvgQual();
		if (bestDNA2[i] != NULL) {
			curQ += bestDNA2[i]->getAvgQual(); curQ /= 2.f;
		}
		uint curL = (uint) bestDNA[i]->length();
		if (bestDNA2[i] != NULL) {
			curL += bestDNA2[i]->length(); 
		}

		if (curL < minL){ minL = curL; }
		if (curL > maxL){ maxL = curL; }
		avgL += curL;
		lengths.push_back(curL);
		
		if (curQ < 1) {//no new Seed found, default seed
			continue;
		}
		float sc = bestDNA[i]->getTempFloat();

		if (curQ < minQ){ minQ = curQ; }
		if (curQ > maxQ){ maxQ = curQ; }
		avgQ += curQ;

		float curA = (float)bestDNA[i]->getAccumError();
		if (bestDNA2[i] != NULL) {
			curA += (float) bestDNA2[i]->getAccumError(); 
		}

		quals.push_back(curQ);
		accums.push_back(curA); 
		if (curA < minA){ minA = curA; }
		if (curA > maxA){ maxA = curA; }
		avgA += curA;	

		sims.push_back(sc);
		if (sc < minS){ minS = sc; }
		if (sc > maxS){ maxS = sc; }
		avgS += sc;
		numCls++; 
	}
	os << "Found " << numCls << " seeds of " << oriKey.size() << " OTU's in " << uclines << " mappings." << endl;

	//sort vectors
	std::sort(lengths.begin(), lengths.end());
	std::sort(quals.begin(), quals.end());
	std::sort(accums.begin(), accums.end());
	std::sort(sims.begin(), sims.end());
	os << "Stats of Seed sequences (0th/10th/50th/90th/100th) percentile:\n";
	if (lengths.size() > 0) { os << "\n     - Seq Length :   " << minL << "/" << calc_median2(lengths, 0.1f) << "/" << calc_median2(lengths, 0.5f) << "/" << calc_median2(lengths, 0.9f) << "/" << maxL; }
	if (quals.size() > 0) { os << "\n     - Quality :      " << minQ << "/" << calc_median2(quals, 0.1f) << "/" << calc_median2(quals, 0.5f) << "/" << calc_median2(quals, 0.9f) << "/" << maxQ; }
	if (accums.size() > 0) { os << "\n     - Accum. Error : " << minA << "/" << calc_median2(accums, 0.1f) << "/" << calc_median2(accums, 0.5f) << "/" << calc_median2(accums, 0.9f) << "/" << maxA;	}
	if (sims.size() > 0) {	os << "\n     - Sim2Consensus: " << minS << "/" << calc_median2(sims, 0.1f) << "/" << calc_median2(sims, 0.5f) << "/" << calc_median2(sims, 0.9f) << "/" << maxS;		}
	os << endl;
}
