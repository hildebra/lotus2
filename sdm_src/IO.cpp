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


#include "IO.h"
//#include <map>
void threadAnalyzeDNA(shared_ptr<DNA> tdn, shared_ptr<MultiDNA> MD,int thrCnt){
	//if (threadActive){
	//	threads[thrCnt].join();
	//}
	//threads[thrdsCnt] = 

	//auto f1 = std::async(&MultiDNA::analyzeDNA,this,tdn);
	int tagIdx(-2);
	MD->analyzeDNA(tdn,thrCnt,-1,tagIdx);
	MD->saveForWrite(tdn);
}
/*void trippleThreadAnalyzeDNA(shared_ptr<MultiDNA> MD, shared_ptr<DNA> tdn,shared_ptr<DNA> tdn2,
							 shared_ptr<DNA> MIDseq,bool changePHead){//,int thrCnt){


	int thrCnt = 0;
	vector<bool> chs = MD->analyzeDNA(tdn,tdn2,MIDseq,changePHead,thrCnt);

	if (chs[0] && chs[1]){
		MD->saveForWrite(tdn,1);
		MD->saveForWrite(tdn2,2);
	} else if (chs[0]){
		MD->saveForWrite(tdn,3);
		MD->getFilters(thrCnt)->colStats[0].singleton++;
	} else if (chs[1]){
		MD->saveForWrite(tdn2,4);
		MD->getFilters(thrCnt)->colStats[1].singleton++;
	}

}*/

void read_single(OptContainer& cmdArgs, shared_ptr<MultiDNA> MD, shared_ptr<InputStreamer> IS){
	//output files
#ifdef _THREADED
	int Nthrds = atoi(cmdArgs["-threads"].c_str()) -2 ;
	int thrCnt = 0;
	MD->setSubfilters(Nthrds);
	bool threadActive(false);
	bool writeThread(false);
	vector<std::thread> threads( 0 );
	if (Nthrds>=0){
		MD->createWriteThread();writeThread=true;
		threads.resize(Nthrds);
	}
#endif
	shared_ptr<Filters> curFil = MD->getFilters();
	bool cont(true); bool sync(false);
	while (cont){
		shared_ptr<DNA> tdn1 = IS->getDNA(cont,0,sync);
		if (tdn1 == NULL) { 
#ifdef DEBUG
			cerr << "NULL read returned" << endl;
#endif
			break; 
		}

		/*if (!tdn1->isPassed()){
		MD->addNoHeadDNA(tdn1);
		tdn1 = tdn2; tdn2 = new DNA("","");
		continue;
		}*/
#ifdef _THREADED
		if (Nthrds>0){
			if (threadActive){
				threads[thrCnt].join();
			}

			//threadAnalyzeDNA(MD,tdn);
			threads[thrCnt] = std::thread(threadAnalyzeDNA,tdn,MD,thrCnt);
			thrCnt++;
			if (thrCnt >= Nthrds){
				thrCnt=0;
				threadActive=true;
			}
		} else { //single Core
			MD->analyzeDNA(tdn);
			MD->saveForWrite(tdn);
		}
#else
		curFil->sTotalPlus(0); 
		int tagIdx(-2);
		MD->analyzeDNA(tdn1,-1,-1,tagIdx);
		//here BC has to be correctly set within DNA object
		if (tagIdx == -1 ) {
			tdn1->setBarcodeDetected(false); 
		}
		MD->depPrep(tdn1,NULL);
		curFil->write2Demulti(tdn1, 0,MD->getfastQoutVer());

		if (!MD->saveForWrite(tdn1)) {
			cont = false;
			break;
		}
		
#endif
		//if (tdn!=NULL && ch1 != tdn->isPassed()){cerr<<"isPassed is != ch1! Aborting..\n";exit(12);}
	}
#ifdef _THREADED
	if (threadActive){
		for (uint i=0; i<threads.size();i++){
			threads[i].join();
		}
	}
#endif
	MD->closeOutStreams();
}


//is called from a while loop, that reads the DNA pairs
bool read_paired_DNAready(shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2, shared_ptr<DNA> MIDseq, bool MIDuse, shared_ptr<MultiDNA> MD, int& revConstellation) {
	
	if (tdn == NULL) { return true; } //|| tdn->length()==0

	shared_ptr<Filters> curFil = MD->getFilters();

	//register read at all with stat counter:
	curFil->sTotalPlus(0); curFil->sTotalPlus(1);

	//prep some variables
	int BCoffs = curFil->getBCoffset();
	bool checkBC2ndRd = curFil->checkBC2ndRd();
	bool dualBCs = curFil->doubleBarcodes();
	bool doBCsAtAll = curFil->doBarcodes();
	bool checkReversedRead = curFil->checkRevRd();


	int tagIdx(-2); int tagIdx2(-2);
	string presentBC(""); int c_err(0);
	bool isReversed(false);//was a reversion detected?

	if (MIDuse && MIDseq!= 0) { 
		tagIdx = curFil->cutTag(MIDseq, presentBC, c_err, true); 
//		delete MIDseq; 
		tdn->setBCnumber(tagIdx, BCoffs);
	}
	if (checkBC2ndRd ) {
		if (!dualBCs) {
			bool revT = false;
			bool Pr1 = curFil->findPrimer(tdn, 0, false, 0);
			bool Pr2 = curFil->findPrimer(tdn2, 0, false, 0);
			tagIdx = curFil->findTag(tdn, presentBC, c_err, true);
			tagIdx2 = curFil->findTag(tdn2, presentBC, c_err, true);
			if ( true &&checkReversedRead && (tagIdx2 < 0 && tagIdx < 0) ) {
				tdn->reverse_transcribe(); tdn2->reverse_transcribe();
				Pr1 = curFil->findPrimer(tdn, 0, false, 0);
				Pr2 = curFil->findPrimer(tdn2, 0, false, 0);
				tagIdx = curFil->findTag(tdn, presentBC, c_err, true);
				tagIdx2 = curFil->findTag(tdn2, presentBC, c_err, true);
				revT = true;
			}
			if ((tagIdx2 >= 0 && tagIdx < 0 && !Pr1) || (Pr2 && !Pr1)) { //swap first & second read
				swap(tdn, tdn2);
				revConstellation++;
			}
			/*else if (tagIdx2 < 0 && tagIdx < 0) {
				int x = 0;
			}*/
			if (revT) {
				tdn->reverse_transcribe(); tdn2->reverse_transcribe();
			}
		}
		tagIdx2 = -2; tagIdx = -2;
		tdn2->setpairREV();		tdn->setpairFWD();
	}

	//tdn->reverse_transcribe();
	MD->analyzeDNA(tdn, -1, 0, tagIdx);
	//tdn->matchSeqRev
	bool ch1(false); if (tdn != NULL) { ch1 = tdn->isPassed(); }
	bool ch2(false); bool ch2n(false);

	//this is all about barcodes..
	if (checkReversedRead  && tdn != NULL && tagIdx < 0) {
		if (!MIDuse) { tagIdx = -2; }
//		curFil->sTotalMinus(0);
		tdn->reverse_transcribe();
		MD->analyzeDNA(tdn, -1, 0, tagIdx);
		ch1 = tdn->isPassed();
		isReversed = ch1;
		if (!isReversed) {//reset
			tdn->reverse_transcribe();
		}
	}

	//test for reverse complemented reads (mohammad samples), when BC not found (NOT dual BC)
	//in that case, this is the first read
	if (false &&checkBC2ndRd && tagIdx < 0 && tdn2 != NULL) {// && !tdn->getBarcodeDetected() ) {
											   //tdn2->reverse_transcribe();
		if (!MIDuse) { tagIdx = -2; }
//		curFil->sTotalMinus(0);
		MD->analyzeDNA(tdn2, -1, 0, tagIdx);
		ch2n = tdn2->isPassed();
		if (!ch2n && checkReversedRead) {
			if (!MIDuse) { tagIdx = -2; }
//			curFil->sTotalMinus(0);
			tdn2->reverse_transcribe();
			MD->analyzeDNA(tdn2, -1, 0, tagIdx);
			ch2n = tdn2->isPassed();
			isReversed = ch2n;
			if (!ch2n) { tdn2->reverse_transcribe(); }//reset to ori
		}
		if (ch2n) {//passed ch2 through BC filter, now really reverse
				   //1st, now 2nd pair
			tdn2->setpairFWD();
			ch1 = ch2n;
			if (tdn != NULL) {
				//tdn->reverse_transcribe(); 
				tdn->setpairREV();
				tdn->reset();
				if (!dualBCs) { tagIdx2 = tdn->getBCnumber(); } // no 2nd BC, thus no BC search in 2nd read
				MD->analyzeDNA(tdn, -1, 1, tagIdx2);
				ch2 = tdn->isPassed();
			}
			swap(tdn, tdn2);
			revConstellation++;
		}

	}



	//if ( ch1 ) {	cerr << cnt << " \n";	}
	//normal case for check 2nd read
	if (!ch2 && tdn2 != NULL) { //ch1&&
								//tdn2->setBCnumber(tdn->getBCnumber());
		if (doBCsAtAll && !dualBCs) { //only check in read1 for BC, if not dual BCing!!
			tagIdx2 = tdn->getBCnumber();  // no 2nd BC, thus no BC search in 2nd read
			if (tagIdx2 >= 0) {
				tagIdx2 -= BCoffs;
			}else if (tagIdx2 < -1) {//something wrong with BCoffs
				cerr << "tagidx2 wrongly truncated to " << tagIdx2 << endl;
			}
		}
		if (isReversed) { tdn2->reverse_transcribe(); }
		MD->analyzeDNA(tdn2, -1, 1, tagIdx2);
		ch2 = tdn2->isPassed();
	}

	//set up BC in DNA header
	//remember that dual BCs are only valid after this step!
	if (dualBCs) {
		//tagIdx2 = -2; //reset just to be sure
		curFil->dblBCeval(tagIdx, tagIdx2, presentBC, tdn, tdn2);
		c_err = -1;

		//check a second time that barcode was correctly identified, just to be double sure...
		if (tagIdx != tagIdx2 || tdn->getBCnumber() != tdn2->getBCnumber()) {
			cerr << "Unequal BC numbers:" << tagIdx << " : " << tagIdx2 << "; in object: " << tdn->getBCnumber() << " : " << tdn2->getBCnumber() << endl;
			cerr << "In read:" << tdn->getID() << endl;
			exit(835);
		}
	}
	else if (tagIdx >= 0) {
		if (MIDuse&&ch1) { curFil->BCintoHead(tagIdx, tdn, presentBC, c_err, true); }
		else { curFil->setBCdna(tagIdx, tdn); }
		if (ch2) { curFil->BCintoHead(tagIdx, tdn2, presentBC, c_err, true); }
	}

	if (tagIdx == -1 || tagIdx2 == -1) {
		if (ch1) {
			tdn->setBarcodeDetected(false);
		}
		if (ch2) {
			tdn2->setBarcodeDetected(false);
		}
	}

	//demultiplex write? do this first before DNA is deleted..
	if (curFil->Demulti2Fls()) {
		curFil->write2Demulti(tdn, 0, MD->getfastQoutVer());
		curFil->write2Demulti(tdn2, 1, MD->getfastQoutVer());
	}


	//at this point the tagIDX *MUST* be correctly set + BCoffset (in the DNA object, tagIDX doesn;t matter)
	MD->depPrep(tdn, tdn2);
	MD->writeNonBCReads(tdn, tdn2);

	int idx1 = 1; int idx2 = 2;
	if (ch1 && !ch2) {
		idx1 = 3; idx2 = 4;
		if (tdn2 != NULL) { tdn2->failed(); }
//		delete tdn2;
	}
	else if (ch2 && !ch1) {
		idx2 = 4; idx1 = 3;
		if (tdn != NULL) { tdn->failed(); }
//		delete tdn;
	}
	else if (!ch1 && !ch2){ //nothing passes
		if (tdn != NULL) { tdn->failed(); }
		if (tdn2 != NULL) { tdn2->failed(); }
//		delete tdn; delete tdn2;
	}

	//save for later .. and collect stats
	if (!MD->saveForWrite(tdn, idx1) ||
		!MD->saveForWrite(tdn2, idx2)) {
		return false;
	}
	return true;
}

bool read_paired(OptContainer& cmdArgs, shared_ptr<MultiDNA> MD, shared_ptr<InputStreamer> IS, bool MIDuse) {
	DNAmap oldMIDs;
	bool fqHeadVer(true);
	shared_ptr<DNA> MIDseq(NULL);
	bool syncedMID(MD->getFilters()->synRdPairs());

	bool sync2pair(true);
	DNAmap pair1rem,pair2rem,MIDrem;

	/*if (sync2pair && MIDuse) {
		cout << "Can not sync read pairs, while explicit MID sequences are being used! (not supported, sorry)\n";
		sync2pair = false;
	} */

	//bool syncedMID = false;
#ifdef DEBUG
	cerr << "Read paired routine" << endl;
#endif	
#ifdef _THREADED
	vector<std::thread> threads( Nthrds );
	int Nthrds = atoi(cmdArgs["-threads"].c_str()) -1 ;
	int thrCnt = 0;
	bool threadActive(false);
	MD->setSubfilters(Nthrds);
	int DNAinMem(0);
#endif
	bool cont(true),cont2(true),cont3(true);
	int revConstellation(0);
	int cnt(0); 
	string tdnSh(""), tdnSh2("");
	bool switching(true); // important to keep track of this, to fix swapped read pairs

	while ( cont ) {
		
		
		bool settdnSh(false);
		shared_ptr<DNA> tdn = IS->getDNA(cont, 0, sync2pair);

		cnt++;
		if ( tdn == NULL && !cont) { break; }
		//tagIdx = -2; tagIdx2 = -2;
		shared_ptr<DNA> tdn2 = IS->getDNA(cont2, 1, sync2pair);//read_fastq_entry(fna2,fastQver,minQScore,lnCnt);
		if (!sync2pair &&  !cont2 && cont ) {
			cerr << "Second provided file has not the same number of entries as first file.\n";
			exit(5);
		}
		//syn 2nd pair
		if (sync2pair) {
			if (!settdnSh) { tdnSh = tdn->getIDshort(); settdnSh = true; }
			if (tdn2 != NULL) { tdnSh2 = tdn2->getIDshort(); }else { tdnSh2 = ""; }
			DNAmap::iterator search;
			if (tdnSh2 != tdnSh) {//something wrong at all? ini search on old pair2
				if (switching) {
					search = pair2rem.find(tdnSh);
					if (search != pair2rem.end()) { pair2rem[tdnSh2] = tdn2; tdn2 = search->second; pair2rem.erase(search); tdnSh2 = tdn2->getIDshort(); switching = !switching;}
				}	else {
					search = pair1rem.find(tdnSh2);
					if (search != pair1rem.end()) { pair1rem[tdnSh] = tdn; tdn = search->second; pair1rem.erase(search); tdnSh = tdn->getIDshort();  switching = !switching;}
				}
			}
			while (tdnSh2 != tdnSh) {// still wrong? search in old unmatched reads, switching between 1/2
				if (switching) {
					search = pair1rem.find(tdnSh2); pair1rem[tdnSh] = tdn;
					if (search != pair1rem.end()) {  
						tdn = search->second; pair1rem.erase(search);
					}else {//nothing? try getting new reads, maybe there is a match here
						tdn = IS->getDNA(cont, 0, sync2pair);
					}
					if (tdn != NULL) { tdnSh = tdn->getIDshort(); } else {tdnSh = "";}
				} else {
					search = pair2rem.find(tdnSh); pair2rem[tdnSh2] = tdn2;
					if (search != pair2rem.end()) {  tdn2 = search->second; pair2rem.erase(search); 
					} else {
						tdn2 = IS->getDNA(cont2, 1, sync2pair);
					}
					if (tdn2 != NULL) { tdnSh2 = tdn2->getIDshort(); }	else { tdnSh2 = ""; }
				}
				
				//read 1 / read 2
				switching = !switching;
			}
		}
		if (cnt % 50 == 0) { switching = !switching; } //add some randomness to the process.
		//security check that pairs are in sync
		if (!sync2pair && cnt % 1000 == 0) {//default check for synced reads, no matter what
			if (!settdnSh) { tdnSh = tdn->getIDshort(); settdnSh = true; }
			if (!tdn2->sameHead(tdnSh) ) { cerr << "WARNING: read pairs out of sync (" << cnt << "): " << tdnSh << " " << tdnSh2 << endl; }
		}
		//sync mid sequence header
		if (MIDuse) {
			if (!settdnSh) {tdnSh = tdn->getIDshort(); settdnSh = true;	}
			//1st try to find in old heap
			if ( !syncedMID && (cont || cont2) ) {
				auto search = oldMIDs.find(tdnSh);
				if ( search != oldMIDs.end() ) {
					//it's in old.. give it to MIDseq, other rountines will del it
					MIDseq = (*search).second;	oldMIDs.erase(search);
				} else {
					//read some new lines, maybe here?
					MIDseq = IS->getDNA(cont3, 2, sync2pair);
					bool SameMIDHead(MIDseq->sameHead(tdnSh));
					while ( !SameMIDHead && cont3 ) {
						//current MID not matching.. stove away
						oldMIDs[MIDseq->getIDshort()] = MIDseq;	MIDseq = IS->getDNA(cont3, 2, sync2pair);
						SameMIDHead = MIDseq->sameHead(tdnSh);
					}
				}
			} else {MIDseq = IS->getDNA(cont3, 2, sync2pair);}
			MIDseq->setMIDseq(true);
			//FQ header version changes have to occur, before the MID tag is labelled
			
		}
		//check if the PE format is right
		if ( fqHeadVer ) { MD->checkFastqHeadVersion(tdn); fqHeadVer = false; }

		cont = read_paired_DNAready(tdn,tdn2, MIDseq, MIDuse, MD, revConstellation);
	}

	//check on remainders in pair2rem / pair1rem
	if (sync2pair) {
		if (pair1rem.size() > 0 && pair2rem.size() > 0) {
			cerr << "Trying to match " << pair1rem.size() << " / " << pair2rem.size() << " out-of-sync read pairs..\n";
			DNAmap::iterator search;
			for (DNAmap::iterator sr = pair1rem.begin(); sr != pair1rem.end(); sr++) {
				search = pair2rem.find(sr->first);
				if (search != pair2rem.end()) {
					cont = read_paired_DNAready(sr->second, search->second, MIDseq, MIDuse, MD, revConstellation);
					pair2rem[search->first] = NULL; pair1rem[sr->first] = NULL;
					pair2rem.erase(search); pair1rem.erase(  sr  );
				}
			}
			cerr << "Writing remaining " << pair1rem.size() << " / " << pair2rem.size() << " out-of-sync reads as singletons.\n";
		} else if (pair1rem.size() > 0 || pair2rem.size() > 0) {
			cerr << "Found " << pair1rem.size() << " / " << pair2rem.size() << " out-of-sync read pairs.\n";
		}
		for (DNAmap::iterator sr = pair1rem.begin(); sr != pair1rem.end(); sr++) {
			cont = read_paired_DNAready(sr->second, NULL, MIDseq, MIDuse, MD, revConstellation);
			pair1rem[sr->first] = NULL;//object needs to remain in mem
		}
		for (DNAmap::iterator sr = pair2rem.begin(); sr != pair2rem.end(); sr++) {
			cont = read_paired_DNAready(NULL, sr->second, MIDseq, MIDuse, MD, revConstellation);
			pair2rem[sr->first] = NULL;
//			pair2rem.erase(sr);
		}

	}
	//close shop
	MD->revConstellationCnts(revConstellation);
	MD->closeOutStreams();
	return true;
}

bool readCmdArgs(int argc, char* argv[],OptContainer& cmdArgs){
	if (argc%2!=1){
		cerr<<"It seems command line arguments were not passed in pairs. Aborting.\n";
		exit(666);
	}
	for (int i=1; i<argc; i+=2){ //parsing of cmdline args
		string theNxtSt = string(argv[i+1]);
		if (theNxtSt[0] != '-'){
			cmdArgs[string(argv[i])] = theNxtSt;
		} else {
			cmdArgs[string(argv[i])] = "T";
		}
	}
	if (cmdArgs.find("-i_MID_fastq") == cmdArgs.end()) {
		cmdArgs["-i_MID_fastq"] = "";
	}
	//set to default (empty)
	if (cmdArgs.find("-OTU_fallback") == cmdArgs.end()){
		cmdArgs["-OTU_fallback"] = "";
	}
	if (cmdArgs.find("-otu_matrix") == cmdArgs.end()) {
		cmdArgs["-otu_matrix"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts") == cmdArgs.end()) {//.ADD
		cmdArgs["-ucAdditionalCounts"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts1") == cmdArgs.end()) {//.REST
		cmdArgs["-ucAdditionalCounts1"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts_refclust") == cmdArgs.end()) {//.ADDREF
		cmdArgs["-ucAdditionalCounts_refclust"] = "";
	}
	if (cmdArgs.find("-optimalRead2Cluster_ref") == cmdArgs.end()) {//.ADDREF
		cmdArgs["-optimalRead2Cluster_ref"] = "";
	}
	//just for debuggin purposes: write out all seqs, where no BC can be detected..
	if (cmdArgs.find("-o_fastq_noBC") == cmdArgs.end()) {
		cmdArgs["-o_fastq_noBC"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts_refclust1") == cmdArgs.end()) {//.RESTREF
		cmdArgs["-ucAdditionalCounts_refclust1"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts_refclust1") == cmdArgs.end()) {//.RESTREF
		cmdArgs["-ucAdditionalCounts_refclust1"] = "";
	}
	if (cmdArgs.find("-XfirstReads") == cmdArgs.end()) {
		cmdArgs["-XfirstReads"] = "";
	}
			//filter sequence file for a specific subset of sequences
	//these arguments can only occur together
	if (cmdArgs.find("-specificReads") == cmdArgs.end()) {
		cmdArgs["-specificReads"] = "";
	} else if (cmdArgs.find("-excludeFile") == cmdArgs.end()) {
		cmdArgs["-excludeFile"] = "";
	}
	if (cmdArgs.find("-onlyPair") == cmdArgs.end()) {
		cmdArgs["-onlyPair"] = "";
	}
	if (cmdArgs.find("-uparseVer") == cmdArgs.end()) {
		cmdArgs["-uparseVer"] = "";
	}


	if (cmdArgs.find("-i_path")  == cmdArgs.end()){ //ok files are not given in mapping file
		//check if dna and qual are passed
		if (cmdArgs.find("-i") != cmdArgs.end()) {
			string fmt = detectSeqFmt(cmdArgs["-i"]);
			if (fmt == "empty") {
				//exit(0);
				cerr << "Only empty input files\n";
			}
			cmdArgs[fmt] = cmdArgs["-i"];
		}
		if (cmdArgs.find("-i_fastq")  == cmdArgs.end()){ // fasta + quality format
			if (cmdArgs.find("-i_fna")  == cmdArgs.end()){
				cerr<<"You did not supply a fasta file. \nPlease give the path to your fasta file as command line argument:\n  -i_fna <yourFastaFile>\n";
				exit(2);
			}
			if (cmdArgs.find("-i_qual")  == cmdArgs.end()){
				string newQ = cmdArgs["-i_fna"];
				int pos = (int)newQ.find_last_of(".");
				newQ = newQ.substr(0,pos);
				newQ += string(".qual");
				fstream fin;
				fin.open(newQ.c_str(),ios::in);
				if( fin.is_open() )	{
					cerr<<"Using quality file: "<<newQ <<endl;
				} else if ((cmdArgs.find("-number")!=  cmdArgs.end() && cmdArgs["-number"] =="T")||
					(cmdArgs["-specificReads"] != "")) {
					cmdArgs["-i_qual"] = "";
				} else {
					cerr<<"You did not supply a quality file. \nPlease give the path to your quality file as command line argument:\n  -i_qual <PathToQualityFile>\n";
					newQ = "";
					//fin.close();	exit(2);
				}
				fin.close();
				cmdArgs["-i_qual"] = newQ;
			}
		}
		//auto create output file name
		if (cmdArgs.find("-o_fna")  == cmdArgs.end()){
			if (cmdArgs.find("-o_fastq")  == cmdArgs.end()){
				//cmdArgs["-o_fna"] = cmdArgs["-i_fna"]+string(".sdm");
				//cerr<<"Writing output fasta into "<<cmdArgs["-o_fna"]<<endl;
				cerr << "No output file will be written\n";
			}
		} else {
			if (cmdArgs.find("-o_fastq")  != cmdArgs.end()){
				cerr<<"\"-o_fna\" was over-writen by \"-o_fastq\"\n";
				cmdArgs["-o_fna"] = "";
			}
		}
	} else {
		if (cmdArgs.find("-o_fna")  == cmdArgs.end() && cmdArgs.find("-o_fastq")  == cmdArgs.end()){
			cerr<<"Please give an output file (\"-o_fna\" || \"-o_fastq\") if you use sdm \"-i_path\" option.\n  Aborting..\n";
			exit(2);
		}
	}

	/*	if (cmdArgs.find("-map")  == cmdArgs.end()){
	cerr<<"You did not supply a mapping file. \nPlease give the path to your mapping file as command line argument:\n  -map <PathToMappingFile>\n";
	exit(2);
	}  */
	if (cmdArgs.find("-o_qual")  == cmdArgs.end()){
		cmdArgs["-o_qual"] = "";
	} else {
		if (cmdArgs.find("-o_fastq")  != cmdArgs.end()){
			cerr<<"\"-o_qual\" was over-writen by \"-o_fastq\"\n";
			cmdArgs["-o_qual"] = "";
		}
	}
	if (cmdArgs.find("-options")  == cmdArgs.end()){
		cmdArgs["-options"] = string("sdm_options.txt");
	}
	if (cmdArgs.find("-threads")  == cmdArgs.end()){
		cmdArgs["-threads"] = "1";
	}
	if (cmdArgs.find("-log")  == cmdArgs.end()){
		string ofile1 = cmdArgs["-o_fna"];
		if (ofile1==""){ofile1 = cmdArgs["-o_fastq"];}
		vector<string> tvec = splitByComma(ofile1,false); 
		ofile1 = tvec[0];
		//remove file ending
		unsigned int pos = (unsigned int) ofile1.find_last_of(".");
		if (pos != string::npos){ofile1 = ofile1.substr(0,pos);	}
		if (tvec.size()==2){
			ofile1+= "_" + getFileNoPath(tvec[1]);
			pos = (unsigned int) ofile1.find_last_of(".");
			if (pos != string::npos){ofile1 = ofile1.substr(0,pos);	}
		}
		cmdArgs["-log"] = ofile1 + string(".log");
	}
	string ofile1 = cmdArgs["-log"];
	ofile1.find_last_of(".log");
	size_t logPos = ofile1.find_last_of(".");
	if (logPos != std::string::npos){
		ofile1 = ofile1.substr(0,logPos);
	}
	if (cmdArgs.find("-length_hist")  == cmdArgs.end()){
		cmdArgs["-length_hist"]  = ofile1 + string("_lenHist.txt");
	}
	if (cmdArgs.find("-qual_hist")  == cmdArgs.end()){
		cmdArgs["-qual_hist"]  = ofile1 + string("_qualHist.txt");
	}
	//-length_hist   -qual_hist

	if (cmdArgs.find("-sample_sep")  == cmdArgs.end()){
		cmdArgs["-sample_sep"] = DEFAULT_BarcodeNameSep;
	} else 	if (cmdArgs["-sample_sep"]==""){
		cerr<<"Invalid sample separator (empty).\nAborting..\n";exit(82);
	}


	if (cmdArgs.find("-o_qual_offset") == cmdArgs.end()) {
		cmdArgs["-o_qual_offset"] = DEFAULT_output_qual_offset;
	}
	
	if (cmdArgs.find("-ignore_IO_errors") == cmdArgs.end()) {
		cmdArgs["-ignore_IO_errors"] = DEFAULT_ignore_IO_errors;
	} else if (cmdArgs["-ignore_IO_errors"] != "0" && cmdArgs["-ignore_IO_errors"] != "1") {
		cerr << "Argument \"ignore_IO_errors\" can only be \"1\" or \"0\". Instead it has value: " << cmdArgs["-ignore_IO_errors"] << endl;
		exit(323);
	}
	if (cmdArgs.find("-o_dereplicate") == cmdArgs.end()) {
		cmdArgs["-o_dereplicate"] = "";
	}
	if (cmdArgs.find("-derep_map") == cmdArgs.end()) {
		cmdArgs["-derep_map"] = "";
	}

	

	//if (cmdArgs.count("-i_fna")==0){}

	return true;
}


/*******************************************
*				read_fasta   			   *
*******************************************

void openOutFiles(string files, string fmt, string xtr){
	ofstream fnaOut;

	vector<string> tfnaout(0);
	if (files.find(",") != string::npos){
		tfnaout = splitByCommas(files);
	} else {
		tfnaout.push_back(files);
	}
	bool multiple = tfnaout.size() > 1;
	string xtr2 = "";
	if (multiple){xtr2 = "paired ";}
	for (uint i =0; i< tfnaout.size(); i++){
		fnaOut.open ( tfnaout[i].c_str(),ios_base::out);
		if (!fnaOut){	cerr<<"Could not open "<<xtr2<<xtr<<fmt<<" output file "<<i<<": "<<tfnaout[i]<<endl;exit(4);	}
		fnaOut.close();
		if (multiple){//also singletonfiles
			string tmp = tfnaout[0]+SingletonFileDescr;
			fnaOut.open(tmp.c_str(),ios_base::out);
			if (!fnaOut){	cerr<<"Could not open Singleton "<<xtr<<fmt<<" output file "<<i<<": "<<tmp<<endl;exit(4);	}
			fnaOut.close();
		}
	}

}

void prepareOutFiles(OptContainer& cmdArgs){
	ofstream fnaOut;
	
	//additional output files (secondary filtering)
	if (cmdArgs.find("-o_fastq2")  != cmdArgs.end() && cmdArgs["-o_fastq2"] != ""){
		openOutFiles(cmdArgs["-o_fastq2"],"fastq","add ");
	}
	if (cmdArgs.find("-o_fna2")  != cmdArgs.end() && cmdArgs["-o_fna2"] != ""){
		openOutFiles(cmdArgs["-o_fna2"],"fna","add ");
	}

	//fastq
	if (cmdArgs["-o_fna"]=="" && cmdArgs["-o_fastq"] != ""){
		openOutFiles(cmdArgs["-o_fastq"],"fastq","");
		return;
	}

	//fasta output
	vector<string> tfnaout = splitByComma(cmdArgs["-o_fna"],false);
	for (unsigned int i=0; i<tfnaout.size();i++){
		fnaOut.open(tfnaout[i].c_str(),ios_base::out);
		if (!fnaOut){	cerr<<"Could not open Fasta "<< i<<" output file "<<tfnaout[0]<<endl;exit(4);	}
		fnaOut.close();
	}
	if (tfnaout.size()==2){//PE - singleton file
		string tmp = tfnaout[0]+SingletonFileDescr;
		fnaOut.open(tmp.c_str(),ios_base::out);
		if (!fnaOut){	cerr<<"Could not open Singleton Fasta output file "<<tmp<<endl;exit(4);	}
		fnaOut.close();
		tmp = tfnaout[1]+SingletonFileDescr;
		fnaOut.open(tmp.c_str(),ios_base::out);
		if (!fnaOut){	cerr<<"Could not open Singleton Fasta output file "<<tmp<<endl;exit(4);	}
		fnaOut.close();
	}
	if (cmdArgs["-o_qual"] != ""){
		vector<string> tqout = splitByComma(cmdArgs["-o_qual"],false);
		for (unsigned int i=0; i<tqout.size();i++){
			fnaOut.open (tqout[0].c_str() ,ios_base::out);
			if (!fnaOut){			cerr<<"Could not open Quality "<<i<<" output file "<<tqout[0]<<endl;			exit(4);		}
			fnaOut.close();
		}
		if (tqout.size()==2){//PE - singleton file
			string tmp = tqout[0]+SingletonFileDescr;
			fnaOut.open(tmp.c_str(),ios_base::out);
			if (!fnaOut){	cerr<<"Could not open Singleton Quality output file "<<tmp<<endl;exit(4);	}
			fnaOut.close();
			tmp = tqout[1]+SingletonFileDescr;
			fnaOut.open(tmp.c_str(),ios_base::out);
			if (!fnaOut){	cerr<<"Could not open Singleton Quality output file "<<tmp<<endl;exit(4);	}
			fnaOut.close();
		}
	}
}
*/

//manages read in of several input files and associated primers / tags to each file
void separateByFile(shared_ptr<Filters> mainFil,OptContainer& cmdArgs){
#ifdef DEBUG
	cerr << "separateByFile"<<endl;
#endif
	vector<string> FastaF = mainFil->getFastaFiles();
	vector<string> QualF = mainFil->getQualFiles();
	vector<string> FastqF = mainFil->getFastqFiles();
	vector<string> MIDfq = mainFil->getMIDfqFiles();
	vector<string> tar;
	vector < vector<int> > idx(0);
	bool bFASTQ = true;
	//prepareOutFiles(cmdArgs);
	string path="";
	if (cmdArgs.find("-i_path")  != cmdArgs.end() && cmdArgs["-i_path"].length() > 2){
		path=cmdArgs["-i_path"] + string("/");
	}


	if (FastaF.size()>0){ //fasta way
		tar = FastaF;
		bFASTQ = false;
	} else { // fastq way
		tar = FastqF;
		if (FastqF.size()==0){
			cerr<<"No FastQ or Fasta file given.\n  Aborting..\n";
			exit(12);
		}
	}

	vector<string> uniqueFas(1,tar[0]);
	idx.push_back(vector<int> (1,0));

	for (unsigned int i=1; i<tar.size(); i++){
		bool suc = false;
		for (unsigned int j=0; j<uniqueFas.size(); j++){
			if (tar[i] == uniqueFas[j]){ //the same
				idx[j].push_back(i);
				suc=true;
				break;
			}
		}
		if (!suc){
			uniqueFas.push_back(tar[i]);
			idx.push_back(vector<int> (1,i));
		}
	}

	//unique Fas files set up.. check for their existence
	shared_ptr<InputStreamer> testFiles = 
		make_shared<InputStreamer>(!bFASTQ, mainFil->getuserReqFastqVer(), "1");
	for (unsigned int i = 0; i < uniqueFas.size(); i++) {
		int tarID = idx[i][0]; string tmp;
		string x = testFiles->setupInput(path, i, tarID, uniqueFas, FastqF, FastaF, QualF, MIDfq, mainFil->isPaired(), cmdArgs["-onlyPair"], tmp, true);
	}
//	delete testFiles;
	
	string mainFile = "", outFile = cmdArgs["-o_fna"];
	
	//prepare for Seed extension or Read subselection, if requested
	UClinks *ucl = NULL; shared_ptr<ReadSubset> RDSset; 
	shared_ptr<Dereplicate> Dere ;
	if (mainFil->doOptimalClusterSeq()){
		ucl = new UClinks(cmdArgs);
		if (cmdArgs.find("-mergedPairs") != cmdArgs.end() && cmdArgs["-mergedPairs"] == "1"){
			ucl->pairedSeqsMerged(mainFil);
		}
		else {
			mainFil->setFloatingEWin(10, 25);
		}
		//are fallback fasta sequences available?
		if (cmdArgs["-OTU_fallback"] != ""){
			shared_ptr<InputStreamer> FALL = make_shared<InputStreamer>(true, mainFil->getuserReqFastqVer(), cmdArgs["-ignore_IO_errors"]);
			FALL->setupFna(cmdArgs["-OTU_fallback"]);
			ucl->setupDefSeeds(FALL,mainFil);
		}
	}
	else if (mainFil->doSubselReads()){
		//this will select a list of reads and distribute these into multiple files
		RDSset = make_shared<ReadSubset>(cmdArgs["-specificReads"],"");
	} else if (mainFil->doDereplicate()) {
		Dere = make_shared<Dereplicate>(cmdArgs);
	}
	//needs to attach to existing file sometimes
	std::ios_base::openmode writeStatus = ios_base::out;
	bool shortStats = false;
	string shrtLog = "";


			// main loop that goes over different files
	int maxRds = mainFil->getXreads();
	int totReadsRead(0);
	for (unsigned int i=0; i<uniqueFas.size();i++ ){
#ifdef DEBUG
		cerr << "new filter in round "<<i << endl;
#endif
		if (maxRds>0 && maxRds - totReadsRead <= 0) { break; }
		shared_ptr<Filters> fil = make_shared<Filters>(mainFil, idx[i][0]);
		unsigned int tarSi = (unsigned int) idx[i].size();
		fil->allResize(tarSi);
		int tarID=-1;
		bool BC2mode = mainFil->doubleBarcodes();
		//int readsRead(0);
		
		
		for (unsigned int j=0; j<tarSi;j++){ //fill in filter
			tarID = idx[i][j];
			if (mainFil->PrimerIdx[tarID]>-1) {
				fil->addPrimerL(mainFil->PrimerL[mainFil->PrimerIdx[tarID]], j);
			}
			if (mainFil->doReversePrimers() && mainFil->PrimerIdxRev[tarID]>-1) {
				fil->addPrimerR(mainFil->PrimerR[mainFil->PrimerIdxRev[tarID]], j);
			}
			fil->Barcode[j] = mainFil->Barcode[tarID];
			if ( BC2mode ) {
				fil->Barcode2[j] = mainFil->Barcode2[tarID];
			}
			fil->SampleID[j] = mainFil->SampleID[tarID];
			fil->SampleID_Combi[j] = mainFil->SampleID_Combi[tarID];
			fil->HeadSmplID[j] = mainFil->HeadSmplID[tarID];
			if (fil->Demulti2Fls()) {
				fil->demultiSinglFiles[j] = mainFil->demultiSinglFiles[tarID];
				fil->demultiSinglFilesF[j] = mainFil->demultiSinglFilesF[tarID];
				
				//closing of ofstreams is only handled on the main object
				//mainFil->demultiSinglFiles[tarID] = vector<ofstream*> (2,NULL);
			}
		}
		fil->checkDoubleBarcode();

		if (tarID==-1){cerr<<"tar == -1. abort.\n";exit(10);}

		//initialize object to handle all input file combinations
		
		shared_ptr<InputStreamer> IS = make_shared<InputStreamer>(!bFASTQ, mainFil->getuserReqFastqVer(), cmdArgs["-ignore_IO_errors"]);
		if (tarSi < 2 && uniqueFas.size() > 1) {
			IS->atFileYofX(i + 1, (uint)uniqueFas.size(), tarSi);
		}
		string mainFileShort = "";
		mainFile = IS->setupInput(path, i, tarID, uniqueFas, FastqF, FastaF, QualF, MIDfq, fil->isPaired(), cmdArgs["-onlyPair"], mainFileShort, false);
		if (!IS->qualityPresent()) {
			fil->deactivateQualFilter();
			cerr << "\n*********\nWarning:: Quality file is not present.\nRecommended to abort demultiplexing.\n*********\n\n";
		}
		fil->BarcodePreStats();
		fil->checkDoubleBarcode();
		fil->checkDoubleSampleIDHead();


		if (mainFil->doOptimalClusterSeq()){
			ucl->findSeq2UCinstruction(IS,bFASTQ,mainFil);
			continue;
		} 


#ifdef DEBUG
		cerr << "Setting up output" << endl;
#endif
		//MultiDNA MD = MultiDNA(&fil, cmdArgs, writeStatus, RDSset);
		shared_ptr<MultiDNA> MD = make_shared<MultiDNA>(fil, cmdArgs, writeStatus, RDSset);
		fil->setMultiDNA(MD);
		if (maxRds > 0) { MD->setReadLimit(maxRds - totReadsRead); }
		writeStatus = ofstream::app;
		//prepare for BC checking (rev/fwd)
		if (fil->doDemultiplex()){
			MD->setBCfixed(false, true);  
			if (MD->isPEseq() == 2) { MD->setBCfixed(false, false); }
		}
		if (cmdArgs.find("-oneLineFastaFormat") != cmdArgs.end() && cmdArgs["-oneLineFastaFormat"] == "1") {
			MD->setOneLinerFastaFmt(true);
		}
		//cout << Dere->Nms_size() << " DEBCs\n";
		MD->attachDereplicator(Dere);
		//only pull out a subset of sequences
		if (mainFil->doSubselReads()) {
			if (cmdArgs.find("-mocatFix") != cmdArgs.end()) {
				cerr << "MOCAT fix appplies\n";
				RDSset->findMatches(IS, MD, true);
			}else{
				RDSset->findMatches(IS, MD, false);
			}
			//delete MD;
			continue;
		}


#ifdef DEBUG
		cerr << "Processing reads" << endl;
#endif





		//**********************
		//heavy reading routine
		//**********************
		if (MD->isPEseq() == 2){
			//read_paired(cmdArgs,MD,IS);
			while ( !read_paired(cmdArgs, MD, IS, IS->hasMIDseqs()) ) {
				//reset output files to previous state
				MD->resetOutFilesAndFilter();
			}
		} else {
			read_single(cmdArgs,MD,IS);
		}




#ifdef DEBUG
		cerr << "All read processed" << endl;
#endif
		outFile = MD->leadOutFile();
//		delete MD;
#ifdef DEBUG
			cerr << "MD deleted" << endl;
#endif

		//stats
		fil->prepStats();
		if (IS->getCurFileN() == 0) {
			fil->printStats(cerr, mainFile, outFile, true);
		} else {
			cerr<<fil->shortStats(""); shortStats = true;
		}

		totReadsRead += fil->totalAccepts();
		
		shrtLog += fil->shortStats( mainFileShort);
//		delete IS;
		//write log file
		if (uniqueFas.size() > 1){//only print sub log if neccessary
			ofstream log;
			string logF = cmdArgs["-log"] + string("0") + itos(i);
			log.open (logF.c_str() ,ios_base::out);
			fil->printStats(log,mainFile,outFile,true);
			log.close();
		}
		mainFil->addStats(fil,idx[i]);
#ifdef DEBUG
		cerr << "Delete tmp filter" << endl;
#endif
		//and cleanup
		//
		//fil;
	}
#ifdef DEBUG
	cerr << "Prep final logging" << endl;
#endif
//write log files
	if (uniqueFas.size() > 1){
		mainFile = "several";
	}
	
	ofstream log; string deLog("");
	string logF = cmdArgs["-log"], logFA = cmdArgs["-log"].substr(0, cmdArgs["-log"].length()-3) + "add.log";

	//different logfile for SEED extension
	if (mainFil->doOptimalClusterSeq()){
		//finish up dereplication file (creating pseudo seeds with counts)
		ucl->finishMAPfile();
		if (cmdArgs["-ucAdditionalCounts"] != ""){
			ucl->set2UC();
			ucl->finishUCfile(mainFil, cmdArgs["-ucAdditionalCounts"], true);//with smplHead (.mid)
			ucl->finishUCfile(mainFil, cmdArgs["-ucAdditionalCounts1"], false);//without smplHead (.rest)
		}
		if (cmdArgs["-ucAdditionalCounts_refclust"] != ""){
			//reference based clustering has some high qual seqs (no replacement with reads..)
			shared_ptr<InputStreamer> FALL = make_shared<InputStreamer>(true, mainFil->getuserReqFastqVer());
			//this reads in the SLV fna's & creates matrix entries for these
			FALL->setupFna(cmdArgs["-OTU_fallback_refclust"]);
			ucl->setRefMode();
			ucl->addDefSeeds(FALL, mainFil);
			ucl->set2UC();
			//mapping from ref OTU clustering
			ucl->finishUCfile(mainFil, cmdArgs["-optimalRead2Cluster_ref"], false);
			//mid / rest mappings
			ucl->finishUCfile(mainFil, cmdArgs["-ucAdditionalCounts_refclust"], true);//with smplHead (.rest)
			ucl->finishUCfile(mainFil, cmdArgs["-ucAdditionalCounts_refclust1"], false);//without smplHead (.rest)

		}
		if (cmdArgs["-log"] != "nolog") {
			log.open(logF.c_str(), ios_base::out);
			ucl->printStats(cerr);
			ucl->printStats(log);
			log.close();
		}
		
		//everything done on DNA? Then write & delete
		if (cmdArgs["-otu_matrix"] != "") {
			ucl->writeOTUmatrix(cmdArgs["-otu_matrix"], mainFil);
		}
		//needs to be written after OTU matrix (renaming scheme)
		shared_ptr<MultiDNA> MD = make_shared<MultiDNA>(mainFil, cmdArgs, ios::out, RDSset);
		mainFil->setMultiDNA(MD);
		ucl->writeNewSeeds(MD, mainFil,false);
		//delete MD;
		//new fastas also need to be written..
		MD.reset(new MultiDNA(mainFil, cmdArgs, ios::out, RDSset, ".ref", 1));//force fna output
		mainFil->setMultiDNA( MD );
		ucl->writeNewSeeds(MD, mainFil,true,true);
		//delete MD;


		return;
	} else if (mainFil->doDereplicate()) {
#ifdef DEBUG
		cerr << "write Dereplicated DNA" << endl;
#endif
		deLog = Dere->writeDereplDNA(mainFil);
#ifdef DEBUG
		cerr << "done write Dere" << endl;
#endif
	}
#ifdef DEBUG
	cerr << "Logging almost finished" << endl;
#endif

	if (cmdArgs["-log"] == "nolog") {
//		delete Dere;
		return;
	}
#ifdef DEBUG
	cerr << "DereLog start" << endl;
#endif
	if (mainFil->doDereplicate()) {
		string dereLog = logF.substr(0,logF.length()-3) + "dere";
		Dere->writeLog(dereLog, deLog);
//		delete Dere;
	}
	
#ifdef DEBUG
	cerr << "DereLog end" << endl;
#endif
	if (shortStats) {
		mainFil->printStats(std::cerr, mainFile, outFile, true);
	}
#ifdef DEBUG
	cerr << "other logs start" << endl;
#endif
	//per sample success rate
	string logPS = logF.substr(0, logF.length() - 3) + "acceptsPerSample.log";
	log.open(logPS.c_str(), ios_base::out);
	mainFil->SmplSpecStats(log);
	log.close();
	log.open(logF.c_str(), ios_base::out);
	mainFil->printStats(log,mainFile,outFile,true);
	log.close();

	string logFs = logF.substr(0, logF.length() - 3) + "acceptsPerFile.log";
	log.open(logFs.c_str(), ios_base::out);
	log << shrtLog;
	log.close();
	
	string logFGC = logF.substr(0, logF.length() - 3) + "GC.txt";
	log.open(logFGC.c_str(), ios_base::out);
	mainFil->printGC(log, mainFil->isPaired());
	log.close();
#ifdef DEBUG
	cerr << "other logs end" << endl;
#endif


//for additional files
	if (mainFil->secondaryOutput()){
		log.open (logFA.c_str() ,ios_base::out);
		mainFil->printStats(log,mainFile,outFile,false);
		log.close();
	}


	//length histogram
	logF = cmdArgs["-length_hist"];
	log.open (logF.c_str() ,ios_base::out);
	mainFil->printHisto(log,0);
	log.close();
	//quality histogram
	logF = cmdArgs["-qual_hist"];
	log.open (logF.c_str() ,ios_base::out);
	mainFil->printHisto(log,1);
	log.close();
	mainFil->close_outFiles_demulti();
#ifdef DEBUG
	cerr << "separateByFile finished" << endl;
#endif

}

void rewriteNumbers(OptContainer& cmdArgs){
	//no renumbering asked for
	if (!(cmdArgs.find("-number")  != cmdArgs.end() && cmdArgs["-number"]=="T")){
		return;
	}
	string prefix="";
	if (cmdArgs.find("-prefix")  != cmdArgs.end()){
		prefix = cmdArgs["-prefix"];
	}
	//read fasta & write with new headers
	int cnt=0;

	string line;
	//ofstream qualOut,fnaOut;
	ifstream fna;
	ofstream ofna;

	//        rerwite input fasta file
	string tname="",tseq="";
	fna.open(cmdArgs["-i_fna"].c_str(),ios::in);
	ofna.open(cmdArgs["-o_fna"].c_str(),ios::out);
	while (getline(fna,line,'\n')){

		if (line[0]=='$'){ //$ marks comment
			continue;
		}
		if(line[0] == '>'){ //fasta description
			if (cnt!=0){
				tname = ">"+prefix+itos(cnt)+"\n";
				ofna << tname << tseq;
			}
			cnt++;tseq="";
			continue;
		}
		tseq += line+"\n";
	}
	tname = ">"+prefix+itos(cnt)+"\n";
	ofna << tname << tseq;
	ofna.close(); fna.close();

	exit(0);
}

void Announce_sdm(){
	cerr << endl << "This is sdm (simple demultiplexer) " << sdm_version << " " << sdm_status << "." << endl << endl;
}
void help_head(){
	cout <<"------------------------------\nThis is sdm version "<<sdm_version <<" "<< sdm_status <<" help print\n------------------------------\n"<<endl;
}
void general_help(){
	help_head();
	cout<<"sdm (simple demultiplexer) is a fast, memory efficient program to demultiplex fasta and fastq files or simply do quality filterings on these.\n";
#ifdef _gzipread
	cout<<"Compiled with gzip support\n";
#else
	cout << "No gzip support compiled\n";
#endif
#ifdef _THREADED
	cout<<"Multithreading not supported"
#else
	cout << "The compiled version does not support multithreading\n";
#endif
	cout<<"Select further help topics by typing:\nsdm -help_options : print help on configuring options files\nsdm -help_commands : help on command arguments for sdm\nsdm -help_map : map files and the keywords for barcodes etc.\n------------------------------\n";
	cout<<"Author: falk.hildebrand@gmail.com"<<endl;

}
void printCmdsHelp(){
	help_head();
	string def_sep = DEFAULT_BarcodeNameSep;
	cout << "Usage:\n./sdm\n  -i_path <path to several fastq / fasta files>\n------OR------\n -i <input sequence file, will autodetect fna/fastq>\n------OR------\n -i_fastq <fastQ file>\n------OR------\n -i_fna <your fasta input file> (required)\n -i_qual <corresponding quality file> (required, unless quality file is \"xx1.qual\" and fasta is \"xx1.yy\")\n\n -map <mapping file in Qiime format> (optional)\n -o_fna <file to write output fasta> (optional)\n -o_qual <file to write corresponding quality values> (optional)\n -o_fastq <fastQ output file (overrides -o_qual & -o_fna)\n";
	cout << " -options <sdm option file>(optional)\n -log <file to save demultiplex log in>(optional). Set to \"nolog\" to deactivate alltogether.\n \n-sample_sep \"X\" string X is used to delimit samples and ID (optional, default:\"" << def_sep << "\")\n -paired 1/2/3 (input is paired end sequenced(2), assumes two input files delimited by \',\'. 1=singleton (default); 3=paired end (R1,R3) + one file with MID (R2))\n";
	cout << " -o_demultiplex [path] write input into single, demultiplexed files\n";
	cout << " -onlyPair [1/2] consider only read pair 1 or 2. Useful when streamlining inputs (LotuS) or considering double barcoding.\n -i_MID_fastq fastq file with only MID sequences; if paired reads are supplied with -i_fna/-i_fastq and the MID identifier via -i_MID_fastq, paired has to be set to 2. If e.g. merged reads are supplied + mids, paired has to be set to 1.\n";
	cout << " -oneLineFastaFormat [0/1] write Fasta and Quality file sequence string in one line, opposed to default 80 characters per line.\n -o_dereplicate <output fasta file> of dereplicated DNA reads (with size in header)\n -dere_size_fmt [0/1] either (0) usearch format \"size=X;\" or (1) \"_X\"\n -min_derep_copies only print seq if at least X copies present. Can be complex terms like \"10:1,3:3\" -> meaning at least 10x in 1 sample or 3x in 3 different samples.\n";
	cout << " -SyncReadPairs [T/F] sdm can check, if read pairs occur in the same (correct) order in the input files, and correct this in case not (T).\n";
	cout << " -maxReadsPerOutput number of filtered reads in output files. If more reads, a new file is created. Only works with -o_fna\n -mergedPairs <1/0> 1: paired sequences were merged externally, important for assumption that read quality is detoriating.\n -OTU_fallback <file>: Fallback fasta sequences for OTU's, only used in SEED extension mode\n";
	cout << " -i_qual_offset [0-64] fastq offset for quality values. Set this to \'0\' or \'auto\' if you are unsure which fastq version is being used (default: read from sdm option file)\n -o_qual_offset [0-64] set quality offset for fastq outfile. Default: 33\n";
	cout << " -ignore_IO_errors [0/1]: 1=Errors in fastq reads are ignored, with sdm trying to sync reads pairs after corrupted single reads (default: 0)\n";
	//-binomialFilterBothPairs [1/0]
	//-count_chimeras [T/F]
	// ucAdditionalCounts_refclust -OTU_fallback_refclust -optimalRead2Cluster_ref
	cout<<"\nMinimal Example:\n./sdm -i test.fna -map mapping.txt (assuming quality file is \"test.qual\")\n";
	// further options (undocumented) : 
	//-length_hist   -qual_hist
	//-suppressOutput[0/1]
	//
}
void printOptionHelp(){
	help_head();
	cout<<"The option file, specified via the \"-options\" argument, provides more specific control over filtering, barcode handling, and sequencing technologies, among others. A reference option file is printed below.\n\n";
	string helpOptionFile="";
	/*helpOptionFile += "minSeqLength - minimal accepted Sequence Length\nmaxSeqLength - maximal Length of Sequence\nminAvgQuality - minimal average Quality\nmaxAmbiguousNT - max number of Ambigous nt's in sequence\nQualWindowThreshhold - Q threshold where seq is rejected\nQualWindowWidth - average quality in this windows is used for QualWindowThreshhold\n";
	helpOptionFile += string("TrimWindowThreshhold - Q value below which sequence is 3' trimmed\nTrimWindowWidth - window size used for TrimWindowThreshhold\nmaxBarcodeErrs - max accepted barcode errors\nmaxPrimerErrs - max accepted Primer errors\nkeepBarcodeSeq - leave Barcode Seq on read? (0/1)\n");
	helpOptionFile += "keepPrimerSeq - keep Primer attached to seq? (0/1)\nmaxHomonucleotide - sequences with a homonucleotide run longer will be rejected\nmaxAccumulatedError - if P is surpassed, sequence is trimmed at that point\nTechnicalAdapter - sequence of the technical adapter (will be removed, if found 5')\nPEheaderPairFmt - ?\nTrimStartNTs - trim X nucleotides from the start of the sequence ";
	helpOptionFile += "fastqVersion - 1 = ";*/

	helpOptionFile += "#--- Example ---\n#copy into new file\n#sequence length refers to sequence length AFTER removal of Primers, Barcodes and trimming. this ensures that downstream analyis tools will have appropiate sequence information\nminSeqLength	250\nmaxSeqLength	1000\nminAvgQuality	25\n\n";
	helpOptionFile += "#Ambiguous bases in Sequence - uclust only supports 0 ambiguous nucleotides\nmaxAmbiguousNT	0\n\n#Homonucleotide Runs.. this should normally be filtered by sequencer software\nmaxHomonucleotide	8\n\n";
	helpOptionFile += "#Filter whole sequence if one window of quality scores is below average\nQualWindowWidth	50\nQualWindowThreshhold	25\n\n#Trim the end of a sequence if a window falls below quality threshhold. Useful for removing low qulaity trailing ends of sequence\n\nTrimWindowWidth	20\nTrimWindowThreshhold	25\n\n#Max number of accumulated P for a mismatch. After this length, the rest of the sequence will be deleted. Complimentary to TrimWindowThreshhold. (-1) deactivates this option.\nmaxAccumulatedError	1\n\n";
	helpOptionFile += "#Barcode Errors - currently this can only be 0; \nmaxBarcodeErrs	0\nmaxPrimerErrs	0\n\n#keep Barcode / Primer Sequence in the output fasta file - in a normal 16S analysis this should be deactivated (0) for Barcode and de-activated (0) for primer\nkeepBarcodeSeq	0\nkeepPrimerSeq	0\n\n";
	helpOptionFile += "#set fastqVersion to 1 if you use Sanger, Illumina 1.8+ or NCBI SRA files. Set fastqVersion to 2, if you use Illumina 1.3+ - 1.7+ or Solexa fastq files.\n\nfastqVersion	1\n\n#if one or more files have a technical adapter still included (e.g. TCAG 454) this can be removed by setting this option\n\nTechnicalAdapter	TCAG\n\n#delete X NTs (e.g. if the first 5 bases are known to have strange biases)\n\nTrimStartNTs	0\n";
	helpOptionFile += "#truncate total Sequence length to X (length after Barcode, Adapter and Primer removals)\nTruncateSequenceLength	200\n";
	helpOptionFile += "#correct PE header format (1/2) this is to accomodate the illumina miSeq paired end annotations 2=\"@XXX 1:0:4\" instead of 1=\"@XXX/1\". Note that the format will be automatically detected\nPEheaderPairFmt	1\n\n#sets if sequences without match to reverse primer will be accepted (T=reject ; F=accept all); default=F\nRejectSeqWithoutRevPrim	F\n";
	helpOptionFile += "#sets if sequences without a forward (LinkerPrimerSequence) primer will be accepted (T=reject ; F=accept all); default=T\nRejectSeqWithoutFwdPrim	T\n\n";
	helpOptionFile += "#checks if the whole amplicon was reverse-transcribed sequenced (Default = F)\nCheckForReversedSeqs	F\n\n";
	helpOptionFile += "#this option should be \"T\" if your amplicons are possibly shorter than a read in a paired end sequencing run (e.g. amplicon of 300 in 250x2 miSeq is \"T\")\nAmpliconShortPE	T\n\n";
	//CheckForMixedPairs CheckForReversedSeqs
	cout<<helpOptionFile<<endl;
}
void printMapHelp(){
	help_head();
	cout<<"The mapping file, specified via the \"-map\" argument, contains all information neccessary to demultiplex a sequencer fasta or fastq output file. If left out, the sequences are only quality checked.\n";
	cout<<"The Mapping file can contain comments, specified by \"#\" at the start of the comment line. Similarly, the header (required) has to start with \"#\"\n. Processed header fields are:\n";
	cout << "SampleID - Sample Identifier, has to be unique for each Barcode\n";
	cout << "CombineSamples - combines Sample Identifiers into one sample\n";
	cout << "BarcodeSequence - The Barcode (MID) tag assigned to each sample. Can contain IUPAC redundant nucleotides\n";
	cout << "Barcode2ndPair - in case of dual indexed reads, use this column to specify the BC on the 2nd read pair\n";
	cout<<"ForwardPrimer (previously LinkerPrimerSequence) - Sequence used for 16S amplification, usually (unless paired end mode) is after the Barcode\n";
	cout<< "ReversePrimer - Reverse Primer Sequence (IUPAC code).\n"; 
	cout<< "fastqFile - if used in -i_path mode, gives relative location of fastq file, such that [-i_path][fastqFile] gives the absolute path to fastq file.\n";
	cout<<"fnaFile - see fastqFile. However, fasta formated file instead of fastq format.\n";
	cout<< "qualFile - see fastqFile. However, quality file corresponding to fasta file instead of fastq format.\n";
	cout<< "MIDfqFile - fastq file containing ONLY the MID sequence, corresponding fastq input file(s)\n";
	cout<< "SampleIDinHead - ID in header of fasta/fastq file, that identifies Sample - replaces Barcode (MID) scanning.\n";
	cout << "derepMin - this column can contain a single number X, only accepting reads in specific sample that have more than X dereplicated copies. Useful if a subset of samples was e.g. sampled 100x as deep as other samples.\nCombined Sample WILL NOT be taken into account.\n";
	cout<<"\n";

}
void printVersion(){
	cout << "sdm " << sdm_version << " " << sdm_status << endl;
}


//bool readCmdArgs(int argc, char* argv[],map<char*, char*, lstr>& cmdArgs);
