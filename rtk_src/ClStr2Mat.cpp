#include "ClStr2Mat.h"


ClStr2Mat::ClStr2Mat(const string inF, const string outF,
	const string mapF, const string basePX, bool covCalc):
	GAs(0), CCH(NULL),smplLoc(0), baseP(0), smplN(0), curr(-1) {
	ifstream incl;
	//set up baseP
	stringstream ss(basePX); string segments;
	while (getline(ss, segments,',') ) { baseP.push_back(segments); }


	ofstream matO, geneNames;
	matO.open(outF + ".mat");
	if (!matO) {
 #ifdef notRpackage
cerr << "Couldn't open matrix output " << outF + ".mat" << endl;
exit(57);
#endif
}
	geneNames.open(outF + ".genes2rows.txt");
	if (!geneNames) {
 #ifdef notRpackage
cerr << "Couldn't open report file " << outF + ".genes2rows.txt" << endl;
exit(56);
#endif
}
	incl.open(inF.c_str());
	if (!incl) {
 #ifdef notRpackage
cerr << "Couldn't open clustering file " << inF << endl;
exit(55);
#endif
}

	//read map(s) and check that
	stringstream ss2(mapF);
	while (getline(ss2, segments, ',')) {
		read_map(segments, covCalc);
	}
	if ( baseP.size() > curr+1) {
	 #ifdef notRpackage
		cerr << "more maps than basePs\n"; 
		exit(72);
	#endif
	}


	//smplnames in out matrix
	vector<string> SmplNmsVec(smplN, "");
	for (auto it = smpls.begin(); it != smpls.end(); it++) {
		vector<int> XX = (*it).second;
		SmplNmsVec[ XX[XX.size()-1] ] = (*it).first;
	}
	matO << "Genes";
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		matO << "\t" << SmplNmsVec[i];
	}
	matO << endl;

	geneNames << "#GID\tCluster\tRepSeq"<<endl;
	CCH = new ContigCrossHit((int)smplN);
	CCH->setSmplNms(SmplNmsVec);
	//SparseMatrix * mat = new SparseMatrix();
	string line;
	int CLidx = 1;
	const string sampleSeq = "__";
	vector<smat_fl> repVec(smplN, 0.f), SmplSum(smplN, 0.f);
	bool repFound = false;
	while (getline(incl, line)) {
		if (line.substr(0, 1) == ">") {//new cluster, add line to matrix
			//mat->newRow();//itos(CLidx)
			//write the vector out
			if (CLidx >= 2) {
				//report part is here
				printVec(matO, repVec, itos(CLidx));
				if (!repFound) {
					geneNames << "\t?\n";

				 #ifdef notRpackage
				cerr << "No RepSeq found for " << line <<" -1"<< endl;
				#endif
				}
			}
			CLidx++; repFound = false;
			geneNames << itos(CLidx) << "\t" << line;
			repVec.clear(); repVec.resize(smplN);
			continue;
		}
		//1 get gene, sample (deparse)
		size_t pos = line.find("nt, >");
		size_t pos2 = line.find("...", pos + 4);
		string gene = line.substr(pos + 5, pos2 -pos -5 );

		if (!repFound && line.back() == '*') {//report representative gene
			geneNames << "\t"<<gene<< endl;
			repFound = true;
		}

		bool geneInAssembl(true);
		pos = gene.find(sampleSeq);
		string sample;
		SmplOccurITmult smNum;
		pos2 = gene.find("_L", pos + 3);
		if (pos != string::npos && pos2 != string::npos) { //has the characteristic "__" sample separator
			sample = gene.substr(0, pos);
			smNum = smpls.find(sample);
			if (smNum == smpls.end()) {
#ifdef notRpackage
				cerr << "incorrect sample name: " << sample << endl << gene << endl;
				exit(55);
#endif
			}
			//2 get abundance
			//Contig + Sample Info will allow to create contig linkage between samples (CCH)
			//int contig = atoi(gene.substr(pos + 3, pos2-pos-3).c_str());
			//can be several samples (combined assembly)
			vector<int> smplLocs = (*smNum).second;
			for (size_t jj = 0; jj < smplLocs.size(); jj++) { //the loop takes account for multiple samples being grouped together in map
															  //CCH->addHit(smplLocs[jj], contig);
				smat_fl abundance = GAs[smplLocs[jj]]->getAbundance(gene);
				//3 add to matrix / output vector
				repVec[smplLocs[jj]] += abundance;
				SmplSum[smplLocs[jj]] += abundance;
				//mat->addCount(sample, CLidx, abundance);
			}
		} else {
			geneInAssembl = false;
		}
	}
	incl.close();
	matO.close(); geneNames.close();
	matO.open(outF + ".mat.sum");
	//print sample sum
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		//cout << SmplNmsVec[i] << "\t" << SmplSum[i] << endl;
		matO << SmplNmsVec[i] << "\t" << SmplSum[i] << endl;
	}
	matO.close();
}
ClStr2Mat::~ClStr2Mat() {
	for (size_t i = 0; i < GAs.size(); i++) { delete GAs[i]; }
	delete CCH;
}
void ClStr2Mat::printVec(ofstream& of,vector<smat_fl>& pr,const string&rowN) {
	of << rowN;
	for (size_t i = 0; i < pr.size(); i++) {
		of << "\t" << pr[i];
	}
	of << "\n";
}
void ClStr2Mat::read_map(const string mapF,bool calcCoverage) {
	ifstream in;
	curr++;//keep track of different maps and inPaths
	uint preMapSize((int)smplLoc.size());
	if (curr > baseP.size()) {
 #ifdef notRpackage
cerr << "more maps than basePs\n";
exit(72);
#endif
}
	in.open(mapF.c_str());
	if (!in) {
 #ifdef notRpackage
cerr << "Couldn't open mapping file " << mapF << endl;
exit(56);
#endif
 }
	#ifdef notRpackage
	cout << "Reading map " << mapF << " on path " << baseP[curr] << endl;
	#endif
	SmplOccurMult CntAssGrps;
	string line; int cnt(-1); int assGrpN(-1); 
	int artiCntAssGrps(0); int skSmplCol(-1);

	while (getline(in, line)) {
		cnt ++; int sbcnt(-1);
		stringstream ss (line); string segments;
		if (line.substr(0, 1) == "#") { //read column position of smpl id, path & assGrps
			if (cnt > 0) { continue; }

			while (getline(ss, segments, '\t')) {
				sbcnt++;
				if (sbcnt==0 && segments != "#SmplID") {
					 #ifdef notRpackage
					cerr << "Map has to start with tag \"#SmplID\"\n";exit(83);
					#endif
				}
				if (sbcnt == 1 && segments != "Path") {
					#ifdef notRpackage
					cerr << "Map has to have tag \"Path\" as second entry\n";
					exit(83);
					#endif
				}
				if (segments == "AssmblGrps") {
					assGrpN = sbcnt; 
					 #ifdef notRpackage
					cout << "Found Assembly groups in map\n";
					#endif
				}
				if (segments == "ExcludeAssembly") {
					skSmplCol = sbcnt; 
					 #ifdef notRpackage
					cout << "Samples can be excluded from assembly\n";
					#endif
				}
			}
			continue;
		}
		vector<string> curLine(0);
		while (getline(ss, segments, '\t')) {
			curLine.push_back(segments);
		}
		if (skSmplCol>-1 && curLine[skSmplCol] == "1") { continue; }


		string smpID = curLine[0];
		if (smpls.find(smpID) != smpls.end()) {
			#ifdef notRpackage
			cerr << "Double sample ID: " << smpID << endl;
			exit(12);
			#endif
		}


		//getline(ss, segments, '\t');
		string locality = curLine[1];

		string assGrp ("");
		if (assGrpN != -1) {
			//handles assembly groups from here
			assGrp  = curLine[assGrpN];
		} else {//simulate CntAssGrps
			assGrp = itos(artiCntAssGrps);
			artiCntAssGrps++;
		}
		if (CntAssGrps.find(assGrp) != CntAssGrps.end()) {
			CntAssGrps[assGrp].push_back( (int)smplLoc.size());
		} else {
			CntAssGrps[assGrp] = vector<int>(1,(int)smplLoc.size());
		}

		if (CntAssGrps[assGrp].size() > 1) {
			string nsmpID = smpID + "M" + std::to_string(CntAssGrps[assGrp].size());
			smpls[nsmpID] = CntAssGrps[assGrp];//(int)smplLoc.size();
		} else {
			smpls[smpID] = vector<int>(1,(int)smplLoc.size());

		}


		smplLoc.push_back(locality);

	}
	in.close();
	smplN = smplLoc.size();
	//read the gene abundances sample-wise in
	for (uint i = preMapSize; i < smplN; i++) {
		string pa2ab = path2counts;
		if (calcCoverage) { pa2ab = path2abundance; }
		GAs.push_back(new GeneAbundance(baseP[curr] + "/" + smplLoc[i], pa2ab));
		#ifdef notRpackage
		cerr << baseP[curr] + "/" + smplLoc[i] << endl;
		#endif
	}
}


///////////////////////////////////////////////////
void ContigCrossHit::addHit(int Smpl, int Ctg) {

}


///////////////////////////////////////////////////

GeneAbundance::GeneAbundance(const string path, const string abunF):
	isPsAss(false){
	ifstream in;
	//first test if this is a pseudoassembly
	in.open((path + pseudoAssMarker).c_str());
	if (in) {
		in.close();
		isPsAss = true;
		return;
	}
	in.close();
	//not? then read abundances
	string newS = path + abunF;
	in.open(newS.c_str());
	if (!in) {
	 #ifdef notRpackage
	cerr << "Couldn't open gene abundance file " << newS << endl;
	exit(36);
	#endif
}
	string line;
	while (getline(in, line)) {
		size_t pos = line.find("\t");
		string gene = line.substr(0, pos);
		GeneAbu[gene] = (smat_fl)atof(line.substr(pos + 1).c_str());
	}
	in.close();
}
smat_fl GeneAbundance::getAbundance(const string x) {
	if (isPsAss) {
		return (smat_fl) 1.f;//return one read count
	}
	SmplAbunIT fnd = GeneAbu.find(x);
	if (fnd == GeneAbu.end()) {
		return (smat_fl) 0;

 #ifdef notRpackage
cerr << "Can't find " << x << endl;
exit(33);
#endif
	}
	return (*fnd).second;
}
