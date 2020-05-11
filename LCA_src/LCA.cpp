// LCA.cpp 

#include "libload.h"
#include "RefTax.h"
#include "LCAimpl.h"
#include "Matrix.h"

const char* LCA_ver = "0.18		alpha";

void helpMsg() {
	cout << "LCA requires at least 3 arguments (-i, -r, -o)\n For more help and options, use \"./LCA -h\"\n";
#ifdef _gzipread
	cout << "Compiled with gzip support\n";
#else
	cout << "No gzip support compiled\n";
#endif
}
void welcomeMsg() {
	cout << "Least common ancestor (LCA) assignments ver " << LCA_ver << endl;
}
int main(int argc, char* argv[])
{
	
	welcomeMsg();
	//measure execution time
	clock_t tStart = clock();

	options* OPT = new options(argc, argv, __default_depth);


	if (argc < 2) {
		helpMsg();
		exit(3);
	}
	/*string RefTaxFile = argv[2];
	string blastres = argv[1];
	string outF = argv[3];*/
	int numThr = OPT->numThr;
	bool highLvl(OPT->calcHighMats);

	//and prep potential higher level mat (vector based, no samples are separated)
	Matrix* mat = new Matrix(OPT->taxDepth, OPT->Taxlvls, OPT->hitRD);
	//parrallel prepare
	int workerNum(numThr - 1);
#ifdef parallel
	std::future<TaxObj*> *parvec = new std::future<TaxObj*>[numThr];
#else
	vector<TaxObj*> sCore(numThr, NULL);
#endif
	//list<TaxObj*> assigns(0);
	unordered_map<string, TaxObj*> assign;

	//ini output stream
	string IDname = "OTU";
	if (OPT->isReads) { IDname = "Reads"; }
	ofstream O(OPT->outF.c_str());
	O << IDname << "\t" << OPT->TaxLvl2string();
	if (OPT->hitRD) { O << "\tHit2DB"; }
	O << endl;// Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tOTU\n";


	//some stats
	int dblSbj(0); int replSbjTax(0);
	//only 1 DB use? Doesn't need to keep hits in mem..
	bool multiDBuse = true; 
	if (OPT->refDBs.size() == 1) { multiDBuse = false; }

	for (int xi = 0; xi < (int) OPT->refDBs.size(); xi++) {
		//read in tax DB
		RefTax* RT = new RefTax( OPT->refDBs[xi], OPT->taxDepth, OPT->nativeSlVdb, OPT->checkTaxoUnkw);
#ifdef DEBUG
		cerr << "Ref tax read\n";
#endif // DEBUG

		RT->setTaxLvls(OPT->Taxlvls);
		//ini & parse blast
		BlastReader* BR = new BlastReader(OPT->blFiles[xi], OPT->input_format);

		//some flags for the parallel execution
		bool burninDone(false);//ini all Cores
		bool allRead(false);//blast file parsed
		int ti(0); int ti_end(-1);

		//loop over all blast results
		while (1) {
			if (allRead && ti_end == ti) {//no more blast results
				cout << "Done Blast File reading\n"; break;
			}
			if (burninDone) {
				//assigns.push_back(parvec[ti].get());
#ifdef parallel
				TaxObj* tmp = parvec[ti].get();
#else
				TaxObj* tmp = sCore[ti];
#endif
				if (!multiDBuse){	//and write the tax out
					O << (tmp)->Subj << "\t" << (tmp)->getWriteString() << endl;
					if (highLvl) {
						mat->add(tmp);
					}
					delete tmp; 
				} else {
					//this part refers to previously determined tax (needs > 1 ref DB)
					auto tf = assign.find(tmp->Subj);
					if (tf == assign.end()) {
						assign[tmp->Subj] = tmp;
					}
					else {
						dblSbj++;
						if (assign[tmp->Subj]->evalAcpyTax(tmp)) {
							replSbjTax++;
						}
						delete tmp;
					}
				}
			}


			//don't forget to destroy this object -> yes, done in LCA algo (better for async operation)
			list<BlastRes*> tmpB = BR->getResBatch();

			//check if this position LCA is finished
			if (!allRead) {
				if (tmpB.size() == 0) {//break routine,write final Results
					allRead = true;
					ti_end = ti;
					//& fetch remaining worker jobs
				} else {
#ifdef parallel
					parvec[ti] = async(std::launch::async, LCA, tmpB, RT, OPT);
#else
					sCore[ti] = LCA(tmpB, RT, OPT);
#endif
				}
			}

			ti++;
			if (ti > workerNum) { //multi core control
				//TODO this is wrong ? if (allRead) { break; } 
				ti = 0; burninDone = true;
			}
		}
		delete RT; delete BR;
	}

//work through remaining hits from multi DBs.. should not go through here in single DB mode
	for (auto it = assign.begin(); it != assign.end(); it++) {
		O << (it->second)->Subj << "\t" << (it->second)->getWriteString() << endl;
		if (highLvl) {
			mat->add(it->second);
		}
		delete it->second; it->second = NULL;
		//into file now
	}
	O.close();
	//clean up 

	if (highLvl) {
		mat->writeAllLevels(OPT->outF);
	}
	delete mat; delete OPT;

	if (dblSbj>0) {
		cout << "Found " << dblSbj << " double subject sequences in "<< OPT->refDBs.size()<< ", reassigned " << replSbjTax << " of these." << endl;
	}
	printf("LCA finished. Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);


}