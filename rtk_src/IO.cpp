#include "IO.h"
std::mutex rarefyMutex;


void lineCntOut(const string inF, const string outF, const string arg4){
    ifstream in(inF.c_str());
    ofstream out(outF.c_str(), ios::out);
    if (!in){
#ifdef notRpackage
        cerr << "Can't open infile " << inF << endl; std::exit(99);
#endif
    }
    if (!out){
#ifdef notRpackage
        cerr << "Can't open outfile " << outF << endl; std::exit(99);
#endif
    }
    //read file that contains nicely ordered all indexes of lines to be extracted
    string line;
    vector<uint> srtTar;
    ifstream idxS(arg4.c_str());
    if (!idxS){
#ifdef notRpackage
        cerr << "Can't open outfile " << arg4 << endl; std::exit(99);
#endif
    }
    while (getline(idxS, line, '\n')) {
        if (line[0] == '>'){
            line.erase(0,1);
        }
        srtTar.push_back(stoi(line));
    }
    idxS.close();
    //sort ascending
    sort(srtTar.begin(), srtTar.end());

    //sort through new file
    if (!out){
#ifdef notRpackage
        cerr << "Can't open outfile " << outF << endl; std::exit(99);
#endif
    }
    uint cnt(1); uint j(0);
    while (getline(in, line, '\n')) {
        if (cnt == srtTar[j]){
            out << line << endl;
            uint cur = srtTar[j];
            while (srtTar[j] == cur){ j++; }
            if (j == srtTar.size()){ break; }
        }
        cnt++;
    }

    in.close(); out.close();
    if (j != srtTar.size()){

#ifdef notRpackage
        cerr << "Missed " << (srtTar.size() - j) << " entries." << endl;
#endif
    }
}
//****************************  smplVec::smplVec ***********
smplVec::smplVec(const vector<mat_fl>& vec, const int nt) :IDs(0),totSum(0),
    num_threads(nt), richness(-1), Shannon(-1.f){
        double cumSum(0.f);
        for (uint i = 0; i < vec.size(); i++){
            cumSum += vec[i];
        }
        if (verbose){
#ifdef notRpackage
            cerr << (long)cumSum << " allocating ";
#endif
        }
        //arr = (int*) malloc((int)cumSum * sizeof(int));
        //arr = new unsigned short[(int)cumSum];
        arr.resize((long)cumSum);
        if (verbose){
#ifdef notRpackage
            cerr << "memory";
#endif
        }
        totSum = cumSum;
        long k(0); uint posInVec(-1);
        //numFeatures = 0;
        for (size_t i = 0; i< vec.size(); i++){

            long maxG = (long)vec[i];
            IDs.push_back( std::to_string(i));

            posInVec++;
            if (maxG == 0){ continue; }//not really a feature, doesnt need ot be counted as cat

            maxG += k; //bring up to current level
            for (; k<maxG; k++){
                arr[k] = posInVec;
            }
            //numFeatures++;
            // some simple numeric id for refernce in chao2 abundance calculations, so it behaves as the swap mode
        }
        posInVec++;
        numFeatures = posInVec;
        if (verbose){
#ifdef notRpackage
            cerr << "..\n";
#endif
        }
    }
smplVec::smplVec(const string inF, const int nt) :IDs(0),totSum(0), num_threads(nt),
    richness(-1),Shannon(-1.f) {

        vector<double> vec;
        ifstream in(inF.c_str());
        string line; double cumSum(0.f);
        while(getline(in,line,'\n')) {
            string ID; float num;
            stringstream ss(line);
            ss>>ID; ss>>num;
            cumSum += num;
            vec.push_back(num); IDs.push_back(ID);
        }
        in.close();

        if (verbose){
#ifdef notRpackage
            cerr<<(long)cumSum<<" allocating ";
#endif
        }
        arr.resize((long)cumSum);
        if (verbose){
#ifdef notRpackage
            cerr<<"memory";
#endif
        }
        totSum = cumSum;
        long k(0); uint posInVec(0);
        for (size_t i = 0; i< vec.size(); i++){

            long maxG = (long)vec[i];
            maxG += k;
            if (maxG == 0){ continue; }//not really a feature, doesnt need ot be counted as cat
            for (; k<maxG; k++){
                arr[k] = posInVec;
            }
            posInVec++;
        }
        numFeatures = posInVec;
        if (verbose){
#ifdef notRpackage
            cerr<<"..\n";
#endif
        }
    }


void smplVec::rarefy(vector<double> depts, string ofile, int rep,
        DivEsts* divs, std::vector<vector<rare_map>> & RareSample,
        vector<string>& retCntsSampleName, string& skippedSample,
        vector<vector<vector<uint>>>* abundInRow, vector<vector<vector<uint>>>* occuencesInRow,
        int writes,bool write, bool fillret){
    bool doShuffle = true;
    long curIdx = 0;
    long dep;
    // resize divvs
    divs->richness.resize(depts.size());
    divs->shannon.resize(depts.size());
    divs->simpson.resize(depts.size());
    divs->invsimpson.resize(depts.size());
    divs->chao1.resize(depts.size());
    divs->eve.resize(depts.size());
    
    for(uint i = 0; i < depts.size(); i++){
        dep = depts[i];

        if (dep>totSum){
            skippedSample = divs->SampleName;
            if (verbose){cout<<"skipped sample, because rowSums < depth \n";}
            return;
        }
        //long curIdx=(long)totSum+1;
        for (int curRep=0;curRep<rep;curRep++){
            if(curIdx+dep >= (long) totSum || doShuffle == true){
                shuffle_singl();	
                curIdx = 0;
                doShuffle = false;
            }


            //count up
            vector<uint> cnts(numFeatures, 0);
            for (long i=(0+curIdx);i<(dep+curIdx);i++){
                cnts[arr[i]]++;
            }


            curIdx += dep;
            string t_out = ofile;
            if (rep!=1){
                std::ostringstream oss;
                oss<<curRep;
                t_out += "_" +oss.str();
            }

            if (curRep < writes && fillret) {
                rare_map cntsMap;
                // fill map:
                for(uint i = 0; i < cnts.size(); i++){
                    if(cnts[i] != 0){
                        cntsMap.insert( std::make_pair(i, cnts[i]) );
                    }
                }
                RareSample[i].push_back(cntsMap);
                if(curRep == 0){
                    retCntsSampleName[i] = divs->SampleName; // safe the sample name as well
                }
            }
            richness = 0;
            divs->richness[i].push_back(this->getRichness(cnts));
            vector<double> three = this->calc_div(cnts, 4);
            divs->shannon[i].push_back(three[0]);
            divs->simpson[i].push_back(three[1]);
            divs->invsimpson[i].push_back(three[2]);
            divs->chao1[i].push_back(this->calc_chao1(cnts,1));
            divs->eve[i].push_back(this->calc_eveness(cnts));
            richness = 0;

            // save abundance for chao2 calculations later
            rarefyMutex.lock();
            for(uint im = 0; im < IDs.size(); im++){
                //sparse convertions in swap mode
                int id = std::stoi(IDs[im]);
                if(cnts[im] != 0){
                    abundInRow->at(i)[curRep][id]++;	
                    occuencesInRow->at(i)[curRep][id] += cnts[im];
                }
            }
            rarefyMutex.unlock();
        }
    }
}


// vector version of diversity fuctions:
long smplVec::getRichness(const vector<unsigned int>& cnts){
    for (size_t i = 0; i<cnts.size(); i++){
        //out<<IDs[i]<<"\t"<<cnts[i]<<endl;
        if (cnts[i]>0){
            richness++;
        }
    }
    return richness;
}


double smplVec::calc_chao1(const vector<uint> & vec,int corrBias){
    double Sobs((double)richness);
    double singl(0); double doubl(0);
    for (size_t i=0;i<vec.size();i++){
        if (vec[i]==1){singl++;}
        else if (vec[i]==2){doubl++;}
    }
    double est=0.0;
    if (corrBias==0){
        est = float( Sobs + (singl*singl)/(2*doubl) );
    } else {
        est = float( Sobs + (singl*(singl-1))/(2*(doubl+1)) );
    }
    /*if (conf.int){
      N = apply(M,2,sum)
      P = exp(-N/Sobs)
      P1 =Sobs/(1-P)
      P2 = 1.96*sqrt((Sobs*P)/(1-P))
      low =  P1 - P2
      idx = low<Sobs
      low[idx] = Sobs[idx]
      hi = P1 + P2
      }*/
    return (est);
}


vector<double> smplVec::calc_div(const vector<uint>& vec,int meth, float base){
    double sum = 0;
    for (size_t i=0; i<vec.size();i++){sum+=(double)vec[i];}
    vector<double> x(vec.begin(),vec.end());
    for (size_t i=0; i<x.size();i++){x[i] /= sum;}
    bool doexp = false;
    if (base <= 2.718284f && base >= 2.718280f){ // account for machine imprecission
        doexp = true;
    }
    vector<double> H(3, 0.0); double H1(0.0), H2(0.0), H3(0.0);
    if (meth == 1 || meth == 4){
        if (doexp){
            for (size_t i=0; i<x.size();i++){if (x[i]>0){H1 += x[i] * -log(x[i])  ;}}
        } else {
            float div = -log10(base);
            for (size_t i = 0; i<x.size(); i++){ if (x[i]>0){ H1 += x[i] * log10(x[i]) / div; } }
        }
        Shannon = H1;
    }
    if (meth == 3 || meth == 4 || meth == 2) {
        for (size_t i = 0; i<x.size(); i++){ H2 += x[i] * x[i]; }
        H3 = H2;
        H2 = 1 - H2;//simpson
        H3 = 1 / H3; //invsimpson
    }
    //for (size_t i=0; i<x.size();i++){H += x[i];}
    //if (meth == (int)2) {		H = 1 - H;	}else if (meth == 3)		H = 1/H;	}
H[0] = H1; H[1] = H2; H[2] = H3;

return(H);
}
double smplVec::calc_eveness(const vector<uint>& vec){
    //double sha = calc_div(vec,1);
    if (Shannon == -1.f){ vector<double> tm = calc_div(vec, 1); }
    return(Shannon / log((double)richness));
}

// map versions of diversity fuctions:
long smplVec::getRichness(rare_map& cnts){
    richness = cnts.size(); // set richness for other functions here
    return richness;
}

double smplVec::calc_chao1(rare_map& cnts,int corrBias){
    double Sobs((double)richness);
    double singl(0); double doubl(0);
    float est(0);
    typedef rare_map::iterator it;
    for(it iterator = cnts.begin(); iterator != cnts.end(); iterator++) {
        if (iterator->second == 1){
            singl++;
        }else if(iterator->second == 2){
            doubl++;
        }
    }

    /*
       for (size_t i=0;i<vec.size();i++){
       if (vec[i]==1){singl++;}
       else if (vec[i]==2){doubl++;}
       }
       */
    if (corrBias==0){
        est = float( Sobs + (singl*singl)/(2*doubl) );
    } else {
        est = float( Sobs + (singl*(singl-1))/(2*(doubl+1)) );
    }

    return est;
}


vector<double> smplVec::calc_div(rare_map& cnts,int meth, float base){
    double sum = 0;
    vector<double> H(3, 0.0); double H1(0.0), H2(0.0), H3(0.0);
    //for (size_t i=0; i<vec.size();i++){sum+=(double)vec[i];}
    sum = std::accumulate(std::begin(cnts), std::end(cnts), 0, [] (int value, const rare_map::value_type& p){ return value + p.second; });

    // copy counts into x
    unordered_map <uint, double> x;
    x.insert(cnts.begin(), cnts.end());
    typedef unordered_map <uint, double>::iterator it;
    for(it iterator = x.begin(); iterator != x.end(); iterator++) {
        // iterator->first = key
        iterator->second /= sum;
    }
    /*
       vector<double> x(vec.begin(),vec.end());
       for (size_t i=0; i<x.size();i++){x[i] /= sum;}*/


    bool doexp = false;
    if (base <= 2.718284f && base >= 2.718280f){ // account for machine imprecission
        doexp = true;
    }

    if (meth == 1 || meth == 4){
        if (doexp){
            typedef unordered_map <uint, double>::iterator it;
            for(it iterator = x.begin(); iterator != x.end(); iterator++) {
                if (iterator->second>0){H1 += iterator->second * -log(iterator->second)  ;}
            }
            //for (size_t i=0; i<x.size();i++){if (x[i]>0){H1 += x[i] * -log(x[i])  ;}}
        } else {
            float div = -log10(base);
            typedef unordered_map <uint, double>::iterator it;
            for(it iterator = x.begin(); iterator != x.end(); iterator++) {
                if (iterator->second>0){H1 += iterator->second * -log10(iterator->second)/ div;}
            }
            //for (size_t i = 0; i<x.size(); i++){ if (x[i]>0){ H1 += x[i] * log10(x[i]) / div; } }
        }
        Shannon = H1;
    }

    if (meth == 3 || meth == 4 || meth == 2) {
        typedef unordered_map <uint, double>::iterator it;
        for(it iterator = x.begin(); iterator != x.end(); iterator++) {
            H2 += iterator->second * iterator->second;
        }
        //for (size_t i = 0; i<x.size(); i++){ H2 += x[i] * x[i]; }
        H3 = H2;
        H2 = 1 - H2;//simpson
        H3 = 1 / H3; //invsimpson
    }
    H[0] = H1; H[1] = H2; H[2] = H3;

    return(H);
}
double smplVec::calc_eveness(rare_map& cnts){
    //double sha = calc_div(vec,1);
    if (Shannon == -1.f){ vector<double> tm = calc_div(cnts, 1); }
    return(Shannon / log((double)richness));
}






void smplVec::print2File(const vector<unsigned int>& cnts,const string t_out){
    richness=0;
    ofstream out(t_out.c_str());
    for (size_t i=0;i<cnts.size();i++){
        //out<<IDs[i]<<"\t"<<cnts[i]<<endl;
        if (cnts[i]>0){
            richness++;
            out<<IDs[i]<<"\t"<<cnts[i]<<endl;
        }
    }
    out.close();
    //cout<<"Richness: "<<richness<<endl;
    //return richness;
}

ulong thr_rng(unsigned long pos,MyRNG& rng) {
    std::uniform_int_distribution<unsigned long> uint_distx(0,pos);
    return uint_distx(rng);
}

/*
   void smplVec::shuffle(){
   time_t seed_val=time(NULL);           // populate somehow
   rng_P.resize(num_threads);
   for (long t=0; t<num_threads; t++){ //initialize N random number generators
   rng_P[t].seed((long)seed_val-(t*13));
   }
   std::vector<std::future<ulong>> futures(num_threads);
//auto *thr = new async[num_threads];
//vector<unsigned long> rngs(num_threads,0);
ulong trng;

unsigned long pos = (unsigned long)totSum- 1;
while ( pos > 0) {
//cout << pos<<" ";
//Launch a group of threads
for (int t = 0; t < num_threads; ++t) {
futures[t] = async(thr_rng,pos-t,rng_P[t]);
}
for (int t = 0; t < num_threads; ++t) {
trng = futures[t].get();
unsigned int temp = arr[pos] ;
arr[pos] = arr[trng];
arr[trng] = temp;
pos--;
}
//std::uniform_int_distribution<unsigned long> uint_distx(0,pos);
//unsigned long j = uint_distx(rng);
//swap(arr[i],arr[j]);
}
//delete [] thr;
if (verbose){
#ifdef notRpackage
cerr<<"fini";
#endif
}
}
*/
/*void smplVec::shuffle_singl(){
  time_t seed_val=time(NULL);           // populate somehow
  rng.seed((long)seed_val);
  unsigned long j; unsigned int temp;
  for (unsigned long i = 0 ; i < (unsigned long)(totSum - 1); i++) {
  std::uniform_int_distribution<unsigned long> uint_distx(0,i);
  j = uint_distx(rng);
  temp = arr[i] ;
  arr[i] = arr[j];
  arr[j] = temp;
//swap(arr[i],arr[j]);
}
if (verbose){
#ifdef notRpackage
//cerr<<"fini";
#endif
}
}*/

void smplVec::shuffle_singl() {
    //auto engine = std::default_random_engine{};
//    std::random_device rd;
//    auto engine = std::mt19937_64{rd()};
	unsigned long long seed = (unsigned long long)chrono::high_resolution_clock::now().time_since_epoch().count();
	auto engine = std::mt19937_64{ seed };
    std::shuffle(std::begin(arr), std::end(arr), engine);	
}


int smplVec::binarySearch( vector<float> vec, const float toFind)
{
    int len =(int) vec.size();
    // Returns index of toFind in sortedArray, or -1 if not found
    int low = 0;
    int high = len - 1;
    int mid;

    float l = vec[low];
    float h = vec[high];

    while (l < toFind && h > toFind) {
        mid = (low + high)>>1;

        float m = vec[mid];

        if (m < toFind) {
            if (vec[mid+1] > toFind){return mid;}
            l = vec[low = mid + 1];
        } else if (m > toFind) {
            if (vec[mid-1] < toFind){return mid-1;}
            h = vec[high = mid - 1];
        } else {
            return mid;
        }
        if (mid==len || mid==0){return mid;}
    }
    if (vec[low] == toFind)
        return low;
    else
        return -1; // Not found
}



/*
void DivEsts::print2file(const string file){
    if (richness.size()<1){return;}
    ofstream out(file.c_str());
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open diversity estimate file " << file << endl; std::exit(99);
#endif
    }
    out<<"Richness\t"<<richness[0];
    for (size_t i=1; i<richness.size();i++){
        out << "\t"<<richness[i];
    }
    out<<"\nShannon\t"<<shannon[0];
    for (size_t i=1; i<shannon.size();i++){
        out << "\t"<<shannon[i];
    }
    out<<"\nSimpson\t"<<simpson[0];
    for (size_t i=1; i<simpson.size();i++){
        out << "\t"<<simpson[i];
    }
    out<<"\nInv. Simpson\t"<<invsimpson[0];
    for (size_t i=1; i<invsimpson.size();i++){
        out << "\t"<<invsimpson[i];
    }
    out<<"\nChao1\t"<<chao1[0];
    for (size_t i=1; i<chao1.size();i++){
        out << "\t"<<chao1[i];
    }
    out<<"\nEveness\t"<<eve[0];
    for (size_t i=1; i<eve.size();i++){
        out << "\t"<<eve[i];
    }
    out.close();
}*/
void printDivMat(const string outF, vector<DivEsts*>& inD, bool printDIV, options* opts){

    string outFmedian = outF + "median_alpha_diversity.tsv";
    ofstream out(outFmedian.c_str());
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open diversity estimate matrix " << outF << endl; std::exit(99);
#endif
    }
    out << "depth";
    for( uint di = 0; di < opts->depth.size(); di++){
        for( uint j = 0; j < 6; j++){
                out << "\t" << opts->depth[di];   
        }
    }
    out << "\n";
    out << "Smpl";
    for( uint di = 0; di < opts->depth.size(); di++){
       out << "\tRichness\tShannon\tSimpson\tInv. Simpson\tChao1\tEveness";
    }
    out << "\n";
    for (size_t i = 0; i < inD.size(); i++){
        if (inD[i] == NULL){
#ifdef notRpackage
            cerr << "Empty vector at index " << i << " in div mat building.\n";
#endif
            out << "-1\t-1\t-1\t-1\t-1\t-1\n";
            continue;
        }
        out << inD[i]->SampleName;
        for( uint di = 0; di < opts->depth.size(); di++){
            out << "\t" << getMedian(inD[i]->richness[di]) ;
            out  << "\t"<< getMedian(inD[i]->shannon[di]);
            out  << "\t"<< getMedian(inD[i]->simpson[di]);
            out  << "\t"<< getMedian(inD[i]->invsimpson[di]) ;
            out  << "\t"<< getMedian(inD[i]->chao1[di]) << "";
            out  << "\t" << getMedian(inD[i]->eve[di]) << "";
        }
        out << "\n";
    }
    out.close();

    // print now each div estimate as well:
    if(printDIV){
        // open all files as streams
        vector<string> divNames;
        divNames.push_back("richness");
        divNames.push_back("shannon");
        divNames.push_back("simpson");
        divNames.push_back("invsimpson");
        divNames.push_back("chao1");
        divNames.push_back("eve");

        vector<ofstream> outFs(divNames.size());

        // open files
        for(uint i = 0; i < divNames.size(); i++){
            string outFdiv = outF + "alpha_" + divNames[i] + ".tsv";
            outFs[i].open(outFdiv.c_str(), ios_base::out);

            // write depth in first line of each file
            // usefull for when more than one depth was
            // requested
            outFs[i] << "depth" ;
            for(uint ii = 0; ii < opts->depth.size(); ii++){
                for(uint ij = 0; ij < opts->repeats ; ij++){
                    outFs[i] << "\t" << opts->depth[ii];
                }
            }
            outFs[i] << std::endl;
            
            outFs[i] << "repeat" ;
            for(uint ii = 0; ii < opts->depth.size(); ii++){
                for(uint ij = 0; ij < opts->repeats ; ij++){
                    outFs[i] << "\t" << ij+1;
                }
            }
            outFs[i] << std::endl;
        }

        // write the divvs to disk
        for (size_t i = 0; i < inD.size(); i++){
            // richness
            uint k = 0;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->richness[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->richness[di][j] ;
                }
            }
            outFs[k] << '\n';

            // shannon
            k = 1;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->shannon[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->shannon[di][j] ;
                }
            }
            outFs[k] << '\n';

            // simpson
            k = 2;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->simpson[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->simpson[di][j] ;
                }
            }
            outFs[k] << '\n';

            // invsimpson
            k = 3;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->invsimpson[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->invsimpson[di][j] ;
                }
            }
            outFs[k] << '\n';

            // chao1
            k = 4;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->chao1[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->chao1[di][j] ;
                }
            }
            outFs[k] << '\n';

            // eve
            k = 5;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->chao1[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->eve[di][j] ;
                }
            }
            outFs[k] << '\n';
        }

        // close streams
        for(uint i = 0; i < divNames.size(); i++){
            outFs[i].close();
        }
    }

}
void printRareMat(const string outF, const vector< rare_map>& rMat, vector< string >& sampleNames, vector < string >& rowId){
    ofstream out(outF.c_str());
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open rarefy matrix file " << outF << endl; std::exit(99);
#endif
    }

    // write the header
    out << "Rarefied";
    for(uint i = 0; i < sampleNames.size(); i++){
        out << "\t"<<sampleNames[i];
    }
    out << "\n";

    // write the tsv body
    for(uint i = 0; i < rowId.size(); i++){
        out << rowId[i] << "\t";
        for(uint j = 0; j < sampleNames.size(); j++){
            auto fnd = rMat[j].find(i);
            if(fnd != rMat[j].end()){
                out << "\t" << fnd->second;
            }else{
                out << "\t0";
            }
        }
        out << "\n";
    }

    out.close();
}



string printSimpleMap(const rare_map & vec, string outF, string id, vector<string> rowNames){
    // takes a map from the rarefaction function and writes the vector
    // to the disk.
    // this way we dont need memory to do
    ofstream out(outF.c_str(),  ios::binary);
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open tmpvec file " << outF << endl; std::exit(99);
#endif
    }
    for(uint i = 0; i < rowNames.size(); i++){
        uint value = 0;
        auto fnd = vec.find(i);
        if(fnd != vec.end()){
            value = fnd->second;
        }
        out.write((char*) &value, sizeof(uint));

    }
    out.close();
    /*
    // write the header
    out << id;
    out << "\n";

    // write the vector
    for(uint i = 0; i < rowNames.size(); i++){
    auto fnd = vec.find(i);
    if(fnd != vec.end()){
    out << fnd->second;
    }else{
    out << "0";
    }
    out << "\n";
    }
    out.close();
    */
    return outF;
}

void reassembleTmpMat(vector<string> inF, vector< string > rowNames, vector< string > colNames, string outF){
    // takes the vectors from printSimpleMap and constrcust a mat from them
    // first open all inF streams
    // iterate through and write line for line
    if(inF.size() == 0){
#ifdef notRpackage
        std::exit(99);
#endif
    }

    vector<std::ifstream*> inFs(inF.size());
    for(uint i = 0; i < inFs.size(); i++){
        std::ifstream* f = new std::ifstream(inF[i].c_str(), std::ios::in | std::ios::binary); // create in free store
        inFs[i] = f;
        inFs[i]->close();
    }

    ofstream out(outF.c_str());
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open tmpvec file " << outF << endl; std::exit(99);
#endif
    }
    out << "Rarefied";
    for(uint i = 0; i < colNames.size(); i++){
        out << '\t' << colNames[i];
    }
    out << '\n';

    string a;
    uint j = 0;
    // bufer for 1000 rows
    uint bj = 0;
    const uint bn = 1000;
    std::vector<vector< uint > > inBuff(bn, std::vector<uint>(inFs.size() ) );

    while(j < rowNames.size()){
        // fill buffer:
        for(uint i = 0; i < inFs.size(); i++){
            uint value;
            inFs[i]->open(inF[i].c_str(), std::ios::in | std::ios::binary);
            long int offset = j * sizeof(value);
            inFs[i]->seekg(offset, ios_base::beg );

            bj = 0;
            while((*inFs[i]) && bj < bn ){
                inFs[i]->read(reinterpret_cast<char*>(&value), sizeof(value));
                inBuff[bj][i] = value;
                bj++;
            }

            inFs[i]->close();
        }
        // write buffers to file
        for(uint ij = 0; ij < bj && j+ij < rowNames.size(); ij++){
            out <<  rowNames[j + ij];
            for(uint i = 0; i < inFs.size(); i++){
                out << '\t' << inBuff[ij][i];
            }
            out << '\n';
        }
        j = j + bj;

    }
    out.close();

}


std::istream& safeGetline(std::istream& is, std::string& t)
{
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



void computeChao2(std::vector<vector<mat_fl>>& chao2, vector<vector<vector<uint>>>& abundInRow){
    for(uint i = 0; i < abundInRow.size(); i++){
        for(uint ii = 0; ii < abundInRow[i].size(); ii++){
            // count No of species
            float NoOfSpec = 0;
            float singletons = 0;
            float doubletons = 0;
            for(uint j = 0; j < abundInRow[i][ii].size(); j++){
                if(abundInRow[i][ii][j] != 0){
                    NoOfSpec++;
                }
                if(abundInRow[i][ii][j] == 1){
                    singletons++;
                }else	if(abundInRow[i][ii][j] == 2){
                    doubletons++;
                }
            }
            // calc chao2
            mat_fl tmpChao2 = 0.0;
            if(doubletons != 0){
                tmpChao2 = float(NoOfSpec + ((singletons*singletons)/(2*doubletons)));
            }
            chao2[i].push_back(tmpChao2);
        }   
    }
}

void computeCE(vector<vector<mat_fl>>& CE, vector<vector<vector<uint>>>& abundInRow){
    // inspired by the ICE implementation in the R package fossil
    // by Matthew Vavrek
    // https://matthewvavrek.com/programs-and-code/fossil/
    // ACE and IE use this functio, one with ror abundance the other with row presence data
    int val;
    for(uint i = 0; i < abundInRow.size(); i++){
        for(uint ii = 0; ii < abundInRow[i].size(); ii++){
            std::vector< int > abundOneToTen(10,0);
            float nr = 0.0, sa = 0.0, sr = 0.0, f1 = 0.0, ca= 0.0,sumf= 0.0, g2a = 0.0;

            for(uint j = 0; j < abundInRow[i].size(); j++){
                val = abundInRow[i][ii][j];
                if(val < 11 && val != 0){
                    nr += val;
                    sr++;
                }else if(val > 10){
                    sa++;
                }
                if(val == 1){
                    f1 += 1;
                }
                if((val < 11) && (val > 0)){
                    abundOneToTen[val-1]++;
                }
            }

            for(uint j = 0; j < abundOneToTen.size(); j++){
                sumf += abundOneToTen[j] * (j+1);
            }
            ca = 1 - (f1)/(nr);
            g2a = (sr/ca) * (sumf/(nr * (nr - 1))) - 1;
            if(g2a < 0){
                g2a = 0;
            }

            if(ca != 0){
                mat_fl tmp = sa + sr/ca + (f1/ca) * g2a;
                CE[i].push_back(tmp);
            }else{
                CE[i].push_back(0);// or compute chao2 here, why would i do that?
            }
        }   }
}


void writeGlobalDiv(options* opts, vector<vector<mat_fl>>& ICE, vector<vector<mat_fl>>& ACE, vector<vector<mat_fl>>& chao2, string outF){
    ofstream out(outF.c_str());
    out << "depth";

    for(uint j = 0; j < opts->depth.size(); j++){
        for(uint i = 0; i < opts->repeats; i++){
            out << "\t" << opts->depth[j];
        }
    }    
    out << '\n';
    
    
    out << "repeat";

    for(uint j = 0; j < opts->depth.size(); j++){
        for(uint i = 0; i < opts->repeats; i++){
            out << "\t" << i+1;
        }
    }    
    out << '\n';

    out << "Chao2";                      
    for(uint i = 0; i < chao2.size(); i++){
        for(uint j = 0; j < chao2[i].size(); j++){
            out << '\t' << chao2[i][j];
        }                            
    }
    out << '\n';

    out << "ICE";
    for(uint i = 0; i < ICE.size(); i++){
        for(uint j = 0; j < ICE.size(); j++){
            out << '\t' << ICE[i][j];
        }   }
    out << '\n';

    out << "ACE";
    for(uint i = 0; i < ACE.size(); i++){
        for(uint j = 0; j < ACE.size(); j++){
            out << '\t' << ACE[i][j];
        }    }
    out << '\n';
    out.close();
}
