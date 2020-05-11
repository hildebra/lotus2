#include "options.h"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

void self_help() {
	cout << "LCA help\nUsage: ./LCA [[ optional_args ]] -i [blast m8 output] -r [taxonomy database] -o [output file]\n";
	cout << "Required arguments:\n";
	cout << "  -i mapping assignments of sequences to ref database in blast .m8 tab delimited format\n";
	cout << "  -r taxonomy file with entries corresponding to sequences in ref database, that was mapped against\n";
	cout << "  -o output file containing the sequence name and the assigned taxonomy against the ref database\n";
	cout << "Optional arguments:\n";
	cout << "  -matHigh       calculate abundance of reads at different taxonomic levels. An extra file (derriving from -o) per tax level is written\n";
	cout << "  -showHitRead   if a hit can be uniquely assigned to a single entry in the ref database, this is reported in the -o file.\n";
	cout << "  -no_bl_filter  use only, if custom scripts were used to pre-filter filter -i file and in-built filter should be deactivated\n";
	//useless and deactivated
	//cout << "  -annotateAll (gives back an annotation for every single read/otu\n";
	cout << "  -readInput     [miTag / OTU] changes the tags attached to single reads\n";
	cout << "  -LCAfrac       [0-1] the fraction of matching taxonomies required to accept this taxonomy on the different levels. Default=\"0.8\"\n";
//	cout << "  -t             number of threads to use, currently deactivated on linux as it requires C++11\n";
	cout << "  -id            comma seperated list of min %identity, to accept a database hit as applicable to this taxonomic level, starting from Species and going to Kingdom. Default=\"97,95,93,91,88,78,0\"\n";
	cout << endl;
	exit(0);
}


options::options(int argc, char **argv,int defDep):
	RefTaxFile(""), blastres(""), outF(""), input_format("bl8"),
	BLfilter(true), calcHighMats(false), hitRD(false), isReads(false),
	annotateAll(false), nativeSlVdb(false), checkTaxoUnkw(true),
	numThr(1), taxDepth(defDep), LCAfract(0.9f), idThr(defDep,0),
	blFiles(0), refDBs(0), Taxlvls(defDep), reportBestHit(false)
{
	idThr[1] = 78; idThr[2] = 88; idThr[3] = 91; idThr[4] = 93;
	idThr[5] = 95; idThr[6] = 97;

	bool newIDthrs = false; string newIDthStr("");

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i"))
			blastres = argv[++i];
		else if (!strcmp(argv[i], "-h"))
			self_help();
		else if (!strcmp(argv[i], "-r"))
			RefTaxFile = argv[++i];
		else if (!strcmp(argv[i], "-f"))//input format
			input_format = argv[++i];
		else if (!strcmp(argv[i], "-o"))
			outF = argv[++i];
		else if (!strcmp(argv[i], "-matHigh"))
			calcHighMats = true;
		else if (!strcmp(argv[i], "-showHitRead"))
			hitRD = true;
		else if (!strcmp(argv[i], "-no_bl_filter"))
			BLfilter = false;
		else if (!strcmp(argv[i], "-no_taxDB_filter"))
			checkTaxoUnkw = false;
		else if (!strcmp(argv[i], "-readInput"))
			isReads = true;
		else if (!strcmp(argv[i], "-reportBestHit"))
			reportBestHit = true;
		else if (!strcmp(argv[i], "-SLVfmt"))
			nativeSlVdb = true;
		else if (!strcmp(argv[i], "-annotateAll"))
			annotateAll = true;
		else if (!strcmp(argv[i], "-t"))
			numThr = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-LCAfrac"))
			LCAfract = atof(argv[++i]);
		else if (!strcmp(argv[i], "-id")) {
			newIDthrs = true; newIDthStr= argv[++i];
		}

	}

	if (hitRD) {//needs to add 1 extra entry to some vectors
		/*taxDepth++;
		idThr.resize(taxDepth, 0);
		Taxlvls.resize(taxDepth,"Read")/*/
	}
	split(blastres, ',', blFiles);
	split(RefTaxFile, ',', refDBs);

	if (blFiles.size() != refDBs.size()) {
		cerr << "Unequal number of blast and refDB files!\n"; exit(24);
	}

	//check that all args were given
	bool isErr(false);
	if (blastres == "") { cerr << "Input file invalid or missing (-i)\n"; isErr = true; }
	if (RefTaxFile == "") { cerr << "RefDb file invalid or missing (-r)\n"; isErr = true; }
	if (outF == "") { cerr << "Output file invalid or missing (-o)\n"; isErr = true; }
	//if (blastres == "") { cerr << "Input file invalid (-f)"; isErr = true; }
	if (isErr) { cerr << "Use \"./LCA -h\" to get full help.\nError in command line args.. exiting\n"; exit(5); }

	if (newIDthrs) {//todo implement
		//cerr << "ID thresh inplement\n"; exit(55);
		vector<string> idthrsrev;
		split(newIDthStr, ',',idthrsrev);
		if (idthrsrev.size() != 7) {
			cerr << "Wrong number of id threshold levels (needs to be 7).\nAborting..\n"; exit(39);
		}
		for (size_t i = 0; i < 7; i++) {
			idThr[i] = atoi(idthrsrev[6 - i].c_str());
		}
	}

	vector<string> defTLvls(8, "");
	defTLvls[0] = "Domain"; defTLvls[1] = "Phylum"; defTLvls[2] = "Class";
	defTLvls[3] = "Order"; defTLvls[4] = "Family";  defTLvls[5] = "Genus";
	defTLvls[6] = "Species"; defTLvls[7] = "Strain";
	for (size_t i = 0; i < (size_t)taxDepth; i++) {
		if (i >= defTLvls.size()) { break; }
		Taxlvls[i] = defTLvls[i];
	}

}
const string options::TaxLvl2string() {
	string ret = Taxlvls[0];
	for (size_t i = 1; i < Taxlvls.size(); i++) {
		ret += __defaultTaxSep + Taxlvls[i];
	}
	return ret;
}


options::~options()
{
}
