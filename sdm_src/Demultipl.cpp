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



int main(int argc, char* argv[])
{
	if (argc<3){
		//help_options,help_map,help_commands
		if (argc==2){
			if (string(argv[1])=="-help_commands"){
				printCmdsHelp();
			} else if (string(argv[1])=="-help_options"){
				printOptionHelp();
			} else if (string(argv[1])=="-help_map"){
				printMapHelp();
			}
			else if (string(argv[1]) == "-version" || string(argv[1]) == "-v"){
				printVersion();
			}
			exit(0);
		}
		general_help();
		exit(0);
	}
#ifdef DEBUG
	cerr << "DEBUG mode"<<endl;
#endif

	ini_DNAconstants();

	Announce_sdm();

	OptContainer cmdArgs;
	readCmdArgs(argc, argv, cmdArgs);
	 
	//rewrites header names
	rewriteNumbers(cmdArgs);

	//reads the sdm_options file
#ifdef DEBUG
	cerr << "Setting up Filter" << endl;
#endif
	shared_ptr<Filters> fil = make_shared<Filters>(cmdArgs);
	bool bReads = fil->readMap(cmdArgs);
#ifdef DEBUG
	cerr << "filter setup & map is read" << endl;
#endif
	if (!bReads){cerr<<"Failed to read Map.\n";exit(3);}
	//cerr<<SAqualP[0]<<" "<<SAqualP[50];
	fil->setcmdArgsFiles(cmdArgs);
	
	clock_t tStart = clock();
	//main function
	
	separateByFile(fil,cmdArgs);
	//cerr << "\nXXXX\n\n";
//	delete fil;

	fprintf(stderr,"Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	return 0;
}



