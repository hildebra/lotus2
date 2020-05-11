#pragma once
#include "ClStr2Mat.h"
#include "options.h"

struct rareStruct{
	int i;
	DivEsts* div;
	vector<string> cntsName;
	vector<vector<rare_map>> cnts;
	string skippedNames;
	vector<string> IDs;
	
};

struct job {
  std::future <rareStruct*> fut;
  bool inUse = false;
};

void binaryStoreSample(options* opts, vector<vector< vector< string >> > & , rareStruct* , 
	vector<string>& , string , vector<vector<string>>& , bool reshapeMap = false);
void memoryStoreSample(options* opts, rareStruct* tmpRS, vector< vector< vector< rare_map >> >& MaRare,  vector<vector<string>>& cntsNames, bool reshapeMap);

void printRarefactionMatrix(options* , vector<vector<vector< string >>>& , string,  
                            vector<vector<string>>& , vector<string>&);
void printRarefactionMatrix(options*, const vector<vector<vector< rare_map>>>& , 
                            string , vector<vector<string>>& , vector<string>& );
