#include "Matrix.h"



Matrix::Matrix(int depth, vector<string> taxs, bool reportRead):mat(depth), 
		colIDs(taxs), rowIDs(depth), readReport(reportRead)
{
	if (readReport){//means that not only species level is reported, but also singular hit (shini table)
		mat.resize(depth + 1);
		rowIDs.resize(depth + 1); colIDs.resize(depth + 1, "Reads");
	}
}
Matrix::~Matrix()
{
}
void Matrix::add(TaxObj* t) {
	string addStr("");// t->get(0));
	for (size_t DL = 0; DL < rowIDs.size(); DL++) {
		string curT = t->get(DL);
		if (curT == __unkwnTax) {curT = __unkwnTaxWR;}
		if (DL == 0) {
			addStr += curT;
		} else {
			addStr += __taxSepMat + curT;
		}
		auto fnd = rowIDs[DL].find(addStr);
		if (fnd == rowIDs[DL].end()) {
			rowIDs[DL][addStr] = 1;
		} else {
			rowIDs[DL][addStr]++;
		}
	}
}
void Matrix::writeAllLevels(const string& outF) {
	
	for (size_t DL = 0; DL < rowIDs.size(); DL++) {
		string outF1 = outF + "_" + colIDs[DL];
		ofstream of(outF1.c_str());
		//for (auto i : rowIDs[DL]) {
		for (auto it = rowIDs[DL].begin(); it != rowIDs[DL].end(); ++it) {
			of << it->first<< __MatSep << it->second << endl;
		}
		of.close();
	}
}
