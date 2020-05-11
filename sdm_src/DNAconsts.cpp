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

#include "DNAconsts.h"
char DNA_trans[256];
short DNA_amb[256];//to count amb chars
short DNA_IUPAC[256 * 256];
short NT_POS[256];


void ini_DNAconstants(){
//static char* DNA_trans[256] = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";	
	//DNA_trans.resize(256,'X');
	for (int i = 0; i<256; i++){ DNA_trans[i] = 'N'; }
	DNA_trans['A'] = 'T';	DNA_trans['T'] = 'A';
	DNA_trans['C'] = 'G'; DNA_trans['G'] = 'C';
	DNA_trans['a'] = 'T';	DNA_trans['t'] = 'A';
	DNA_trans['c'] = 'G'; DNA_trans['g'] = 'C';
	for (int i = 0; i<256; i++){ DNA_amb[i] = 0; }
	DNA_amb['A'] = 1;	DNA_amb['T'] = 1;
	DNA_amb['C'] = 1; DNA_amb['G'] = 1;
	DNA_amb['a'] = 1;	DNA_amb['t'] = 1;
	DNA_amb['c'] = 1; DNA_amb['g'] = 1;


	for ( int i = 0; i < 256 ; i++ ) {  NT_POS[i] = 5; }
	NT_POS['A'] = 0; NT_POS['T'] = 1; NT_POS['G'] = 2; NT_POS['C'] = 3;
	NT_POS['N'] = 4;
	NT_POS['a'] = 0; NT_POS['t'] = 1; NT_POS['g'] = 2; NT_POS['c'] = 3;
	NT_POS['n'] = 4;
	for ( int i = 0; i < 256 * 256; i++ ) { DNA_IUPAC[i] = 1;  }
	for ( int i = 0; i<14; i++ ) {//first: N is always a hit
		DNA_IUPAC['N'+256*DNA_SPACE[i]]= 0;
		DNA_IUPAC[256*'N'+DNA_SPACE[i]]= 0;
	}
	for ( int i = 0; i<5; i++ ) {//first: N is always a hit
		DNA_IUPAC['B'+256*DNA_SPACE[i]]=0;
		DNA_IUPAC[256*'B'+DNA_SPACE[i]]=0;
		DNA_IUPAC['H'+256*DNA_SPACE[i]]=0;
		DNA_IUPAC[256*'H'+DNA_SPACE[i]]=0;
		DNA_IUPAC['D'+256*DNA_SPACE[i]]=0;
		DNA_IUPAC[256*'D'+DNA_SPACE[i]]=0;
		DNA_IUPAC['V'+256*DNA_SPACE[i]]=0;
		DNA_IUPAC[256*'V'+DNA_SPACE[i]]=0;
	}

	DNA_IUPAC['B'+256*'A']=1;DNA_IUPAC[256*'B'+'A']=1;
	DNA_IUPAC[('D'+256*'C')]=1;DNA_IUPAC[256*'D'+'C']=1;
	DNA_IUPAC['H'+256*'G']=1;DNA_IUPAC[256*'H'+'G']=1;
	DNA_IUPAC['V'+256*'T']=1;DNA_IUPAC[256*'V'+'T']=1;



	DNA_IUPAC['R'+256*'A']=0;DNA_IUPAC[256*'R'+'A']=0;
	DNA_IUPAC['R'+256*'G']=0;DNA_IUPAC[256*'R'+'G']=0;
	DNA_IUPAC['M'+256*'C']=0;DNA_IUPAC[256*'M'+'C']=0;
	DNA_IUPAC['M'+256*'A']=0;DNA_IUPAC[256*'M'+'A']=0;
	DNA_IUPAC['Y'+256*'C']=0;DNA_IUPAC[256*'Y'+'C']=0;
	DNA_IUPAC['Y'+256*'T']=0;DNA_IUPAC[256*'Y'+'T']=0;
	DNA_IUPAC['K'+256*'G']=0;DNA_IUPAC[256*'K'+'G']=0;
	DNA_IUPAC['K'+256*'T']=0;DNA_IUPAC[256*'K'+'T']=0;
	DNA_IUPAC['W'+256*'A']=0;DNA_IUPAC[256*'W'+'A']=0;
	DNA_IUPAC['W'+256*'T']=0;DNA_IUPAC[256*'W'+'T']=0;
	DNA_IUPAC['S'+256*'C']=0;DNA_IUPAC[256*'S'+'C']=0;
	DNA_IUPAC['S'+256*'G']=0;DNA_IUPAC[256*'S'+'G']=0;

	DNA_IUPAC['A'+256*'A']=0;DNA_IUPAC[256*'A'+'A']=0;
	DNA_IUPAC['T'+256*'T']=0;DNA_IUPAC[256*'T'+'T']=0;
	DNA_IUPAC['G'+256*'G']=0;DNA_IUPAC[256*'G'+'G']=0;
	DNA_IUPAC['C'+256*'C']=0;DNA_IUPAC[256*'C'+'C']=0;
	//DEBUG
	
	//fake use to surpress compiler warnings
	//if (sdm_version == 0.f){ cout << "too low version" << sdm_status; }
	//if (strcmp(sdm_status, "XX")){ cout << "Some"; }

}