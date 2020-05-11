#pragma once
#include "RefTax.h"

//routine that calls various subfunctions like blast filter, LCA etc
TaxObj* LCA(list<BlastRes*>, RefTax*, options*);

//filtering of the blast results
double filterBlastPrimary(list<BlastRes*>&, bool maxHitOnly=false);
//routine that performs actual LCA matching etc
TaxObj* LCAcore(list<TaxObj*>, bool &hitRd, double LCAfrac=0.9f, int tdepth=__default_depth);


list<TaxObj*> BlastToTax(list<BlastRes*>& BR, RefTax* RT, options*, float& consPerID);