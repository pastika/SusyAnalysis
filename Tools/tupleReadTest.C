#include "StopTupleReader.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <iostream>

int main()
{
    TFile *f = new TFile("/eos/uscms/store/user/pastika/DYJetsToLL_M-50_13TeV-madgraph-pythia8/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141227_223539/0000/stopFlatNtuples_4.root");
    
    TTree *t = (TTree*)f->Get("stopTreeMaker/AUX");

    NTupleReader tr(t);

    while(tr.getNextEvent())
    {
        //std::cout << tr.getNEntries() << "\t" << tr.run << "\t" << tr.lumi << "\t" << tr.event << std::endl;
    }
}
