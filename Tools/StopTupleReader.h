#ifndef STOP_NTUPLE_READER_H
#define STOP_NTUPLE_READER_H

#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"

#include <vector>
#include <cstdio>

class NTupleReader
{

public:
    // List of all variables used in tuple
    int run, lumi, event;

    NTupleReader(TTree * tree);

    int getNEntries()
    {
        if(tree_) return tree_->GetEntries();
        else 
        {
            printf("NTupleReader::getNEntries() - NO tree defined yet!!!\n");
            return -1;
        }
    }

    bool getNextEvent();

private:
    // private variabals for internal use
    TTree *tree_;
    int nevt_;

    void activateBranches();
    void clearTuple();
};

#endif
