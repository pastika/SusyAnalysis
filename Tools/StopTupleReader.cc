#include "StopTupleReader.h"
#include <iostream>

NTupleReader::NTupleReader(TTree * tree)
{
    tree_ = tree;
    nevt_ = 0;
    clearTuple();
    activateBranches();
    //getNextEvent();
}

void NTupleReader::activateBranches()
{
    tree_->SetBranchStatus("*", 0);

    tree_->SetBranchStatus("run", 1);   tree_->SetBranchAddress("run", &run);
    tree_->SetBranchStatus("event", 1); tree_->SetBranchAddress("event", &event);
    tree_->SetBranchStatus("lumi", 1);  tree_->SetBranchAddress("lumi", &lumi);
}

bool NTupleReader::getNextEvent()
{
    if(nevt_ >= getNEntries())
    {
        return false;
    }
    clearTuple();
    int status = tree_->GetEntry(nevt_);
    std::cout << status << "\t" << run << std::endl;
    nevt_++;
    return status > 0;
}

void NTupleReader::clearTuple()
{
    run = lumi = event = 0;
}

