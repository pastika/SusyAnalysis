#!/bin/tcsh

date

source /uscmst1/prod/sw/cms/cshrc prod
setenv SCRAM_ARCH slc5_amd64_gcc462

          #shown for Fermilab UAF users as example for remote site
cd  /uscms_data/d2/vasurang/Run2012/SUSYSTop/CMSSW_5_2_5/src/
echo "Using  /uscms_data/d2/vasurang/Run2012/SUSYSTop/CMSSW_5_2_5/src/"
cmsenv                                  #shown for C shell
cmscvsroot CMSSW
echo "CVSROOT is"
echo $CVSROOT
echo "PWD is: "
pwd


cd CURDIR 


cmsRun CONFIGFILE

