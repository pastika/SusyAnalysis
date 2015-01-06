#!/bin/tcsh

set CMSRELDIR = /uscms/home/javiert/work/CMSSW_5_2_5
set workDir = /uscms/home/javiert/work/CMSSW_5_2_5/src/SusyAnalysis/LostLepton/test
set pythonConfFile = /uscms/home/javiert/work/CMSSW_5_2_5/src/SusyAnalysis/LostLepton/test/default_cfg.py
set listDir = /uscms/home/javiert/work/CMSSW_5_2_5/src/SusyAnalysis/LostLepton/test/lists
set outBaseDir = /uscms_data/d3/javiert

echo =============================
date
echo Using: $CMSRELDIR
echo jobID: 
echo CMSRELDIR: $CMSRELDIR
echo workDir: $workDir
echo pythonConfFile: $pythonConfFile
echo =============================
echo
echo Setting Fermilab CMS environment..
source /uscmst1/prod/sw/cms/cshrc prod
cd  $CMSRELDIR
pwd
cmsenv                                 
cmscvsroot CMSSW
echo CVSROOT: $CVSROOT
echo PWD `pwd`
echo Done setting Fermilab CMS environment.

echo
echo Actual job begins..
echo

cd $workDir
echo ARG1: $1
echo ARG2: $2
echo ARG3: $3
echo ARG4: $4
echo
#cmsRun $pythonConfFile inputFiles=dcap:///pnfs/cms/WAX/11/store/user/lpcsusyhad/dhare/TTJets_TuneZ2star_8TeV-madgraph-tauola/SUSYPAT_76_1_mgR.root outputFile=out.root 

cmsRun $pythonConfFile fileList=$listDir/$1 outputFile=$outBaseDir/$2_$1.root

echo
echo End of the Job!
echo
date
echo
