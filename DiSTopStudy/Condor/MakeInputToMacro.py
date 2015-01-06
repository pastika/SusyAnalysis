import parser
import os
import shutil
import ROOT


mcStr= [   "TTbar", "WJetsToLNu",   "Zinv_400HTinf",   "Zinv_200HT400",   "DYJets",    "QCDFlat",    "T_t",    "T_tW",
                                "Signal_x500_y100",       "Signal_x450_y150",        "Signal_x350_y50",          "Signal_x500_y50",
                                "Signal_mStop250_mLSP0",  "Signal_mStop250_mLSP50",
                                "Signal_mStop300_mLSP0",  "Signal_mStop300_mLSP50",  "Signal_mStop300_mLSP100",
                                "Signal_mStop350_mLSP0",  "Signal_mStop350_mLSP50",  "Signal_mStop350_mLSP100",  "Signal_mStop350_mLSP150",
                                "Signal_mStop400_mLSP150","Signal_mStop400_mLSP200",
                                "Signal_mStop450_mLSP0",  "Signal_mStop450_mLSP50",  "Signal_mStop450_mLSP100",  "Signal_mStop450_mLSP150",  "Signal_mStop450_mLSP200",  "Signal_mStop450_mLSP250",
                                "Signal_mStop500_mLSP0",  "Signal_mStop500_mLSP50",  "Signal_mStop500_mLSP100",  "Signal_mStop500_mLSP150",  "Signal_mStop500_mLSP200",  "Signal_mStop500_mLSP250",
                                "Signal_mStop550_mLSP0",  "Signal_mStop550_mLSP50",  "Signal_mStop550_mLSP100",  "Signal_mStop550_mLSP150",  "Signal_mStop550_mLSP200",  "Signal_mStop550_mLSP250",
                                "Signal_mStop600_mLSP0",  "Signal_mStop600_mLSP50",  "Signal_mStop600_mLSP100",  "Signal_mStop600_mLSP150",  "Signal_mStop600_mLSP200",  "Signal_mStop600_mLSP250",
                                "Signal_mStop650_mLSP0",  "Signal_mStop650_mLSP50",  "Signal_mStop650_mLSP100",  "Signal_mStop650_mLSP150",  "Signal_mStop650_mLSP200",  "Signal_mStop650_mLSP250",
                                "Signal_mStop700_mLSP0",  "Signal_mStop700_mLSP50",  "Signal_mStop700_mLSP100",  "Signal_mStop700_mLSP150",  "Signal_mStop700_mLSP200",  "Signal_mStop700_mLSP250",
                                "Signal_mStop750_mLSP0",  "Signal_mStop750_mLSP50",  "Signal_mStop750_mLSP100",  "Signal_mStop750_mLSP150",  "Signal_mStop750_mLSP200",  "Signal_mStop750_mLSP250",
                                "Signal_mStop850_mLSP0",  "Signal_mStop850_mLSP50",  "Signal_mStop850_mLSP100",  "Signal_mStop850_mLSP150",  "Signal_mStop850_mLSP200",  "Signal_mStop850_mLSP250" ]

mcTyp= [   "TTbar", "WJetsToLNu",   "Zinv_400HTinf",   "Zinv_200HT400",   "DYJets",    "QCDFlat",    "T_t",    "T_tW",
                                "Signal",       "Signal",        "Signal",          "Signal",
                                "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",  "Signal",  "Signal",
                                "Signal",  "Signal",  "Signal",  "Signal",  "Signal",  "Signal" ]

mcDir= [   "TTJets_TuneZ2star_8TeV-madgraph-tauola", "WJetsToLNu",   "ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph",   "ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgrap",   "DYJets",    "QCD_selTag_Pt-15to3000_TuneZ2_Flat",    "T_t",    "T_tW",
           
                                "T2tt_x500_y100_NOCUTS_09Aug2012V1",       "T2tt_x450_y150_NOCUTS_09Aug2012V1",        "T2tt_x350_y50_NOCUTS_09Aug2012V1",          "T2tt_x500_y50_NOCUTS_09Aug2012V1",
                                "T2tt_mStop250_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop250_mLSP50_NOCUTS_09Aug2012V1",
                                "T2tt_mStop300_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop300_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop300_mLSP100_NOCUTS_09Aug2012V1",
                                "T2tt_mStop350_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop350_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop350_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop350_mLSP150_NOCUTS_09Aug2012V1",
                                "T2tt_mStop400_mLSP150_NOCUTS_09Aug2012V1","T2tt_mStop400_mLSP200_NOCUTS_09Aug2012V1",
                                "T2tt_mStop450_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop450_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop450_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop450_mLSP150_NOCUTS_09Aug2012V1",  "T2tt_mStop450_mLSP200_NOCUTS_09Aug2012V1",  "T2tt_mStop450_mLSP250_NOCUTS_09Aug2012V1",
                                "T2tt_mStop500_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop500_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop500_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop500_mLSP150_NOCUTS_09Aug2012V1",  "T2tt_mStop500_mLSP200_NOCUTS_09Aug2012V1",  "T2tt_mStop500_mLSP250_NOCUTS_09Aug2012V1",
                                "T2tt_mStop550_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop550_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop550_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop550_mLSP150_NOCUTS_09Aug2012V1",  "T2tt_mStop550_mLSP200_NOCUTS_09Aug2012V1",  "T2tt_mStop550_mLSP250_NOCUTS_09Aug2012V1",
                                "T2tt_mStop600_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop600_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop600_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop600_mLSP150_NOCUTS_09Aug2012V1",  "T2tt_mStop600_mLSP200_NOCUTS_09Aug2012V1",  "T2tt_mStop600_mLSP250_NOCUTS_09Aug2012V1",
                                "T2tt_mStop650_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop650_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop650_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop650_mLSP150_NOCUTS_09Aug2012V1",  "T2tt_mStop650_mLSP200_NOCUTS_09Aug2012V1",  "T2tt_mStop650_mLSP250_NOCUTS_09Aug2012V1",
                                "T2tt_mStop700_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop700_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop700_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop700_mLSP150_NOCUTS_09Aug2012V1",  "T2tt_mStop700_mLSP200_NOCUTS_09Aug2012V1",  "T2tt_mStop700_mLSP250_NOCUTS_09Aug2012V1",
                                "T2tt_mStop750_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop750_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop750_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop750_mLSP150_NOCUTS_09Aug2012V1",  "T2tt_mStop750_mLSP200_NOCUTS_09Aug2012V1",  "T2tt_mStop750_mLSP250_NOCUTS_09Aug2012V1",
                                "T2tt_mStop850_mLSP0_NOCUTS_09Aug2012V1",  "T2tt_mStop850_mLSP50_NOCUTS_09Aug2012V1",  "T2tt_mStop850_mLSP100_NOCUTS_09Aug2012V1",  "T2tt_mStop850_mLSP150_NOCUTS_09Aug2012V1",  "T2tt_mStop850_mLSP200_NOCUTS_09Aug2012V1",  "T2tt_mStop850_mLSP250" ]

xSecArr = [  225.1967,     30.08,                6.26,           49.28,     3532.8149,  2.99913994E10,   47.0,  11.1773*2.0,
                                        0.0855847,                0.169668,                   0.807323,                0.0855847,
                            5.57596,                  5.57596,
                            1.99608,                  1.99608,                    1.99608,
                            0.807323,                 0.807323,                   0.807323,                0.807323,
                            0.35683,                  0.35683,
                          0.169668,                 0.169668,                   0.169668,                0.169668,                    0.169668,                     0.169668,
                            0.0855847,                0.0855847,                  0.0855847,               0.0855847,                   0.0855847,                    0.0855847,
                            0.0452067,                0.0452067,                  0.0452067,               0.0452067,                   0.0452067,                    0.0452067,
                            0.0248009,                0.0248009,                  0.0248009,               0.0248009,                   0.0248009,                    0.0248009,
                            0.0139566,                0.0139566,                  0.0139566,               0.0139566,                   0.0139566,                    0.0139566,
                            0.0081141,                0.0081141,                  0.0081141,               0.0081141,                   0.0081141,                    0.0081141,
                            0.00480639,               0.00480639,                 0.00480639,              0.00480639,                  0.00480639,                   0.00480639,
                            0.00176742,               0.00176742,                 0.00176742,              0.00176742,                  0.00176742,                   0.00176742 ]
                                                                                                             


inputtomacro = open("FileHeader.h","w")
inputtomacro.write("Char_t *cFileNames[][50] = {");


list_dir = []
index_dir = []
list_numfiles = []
for dirname in os.listdir('.'):
    if os.path.isdir(dirname):
         print ("Working in ------------------------------ "+ dirname)
         temp_files = []
         for files in os.listdir(dirname):
             if files.find(".root") > -1:
                 temp_files.append(files)

         if len(temp_files) > 0:
             list_numfiles.append(len(temp_files))
             inputtomacro.write("{")
             for files in temp_files:
                 inputtomacro.write("\""+files+"\",")
             for i in range(len(temp_files),49):
                 inputtomacro.write("\"-\",")
             inputtomacro.write("\"-\"")
             inputtomacro.write("},"+"\n")
             list_dir.append(dirname)
             index = 0
             for tempdir in mcDir:
                 if tempdir.find(dirname) > -1 :
                     print str(index) + " : " +tempdir + " : " + dirname
                     break
                 index = index + 1
             if index < len(mcDir):
                 index_dir.append(index)
             else:
                 index_dir.append(-1)
inputtomacro.write("};\n");



for indicies in index_dir:
    if indicies > -1:
        print mcStr[indicies] + " : " + mcDir[indicies] + " : " + mcTyp[indicies] + " : " +str(xSecArr[indicies])
        

inputtomacro.write("Char_t *cFileStr[] = {");
tempindex = 0
for indicies in index_dir:
    if indicies > -1:
        inputtomacro.write("\""+mcStr[indicies]+"\"")
    else :
        inputtomacro.write("\" - \"")
    if tempindex < len(index_dir) - 1:
        inputtomacro.write(",")
    tempindex =  tempindex + 1
    
inputtomacro.write("};\n\n");



inputtomacro.write("Char_t *cFileDir[] = {");
tempindex = 0
for indicies in index_dir:
    if indicies > -1:
        inputtomacro.write("\""+mcDir[indicies]+"\"")
    else :
        inputtomacro.write("\" - \"")
    if tempindex < len(index_dir) - 1:
        inputtomacro.write(",")
    tempindex =  tempindex + 1
    
inputtomacro.write("};\n\n");




inputtomacro.write("Int_t iFileType[] = {");
tempindex = 0
for indicies in index_dir:
    if indicies > -1:
        if mcTyp[indicies].find("Signal") > -1:
            inputtomacro.write("0")
        else:
            if mcTyp[indicies].find("QCD") > -1:
                inputtomacro.write("2")
            else :
                  inputtomacro.write("1")
    else :
        inputtomacro.write("-1")
    if tempindex < len(index_dir) - 1:
        inputtomacro.write(",")
    tempindex =  tempindex + 1
    
inputtomacro.write("};\n\n");




inputtomacro.write("Double_t  dCrossSection[] = {");
tempindex = 0
for indicies in index_dir:
    if indicies > -1:
        inputtomacro.write(str(xSecArr[indicies]))
    else :
        inputtomacro.write("-1")
    if tempindex < len(index_dir) - 1:
        inputtomacro.write(",")
    tempindex =  tempindex + 1
    
inputtomacro.write("};\n\n");

inputtomacro.write("Int_t iFileLen[] = {");
tempindex = 0
for indicies in index_dir:
    print str(list_numfiles[tempindex])
    inputtomacro.write(str(list_numfiles[tempindex]))
    if tempindex < len(index_dir) - 1:
        inputtomacro.write(",")
    tempindex =  tempindex + 1
inputtomacro.write("};\n\n")


inputtomacro.write("Int_t iSamples="+str(len(index_dir))+";\n\n")

