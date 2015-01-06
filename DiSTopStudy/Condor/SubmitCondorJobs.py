import parser
import os
import shutil

numberoffilesperjob = 100

f = open('InputFiles.txt')
for line in iter(f):
    print "******************************************************************************************************"
    print line
    linesplit = line.split(" ")
    print linesplit[0]
    print linesplit[1]
    dirname = linesplit[1].replace("\n","");
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    print "In Dir" + os.getcwd()
    #Code creating the list of files
    fileindir =  os.listdir(linesplit[0])
    numberoffiles = len(fileindir)
    print "Number of files in this dir  is :" + str(numberoffiles)
    icurcnt = 1
    ijobcount = 0
    strfilenames = ""
    fulldirname = linesplit[0].replace("/pnfs/","dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/")
    for filename in fileindir:
        if icurcnt > numberoffilesperjob:
            icurcnt = 1
        if icurcnt == 1:
             ijobcount = ijobcount + 1
             strfilenames = "inputfiles_" + str(ijobcount) + ".txt"
             print strfilenames
             ffilenames = open(strfilenames,'w'); 
        print filename + " : " + str(ijobcount)
        
        tempstringname = "       \""+fulldirname+"/"+filename + "\", \n" 
        ffilenames.write(tempstringname)
        icurcnt = icurcnt + 1

    #define a submit file
    submitfilename =  "Submitfile.sh"
    fSubmitfile = open(submitfilename,'w')

    for ifile in range (1,ijobcount+1):
        strtempfilename = "StopAnalysis_job_"+str(ifile)+"_cfg.py"
        fconfig =  open('../StopAnalysis_cfg.py')
        fconfignew =  open(strtempfilename,'w')
        print "Creating : " + strtempfilename
        for configcontent in iter(fconfig):
           # print configcontent
            if configcontent.find('REPLACEWITHFILENAMES') > -1:
                configcontent = ''
                strfilenames = "inputfiles_" + str(ijobcount) + ".txt"
                ffilenames = open(strfilenames); 
                for filenamecontent in iter(ffilenames):
                    fconfignew.write(filenamecontent)
                ffilenames.close()
            if configcontent.find('OUTFILENAME') > -1:
                outfilename = "out_"+dirname+"_"+str(ifile)
                print outfilename
                print configcontent
                configcontent = configcontent.replace('OUTFILENAME',outfilename)
            fconfignew.write(configcontent)
        #close the for loop for config file
        fconfignew.close()    
        fconfig.close()
        strtemprsfilename = "cmsRunScript_"+str(ifile)+".sh"
        frunscript =  open('../cmsRunScript.sh')
        frunscriptnew =  open(strtemprsfilename,'w')
        for runscriptcontent in iter(frunscript):
            if runscriptcontent.find('CONFIGFILE') > -1:
                runscriptcontent = runscriptcontent.replace('CONFIGFILE',strtempfilename)
            if runscriptcontent.find('CURDIR') > -1:
                runscriptcontent = runscriptcontent.replace('CURDIR',os.getcwd())
            frunscriptnew.write(runscriptcontent)
        #close the for loop for run script
        os.system("chmod 777 " + strtemprsfilename)
        frunscriptnew.close()
        frunscript.close()
        strtempcondorfilename = "CondorScript"+str(ifile)
        fcondorscript =  open('../MakePythonScript.sh')
        fcondorscriptnew =  open(strtempcondorfilename,'w')
        for condorscriptcontent in iter(fcondorscript):
            if condorscriptcontent.find('CMSRUNSCRIPT') > -1:
                condorscriptcontent = condorscriptcontent.replace('CMSRUNSCRIPT',strtemprsfilename)
            if condorscriptcontent.find('FILENAME') > -1:
                condorscriptcontent = condorscriptcontent.replace('FILENAME',dirname+"_"+str(ijobcount))
            fcondorscriptnew.write(condorscriptcontent)
        fcondorscriptnew.close()
        fcondorscript.close()
        fSubmitfile.write("condor_submit "+strtempcondorfilename+"\n")
    #Close the for loop that creates the files and jobs
    os.system("chmod 777 " + submitfilename)
    
    fSubmitfile.close()
    os.system("./" + submitfilename)
    os.chdir("../")
    print "In Dir" +os.getcwd()
f.close()
