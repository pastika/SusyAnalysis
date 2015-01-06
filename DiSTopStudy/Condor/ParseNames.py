import parser
import os
import shutil

f = open('InputFiles.txt_BK')
fout = open('InputFiles.txt_out',"w")
for line in iter(f):
    print line
    line = line.replace("\n","")
    parsedstrings = line.split("/")
    numberofsplis = len(parsedstrings)
    nameofthedir = ""
    for ifiles in range(0,numberofsplis):
        #print parsedstrings[ifiles]
        if parsedstrings[ifiles].find("T2tt_mS")>-1:
            print parsedstrings[ifiles]
            nameofthedir = parsedstrings[ifiles]
    fout.write(line+" "+nameofthedir+"\n")
fout.close()    
f.close()
