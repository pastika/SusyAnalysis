#!/bin/sh

for f in  $(cat InputFiles.txt)
do
echo $f
ls -1 $f
done 
