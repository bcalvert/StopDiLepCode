#! /bin/bash

for indFile in `/bin/ls sortOutTree_*_1.txt | sed "s/sortOut//" | sed "s/_1.txt//"`
  do
  rm -rf sortOut${indFile}_CompSort.txt
  rm -rf sortOut${indFile}_CompUnSort.txt
  echo $indFile
  for file in `/bin/ls sortOut${indFile}_*.txt`
    do
    echo $file
    cat $file >> sortOut${indFile}_CompUnSort.txt
  done
  cat sortOut${indFile}_CompUnSort.txt | sort -s -n -k 1,1 >> sortOut${indFile}_CompSort.txt
#  echo sortOut${indFile}_*.txt
#cat notSortOut${name}.txt | sort -s -n -k 1,1 >> sortOut${name}.txt
done