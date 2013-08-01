#! /bin/bash
for sortFile in `/bin/ls sortOutTree_*_CompSort.txt* | sed "s/sortOut//" | sed "s/_CompSort.txt//"`
  do
    rm -rf ${sortFile}_Comp_UniqMassPts.txt
    echo ${sortFile}
    cat sortOut${sortFile}_CompSort.txt | cut -d\: -f2- | sort -u >> ${sortFile}_Comp_UniqMassPts.txt
done