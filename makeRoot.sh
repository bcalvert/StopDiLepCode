#! /bin/bash

make_root() {
    cflags=`root-config --cflags`
    libs=`root-config --libs`
#    libs="$libs mt2bisect.cc"
    name=`echo $1 | sed 's,\..*,,g'`
    if [[ $2 ]] ; then name=$2; fi
    g++ -g -Wall -fPIC -pthread -m64 $cflags $libs -o $name $1
}

make_root $1 $2
