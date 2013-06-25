#! /bin/bash

# Usage:  runCondor.sh 'command to run'
pwd
DATE=`date +%s%N | sed -e 's/0*$//'`
NAME=`echo $* | tr -s ' ' '_' | tr -s ':' '_' | tr -s '/' '_'`
SCRIPTNAME=${NAME}_${DATE}.sh
SUBNAME=${NAME}_${DATE}.submit
echo $SCRIPTNAME
echo $SUBNAME
cat << EOF >> $SCRIPTNAME
#!/bin/sh
$*
EOF

chmod +x $SCRIPTNAME

cat << EOF >> $SUBNAME
universe = vanilla
getenv = True
executable = $SCRIPTNAME
output = ${SCRIPTNAME}.stdout
error = ${SCRIPTNAME}.stderr
log = ${SCRIPTNAME}.logfile
notification=never
queue
EOF

condor_submit $SUBNAME
#requirements = regexp("R510*",machine,"i") != True


#+IsLocalJob = true
#requirements = TARGET.IsLocalSlot