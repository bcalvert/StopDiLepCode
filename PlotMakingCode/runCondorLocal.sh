#! /bin/bash

# Usage:  runCondor.sh 'command to run'

DATE=`date +%s%N | sed -e 's/0*$//'`
NAME=`echo $* | tr -s ' ' '_' | tr -s ':' '_' | tr -s '/' '_' | cut -c 1-24`
SCRIPTNAME=${NAME}_${DATE}.sh
SUBNAME=${NAME}_${DATE}.submit

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
+IsLocalJob=true
requirements=TARGET.IsLocalSlot
queue
EOF

condor_submit $SUBNAME
#requirements = regexp("R510-0-4*",machine,"i") != True
