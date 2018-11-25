#!/bin/bash

LHEINTERFACEDIR="/afs/cern.ch/work/s/sosen/ChongbinTop/SAB";
SABHOME="${LHEINTERFACEDIR}/../common-framework/Framework/TTHbbAnalysis/BucketofTops";

mkdir build;
cd build;
setupATLAS;
asetup Athena,21.0.75,here; #for c++11

g++ -fPIC -shared ${SABHOME}/Root/BucketofTops.cxx -I ${SABHOME}/ -o libBucketofTops.so $(root-config --cflags --libs);

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH;

g++ -Wall -L. -I ${SABHOME}/ ${LHEINTERFACEDIR}/LHEinterface/main.cxx -lBucketofTops $(root-config --cflags --libs);
#g++ -Wall -L ${LHEINTERFACEDIR}/lib -WI,-rpath=${LHEINTERFACEDIR}/lib ${LHEINTERFACEDIR}/LHEinterface/main.cxx -lBucketofTops $(root-config --cflags --libs);

#g++ ${LHEINTERFACEDIR}/LHEinterface/main.cxx -I ${SABHOME}/ -o runbuckets $(root-config --cflags --libs);
#g++ ${LHEINTERFACEDIR}/LHEinterface/main.cxx ${SABHOME}/Root/BucketofTops.cxx -I ${SABHOME}/ -o runbuckets $(root-config --cflags --libs);
#g++ -c ${SABHOME}/Root/BucketofTops.cxx -I ${SABHOME}/ $(root-config --cflags --libs) -fPIC -o BucketofTops.o;
cd ..;
