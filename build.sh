#!/bin/bash
MADPATH=/Users/antonk/src/madness

mpic++ -std=c++11 -g3 -O0 ./fpatom.cc  -I/Users/antonk/src/madness/src $MADPATH/src/apps/chem/libMADchem.a $MADPATH/src/madness/mra/libMADmra.a \
$MADPATH/src/madness/tensor/libMADlinalg.a $MADPATH/src/madness/tensor/libMADtensor.a $MADPATH/src/madness/misc/libMADmisc.a \
$MADPATH/src/madness/muParser/libMADmuparser.a $MADPATH/src/madness/tinyxml/libMADtinyxml.a $MADPATH/src/madness/world/libMADworld.a  -framework Accelerate -lxc
