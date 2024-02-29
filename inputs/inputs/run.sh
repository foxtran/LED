#!/bin/bash

export OMP_NUM_THREADS=32
export MRCC_BASIS_PATH=/home/lokalgi/lokalkorr/ger/MRCC.pipeline/install/share/MRCC/BASIS/
export PATH=/home/lokalgi/lokalkorr/ger/MRCC.pipeline/install/bin:$PATH

cd HF
dmrcc | tee ../HF.log
cp MOCOEF ../CC
cp MOCOEF ../LED
cd ../CC
dmrcc | tee ../CC.log
cp CFFILE.tadump ../LED
cp CFFILE.tadump ../
cd ../LED
dmrcc | tee ../LED.log
dftrun | tee ../LED.log
cp DFTGRID.tadump ../
cd ..
