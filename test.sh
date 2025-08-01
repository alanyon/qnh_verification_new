#!/bin/bash -l
CODE_DIR=/home/users/andre.lanyon/qnh_verification_new
export CYCLE_POINT=20240703T0000Z
export MOO_DIR=moose:/adhoc/users/ppdev
export DATA_DIR=/data/users/andre.lanyon/QNH

cd ${CODE_DIR}
module load scitools
python ${CODE_DIR}/qnh_imp_glb_test.py