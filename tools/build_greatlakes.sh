#! /bin/bash

cmake .. -DTASK=exo3 -DCMAKE_PREFIX_PATH="$NC_ROOT;$HOME/opt/" -DCMAKE_CXX_COMPILER="icc" -DCMAKE_C_COMPILER="icc"
