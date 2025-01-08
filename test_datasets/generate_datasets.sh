#!/bin/bash

/home/msayos/Documents/smarco/WFA2-lib/bin/generate_dataset -o ./v3/simulated/sim_100b_1pc_100000s.seq -n 100000 -l 100 -e 1
/home/msayos/Documents/smarco/WFA2-lib/bin/generate_dataset -o ./v3/simulated/sim_100b_10pc_100000s.seq -n 100000 -l 100 -e 10
/home/msayos/Documents/smarco/WFA2-lib/bin/generate_dataset -o ./v3/simulated/sim_1000b_1pc_1000s.seq -n 1000 -l 1000 -e 1
/home/msayos/Documents/smarco/WFA2-lib/bin/generate_dataset -o ./v3/simulated/sim_1000b_10pc_1000s.seq -n 1000 -l 1000 -e 10
/home/msayos/Documents/smarco/WFA2-lib/bin/generate_dataset -o ./v3/simulated/sim_10000b_1pc_10s.seq -n 10 -l 10000 -e 1
/home/msayos/Documents/smarco/WFA2-lib/bin/generate_dataset -o ./v3/simulated/sim_10000b_10pc_10s.seq -n 10 -l 10000 -e 10
