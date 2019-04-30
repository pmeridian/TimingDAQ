#!/bin/bash

run_number=$1
config_version=$2

code_dir=/afs/cern.ch/user/m/meridian/work/MTD/MTDTB_FNAL_201811/TimingDAQ/
data_dir=/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Nov2018/raw/
out_dir=${data_dir}

echo "============= Job for run $1 ===================="
source /cvmfs/sft.cern.ch/lcg/releases/LCG_94/ROOT/6.14.04/x86_64-slc6-gcc62-opt/ROOT-env.sh
cd $code_dir
echo $PWD

ls -ltrh ${data_dir}

logfile=automation/lsfbatch/log/$run_number.log

python automation/DecodeData.py --vVME $config_version -R $run_number --data_dir=$data_dir --code_dir=$code_dir --out_dir=$out_dir -v
