#!/bin/bash
run_number=$1
config_version=$2

dir=/afs/cern.ch/user/m/meridian/work/MTD/MTDTB_FNAL_201811/TimingDAQ/automation/lsfbatch

cd $dir

bsub -q cmscaf1nh -o ${dir}/log/${run_number}.log -e ${dir}/log/${run_number}.err -J decode_$run_number ${dir}/job.sh $run_number $config_version 

cd -
echo
echo '======= Done =========='
