#!/bin/bash
path=$(pwd)

nmpi=4
input_path=$path/../condition
initial_POSCAR_file=$path/data/POSCAR_$1

install_dir="/opt/intel/oneapi"
source ${install_dir}/setvars.sh intel64 --force
vasp=~/programs/vasp.5.4.4_vtst_k/bin/vasp_std

cd $path/VASP_exe
if [ -e OUTCAR ];then
    rm ./*
fi

files=(INCAR POTCAR KPOINTS)
cp $initial_POSCAR_file POSCAR
for input in ${files[@]};do
    cp $input_path/$input ./ 
done

OMP_NUM_THREADS=1
mpirun --bind-to none -np ${nmpi} ${vasp} > log
cat OUTCAR | grep 'free  energy   TOTEN  =' | tail -n 1 | awk '{print $(NF -1)}' > $path/data/TOTEN_$1.csv
cd $path
