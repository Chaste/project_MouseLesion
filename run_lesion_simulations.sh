#!/bin/bash

BUILD_TYPE=intel_production
NUM_PROCS=1
TEST=TestFibroblasts

cd ../..
scons cl=1 -j12 b=IntelProduction projects/MouseLesion/test/${TEST}.hpp
cd projects/MouseLesion

# Run main control
mpirun -np ${NUM_PROCS} ./build/${BUILD_TYPE}/${TEST}Runner --boundary-conduction-scaling 1 --end-time 80 --scar-conduction-scaling 1 --scar-membrane-area-scaling 1 true 2>&1 &

declare -a arr=("--use-neutral-cells" "--use-fibroblasts")

## now loop through the above array
for i in "${arr[@]}"
do
	# Run 'controls'
	mpirun -np ${NUM_PROCS} ./build/${BUILD_TYPE}/${TEST}Runner --boundary-conduction-scaling 1 --end-time 80 --scar-conduction-scaling 1 --scar-membrane-area-scaling 1 ${i} true 2>&1 &

	# Immediately after injury

	mpirun -np ${NUM_PROCS} ./build/${BUILD_TYPE}/${TEST}Runner --boundary-conduction-scaling 1 --end-time 80 --scar-conduction-scaling 0.001 --scar-membrane-area-scaling 1 ${i} true 2>&1 &

	# Results for 'capacitance' study

	mpirun -np ${NUM_PROCS} ./build/${BUILD_TYPE}/${TEST}Runner --boundary-conduction-scaling 1 --end-time 80 --scar-conduction-scaling 0.1 --scar-membrane-area-scaling 1 ${i} true 2>&1 &

	mpirun -np ${NUM_PROCS} ./build/${BUILD_TYPE}/${TEST}Runner --boundary-conduction-scaling 1 --end-time 80 --scar-conduction-scaling 0.1 --scar-membrane-area-scaling 0.1 ${i} true 2>&1 &

	mpirun -np ${NUM_PROCS} ./build/${BUILD_TYPE}/${TEST}Runner --boundary-conduction-scaling 1 --end-time 80 --scar-conduction-scaling 0.1 --scar-membrane-area-scaling 0.01 ${i} true 2>&1 &

	mpirun -np ${NUM_PROCS} ./build/${BUILD_TYPE}/${TEST}Runner --boundary-conduction-scaling 1 --end-time 80 --scar-conduction-scaling 0.1 --scar-membrane-area-scaling 0.001 ${i} true 2>&1 &

	mpirun -np ${NUM_PROCS} ./build/${BUILD_TYPE}/${TEST}Runner --boundary-conduction-scaling 1 --end-time 80 --scar-conduction-scaling 0.1 --scar-membrane-area-scaling 0.0001 ${i} true 2>&1 &

done
