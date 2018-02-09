#!/bin/bash

MESH_SIZE=(510 632 1600 6840)
COMPILER_PARAMS=(\  gpu)
MODE_STR=(cpu gpu)

NUM_WARMUPS=1
NUM_EXECUTIONS=3

NUM_MODES=${#COMPILER_PARAMS[@]}
NUM_MESHES=${#MESH_SIZE[@]}

for ((i=0; i<${NUM_MODES}; i++)); do
	compiler_param=${COMPILER_PARAMS[$i]}
	mode_str=${MODE_STR[$i]}
	extra_flag=${EXTRA_FLAGS[$i]}

	make clean
	FPREC=double make ${compiler_param}

	for thread in {4,4}; do
		export OMP_NUM_THREADS=${thread}

		for ((k=0;k<${MESH_SIZE};k++)); do
			echo $k
			mesh=${MESH_SIZE[$k]}
			file_sta=E${mesh}e.dat
			file_dyn=E${mesh}d.dat
			sol_sta=S${mesh}e_${MODE_STR[$i]}_${thread}.dat 
			sol_dyn=S${mesh}d_${MODE_STR[$i]}_${thread}.dat
		
			echo "Executando para $mesh em modo $mode_str com $thread threads"
			echo ""
			for ((l=1;l<=NUM_WARMUPS;l++)); do
				echo "Warmup ${l}: "
				./main $file_sta $file_dyn $sol_sta $sol_dyn
			done
			for ((l=1;l<=NUM_EXECUTIONS;l++)); do
				echo "Execução número $l"
				output_file="results/results_${mode_str}_${mesh}_${thread}_${l}.txt"
				./main $file_sta $file_dyn $sol_sta $sol_dyn $extra_flag > ${output_file}
			done
		done
	done
done

