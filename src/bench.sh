#!/bin/bash

MESH_SIZE=(510 1820 6840 26480)
COMPILER_PARAMS=(\  gpu)
MODE_STR=(cpu gpu)

NUM_MODES=${#COMPILER_PARAMS[@]}
NUM_MESHES=${#MESH_SIZE[@]}

mkdir -p results
rm -f results/*.*

for ((i=0; i<${NUM_MODES}; i++)); do
	compiler_param=${COMPILER_PARAMS[$i]}
	mode_str=${MODE_STR[$i]}
	extra_flag=${EXTRA_FLAGS[$i]}

	make clean
	FPREC=double make ${compiler_param}

	for thread in {1,4}; do
		export OMP_NUM_THREADS=${thread}

		for ((k=0;k<${NUM_MESHES};k++)); do
			mesh=${MESH_SIZE[$k]}
			file_sta=E${mesh}e.dat
			file_dyn=E${mesh}d.dat
			sol_sta=S${mesh}e_${MODE_STR[$i]}_${thread}.dat 
			sol_dyn=S${mesh}d_${MODE_STR[$i]}_${thread}.dat
			
			warmups=3
			executions=30

			if (( $thread == 1 )); then
				warmups=1
				executions=5
			fi

			echo "Executando para $mesh em modo $mode_str com $thread threads"
			echo ""
			for ((l=1;l<=warmups;l++)); do
				echo "Warmup ${l}: "
				./main $file_sta $file_dyn $sol_sta $sol_dyn
			done
			for ((l=1;l<=executions;l++)); do
				echo "Execução número $l"
				output_file="results/results_${mode_str}_${mesh}_${thread}_${l}.txt"
				./main $file_sta $file_dyn $sol_sta $sol_dyn ${extra_flag} >> ${output_file}
			done
		done
	done
done

