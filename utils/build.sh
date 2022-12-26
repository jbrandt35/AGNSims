
BBH_separation=$1
perturber_separation=$2
i=$3

mkdir -p ../runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}

cp code/* ../runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}
cp submit.sbatch ../runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}

cd ../runs

sed -i "s|NAME|B${BBH_separation}_p${perturber_separation}_${i}|g" BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.sbatch
sed -i "s|WORKDIR|../runs/BBH_separation_${BBH_separation}/perturber_separation_${perturber_separation}/${i}|g" BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.sbatch

sed -i "s|PERTSEPARATION|${perturber_separation}|g" BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/Simulation.py
sed -i "s|BINSEPARATION|${BBH_separation}|g" BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/Simulation.py

sbatch BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.sbatch

