
BBH_separation=$1
perturber_separation=$2
i=$3

mkdir -p ../runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}

cp code/* ../runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}
cp submit.pbs ../runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}

cd ../runs

sed -i "s|NAME|B_${BBH_separation}_p_${perturber_separation}_${i}|g" BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.pbs
sed -i "s|WORKDIR|/storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_${BBH_separation}/perturber_separation_${perturber_separation}/${i}|g" BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.pbs

sed -i "s|PERTSEPARATION|${perturber_separation}|g" BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/Simulation.py
sed -i "s|BINSEPARATION|${BBH_separation}|g" BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/Simulation.py

qsub BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.pbs

