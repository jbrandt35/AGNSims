
BBH_separation=$1
perturber_separation=$2
i=$3

mkdir -p /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}

cp /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/utils/code/* /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}
cp /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/utils/submit.pbs /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}

sed -i "s|NAME|B_${BBH_separation}_p_${perturber_separation}_${i}|g" /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.pbs
sed -i "s|WORKDIR|/storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_${BBH_separation}/perturber_separation_${perturber_separation}/${i}|g" /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.pbs

sed -i "s|PERTSEPARATION|${perturber_separation}|g" /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/Simulation.py
sed -i "s|BINSEPARATION|${BBH_separation}|g" /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/Simulation.py

qsub /storage/home/hhive1/jbrandt35/data/AGN_sims/heatmap/runs/BBH_separation_$BBH_separation/perturber_separation_$perturber_separation/${i}/submit.pbs

