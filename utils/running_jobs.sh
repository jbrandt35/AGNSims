squeue -u jbrandt35 --noheader > temp.txt
number_running=$(grep R -o temp.txt | wc -l)
number_queued=$(grep PD -o temp.txt | wc -l)
number=$(($number_running + $number_queued))
rm temp.txt
echo $number