qstat -u jbrandt35 > temp.txt
sed -i '1,4d' temp.txt
number_running=$(grep R -o temp.txt | wc -l)
number_queued=$(grep Q -o temp.txt | wc -l)
number=$(($number_running + $number_queued))
echo $number
rm temp.txt