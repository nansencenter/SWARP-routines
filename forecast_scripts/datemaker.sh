# Creates a datelist for testing out process_FCresults.sh

datelist=datelist.txt
touch $datelist
echo $(date +%Y%m%d) >> $datelist
echo $(date +%Y-%m-%d) >> $datelist
echo $(date +%Y) >> $datelist
echo $(date +%m) >> $datelist
echo $(date +%d) >> $datelist
echo $(date +%j) >> $datelist
