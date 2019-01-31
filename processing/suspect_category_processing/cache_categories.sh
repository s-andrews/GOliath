# to run the script ./cache_categories.sh suspect_categories/*
# the description is taken from the file name
# the file must contain GO (or other) ids in the firt column - this doesn't check
outfile="../../suspect_GO_categories/suspect_categories.txt"
echo -e "GO_ID\tbias_source" > $outfile
echo -e "writing to file $outfile"

all_files=$@

for file in $all_files; do 
	filename=$(basename $file);
	description=${filename%".txt"}; # remove .txt from filename
	first_col=$(cut -f1 $file)
	
	# if the whole description etc is present e.g. HABITUATION%GOBP%GO:0046959, extract GO id
	if [[ $first_col == *"%"* ]]; then
		cut -f1 $file | while IFS='%' read -ra array; do
			for i in "${array[@]}"; do
				if [[ $i =~ .*GO:.* ]]; then
					# write out to file						
					echo -e "${i}\t${description}" >> $outfile;
				fi  
			done
		done
	else
		cut -f1 $file | while IFS='\n' read -ra array; do			
			for i in "${array[@]}"; do			
				echo -e "${i}\t${description}" >> $outfile;
			done
		done	
	fi
	echo -e "added $description categories to file"
done
