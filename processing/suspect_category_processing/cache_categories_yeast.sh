# to run the script ./cache_categories_yeast.sh
# the description is taken from the file name
# the file must contain GO (or other) ids in the first column - this doesn't check
outfile="../../suspect_GO_categories/yeast/suspect_categories.txt"
echo -e "GO_ID\tbias_source" > $outfile
echo -e "writing to file $outfile"

#all_files=$@
all_files="suspect_categories/yeast/*";

for file in $all_files; do 
	filename=$(basename $file);
	description=${filename%".txt"}; # remove .txt from filename
	first_col=$(cut -f1 $file)

	cut -f1 $file | while read -r line; do
		if [[ $line =~ .*%GO:.* ]]; then
			IFS='%' read -ra category_info <<< $line					
			for i in "${category_info[@]}"; do
				if [[ $i =~ .*GO:.* ]]; then
					echo -e "${i}\t${description}" >> $outfile;
				fi
			done
		else
		# some descriptions are split by % but are not GO IDs - we'll just use the whole description as sometimes the IDs
		# can be difficult to extract
		# some descriptions do not contain %
			echo -e "$line\t${description}" >> $outfile;
		fi
	done	
done					
