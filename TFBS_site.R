Hi Rajesh,
 
This is my code, based on a couple of websites/message boards I found:
eg starting with bedGraph data from GEO
 
# Filter out zeros
awk '$4 > 0' bedGraph_file > filtered_bedGraph_file
 
# sort
sort -k1,1 -k2,2n > sorted_bedGraph
 
# Convert to BED
awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' sorted_bedGraph > bed_file
 
# Call peaks 
macs2 callpeak -t bed_file -f BED -n output_prefix -g hs --outdir output_dir --call-summits --qvalue 0.05 --extsize 200
 
As for filtering for specific genes, as you know it’s not easy, especially for distal enhancers! nearest gene is a safe bet.

