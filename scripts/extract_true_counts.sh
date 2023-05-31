input_reads=$1
output_quants=$2

gunzip -c $input_reads | \
  awk 'NR % 4 == 1 { split($0,arr,"_"); if (arr[3] != "unaligned") { ++count[substr(arr[1], 2)]; } } END { for (key in count) { print key"\t"count[key] } }' > $output_quants
