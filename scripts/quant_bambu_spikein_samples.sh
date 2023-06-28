samples_drna=("SGNex_Hct116_directRNA_replicate3_run1" "SGNex_K562_directRNA_replicate4_run1" "SGNex_K562_directRNA_replicate5_run1" "SGNex_MCF7_directRNA_replicate4_run1")
samples_cdna=("SGNex_Hct116_cDNA_replicate3_run3" "SGNex_Hct116_cDNA_replicate4_run2" "SGNex_Hct116_directcDNA_replicate3_run2" "SGNex_Hct116_directcDNA_replicate4_run1" "SGNex_Hct116_directcDNA_replicate5_run1" "SGNex_K562_directcDNA_replicate4_run2" "SGNex_MCF7_cDNAStranded_replicate3_run2" "SGNex_MCF7_directcDNA_replicate3_run3" "SGNex_MCF7_directcDNA_replicate4_run2")

## using GNU parallel
#parallel --progress -j 16 '../target/release/oarfish -a ../../data/bambu_paper/aln/{}.bam -o ../../data/bambu_paper/quants/oarfish/{}.tsv' ::: ${samples_cdna[@]}
#parallel --progress -j 16 '../target/release/oarfish -t 50 -a ../../data/bambu_paper/aln/{}.bam -o ../../data/bambu_paper/quants/oarfish/{}.tsv' ::: ${samples_drna[@]}

## using rush (Go program)
RUSH=~/go/bin/rush
echo ${samples_drna[@]} | $RUSH -j 16 --eta -T "r" -D " " "../target/release/oarfish -t 50 -a ../../data/bambu_paper/aln/{}.bam -o ../../data/bambu_paper/quants/oarfish/{}.tsv"
echo ${samples_cdna[@]} | $RUSH -j 16 --eta -T "r" -D " " "../target/release/oarfish -n -a ../../data/bambu_paper/aln/{}.bam -o ../../data/bambu_paper/quants/oarfish/{}.tsv"

