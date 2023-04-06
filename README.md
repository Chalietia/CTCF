#Code and pipelines I used in this project

## SON C&R data processing

1. For alignment with Bowtie2, I use
`--local --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 1000- x hg19`
2. Reads filtered with:
`samtools view -q 5`
3. Dedup with Picard:
`MarkDuplicates REMOVE_DUPLICATES=True ASSUME_SORT_ORDER=queryname`
4. After generating final bam files, all noise peaks are called with macs3. The Cut&Run noise peaks are very sharp, clearly differs from the broad islands of SON&igg.
`macs3 callpeak -f BAM -g 2.7e9 -q 0.1`
5. The called narrowpeak files are merged with each other then with encode hg19 blacklist file with `bedops -m`
6. bw files generated with
`bamCompare -b1 $SON -b2 $IGG -bs 5000 --smoothLength 15000 --effectiveGenomeSize 2864785220 --exactScaling --normalizeUsing RPKM --scaleFactorsMethod None -bl blacklistPEAK.merge.bed`
7. Speckle associated are called with sicer:
`sicer -t $SON -c $IGG -s hg19 -fdr 0.1 -w 2000 -g 20000`
8. All called SON islands (SPAD) from control samples are intersected with `bedtools interesect`
9. Finally, heatmap matrix calculated by deeptools:
`computeMatrix scale-regions -S $(find . -name '*.bw' -exec echo {} +) -R SPAD.bed -b 500 -a 500 --missingDataAsZero`

