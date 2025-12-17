### Nanopore short tandem repeats analysis ###
The pipeline focuses on analysing PCR amplified repeat regions <br />

1.) BWA alignment with -K 500M -r 50000,50000 -end-bonus=10000 -no-end-flt parameter (Miller, D. E. et al. Targeted long-read sequencing identifies missing disease-causing variation. The American Journal of Human Genetics 108, 1436–1449 (2021)) <br />
<br />
2.) TRGT(Dolzhenko, E. et al. Characterization and visualization of tandem repeats at genome scale. Nat Biotechnol 42, 1606–1614 (2024)) was used for calling and plotting the repeats. <br />
<br />
3.) The repeats were also annotated using stranger(https://github.com/Clinical-Genomics/stranger.git) <br />

<br />
Note : The bam files generated can be used to run epi2me workflow which uses straglr for calling the repeats.
