../eagle \
    --bfile=EUR_test \
    --geneticMapFile=USE_BIM \
    --chrom=21 \
    --outPrefix=phased \
    --numThreads=4 \
    2>&1 | tee example.log

### run eagle without any parameters to list options

### typical options for phasing without a reference:
# to import genetic map coordinates: --geneticMapFile=tables/genetic_map_hg##.txt.gz
# to remove indivs or exclude SNPs: --remove, --exclude
# to perform QC on missingness:  --maxMissingPerIndiv, --maxMissingPerSnp
# to select a region to phase: --bpStart, --bpEnd, --bpFlanking

### old:
# to use Eagle1 algorithm: --v1
# to use Eagle1 fast mode: --v1fast
