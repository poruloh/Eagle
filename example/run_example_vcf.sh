../eagle \
    --vcf=EUR_test.vcf.gz \
    --geneticMapFile=../tables/genetic_map_hg19_withX.txt.gz \
    --chrom=21 \
    --outPrefix=phased \
    --numThreads=4 \
    2>&1 | tee example_vcf.log

### run eagle without any parameters to list options

### typical options for phasing without a reference in VCF/BCF mode:
# to import genetic map coordinates: --geneticMapFile=tables/genetic_map_hg##.txt.gz
# to select a region to phase: --bpStart, --bpEnd

### old:
# to use Eagle1 algorithm: --v1
# to use Eagle1 fast mode: --v1fast
