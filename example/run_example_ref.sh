../eagle \
    --vcfRef=ref.bcf \
    --vcfTarget=target.vcf.gz \
    --geneticMapFile=../tables/genetic_map_hg19_example.txt.gz \
    --outPrefix=target.phased \
    2>&1 | tee example_ref.log
