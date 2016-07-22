echo 'If pulling from github, change --geneticMapFile to ../tables/genetic_map_hg19_example.txt.gz'
echo

../eagle \
    --vcfRef=ref.bcf \
    --vcfTarget=target.vcf.gz \
    --geneticMapFile=../tables/genetic_map_hg19_withX.txt.gz \
    --outPrefix=target.phased \
    2>&1 | tee example_ref.log
