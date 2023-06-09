set -e

echo -e "\n\n1.hapsim_a_population"
sleep 1
../1.hapsim_a_population.py --nt 3 --maf-min 0.05 --n-gwa 1000 --n-ld 100 --n-rep 2 \
    --vcf-ref 1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.head10k.gz \
    --hm3-snp hapmap3.snp --gene-bed chr1.bed.head3 --out-dir /app/user/ml/h2_simulate/test_runs/test_quantitative

echo -e "\n\n2.realize_beta_and_gwa"
sleep 1
../2.realize_beta_and_gwa.py --nt 3 --h2g-vals 0.1 --pqtl-vals 0.1 \
    --neg-alpha-vals 1.0 --prevalence-vals 0.1 --n-gwa 1000 --n-rep 2 --out-dir /app/user/ml/h2_simulate/test_runs/test_quantitative

echo -e "\n\n3.makefile_step1"
sleep 1
../3.makefile_step1.py --n-ld 100 --n-rep 2 --out-dir /app/user/ml/h2_simulate/test_runs/test_quantitative


for method in kggsee ldak_sumher lder ldsc gcta; do
    echo -e "\n\nmakefile.step1.$method"
    sleep 1
    make -j 6 -f test_quantitative/makefile.step1.$method
done


echo -e "\n\n4.makefile_step2"
sleep 1
../4.makefile_step2.py --chrom 1 --n-gwa 1000 --n-ld 100 --n-rep 2 --h2g-vals 0.1 \
    --pqtl-vals 0.1 --neg-alpha-vals 1.0 --maf-min 0.05 --out-dir /app/user/ml/h2_simulate/test_runs/test_quantitative \
    --kggsee --gcta --ldsc --lder --ldak-gbat --ldak-sumher --hess --jobs 6

for method in kggsee ldak_sumher lder ldsc gcta hess ldak_gbat; do
    echo -e "\n\nmakefile.step2.$method"
    sleep 1
    sh test_quantitative/makefile.step2.$method.sh
done


echo -e "\n\n5.harvest_outputs"
sleep 1
../5.harvest_outputs.py --nt 3 --n-ld 100 --n-rep 2 --h2g-vals 0.1 \
   --pqtl-vals 0.1 --neg-alpha-vals 1.0 --out-dir /app/user/ml/h2_simulate/test_runs/test_quantitative \
   --kggsee --gcta --ldsc --lder --ldak-gbat --ldak-sumher --hess --out-suffix all.tsv


echo -e "\n\nAll done succesfully.\n"
