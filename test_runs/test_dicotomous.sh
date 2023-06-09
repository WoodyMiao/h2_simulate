set -e

echo -e "\n\n1.hapsim_and_realize_binary"
sleep 1
../binary_traits/1.hapsim_and_realize_binary.py --nt 3 --maf-min 0.05 --n-gwa 1000 --n-ld 100 --n-rep 2 --n-pop-beta 1e5 \
    --h2g 0.1 --pqtl-vals 0.1 --neg-alpha 1.0 --prevalence-vals 0.1 \
    --vcf-ref 1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.head10k.gz \
    --hm3-snp hapmap3.snp --gene-bed chr1.bed.head3 --out-dir /app/user/ml/h2_simulate/test_runs/test_dicotomous


echo -e "\n\n3.makefile_step1"
sleep 1
../3.makefile_step1.py --n-ld 100 --n-rep 2 --out-dir /app/user/ml/h2_simulate/test_runs/test_dicotomous --skip-gcta --skip-assoc-ld

for method in lder ldak_sumher kggsee ldsc; do
    make -j 6 -f test_dicotomous/makefile.step1.$method
done


for assoc in logit linear; do
    echo -e "\n\n4.makefile_step2 $assoc"
    sleep 1
    ../binary_traits/4.makefile_step2_binary.py --chrom 1 --n-gwa 1000 --n-ld 100 --n-rep 2 --h2g 0.1 --pqtl-vals 0.1 --neg-alpha 1.0 --prevalence-vals 0.1 \
        --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess --jobs 6 --assoc-test $assoc --out-dir /app/user/ml/h2_simulate/test_runs/test_dicotomous

    for method in lder ldak_sumher hess kggsee ldak_gbat ldsc; do
        echo -e "\n\nmakefile.step2.$method"
        sleep 1
        sh test_dicotomous/makefile.step2.$assoc.$method.sh
    done

    echo -e "\n\n5.harvest_outputs $assoc"
    sleep 1
    ../5.harvest_outputs.py --nt 3 --binary-assoc-test $assoc --n-ld 100 --n-rep 2 --h2g-vals 0.1 --prevalence-vals 0.1 \
        --pqtl-vals 0.1 --neg-alpha-vals 1.0 --out-dir /app/user/ml/h2_simulate/test_runs/test_dicotomous \
        --skip-assoc-ld --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess --out-suffix all.tsv
done


echo -e "\n\nAll done succesfully.\n"
