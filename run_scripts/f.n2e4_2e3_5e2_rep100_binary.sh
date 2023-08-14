#PBS -N f.binary
#PBS -l nodes=pan03:ppn=64
#PBS -l mem=500g
#PBS -o /dev/null
#PBS -e /dev/null

set -e

cd /home/ml/work/q.EHE_paper/h.hapsim_h2

export MKL_NUM_THREADS=20
export NUMEXPR_NUM_THREADS=20
export OMP_NUM_THREADS=20

# 1.hapsim_and_realize_binary
/app/user/ml/h2_simulate/binary_traits/1.hapsim_and_realize_binary.py --nt 5 --n-gwa 2e4 --n-ld 5e2,2e3 --n-rep 100 \
    --hm3-snp /app/resources/HapMap3_b36/hapmap3_all_pop_union.snp --gene-bed a.coding_gene_bed/chr1.one4th.bed \
    --vcf-ref /home/ml/resources/1kg_eur_hg19/a.original_vcf_bychr/1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.gz \
    --out-dir f.n2e4_2e3_5e2_rep100_binary >f.n2e4_2e3_5e2_rep100_binary.log1 2>f.n2e4_2e3_5e2_rep100_binary.log2

# 1.realize large genes
for gene in gene0184 gene1860 gene0778 gene0718 gene1830; do
    /app/user/ml/h2_simulate/binary_traits/1.hapsim_and_realize_binary_one_gene.py --n-gwa 2e4 --n-ld 5e2,2e3 --n-rep 100 --gene $gene \
        --hm3-snp /app/resources/HapMap3_b36/hapmap3_all_pop_union.snp --gene-bed a.coding_gene_bed/chr1.one4th.bed \
        --vcf-ref /home/ml/resources/1kg_eur_hg19/a.original_vcf_bychr/1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.gz \
        --out-dir f.n2e4_2e3_5e2_rep100_binary 2>>f.n2e4_2e3_5e2_rep100_binary.log2.large2 >>f.n2e4_2e3_5e2_rep100_binary.log1.large2 &
    sleep 10000
done
wait


# 3.makefile_step1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

/app/user/ml/h2_simulate/3.makefile_step1.py --n-ld 5e2,2e3 --n-rep 100 --out-dir f.n2e4_2e3_5e2_rep100_binary \
    --skip-gcta --skip-assoc-ld --gene-list a.coding_gene_bed/chr1.one4th.lst

for method in lder ldak_sumher kggsee ldsc; do
    echo "Making step1.$method" >>f.n2e4_2e3_5e2_rep100_binary.log
    make -kj64 -f f.n2e4_2e3_5e2_rep100_binary/makefile.step1.$method &>f.n2e4_2e3_5e2_rep100_binary/makefile.step1.$method.log
done


# 4.makefile_step2
for assoc in chi2 logit linear; do
    /app/user/ml/h2_simulate/binary_traits/4.makefile_step2_binary.py \
        --assoc-test $assoc --chrom 1 --n-gwa 2e4 --n-ld 5e2,2e3 --n-rep 100 \
        --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess --jobs 64 \
        --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/f.n2e4_2e3_5e2_rep100_binary --gene-list a.coding_gene_bed/chr1.one4th.lst
done


# 5.harvest_outputs
for assoc in logit linear chi2; do
    for q in q0-q1 q1-q2 q2-q3 q3-q4; do
        /app/user/ml/h2_simulate/5.harvest_outputs.py --nt 5 --binary-assoc-test $assoc \
            --n-ld 2e3,5e2 --n-rep 100 --h2g-vals 1e-3 --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25 \
            --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess --skip-assoc-ld \
            --gene-list a.coding_gene_bed/chr1.one4th_$q.lst --out-suffix one4th_$q.tsv \
            --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/f.n2e4_2e3_5e2_rep100_binary &
    done
    wait
done
