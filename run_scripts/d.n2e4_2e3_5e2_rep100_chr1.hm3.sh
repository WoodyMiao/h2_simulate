#PBS -N d.hm3
#PBS -l nodes=pan03:ppn=64
#PBS -l mem=500g
#PBS -o /dev/null
#PBS -e /dev/null

set -e

cd /home/ml/work/q.EHE_paper/h.hapsim_h2

# Use summary statistics calculated from genetic models containing all SNPs but estimate h2 with only HapMap3 SNPs
# 3.makefile_step1
/app/user/ml/h2_simulate/3.makefile_step1.py --extract-hm3 --skip-mkdir \
    --n-ld 2e3,5e2 --n-rep 100 --gene-list a.coding_gene_bed/chr1.one4th.3more_hm3snps.lst \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1

for method in kggsee lder ldak_sumher ldsc; do
    make -j 64 -f d.n2e4_2e3_5e2_rep100_chr1/makefile.hm3_step1.$method
done


# 4.makefile_step2
/app/user/ml/h2_simulate/4.makefile_step2.py --hm3-h2 --skip-mkdir \
    --chrom 1 --n-gwa 2e4 --n-ld 2e3,5e2 --n-rep 100 --h2g-vals 1e-3 \
    --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25 --maf-min 0.01 \
    --gene-list a.coding_gene_bed/chr1.one4th.3more_hm3snps.lst \
    --jobs 64 --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1

for method in lder ldak_sumher hess kggsee ldak_gbat ldsc; do
    sh /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1/makefile.hm3_step2.$method.sh
done


# 5.harvest_outputs
for q in q0-q1 q1-q2 q2-q3 q3-q4; do
    /app/user/ml/h2_simulate/5.harvest_outputs.py --nt 20 --hm3-h2 \
        --n-ld 2e3,5e2 --n-rep 100 --h2g-vals 1e-3 --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25 \
        --gene-list a.coding_gene_bed/chr1.one4th_$q.lst --out-suffix one4th_hm3_$q.tsv \
        --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess \
        --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1
done


# Use only HapMap3 SNPs
export MKL_NUM_THREADS=10
export NUMEXPR_NUM_THREADS=10
export OMP_NUM_THREADS=10

# 2.realize_beta_and_gwa
/app/user/ml/h2_simulate/2.realize_beta_and_gwa.py --nt 10 --hm3-only \
    --gene-list a.coding_gene_bed/chr1.one4th.3more_hm3snps.lst --n-gwa 2e4 --n-rep 100 \
    --h2g-vals 1e-3 --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25 --prevalence-vals 0.1 \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1


# 4.makefile_step2
/app/user/ml/h2_simulate/4.makefile_step2.py --hm3-only \
    --chrom 1 --n-ld 2e3,5e2 --n-gwa 2e4 --n-rep 100 --h2g-vals 1e-3 \
    --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25 --maf-min 0.01 \
    --gene-list a.coding_gene_bed/chr1.one4th.3more_hm3snps.lst \
    --jobs 64 --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1

for method in lder ldak_sumher hess kggsee ldak_gbat ldsc; do
    sh /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1/makefile.step2.$method.sh
done


# 5.harvest_outputs
for q in q0-q1 q1-q2 q2-q3 q3-q4; do
    /app/user/ml/h2_simulate/5.harvest_outputs.py --nt 20 --hm3-only --n-ld 2e3,5e2 \
        --n-rep 100 --h2g-vals 1e-3 --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25 \
        --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess \
        --gene-list a.coding_gene_bed/chr1.one4th_$q.lst \
        --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1 \
        --out-suffix one4th_$q.tsv
done
