#PBS -N d.n2e4
#PBS -l nodes=pan03:ppn=64
#PBS -l mem=500g
#PBS -o /dev/null
#PBS -e /dev/null

set -e

cd /home/ml/work/q.EHE_paper/h.hapsim_h2

export MKL_NUM_THREADS=10
export NUMEXPR_NUM_THREADS=10
export OMP_NUM_THREADS=10

# 1.hapsim_a_population
/app/user/ml/h2_simulate/1.hapsim_a_population.py --nt 10 \
    --maf-min 0.01 --n-gwa 2e4 --n-ld 2e3,5e2 --n-rep 100 \
    --hm3-snp /app/resources/HapMap3_b36/hapmap3_all_pop_union.snp \
    --gene-bed /home/ml/work/q.EHE_paper/h.hapsim_h2/a.coding_gene_bed/chr1.bed \
    --vcf-ref /home/ml/resources/1kg_eur_hg19/a.original_vcf_bychr/1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.gz \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1


# 2.realize_beta_and_gwa
/app/user/ml/h2_simulate/2.realize_beta_and_gwa.py --nt 10 \
    --gene-list /home/ml/work/q.EHE_paper/h.hapsim_h2/a.coding_gene_bed/chr1.one4th.lst --n-gwa 2e4 --n-rep 100 \
    --h2g-vals 1e-2,1e-3,1e-4 --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25,1.0 --prevalence-vals 0.1,0.01 \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1


# 3.makefile_step1
/app/user/ml/h2_simulate/3.makefile_step1.py --skip-gcta \
    --gene-list /home/ml/work/q.EHE_paper/h.hapsim_h2/a.coding_gene_bed/chr1.one4th.lst --n-ld 2e3,5e2 --n-rep 100 \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1

for method in lder ldak_sumher kggsee ldsc; do
    make -j 64 -f /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1/makefile.step1.$method
done


# 4.makefile_step2
/app/user/ml/h2_simulate/4.makefile_step2.py \
    --chrom 1 --n-gwa 2e4 --n-ld 2e3,5e2 --n-rep 100 --h2g-vals 1e-2,1e-3,1e-4 \
    --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25,1.0 --maf-min 0.01 \
    --gene-list /home/ml/work/q.EHE_paper/h.hapsim_h2/a.coding_gene_bed/chr1.one4th.lst \
    --jobs 64 --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1

for method in lder ldak_sumher hess kggsee ldak_gbat ldsc; do
    sh /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1/makefile.step2.$method.sh
done


# 5.harvest_outputs
for q in q0-q1 q1-q2 q2-q3 q3-q4; do
    /app/user/ml/h2_simulate/5.harvest_outputs.py --nt 20 \
        --n-ld 2e3,5e2 --n-rep 100 --h2g-vals 1e-2,1e-3,1e-4 --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 0.25,1.0 \
        --gene-list /home/ml/work/q.EHE_paper/h.hapsim_h2/a.coding_gene_bed/chr1.one4th_$q.lst \
        --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess --out-suffix one4th_$q.tsv \
        --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1
done


# Run GCTA
# 3.makefile_step1
/app/user/ml/h2_simulate/3.makefile_step1.py --n-ld 2e3,5e2 --n-rep 50 \
    --skip-mkdir --gene-list a.coding_gene_bed/chr1.one32th.lst \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1

make -j 10 -f d.n2e4_2e3_5e2_rep100_chr1/makefile.step1.gcta


# 4.makefile_step2
/app/user/ml/h2_simulate/4.makefile_step2.py \
    --chrom 1 --n-gwa 2e4 --n-ld 2e3,5e2 --n-rep 50 --h2g-vals 1e-3 \
    --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 1.0 --maf-min 0.01 \
    --gene-list a.coding_gene_bed/chr1.one32th.lst --jobs 40 --gcta --skip-mkdir \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1

make -j 40 -f d.n2e4_2e3_5e2_rep100_chr1/makefile.pqtl1_alpha-1.0_h2g1e-3.gcta
make -j 40 -f d.n2e4_2e3_5e2_rep100_chr1/makefile.pqtl0.1_alpha-1.0_h2g1e-3.gcta
make -j 40 -f d.n2e4_2e3_5e2_rep100_chr1/makefile.pqtl0.01_alpha-1.0_h2g1e-3.gcta


# 5.harvest_outputs
/app/user/ml/h2_simulate/5.harvest_outputs.py --nt 20 \
    --n-ld 2e3,5e2 --n-rep 50 --h2g-vals 1e-3 --pqtl-vals 1,0.1,0.01 --neg-alpha-vals 1.0 \
    --gene-list a.coding_gene_bed/chr1.one32th.lst --out-suffix one32th.tsv \
    --gcta --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess \
    --out-dir /home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1 \
