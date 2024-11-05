params.bed_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_step_all/08.clean_data_step1'
params.covariate_file = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_pca/05.pca_calc_projection/pca_calc_projection.sscore'
params.script_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/scripts'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/gwas_plink2'

process pheno_prepare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/01.pheno_prepare", mode: 'symlink'

    input:
    path(covariate_file) from params.covariate_file

    output:
    file("*.txt") into pheno_file_ch_1

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/gwas_pheno_prepare.py --covariate_file ${covariate_file}
    """
}

process gwas_plink2 { 
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/02.gwas_plink2", mode: 'symlink'

    input:
    path(bed_path) from params.bed_path
    path(covariate_file) from params.covariate_file
    file(pheno_file) from pheno_file_ch_1

    output:
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'cteph_gwas_plink2'
    covariate_cols='6-10'
    col_name="PHENO1"
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 \\
    --bfile ${bed_prefix} \\
    --pheno ${pheno_file} \\
    --pheno-name ${col_name} \\
    --covar ${covariate_file} \\
    --covar-col-nums ${covariate_cols} \\
    --glm hide-covar firth single-prec-cc \\
    --threads 8 \\
    --out ${out_prefix}
    """
}