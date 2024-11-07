params.bed_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_step_all/08.clean_data_step1'
params.covariate_file = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_pca/05.pca_calc_projection/pca_calc_projection.sscore'
params.age_sex_file = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/sample_ls/plink_covar.txt'
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
    file("*.txt") into pheno_file_ch, pheno_file_ch_2, pheno_file_ch_3, pheno_file_ch_4

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/gwas_pheno_prepare.py --covariate_file ${covariate_file}
    """
}

process total_covar_prepare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/02.total_covar_prepare", mode: 'symlink'

    input:
    path(covariate_file) from params.covariate_file
    path(age_sex_file) from params.age_sex_file

    output:
    file("*.txt") into total_covar_file_ch, total_covar_file_ch_2, total_covar_file_ch_3, total_covar_file_ch_4

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/gwas_total_covar_prepare.py --covariate_file ${covariate_file} --age_sex_file ${age_sex_file}
    """
}

process gwas_plink2_firth { 
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/03.gwas_plink2_firth", mode: 'symlink'

    input:
    path(bed_path) from params.bed_path
    file(covariate_file) from total_covar_file_ch
    file(pheno_file) from pheno_file_ch

    output:
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'cteph_gwas_plink2'
    covariate_cols='6-7'
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

process gwas_plink2_glm { 
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/04.gwas_plink2_glm", mode: 'symlink'

    input:
    path(bed_path) from params.bed_path
    file(covariate_file) from total_covar_file_ch_2
    file(pheno_file) from pheno_file_ch_2

    output:
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'cteph_gwas_plink2'
    covariate_cols='6-7'
    col_name="PHENO1"
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 \\
    --bfile ${bed_prefix} \\
    --pheno ${pheno_file} \\
    --pheno-name ${col_name} \\
    --covar ${covariate_file} \\
    --covar-col-nums ${covariate_cols} \\
    --glm hide-covar \\
    --threads 8 \\
    --out ${out_prefix}
    """
}

process gwas_plink2_firth_PCs { 
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/05.gwas_plink2_firth_PCs", mode: 'symlink'

    input:
    path(bed_path) from params.bed_path
    file(covariate_file) from total_covar_file_ch_3
    file(pheno_file) from pheno_file_ch_3

    output:
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'cteph_gwas_plink2'
    covariate_cols='6-12'
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

process gwas_plink2_glm_PCs { 
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/06.gwas_plink2_glm_PCs", mode: 'symlink'

    input:
    path(bed_path) from params.bed_path
    file(covariate_file) from total_covar_file_ch_4
    file(pheno_file) from pheno_file_ch_4

    output:
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'cteph_gwas_plink2'
    covariate_cols='6-12'
    col_name="PHENO1"
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 \\
    --bfile ${bed_prefix} \\
    --pheno ${pheno_file} \\
    --pheno-name ${col_name} \\
    --covar ${covariate_file} \\
    --covar-col-nums ${covariate_cols} \\
    --glm hide-covar \\
    --threads 8 \\
    --out ${out_prefix}
    """
}