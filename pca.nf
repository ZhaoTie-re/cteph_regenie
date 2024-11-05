params.bed_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_step_all/08.clean_data_step1'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_pca'
params.high_ldfile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/plinkQC_extdata/high-LD-regions-hg38-GRCh38.txt'
params.script_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/scripts'

process hild_set {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/01.hild_set", mode: 'symlink'

    input:
    path(bed_path) from params.bed_path
    path(high_ldfile) from params.high_ldfile

    output:
    file("hild.set") into hild_set_ch
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --make-set ${high_ldfile} --write-set --out hild
    """
}

process ld_prune {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/02.ld_prune", mode: 'symlink'

    input:
    file(hild_set) from hild_set_ch
    path(bed_path) from params.bed_path

    output:
    file("${out_prefix}.prune.in") into ld_prune_out, retain_snps
    file("${out_prefix}.prune.out")
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'rm_hild.ld_prune'
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 --bfile ${bed_prefix} --maf 0.01 --threads 8 --exclude ${hild_set} --indep-pairwise 500 50 0.2 --out ${out_prefix}
    """
}

process king_cutoff {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/03.king_cutoff", mode: 'symlink'

    input:
    file(prune_in) from ld_prune_out
    path(bed_path) from params.bed_path

    output:
    file("*.king.cutoff.in.id") into retain_id_king
    file("*.king.cutoff.out.id")
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'king_cutoff'
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 --bfile ${bed_prefix} --extract ${prune_in} --king-cutoff 0.0884 --threads 8 --out ${out_prefix}
    """
}

process pca_calc {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/04.pca_calc", mode: 'symlink'

    input:
    file(retain_id_king) from retain_id_king
    file(retain_snps) from retain_snps
    path(bed_path) from params.bed_path

    output:
    tuple file("*.acount"), file("*.eigenvec.allele") into pca_calc_projection
    file("*.log")
    file("*.eigenvec")
    file("*.eigenval") into eigenval_ch

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'pca_calc'
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 --bfile ${bed_prefix} --keep ${retain_id_king} --extract ${retain_snps} --freq counts --threads 8 --pca allele-wts 10 --out ${out_prefix}
    """
}

process pca_calc_projection {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/05.pca_calc_projection", mode: 'symlink'

    input:
    tuple file(acount), file(eigenvec_allele) from pca_calc_projection
    path(bed_path) from params.bed_path

    output:
    file("*.sscore") into pca_calc_projection_out
    file("*.log")

    script:
    bed_prefix = bed_path + '/allchr.snp.qc'
    out_prefix = 'pca_calc_projection'
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 --bfile ${bed_prefix} --threads 8 --read-freq ${acount} --score ${eigenvec_allele} 2 6 header-read no-mean-imputation variance-standardize --score-col-nums 7-16 --out ${out_prefix}
    """
}

process plot_pca_calc_projection {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/05.pca_calc_projection", mode: 'symlink'

    input:
    file(sscore) from pca_calc_projection_out
    file(eigenval) from eigenval_ch

    output:
    file("*.html")

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/plot_pca.py --pca_path ${sscore} --eigenval_path ${eigenval}
    """
}




