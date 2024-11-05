params.snp_dir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_step/01.snp_select'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_step_all'
params.script_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/scripts'
params.phenotype_lst = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/sample_ls/plink_pheno.txt'

Channel
    .from(1..22)
    .map { chr ->
        "${params.snp_dir}/${chr}.tommo.snp.vcf.gz"
    }
    .toList()
    .map { paths ->
        paths.join('\n')
    }
    .set { merge_list_content }

process writeMergeList {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/01.vcf_merge", mode: 'symlink'

    input:
    val(merge_list) from merge_list_content

    output:
    file("merge_list.txt") into writeMergeList

    script:
    """
    echo "$merge_list" > merge_list.txt
    """
}

process MergeChromosomes {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/01.vcf_merge", mode: 'symlink'

    input:
    file(merge_list) from writeMergeList

    output:
    tuple file(merged_vcf), file(merged_vcf_tbi) into merged_vcf_ch

    script:
    merged_vcf = 'allchr.tommo.snp.vcf.gz'
    merged_vcf_tbi = 'allchr.tommo.snp.vcf.gz.tbi'
    """
    bcftools concat -f ${merge_list} -Oz -o ${merged_vcf} --threads 8
    bcftools index -t ${merged_vcf} --threads 8
    """
}

process bed_prepare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/02.bed_prepare", mode: 'symlink'

    input:
    tuple file(snp_vcf), file(snp_vcf_tbi) from merged_vcf_ch
    path(phenotype) from params.phenotype_lst

    output:
    tuple file("${bed_prefix}.bed"), file("${bed_prefix}.bim"), file("${bed_prefix}.fam") into qc_call_rate_ch, qc_maf_ch, qc_hwe_ch, qc_filtering_ch, qc_fcoef_ch_1, clean_data_ch

    script:
    bed_prefix = 'allchr.snp'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --vcf ${snp_vcf} --make-bed --out ${bed_prefix} --pheno ${phenotype}
    """
}

process qc_call_rate {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/03.qc_call_rate", mode: 'symlink'

    input:
    tuple file(bed_file), file(bim_file), file(fam_file) from qc_call_rate_ch

    output:
    file("*.imiss") into missing_rate_sample_plot
    file("*.lmiss") into missing_rate_snp_plot

    script:
    bed_prefix = bed_file.baseName
    out_prefix = 'allchr.snp.qc.call_rate'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --missing --out ${out_prefix}
    """
}

process plot_qc_call_rate_sample {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/03.qc_call_rate/missing_rate_sample", mode: 'symlink'

    input:
    file(imiss_file) from missing_rate_sample_plot

    output:
    file("*.pdf")

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/qc_plot_call_rate_sample.py --chr all --sample_miss_path ${imiss_file}
    """
}

process plot_qc_call_rate_snp {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/03.qc_call_rate/missing_rate_snp", mode: 'symlink'

    input:
    file(lmiss_file) from missing_rate_snp_plot

    output:
    file("*.pdf")

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/qc_plot_call_rate_snp.py --chr all --snp_miss_path ${lmiss_file}
    """
}

process qc_maf {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/04.qc_maf", mode: 'symlink'

    input:
    tuple file(bed_file), file(bim_file), file(fam_file) from qc_maf_ch

    output:
    file("*.frq") into maf_plot

    script:
    bed_prefix = bed_file.baseName
    out_prefix = 'allchr.snp.qc.maf'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --freq --out ${out_prefix}
    """
}

process plot_qc_maf {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/04.qc_maf", mode: 'symlink'

    input:
    file(maf_file) from maf_plot

    output:
    file("*.pdf")

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/qc_plot_maf.py --frq_path ${maf_file}
    """
}

process qc_hwe {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/05.qc_hwe", mode: 'symlink'

    input:
    tuple file(bed_file), file(bim_file), file(fam_file) from qc_hwe_ch

    output:
    file("*.hwe") into hwe_plot

    script:
    bed_prefix = bed_file.baseName
    out_prefix = 'allchr.snp..qc.hwe'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --hardy --out ${out_prefix}
    """
}

process plot_qc_hwe {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/05.qc_hwe", mode: 'symlink'

    input:
    file(hwe_file) from hwe_plot

    output:
    file("*.pdf")

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/qc_plot_hwe.py --hwe_path ${hwe_file}
    """
}

process qc_filtering {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/06.qc_filtering", mode: 'symlink'

    input:
    tuple file(bed_file), file(bim_file), file(fam_file) from qc_filtering_ch

    output:
    file("*.prune.in") into qc_fcoef_ch_2
    file("*.irem")
    file("*.log")
    file("*.prune.out")

    script:
    bed_prefix = bed_file.baseName
    out_prefix = 'allchr.snp.qc.filtering'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --maf 0.01 --geno 0.4 --mind 0.4 --hwe 1e-6 --indep-pairwise 50 5 0.2 --out ${out_prefix} --threads 8
    """
}

process qc_fcoef {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/07.qc_fcoef", mode: 'symlink'

    input:
    tuple file(bed_file), file(bim_file), file(fam_file)from qc_fcoef_ch_1
    file(select_snp) from qc_fcoef_ch_2

    output:
    file("*.het") into plot_qc_fcoef

    script:
    bed_prefix = bed_file.baseName
    out_prefix = 'allchr.snp.qc.fcoef'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --extract ${select_snp} --het --out ${out_prefix}
    """
}

process plot_qc_fcoef {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/07.qc_fcoef", mode: 'symlink'

    input:
    file(het_file) from plot_qc_fcoef

    output:
    file("*.pdf")
    file("*.sample") into high_het_sample

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/qc_plot_fcoef.py --chr all --het_path ${het_file}
    """
}

process clean_data_step1 {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outdir}/08.clean_data_step1", mode: 'symlink'

    input:
    tuple file(bed_file), file(bim_file), file(fam_file) from clean_data_ch
    file(high_het) from high_het_sample

    output:
    tuple file("${out_prefix}.bed"), file("${out_prefix}.bim"), file("${out_prefix}.fam") into clean_data_ch_2
    file("*.log")

    script:
    bed_prefix = bed_file.baseName
    out_prefix = 'allchr.snp.qc'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --maf 0.01 --geno 0.4 --mind 0.4 --hwe 1e-6 --remove ${high_het} --keep-allele-order --make-bed --out ${out_prefix}
    """
}






