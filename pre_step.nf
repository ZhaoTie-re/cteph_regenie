params.cteph_vcf_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/CTEPH/08.addTommo_HGVD_AF_vcf/'
params.naga_vcf_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/NAGAHAMA/08.addTommo_HGVD_AF_vcf/'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/pre_step'
params.script_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie/scripts'
params.phenotype_lst = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/sample_ls/phenotype.txt'


chromosomes = (1..22)

Channel
    .from(chromosomes)
    .map { chr ->
        def cteph_vcf = "${params.cteph_vcf_path}${chr}.tommo.vcf.gz"
        def naga_vcf = "${params.naga_vcf_path}${chr}.tommo.vcf.gz"
        return [chr, cteph_vcf, naga_vcf]
    }
    .set { vcf_files_ch }

process merge_vcf {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/00.merge_vcf", mode: 'symlink'

    input:
    tuple val(chr), val(cteph_vcf), val(naga_vcf) from vcf_files_ch

    output:
    tuple val(chr), file(merged_vcf), file(merged_vcf_tbi) into vcf_snp_ch

    script:
    merged_vcf = chr + '.tommo.merged.vcf.gz'
    merged_vcf_tbi = chr + '.tommo.merged.vcf.gz.tbi'
    """
    bcftools merge ${cteph_vcf} ${naga_vcf} -Oz -o ${merged_vcf} --threads 2
    bcftools index -t ${merged_vcf} --threads 2
    """
}

process snp_select {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/01.snp_select", mode: 'symlink'

    input:
    tuple val(chr), file(merged_vcf), file(merged_vcf_tbi) from vcf_snp_ch

    output:
    tuple val(chr), file(snp_vcf), file(snp_vcf_tbi) into vcf_info_ch

    script:
    snp_vcf = chr + '.tommo.snp.vcf.gz'
    snp_vcf_tbi = chr + '.tommo.snp.vcf.gz.tbi'
    """
    bcftools view -m2 -M2 -v snps ${merged_vcf} -Oz -o ${snp_vcf} --threads 2
    bcftools index -t ${snp_vcf} --threads 2
    """
}

process info_recalculate {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/02.info_recalculate", mode: 'symlink'

    input:
    tuple val(chr), file(snp_vcf), file(snp_vcf_tbi) from vcf_info_ch

    output:
    tuple val(chr), file(info_vcf), file(info_vcf_tbi) into bed_prepare_ch

    script:
    info_vcf = chr + '.tommo.snp.reinfo.vcf.gz'
    info_vcf_tbi = chr + '.tommo.snp.reinfo.vcf.gz.tbi'
    """
    bcftools +fill-tags ${snp_vcf} -Oz -o ${info_vcf} --threads 2 -- -t 'AF,AC,AN,DP:1=int(sum(FORMAT/DP))'
    bcftools index -t ${info_vcf} --threads 2
    """
}

process bed_prepare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/03.bed_prepare", mode: 'symlink'

    input:
    tuple val(chr), file(info_vcf), file(info_vcf_tbi) from bed_prepare_ch
    path(phenotype) from params.phenotype_lst

    output:
    tuple val(chr), file("${bed_prefix}.bed"), file("${bed_prefix}.bim"), file("${bed_prefix}.fam") into qc_call_rate_ch, qc_maf_ch, qc_hwe_ch, qc_filtering_ch

    script:
    bed_prefix = chr + '.cteph.naga.tommo'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --vcf ${info_vcf} --make-bed --out ${bed_prefix} --pheno ${phenotype}
    """
}

process qc_call_rate {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/04.qc_call_rate", mode: 'symlink'

    input:
    tuple val(chr), file(bed_file), file(bim_file), file(fam_file) from qc_call_rate_ch

    output:
    tuple val(chr), file("*.imiss") into missing_rate_sample_plot
    tuple val(chr), file("*.lmiss") into missing_rate_snp_plot

    script:
    bed_prefix = bed_file.baseName
    out_prefix = chr + '.qc.call_rate'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --missing --out ${out_prefix}
    """
}

process qc_maf {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/05.qc_maf", mode: 'symlink'

    input:
    tuple val(chr), file(bed_file), file(bim_file), file(fam_file) from qc_maf_ch

    output:
    tuple val(chr), file("*.frq") into maf_plot

    script:
    bed_prefix = bed_file.baseName
    out_prefix = chr + '.qc.maf'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --freq --out ${out_prefix}
    """
}

process qc_hwe {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/06.qc_hwe", mode: 'symlink'

    input:
    tuple val(chr), file(bed_file), file(bim_file), file(fam_file) from qc_hwe_ch

    output:
    tuple val(chr), file("*.hwe") into hwe_plot

    script:
    bed_prefix = bed_file.baseName
    out_prefix = chr + '.qc.hwe'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --hardy --out ${out_prefix}
    """
}

process plot_qc_call_rate_sample {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/04.qc_call_rate/missing_rate_sample", mode: 'symlink'

    input:
    tuple val(chr), file(imiss_file) from missing_rate_sample_plot

    output:
    file("*.pdf")

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/qc_plot_call_rate_sample.py --chr ${chr} --sample_miss_path ${imiss_file}
    """
}

process plot_qc_call_rate_snp {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/04.qc_call_rate/missing_rate_snp", mode: 'symlink'

    input:
    tuple val(chr), file(lmiss_file) from missing_rate_snp_plot

    output:
    file("*.pdf")

    script:
    """
    source activate cteph_geno_pro
    python ${params.script_path}/qc_plot_call_rate_snp.py --chr ${chr} --snp_miss_path ${lmiss_file}
    """
}

process qc_filtering {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/07.qc_filtering", mode: 'symlink'

    input:
    tuple val(chr), file(bed_file), file(bim_file), file(fam_file) from qc_filtering_ch

    output:
    file("*.log")

    script:
    bed_prefix = bed_file.baseName
    out_prefix = chr + '.qc.geno.mind.maf.hwe'
    """
    export PATH=/home/b/b37974/:$PATH
    plink --bfile ${bed_prefix} --geno 0.8 --mind 0.8 --maf 0.05 --hwe 1e-6 --indep-pairwise 50 5 0.2 --out ${out_prefix}
    """
}
