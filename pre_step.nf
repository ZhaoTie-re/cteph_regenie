params.cteph_vcf_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/CTEPH/08.addTommo_HGVD_AF_vcf/'
params.naga_vcf_path = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/NAGAHAMA/08.addTommo_HGVD_AF_vcf/'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_regenie'

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
    tuple val(chr), file(snp_vcf), file(snp_vcf_tbi) into snp_vcf_ch

    script:
    snp_vcf = chr + '.tommo.snp.vcf.gz'
    snp_vcf_tbi = chr + '.tommo.snp.vcf.gz.tbi'
    """
    bcftools view -m2 -M2 -v snps ${merged_vcf} -Oz -o ${snp_vcf} --threads 2
    bcftools index -t ${snp_vcf} --threads 2
    """
}



// process plink_prepare_cteph {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${group}:${chr}"

//     publishDir "${params.outdir}/01.plink_prepare/${group}", mode: 'symlink'

//     input:
//     tuple val(group), val(chr) from cteph_chr_ch
//     val(filePath) from params.cteph_vcf_path

//     output:
//     tuple val(group), val(chr), file("*.bed"), file("*.bim"), file("*.fam") into bed_cteph_ch

//     script:
//     vcf = filePath + chr + '.tommo.vcf.gz'
//     bed_prefix = group + '.chr' + chr
//     """
//     export PATH=/home/b/b37974/plink:\$PATH
//     plink --vcf ${vcf} --make-bed --out ${bed_prefix}
//     """
// }

// process plink_prepare_naga {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${group}:${chr}"

//     publishDir "${params.outdir}/01.plink_prepare/${group}", mode: 'symlink'

//     input:
//     tuple val(group), val(chr) from naga_chr_ch
//     val(filePath) from params.naga_vcf_path

//     output:
//     tuple val(group), val(chr), file("*.bed"), file("*.bim"), file("*.fam") into bed_naga_ch

//     script:
//     vcf = filePath + chr + '.tommo.vcf.gz'
//     bed_prefix = group + '.chr' + chr
//     """
//     export PATH=/home/b/b37974/plink:\$PATH
//     plink --vcf ${vcf} --make-bed --out ${bed_prefix}
//     """
// }


