from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from Bio.Seq import Seq
HTTP = HTTPRemoteProvider()

configfile: "config.json"

def get_samples():
    samples,= glob_wildcards(config['FASTQ-path'] + '/{ID}_1.fq.gz')
    return [sample for sample in samples if not sample.startswith('._')]

def get_samples_from_BAM():
    samples,= glob_wildcards(config['bam-path'] + '/{ID}.bam')
    return [sample for sample in samples if not sample.startswith('._')]

rule list_samples:
    run:
        for sample in get_samples():
            print(sample)

rule download_genome:
    input:
        HTTP.remote(config['SC5314-genome-url'],keep_local=False)
    output:
        config['SC5314-genome-path'] + '/C_albicans_SC5314.fasta'
    shell:
        'gzcat {input} > {output}'

rule download_GFF:
    input:
        HTTP.remote(config['SC5314-GFF-url'],keep_local=False)
    output:
        config['SC5314-GFF-path'] + '/C_albicans_SC5314_features.gff'
    shell:
        'cp {input} {output}'

rule remove_chrom_features:
    input:
        config['SC5314-GFF-path'] + '/C_albicans_SC5314_features.gff'
    output:
        config['SC5314-GFF-path'] + '/C_albicans_SC5314_features_no_chroms.gff'
    shell:
        'cat {input} | grep -v "C_albicans_SC5314	CGD	chromosome	1" > {output}'

rule split_genome:
    input:
        config['SC5314-genome-path'] + '/C_albicans_SC5314.fasta'
    output:
        [config['SC5314-genome-path'] + '/C_albicans_SC5314.fasta.split/C_albicans_SC5314.id_' + chrom + '.fasta' for chrom in config['chromosomes']]
    shell:
        "seqkit split -i {input}"

rule haplotype_A_genome:
    input:
        [config['SC5314-genome-path'] + '/C_albicans_SC5314.fasta.split/C_albicans_SC5314.id_' + chrom + '.fasta' for chrom in config['chromosomes_A']]
    output:
        config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta'
    shell:
        "cat {input} > {output}"

rule haplotype_B_genome:
    input:
        [config['SC5314-genome-path'] + '/C_albicans_SC5314.fasta.split/C_albicans_SC5314.id_' + chrom + '.fasta' for chrom in config['chromosomes_B']]
    output:
        config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_B.fasta'
    shell:
        "cat {input} > {output}"

rule bwa_index_haplotype_A:
    input:
        config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta'
    output:
        [config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta' + ending for ending in ['.amb','.ann','.bwt','.pac','.sa']]
    shell:
        "bwa index {input}"

rule bwa_single:
    input:
        ref=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta',
        index=[config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta' + ending for ending in ['.amb','.ann','.bwt','.pac','.sa']],
        forward=config['FASTQ-path'] + '/{sample}_1.fq.gz',
        reverse=config['FASTQ-path'] + '/{sample}_2.fq.gz'
    threads: 8
    output:
        temp(config['sam-path'] + '/{sample}.sam')
    shell:
        "bwa mem -t {threads} {input.ref} {input.forward} {input.reverse} | samtools view -h -F 4 > {output}"

rule bam_single:
    input:
        config['sam-path'] + '/{sample}.sam'
    threads: 1
    output:
        config['bam-path'] + '/{sample}.bam'
    shell:
        "samtools sort {input} > {output}"

rule bam:
    input:
        [config['bam-path'] + '/' + sample + '.bam' for sample in get_samples()]

rule index_bam_single:
    input:
        config['bam-path'] + '/{sample}.bam'
    output:
        config['bam-path'] + '/{sample}.bam.bai'
    shell:
        "samtools index {input}"

rule index_bam:
    input:
        [config['bam-path'] + '/' + sample + '.bam.bai' for sample in get_samples_from_BAM()]

