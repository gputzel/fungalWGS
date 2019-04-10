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

def get_candida_albicans_samples():
    with open('sample_data/candida_albicans_samples.txt') as fi:
        return [l.rstrip('\n') for l in fi.readlines()]

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

rule get_known_snps:
    input:
        HTTP.remote(config['known-snps-url'],keep_local=False)
    output:
        "resources/Candida_SC5314_genome/VCF/A22_Jones_PMID_15123810_Polymorphisms.vcf"
    shell:
        "cp {input} {output}"

rule mark_duplicates_single:
    input:
        config['bam-path'] + '/{sample}.bam'
    output:
        bam=config['bam-markdups-path'] + '/{sample}.bam',
        metrics=config['bam-markdups-path'] + '/mark-duplicates-metrics/{sample}-metrics.txt'
    shell:
        "gatk MarkDuplicates --INPUT {input} --OUTPUT {output.bam} " +
            "--METRICS_FILE {output.metrics} --VALIDATION_STRINGENCY SILENT " +
            "--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 " +
            '--ASSUME_SORT_ORDER "coordinate" '

rule mark_duplicates:
    input:
        [config['bam-markdups-path'] + '/' + sample + '.bam' for sample in get_samples_from_BAM()]

rule reference_index:
    input:
        config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta'
    output:
        config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta.fai'
    shell:
        "samtools faidx {input}"

rule dictionary_reference:
    input:
        config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta'
    output:
        config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.dict'
    shell:
        "gatk CreateSequenceDictionary -R {input}"

rule add_read_groups_single:
    input:
        config['bam-markdups-path'] + '/{sample}.bam'
    output:
        config['bam-readgroups-path'] + '/{sample}.bam'
    shell:
        'gatk AddOrReplaceReadGroups -I {input} -O {output} -LB "{wildcards.sample}" -PL "illumina" -SM "{wildcards.sample}" -PU "{wildcards.sample}"'

rule add_read_groups:
    input:
        [config['bam-readgroups-path'] + '/' + sample + '.bam' for sample in get_samples_from_BAM()]

rule base_recalibrate_single:
    input:
        index=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta.fai',
        dictionary=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.dict',
        reference=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta',
        bam=config['bam-readgroups-path'] + '/{sample}.bam',
        vcf="resources/Candida_SC5314_genome/VCF/A22_Jones_PMID_15123810_Polymorphisms.vcf"
    output:
        bqsr=config['bam-readgroups-path'] + '/bqsr-fits/bqsr-{sample}.txt'
    shell:
        'gatk BaseRecalibrator -R {input.reference} -I {input.bam} ' +
            "--use-original-qualities -O {output.bqsr} " +
            "--known-sites {input.vcf}"

rule base_recalibrate:
    input:
        [config['bam-readgroups-path'] + '/bqsr-fits/bqsr-' + sample + '.txt' for sample in get_samples_from_BAM()]

rule apply_bqsr_single:
    input:
        reference=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta',
        bam=config['bam-readgroups-path'] + '/{sample}.bam',
        bqsr=config['bam-readgroups-path'] + '/bqsr-fits/bqsr-{sample}.txt'
    output:
        bam=config['bam-recalibrated-path'] + '/{sample}.bam'
    shell:
        'gatk ApplyBQSR -R {input.reference} -bqsr {input.bqsr} -I {input.bam} -O {output.bam} ' +
            "--static-quantized-quals 10 --static-quantized-quals 20 " +
            "--static-quantized-quals 30 --add-output-sam-program-record --create-output-bam-md5 --use-original-qualities"

rule apply_bqsr:
    input:
        [config['bam-recalibrated-path'] + '/' + sample + '.bam' for sample in get_samples_from_BAM()]

rule gvcf_single:
    input:
        reference=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta',
        bam=config['bam-recalibrated-path'] + '/{sample}.bam'
    output:
        config["gvcf-path"] + "/{sample}.g.vcf"
    shell:
        'gatk HaplotypeCaller -I {input.bam} -R {input.reference} -O {output} -ERC GVCF'

rule gvcf:
    input:
        [config["gvcf-path"] + "/" + sample + ".g.vcf" for sample in get_samples_from_BAM()]

rule combine_gvcfs:
    input:
        c_albicans_list="sample_data/candida_albicans_samples.txt",
        reference=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta',
        gvcfs=[config["gvcf-path"] + "/" + sample + ".g.vcf" for sample in get_samples_from_BAM()]
    output:
        config["gvcf-combined-path"] + "/combined.g.vcf"
    shell:
        'gatk CombineGVCFs --reference {input.reference} ' +
            " ".join(['-V ' + config["gvcf-path"] + "/" + sample + ".g.vcf " for sample in get_candida_albicans_samples()]) +
            " -O {output}"

rule combined_vcf:
    input:
        reference=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta',
        gvcf=config["gvcf-combined-path"] + "/combined.g.vcf"
    output:
        config["vcf-path"] + "/combined.vcf"
    shell:
        "gatk GenotypeGVCFs -R {input.reference} -V {input.gvcf} -O {output}"

rule zip:
    input:
        config["vcf-path"] + "/combined.vcf"
    output:
        config["vcf-path"] + "/combined.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

rule tabix_index:
    input:
        config["vcf-path"] + "/combined.vcf.gz"
    output:
        config["vcf-path"] + "/combined.vcf.gz.tbi"
    shell:
        "tabix {input}"
