import os
import sys
import json
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from Bio.Seq import Seq
HTTP = HTTPRemoteProvider()

configfile: "config.json"

projectName=os.environ.get("PROJECT",None)
if projectName is not None:
    project_config_file = "projects/" + projectName + ".json"
    with open(project_config_file) as json_file:
        project = json.load(json_file)
else:
    print("Please specify a project as a command line environment variable, as in:")
    print('PROJECT="myproject" snakemake ...')
    sys.exit()

def get_samples():
    return [sample for sample in project["samples"].keys()]

rule list_samples:
    run:
        for sample in get_samples():
            print(sample)

rule download_genome:
    output:
        "resources/" + projectName + "/genome.fasta.gz"
    shell:
        'wget ' + project["genome-url"] + ' -O {output}'

rule decompress_genome:
    input:
        "resources/" + projectName + "/genome.fasta.gz"
    output:
        "resources/" + projectName + "/genome.fasta"
    shell:
        "gzcat {input} > {output}"

rule download_GFF:
    output:
        "resources/" + projectName + "/annotation.gff.gz"
    shell:
        'wget ' + project["GFF-url"] + ' -O {output}'

rule decompress_GFF:
    input:
        "resources/" + projectName + "/annotation.gff.gz"
    output:
        "resources/" + projectName + "/annotation.gff"
    shell:
        "gzcat {input} > {output}"

rule bwa_index:
    input:
        "resources/" + projectName + "/genome.fasta"
    output:
        ["resources/" + projectName + "/genome.fasta" + ending for ending in ['.amb','.ann','.bwt','.pac','.sa']]
    shell:
        "bwa index {input}"

rule download_SRA:
    output:
        forward="output/" + projectName + "/FASTQ_SRA/{SRA}_pass_1.fastq.gz",
        reverse="output/" + projectName + "/FASTQ_SRA/{SRA}_pass_2.fastq.gz"
    threads: 4
    shell:
        "fastq-dump --outdir output/" + projectName + "/FASTQ_SRA --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip {wildcards.SRA}"

def merged_FASTQ_forward_input(wildcards):
    sample_data=project["samples"][wildcards.sample]
    if "FASTQ-pairs" in sample_data.keys():
        pairs = sample_data["FASTQ-pairs"]
        return [p["forward"] for p in pairs]
    if "SRA" in sample_data.keys():
        SRA = sample_data["SRA"]
        fname = "output/" + projectName + "/FASTQ_SRA/" + SRA + "_pass_1.fastq.gz"
        return [fname]

rule merged_FASTQ_forward:
    input:
        unpack(merged_FASTQ_forward_input)
    output:
        temp("output/" + projectName + "/FASTQ_merged/{sample}_R1.fastq.gz")
    shell:
        "cat {input} > {output}"

def merged_FASTQ_reverse_input(wildcards):
    sample_data=project["samples"][wildcards.sample]
    if "FASTQ-pairs" in sample_data.keys():
        pairs = sample_data["FASTQ-pairs"]
        return [p["reverse"] for p in pairs]
    if "SRA" in sample_data.keys():
        SRA = sample_data["SRA"]
        fname = "output/" + projectName + "/FASTQ_SRA/" + SRA + "_pass_2.fastq.gz"
        return [fname]


rule merged_FASTQ_reverse:
    input:
        unpack(merged_FASTQ_reverse_input)
    output:
        temp("output/" + projectName + "/FASTQ_merged/{sample}_R2.fastq.gz")
    shell:
        "cat {input} > {output}"

rule bwa_single:
    input:
        ref="resources/" + projectName + "/genome.fasta",
        index=["resources/" + projectName + "/genome.fasta" + ending for ending in ['.amb','.ann','.bwt','.pac','.sa']],
        forward="output/" + projectName + "/FASTQ_merged/{sample}_R1.fastq.gz",
        reverse="output/" + projectName + "/FASTQ_merged/{sample}_R2.fastq.gz"
    threads: 4
    output:
        temp("output/" + projectName + "/SAM/{sample}.sam")
    shell:
        "bwa mem -t {threads} {input.ref} {input.forward} {input.reverse} | samtools view -h -F 4 > {output}"

rule bam_single:
    input:
        "output/" + projectName + "/SAM/{sample}.sam"
    threads: 4
    output:
        "output/" + projectName + "/BAM/{sample}.bam"
    shell:
        "samtools sort -@ 3 -m 4G {input} > {output}"

rule index_bam_single:
    input:
        "output/" + projectName + "/BAM/{sample}.bam"
    output:
        "output/" + projectName + "/BAM/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule index_bam_all:
    input:
        ["output/" + projectName + "/BAM/" + sample + ".bam.bai" for sample in get_samples()]

rule mark_duplicates_single:
    conda:
        "envs/gatk4.yml"
    input:
        "output/" + projectName + "/BAM/{sample}.bam"
    output:
        bam="output/" + projectName + "/BAM-markdups/{sample}.bam",
        metrics="output/" + projectName + "/BAM-markdups/{sample}-metrics.txt"
    run:
        batch = project['samples'][wildcards.sample]['batch']
        optical_distance=str(project['batches'][batch]['optical-duplicate-distance'])
        cmd = "gatk MarkDuplicates --INPUT {input} --OUTPUT {output.bam} --METRICS_FILE {output.metrics} --VALIDATION_STRINGENCY SILENT " + "--OPTICAL_DUPLICATE_PIXEL_DISTANCE " + optical_distance + " " + '--ASSUME_SORT_ORDER "coordinate" '
        shell(cmd)

rule reference_index:
    input:
        "resources/" + projectName + "/genome.fasta"
    output:
        "resources/" + projectName + "/genome.fasta.fai"
    shell:
        "samtools faidx {input}"

rule dictionary_reference:
    conda:
        "envs/gatk4.yml"
    input:
        "resources/" + projectName + "/genome.fasta"
    output:
        "resources/" + projectName + "/genome.dict"
    shell:
        "gatk CreateSequenceDictionary -R {input}"

rule add_read_groups_single:
    conda:
        "envs/gatk4.yml"
    input:
        "output/" + projectName + "/BAM-markdups/{sample}.bam"
    output:
        "output/" + projectName + "/BAM-readgroups/{sample}.bam"
    run:
        sample=wildcards.sample
        batch=project['samples'][sample]['batch']
        platform=project['batches'][batch]['platform']
        cmd = "gatk AddOrReplaceReadGroups -I {input} -O {output} -LB " + sample + " -PL " + platform + " -SM " + sample + " -PU " + batch
        shell(cmd)
#'gatk AddOrReplaceReadGroups -I {input} -O {output} -LB "{wildcards.sample}" -PL "illumina" -SM "{wildcards.sample}" -PU "{wildcards.sample}"'

rule index_RG_bam:
    input:
        "output/" + projectName + "/BAM-readgroups/{sample}.bam"
    output:
        "output/" + projectName + "/BAM-readgroups/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule index_RG_bam_all:
    input:
        ["output/" + projectName + "/BAM-readgroups/" + sample + ".bam.bai" for sample in get_samples()]

rule coverage:
    conda:
        "envs/bedtools.yml"
    input:
        "output/" + projectName + "/BAM-readgroups/{sample}.bam"
    output:
        "output/" + projectName + "/coverage/{sample}.txt"
    threads: 4
    shell:
        "genomeCoverageBed -ibam {input} > {output}"

rule coverage_all:
    input:
        ["output/" + projectName + "/coverage/" + sample + ".txt" for sample in get_samples()]

rule base_recalibrate_single:
    conda:
        "envs/gatk4.yml"
    input:
        index="resources/" + projectName + "/genome.fasta.fai",
        dictionary="resources/" + projectName + "/genome.dict",
        reference="resources/" + projectName + "/genome.fasta",
        bam="output/" + projectName + "/BAM-readgroups/{sample}.bam",
        vcf="resources/Candida_SC5314_genome/VCF/A22_Jones_PMID_15123810_Polymorphisms.vcf"
    output:
        bqsr=config['bam-readgroups-path'] + '/bqsr-fits/bqsr-{sample}.txt'
    shell:
        'gatk BaseRecalibrator -R {input.reference} -I {input.bam} ' +
            "--use-original-qualities -O {output.bqsr} " +
            "--known-sites {input.vcf}"

rule gvcf_single:
    conda:
        "envs/gatk4.yml"
    input:
        reference="resources/" + projectName + "/genome.fasta",
        index="resources/" + projectName + "/genome.fasta.fai",
        dictionary="resources/" + projectName + "/genome.dict",
        bam="output/" + projectName + "/BAM-readgroups/{sample}.bam",
        bam_index="output/" + projectName + "/BAM-readgroups/{sample}.bam.bai"
    output:
        #config["gvcf-path"] + "/{sample}.g.vcf"
        "output/" + projectName + "/GVCF/{sample}.g.vcf"
    shell:
        'gatk HaplotypeCaller -I {input.bam} -R {input.reference} -O {output} -ERC GVCF'

rule combine_gvcfs:
    conda:
        "envs/gatk4.yml"
    input:
        reference="resources/" + projectName + "/genome.fasta",
        gvcfs=["output/" + projectName + "/GVCF/" + sample + ".g.vcf" for sample in get_samples() if not sample in project["exclude_from_VCF"]]
    output:
        "output/" + projectName + "/GVCF_combined/combined.g.vcf"
    shell:
        'gatk CombineGVCFs --reference {input.reference} ' +
            " ".join(['-V ' + 'output/' + projectName + "/GVCF/" + sample + ".g.vcf " for sample in get_samples() if not sample in project["exclude_from_VCF"]]) +
            " -O {output}"

rule combined_vcf:
    conda:
        "envs/gatk4.yml"
    input:
        reference="resources/" + projectName + "/genome.fasta",
        gvcf="output/" + projectName + "/GVCF_combined/combined.g.vcf"
    output:
        "output/" + projectName + "/VCF/combined.vcf"
    shell:
        "gatk GenotypeGVCFs -R {input.reference} -V {input.gvcf} -O {output}"

rule zip:
    input:
        "output/" + projectName + "/VCF/combined.vcf"
    output:
        "output/" + projectName + "/VCF/combined.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

rule tabix_index:
    input:
        "output/" + projectName + "/VCF/combined.vcf.gz"
    output:
        "output/" + projectName + "/VCF/combined.vcf.gz.tbi"
    shell:
        "tabix {input}"

rule variant_filtration:
    input:
        "output/" + projectName + "/VCF/combined.vcf"
    output:
        "output/" + projectName + "/VCF_filter_labels/combined.vcf"
    shell:
        "gatk VariantFiltration -V {input} -O {output} " +
            '--filter-expression "QD < 2.0" --filter-name "LowQD" ' +
            '--filter-expression "ReadPosRankSum < -8.0" --filter-name "LowRankSum" ' +
            '--filter-expression "FS > 60.0" --filter-name "HighFS" ' +
            '--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" ' +
            '--filter-expression "MQ < 40.0" --filter-name "LowMQ" '

rule filter_VCF:
    input:
        "output/" + projectName + "/VCF_filter_labels/combined.vcf"
    output:
        "output/" + projectName + "/VCF_PASS/combined.vcf"
    shell:
        "gatk SelectVariants -V {input} --exclude-filtered -O {output}"

rule gene_interval_bam:
    input:
        "output/" + projectName + "/BAM-readgroups/{sample}.bam"
    output:
        "output/" + projectName + "/region_BAM/{region}/{sample}.bam"
    run:
        region = wildcards.region
        sample = wildcards.sample
        chrom = project["regions"][region]["chromosome"]
        start = str(project["regions"][region]["start"]) 
        end = str(project["regions"][region]["end"]) 
        interval = chrom + ":" + start + "-" + end
        cmd = "samtools view -bh {input} " + interval + " > {output}"
        shell(cmd)

rule gene_interval_bam_index:
    input:
        "output/" + projectName + "/region_BAM/{region}/{sample}.bam"
    output:
        "output/" + projectName + "/region_BAM/{region}/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule all_gene_bam:
    input:
        ["output/" + projectName + "/region_BAM/" + region + "/" + sample + ".bam.bai" for region in project["regions"].keys() for sample in get_samples()]

rule gene_interval_vcf:
    input:
        ref="resources/" + projectName + "/genome.fasta",
        vcf="output/" + projectName + "/VCF_PASS/combined.vcf"
    output:
        "output/" + projectName + "/region_VCF/{region}/{sample}.vcf"
    run:
        region = wildcards.region
        sample = wildcards.sample
        chrom = project["regions"][region]["chromosome"]
        start = str(project["regions"][region]["start"]) 
        end = str(project["regions"][region]["end"]) 
        interval = chrom + ":" + start + "-" + end
        cmd = "gatk SelectVariants -R {input.ref} -V {input.vcf} -sn " + sample + " -O {output} -L " + interval
        shell(cmd)

#extractHAIRS doesn't like it where you have more than two alleles - even if only two of them occur in each sample!
rule trim_gene_vcf:
    input:
        "output/" + projectName + "/region_VCF/{region}/{sample}.vcf"
    output:
        "output/" + projectName + "/region_VCF_trimmed/{region}/{sample}.vcf"
    shell:
        "bcftools view --trim-alt-alleles {input} > {output}"

rule exclude_non_variants:
    conda:
        "envs/gatk4.yml"
    input:
        ref="resources/" + projectName + "/genome.fasta",
        vcf="output/" + projectName + "/region_VCF_trimmed/{region}/{sample}.vcf"
    output:
        "output/" + projectName + "/region_VCF_exclude_non_variant/{region}/{sample}.vcf"
    shell:
        "gatk SelectVariants -R {input.ref} -V {input.vcf} --exclude-non-variants -O {output}"

rule fragment_file:
    input:
        bam="output/" + projectName + "/BAM-readgroups/{sample}.bam",
        vcf="output/" + projectName + "/region_VCF_exclude_non_variant/{region}/{sample}.vcf"
    output:
        "output/" + projectName + "/fragment_files/{region}/{sample}_fragment_file"
    shell:
        config["hapcut-path"] + "/extractHAIRS --bam {input.bam} --VCF {input.vcf} --out {output}"

rule hapcut:
    input:
        fragment="output/" + projectName + "/fragment_files/{region}/{sample}_fragment_file",
        vcf="output/" + projectName + "/region_VCF_exclude_non_variant/{region}/{sample}.vcf"
    output:
        "output/" + projectName + "/haplocut_output/{region}/{sample}_haplocut_output"
    shell:
        config["hapcut-path"] + "/HAPCUT2 --fragments {input.fragment} --vcf {input.vcf} --output {output}"

rule phased_vcf:
    conda:
        "envs/fgbio.yml"
    input:
        hapcut="output/" + projectName + "/haplocut_output/{region}/{sample}_haplocut_output",
        vcf="output/" + projectName + "/region_VCF_exclude_non_variant/{region}/{sample}.vcf"
    output:
        "output/" + projectName + "/region_phased_VCF/{region}/{sample}.vcf"
    shell:
        "fgbio HapCutToVcf -i {input.hapcut} -v {input.vcf} -o {output}" 

rule phased_vcf_zip:
    input:
        "output/" + projectName + "/region_phased_VCF/{region}/{sample}.vcf"
    output:
        "output/" + projectName + "/region_phased_VCF/{region}/{sample}.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

rule phased_vcf_index:
    input:
        "output/" + projectName + "/region_phased_VCF/{region}/{sample}.vcf.gz"
    output:
        "output/" + projectName + "/region_phased_VCF/{region}/{sample}.vcf.gz.tbi"
    shell:
        "tabix {input}"

rule haplotypes_full:
    input:
        ref="resources/" + projectName + "/genome.fasta",
        vcf="output/" + projectName + "/region_phased_VCF/{region}/{sample}.vcf.gz",
        index="output/" + projectName + "/region_phased_VCF/{region}/{sample}.vcf.gz.tbi"
    output:
        temp("output/" + projectName + "/haplotype_sequences_full/{region}/haplotype_{hap}/{sample}.fasta")
    shell:
        'bcftools consensus -s {wildcards.sample} -H {wildcards.hap} ' + 
        '-f {input.ref} {input.vcf} > {output}'

rule haplotypes:
    input:
        "output/" + projectName + "/haplotype_sequences_full/{region}/haplotype_{hap}/{sample}.fasta"
    output:
        "output/" + projectName + "/haplotype_sequences/{region}/haplotype_{hap}/{sample}.fasta"
    run:
        region = wildcards.region
        sample = wildcards.sample
        chromosome = project["regions"][region]["chromosome"]
        start = str(project["regions"][region]["start"]) 
        end = str(project["regions"][region]["end"]) 
        interval = chromosome + ":" + start + "-" + end
        cmd = "samtools faidx {input} " + interval + ' | seqkit replace -p ".*" -r "{wildcards.sample}-hap-{wildcards.hap}" > {output}'
        shell(cmd)

rule all_haplotypes:
    input:
        ["output/" + projectName + "/haplotype_sequences/" + region + "/haplotype_" + hap + "/" + sample + ".fasta" for region in project["regions"].keys() for sample in get_samples() for hap in ["1","2"] if not sample in project["exclude_from_VCF"]]

rule combined_haplotype_fasta:
    input:
        ["output/" + projectName + "/haplotype_sequences/{region}/haplotype_" + hap + "/" + sample + ".fasta" for sample in get_samples() for hap in ["1","2"] if not sample in project["exclude_from_VCF"]]
    output:
        "output/" + projectName + "/haplotype_sequences_combined/{region}.fasta" 
    shell:
        "cat {input} > {output}"

rule genetree:
    input:
        "output/" + projectName + "/haplotype_sequences_combined/{region}.fasta"
    output:
        outdir=directory("output/" + projectName + "/gene_trees/{region}/{workflow}")
    shell:
        "ete3 build -w {wildcards.workflow} -n {input} -o {output.outdir}/ --clearall"

rule all_genetrees:
    input:
        ["output/" + projectName + "/gene_trees/" + region + "/" + workflow for region in project["regions"].keys() for workflow in config["ete3-genetree-workflows"].keys()]

def genetree_graphic_input(wildcards):
    region = wildcards["region"]
    workflow = wildcards["workflow"]
    treedir = config["ete3-genetree-workflows"][workflow]["tree-directory"]
    filename = "output/" + projectName + "/gene_trees/" + region + "/" + workflow + "/" + treedir + "/" + region + ".fasta.final_tree.nw"
    return {"tree":filename,"project_JSON_file":project_config_file}

rule genetree_graphic:
    input:
        unpack(genetree_graphic_input)
    output:
        pdf="output/" + projectName + "/gene_tree_graphics/{workflow}/{region}.pdf"
    script:
        "scripts/tree_graphic.R"

rule sort_gff:
    input:
        "resources/Candida_SC5314_genome/GFF/C_albicans_SC5314_features_no_chroms.gff"
    output:
        "resources/Candida_SC5314_genome/GFF/C_albicans_SC5314_features_no_chroms_sorted.gff.gz"
    shell:
        """grep -v "^#" {input} | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > {output}"""

rule index_sorted_gff:
    input:
        "resources/Candida_SC5314_genome/GFF/C_albicans_SC5314_features_no_chroms_sorted.gff.gz"
    output:
        "resources/Candida_SC5314_genome/GFF/C_albicans_SC5314_features_no_chroms_sorted.gff.gz.tbi"
    shell:
        "tabix {input}"

rule vep:
    input:
        gff_index="resources/Candida_SC5314_genome/GFF/C_albicans_SC5314_features_no_chroms_sorted.gff.gz.tbi",
        gff="resources/Candida_SC5314_genome/GFF/C_albicans_SC5314_features_no_chroms_sorted.gff.gz",
        vcf="output/gene_sequences_phased_vcf/{gene}/{sample}.vcf.gz",
        vcf_index="output/gene_sequences_phased_vcf/{gene}/{sample}.vcf.gz.tbi",
        ref=config['SC5314-genome-path'] + '/C_albicans_SC5314_haplotype_A.fasta'
    output:
        vep="output/vep_output/{gene}/{sample}_vep",
        html="output/vep_output/{gene}/{sample}_vep_summary.html"
    shell:
        """vep --species "Candida albicans" --gff {input.gff} --fasta {input.ref} --format vcf -i {input.vcf} -o {output.vep}"""

#rule all_vep:
#    input:
#        ["output/vep_output/" + gene + "/" + sample + "_vep" for gene in config["genes"].keys() for sample in get_candida_albicans_samples()]
