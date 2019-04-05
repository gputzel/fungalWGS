from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from Bio.Seq import Seq
HTTP = HTTPRemoteProvider()

configfile: "config.json"

def get_samples():
    samples,= glob_wildcards(config['FASTQ-path'] + '/{ID}_1.fq.gz')
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
