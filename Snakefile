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
