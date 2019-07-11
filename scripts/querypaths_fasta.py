#print(snakemake.input)
#print(snakemake.wildcards.region)

sample=snakemake.wildcards.sample

with open(snakemake.input.tsv,'r') as fi:
    _ = fi.readline()
    with open(snakemake.output.fasta,'w') as fo:
        for i,l in enumerate(fi.readlines()):
            fo.write(">" + sample + "-path-" + str(i+1) + '\n')
            fo.write(l.rstrip('\n').split('\t')[11] + '\n')
