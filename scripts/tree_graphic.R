library(jsonlite)
library(ggtree)
library(treeio)
library(tidyverse)

treefile <- snakemake@input[["tree"]]
project <- jsonlite::fromJSON(snakemake@input[["project_JSON_file"]])

df.strains <- data.frame(
    stringsAsFactors = FALSE,
    strain = names(project$samples),
    Site = sapply(project$samples,FUN=function(s) s[["Site"]]),
    Region = sapply(project$samples,FUN=function(s) s[["Region"]])
)

mytree <- treeio::read.tree(treefile)

mytree$tip.label %>%
    enframe(value='sequence') %>%
    select(-name) %>%
    extract(sequence,into=c('strain','path'),regex='([^-]*)-path-([^-]*)',remove=F) %>%
    select(-path)->
    df.sequence

df <- left_join(df.sequence,df.strains,by='strain')

ggtree(mytree) %<+% df +
    geom_tiplab(aes(fill=Site),color='black',geom='label',label.padding=unit(0.15,"lines"),label.size=0) +
    theme(legend.position='right') -> g.mytree.labelled

split(df,df$strain) %>%
    lapply(FUN=function(d) d$sequence) -> pairs.list

for (strain in names(pairs.list)){
    pair <- pairs.list[[strain]]
    g.mytree.labelled <- g.mytree.labelled +
        geom_taxalink(pair[1],pair[2],color='red',alpha=0.3)
}

ggsave(filename=snakemake@output[['pdf']],g.mytree.labelled,height=20,width=15)