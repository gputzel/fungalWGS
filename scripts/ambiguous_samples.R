library(tidyverse)

read_table(snakemake@input[['distmat']],skip = 1,col_names = F) %>%
    column_to_rownames('X1') %>%
    as.matrix -> dist.m
colnames(dist.m) <- rownames(dist.m)

rownames(dist.m) %>%
    enframe(value='path') %>%
    select(-name) %>%
    separate(col=path,into=c('sample','X','pathnumber'),sep='-') %>%
    select(-X) %>%
    group_by(sample) %>%
    summarize(num.paths = n()) %>%
    arrange(-num.paths) -> num.paths.df

num.paths.df %>%
    filter(num.paths > 2) -> ambiguous.samples.df

num.paths.df %>%
    filter(num.paths < 3) -> unambiguous.samples.df

ambiguous.samples.df %>%
    write_tsv(snakemake@output[['txt']])

ambiguous.samples <- ambiguous.samples.df$sample

best.paths <- function(sample,distance.matrix,ambiguous.samples){
    sample.paths <- grep(pattern=paste0(sample,"-"),rownames(distance.matrix),value=T)
    ambiguous.paths.ind <- sapply(strsplit(rownames(dist.m),split="-"),'[[',1) %in% ambiguous.samples
    submatrix <- distance.matrix[sample.paths,!ambiguous.paths.ind]
    min.distances <- apply(submatrix,MARGIN=1,min)
    names(sort(min.distances)[1:2])
}

filtered.ambiguous.paths <- do.call(c,
    lapply(ambiguous.samples,
       function(s) best.paths(s,distance.matrix = dist.m,
                              ambiguous.samples = ambiguous.samples)
    )
)

rownames(dist.m) %>%
    enframe(value='path') %>%
    select(-name) %>%
    separate(col=path,into=c('sample','X','pathnumber'),sep='-',remove=F) %>%
    select(-X) %>%
    filter(sample %in% unambiguous.samples.df$sample) %>%
    .$path -> unambiguous.paths

c(unambiguous.paths,filtered.ambiguous.paths) %>%
    enframe(value='path') %>%
    select(-name) %>%
    write_tsv(snakemake@output[['filtered']],col_names = F)