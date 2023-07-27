#!/usr/bin/Rscript

suppressWarnings(library(ggplot2, quietly = T))
suppressWarnings(library(viridis, quietly = T))
suppressWarnings(library(seqinr, quietly = T))
suppressMessages(library(dplyr, quietly = T))

pdf(NULL)

parser <- argparse::ArgumentParser(description = 'Draw common pairwise k-mers.')
parser$add_argument('fasta', metavar = 'FILE',
    help = 'Input FASTA file.')
parser$add_argument('-n', '--names', metavar = 'STR', required = T,
    help = 'Two sequences names, through a comma.')
parser$add_argument('-k', '--kmer', metavar = 'STR', default = 25, type = 'integer',
    help = 'k-mer size [%(default)s].')
parser$add_argument('-f', '--max-freq', metavar = 'INT', default = 10, type = 'integer',
    help = 'Ignore k-mers with frequency over INT [%(default)s].')
parser$add_argument('--no-canon', action = 'store_false', dest = 'canon',
    help = 'Do not canonize k-mers (forward and reverse-compl will not match).')
parser$add_argument('-o', '--out', metavar = 'FILE', required = T,
    help = 'Output plot file (PNG, PDF, etc.). Literal {} is replaced with --names STR.')
args <- parser$parse_args()

###################

fasta <- read.fasta(args$fasta)
names <- unlist(strsplit(args$names, ',', fixed = T))
if (length(names) != 2) {
    stop('There must be exactly two names.')
}
name1 <- names[1]
name2 <- names[2]
rm(names)
seq1 <- toupper(paste(fasta[[name1]], collapse = ''))
seq2 <- toupper(paste(fasta[[name2]], collapse = ''))

###################

rcomp <- function(seq) {
    stringi::stri_replace_all_fixed(seq,
        pattern = c('A', 'C', 'G', 'T'), replacement = c('t', 'g', 'c', 'a'),
        vectorize_all = FALSE) |>
    sapply(function(x) rev(toupper(x)))
}

get_kmers <- function(seq, k, canon = T) {
    n <- nchar(seq)
    stopifnot(k <= n)
    fkmers <- sapply(1:(n + 1 - k), function(i) substr(seq, i, i + k - 1))
    if (canon) {
        rkmers <- rcomp(fkmers)
        pmin(fkmers, rkmers)
    } else {
        fkmers
    }
}

####################

kmers1 <- get_kmers(seq1, args$k, args$canon)
kmers2 <- get_kmers(seq2, args$k, args$canon)
freq <- table(c(kmers1, kmers2))

kmers1 <- data.frame(kmer = kmers1) |>
    mutate(pos = 1:n(), freq = as.vector(freq[kmer])) |>
    filter(freq <= args$max_freq)
kmers2 <- data.frame(kmer = kmers2) |>
    mutate(pos = 1:n(), freq = as.vector(freq[kmer])) |>
    filter(freq <= args$max_freq)
kmer_pairs <- inner_join(
    kmers1, kmers2, by = c('kmer', 'freq'), suffix = c('1', '2'), relationship = 'many-to-many')

n1 <- length(unique(kmers1$kmer))
n2 <- length(unique(kmers2$kmer))
n12 <- length(unique(kmer_pairs$kmer))
jaccard <- n12 / (n1 + n2 - n12)

max_freq_legend <- min(20, args$max_freq)
freq_limits <- c(2, max_freq_legend)

ggplot(kmer_pairs) +
    geom_point(aes(pos1, pos2, color = pmax(2, pmin(freq, max_freq_legend))),
        size = 0.1) +
    labs(title = sprintf('%s - %s', name1, name2),
        subtitle = sprintf('Jaccard index: %.5f', jaccard)) +
    scale_x_continuous(sprintf('Position (%s)', name1), expand = expansion(mult = 0.02)) +
    scale_y_continuous(sprintf('Position (%s)', name2), expand = expansion(mult = 0.02)) +
    coord_fixed() +
    scale_color_viridis('k-mer frequency',
        limits = freq_limits) +
    theme_bw() +
    theme(
        legend.position = 'bottom',
        legend.key.height = unit(0.8, 'lines'),
    )

out_filename <- gsub('{}', args$names, args$out, fixed = T)
ggsave(out_filename, width = 8, height = 8.7, dpi = 400)
