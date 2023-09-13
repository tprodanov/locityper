#!/usr/bin/Rscript

pdf(NULL)
suppressWarnings(library(ggplot2, quietly = T))
suppressMessages(library(ggnewscale, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(stringi, quietly = T))

parser <- argparse::ArgumentParser(description = 'Draw read assignment.')
parser$add_argument('dir', metavar = 'DIR',
    help = 'Input directory with locus analysis.')
parser$add_argument('-g', '--genotype', metavar = 'STR',
    help = paste('Genotype to draw (through a comma).',
    'Required, unless there is only one genotype in the solution.'))
parser$add_argument('-s', '--stage', metavar = 'INT', type = 'integer',
    help = 'Draw this stage. Default: stage with highest likelihood.')
parser$add_argument('-a', '--attempt', metavar = 'INT', type = 'integer',
    help = 'Draw this stage iteration. Default: iteration with highest likelihood.')
parser$add_argument('-d', '--depth', metavar = 'FLOAT', type = 'double',
    help = 'Read depth limit (inferred by default).')
parser$add_argument('-l', '--lik', metavar = 'FLOAT', type = 'double',
    help = 'Likelihood limit (inferred by default).')
parser$add_argument('--no-lik', action = 'store_true',
    help = 'Do not draw likelihood axis.')
parser$add_argument('-o', '--out', metavar = 'FILE',
    help = paste('Output plot file (PNG, PDF, etc.).',
        'Literal {%s} is replaced with the value for keys "gt", "stage", "attempt".',
        'Default: "<dir>/plots/{gt}.png".'))
args <- parser$parse_args()

# Load data.

dir <- args$dir
sol <- read.csv(file.path(dir, 'depth.csv.gz'), sep = '\t', comment = '#')

# Filter read depth according to genotype and the best likelihood.

if (is.null(args$genotype)) {
    if (length(unique(sol$genotype)) > 1) {
        stop(paste('There is more than 1 genotype in the solution,',
            'must provide -g/--genotype <str>'))
    }
    args$genotype <- sol$genotype[1]
} else {
    sol <- filter(sol, genotype == args$genotype)
}
if (!is.null(args$stage)) {
    sol <- filter(sol, stage == args$stage)
}
if (!is.null(args$attempt)) {
    sol <- filter(sol, attempt == args$attempt)
}
if (nrow(sol) == 0) {
    stop('Empty dataframe after filtering!')
}
gt <- strsplit(args$genotype, ',', fixed = T) |> unlist()

summary_info <- filter(sol, contig == 'summary') |>
    rename(info = window) |>
    mutate(depth = NULL, lik = NULL,
        lik = as.numeric(sub(' .*', '', sub('.*lik=', '', info)))) |>
    arrange(-lik) |> head(n = 1)
sol <- filter(sol,
    stage == summary_info$stage & attempt == summary_info$attempt)
if (sum(sol$contig == 'summary') != 1) {
    stop('Unexpected number of solutions!')
}
sol <- filter(sol, contig != 'summary') |>
    mutate(
        contig_ix = as.numeric(contig),
        contig = gt[contig_ix],
        contig_ext = sprintf('%s-%d', contig, contig_ix),
        window = as.numeric(window))

# Parse information.

info <- unlist(strsplit(summary_info$info, ' ')) |>
    strsplit('=')
info <- setNames(as.numeric(sapply(info, `[`, 2)), sapply(info, `[`, 1))

# Create title and subtitle.

comma_int <- scales::label_comma(accuracy = 1)
comma_prec2 <- scales::label_comma(accuracy = 0.01, style_positive = 'plus')
title <- paste(gt, collapse = ', ')
subtitle <- paste(
    sprintf('Reads: %s,  Unmapped: %s,  On boundary: %s.',
        comma_int(info['reads']),
        comma_int(info['unmapped']),
        comma_int(info['boundary'])),
    sprintf('log10-likelihood: %s  =  %s (alns)   %s (depth)',
        comma_prec2(info['lik']),
        comma_prec2(info['aln_lik']),
        comma_prec2(info['depth_lik'])),
    sep = '\n')

# Select scale limits.

no_lik <- args$no_lik
max_depth <- (if (is.null(args$depth)) { max(sol$depth) } else { args$depth })
min_lik <- (if (is.null(args$lik)) { min(sol$lik) * 1.00001 } else { -abs(args$lik) })

depth_breaks <- scales::breaks_extended(3)(c(0, max_depth))
lik_breaks <- scales::breaks_extended(3)(c(min_lik, 0))

lik_axis_mult <- -0.5 * (max_depth / min_lik)
breaks <- c(lik_axis_mult * lik_breaks, depth_breaks)
labels_left <- c(rep('', length(lik_breaks)), depth_breaks)
labels_right <- c(lik_breaks, rep('', length(depth_breaks)))
ylim <- (if (no_lik) { c(0, max_depth) } else { c(lik_axis_mult * min_lik, max_depth) }) * 1.01

# Drawing and saving.

fill_colors <- rev(c('#355070', '#6D597A', '#B56576', '#E56B6F'))
main_color <- tail(fill_colors, n = 1)
fill_rescale <- -c(seq(min_lik, -1, length.out = length(fill_colors)), 0) /
    min_lik + 1
nwindows <- max(sol$window)

ggplot(sol) +
    # Window weights behind likelihoods.
    geom_rect(aes(xmin = window - 0.5, xmax = window + 0.5,
        ymin = ifelse(no_lik, 0, -Inf),
        ymax = ifelse(no_lik, Inf, 0),
        fill = weight)) +
    scale_fill_gradientn('Window weight', limits = c(0, 1),
        colors = c('#ff000055', '#ffffff00')) +
    new_scale('fill') +

    # 0-line.
    geom_hline(yintercept = 0, color = main_color) +
    # Read depth bars.
    geom_bar(aes(window, pmin(depth, max_depth), fill = pmax(lik, min_lik)),
        stat = 'identity', width = 1) +
    # Likelihood points.
    (if (no_lik) {
        list()
    } else {
        geom_point(aes(window, lik_axis_mult * pmax(lik, min_lik)),
          size = 0.5, color = main_color)
    }) +

    facet_wrap(~ contig_ext, ncol = 1, strip.position = 'top',
        labeller = as_labeller(function(x) sub('-[0-9]$', '', x))) +
    ggtitle(title, subtitle) +
    scale_x_continuous('Window',
        expand = expansion(mult = 0.005),
        breaks = seq(0, nwindows, 100),
        minor_breaks = seq(0, nwindows, 10),
        ) +
    scale_y_continuous('Read depth',
        expand = c(0, 0),
        limits = ylim,
        breaks = breaks,
        minor_breaks = NULL,
        labels = labels_left,

        sec.axis = (if (no_lik) {
            waiver()
        } else {
            dup_axis(
            name = 'Depth log10-likelihood',
            labels = labels_right)
        })) +
    scale_fill_gradientn('Depth log10-likelihood',
        colors = c(fill_colors, main_color),
        values = fill_rescale,
        limits = c(min_lik, 0),
        breaks = lik_breaks,
        minor_breaks = NULL,
        expand = c(0, 0)) +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = 'gray95', color = NA),
        strip.text = element_text(
            face = 'bold',
            margin = margin(t = 2, b = 2)),

        legend.position = 'bottom',
        legend.box.margin = margin(t = -10, b = -5),
        legend.key.height = unit(0.8, 'lines'),
        panel.grid.minor = element_line(linetype = 'dashed'),
        panel.background = element_rect(fill = NA),

        axis.title.x = element_text(margin = margin(t = 0)),
        axis.title.y.right = element_text(margin = margin(l = 7)),
        axis.text.y.right = element_text(hjust = 1),
        axis.ticks.y = element_blank(),
    )

if (is.null(args$out)) {
    args$out <- file.path(dir, 'plots/{gt}.png')
}

out_filename <- stri_replace_all_fixed(
    args$out,
    pattern = c('{gt}', '{stage}', '{attempt}'),
    replacement = c(summary_info$genotype, summary_info$stage, summary_info$attempt),
    vectorize_all = F
)
dir.create(dirname(out_filename), showWarnings = F)
ggsave(out_filename, width = 9, height = 6, dpi = 400, scale = 1.2)
warnings()
