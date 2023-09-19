#!/usr/bin/Rscript

pdf(NULL)
suppressWarnings(library(ggplot2, quietly = T))
suppressMessages(library(ggnewscale, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(stringi, quietly = T))

msg <- function(...) cat(sprintf(...), sep='', file=stderr())
int_fmt <- scales::label_comma(accuracy = 1)
float_fmt2 <- scales::label_comma(accuracy = 0.01)

parser <- argparse::ArgumentParser(description = 'Draw read assignment.')
parser$add_argument('dir', metavar = 'DIR',
    help = 'Input directory with locus analysis.')
parser$add_argument('-g', '--genotypes', metavar = 'STR', nargs='+',
    help = 'Genotype(s) to draw.')
parser$add_argument('--no-lik', action = 'store_true',
    help = 'Do not draw likelihood axis.')
parser$add_argument('-o', '--out', metavar = 'FILE', default = '{dir}/plots/{gt}.png',
    help = paste('Output plot file (PNG, PDF, etc.).',
        'Literals {dir} and {gt} are replaced with input directory and genotype.',
        'Default: %(default)s.'))
args <- parser$parse_args()

if (length(args$genotype) > 1 && !grepl('{gt}', args$out, fixed = T)) {
    stop('More than one genotype provided and output file does not contain `{gt}`')
}

# Load data.

msg('Loading data\n')
dir <- args$dir
full_sol <- read.csv(file.path(dir, 'depth.csv.gz'), sep = '\t', comment = '#') |>
    filter(genotype %in% args$genotypes)
if (nrow(full_sol) == 0) {
    stop('Cannot find any of the genotypes!')
}

# Select scale limits.

msg('Processing dataset\n')
no_lik <- args$no_lik
wo_summary <- filter(full_sol, contig != 'summary')
max_depth <- max(wo_summary$depth)
lik_range <- if (no_lik) { c(0, 0) } else { range(wo_summary$lik) }
min_lik <- lik_range[1] * 1.00001

depth_breaks <- scales::breaks_extended(3)(c(0, max_depth))
lik_breaks <- scales::breaks_extended(3)(c(min_lik, 0))

lik_axis_mult <- -0.5 * (max_depth / min_lik)
breaks <- c(lik_axis_mult * lik_breaks, depth_breaks)
labels_left <- c(rep('', length(lik_breaks)), depth_breaks)
labels_right <- c(lik_breaks, rep('', length(depth_breaks)))
ylim <- (if (no_lik) { c(0, max_depth) } else { c(lik_axis_mult * min_lik, max_depth) }) * 1.03

# Colors

fill_colors <- rev(c('#355070', '#6D597A', '#B56576', '#E56B6F'))
main_color <- tail(fill_colors, n = 1)
max_lik = max(c(min_lik, lik_range[2], -1))
fill_rescale <- 1 - c(seq(min_lik, max_lik, length.out = length(fill_colors)), 0) / min_lik

for (gt_str in args$genotype) {
    msg('    Processing genotype %s\n', gt_str)
    sol <- filter(full_sol, genotype == gt_str)
    if (nrow(sol) == 0) {
        msg('ERROR: Cannot find genotype %s\n', gt_str)
        next
    }

    out_filename <- stri_replace_all_fixed(args$out,
        pattern = c('{dir}', '{gt}'), replacement = c(dir, gt_str), vectorize_all = F)
    dir.create(dirname(out_filename), showWarnings = FALSE)

    gt <- strsplit(gt_str, ',', fixed = T)[[1]]
    summary_info <- filter(sol, contig == 'summary') |>
        rename(info = window) |>
        mutate(depth = NULL, lik = NULL, lik = as.numeric(sub(' .*', '', sub('.*lik=', '', info)))) |>
        arrange(-lik) |> head(n = 1)
    sol <- filter(sol, stage == summary_info$stage & attempt == summary_info$attempt)
    if (sum(sol$contig == 'summary') != 1) {
        stop('Unexpected number of solutions!')
    }
    sol <- filter(sol, contig != 'summary') |>
        mutate(contig_ix = as.numeric(contig),
            contig = gt[contig_ix], contig_ext = sprintf('%s-%d', contig, contig_ix), window = as.numeric(window))

    # Parse information.

    info <- unlist(strsplit(summary_info$info, ' ')) |>
        strsplit('=')
    info <- setNames(as.numeric(sapply(info, `[`, 2)), sapply(info, `[`, 1))

    # Create title and subtitle.

    title <- paste(gt, collapse = ', ')
    subtitle <- paste(
        sprintf('Total reads: %s  (unmapped: %s,  out of bounds: %s).',
            int_fmt(info['reads']), int_fmt(info['unmapped']), int_fmt(info['boundary'])),
        sprintf('log₁₀-likelihood: %s  =  [alignment] %s × %s  +  [depth] %s × %s',
            float_fmt2(info['lik']), float_fmt2(info['aln_contrib']), float_fmt2(info['aln_lik']),
            float_fmt2(info['depth_contrib']), float_fmt2(info['depth_lik'])),
        sep = '\n')

    # Drawing and saving.
    msg('    Drawing read depth for %s to %s\n', gt_str, out_filename)
    nwindows <- max(sol$window)
    ggplot(sol) +
        # Window weights behind likelihoods.
        geom_rect(aes(xmin = window - 0.5, xmax = window + 0.5,
            ymin = ifelse(no_lik, 0, -Inf), ymax = ifelse(no_lik, Inf, 0), fill = weight)) +
        scale_fill_gradientn('Window weight ', limits = c(0, 1), colors = c('#ff000055', '#ffffff00')) +
        new_scale('fill') +

        # 0-line.
        geom_hline(yintercept = 0, color = main_color) +
        # Read depth bars.
        geom_bar(aes(window, pmin(depth, max_depth), fill = pmax(lik, min_lik)), stat = 'identity', width = 1) +
        # Likelihood points.
        (if (no_lik) {
            list()
        } else {
            geom_point(aes(window, lik_axis_mult * pmax(lik, min_lik)), size = 0.5, color = main_color)
        }) +

        facet_wrap(~ contig_ext, ncol = 1, strip.position = 'top',
            labeller = as_labeller(function(x) sub('-[0-9]$', '', x))) +
        ggtitle(title, subtitle) +
        coord_cartesian(ylim = ylim) +
        scale_x_continuous('Window', expand = expansion(mult = 0.005),
            breaks = seq(0, nwindows, 100), minor_breaks = seq(0, nwindows, 10)) +
        scale_y_continuous('Read depth',
            expand = c(0, 0),
            breaks = breaks,
            minor_breaks = NULL,
            labels = labels_left,

            sec.axis = (if (no_lik) {
                waiver()
            } else {
                dup_axis(
                name = 'Depth log₁₀-likelihood',
                labels = labels_right)
            })) +
        scale_fill_gradientn('log₁₀-likelihood',
            colors = c(fill_colors, main_color),
            values = fill_rescale,
            limits = c(min_lik, 0),
            breaks = lik_breaks,
            minor_breaks = NULL,
            expand = c(0, 0)) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = 'gray95', color = NA),
            strip.text = element_text(face = 'bold', margin = margin(t = 2, b = 2)),

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
    ggsave(out_filename, width = 9, height = 6, dpi = 400, scale = 1.2)
}
warnings()
msg('Successfully finished\n')
