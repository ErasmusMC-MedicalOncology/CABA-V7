# Author:                      Job van Riet
# Date:                        21-04-2021
# Function:                    Determine mFASTSeqS deviations between T1 and T2.


# Import libraries --------------------------------------------------------

library(dplyr)
library(ggplot2)
library(multcomp)

# Import data -------------------------------------------------------------

## Import patient metadata ----

overviewPatients <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 1)

## Import mFAST-Seq - T1 ----
dataFastSeq.T1 <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 6) %>%
    reshape2::melt(id = 'Chromosomal arm', variable.name = 'L-code', value.name = 'value.T1') %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        value.T1 = as.numeric(value.T1),
        `Chromosomal arm` = gsub('chr', '', `Chromosomal arm`)
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `L-code`))

## Import mFAST-Seq - T2 ----
dataFastSeq.T2 <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 7) %>%
    reshape2::melt(id = 'Chromosomal arm', variable.name = 'Follow up L-code', value.name = 'value.T2') %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        value.T2 = as.numeric(value.T2),
        `Chromosomal arm` = gsub('chr', '', `Chromosomal arm`)
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `Follow up L-code`))

## Merge mFAST-Seq T1 and T2 ----
dataFastSeq.Combined <- dataFastSeq.T1 %>%
    dplyr::inner_join(dataFastSeq.T2, by = c('Chromosomal arm', 'Subject Number')) %>%
    dplyr::mutate(
        chr = gsub('p|q', '', `Chromosomal arm`),
        `Chromosomal arm` = factor(`Chromosomal arm`, levels = gtools::mixedsort(unique(`Chromosomal arm`))),
        value.deltaT2vsT1 = value.T2 - value.T1
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `AR-V7 (Baseline)`, `AR-V7 Conversion`, `Genome-Wide Z Score (Baseline)`, `Genome-Wide Z Score (T2)`), by = 'Subject Number') %>%

    # Remove samples with no first and second AR-V7 determination.
    dplyr::filter(!is.na(`AR-V7 Conversion`), !is.na(`AR-V7 (Baseline)`)) %>%

    # Add nr. of patients per conversion group.
    dplyr::group_by(`AR-V7 Conversion`) %>%
    dplyr::mutate(`AR-V7 Conversion` = sprintf('%s (<i>n</i> = %s)', ifelse(`AR-V7 Conversion` == 'Neg.', 'AR-V7<sup>Conv.</sup>', ifelse(`AR-V7 Conversion` == 'Pos.', 'AR-V7<sup>Pos.</sup>', 'Und.')), dplyr::n_distinct(`Subject Number`))) %>%
    dplyr::ungroup() %>%

    # Remove acrocentric chromosomes which we cannot measure
    dplyr::filter(!is.na(value.deltaT2vsT1))


# Test chromosomal disparities --------------------------------------------

# Make model based on chromosomal arm vs. AR-V7 status (T2) and use deltaT as measurement.
fit <- aov(value.deltaT2vsT1 ~ `Chromosomal arm` * `AR-V7 Conversion`, data = dataFastSeq.Combined)
multiFit <- multcomp::glht(fit, alternative = 'two.sided')
multiFit <- summary(multiFit, test = adjusted('BH'))
multiFit.df <- tibble::tibble(test = names(multiFit$test$sigma), chr = gsub('.*`', '', gsub(':.*', '', names(multiFit$test$sigma))), group = ifelse(grepl('Pos.', names(multiFit$test$sigma)), 'AR-V7 Retainers', 'AR-V7 Loss'), p.adj = multiFit$test$pvalues, tstat = multiFit$test$tstat, coef = multiFit$test$coefficients, sigma = multiFit$test$sigma)

dataFastSeq.Combined %>%
    dplyr::group_by(`Chromosomal arm`, `AR-V7 Conversion`) %>%
    dplyr::mutate(medianDelta = round(mean(value.deltaT2vsT1), 1)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(., ggplot2::aes(x = `Chromosomal arm`, y = value.deltaT2vsT1, fill = `AR-V7 Conversion`, color = `AR-V7 Conversion`, label = medianDelta)) +
    ggplot2::geom_hline(yintercept = 0, color = 'red', lty = 'dashed') +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = .25, color = 'black') +
    ggbeeswarm::geom_quasirandom(size = .5) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(-150, -100, -50, -25, -10, -5, -2.5, 0, 2.5, 5, 10, 25, 50, 150), limits = c(-25, 25)) +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = F) +
    ggplot2::scale_color_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = F) +
    ggplot2::stat_summary(ggplot2::aes(y = stage(value.deltaT2vsT1, after_stat = -25)), fun = min, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, vjust = .5, angle = 0) +
    ggplot2::labs(x = NULL, y = 'Delta Z-score (T2 - T1)') +
    ggplot2::facet_wrap(. ~`AR-V7 Conversion`, nrow = 2) +
    theme_Job + ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank())
