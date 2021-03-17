# Author:                      Job van Riet
# Date:                        14-01-2021
# Function:                    Import unfiltered QIAseq mutations, add the CPCT-02 mPRAD incidence and perform filtering.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)


# Import data -------------------------------------------------------------

# Load the CPCT-02 DR-071 cohort.
load('/mnt/data2/hartwig/DR71/Oct2020/RData/data.Cohort.RData')

# Import unfiltered QIAseq mutations.
unfilteredMuts <- readxl::read_xlsx('Misc./QIAseq_UnfilteredMuts.xlsx', trim_ws = T)

# Remove mutations with a known location in hg19.
unfilteredMuts <- unfilteredMuts %>% dplyr::filter(`[1]LOC_hg19` != 'not found')


# Add Mut. Frequencies ----------------------------------------------------

# Convert mutations into GRanges to check overlap.
muts.Gr <- GenomicRanges::sort(GenomicRanges::GRanges(seqnames = gsub(':.*', '', unfilteredMuts$`[1]LOC_hg19`), ranges = IRanges::IRanges(start = as.numeric(gsub('.*:', '', unfilteredMuts$`[1]LOC_hg19`)), end = as.numeric(gsub('.*:', '', unfilteredMuts$`[1]LOC_hg19`)))))
muts.Gr$`[1]LOC_hg19` <- unfilteredMuts$`[1]LOC_hg19`
muts.Gr <- unique(muts.Gr)

# Perform overlap.
muts.CPCT <- IRanges::subsetByOverlaps(data.Cohort$somaticVariants, muts.Gr, type = 'within')
muts.CPCT$LOC <- sprintf('%s:%s', seqnames(muts.CPCT), start(muts.CPCT))

muts.CPCT <- tibble::as_tibble(muts.CPCT) %>%
    dplyr::group_by(LOC) %>%
    dplyr::summarise(totalCPCT = dplyr::n_distinct(sample)) %>%
    dplyr::ungroup()

# Add back to QIAseq mutations.
unfilteredMuts <- unfilteredMuts %>% dplyr::left_join(muts.CPCT, by = c('[1]LOC_hg19' = 'LOC'))


# Perform filtering. ------------------------------------------------------

# Determine whether mutations overlap the panel-design.
regions.panel <- c(rtracklayer::import.bed('Misc./QIAseq_DNA_panel.CDHS-25058Z-3602.primers-150bphg19.bed'), rtracklayer::import.bed('Misc./QIAseq_DNA_panel.CDHS-25058Z-3602.roihg19.bed'))
regions.panel <- GenomicRanges::reduce(regions.panel)

overlappingPanel <- IRanges::subsetByOverlaps(muts.Gr, regions.panel, type = 'within')

# Add within panel information back to the table.
unfilteredMuts$withinPanel <- unfilteredMuts$`[1]LOC_hg19` %in% overlappingPanel$`[1]LOC_hg19`

unfilteredMuts <- unfilteredMuts %>%
    dplyr::group_by(`[1]LOC_hg19`) %>%
    dplyr::mutate(totalCaba = dplyr::n_distinct(lcode)) %>%
    dplyr::ungroup()

# Perform filtering.
filteredMuts <- unfilteredMuts %>%
    dplyr::filter(
        # Filter mutations within panel-design.
        (withinPanel | totalCPCT >= 5),
        # Filter on incidence within CABA-V7 samples.
        (totalCaba <= 30 | totalCPCT >= 5),
        # Filter out SNPs.
        (SNP == 'n' | totalCPCT >= 5),
        # Filter on DP.
        (`[35]consensus-5:DP` >= 100 | totalCPCT >= 5),
        # Filter on min. VAF.
        ((as.numeric(`[34]consensus-5:AF`) >= 0.05 & !is.na(as.numeric(`[34]consensus-5:AF`))) | totalCPCT >= 5),
        # Filter on max. VAF.
        ((as.numeric(`[34]consensus-5:AF`) <= 0.95 & !is.na(as.numeric(`[34]consensus-5:AF`))) | totalCPCT >= 5),
        # Filter on coding.
        `[19]ann_hgvs_p` != '.'
    ) %>%
    dplyr::group_by(`[1]LOC_hg19`) %>%
    dplyr::mutate(totalCaba = dplyr::n_distinct(lcode)) %>%
    dplyr::ungroup()


# Output mutations --------------------------------------------------------

write.table(filteredMuts, file = 'asd.txt', quote = F, sep = '\t', row.names = F)
