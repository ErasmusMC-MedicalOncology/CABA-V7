# Author:                      Job van Riet
# Date:                        29-04-2021
# Function:                    Import unfiltered QIAseq mutations, add the CPCT-02 mPRAD incidence and perform filtering.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)


# Import data -------------------------------------------------------------

# Import unfiltered QIAseq mutations.
unfilteredMuts <- dplyr::bind_rows(pbapply::pblapply(list.files('/mnt/data2/hartwig/DR71/Misc/Annotated/', pattern = 'tsv.gz$', full.names = T), function(x){

    # Import.
    data <- readr::read_tsv(x, guess_max = 299999, trim_ws = T)
    data$lcode <- gsub('.ann.tsv.gz', '', basename(x))

    # Clean-up columns.
    colnames(data) <- gsub('# ', '', gsub('-|:', '', gsub(unique(data$lcode), '', colnames(data))))

    # Filter non-standard chromosomes.
    data <- data %>% dplyr::filter(`[1]CHROM` %in% c(1:22, 'X', 'Y'))

    # Fix nomenclature.
    data <- data %>% dplyr::mutate(`[1]CHROM` = paste0('chr', `[1]CHROM`))

    # Add position-based identifier.
    data <- data %>% dplyr::mutate(ID_Job = sprintf('%s.%s:%s>%s', `[1]CHROM`, `[2]POS`, `[4]REF`, `[5]ALT`))

    # Retrieve ALT AF.
    data$`[28]consensus5AF_Job` <- unlist(lapply(strsplit(data$`[28]consensus5AF`, ','), function(x) as.numeric(x)[2]))
    data$`[28]consensus5AF_Job` <- ifelse(is.infinite(data$`[28]consensus5AF_Job`), NA, data$`[28]consensus5AF_Job`)

    # Fix gene names.
    data <- data %>% dplyr::mutate(`[17]ann_gene_name_Job` = gsub('\\|.*', '', `[17]ann_gene_name`))

    # Fix annotations.
    data <- data %>% dplyr::mutate(`[15]ann_annotation_Job` = gsub('\\|.*', '', `[15]ann_annotation`))

    # Remove duplicates (if any)
    data <- data %>% dplyr::distinct()

    # Return cleaned data.
    return(data)

}, cl = 5))


# Perform filtering. ------------------------------------------------------

## Determine whether mutations overlap the panel-design. ----

regions.panel <- c(rtracklayer::import.bed('Misc./QIAseq_DNA_panel.CDHS-25058Z-3602.primers-150bphg38.bed'), rtracklayer::import.bed('Misc./QIAseq_DNA_panel.CDHS-25058Z-3602.roihg38.bed'))
BiocGenerics::strand(regions.panel) <- '*'

regions.panel <- GenomicRanges::reduce(regions.panel)

muts.Gr <- GenomicRanges::makeGRangesFromDataFrame(unfilteredMuts, keep.extra.columns = F, seqnames.field = '[1]CHROM', start.field = '[2]POS', end.field = '[2]POS')
muts.Gr$ID_Job <- unfilteredMuts$ID_Job
overlappingPanel <- IRanges::subsetByOverlaps(muts.Gr, regions.panel)

# Add within panel information back to the table.
unfilteredMuts$withinPanel <- unfilteredMuts$ID_Job %in% overlappingPanel$ID_Job


## Determine presence within PON ----

HMF.PON <- VariantAnnotation::readVcfAsVRanges('/mnt/data2/hartwig/DR71/Misc/SageGermlinePon.98x.38.vcf.gz')
HMF.PON.Filtered <- IRanges::subsetByOverlaps(HMF.PON, muts.Gr)
HMF.PON.Filtered <- HMF.PON.Filtered[HMF.PON.Filtered$PON_COUNT >= 10 & HMF.PON.Filtered$PON_MAX >= 5]
HMF.PON.Filtered$ID_Job <- sprintf('%s.%s:%s>%s', GenomeInfoDb::seqnames(HMF.PON.Filtered), BiocGenerics::start(HMF.PON.Filtered), VariantAnnotation::ref(HMF.PON.Filtered), VariantAnnotation::alt(HMF.PON.Filtered))

# Add to table.
unfilteredMuts$inPON <- unfilteredMuts$ID_Job %in% HMF.PON.Filtered$ID_Job


## Perform filtering. ----

filteredMuts <- unfilteredMuts %>%
    dplyr::filter(
        # Filter mutations within panel-design.
        withinPanel,
        # Filter on HMF PON.
        !inPON,
        # Filter on DP.
        `[29]consensus5DP` >= 100,
        # Filter on min. VAF.
        (`[28]consensus5AF_Job` >= 0.05 & !is.na(`[28]consensus5AF_Job`)),
        # Filter on max. VAF.
        (`[28]consensus5AF_Job` <= 0.95 & !is.na(`[28]consensus5AF_Job`)),
        # Filter on coding.
        `[18]ann_hgvs_p` != '.',
        # Only take known COSMIC mutations.
        grepl('COS', `[3]ID`)
    ) %>%
    # Filter on incidence within CABA-V7 samples.
    dplyr::group_by(ID_Job) %>%
    dplyr::mutate(totalCaba = dplyr::n_distinct(lcode)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(totalCaba <= 30)


# Output mutations --------------------------------------------------------

filteredMuts <- filteredMuts %>%
    dplyr::group_by(ID_Job) %>%
    dplyr::mutate(totalCaba = dplyr::n_distinct(lcode)) %>%
    dplyr::ungroup()

write.table(filteredMuts, file = 'filteredMuts.txt', quote = F, sep = '\t', row.names = F)
