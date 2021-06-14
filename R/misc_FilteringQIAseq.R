# Author:                      Job van Riet
# Date:                        14-06-2021
# Function:                    Import unfiltered QIAseq mutations and perform (germline) filtering.


# Import libraries --------------------------------------------------------

library(dplyr)
library(plyr)
library(ggplot2)


# Helper functions --------------------------------------------------------

readVCF.Biomics <- function(x){

    # Import.
    data <- readr::read_tsv(x, guess_max = 299999, trim_ws = T)
    data$biomicsID <- gsub('.ann.tsv.gz', '', basename(x))

    # Clean-up columns.
    colnames(data) <- gsub('# ', '', gsub('-|:', '', gsub(unique(data$biomicsID), '', colnames(data))))

    # Filter non-standard chromosomes.
    data <- data %>% dplyr::filter(`[1]CHROM` %in% c(1:22, 'X', 'Y'))

    # Fix nomenclature.
    data <- data %>% dplyr::mutate(`[1]CHROM` = paste0('chr', `[1]CHROM`))

    # Add position-based identifier.
    data <- data %>% dplyr::mutate(ID_Job = sprintf('%s.%s:%s>%s', `[1]CHROM`, `[2]POS`, `[4]REF`, `[5]ALT`))

    # Retrieve ALT AF.
    data$`[28]consensus5AF_Job` <- unlist(lapply(strsplit(data$`[28]consensus5AF`, ','), function(x) as.numeric(x)[2]))
    data$`[28]consensus5AF_Job` <- ifelse(is.infinite(data$`[28]consensus5AF_Job`), NA, data$`[28]consensus5AF_Job`)

    # Retrieve ALT AD.
    data$consensus5AltDepth <- unlist(lapply(strsplit(data$`[30]consensus5AD`, ','), function(x) as.numeric(x)[2]))
    data$consensus5AltDepth <- ifelse(is.infinite(data$consensus5AltDepth), NA, data$consensus5AltDepth)

    # Fix gene names.
    data <- data %>% dplyr::mutate(`[17]ann_gene_name_Job` = gsub('\\|.*', '', `[17]ann_gene_name`))

    # Fix annotations.
    data <- data %>% dplyr::mutate(`[15]ann_annotation_Job` = gsub('\\|.*', '', `[15]ann_annotation`))

    # Fix gnoMAD AF.
    data <- data %>% dplyr::mutate(
        `[11]non_cancer_AF` = as.numeric(`[11]non_cancer_AF`),
        `[11]non_cancer_AF` = ifelse(is.na(`[11]non_cancer_AF`), 0, `[11]non_cancer_AF`)
    )

    # Remove duplicates (if any)
    data <- data %>% dplyr::distinct()

    # Return cleaned data.
    return(data)

}


# Import data -------------------------------------------------------------

biomicsIDs <- readr::read_tsv('Misc./biomicsIDs.txt')

# Import germline and somatic variants.
vcfFiles <- list.files('//mnt/data2/hartwig/DR71/Misc/CABAV7', pattern = 'tsv.gz$', full.names = T)
vcfFiles <- vcfFiles[grepl(paste(biomicsIDs$biomicsID, collapse = '|'), vcfFiles)]

# Import unfiltered QIAseq mutations.
unfilteredMuts <- dplyr::bind_rows(pbapply::pblapply(vcfFiles, readVCF.Biomics, cl = 10))
unfilteredMuts <- unfilteredMuts %>% dplyr::left_join(biomicsIDs)

# Only retain included samples.
data.Patients <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Sample Overview')
sampleLinks <- data.Patients %>%
    dplyr::select(`L-code (T1)`, `L-code (Germline)`, `Subject Number`, `AR-V7 (Baseline)`) %>%
    dplyr::mutate(hasGermline = (`L-code (T1)` != 'N/A' & `L-code (Germline)` != 'N/A'))

# Filter.
unfilteredMuts <- unfilteredMuts %>% dplyr::filter(!is.na(`L-code`), `L-code` %in% c(sampleLinks$`L-code (T1)`, sampleLinks$`L-code (Germline)`))

# Add subject #.
unfilteredMuts <- unfilteredMuts %>% dplyr::left_join(rbind(sampleLinks %>% dplyr::select(`L-code` = `L-code (T1)`, `Subject Number`, hasGermline), sampleLinks %>% dplyr::select(`L-code` = `L-code (Germline)`, `Subject Number`, hasGermline)))


# Generate PON for the missing germline samples ---------------------------

# Germline variants present in 5 or more sample are used within the PON for sample lacking matching germline data.
PON.CABA <- unfilteredMuts %>%
    dplyr::filter(
        Type == 'Germline',
        `[25]consensusDP` >= 5
    ) %>%
    dplyr::group_by(ID_Job) %>%
    dplyr::summarise(totalN = dplyr::n_distinct(`Subject Number`)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(totalN >= 5)


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


## Add germline status ----

determineGermline <- function(x, PON){

    # Select germline variants.
    muts.Germline <- x %>% dplyr::filter(Type == 'Germline') %>% dplyr::filter(`[25]consensusDP` >= 5) %>% dplyr::pull(ID_Job)

    # Germline status.
    germlineStatus <- 'Germline (Matched)'
    if(base::length(muts.Germline) == 0) germlineStatus <- 'Germline (PON)'

    # If sample had no germline, use the PON.
    if(base::length(muts.Germline) == 0) muts.Germline <- PON %>% dplyr::pull(ID_Job)



    # Filter germline variants (or PON)
    x <- x %>% dplyr::mutate(germlineStatus = ifelse(Type == 'Somatic' & !ID_Job %in% muts.Germline, 'Somatic', germlineStatus))

    return(x)

}

unfilteredMuts.perSubject <- base::split(unfilteredMuts, unfilteredMuts$`Subject Number`)

unfilteredMuts <- dplyr::bind_rows(pbapply::pblapply(unfilteredMuts.perSubject, determineGermline, PON = PON.CABA, cl = 10))


## Determine whether COSMIC has listed the mutations as SNP ----

cosmicSNPs <- readr::read_tsv('/mnt/data2/hartwig/DR71/Misc/cosmicSNPs.txt', col_names = 'ID', skip = 1)
unfilteredMuts$cosmicSNP <- unfilteredMuts$`[3]ID` %in% cosmicSNPs$ID


## Perform filtering. ----

filteredMuts <- unfilteredMuts %>%
    dplyr::filter(
        # Filter mutations within panel-design.
        withinPanel,
        # Filter on DP.
        `[29]consensus5DP` >= 100,
        # Filter on min. VAF.
        !is.na(`[28]consensus5AF_Job`),
        (`[28]consensus5AF_Job` >= 0.1 | (`[28]consensus5AF_Job` >= 0.05 & consensus5AltDepth >= 10)),
        # Filter on max. VAF.
        (`[28]consensus5AF_Job` <= 0.95 & !is.na(`[28]consensus5AF_Job`)),
        # Filter on coding.
        `[18]ann_hgvs_p` != '.',
        # Filter germline variants.
        germlineStatus == 'Somatic',
        # Remove gnomad non-cancer variants.
        `[11]non_cancer_AF` <= 0.025,
        # Remove COSMIC SNPs.
        !cosmicSNP,
        # Remove repeat-like mis-alignments.
        !grepl('conservative_inframe_insertion', `[15]ann_annotation`)
    ) %>%
    # Filter on incidence within CABA-V7 samples.
    dplyr::group_by(ID_Job) %>%
    dplyr::mutate(totalCaba = dplyr::n_distinct(`L-code`)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(totalCaba <= 30)


# Output mutations --------------------------------------------------------

write.table(filteredMuts, file = 'filteredMuts.txt', quote = F, sep = '\t', row.names = F)
