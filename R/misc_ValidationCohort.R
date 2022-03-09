
library(plyr)
library(dplyr)
library(ggplot2)
library(tableone)
library(patchwork)

data.Patient <- list()
data.Patient$Overview <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Sample Overview') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))

# Selection of validation cohort ----

data.Patient$Overview %>%
    dplyr::group_by(`AR-V7 (Baseline)`) %>%
    dplyr::summarise(
        q1CTC = summary(`CTC Count (Baseline – 7.5mL)`)[2],
        mCTC = median(`CTC Count (Baseline – 7.5mL)`),
        q3CTC = summary(`CTC Count (Baseline – 7.5mL)`)[5],
        total = dplyr::n()
    )


newCohort <- readxl::read_xlsx('Misc./newCohort.xlsx') %>%
    dplyr::filter(`Plasma ml/ampul` >= .5) %>%
    dplyr::filter(!is.na(CTC), !is.na(`Plasma ampul`), `Plasma ampul` >= 1, !is.na(`AR-V7`)) %>%
    dplyr::mutate(ARV7 = ifelse(grepl('No res', `AR-V7`), 'Und.', `AR-V7`))

pats <- list()

#pats$und <- newCohort[newCohort$ARV7 == 'Und.',] %>% dplyr::filter(CTC < 15) %>%  dplyr::sample_n(20)
#pats$pos <- newCohort[newCohort$ARV7 == 'AR-V7 positive',] %>% dplyr::filter(CTC < 1500) %>% dplyr::sample_n(23)
#pats$neg <- newCohort[newCohort$ARV7 == 'AR-V7 negative',] %>% dplyr::filter(CTC < 500) %>% dplyr::sample_n(37)

pats$neg %>%
    dplyr::summarise(
        q1CTC = summary(CTC)[2],
        mCTC = median(CTC),
        q3CTC = summary(CTC)[5],
        total = dplyr::n()
    )

wilcox.test(data.Patient$Overview %>% dplyr::filter(`AR-V7 (Baseline)` == 'Und.') %>% dplyr::pull(`CTC Count (Baseline – 7.5mL)`), pats$und$CTC)
wilcox.test(data.Patient$Overview %>% dplyr::filter(`AR-V7 (Baseline)` == 'Pos.') %>% dplyr::pull(`CTC Count (Baseline – 7.5mL)`), pats$pos$CTC)
wilcox.test(data.Patient$Overview %>% dplyr::filter(`AR-V7 (Baseline)` == 'Neg.') %>% dplyr::pull(`CTC Count (Baseline – 7.5mL)`), pats$neg$CTC)