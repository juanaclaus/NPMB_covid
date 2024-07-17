##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##Script for Microbiome analysis of BCG-Corona and MB-COVID patients
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##this is for the clean data
##A lot of data formatting and creation of new variables was done prior
##Examples of visualizations are shown

####load packages####
library(reshape2)
library(stringi)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(openxlsx)
#library(phyloseq)
library(plotly)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(patchwork)
library(gtools) 
library(gdata)
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
#BiocManager::install("ANCOMBC")
#BiocManager::install("microbiomeutilities")
#BiocManager::install("microbiome")
library(ANCOMBC) 
library(microbiome)
library(microbiomeutilities)
library(readr)
library(lattice)
library(vegan)  
library(tidyverse)
library(hrbrthemes)
library(viridisLite)
library(viridis)
library(knitr)
library(janitor)
library(ggpubr)
library(ggsci)
library(readxl)
library(writexl)
library(BiocManager)
library(Biostrings)
library(mia)
library(phyloseq)
library(miaViz) #BiocManager::install("miaViz")
library(scater)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##MIA creates the TreeSummarized experiments
#generic and highly optimized container for complex data structures
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####loading data####

samdf<-read_excel("231113_clinical_data_merged.xlsx")
#row names of sample data have to be the sample names (from sequencing)
rownames(samdf)
samdf<-samdf %>% remove_rownames %>% column_to_rownames(var="sample_name")
sample_names<-rownames(samdf)

counts<-read_excel("231120_ASV_BLAST_BCG-MB_merged.xlsx")
rownames(counts)
counts<-counts %>% remove_rownames %>% column_to_rownames(var="OTUID")
counts<-counts[order(row.names(counts)), ]

all(rownames(samdf) == colnames(counts)) #should be TRUE

taxa<-read_excel("231130_taxonomy_BLAST_BCG-MB_merged.xlsx")
rownames(taxa)
taxa<-taxa %>% remove_rownames %>% column_to_rownames(var="OTUID")
taxa_rownames<-rownames(taxa)
#row name for taxa should be ASV
#have to be in same order as count
taxa<-taxa[order(row.names(taxa)), ]

# rowdata rownames match assay rownames
all(rownames(taxa) == rownames(counts)) # should be TRUE

# Create a list that contains assays
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counts <- as.matrix(counts)
assays <- SimpleList(counts = counts)

# Convert colData and rowData into DataFrame
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samdf <- DataFrame(samdf)
class(samdf)
taxa <- DataFrame(taxa)
class(taxa)

#### Create TreeSE ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tse <- TreeSummarizedExperiment(assays = assays,
                                colData = samdf,
                                rowData = taxa
)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Remove mock sample 
tse <- tse[ , colnames(tse) != "mock"]

# Convert rownames into ASV_number format
rownames(tse) <- paste0("ASV", seq( nrow(tse) ))
tse

tree_genus <- agglomerateByRank(tse, rank = "Genus", onRankOnly=TRUE)
nrow(tree_genus)
dim(tree_genus)
# 413 taxa and 400 samples
###relative abundance
tree_genus<-transformAssay(tree_genus, method = "relabundance")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###Create BCG tse and remove the singletons ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filtered_tree_genus<-tree_genus
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##only for BCG-C
tse_genus_BCG<-filtered_tree_genus[, filtered_tree_genus$study %in% c("BCG-Corona")]

# Check taxa abundance and identify taxa with zero reads and remove
taxa_abundance_bcg <- rowSums(assay(tse_genus_BCG, "counts"))
zero_read_taxa_bcg <- names(taxa_abundance_bcg[taxa_abundance_bcg == 0])
non_zero_read_taxa_bcg <-names(taxa_abundance_bcg[!(taxa_abundance_bcg == 0)])
tse_genus_BCG <- tse_genus_BCG[rowData(tse_genus_BCG)$Genus %in% non_zero_read_taxa_bcg, ]

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###now remove the singletons
tse_genus_BCG_single<-tse_genus_BCG

# Step 1: Extract abundance data
abundance_matrix <- assay(tse_genus_BCG_single)

# Step 2: Identify features present in only one sample
singleton_features <- rownames(abundance_matrix)[rowSums(abundance_matrix > 0) == 1]
singleton_features<-base::as.data.frame(singleton_features)

non_singletons <-rownames(abundance_matrix)[!(rowSums(abundance_matrix > 0) == 1)]

#Filter them out
tse_genus_BCG_single <- tse_genus_BCG_single[rowData(tse_genus_BCG_single)$Genus %in% non_singletons, ]
tse_genus_BCG_filtered<-tse_genus_BCG_single

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###Create MB tse and remove the singletons ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tse_genus_MB<-filtered_tree_genus[, filtered_tree_genus$study %in% c("MB-COVID")]

# Check taxa abundance and identify taxa with zero reads and remove
taxa_abundance_mb <- rowSums(assay(tse_genus_MB, "counts"))
zero_read_taxa_mb <- names(taxa_abundance_mb[taxa_abundance_mb == 0])
non_zero_read_taxa_mb <-names(taxa_abundance_mb[!(taxa_abundance_mb == 0)])
tse_genus_MB <- tse_genus_MB[rowData(tse_genus_MB)$Genus %in% non_zero_read_taxa_mb, ]

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###now remove the singletons
tse_genus_MB_single<-tse_genus_MB


# Step 1: Extract abundance data
abundance_matrix_mb <- assay(tse_genus_MB_single)

# Step 2: Identify features present in only one sample
singleton_features_mb <- rownames(abundance_matrix_mb)[rowSums(abundance_matrix_mb > 0) == 1]
singleton_features_mb<-base::as.data.frame(singleton_features_mb)
#write.xlsx(singleton_features_mb, "~/surfdrive/PhD/MMB/Merge BCG-MB/ANALYSIS-Sept 2023/Controls/singletons/231127_MB-COVID_genus_singleton.xlsx")

non_singleton_mbs <-rownames(abundance_matrix_mb)[!(rowSums(abundance_matrix_mb > 0) == 1)]

#Filter them out
tse_genus_MB_single <- tse_genus_MB_single[rowData(tse_genus_MB_single)$Genus %in% non_singleton_mbs, ]
tse_genus_MB_filtered<-tse_genus_MB_single

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####top 10 abundance-MB-COVID participants ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##143 samples

##by infection type genus ##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tse_genus_MB_10<-tse_genus_MB_filtered
top_taxa <- getTopFeatures(tse_genus_MB_10,top = 10, assay.type = "relabundance")
rowData(tse_genus_MB_10)[top_taxa, taxonomyRanks(tse_genus_MB_10)]
genus_renamed <- lapply(rowData(tse_genus_MB_10)$Genus,
                        function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_genus_MB_10)$Genus <- as.character(genus_renamed)
colData(tse_genus_MB_10)$covid_infection<-factor(colData(tse_genus_MB_10)$covid_infection,
                                                 levels = c("covid", "other"))

##sample code of visualizations
plot<-plotAbundance(
  tse_genus_MB_10,
  rank = "Genus",
  features = "covid_infection",
  order_rank_by = "name",
  order_sample_by = "covid_infection",
  decreasing = TRUE,
  use_relative = TRUE,
  layout = "bar"
)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####top 10 abundance-BCG-Corona participants ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##by infection type genus ##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tse_genus_BCG_10<-tse_genus_BCG_filtered

## Compositional barplot with top 10 taxa 
# Renaming the "genus" rank to keep only top taxa and the rest to "Other"

top_taxa <- getTopFeatures(tse_genus_BCG_10,top = 10, assay.type = "relabundance")
rowData(tse_genus_BCG_10)[top_taxa, taxonomyRanks(tse_genus_BCG_10)]
genus_renamed <- lapply(rowData(tse_genus_BCG_10)$Genus,
                        function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_genus_BCG_10)$Genus <- as.character(genus_renamed)
colData(tse_genus_BCG_10)$covid_infection<-factor(colData(tse_genus_BCG_10)$covid_infection,
                                                  levels = c("covid", "other", "No"))

##sample code of visualizations
plot<-plotAbundance(
  tse_genus_BCG_10,
  rank = "Genus",
  features = "covid_infection",
  order_rank_by = "name",
  order_sample_by = "covid_infection",
  decreasing = TRUE,
  use_relative = TRUE,
  layout = "bar"
)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####Alpha diversities MB-COVID#### 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tse_genus_MB_filtered <- mia::estimateDiversity(tse_genus_MB_filtered, 
                                                assay.type = "counts",
                                                index = c("shannon","inverse_simpson", "coverage"), 
                                                name = c("shannon","inverse_simpson", "coverage") )


alpha_mb<-as.data.frame(colData(tse_genus_MB_filtered))
alpha_mb_co<-alpha_mb%>%mutate(coinfection=covid_indicator)
alpha_mb_co%>%tabyl(coinfection)
alpha_mb_co%>%tabyl(coinfection, other_indicator)

alpha_mb_co$coinfection<-as.character(alpha_mb_co$coinfection)
alpha_mb_co<-alpha_mb_co%>%mutate(coinfection=ifelse(
  coinfection=="1" & other_indicator=="1", "co", coinfection
))
alpha_mb_co<-alpha_mb_co[!(alpha_mb_co$coinfection=="co"),]
alpha_mb_co_cf<-alpha_mb_co[!(alpha_mb_co$cystic_fibrosis=="YES"),]

#tables
wilcox.test(shannon ~ covid_infection, data = colData(tse_genus_MB_filtered))
wilcox.test(inverse_simpson ~ covid_infection, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ severity_ab, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ severity_ab, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ antibiotic_users, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ antibiotic_users, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ fever, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ fever, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ respiratory_symptoms, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ respiratory_symptoms, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ Dyspnea, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ Dyspnea, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ non_resp, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ non_resp, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ covid_indicator, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ covid_indicator, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ other_indicator, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ other_indicator, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ comorbidity_summary_cat, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ comorbidity_summary_cat, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ CT_covid_cat, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ CT_covid_cat, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ in_out_patient, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ in_out_patient, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ age_cat, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ age_cat, data = colData(tse_genus_MB_filtered))

kruskal.test(shannon ~ sex, data = colData(tse_genus_MB_filtered))
kruskal.test(inverse_simpson ~ sex, data = colData(tse_genus_MB_filtered))

#tables mean and sd

alpha_mb %>%
  group_by(covid_infection) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(severity_ab) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(antibiotic_users) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(fever) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(respiratory_symptoms) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(Dyspnea) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(non_resp) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(covid_indicator) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(other_indicator) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(comorbidity_summary_cat) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(CT_covid_cat) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(in_out_patient) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(age_cat) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_mb %>%
  group_by(sex) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)

alpha_mb_co %>%
  group_by(covid_infection) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
kruskal.test(shannon ~ covid_infection, data = alpha_mb_co)
kruskal.test(inverse_simpson ~ covid_infection, data = alpha_mb_co)


alpha_mb_co %>%
  group_by(severity_ab) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
kruskal.test(shannon ~ severity_ab, data = alpha_mb_co)
kruskal.test(inverse_simpson ~ severity_ab, data = alpha_mb_co)


#linear regression crude
tabyl(alpha_mb_co$covid_infection)
alpha_mb_co_severity<-alpha_mb_co[!(alpha_mb_co$covid_infection=="other"),]
tabyl(alpha_mb_co_severity$WHO_severity)

lmshannon <- lm(shannon~covid_infection, data = alpha_mb_co_cf)
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~WHO_severity, data = alpha_mb_co_severity) 
summary(lmshannon) #Review the results
confint(lmshannon)


lmshannon <- lm(shannon~comorbidity_summary_cat, data = alpha_mb_co) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~CT_covid_cat, data = alpha_mb_co) 
summary(lmshannon) #Review the results
confint(lmshannon)


lmsimp <- lm(inverse_simpson~covid_infection, data = alpha_mb_co_cf) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~WHO_severity, data = alpha_mb_co_severity) 
summary(lmsimp) #Review the results
confint(lmsimp)


lmsimp <- lm(inverse_simpson~comorbidity_summary_cat, data = alpha_mb_co) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~CT_covid_cat, data = alpha_mb_co) 
summary(lmsimp) #Review the results
confint(lmsimp)

#linear regression adjusted
colData(tse_genus_MB_filtered)$age_yrs<-as.numeric(colData(tse_genus_MB_filtered)$age_yrs)
colData(tse_genus_MB_filtered)$comorbidity_summary<-as.numeric(colData(tse_genus_MB_filtered)$comorbidity_summary)

lmshannon <- lm(shannon~covid_infection +sex +age_yrs +comorbidity_summary+antibiotic_users, data = alpha_mb_co_cf) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~WHO_severity+sex +age_yrs +comorbidity_summary+antibiotic_users, data = alpha_mb_co_severity) 
summary(lmshannon) #Review the results
confint(lmshannon)


lmshannon <- lm(shannon~CT_covid_cat+sex +age_yrs +comorbidity_summary+antibiotic_users, data = alpha_mb_co) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~comorbidity_summary_cat+sex +age_yrs, data = alpha_mb_co) 
summary(lmshannon) #Review the results
confint(lmshannon)


lmsimp <- lm(inverse_simpson~covid_infection+sex +
               age_yrs +comorbidity_summary+antibiotic_users, data = alpha_mb_co_cf) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~WHO_severity+sex +
               age_yrs +comorbidity_summary+antibiotic_users, data = alpha_mb_co_severity) 
summary(lmsimp) #Review the results
confint(lmsimp)


lmsimp <- lm(inverse_simpson~CT_covid_cat+sex +age_yrs +comorbidity_summary+antibiotic_users, data = alpha_mb_co) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~comorbidity_summary_cat+sex +age_yrs, data = alpha_mb_co) 
summary(lmsimp) #Review the results
confint(lmsimp)

##sample code of visualizations
ggplot(as.data.frame(colData(tse_genus_MB_filtered_ab)),
       aes(x = covid_infection, 
           y = inverse_simpson, 
           fill=covid_infection)) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA, width = 0.50)+
  geom_jitter(aes(color=covid_infection), width = 0.30)+
  #stat_compare_means(label = "p.format",label.y = 3.1)+
  stat_compare_means(label.x.npc = 0.01, size = 5, label="p.format")+
  theme_classic()+
  theme(legend.position = "none")+
  labs(
    y="Inverse Simpson")+
  ggtitle(bquote(~ bold('e)') ~'Alpha diversity by virus type regardless of'~symptoms^"2,3,4"))+
  #ggtitle(bquote('E) Alpha diversity by virus type regardless of'~symptoms^"2,3,4"))+
  theme(axis.text.x = element_text( size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        title = element_text(size = 15),
        axis.text.y = element_text(size=15), 
        legend.text =element_text(size = 15) )+
  scale_fill_lancet()+
  scale_colour_lancet()+
  scale_x_discrete(
    labels=c("SARS-CoV-2 \n n=85", "Other virus \n n=30"))+
  stat_summary(geom="pointrange", color="black")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####Alpha diversities BCG-Corona#### 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tse_genus_BCG_filtered <- mia::estimateDiversity(tse_genus_BCG_filtered, 
                                                 assay.type = "counts",
                                                 index = c("shannon","inverse_simpson", "coverage"), 
                                                 name = c("shannon","inverse_simpson", "coverage") )


#tables
kruskal.test(shannon ~ covid_infection, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ covid_infection, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ WHO_severity, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ WHO_severity, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ overall_severity, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ overall_severity, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ symptomatic_swab, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ symptomatic_swab, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ hospital, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ hospital, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ antibiotic_users, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ antibiotic_users, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ covid_indicator, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ covid_indicator, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ other_indicator, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ other_indicator, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ comorbidity_summary_cat, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ comorbidity_summary_cat, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ sex, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ sex, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ age_cat, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ age_cat, data = colData(tse_genus_BCG_filtered))

kruskal.test(shannon ~ severity, data = colData(tse_genus_BCG_filtered))
kruskal.test(inverse_simpson ~ severity, data = colData(tse_genus_BCG_filtered))

#tables mean and sd
alpha_bcg<-as.data.frame(colData(tse_genus_BCG_filtered))
alpha_bcg_severity<-alpha_bcg[(alpha_bcg$covid_infection)=="covid",]
alpha_bcg_severity<-alpha_bcg_severity[!(alpha_bcg_severity$virus_other=="Adenovirus"),]
alpha_bcg_severity<-alpha_bcg_severity[!(alpha_bcg_severity$virus_other=="BOCA"),]
alpha_bcg_severity<-alpha_bcg_severity[!(alpha_bcg_severity$virus_other=="Rhino"),]



alpha_bcg_covid<-alpha_bcg[!(alpha_bcg$covid_infection=="other"),]
alpha_bcg_covid<-alpha_bcg_covid[!(alpha_bcg_covid$virus_other=="Adenovirus"),]
alpha_bcg_covid<-alpha_bcg_covid[!(alpha_bcg_covid$virus_other=="BOCA"),]
alpha_bcg_covid<-alpha_bcg_covid[!(alpha_bcg_covid$virus_other=="Rhino"),]

alpha_bcg_other<-alpha_bcg[!(alpha_bcg$covid_infection=="covid"),]

alpha_bcg_covid_other<-alpha_bcg[!(alpha_bcg$covid_infection=="No"),]
alpha_bcg_covid_other<-alpha_bcg_covid_other%>%mutate(coinfection=covid_indicator)
alpha_bcg_covid_other%>%tabyl(coinfection)
alpha_bcg_covid_other%>%tabyl(coinfection, other_indicator)

alpha_bcg_covid_other$coinfection<-as.character(alpha_bcg_covid_other$coinfection)
alpha_bcg_covid_other<-alpha_bcg_covid_other%>%mutate(coinfection=ifelse(
  coinfection=="1" & other_indicator=="1", "co", coinfection
))
alpha_bcg_covid_other<-alpha_bcg_covid_other[!(alpha_bcg_covid_other$coinfection=="co"),]

alpha_bcg_covid
alpha_bcg_covid$covid_infection<-factor(alpha_bcg_covid$covid_infection, levels=c("No", "covid"))
alpha_bcg_other
alpha_bcg_other$covid_infection<-factor(alpha_bcg_other$covid_infection, levels=c( "No", "other"))
alpha_bcg_covid_other
alpha_bcg_covid_other$covid_infection<-factor(alpha_bcg_covid_other$covid_infection, levels=c("other" ,"covid"))

alpha_bcg %>%
  group_by(covid_infection) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(WHO_severity) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(overall_severity) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(symptomatic_swab) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(hospital) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(antibiotic_users) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(covid_indicator) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(other_indicator) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(comorbidity_summary_cat) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(age_cat) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
alpha_bcg %>%
  group_by(sex) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)

alpha_bcg %>%
  group_by(severity) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)

alpha_bcg_severity %>%
  group_by(severity) %>%
  summarise_at(vars(shannon, inverse_simpson), list(mean,sd), na.rm = TRUE)
kruskal.test(shannon ~ severity, data = alpha_bcg_severity)
kruskal.test(inverse_simpson ~ severity, data = alpha_bcg_severity)

#linear regression crude
lmshannon <- lm(shannon~covid_infection, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~WHO_severity, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~overall_severity, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~symptomatic_swab, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~covid_indicator, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~other_indicator, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~severity, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

alpha_bcg_severity<-alpha_bcg[(alpha_bcg$covid_infection)=="covid",]

lmshannon <- lm(shannon~severity, data = alpha_bcg_severity) 
summary(lmshannon) #Review the results
confint(lmshannon)

alpha_bcg_covid
alpha_bcg_other
alpha_bcg_covid_other

lmshannon <- lm(shannon~covid_infection, data = alpha_bcg_covid) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~covid_infection, data = alpha_bcg_covid_other) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~covid_infection, data = alpha_bcg_other) 
summary(lmshannon) #Review the results
confint(lmshannon)



lmsimp <- lm(inverse_simpson~covid_infection, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~WHO_severity, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~overall_severity, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~symptomatic_swab, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~covid_indicator, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~other_indicator, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~severity, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)
lmsimp <- lm(inverse_simpson~severity, data = alpha_bcg_severity) 
summary(lmsimp) #Review the results
confint(lmsimp)


lmsimp <- lm(inverse_simpson~covid_infection, data = alpha_bcg_covid) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~covid_infection, data = alpha_bcg_covid_other) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~covid_infection, data = alpha_bcg_other) 
summary(lmsimp) #Review the results
confint(lmsimp)


#linear regression adjusted
colData(tse_genus_BCG_filtered)$age_yrs<-as.numeric(colData(tse_genus_BCG_filtered)$age_yrs)

lmshannon <- lm(shannon~covid_infection +sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~WHO_severity+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~overall_severity+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~symptomatic_swab+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~covid_indicator+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~other_indicator+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~severity+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~severity+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data = alpha_bcg_severity) 
summary(lmshannon) #Review the results
confint(lmshannon)


lmshannon <- lm(shannon~covid_infection+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data =alpha_bcg_covid) 
summary(lmshannon) #Review the results
confint(lmshannon)

lmshannon <- lm(shannon~covid_infection+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data =alpha_bcg_covid_other) 
summary(lmshannon) #Review the results
confint(lmshannon)


lmshannon <- lm(shannon~covid_infection+sex +age_yrs +comorbidity_summary+
                  antibiotic_users, data  = alpha_bcg_other) 
summary(lmshannon) #Review the results
confint(lmshannon)




lmsimp <- lm(inverse_simpson~covid_infection+sex +age_yrs +comorbidity_summary, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~WHO_severity+sex +age_yrs +comorbidity_summary, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~overall_severity+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~symptomatic_swab+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~covid_indicator+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~other_indicator+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~severity+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data = colData(tse_genus_BCG_filtered)) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~severity+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data = alpha_bcg_severity) 
summary(lmsimp) #Review the results
confint(lmsimp)


lmsimp <- lm(inverse_simpson~covid_infection+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data =alpha_bcg_covid) 
summary(lmsimp) #Review the results
confint(lmsimp)

lmsimp <- lm(inverse_simpson~covid_infection+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data =alpha_bcg_covid_other) 
summary(lmsimp) #Review the results
confint(lmsimp)


lmsimp <- lm(inverse_simpson~covid_infection+sex +age_yrs +comorbidity_summary+
               antibiotic_users, data  = alpha_bcg_other) 
summary(lmsimp) #Review the results
confint(lmsimp)


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####Beta diversites MB-COVID#### 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Estimate beta diversity
tse_genus_MB_filtered <- runMDS(tse_genus_MB_filtered, FUN = vegan::vegdist,
                                name = "MDS_BC", exprs_values = "counts")

# Create ggplot object
p <- plotReducedDim(tse_genus_MB_filtered, "MDS_BC", colour_by = "severity_ab")

# Add explained variance for each axis
e <- attr(reducedDim(tse_genus_MB_filtered, "MDS_BC"), "eig");
rel_eig <- e/sum(e[e>0])          
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""),
              color="Severity")

p <- p + scale_color_lancet(labels=c("Other viral infection", "Mild", "Moderate", 
                                     "Severe", "Fatal", "SARS-CoV-2 AB",
                                     "Other Virus AB/CF"))
print(p)

###Permanova for testing differences 
permanova <- vegan::adonis2(t(assay(tse_genus_MB_filtered,"relabundance")) ~ severity_ab,
                            data = colData(tse_genus_MB_filtered),
                            permutations = 9999)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####Beta diversites BCG-Corona#### 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Estimate beta diversity
tse_genus_BCG_filtered <- runMDS(tse_genus_BCG_filtered, FUN = vegan::vegdist, name = "MDS_BC", exprs_values = "counts",
                                 method = "bray")


vegan::adonis2(t(assay(tse_genus_BCG_filtered_alpha,"relabundance")) ~ covid_infection,
               data = colData(tse_genus_BCG_filtered_alpha),
               permutations = 9999)

plot_pathogen_b <- plotReducedDim(tse_genus_BCG_filtered_alpha, "MDS_BC", colour_by = "covid_infection",
                                  point_size=3)
tabyl(colData(tse_genus_BCG_filtered)$covid_infection)
plot_pathogen_b <- plot_pathogen_b +
  scale_color_lancet(limits=c("No", "covid", "other") ,
                     labels=c("None detected","SARS-CoV-2", "Other virus"))

plot_pathogen_b <- plot_pathogen_b +labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
                                         y = paste("PCoA 2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""),
                                         
                                         color="",
                                         subtitle = "Pr(>F) = 0.16")

plot_pathogen_b <- plot_pathogen_b+ 
  ggtitle(bquote(~ bold('g)') ~'MDS Bray-Curtis by virus type regardless of'~symptoms^"2,6"))
# ggtitle(bquote('G) MDS Bray-Curtis by virus type regardless of'~symptoms^"2,6"))

plot_pathogen_b<-plot_pathogen_b+theme_classic()

plot_pathogen_b<-plot_pathogen_b +theme(axis.title.x = element_text(size=15),
                                        legend.text = element_text(size=15),
                                        axis.title.y = element_text(size=15),
                                        strip.text.x = element_text(size = 15),
                                        title = element_text(size=15), 
                                        axis.text = element_text(size=15))


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####ANCOMBC2 MB-COVID####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#https://microbiome.github.io/OMA/differential-abundance.html#daa-with-confounding
#correct for confounding 

#first run the functions by Bart
phy_mb_an<-makePhyloseqFromTreeSummarizedExperiment(tse_genus_MB_filtered) 


#Here, confounders can be added to the formula along withthe main outcome variable.
#This way, the model evaluates whether differentially abundant taxa are associated 
#with one of the variables when 'age_cat' , 'sex', 'antibiotic_users' are kept constant.


sample_data(phy_mb_an)$age_cat<-factor(sample_data(phy_mb_an)$age_cat,
                                       levels=c("≤50","50-70" ,"≥70"))

sample_data(phy_mb_an)$sex<-factor(sample_data(phy_mb_an)$sex,
                                   levels=c("MALE","FEMALE"))

sample_data(phy_mb_an)$comorbidity_summary_cat
sample_data(phy_mb_an)$antibiotic_users


#pathogen type
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mb_pANC_run_covid <- anc_run(phylo_obj = phy_mb_an, 
                             var = c('covid_infection', 'age_cat' , 'sex', 'antibiotic_users', "comorbidity_summary_cat"),
                             var_levels = list(c( 'other','covid'), 
                                               c("≤50","50-70" ,"≥70"),
                                               c("MALE","FEMALE"),
                                               c("NO", "YES", "possible"),
                                               c("0", "1-5", "6-10", "11-15")),
                             level = 'Genus', 
                             pv_cut = 0.10,
                             confounder = NULL,
                             struc_zero = T,             
                             neg_lb = T,                 
                             pairwise = F,               
                             dunnett = T,                
                             global = T,
                             trend = F,
                             padj_method = 'BH',
                             padj_m = 'mdFDR',
                             subset_column = NULL,         
                             subset_value = NULL)


pANC_data_covid <- anc_data(mb_pANC_run_covid,                   # Output from anc_run
                            var = 'covid_infection',
                            var_levels = c( 'other','covid'),                  
                            val = 'padj',                 
                            sig_cut_off = 0.05,            
                            lfc_cut_off = 0,               
                            test = 'primary')


#structural zeros
tab_zero = mb_pANC_run_covid$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####ANCOMBC2 BCG-Corona####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phy_bcg_an<-makePhyloseqFromTreeSummarizedExperiment(tse_genus_BCG_filtered) 
sample_data(phy_bcg_an)$age_cat<-factor(sample_data(phy_bcg_an)$age_cat,
                                        levels=c("≤50","50-70" ,"≥70"))

sample_data(phy_bcg_an)$sex<-factor(sample_data(phy_bcg_an)$sex,
                                    levels=c("MALE","FEMALE"))

sample_data(phy_bcg_an)$overall_severity
sample_data(phy_bcg_an)$symptomatic_swab
sample_data(phy_bcg_an)$comorbidity_summary

pANC_run_covid_aj <- anc_run(phylo_obj = phy_bcg_an, 
                             var = c('covid_infection', 'age_cat' , 'sex', 'antibiotic_users'),
                             var_levels = list(c('No', 'other','covid'), 
                                               c("≤50","50-70"),
                                               c("MALE","FEMALE"),
                                               c("NO",  "possible")),
                             level = 'Genus', 
                             pv_cut = 0.10,
                             confounder = NULL,
                             struc_zero = T,             
                             neg_lb = T,                 
                             pairwise = F,               
                             dunnett = T,                
                             global = T,
                             trend = F,
                             padj_method = 'BH',
                             padj_m = 'mdFDR',
                             subset_column = NULL,         
                             subset_value = NULL)


#structural zeros
tab_zero = pANC_run_covid_aj$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

pANC_data_covid <- anc_data(pANC_run_covid_aj,                   # Output from anc_run
                            var = 'covid_infection',
                            var_levels = c('No', 'other','covid'),                  
                            val = 'padj',                 
                            sig_cut_off = 0.05,            
                            lfc_cut_off = 0,               
                            test = 'primary')



##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####    Linear regression models     ####
#      EC concentrations of the bacterial groups
#      Removing the co-infections and the 5CF patients 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(gtsummary)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(janitor)
library(lubridate) #to work with dates
library(linelist)
library(tidyr)
library(rstatix)
library(flextable)
library(reshape2)
library(tableone)
library(nnet)
library(broom)
library(ggsci)
library(cowplot)
library(ggvis)
library(gridExtra)
library(devtools)



mem_bcg<-read.xlsx("~/surfdrive/PhD/MMB/Merge BCG-MB/ANALYSIS-Sept 2023/Blasted databases/Mixed_models_dataset/240216_EC_BCG_bacterial_groups.xlsx")
#row names of sample data have to be the sample names (from sequencing)
rownames(mem_bcg)
mem_bcg<-mem_bcg %>% remove_rownames %>% column_to_rownames(var="X1")

mem_bcg$g__Moraxella<-as.numeric(mem_bcg$g__Moraxella)
mem_bcg$g__Corynebacterium<-as.numeric(mem_bcg$g__Corynebacterium)
mem_bcg$g__Dolosigranulum<-as.numeric(mem_bcg$g__Dolosigranulum)
mem_bcg$g__Streptococcus<-as.numeric(mem_bcg$g__Streptococcus)
mem_bcg$g__Staphylococcus<-as.numeric(mem_bcg$g__Staphylococcus)
mem_bcg$g__Actinomyces<-as.numeric(mem_bcg$g__Actinomyces)
mem_bcg$g__Haemophilus<-as.numeric(mem_bcg$g__Haemophilus)
mem_bcg$g__Neisseria<-as.numeric(mem_bcg$g__Neisseria)
mem_bcg$g__Pseudomonas<-as.numeric(mem_bcg$g__Pseudomonas)
mem_bcg$oral_cavity_bacteria<-as.numeric(mem_bcg$oral_cavity_bacteria)
mem_bcg$other.minority.species<-as.numeric(mem_bcg$other.minority.species)


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##loading estimated concentations 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mem_bcg<-read.xlsx("240216_EC_BCG_bacterial_groups.xlsx")
#row names of sample data have to be the sample names (from sequencing)
rownames(mem_bcg)
mem_bcg<-mem_bcg %>% remove_rownames %>% column_to_rownames(var="X1")

mem_bcg$g__Moraxella<-as.numeric(mem_bcg$g__Moraxella)
mem_bcg$g__Corynebacterium<-as.numeric(mem_bcg$g__Corynebacterium)
mem_bcg$g__Dolosigranulum<-as.numeric(mem_bcg$g__Dolosigranulum)
mem_bcg$g__Streptococcus<-as.numeric(mem_bcg$g__Streptococcus)
mem_bcg$g__Staphylococcus<-as.numeric(mem_bcg$g__Staphylococcus)
mem_bcg$g__Actinomyces<-as.numeric(mem_bcg$g__Actinomyces)
mem_bcg$g__Haemophilus<-as.numeric(mem_bcg$g__Haemophilus)
mem_bcg$g__Neisseria<-as.numeric(mem_bcg$g__Neisseria)
mem_bcg$g__Pseudomonas<-as.numeric(mem_bcg$g__Pseudomonas)
mem_bcg$oral_cavity_bacteria<-as.numeric(mem_bcg$oral_cavity_bacteria)
mem_bcg$other_minority_species<-as.numeric(mem_bcg$other_minority_species)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mem_mb<-read.xlsx("240130_EC_bacterial_groups_MB.xlsx")
#row names of sample data have to be the sample names (from sequencing)
rownames(mem_mb)
mem_mb<-mem_mb %>% remove_rownames %>% column_to_rownames(var="X1")

mem_mb$g__Moraxella<-as.numeric(mem_mb$g__Moraxella)
mem_mb$g__Corynebacterium<-as.numeric(mem_mb$g__Corynebacterium)
mem_mb$g__Dolosigranulum<-as.numeric(mem_mb$g__Dolosigranulum)
mem_mb$g__Streptococcus<-as.numeric(mem_mb$g__Streptococcus)
mem_mb$g__Staphylococcus<-as.numeric(mem_mb$g__Staphylococcus)
mem_mb$g__Actinomyces<-as.numeric(mem_mb$g__Actinomyces)
mem_mb$g__Haemophilus<-as.numeric(mem_mb$g__Haemophilus)
mem_mb$g__Neisseria<-as.numeric(mem_mb$g__Neisseria)
mem_mb$g__Pseudomonas<-as.numeric(mem_mb$g__Pseudomonas)
mem_mb$oral_cavity_bacteria<-as.numeric(mem_mb$oral_cavity_bacteria)
mem_mb$other_minority_species<-as.numeric(mem_mb$other_minority_species)



MB_data<-as.data.frame(colData(tse_genus_MB_filtered))

MB_data_EC<-merge(MB_data,mem_mb, by = 'row.names', all = TRUE )

MB_data_EC$Moraxella_log<-log10(MB_data_EC$g__Moraxella+1)
MB_data_EC$Streptococcus_log<-log10(MB_data_EC$g__Streptococcus+1)
MB_data_EC$Corynebacterium_log<-log10(MB_data_EC$g__Corynebacterium+1)
MB_data_EC$Dolosigranulum_log<-log10(MB_data_EC$g__Dolosigranulum+1)
MB_data_EC$Staphylococcus_log<-log10(MB_data_EC$g__Staphylococcus+1)
MB_data_EC$Actinomyces_log<-log10(MB_data_EC$g__Actinomyces+1)
MB_data_EC$Streptococcus_log<-log10(MB_data_EC$g__Streptococcus+1)
MB_data_EC$Haemophilus_log<-log10(MB_data_EC$g__Haemophilus+1)
MB_data_EC$Neisseria_log<-log10(MB_data_EC$g__Neisseria+1)
MB_data_EC$Pseudomonas_log<-log10(MB_data_EC$g__Pseudomonas+1)
MB_data_EC$oral_cavity_bacteria_log<-log10(MB_data_EC$oral_cavity_bacteria+1)
MB_data_EC$other_minority_species<-log10(MB_data_EC$other_minority_species+1)

##Filtering and removing 2 co-infections and 5CF cases
MB_data_EC$cystic_fibrosis
MB_data_EC_cf<-MB_data_EC[MB_data_EC$cystic_fibrosis=="NO",]

MB_data_EC_cf<-MB_data_EC_cf%>%mutate(coinfection=covid_indicator)
MB_data_EC_cf%>%tabyl(coinfection)
MB_data_EC_cf%>%tabyl(coinfection, other_indicator)

MB_data_EC_cf$coinfection<-as.character(MB_data_EC_cf$coinfection)
MB_data_EC_cf<-MB_data_EC_cf%>%mutate(coinfection=ifelse(
  coinfection=="1" & other_indicator=="1", "co", coinfection
))
MB_data_EC_final<-MB_data_EC_cf[!(MB_data_EC_cf$coinfection=="co"),]

MB_data_EC_final%>%tabyl(covid_indicator, other_indicator)

MB_data_EC_final%>%tabyl(comorbidity_summary_cat)
MB_data_EC_final$comorbidity_summary_cat_2<-MB_data_EC_final$comorbidity_summary_cat
MB_data_EC_final$comorbidity_summary_cat_2<-as.character(MB_data_EC_final$comorbidity_summary_cat_2)
MB_data_EC_final<-MB_data_EC_final%>% mutate(comorbidity_summary_cat_2=ifelse(
  comorbidity_summary_cat_2=="11-15", "6-10", comorbidity_summary_cat_2))
MB_data_EC_final$comorbidity_summary_cat_2<-as.factor(MB_data_EC_final$comorbidity_summary_cat_2)
MB_data_EC_final%>%tabyl(comorbidity_summary_cat_2)


MB_data_EC_final$WHO_severity

MB_data_EC_final_severity<-MB_data_EC_final[!(MB_data_EC_final$WHO_severity=="COVID_negative_infection"),]
MB_data_EC_final_severity$WHO_severity
MB_data_EC_final_severity%>%tabyl(WHO_severity)
MB_data_EC_final_severity$WHO_severity<-factor(MB_data_EC_final_severity$WHO_severity, 
                                               levels=c("WHO_Mild", "WHO_Moderate",
                                                        "WHO_Severe", "WHO_Fatal"))

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#creating BCG dataset
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCG_data<-as.data.frame(colData(tse_genus_BCG_filtered))

BCG_data_EC<-merge(BCG_data,mem_bcg, by = 'row.names', all = TRUE )

BCG_data_EC$Moraxella_log<-log10(BCG_data_EC$g__Moraxella+1)
BCG_data_EC$Streptococcus_log_10<-log10(BCG_data_EC$g__Streptococcus+1)
BCG_data_EC$Corynebacterium_log<-log10(BCG_data_EC$g__Corynebacterium+1)
BCG_data_EC$Dolosigranulum_log<-log10(BCG_data_EC$g__Dolosigranulum+1)
BCG_data_EC$Staphylococcus_log<-log10(BCG_data_EC$g__Staphylococcus+1)
BCG_data_EC$Actinomyces_log<-log10(BCG_data_EC$g__Actinomyces+1)
BCG_data_EC$Streptococcus_log<-log10(BCG_data_EC$g__Streptococcus+1)
BCG_data_EC$Haemophilus_log<-log10(BCG_data_EC$g__Haemophilus+1)
BCG_data_EC$Neisseria_log<-log10(BCG_data_EC$g__Neisseria+1)
BCG_data_EC$Pseudomonas_log<-log10(BCG_data_EC$g__Pseudomonas+1)
BCG_data_EC$oral_cavity_bacteria_log<-log10(BCG_data_EC$oral_cavity_bacteria+1)
BCG_data_EC$other_minority_species<-log10(BCG_data_EC$other_minority_species+1)

BCG_data_EC$WHO_severity
BCG_data_EC%>%tabyl(severity)
BCG_data_EC<-BCG_data_EC%>%mutate(covid_severity=severity)
BCG_data_EC%>%tabyl(covid_severity)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCG_data_EC$covid_severity<-as.character(BCG_data_EC$covid_severity)
BCG_data_EC_severity<-BCG_data_EC[!(BCG_data_EC$covid_severity=="COVID_negative_infection"),]
BCG_data_EC_severity%>%tabyl(covid_severity)

BCG_data_EC_severity$covid_severity<-factor(BCG_data_EC_severity$covid_severity, 
                                            levels = c("no infection", "Very Mild", "Mild"))

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCG_data_EC_severity_final<-BCG_data_EC_severity[!(BCG_data_EC_severity$covid_severity=="no infection"),]
BCG_data_EC_severity_final%>%tabyl(covid_severity)

BCG_data_EC_severity_final$covid_severity<-factor(BCG_data_EC_severity_final$covid_severity, 
                                                  levels = c("Very Mild", "Mild"))


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCG_data_EC<-BCG_data_EC%>%mutate(coinfection=covid_indicator)
BCG_data_EC%>%tabyl(coinfection)
BCG_data_EC%>%tabyl(coinfection, other_indicator)

BCG_data_EC$coinfection<-as.character(BCG_data_EC$coinfection)
BCG_data_EC<-BCG_data_EC%>%mutate(coinfection=ifelse(
  coinfection=="1" & other_indicator=="1", "co", coinfection
))

BCG_data_EC_final<-BCG_data_EC[!(BCG_data_EC$coinfection=="co"),]

BCG_data_EC_final%>%tabyl(covid_indicator, other_indicator)
BCG_data_EC_final%>%tabyl(covid_infection)


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#create variable for dividing negative infections as well
BCG_data_EC_final<-BCG_data_EC_final%>%mutate(severity_modified=severity)
BCG_data_EC_final%>%tabyl(severity_modified, overall_severity)
BCG_data_EC_final$severity_modified<-as.character(BCG_data_EC_final$severity_modified)

BCG_data_EC_final<-BCG_data_EC_final%>%mutate(severity_modified=ifelse(
  severity_modified=="no infection" & overall_severity=="Asymp", "None_Asymp", severity_modified
))
BCG_data_EC_final<-BCG_data_EC_final%>%mutate(severity_modified=ifelse(
  severity_modified=="no infection" & overall_severity=="Very Mild", "None_very_mild", severity_modified
))

BCG_data_EC_final<-BCG_data_EC_final%>%mutate(severity_modified=ifelse(
  severity_modified=="no infection" & overall_severity=="Mild", "None_mild", severity_modified
))

BCG_data_EC_final<-BCG_data_EC_final%>%mutate(severity_modified=ifelse(
  severity_modified=="no infection" & overall_severity=="UNK", "None_UNK", severity_modified
))

BCG_data_EC_final%>%tabyl(severity_modified)

BCG_data_EC_final_sev_mod<-BCG_data_EC_final[!(BCG_data_EC_final$severity_modified=="None_mild"),]
BCG_data_EC_final_sev_mod<-BCG_data_EC_final_sev_mod[!(BCG_data_EC_final_sev_mod$severity_modified=="None_UNK"),]
BCG_data_EC_final_sev_mod<-BCG_data_EC_final_sev_mod[!(BCG_data_EC_final_sev_mod$severity_modified=="COVID_negative_infection"),]
BCG_data_EC_final_sev_mod%>%tabyl(severity_modified)

BCG_data_EC_final_sev_mod$severity_modified<-factor(BCG_data_EC_final_sev_mod$severity_modified,
                                                    levels = c("None_Asymp", "None_very_mild", "Very Mild", "Mild"))

BCG_data_EC_final<-BCG_data_EC_final%>%mutate(severity_modified_2=severity_modified)
BCG_data_EC_final%>%tabyl(severity_modified_2, severity_modified)
BCG_data_EC_final$severity_modified_2<-as.character(BCG_data_EC_final$severity_modified_2)

BCG_data_EC_final<-BCG_data_EC_final%>%mutate(severity_modified_2=ifelse(
  severity_modified_2=="None_Asymp" , "None_Asymp_verymild", severity_modified_2
))
BCG_data_EC_final<-BCG_data_EC_final%>%mutate(severity_modified_2=ifelse(
  severity_modified_2=="None_very_mild" , "None_Asymp_verymild", severity_modified_2
))
BCG_data_EC_final%>%tabyl(severity_modified_2)


BCG_data_EC_final_sev_mod2<-BCG_data_EC_final[!(BCG_data_EC_final$severity_modified=="None_mild"),]
BCG_data_EC_final_sev_mod2<-BCG_data_EC_final_sev_mod2[!(BCG_data_EC_final_sev_mod2$severity_modified=="None_UNK"),]
BCG_data_EC_final_sev_mod2<-BCG_data_EC_final_sev_mod2[!(BCG_data_EC_final_sev_mod2$severity_modified=="COVID_negative_infection"),]
BCG_data_EC_final_sev_mod2%>%tabyl(severity_modified_2)

BCG_data_EC_final_sev_mod2$severity_modified_2<-factor(BCG_data_EC_final_sev_mod2$severity_modified_2,
                                                       levels = c("None_Asymp_verymild","Very Mild", "Mild"))
BCG_data_EC_final_sev_mod2%>%tabyl(severity_modified_2)
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BCG_data_EC_final_covidvsnone<-BCG_data_EC_final[!(BCG_data_EC_final$covid_infection=="other"),]
BCG_data_EC_final_covidvsnone%>%tabyl(covid_infection)
BCG_data_EC_final_covidvsnone$covid_infection<-factor(BCG_data_EC_final_covidvsnone$covid_infection,
                                                      levels = c("No", "covid"))

BCG_data_EC_final_covidvsother<-BCG_data_EC_final[!(BCG_data_EC_final$covid_infection=="No"),]
BCG_data_EC_final_covidvsother%>%tabyl(covid_infection)
BCG_data_EC_final_covidvsother$covid_infection<-factor(BCG_data_EC_final_covidvsother$covid_infection,
                                                       levels = c( "other", "covid"))

BCG_data_EC_final_othervsnone<-BCG_data_EC_final[!(BCG_data_EC_final$covid_infection=="covid"),]
BCG_data_EC_final_othervsnone%>%tabyl(covid_infection)
BCG_data_EC_final_othervsnone$covid_infection<-factor(BCG_data_EC_final_othervsnone$covid_infection,
                                                      levels = c("No", "other"))


BCG_data_EC_final_covidvsnone
BCG_data_EC_final_covidvsother
BCG_data_EC_final_othervsnone


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#####Creating distribution of concentration ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

co<-ggplot(BCG_data_EC, aes(x=Corynebacterium_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()
st<-ggplot(BCG_data_EC, aes(x=Streptococcus_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()
mor<-ggplot(BCG_data_EC, aes(x=Moraxella_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()
dol<-ggplot(BCG_data_EC, aes(x=Dolosigranulum_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()
sta<-ggplot(BCG_data_EC, aes(x=Staphylococcus_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()
act<-ggplot(BCG_data_EC, aes(x=Actinomyces_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()
hae<-ggplot(BCG_data_EC, aes(x=Haemophilus_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()
nei<-ggplot(BCG_data_EC, aes(x=Neisseria_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()
ocb<-ggplot(BCG_data_EC, aes(x=oral_cavity_bacteria_log)) + geom_histogram(bins=100,color="black", fill="gold")+
  theme_bw()

grid.arrange(co, st, mor, dol, sta, act, hae, nei, ocb, ncol=2)



co<-ggplot(MB_data_EC_final, aes(x=Corynebacterium_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()
st<-ggplot(MB_data_EC_final, aes(x=Streptococcus_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()
mor<-ggplot(MB_data_EC_final, aes(x=Moraxella_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()
dol<-ggplot(MB_data_EC_final, aes(x=Dolosigranulum_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()
sta<-ggplot(MB_data_EC_final, aes(x=Staphylococcus_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()
act<-ggplot(MB_data_EC_final, aes(x=Actinomyces_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()
hae<-ggplot(MB_data_EC_final, aes(x=Haemophilus_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()
nei<-ggplot(MB_data_EC_final, aes(x=Neisseria_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()
ocb<-ggplot(MB_data_EC_final, aes(x=oral_cavity_bacteria_log)) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()

grid.arrange(co, st, mor, dol, sta, act, hae, nei, ocb, ncol=2)

#variability 
#log10(MB_data_EC$g__Actinomyces+1)
ggplot(MB_data_EC, aes(x=log10(as.numeric(BQ_copy_per_ul)+1))) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()+labs(title = "UMCU patients")

summary(as.numeric(MB_data_EC$BQ_copy_per_ul))

#variability 
ggplot(BCG_data_EC, aes(x=log10(as.numeric(BQ_copy_per_ul)+1))) + geom_histogram(bins=100,color="black", fill="steelblue")+
  theme_bw()+labs(title = "BCG-Corona")

summary(as.numeric(BCG_data_EC$BQ_copy_per_ul))


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Crude MB-COVID####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Actinomyces
Actinomyces<-lm(Actinomyces_log~covid_indicator, data=MB_data_EC_final)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~WHO_severity, data=MB_data_EC_final)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)


Actinomyces<-lm(Actinomyces_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)
MB_data_EC_final_severity

#Corynebacterium
Corynebacterium<-lm(Corynebacterium_log~covid_indicator, data=MB_data_EC_final)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)


Corynebacterium<-lm(Corynebacterium_log~WHO_severity, data=MB_data_EC_final)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)


#Dolosigranulum
Dolosigranulum<-lm(Dolosigranulum_log~covid_indicator, data=MB_data_EC_final)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~WHO_severity, data=MB_data_EC_final)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

#Haemophilus
Haemophilus<-lm(Haemophilus_log~covid_indicator, data=MB_data_EC_final)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~WHO_severity, data=MB_data_EC_final)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

#Moraxella
Moraxella<-lm(Moraxella_log~covid_indicator, data=MB_data_EC_final)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~WHO_severity, data=MB_data_EC_final)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

#Neisseria
Neisseria<-lm(Neisseria_log~covid_indicator, data=MB_data_EC_final)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~WHO_severity, data=MB_data_EC_final)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

#Pseudomonas
Pseudomonas<-lm(Pseudomonas_log~covid_indicator, data=MB_data_EC_final)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~WHO_severity, data=MB_data_EC_final)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)


#Staphylococcus
Staphylococcus<-lm(Staphylococcus_log~covid_indicator, data=MB_data_EC_final)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)


Staphylococcus<-lm(Staphylococcus_log~WHO_severity, data=MB_data_EC_final)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

#Streptococcus
Streptococcus<-lm(Streptococcus_log~covid_indicator, data=MB_data_EC_final)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~WHO_severity, data=MB_data_EC_final)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)


#Oral_cavity
Oral_cavity<-lm(oral_cavity_bacteria_log~covid_indicator, data=MB_data_EC_final)
summary(Oral_cavity)
confint(Oral_cavity, level=0.95)
tbl_regression(Oral_cavity, exponentiate = FALSE, intercept = TRUE)

Oral_cavity<-lm(oral_cavity_bacteria_log~WHO_severity, data=MB_data_EC_final)
summary(Oral_cavity)
confint(Oral_cavity, level=0.95)
tbl_regression(Oral_cavity, exponentiate = FALSE, intercept = TRUE)

Oral_cavity<-lm(oral_cavity_bacteria_log~WHO_severity, data=MB_data_EC_final_severity)
summary(Oral_cavity)
confint(Oral_cavity, level=0.95)
tbl_regression(Oral_cavity, exponentiate = FALSE, intercept = TRUE)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#####Crude BCG-Corona####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCG_data_EC_final_covidvsnone
BCG_data_EC_final_covidvsother
BCG_data_EC_final_othervsnone

#Actinomyces
Actinomyces<-lm(Actinomyces_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~covid_severity, data=BCG_data_EC_severity)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

BCG_data_EC_final_sev_mod

Actinomyces<-lm(Actinomyces_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

#Corynebacterium

Corynebacterium<-lm(Corynebacterium_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~covid_severity, data=BCG_data_EC_severity)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

#Dolosigranulum
Dolosigranulum<-lm(Dolosigranulum_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~covid_severity, data=BCG_data_EC_severity)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

#Haemophilus

Haemophilus<-lm(Haemophilus_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~covid_severity, data=BCG_data_EC_severity)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

#Moraxella
Moraxella<-lm(Moraxella_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~covid_severity, data=BCG_data_EC_severity)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

#Neisseria
Neisseria<-lm(Neisseria_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~covid_severity, data=BCG_data_EC_severity)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

#Pseudomonas
Pseudomonas<-lm(Pseudomonas_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~covid_severity, data=BCG_data_EC_severity)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

#Staphylococcus
Staphylococcus<-lm(Staphylococcus_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~covid_severity, data=BCG_data_EC_severity)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)


Staphylococcus<-lm(Staphylococcus_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

#Streptococcus

Streptococcus<-lm(Streptococcus_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~covid_severity, data=BCG_data_EC_severity)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~covid_severity, data=BCG_data_EC_severity_final)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

#oral_bacteria

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_infection, data=BCG_data_EC_final_covidvsnone)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_infection, data=BCG_data_EC_final_covidvsother)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_infection, data=BCG_data_EC_final_othervsnone)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_severity, data=BCG_data_EC_severity)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_severity, data=BCG_data_EC_severity_final)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~severity_modified_2, data=BCG_data_EC_final_sev_mod2)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~severity_modified, data=BCG_data_EC_final_sev_mod)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

BCG_data_EC_final_sev_mod%>%tabyl(covid_infection, overall_severity)
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Adjusted MB-COVID####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Actinomyces
Actinomyces<-lm(Actinomyces_log~covid_indicator+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~WHO_severity+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~WHO_severity+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final_severity)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)


#Corynebacterium
Corynebacterium<-lm(Corynebacterium_log~covid_indicator+
                      sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)


Corynebacterium<-lm(Corynebacterium_log~WHO_severity+
                      sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)


Corynebacterium<-lm(Corynebacterium_log~WHO_severity+
                      sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final_severity)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

MB_data_EC_final_severity

#Dolosigranulum
Dolosigranulum<-lm(Dolosigranulum_log~covid_indicator+
                     sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~WHO_severity+
                     sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~WHO_severity+
                     sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final_severity)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

#Haemophilus
Haemophilus<-lm(Haemophilus_log~covid_indicator+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~WHO_severity+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~WHO_severity+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final_severity)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

#Moraxella
Moraxella<-lm(Moraxella_log~covid_indicator+
                sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~WHO_severity+
                sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)


Moraxella<-lm(Moraxella_log~WHO_severity+
                sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final_severity)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

#Neisseria
Neisseria<-lm(Neisseria_log~covid_indicator+
                sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~WHO_severity+
                sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)
MB_data_EC_final_severity

Neisseria<-lm(Neisseria_log~WHO_severity+
                sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final_severity)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

#Pseudomonas
Pseudomonas<-lm(Pseudomonas_log~covid_indicator+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~WHO_severity+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~WHO_severity+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, 
                data=MB_data_EC_final_severity)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

#Staphylococcus
Staphylococcus<-lm(Staphylococcus_log~covid_indicator+
                     sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)


Staphylococcus<-lm(Staphylococcus_log~WHO_severity+
                     sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~WHO_severity+
                     sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users,
                   data=MB_data_EC_final_severity)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

#Streptococcus
Streptococcus<-lm(Streptococcus_log~covid_indicator+
                    sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~WHO_severity+
                    sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~WHO_severity+
                    sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final_severity)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)


#Oral_cavity
Oral_cavity<-lm(oral_cavity_bacteria_log~covid_indicator+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Oral_cavity)
confint(Oral_cavity, level=0.95)
tbl_regression(Oral_cavity, exponentiate = FALSE, intercept = TRUE)

Oral_cavity<-lm(oral_cavity_bacteria_log~WHO_severity+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, data=MB_data_EC_final)
summary(Oral_cavity)
confint(Oral_cavity, level=0.95)
tbl_regression(Oral_cavity, exponentiate = FALSE, intercept = TRUE)


Oral_cavity<-lm(oral_cavity_bacteria_log~WHO_severity+
                  sex+age_yrs+comorbidity_summary_cat_2+antibiotic_users, 
                data=MB_data_EC_final_severity)
summary(Oral_cavity)
confint(Oral_cavity, level=0.95)
tbl_regression(Oral_cavity, exponentiate = FALSE, intercept = TRUE)
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#####Adjusted BCG-Corona####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Actinomyces
Actinomyces<-lm(Actinomyces_log~covid_infection +
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~covid_infection+
                  sex+ age_yrs+antibiotic_users, data=BCG_data_EC_final_covidvsother)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~covid_infection+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

Actinomyces<-lm(Actinomyces_log~covid_severity+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)



Actinomyces<-lm(Actinomyces_log~covid_severity+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)


Actinomyces<-lm(Actinomyces_log~severity_modified_2+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

BCG_data_EC_final_sev_mod

Actinomyces<-lm(Actinomyces_log~severity_modified+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Actinomyces)
confint(Actinomyces, level=0.95)
tbl_regression(Actinomyces, exponentiate = FALSE, intercept = TRUE)

#Corynebacterium

Corynebacterium<-lm(Corynebacterium_log~covid_infection+
                      sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~covid_infection+
                      sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~covid_infection +
                      sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~covid_severity +
                      sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)



Corynebacterium<-lm(Corynebacterium_log~covid_severity +
                      sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~severity_modified_2+
                      sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

Corynebacterium<-lm(Corynebacterium_log~severity_modified+
                      sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Corynebacterium)
confint(Corynebacterium, level=0.95)
tbl_regression(Corynebacterium, exponentiate = FALSE, intercept = TRUE)

#Dolosigranulum
Dolosigranulum<-lm(Dolosigranulum_log~covid_infection +
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~covid_infection +
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~covid_infection +
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~covid_severity +
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)



Dolosigranulum<-lm(Dolosigranulum_log~covid_severity +
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~severity_modified_2+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

Dolosigranulum<-lm(Dolosigranulum_log~severity_modified+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Dolosigranulum)
confint(Dolosigranulum, level=0.95)
tbl_regression(Dolosigranulum, exponentiate = FALSE, intercept = TRUE)

#Haemophilus

Haemophilus<-lm(Haemophilus_log~covid_infection+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~covid_infection+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~covid_infection+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~covid_severity+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)


Haemophilus<-lm(Haemophilus_log~covid_severity+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

Haemophilus<-lm(Haemophilus_log~severity_modified_2+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)


Haemophilus<-lm(Haemophilus_log~severity_modified+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Haemophilus)
confint(Haemophilus, level=0.95)
tbl_regression(Haemophilus, exponentiate = FALSE, intercept = TRUE)

#Moraxella
Moraxella<-lm(Moraxella_log~covid_infection+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~covid_infection+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~covid_infection+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~covid_severity+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)


Moraxella<-lm(Moraxella_log~covid_severity+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~severity_modified_2+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

Moraxella<-lm(Moraxella_log~severity_modified+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Moraxella)
confint(Moraxella, level=0.95)
tbl_regression(Moraxella, exponentiate = FALSE, intercept = TRUE)

#Neisseria
Neisseria<-lm(Neisseria_log~covid_infection+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~covid_infection+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~covid_infection+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~covid_severity+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)


Neisseria<-lm(Neisseria_log~covid_severity+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~severity_modified_2+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

Neisseria<-lm(Neisseria_log~severity_modified+
                sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Neisseria)
confint(Neisseria, level=0.95)
tbl_regression(Neisseria, exponentiate = FALSE, intercept = TRUE)

#Pseudomonas
Pseudomonas<-lm(Pseudomonas_log~covid_infection+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~covid_infection+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~covid_infection+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~covid_severity+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)



Pseudomonas<-lm(Pseudomonas_log~covid_severity+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~severity_modified_2+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

Pseudomonas<-lm(Pseudomonas_log~severity_modified+
                  sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Pseudomonas)
confint(Pseudomonas, level=0.95)
tbl_regression(Pseudomonas, exponentiate = FALSE, intercept = TRUE)

#Staphylococcus
Staphylococcus<-lm(Staphylococcus_log~covid_infection+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~covid_infection+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~covid_infection+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~covid_severity+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~covid_severity+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~severity_modified_2+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

Staphylococcus<-lm(Staphylococcus_log~severity_modified+
                     sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Staphylococcus)
confint(Staphylococcus, level=0.95)
tbl_regression(Staphylococcus, exponentiate = FALSE, intercept = TRUE)

#Streptococcus

Streptococcus<-lm(Streptococcus_log~covid_infection+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~covid_infection+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~covid_infection+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~covid_severity+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~covid_severity+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~severity_modified_2+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

Streptococcus<-lm(Streptococcus_log~severity_modified+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(Streptococcus)
confint(Streptococcus, level=0.95)
tbl_regression(Streptococcus, exponentiate = FALSE, intercept = TRUE)

#oral_bacteria

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_infection+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsnone)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_infection+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_covidvsother)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_infection+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_othervsnone)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_severity+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~covid_severity+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_severity_final)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~severity_modified_2+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod2)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

oral_bacteria<-lm(oral_cavity_bacteria_log~severity_modified+
                    sex+ age_yrs+antibiotic_users , data=BCG_data_EC_final_sev_mod)
summary(oral_bacteria)
confint(oral_bacteria, level=0.95)
tbl_regression(oral_bacteria, exponentiate = FALSE, intercept = TRUE)

