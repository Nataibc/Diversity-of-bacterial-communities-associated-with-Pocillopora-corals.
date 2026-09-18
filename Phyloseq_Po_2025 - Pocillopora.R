#-------------------------------------------------------------------------------
#ANALISIS DEL MICROBIOMA BASADO EN SECUENCIACION AMPLICON DEL GEN 16S ARNr
#ANALISIS BASADO EN PHYLOSEQ Y OTROS
#Instalar y cargar los paquetes necesarios

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.20")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "seqLogo"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(BiocManager)
BiocManager::install("microbiome")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

remove.packages("globals")
install.packages("globals", dependencies = TRUE)

library(dada2); packageVersion("dada2")
library(BiocManager)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(tidyr)
library(stringr)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(vegan)
library(lessR)
library(ape)
library(data.table)
library(breakaway)
library(microbiome)
library(ggpubr)
library(kableExtra)
library (DESeq2)
library(ComplexHeatmap)
library(ANCOMBC)
library(DECIPHER)
library(phangorn)
#-------------------------------------------------------------------------------
#PHYLOSEQ pipeline (modificado por Jordan Ruiz, 2024)

#Hacer los anlisis usando el paquete de Phyloseq (McMurdie y Holmes, 2013)
#El obejto de phyloseq tiene 4 componentes: la tabla de otus, los datos de las muestras, la tabla de taxa, y el arbol filogenetico

#Datos de las muestras (Metadata)
#Importar los metadatos y fusionar las columnas codigo y muestra
Metadata <- read.csv("/Users/Usuario/Desktop/Analisis_POC_16S/Metadatos.csv", sep = ",", header = TRUE, check.names = FALSE)
Metadata$Sample <- paste(Metadata$Code, Metadata$Sample, sep = "-")
colnames(Metadata)
head(Metadata)
Metadata2 <- select(Metadata, -Code)
head(Metadata2)
Locality <- as.factor(Metadata2$Locality)
ORF_Type <- as.factor(Metadata2$ORF_Type)
Year_Month <- as.factor(Metadata2$Year_Month)
Upwelling <- as.factor(Metadata2$Upwelling)

#Tabla de otus y tabla de taxa
#Crear las matrices que van a hacer parte del objeto de phyloseq
asv_mat<- read.csv("/Users/Usuario/Desktop/Analisis_POC_16S/Secuencias/Pocillopora_sample.asvs_taxa_2024.csv", sep = ",", header = TRUE, check.names = FALSE)
tax_mat<- read.csv("/Users/Usuario/Desktop/Analisis_POC_16S/Secuencias/Pocillopora_asvs_taxa_2024.csv", sep = ",", header = TRUE, check.names = FALSE)

head(asv_mat)
head(tax_mat)

asv_mat <- asv_mat %>% 
  tibble::column_to_rownames("asv")

tax_mat <- tax_mat %>%
  tibble::column_to_rownames("asv")

Metadata2 <- Metadata2 %>%
  tibble::column_to_rownames("Sample")

asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
Metadata = sample_data(Metadata2)

#OBJETO DE PHYLOSEQ Y EXPLORACION DE LOS DATOS
#Crear el objeto de phyloseq
Po <- phyloseq(OTU, TAX, Metadata)
Po

sample_names(Po)
rank_names(Po)
sample_variables(Po)

#Revisar el numero de lecturas por muestras, si hay alguna muestra con muy pocas secuencias, se debe remover
#Todas las muestras tienen entre 17000 y 17400 secuencias por lo que no se debe remover ninguna
readcount = data.table(as(sample_data(Po), "data.frame"),
                       TotalReads = sample_sums(Po), 
                       keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")
ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
head(readcount[order(readcount$TotalReads), c("SampleID", "TotalReads")])

write.csv(readcount, 
          file="/Users/Usuario/Desktop/Analisis_POC_16S/read_count_table.csv")

#Calcular el tamano de las librerias por muestra
tab <- otu_table(asv_mat, taxa_are_rows = TRUE)
class(tab) <- "matrix"
tab <- t(tab) # transpose observations to rows
sum_seq <- rowSums(tab)
plot(sum_seq, ylim=c(0,65000), main=c("Number of counts per sample"), xlab=c("Samples"), ylab=c("Number of reads"))
sum_seq
mean(sum_seq)
sd(sum_seq)
min(sum_seq)
max(sum_seq)

#Hacer la curva de rarefraccion para determinar la profundidad de la secuenciacon
#La profundidad nos dice si se secuencio lo suficiente para representar la diversidad de la comunidad


source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")

rare1 <- ggrare(Po, step = 100, color = "Locality", plot = TRUE, parallel = TRUE, se = FALSE)
figure_rare <- rare1 + theme_classic() +
  labs(x = "\nSample size", y = "Number of ASV observed\n") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "gray")) +
  scale_x_continuous(labels=comma) +
  scale_color_manual(values=c("#404788FF",  "#238A8D",  "#FDE725FF" )) +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"))
figure_rare


ggsave("/Users/Usuario/Desktop/Analisis_POC_16S/Plot_rarefraction_curve_2025.png", 
       figure_rare, width = 8, height = 6, dpi = 600)


#Filtrar la taxonomia, remover filos que tienen solo 1 (que solo estan en el dominio bacteria) y el NA (filos no asignados) los cuales son artifactos
table(tax_table(Po) [, "Phylum"], exclude = NULL)
Po <- subset_taxa(Po, !is.na(Phylum) & !Phylum %in% c("", "Abditibacteriota", "Armatimonadota", "Dadabacteria", 
                                                      "Fibrobacterota", "Hydrogenedentes", "Nitrospirota", "AncK6", "Deinococcota", "Latescibacterota", "Nitrospinota", "<NA>"))
Po

plot_bar(Po,  fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

ggsave("/Users/Usuario/Desktop/Analisis_POC_16S/bar_plot_phylum_filtered_2024.png", 
       plot = (plot_bar(Po, fill = "Phylum")), width = 8, height = 6, dpi = 600)

#Hacer un analisis de la prevalencia de los filos en las muestras
prevdf <- apply(X = otu_table(Po),
                MARGIN = ifelse(taxa_are_rows(Po), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(Po),
                     tax_table(Po))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(Po, "Phylum"))
fig_prev <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(Po),color=Phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

fig_prev

ggsave("/Users/Usuario/Desktop/Analisis_POC_16S/plot_phylum_prevalence_2025.png", 
       fig_prev, width = 8, height = 6, dpi = 600)

#Remover las ASVs que estan por debajo de una prevalencia del 5% (todos lo que este en menos de 3 muestras lo va a quitar)
#Forma de quitar artefactos y contaminantes 
prevalenceThreshold <- 0.05 * nsamples(Po)
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
Po1 <- prune_taxa(keepTaxa, Po)
saveRDS(Po1, "Po1.rds")

Po1
class(Po)
str(Po)

tax_table_df <- as.data.frame(tax_table(Po1))
n_phylum  <- n_distinct(tax_table_df$Phylum, na.rm = TRUE)
n_class   <- n_distinct(tax_table_df$Class, na.rm = TRUE)
n_order   <- n_distinct(tax_table_df$Order, na.rm = TRUE)
n_family  <- n_distinct(tax_table_df$Family, na.rm = TRUE)
n_genus   <- n_distinct(tax_table_df$Genus, na.rm = TRUE)
install.packages("writexl")  # solo si no lo tienes
library(writexl)

write_xlsx(tax_table_df, path = "C:/Users/Usuario/Desktop/Analisis_POC_16S/tax_table_Po1.xlsx")

write.csv(cbind(data.frame(otu_table(Po1)),
                tax_table(Po1)), 
          file="/Users/Usuario/Desktop/Analisis_POC_16S/filtered_otu_table.csv")


#Hacer la rarefraccion (normalizar) de los datos para los analisis de diversidad
set.seed(111); .Random.seed
Po_rarefied <- rarefy_even_depth(Po1, 
                               sample.size = min(colSums(otu_table(Po1))), 
                               rngseed = 63)
Po_rarefied

head(phyloseq::sample_sums(Po_rarefied))

write.csv(cbind(data.frame(otu_table(Po_rarefied)),
                tax_table(Po_rarefied)), 
          file="/Users/Usuario/Desktop/Analisis_POC_16S/filtered_rarefied_otu_table.csv")

saveRDS(Po_rarefied, file = "Po_rarefied.rds")
saveRDS(Po_rarefied, file = "/Users/Usuario/Desktop/Analisis_POC_16S/Po_rarefied.rds")

#Subset de datos por tipo y localidad(nat)-----

#tipo1

Po_rarefiedtipo1 <- subset_samples(Po_rarefied, ORF_Type == "Type_1")
Po_rarefiedtipo3 <- subset_samples(Po_rarefied, ORF_Type == "Type_3")

#Localidades
Po_rarefie_Cupica <- subset_samples(Po_rarefied, Locality == "Golfo Cupica")
Po_rarefie_Utria <- subset_samples(Po_rarefied, Locality == "Utria")
Po_rarefie_Gorgona <- subset_samples(Po_rarefied, Locality == "Gorgona")
#-------------------------------------------------------------------------------
#ANALISIS DE LA DIVERSIDAD DE LA COMUNIDAD MICROBIANA parte 1
#DIVERSIDAD ALFA
#Hacer el analisis de diversidad alfa (diversidad dentro de la muestra)
#Usar el objeto filtrado y normalizado (rarefraccion) para el analisis
#Este es el analisis univariado de la diversidad

#Revisar la correlacion entre las lecturas totales y la riqueza observada
#La diversidad aumenta a medida que se obtienen mas lecturas como es de esperarse (en mi caso no)
#hacerlo sobre el objeto original de phyloseq (Po)
figure_obs_rich <- ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(Po),
                         "observed" = phyloseq::estimate_richness(Po, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = TRUE) +
  labs(x = "\nTotal Reads", y = "Number of ASV observed\n") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm")) +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches")) +
  scale_x_continuous(labels=comma)
figure_obs_rich

ggsave("/Users/Usuario/Desktop/Analisis_POC_16S/plot_observed_richness_2025.png", 
       figure_obs_rich, width = 8, height = 6, dpi = 600)

#Determinar los indices de diversidad, riqueza y dominancia para cada factor y exportar la tabla
plot_richness(Po_rarefied) 

adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(Po_rarefied, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(Po_rarefied, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(Po_rarefied, measures = "Simpson"),
  "Chao1" = phyloseq::estimate_richness(Po_rarefied, measures = "Chao1"),
  "InvSimpson" = phyloseq::estimate_richness(Po_rarefied, measures = "InvSimpson"),
  "Locality" = phyloseq::sample_data(Po_rarefied)$Locality, 
  "ORF_Type" = phyloseq::sample_data(Po_rarefied)$ORF_Type,
  "Upwelling" = phyloseq::sample_data(Po_rarefied)$Upwelling,
  "Year_Month" = phyloseq::sample_data(Po_rarefied)$Year_Month)

head(adiv)

write.csv(adiv, file="/Users/Usuario/Desktop/Analisis_POC_16S/alpha_div.csv")

#-------------teniendo en cuenta el indice de Pielou (nat)----------------

observed <- phyloseq::estimate_richness(Po_rarefied, measures = "Observed")$Observed
shannon <- phyloseq::estimate_richness(Po_rarefied, measures = "Shannon")$Shannon
simpson <- phyloseq::estimate_richness(Po_rarefied, measures = "Simpson")$Simpson
chao1 <- phyloseq::estimate_richness(Po_rarefied, measures = "Chao1")$Chao1
invsimpson <- phyloseq::estimate_richness(Po_rarefied, measures = "InvSimpson")$InvSimpson

adivpielou <- data.frame(
  "Observed" = observed,
  "Shannon" = shannon,
  "Simpson" = simpson,
  "Chao1" = chao1,
  "InvSimpson" = invsimpson,
  "Pielou" = shannon / log(observed),
  "Locality" = phyloseq::sample_data(Po_rarefied)$Locality, 
  "ORF_Type" = phyloseq::sample_data(Po_rarefied)$ORF_Type,
  "Upwelling" = phyloseq::sample_data(Po_rarefied)$Upwelling,
  "Year_Month" = phyloseq::sample_data(Po_rarefied)$Year_Month
)
head(adivpielou)

write.csv(adivpielou, file="/Users/Usuario/Desktop/Analisis_POC_16S/alpha_div_Pielou.csv")

#-----no va ya porque no se cumplieron los supuestos----------------
#Hacer analisis estadistico (univariado) no parametrico para la diversidad, obvervada y shannon con base en cada factor
#Homegeneidad de varianzas (si la varianzas son homogeneas, probbar la normalidad para decidir entre ANOVA o Kruskal-Wallis) con test de Levene y diferencias con test de Wilcox
library(car)
leveneTest(y = adivpielou$Observed, group = adiv$ORF_Type, center = "median")
by(adivpielou$Observed, adiv$ORF_Type, shapiro.test)

#Locality
leveneTest(y = adiv$Shannon, group = adiv$Locality, center = "median")
by(adiv$Shannon, adiv$Locality, shapiro.test)
shapiro.test(adiv$Shannon, adiv$Locality) 
hist(adiv$Shannon)
kruskal.test(Shannon ~ Locality, data = adiv)
kruskal.test(Observed ~ Locality, data = adiv)
TukeyHSD_Shannon <- TukeyHSD(aov(Shannon ~ Locality, data =  adiv))
TukeyHSD_Shannon_df <- data.frame(TukeyHSD_Shannon$Locality)
TukeyHSD_Shannon_df$measure = "Shannon"
TukeyHSD_Shannon_df$shapiro_test_pval = (shapiro.test(residuals(aov(Shannon ~ Locality, data =  adiv))))$p.value
TukeyHSD_Shannon_df

TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ Locality, data =  adiv))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$Locality)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ Locality, data =  adiv))))$p.value
TukeyHSD_Observed_df

#ORF_ Type
leveneTest(y = adiv$Shannon, group = adiv$ORF_Type, center = "median")
by(adiv$Shannon, adiv$ORF_Type, shapiro.test)
shapiro.test(adiv$Shannon) 
hist(adiv$Shannon)
kruskal.test(Shannon ~ ORF_Type, data = adiv)
kruskal.test(Observed ~ ORF_Type, data = adiv)
TukeyHSD_Shannon <- TukeyHSD(aov(Shannon ~ ORF_Type, data =  adiv))
TukeyHSD_Shannon_df <- data.frame(TukeyHSD_Shannon$ORF_Type)
TukeyHSD_Shannon_df$measure = "Shannon"
TukeyHSD_Shannon_df$shapiro_test_pval = (shapiro.test(residuals(aov(Shannon ~ ORF_Type, data =  adiv))))$p.value
TukeyHSD_Shannon_df
TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ ORF_Type, data =  adiv))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$ORF_Type)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ ORF_Type, data =  adiv))))$p.value
TukeyHSD_Observed_df

#Hacer otros analisis y graficas con la significancia----
# Generamos un objeto `phyloseq` sin taxa que sume 0 reads
Po2 <- prune_taxa(taxa_sums(Po_rarefied) > 0, Po_rarefied)
Po2
# Calculamos los índices de diversidad
#(Asegurarse de cargar los paquetes global para que funcione)

tab <- microbiome::diversity(Po2, index = "all")

# Y finalmente visualizamos la tabla de resultados
head(tab) %>%
  kable(format = "html", col.names = colnames(tab), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "310px")

write.csv(tab, file="/Users/Usuario/Desktop/Analisis_POC_16S/alpha_div2.csv")


Po2_meta <- meta(Po2)

head(Po2_meta) %>%
  kable(format = "html", col.names = colnames(Po2_meta), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "500px")

Po2_meta$Shannon <- tab$diversity_shannon
Po2_meta$Simpson <- tab$diversity_gini_simpson

Shannon <- tab$shannon
Simpson <- tab$gini_simpson

metadf <- Po2_meta
metadf$Locality.fact <- as.factor(metadf$Locality)
metadf$ORF_Type.fact <- as.factor(metadf$ORF_Type)
metadf$Year_Month.fact <- as.factor(metadf$Year_Month)
metadf$Upwelling.fact <- as.factor(metadf$Upwelling)

# Obtenemos las variables desde nuestro objeto `phyloseq`
Locality <- levels(metadf$Locality.fact)
ORF_Type <- levels(metadf$ORF_Type.fact)
Year_Month <- levels(metadf$Year_Month.fact)
Upwelling <- levels(metadf$Upwelling.fact)

# Creamos una lista de lo que queremos comparar
pares.Locality <- combn(seq_along(Locality), 2, simplify = FALSE, FUN = function(i)Locality[i])
pares.ORF_Type <- combn(seq_along(ORF_Type), 2, simplify = FALSE, FUN = function(i)ORF_Type[i])
pares.Year_Month <- combn(seq_along(Year_Month), 2, simplify = FALSE, FUN = function(i)Year_Month[i])
pares.Upwelling <- combn(seq_along(Upwelling), 2, simplify = FALSE, FUN = function(i)Upwelling[i])

# Imprimimos en pantalla el resultado
print(pares.Locality)
print(pares.ORF_Type)
print(pares.Year_Month)
print(pares.Upwelling)


boxplot1 <- ggboxplot(Po2_meta, x = "Locality", y = "Shannon",
                      add = "jitter", fill = "Locality", palette = c( "#fdae6b", "#bcbddc", "#7fcdbb")) +
  labs(x = "\nLocality", y = "Shannon\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = ))) +
  ylim(3.0, 6.0)

boxplot2 <- ggboxplot(Po2_meta, x = "ORF_Type", y = "Shannon",
                      add = "jitter", fill = "ORF_Type", palette = c("#ffeda0","#9ecae1")) +
  labs(x = "\nORF_Type", y = "Shannon\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(3.2, 6.0)


b1 <- boxplot1 + stat_compare_means(comparisons = pares.Locality)
print(b1)

b2 <- boxplot2 + stat_compare_means(comparisons = pares.ORF_Type)
print(b2)


figure_diversities <- ggarrange(b5,b4,b1,b8,b7,b6,
                                labels=c("A","B","C","D","E","F"),
                                ncol=3,nrow=2)

figure_diversities_index <- figure_diversities + theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"))
figure_diversities_index

ggsave("/Users/Usuario/Desktop/Analisis_POC_16S/plot_diversities_2024.png", 
       figure_diversities_index, width = 18, height = 12, dpi = 600)
#-----------------------
#Modelos lineales generalizados con distribución de Poisson (variables discretas)

#install.packages("broom")
#install.packages("gt")
#install.packages("ggeffects")
#install.packages("emmeans") 
library(ggeffects)
library(broom)
library(dplyr)
library(gt)
library(emmeans)
library(car) #Anova
library(MASS) #binomial negativo

#Riqueza observada 

# Efecto ppal de la localidad
#glm_obs_inter_l<- glm(Observed ~ Locality, data = adiv, family = quasipoisson(link = "log"))
#summary(glm_obs_inter_l)
#exp(confint(glm_obs_inter_l))
#     anova II
#anova(glm_obs_inter_l, test = "Chisq")

#         Post-hoc
#emm_locality <- emmeans(glm_obs_inter_l, ~ Locality)
#pairs(emm_locality, adjust = "tukey")  # Comparaciones entre localidades

#    Boxplot localidad
#adivpielou$Locality <- factor(adivpielou$Locality, levels = c("Golfo Cupica", "Utria", "Gorgona"))
#boxplot_glm_obs_inter_l<- ggboxplot(adivpielou, x = "Locality", y = "Observed",
 #                                   add = "jitter", fill = "Locality", palette = c( "#bcbddc","#fdae6b","#7fcdbb")) +
  #labs(x = "\nLocality", y = "Observed\n") +
  #theme(axis.title = element_text(size = 16),
   #     axis.title.x = element_blank(),
#        axis.text = element_text(size = 14),
#       axis.ticks.length = unit(.25, "cm"),
#       axis.text.x = element_text(angle = 90)) +
# theme(legend.position = "none") +
# theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
#       panel.border = element_rect(colour = "black", fill=NA, size=1),
#       axis.line = element_blank()) +
# guides(color = guide_legend(override.aes = list(size = 1))) +
# ylim(10, 250)

#boxplot_glm_obs_inter_l

## Efecto ppal del linaje evolutivo
#glm_obs_inter_t<- glm(Observed ~ ORF_Type, data = adiv, family = poisson(link = "log"))
#summary(glm_obs_inter_t)
#anova(glm_obs_inter_t, test = "Chisq")
#Boxplot
#boxplot_glm_obs_inter_t<- ggboxplot(adivpielou, x = "ORF_Type", y = "Observed",
#                                   add = "jitter", fill = "ORF_Type", palette = c("#ffeda0","#9ecae1")) +
# labs(x = "\nORF_Type", y = "Observed\n") +
# theme(axis.title = element_text(size = 16),
#       axis.title.x = element_blank(),
#       axis.text = element_text(size = 14),
#       axis.ticks.length = unit(.25, "cm"),
#       axis.text.x = element_text(angle = 90)) +
# theme(legend.position = "none") +
# theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
#       panel.border = element_rect(colour = "black", fill=NA, size=1),
  #       axis.line = element_blank()) +
#guides(color = guide_legend(override.aes = list(size = 1))) +
# ylim(10, 250)

#boxplot_glm_obs_inter_t

#Interacción
#glm Poisson
glm_obs_inter<- glm(Observed ~ Locality * ORF_Type, data = adiv, family = poisson(link = "log"))
deviance(glm_obs_inter)/df.residual(glm_obs_inter)
summary(glm_obs_inter)
anova(glm_obs_inter, test = "Chisq")
exp(confint(glm_obs_inter))

#glm binomial negativo
glm_obs_inter_nb<- glm.nb(Observed ~ Locality * ORF_Type, data = adiv, link = "log")
summary(glm_obs_inter_nb)

AIC(glm_obs_inter, glm_obs_inter_nb)
BIC(glm_obs_inter, glm_obs_inter_nb)

Anova(glm_obs_inter_nb)

#Grafico de interacción
eff_obs <- ggemmeans(glm_obs_inter_nb, terms = c("Locality", "ORF_Type"))
eff_obs$x <- factor(eff_obs$x, levels = c("Golfo Cupica", "Utria", "Gorgona")) #cambia el orden de las localidades

Riqueza_obs_inter <- ggplot(eff_obs, aes(x = x, y = predicted, color = group, group = group)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + 
  labs(
    x = "Localidad",
    y = "Riqueza Observada (media ± IC 95%)",
    color = "ORF_Type",
    title = "Interacción entre Localidad y Tipo ORF sobre la Riqueza Observada"
  ) +
  scale_color_manual(
    values = c("Type_1" = "#ffeda0", 
               "Type_3" = "#9ecae1")
  ) +
  scale_y_continuous(
    limits = c(105, 200),     # Cambia estos valores según tus datos
    breaks = seq(105, 200, 10) # Controla los cortes del eje Y
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),        # Quita líneas verticales
    panel.grid.minor.x = element_blank(),        
    panel.grid.major.y = element_line(color = "grey98"), # Deja solo horizontales
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),   # Añade líneas de los ejes
    axis.ticks = element_line(color = "black")   # Añade ticks en los ejes
  )

Riqueza_obs_inter

#Modelo lineal generalizado con distribución de Gamma (variables continuas)
#Shannon Localidad

lm1 <- lm(Shannon ~ Locality*ORF_Type, data = adiv)
summary(lm1)

res1 <- residuals(lm1)
qqPlot(res1)

#install.packages("nortest")
library(nortest)
ad.test(res1)

adiv$inter <- interaction(adiv$Locality, adiv$ORF_Type)
leveneTest(res1 ~ Locality, data = adiv)
leveneTest(res1 ~ ORF_Type, data = adiv)
leveneTest(res1 ~ inter, data = adiv)

hist(res1)

#Shannon interacción 
glm_shannon_inter <- glm(Shannon ~ Locality * ORF_Type, data = adiv, family = Gamma(link = "log"))
summary(glm_shannon_inter)
Anova(glm_shannon_inter)

# Generar los efectos estimados del modelo
eff_shannon <- ggemmeans(glm_shannon_inter, terms = c("Locality", "ORF_Type"))
eff_shannon$x <- factor(eff_shannon$x, levels = c("Golfo Cupica", "Utria", "Gorgona"))

Shannon_inter <- ggplot(eff_shannon, aes(x = x, y = predicted, color = group, group = group)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  labs(
    x = "Localidad",
    y = "Índice de Shannon (media ± IC 95%)",
    color = "ORF_Type",
    title = "Efecto del linaje y localidad sobre la diversidad (Shannon)"
  ) +
  scale_color_manual(
    values = c("Type_1" = "#ffeda0", 
               "Type_3" = "#9ecae1")
  ) +
  scale_y_continuous(
    limits = c(4, 5.1),     # Cambia estos valores según tus datos
    breaks = seq(4, 5, 0.5) # Controla los cortes del eje Y
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),        # Quita líneas verticales
    panel.grid.minor.x = element_blank(),        
    panel.grid.major.y = element_line(color = "grey98"), # Deja solo horizontales
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),   # Añade líneas de los ejes
    axis.ticks = element_line(color = "black")   # Añade ticks en los ejes
  )
Shannon_inter

#Modelo lineal generalizado con Beta-dispersión(variables entre 0 - 1)
#install.packages("betareg")
library(betareg)

##Equitabilidad de Pielou

glm_pielou_inter <- betareg(Pielou ~ Locality * ORF_Type, data = adivpielou, link = "logit")
summary(glm_pielou_inter)
library(car)
Anova(glm_pielou_inter, type = 2)

# Generar los efectos estimados del modelo
eff_pielou<- ggemmeans(glm_pielou_inter, terms = c("Locality", "ORF_Type"))
eff_pielou$x <- factor(eff_pielou$x, levels = c("Golfo Cupica", "Utria", "Gorgona"))

Pielou_inter <- ggplot(eff_pielou, aes(x = x, y = predicted, color = group, group = group)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  #geom_text(aes(label = .group, y = predicted + 2), color = "black", size = 5) + 
  labs(
    x = "Localidad",
    y = "Equitatividad de Pielou (media ± SD)",
    color = "ORF_Type",
    title = "Efecto del linaje y localidad sobre la equitabilidad de (Pielou)"
  ) +
  scale_color_manual(
    values = c("Type_1" = "#ffeda0", 
               "Type_3" = "#9ecae1")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),        # Quita líneas verticales
    panel.grid.minor.x = element_blank(),        
    panel.grid.major.y = element_line(color = "grey98"), # Deja solo horizontales
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),   # Añade líneas de los ejes
    axis.ticks = element_line(color = "black")   # Añade ticks en los ejes
  )
Pielou_inter


#-------------------------------------------------------------------------------
#ANALISIS DE LA DIVERSIDAD DE LA COMUNIDAD MICROBIANA parte 2
#DIVERSIDAD BETA
#Hacer el analisis de diversidad beta (diversidad entre muestras)
#Usar el objeto filtrado y normalizado (rarefraccion) para el analisis
#Este es el analisis multivariado de la diversidad
library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)

#ANALISIS Y PCoA
#Determinar la disimilutud como medida de distancia
Po_rarefied %>% transform_sample_counts(function(x) x/sum(x)) %>%
  otu_table() %>%
 t() %>%
 sqrt() %>%
 as.data.frame() %>%
 vegdist(binary=F, method = "bray") -> dist

#Correr analisis de ordenacion PCoA con la distancia
#Bray Curtis
#Modifcado para ver tipos y localidaes (se quitó shannon)
ord <- ordinate(Po_rarefied,"PCoA",dist)
ord$vectors

#pcoa <- plot_ordination(Po_rarefied, 
                ord,
                color = "Locality", 
                label= "Locality") +
  geom_point(aes(shape = ORF_Type), size = 4) + 
  stat_ellipse(geom = "polygon", aes(group = Locality, label = Locality, fill = Locality), 
               alpha = 0.1, size = 1, linetype = 1) +
  geom_vline(xintercept = c(0), color = "grey", linetype = 1) +
  geom_hline(yintercept = c(0), color = "grey", linetype = 1) +
  labs(x = "PCo1 [33.7%]", y = "PCo2 [15.4%]") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #geom_point(aes(size=Shannon)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values=c("#bcbddc", "#7fcdbb","#fdae6b")) +
  scale_fill_manual(values=c("#bcbddc", "#7fcdbb","#fdae6b"))

print(pcoa)

#pcoa1 <- pcoa +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm"),
        strip.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ggtitle("ADONIS: p=0.001, R2=0.5865 \nNPMANOVA: p=0.001") 

  pcoa1

  #ggsave("/Users/Usuario/Desktop/Analisis_POC_16S/plot_PCoA_beta_2024.png", 
       pcoa1, width = 10, height = 10, dpi = 600)

###Tipo
#pcoa2 <- plot_ordination(Po_rarefied, 
                         ord,
                         color = "ORF_Type", 
                         label= "ORF_Type") +
   geom_point(aes(shape = Locality), size = 4) + 
  stat_ellipse(geom = "polygon", aes(group = ORF_Type, label = ORF_Type, fill = ORF_Type), 
               alpha = 0.1, size = 1, linetype = 1) +
  geom_vline(xintercept = c(0), color = "grey", linetype = 1) +
  geom_hline(yintercept = c(0), color = "grey", linetype = 1) +
  #labs(x = "PCo1 [62.0%]", y = "PCo2 [11.1%]") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #geom_point(aes(size=Shannon)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values=c("#ffeda0","#9ecae1")) +
  scale_fill_manual(values=c("#ffeda0","#9ecae1"))

print(pcoa2)

#pcoa3 <- pcoa2 +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm"),
        strip.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ggtitle("ADONIS: p=0.001, R2=0.5865 \nNPMANOVA: p=0.001") 

 print (pcoa3)

  #ggsave("/Users/Usuario/Desktop/Analisis_POC_16S/plot_PCoA_beta_2024.png", 
       pcoa1, width = 10, height = 10, dpi = 600)

###con los 6 grupos separados
sample_data(Po_rarefied)$Group <- interaction(sample_data(Po_rarefied)$Locality, 
                                              sample_data(Po_rarefied)$ORF_Type)
pcoagrp <- plot_ordination(Po_rarefied, 
                           ord,
                           color = "Locality", 
                           label = "Locality") +
  geom_point(aes(shape = ORF_Type), size = 4) + 
  stat_ellipse(geom = "polygon", 
               aes(group = Group, fill = Locality), 
               alpha = 0.2, linetype = 1) +
  geom_vline(xintercept = 0, color = "grey", linetype = 1) +
  geom_hline(yintercept = 0, color = "grey", linetype = 1) +
  labs(x = "PCo1 [33.7%]", y = "PCo2 [15.4%]") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = c("#bcbddc", "#7fcdbb", "#fdae6b")) +
  scale_fill_manual(values = c("#bcbddc", "#7fcdbb", "#fdae6b"))

print(pcoagrp)

pcoagrp1<-pcoagrp +
theme(axis.title = element_text(size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.ticks.length = unit(.25, "cm"),
      strip.text = element_text(size = 14),
      legend.key.size = unit(1, "cm"),
      plot.title = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
  #+ ggtitle("ADONIS: p=0.001, R2=0.5865 \nNPMANOVA: p=0.001") 

print(pcoagrp1)


#-----------

#----Hacer prueba estadistica de significancia con PERMANOVA (Adonis)---
#(actualizar o recargar los paquetes vegan y phyloseq por si no funciona)
#Hacer analisis de heterostasidad betadisper

adonis(dist ~ get_variable(Po_rarefied, "Locality"), permutations = 1000)$aov.tab
adonis(dist ~ get_variable(Po_rarefied, "ORF_Type"), permutations = 1000)$aov.tab

##----  #interaction Locality ORF type (nat)
all_test <- adonis(dist ~ get_variable(Po_rarefied, "Locality")+get_variable(Po_rarefied, "ORF_Type") + get_variable(Po_rarefied, "Locality"):get_variable(Po_rarefied, "ORF_Type"), permutations=1000)$aov.tab
all_test

#Pair-wise PERMANOVA para localidad 
cbn <- combn(x=unique(Metadata2$Locality), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  Po_subs <- subset_samples(Po_rarefied, Locality %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(Po_subs))
  permanova_pairwise <- adonis(phyloseq::distance(Po_subs, method = "bray") ~ Locality, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#Pair-wise PERMANOVA para Tipo ORF (nat)
cbn <- combn(x=unique(Metadata2$ORF_Type), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  Po_subs <- subset_samples(Po_rarefied, ORF_Type %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(Po_subs))
  permanova_pairwise <- adonis(phyloseq::distance(Po_subs, method = "bray") ~ ORF_Type, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#-----Pair-wise separando cada tipo (nat) ---------
#Pair-wise PERMANOVA para localidad (fuente)
#TIPO 1
cbn1 <- combn(x=unique(Metadata2$Locality), m = 2)
p1 <- c()

for(i in 1:ncol(cbn1)){
  Po_subs1 <- subset_samples(Po_rarefiedtipo1, Locality %in% cbn1[,i])
  metadata_sub1 <- data.frame(sample_data(Po_subs))
  permanova_pairwise1 <- adonis(phyloseq::distance(Po_subs1, method = "bray") ~ Locality, 
                               data = metadata_sub1)
  p1 <- c(p1, permanova_pairwise1$aov.tab$`Pr(>F)`[1])
}

p.adj1 <- p.adjust(p1, method = "BH")
p.table1 <- cbind.data.frame(t(cbn1), p=p1, p.adj=p.adj1)
p.table1
#--------------
#TIPO 3
cbn3 <- combn(x=unique(Metadata2$Locality), m = 2)
p3 <- c()

for(i in 1:ncol(cbn3)){
  Po_subs3 <- subset_samples(Po_rarefiedtipo3, Locality %in% cbn3[,i])
  metadata_sub3 <- data.frame(sample_data(Po_subs3))
  permanova_pairwise3 <- adonis(phyloseq::distance(Po_subs3, method = "bray") ~ Locality, 
                               data = metadata_sub3)
  p3 <- c(p3, permanova_pairwise3$aov.tab$`Pr(>F)`[1])
}

p.adj3 <- p.adjust(p3, method = "BH")
p.table3 <- cbind.data.frame(t(cbn3), p=p3, p.adj=p.adj3)
p.table3

#------------ separar cada localidad para ver la diferencias de los tipos
Po_rarefie_Gorgona
Po_rarefie_Utria
Po_rarefie_Cupica 
#-------------------Golfo cupica
cbn2 <- combn(x=unique(Metadata2$ORF_Type), m = 2)
p2 <- c()

for(i in 1:ncol(cbn2)){
  Po_subs2 <- subset_samples(Po_rarefie_Cupica, ORF_Type %in% cbn2[,i])
  metadata_sub2 <- data.frame(sample_data(Po_subs2))
  permanova_pairwise2 <- adonis(phyloseq::distance(Po_subs2, method = "bray") ~ ORF_Type, 
                                data = metadata_sub2)
  p2 <- c(p2, permanova_pairwise2$aov.tab$`Pr(>F)`[1])
}

p.adj2 <- p.adjust(p2, method = "BH")
p.table2 <- cbind.data.frame(t(cbn2), p=p2, p.adj=p.adj2)
p.table2
#-----------Utria
cbn4 <- combn(x=unique(Metadata2$ORF_Type), m = 2)
p4 <- c()

for(i in 1:ncol(cbn4)){
  Po_subs4 <- subset_samples(Po_rarefie_Utria, ORF_Type %in% cbn4[,i])
  metadata_sub4 <- data.frame(sample_data(Po_subs4))
  permanova_pairwise4 <- adonis(phyloseq::distance(Po_subs4, method = "bray") ~ ORF_Type, 
                                data = metadata_sub4)
  p4 <- c(p4, permanova_pairwise4$aov.tab$`Pr(>F)`[1])
}

p.adj4 <- p.adjust(p4, method = "BH")
p.table4 <- cbind.data.frame(t(cbn4), p=p4, p.adj=p.adj4)
p.table4
#-------------Gorgona
cbn5 <- combn(x=unique(Metadata2$ORF_Type), m = 2)
p5 <- c()

for(i in 1:ncol(cbn5)){
  Po_subs5 <- subset_samples(Po_rarefie_Gorgona, ORF_Type %in% cbn5[,i])
  metadata_sub5 <- data.frame(sample_data(Po_subs5))
  permanova_pairwise5 <- adonis(phyloseq::distance(Po_subs5, method = "bray") ~ ORF_Type, 
                                data = metadata_sub5)
  p5 <- c(p5, permanova_pairwise5$aov.tab$`Pr(>F)`[1])
}

p.adj5 <- p.adjust(p5, method = "BH")
p.table5 <- cbind.data.frame(t(cbn5), p=p5, p.adj=p.adj5)
p.table5


#Dispersion
disper1 <- boxplot(betadisper(dist, 
                   get_variable(Po_rarefied, "Locality")),las=2, 
        main=paste0("Multivariate Dispersion Test Bray-Curtis "," pvalue = ", 
                    permutest(betadisper(dist, get_variable(Po_rarefied, "Locality")))$tab$`Pr(>F)`[1]))
disper1

disper1reuslt <- permutest(betadisper(dist, get_variable(Po_rarefied, "Locality")))
disper1reuslt 



disper2 <- boxplot(betadisper(dist, 
                   get_variable(Po_rarefied, "ORF_Type")),las=2, 
        main=paste0("Multivariate Dispersion Test Bray-Curtis "," pvalue = ", 
                    permutest(betadisper(dist, get_variable(Po_rarefied, "ORF_Type")))$tab$`Pr(>F)`[1]))
disper2
disper2reuslt <- permutest(betadisper(dist, get_variable(Po_rarefied, "ORF_Type")))
disper2reuslt 

disper3 <- boxplot(betadisper(dist, 
                              get_variable(Po_rarefied, "Group")),las=2, 
                   main=paste0("Multivariate Dispersion Test Bray-Curtis "," pvalue = ", 
                               permutest(betadisper(dist, get_variable(Po_rarefied, "Group")))$tab$`Pr(>F)`[1]))
disper3
disper3reuslt <- permutest(betadisper(dist, get_variable(Po_rarefied, "Group")))
disper3reuslt 


#ANOSIM
plot(anosim(dist, get_variable(Po_rarefied, "Locality"))
     ,main="ANOSIM Bray-Curtis "
     ,las=2)
plot(anosim(dist, get_variable(Po_rarefied, "ORF_Type"))
     ,main="ANOSIM Bray-Curtis "
     ,las=2)
plot(anosim(dist, get_variable(Po_rarefied, "Group"))
     ,main="ANOSIM Bray-Curtis "
     ,las=2)


#-------------------------------------------------------------------------------
#ABUNDANCIA RELATIVA Y VISUALIZACION DE LA COMPOSICION DE LA COMUNIDAD MICROBIANA

#-------Visualizar y hacer analisis de la diversidad con base en la abundancia relativa
#Combinar todo en el nivel filo y hallar la abundancia relativa
#Hacer bar plot y box plot con la abundancia relativa para el filo y el genero

#BARRAS APILADAS, BOXPLOTS Y % DE ABUNDANCIA RELATIVA (TAXA RARA Y POR GRUPOS)
#Otra forma de hacer el analisis de la abundancia relativa basada en barras apiladas
#install.packages("ggsci")
library(ggsci)

Po1.rel = transform_sample_counts(Po1, function(x) x/sum(x)*100)

#Phylum
# agglomerate taxa
glom <- tax_glom(Po1.rel, taxrank = 'Phylum', NArm = FALSE)
Po1.melt <- psmelt(glom)


#Organizar los resultados en el orden que yo quiero 
Po1.melt$Locality <- factor(Po1.melt$Locality, levels = c("Golfo Cupica", "Utria", "Gorgona"))

# change to character for easy-adjusted level
Po1.melt$Phylum <- as.character(Po1.melt$Phylum)

Po1.melt <- Po1.melt %>%
  group_by(Locality, ORF_Type, Phylum) %>%
  mutate(median=median(Abundance))

#to get the same rows together
Po1.melt_sum <- Po1.melt %>%
  group_by(Sample, Locality, ORF_Type, Phylum) %>%
  summarise(Abundance=sum(Abundance))

###Para organizar en orden descendente (del filo más abundante)
# 1. Extraer el orden de las muestras según la abundancia de Endozoicomonadaceae
orden_muestras <- Po1.melt_sum %>%
  filter(Phylum == "Proteobacteria") %>%
  arrange(desc(Abundance)) %>%
  pull(Sample)

# 2. Reordenar el factor Sample en ese orden
Po1.melt_sum$Sample <- factor(Po1.melt_sum$Sample, levels = orden_muestras)

bar1 <- ggplot(Po1.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Locality + ORF_Type, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90)) +
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size= 14, color = "black"), 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 16),
        #axis.text.x.bottom = element_blank(),
        strip.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#9ecae1", "#636386", "#D89080",
                               "#bcbddc", "#c51b7d", "#fdae6b","#227487","#e0f3db", "#B1776E" ))
guides(fill = guide_legend(ncol=1))
bar1

#Phylum localidades
BoxPlot_Phylum1 <- Po1.melt_sum %>% 
  ggplot(aes(x = Phylum, y = Abundance, fill = Phylum)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~Locality, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#9ecae1", "#636386", "#D89080",
                               "#bcbddc", "#c51b7d", "#fdae6b","#227487","#e0f3db", "#B1776E"))
BoxPlot_Phylum1
#Phylum tipos
BoxPlot_Phylum2 <- Po1.melt_sum %>% 
  ggplot(aes(x = Phylum, y = Abundance, fill = Phylum)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~ORF_Type, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#9ecae1", "#636386", "#D89080","#bcbddc", "#c51b7d", "#fdae6b","#227487","#e0f3db", "#B1776E"))
BoxPlot_Phylum2

#Tabla Phylum
# Calcular media y error estándar de la abundancia relativa
tabla_phylum <- Po1.melt_sum %>%
  group_by(Locality, ORF_Type, Phylum) %>%
  summarise(
    Media = mean(Abundance),
    SE = sd(Abundance) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Media_SE = paste0(round(Media, 2), " ± ", round(SE, 2)))  # Formato para tabla

# Calcular media y error estándar
tabla_phylum <- Po1.melt_sum %>%
  group_by(Locality, ORF_Type, Phylum) %>%
  summarise(
    Media = mean(Abundance),
    SE = sd(Abundance) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Media_SE = paste0(round(Media, 2), " ± ", round(SE, 2)))

# Crear una columna combinada "Localidad_ORF"
tabla_phylum$Locality_ORF <- paste(tabla_phylum$Locality, tabla_phylum$ORF_Type, sep = "_")
#organizar en el orden que quiero ver los resutlados
tabla_phylum$Locality <- factor(tabla_phylum$Locality, levels = c("Golfo Cupica", "Utria", "Gorgona"))

# Reorganizar la tabla para que las combinaciones Locality_ORF sean columnas
tabla_phylum_final <- tabla_phylum %>%
  select(Phylum, Locality_ORF, Media_SE) %>%
  pivot_wider(names_from = Locality_ORF, values_from = Media_SE)

print(tabla_phylum_final)
write.csv(tabla_phylum_final, "tabla_abundancia_phylum_formato_publicacion.csv", row.names = FALSE)

#Familia
# agglomerate taxa
glom3 <- tax_glom(Po1.rel, taxrank = 'Family', NArm = FALSE)
Po1.melt3 <- psmelt(glom3)

#Organizar los resultados en el orden que yo quiero 
Po1.melt3$Locality <- factor(Po1.melt3$Locality, levels = c("Golfo Cupica", "Utria", "Gorgona"))

# change to character for easy-adjusted level
Po1.melt3$Family <- as.character(Po1.melt3$Family)

Po1.melt3 <- Po1.melt3 %>%
  group_by(Locality, ORF_Type, Family) %>%
  mutate(median=median(Abundance))

# select group median > 1.7
keep <- unique(Po1.melt3$Family[Po1.melt3$median >1.5])
Po1.melt3$Family[!(Po1.melt3$Family %in% keep)] <- "< 1.5"

#to get the same rows together
Po1.melt3_sum <- Po1.melt3 %>%
  group_by(Sample, Locality, ORF_Type, Family) %>%
  summarise(Abundance=sum(Abundance))
###Para organizar en orden descendente (primero Endozocomonadaceae)
# 1. Extraer el orden de las muestras según la abundancia de Endozoicomonadaceae
orden_muestras <- Po1.melt3_sum %>%
  filter(Family == "Endozoicomonadaceae") %>%
  arrange(desc(Abundance)) %>%
  pull(Sample)

# 2. Reordenar el factor Sample en ese orden
Po1.melt3_sum$Sample <- factor(Po1.melt3_sum$Sample, levels = orden_muestras)


bar3 <- ggplot(Po1.melt3_sum, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", aes(fill=Family)) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Locality + ORF_Type, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90)) +
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=14, color = "black"), 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 16),
        #axis.text.x.bottom = element_blank(),
        strip.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#227487", "#636386", "#D89080",
                    "#bcbddc", "#c51b7d", "#fdae6b", "#B1776E","#e0f3db", "#9ecae1"))
guides(fill = guide_legend(ncol=1))
bar3

#familia
BoxPlot_Family1 <- Po1.melt3_sum %>% 
  ggplot(aes(x = Family, y = Abundance, fill = Family)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~Locality, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#227487", "#636386", "#D89080",
                               "#bcbddc", "#c51b7d", "#fdae6b", "#B1776E","#e0f3db", "#9ecae1"))
BoxPlot_Family1

BoxPlot_Family2 <- Po1.melt3_sum %>% 
  ggplot(aes(x = Family, y = Abundance, fill = Family)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~ORF_Type, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#227487", "#636386", "#D89080",
                               "#bcbddc", "#c51b7d", "#fdae6b", "#B1776E","#e0f3db", "#9ecae1"))
BoxPlot_Family2

#Tabla Family
# Calcular media y error estándar de la abundancia relativa
tabla_Family <- Po1.melt3_sum %>%
  group_by(Locality, ORF_Type, Family) %>%
  summarise(
    Media = mean(Abundance),
    SE = sd(Abundance) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Media_SE = paste0(round(Media, 2), " ± ", round(SE, 2)))  # Formato para tabla

# Calcular media y error estándar
tabla_Family <- Po1.melt3_sum %>%
  group_by(Locality, ORF_Type, Family) %>%
  summarise(
    Media = mean(Abundance),
    SE = sd(Abundance) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Media_SE = paste0(round(Media, 2), " ± ", round(SE, 2)))

# Crear una columna combinada "Localidad_ORF"
tabla_Family$Locality_ORF <- paste(tabla_Family$Locality, tabla_Family$ORF_Type, sep = "_")
#organizar en el orden que quiero ver los resutlados
tabla_Family$Locality <- factor(tabla_Family$Locality, levels = c("Golfo Cupica", "Utria", "Gorgona"))

# Reorganizar la tabla para que las combinaciones Locality_ORF sean columnas
tabla_Family_final <- tabla_Family %>%
  select(Family, Locality_ORF, Media_SE) %>%
  pivot_wider(names_from = Locality_ORF, values_from = Media_SE)

print(tabla_Family_final)

#Genero
# agglomerate taxa
glom4 <- tax_glom(Po1.rel, taxrank = 'Genus', NArm = FALSE)
Po1.melt4 <- psmelt(glom4)

#Organizar los resultados en el orden que quiero ver los resultados
Po1.melt4$Locality <- factor(Po1.melt4$Locality, levels = c("Golfo Cupica", "Utria", "Gorgona"))

# change to character for easy-adjusted level
Po1.melt4$Genus <- as.character(Po1.melt4$Genus)

Po1.melt4 <- Po1.melt4 %>%
  group_by(Locality, ORF_Type, Genus) %>%
  mutate(median=median(Abundance))

# select group median > 1.7
keep <- unique(Po1.melt4$Genus[Po1.melt4$median > 0])
Po1.melt4$Genus[!(Po1.melt4$Genus %in% keep)] <- "<0%"

#to get the same rows together
Po1.melt4_sum <- Po1.melt4 %>%
  group_by(Sample, Locality, ORF_Type, Genus) %>%
  summarise(Abundance=sum(Abundance))

###Para organizar en orden descendente (primero Endozoicomonas)
# 1. Extraer el orden de las muestras según la abundancia de Endozoicomonadaceae
orden_muestras <- Po1.melt4_sum %>%
  filter(Genus == "Endozoicomonas") %>%
  arrange(desc(Abundance)) %>%
  pull(Sample)

# 2. Reordenar el factor Sample en ese orden
Po1.melt4_sum$Sample <- factor(Po1.melt4_sum$Sample, levels = orden_muestras)

bar4 <- ggplot(Po1.melt4_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", aes(fill=Genus)) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Locality + ORF_Type, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90)) +
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=14, color = "black"), 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 16),
        #axis.text.x.bottom = element_blank(),
        strip.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#636386","#227487", "#D89080",
                               "#bcbddc", "#c51b7d", "#fdae6b", "#B1776E","#e0f3db", "#bc80bd", "#9ecae1"))
guides(fill = guide_legend(ncol=1))
bar4

#Genero
BoxPlot_Genus1 <- Po1.melt4_sum %>% 
  ggplot(aes(x = Genus, y = Abundance, fill = Genus)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~Locality, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, , face = "italic"),
        axis.text.y = element_text(size = 15),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 17, face = "italic"),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 17),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#636386","#227487", "#D89080",
                               "#bcbddc", "#c51b7d", "#fdae6b", "#B1776E","#e0f3db", "#bc80bd", "#9ecae1"))
BoxPlot_Genus1

BoxPlot_Genus2 <- Po1.melt4_sum %>% 
  ggplot(aes(x = Genus, y = Abundance, fill = Genus)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~ORF_Type, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, , face = "italic"),
        axis.text.y = element_text(size = 15),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 17, face = "italic"),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 17),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_manual(values = c("#E64B35", "#ffeda0","#37B3BE", "#99d8c9","#636386","#227487", "#D89080",
                               "#bcbddc", "#c51b7d", "#fdae6b", "#B1776E","#e0f3db", "#bc80bd", "#9ecae1"))
BoxPlot_Genus2

#Tabla Genus
# Calcular media y error estándar de la abundancia relativa
tabla_Genus <- Po1.melt4_sum %>%
  group_by(Locality, ORF_Type, Genus) %>%
  summarise(
    Media = mean(Abundance),
    SE = sd(Abundance) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Media_SE = paste0(round(Media, 2), " ± ", round(SE, 2)))  # Formato para tabla

# Calcular media y error estándar
tabla_Genus <- Po1.melt4_sum %>%
  group_by(Locality, ORF_Type, Genus) %>%
  summarise(
    Media = mean(Abundance),
    SE = sd(Abundance) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Media_SE = paste0(round(Media, 2), " ± ", round(SE, 2)))

# Crear una columna combinada "Localidad_ORF"
tabla_Genus$Locality_ORF <- paste(tabla_Genus$Locality, tabla_Genus$ORF_Type, sep = "_")
#organizar en el orden que quiero ver los resutlados
tabla_Genus$Locality <- factor(tabla_Genus$Locality, levels = c("Golfo Cupica", "Utria", "Gorgona"))

# Reorganizar la tabla para que las combinaciones Locality_ORF sean columnas
tabla_Genus_final <- tabla_Genus %>%
  select(Genus, Locality_ORF, Media_SE) %>%
  pivot_wider(names_from = Locality_ORF, values_from = Media_SE)

print(tabla_Genus_final)

###-----MATRIZ DE TAXA X MUESTRAS PARA EL DIAGRAMA DE VENN---
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)

# 1. Crear variable combinada Localidad_Tipo en el objeto phyloseq original
sample_data(Po1)$Grupo <- paste(sample_data(Po1)$Locality, sample_data(Po1)$ORF_Type, sep = "_")

# 2. Agrupar el objeto phyloseq a nivel de Family
Po1_Fam <- tax_glom(Po1, taxrank = "Family")

# 3. Obtener tabla de conteos
otu_df <- as.data.frame(otu_table(Po1_Fam))
otu_df$ASV <- rownames(otu_df)

# 4. Extraer taxonomía
tax_df <- as.data.frame(tax_table(Po1_Fam))
tax_df$ASV <- rownames(tax_df)

# 5. Agregar el Family a la tabla OTU
otu_tax <- left_join(otu_df, tax_df[, c("ASV", "Family")], by = "ASV")

# 6. Convertir a formato largo
otu_long <- otu_tax %>%
  pivot_longer(cols = -c(ASV, Family), names_to = "SampleID", values_to = "Abundancia")

# 7. Crear tabla de metadatos con columna 'Grupo'
meta_df <- as(sample_data(Po1_Fam), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(Grupo = paste(Locality, ORF_Type, sep = "_"))

# 8. Asociar cada SampleID a su Grupo
otu_long <- left_join(otu_long, meta_df[, c("SampleID", "Grupo")], by = "SampleID")

# 9. Sumar abundancias por Family y Grupo
tabla_absoluta <- otu_long %>%
  group_by(Family, Grupo) %>%
  summarise(Abundancia = sum(Abundancia), .groups = "drop") %>%
  pivot_wider(names_from = Grupo, values_from = Abundancia, values_fill = 0)

# 10. Visualizar la tabla
View(tabla_absoluta)

###Diagrama de Venn sin necesitar la tabla anterior xd
#install.packages("VennDiagram")
library(VennDiagram)  # o VennDiagram si prefieres base
library(ggVennDiagram)
Po1_Genus <- tax_glom(Po1, taxrank = "Genus")

# 1. Crear variable combinada Localidad_Tipo en el objeto phyloseq original
sample_data(Po1)$Grupo <- paste(sample_data(Po1)$Locality, sample_data(Po1)$ORF_Type, sep = "_")

otup <- as.data.frame(otu_table(Po1_Genus))
taxp <- as.data.frame(tax_table(Po1_Genus))
metap <- as.data.frame(sample_data(Po1))
otup$Genus <- taxp$Genus
metap$SampleID <- rownames(metap)
otup$ASV <- rownames(otup)

otu_long <- otup %>%
  pivot_longer(cols = -c(ASV, Genus),
               names_to = "Grupo",
               values_to = "Abundancia") %>%
  left_join(meta %>% rownames_to_column("Grupo") %>% select(Grupo),
            by = "Grupo")

otu_long <- otu_long %>%
  mutate(Presencia = ifelse(Abundancia > 0, 1, 0))

venn_list <- otu_long %>%
  group_by(Genus, Grupo) %>%
  summarise(Presente = any(Presencia == 1), .groups = "drop") %>%
  filter(Presente == TRUE) %>%
  group_by(Grupo) %>%
  summarise(Phyla = list(unique(Genus))) %>%
  deframe()  # Convierte a lista nombrada

ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_void()

library(VennDiagram)
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = names(venn_list),
  filename = NULL,
  fill = c("#E69F00", "#56B4E9", "#009E73"),
  alpha = 0.5,
  cex = 1.5
)
grid::grid.draw(venn.plot)

# Encontrar géneros únicos por localidad
unique_genus <- venn_list

genus_unicos <- lapply(names(unique_genus), function(loc) {
  otros <- setdiff(names(unique_genus), loc)
  setdiff(unique_genus[[loc]], unlist(unique_genus[otros]))
})

names(genus_unicos) <- names(unique_genus)

# Ver en consola
genus_unicos

#-----MATRIZ DE MUESTRAS X genero PARA EL ANÁLISIS AMBIENTAL NAT----------
# Cargar paquetes necesarios
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
# 1. Agregar Localidad al nombre de muestra
#sample_data(Po1)$Sample_Loc <- paste0(sample_names(Po1), "_", sample_data(Po1)$Locality)
#sample_names(Po1) <- sample_data(Po1)$Sample_Loc

# 2. Agrupar por genero (sin transformar a abundancia relativa)
Po1.Genus <- tax_glom(Po1, taxrank = "Genus", NArm = TRUE)

# 3. Crear objeto en abundancia relativa para filtrar generos dominantes
Po1.Genus.rel <- transform_sample_counts(Po1.Genus, function(x) x / sum(x) * 100)

# 4. Obtener generos con abundancia relativa media > 1.5%
df_rel <- psmelt(Po1.Genus.rel)
generos_dominantes <- df_rel %>%
  group_by(Genus) %>%
  summarise(media_abund = mean(Abundance, na.rm = TRUE)) %>%
  filter(media_abund > 1.5) %>%
  pull(Genus)

# 5. Extraer matriz de conteo absoluto (OTU table)
asv_abs <- otu_table(Po1.Genus)
if (!taxa_are_rows(Po1.Genus)) {
  asv_abs <- t(asv_abs)
}
asv_abs <- as.data.frame(asv_abs)

# 6. Agregar la asignación taxonómica (genero)
taxa_info <- tax_table(Po1.Genus)[, "Genus"]
taxa_info <- as.data.frame(taxa_info)
taxa_info$ASV <- rownames(taxa_info)

# 7. Combinar con tabla de abundancia
asv_abs$ASV <- rownames(asv_abs)
asv_long <- pivot_longer(asv_abs, -ASV, names_to = "Sample", values_to = "Abundance")
asv_long <- left_join(asv_long, taxa_info, by = "ASV")

# 8. Filtrar solo las generos dominantes
asv_long_filt <- asv_long %>%
  filter(Genus %in% generos_dominantes)

# 9. Sumar por muestra y genero
matriz_abund_genero <- asv_long_filt %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

# 10. Ver resultado
head(matriz_abund)
matriz_abund_genero

write.csv(matriz_abund, "tabla_abundancia_generos_ambiente.csv", row.names = FALSE)
dca_result <- decorana(matriz_abund_genero)

# Ver el resumen
summary(dca_result)
#-------------------------------------------------------------------------------
#ANALISIS DIFERENCIAL DE LA ABUNDANCIA (DESeq2)

#Hacer un analisis diferencial de la abundancia con DEseq2 para detectar los grupos bacterianos por condicion
#Factorizar para DESeq2
sample_data(Po1)$Locality <- as.factor(sample_data(Po1)$Locality)
Po1.taxa <- tax_glom(Po1, taxrank = 'Family', NArm = FALSE)
#Comparacion pairwise
Po1.taxa.sub <- subset_samples(Po1.taxa, Locality %in% c("Golfo Cupica", "Utria", "Gorgona"))
#filtrar los caracteres con >90% de ceros
Po1.taxa.pse.sub <- prune_taxa(rowSums(otu_table(Po1.taxa.sub) == 0) < ncol(otu_table(Po1.taxa.sub)) * 0.9, Po1.taxa.sub)
Po1_ds = phyloseq_to_deseq2(Po1.taxa.pse.sub, ~ Locality)
#Estimador alternativo "cada gen contiene una muestra con un cero"
ds <- estimateSizeFactors(Po1_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05 
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res
taxa_sig = rownames(res[1:13, ]) # seleccionar 13 con los p.adj valores mas bajos
Po1.taxa.rel <- transform_sample_counts(Po1, function(x) x/sum(x)*100)
Po1.taxa.rel.sig <- prune_taxa(taxa_sig, Po1.taxa.rel)
#Hacer la matriz
Po1.taxa.rel.sig <- prune_samples(colnames(otu_table(Po1.taxa.pse.sub)), Po1.taxa.rel.sig)
matrix <- as.matrix(data.frame(otu_table(Po1.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(Po1.taxa.rel.sig)[, "Family"])
metadata_sub <- data.frame(sample_data(Po1.taxa.rel.sig))
colnames(matrix) = NULL

#Definir la anotacion
annotation_col = data.frame(
  Locality = as.factor(metadata_sub$Locality),
  ORF_Type = as.factor(metadata_sub$ORF_Type),
  check.names = FALSE)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(Phylum = as.factor(tax_table(Po1.taxa.rel.sig)[, "Phylum"]))
rownames(annotation_row) = rownames(matrix)

#Definir la anotacion de los colores
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
  Compartment = c(`Mucus` = "blue", `Tissue` = "red"),
  `Site` = c(`Protected` = "blue", `Urban` = "red"),
  `Season` = c(`Dry` = "blue", `Rainy` = "red"),
  `Phylum` = c(`Bacteroidota` = "red", `Cyanobacteria` = "white", `Firmicutes` = "black",
               `Proteobacteria` = "blue", `Spirochaetota` = "grey"))
library(circlize)
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
col_fun(seq(-3, 3))

Deseq_heatmap <- ComplexHeatmap::pheatmap(matrix, scale= "row", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors, 
                         col = col_fun, 
                         heatmap_legend_param = list( 
                           title = "Abundance",
                           legend_height = unit(4, "cm"),
                           title_position = "leftcenter-rot"))
Deseq_heatmap
                         
png(filename = "Deseq_Heatmap_family2.png", width = 10, height = 10, res = 1200, units = "in")
Heatmap <- draw(Deseq_heatmap, heatmap_legend_side = "left", annotation_legend_side = "right")
dev.off()

#Hacer un ANCOM-CB Analisis basado en la abundancia relativa para datos de microbioma
#el ANCOM-CB detecta que taxones cambiaron su abundancia relativa significativamente
#ancombc
library(tidyverse)
library(DT)
# Convertir abundancias a porcentajes para calcular las más abundantes
Po1.rel <- transform_sample_counts(Po1, function(x) x / sum(x) * 100)

# Agregar por Genus
glom_genus <- tax_glom(Po1.rel, taxrank = "Genus", NArm = TRUE)

# Derretir para tabla larga
df_genus <- psmelt(glom_genus)

df_genus <- as.data.frame(df_genus)
df_genus$Genus <- as.character(df_genus$Genus)
# Calcular abundancia media por género
library(dplyr)

str(df_genus)
top_genera <- df_genus %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
  filter(mean_abund > 0.7) %>%
  pull(Genus)
top_genera
keep_taxa <- taxa_names(Po1)[tax_table(Po1)[, "Genus"] %in% top_genera]

Po1.genus.filt <- prune_taxa(keep_taxa, Po1)

# Convertir Locality en factor y cambiar el nivel de referencia a "Gorgona"
sample_data(Po1.genus.filt)$Locality <- relevel(factor(sample_data(Po1.genus.filt)$Locality), ref = "Gorgona")

Po1.taxa.sub<- subset_samples(Po1.genus.filt, Locality %in% c("Golfo Cupica", "Utria", "Gorgona") &
                                                       ORF_Type %in% c("Type_1", "Type_3"))
out = ancombc(data = Po1.taxa.sub, assay_name = "counts", 
              tax_level ="Genus",  
              formula = "Locality + ORF_Type", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "Locality", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)
res <- out$res
res
res_global = out$res_global

#---
# 1. Extraer resultados
lfc <- as.data.frame(out$res$lfc)
pval <- as.data.frame(out$res$p_val)

# 2. Dar nombres a las columnas
colnames(lfc) <- paste0(colnames(lfc), "_logFC")
colnames(pval) <- paste0(colnames(pval), "_pval")

#3. Unir ambas tablas
resultados <- cbind(Taxon = rownames(lfc), lfc, pval)

# 4. Formatear p-valor significativo
resultados_fmt <- resultados %>%
  mutate(across(ends_with("_pval"), 
                ~ifelse(. < 0.05, paste0(formatC(., format = "f", digits = 2), " **"), formatC(., format = "f", digits = 2))))
# 5. Unir logFC con p-valores en columnas lado a lado
cols <- grep("_logFC", colnames(resultados_fmt), value = TRUE)
tabla_final <- resultados_fmt %>%
  select(Taxon, all_of(cols)) %>%
  mutate(across(starts_with("Taxon"), as.character))

for (col in cols) {
  var <- gsub("_logFC", "", col)
  tabla_final[[var]] <- paste0(
    formatC(resultados[[paste0(var, "_logFC")]], format = "f", digits = 2), " (",
    resultados_fmt[[paste0(var, "_pval")]], ")"
  )
}
# 7. Mostrar como tabla DT interactiva
library(DT)
datatable(tabla_final, 
          caption = "Diferencias en abundancia relativa (logFC) y valores p (entre paréntesis).",
          options = list(pageLength = 20)) 

#--
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "Utria vs Cupica", "Gorgona vs Cupica", "Type 3 vs Type 1")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")


#que bacterias contribuyen mas a las diferencias en la composicion de las comunidades entre los grupos (beta-diversidad)
#hacer un analisis PERMANOVA
pseq <- Po1

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

# samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ Locality,
                    data = meta, permutations=999, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["Locality", "Pr(>F)"])
#revisar la homogeneidad de condicion
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$Source))
permutest(betadisper(dist, meta$Source), pairwise = TRUE)
#taxa que contribuyen mas a las diferencias entre compartimentos o fuente (source)
coef <- coefficients(permanova)["Source1",]
top.coef <- coef[rev(order(abs(coef)))[1:40]]
par(mar = c(3, 6, 3, 6))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

#-------------------------------------------------------------------------------
#ANALISIS DEL CORE MICROBIOME (género)

#Definir el CORE microbioma (core-microbiome) y hacer analisis de taxones compartidos (Venn)
#El "core" se define como un grupo de taxones que se detectan en una marcada fraccion de la comunidad a un umbral de abundancia dada
library(microbiomeutilities)
library(RColorBrewer)
Po.rel <- microbiome::transform(Po1, "compositional")
Po.rel.f <- format_to_besthit(Po.rel)
#Set different detection levels and prevalence
prevalences <- seq(.05, 1, .05) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-2), log10(.2), length = 10)
#(1e-3) = 0.001% abundance; change "-3" to -2 to increase to 0.01%
core_plot <- plot_core(Po.rel.f, plot.type = "heatmap", 
               colours = rev(brewer.pal(10, "Spectral")),
               min.prevalence = 0.8, 
               prevalences = prevalences, 
               detections = detections) +
  xlab("Detection Threshold (Relative Abundance (%))") +
  theme_classic()

print(core_plot)

#Familia
Po.rel.f.fam <- aggregate_taxa(Po.rel.f, "Family")
prevalences <- seq(.05, 1, .05)
#detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)
detections <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)

core_plot1 <- plot_core(Po.rel.f.fam, 
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences, 
                detections = detections, min.prevalence = 0.75) +
  xlab("Relative Abundance (%)") 
core_plot1 <- core_plot1 + theme_classic() + ylab("Family") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14, color = "black"), 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 16),
        legend.key.height = unit(1, "cm"),
        axis.line = element_blank(),
        axis.ticks.length = unit(.25, "cm")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
core_plot1

#Genero
Po.rel.f.gen <- aggregate_taxa(Po.rel.f, "Genus")
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)

core_plot2 <- plot_core(Po.rel.f.gen, 
                plot.type = "heatmap", 
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences, 
                detections = detections, min.prevalence = .75) +
  xlab("Relative Abundance (%)")
core_plot2 <- core_plot2 + theme_classic() + ylab("Genus") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14, colour = "black"), 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 16),
        legend.key.height = unit(1, "cm"),
        axis.line = element_blank(),
        axis.ticks.length = unit(.25, "cm")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
core_plot2

figure_core <- ggarrange(core_plot1, core_plot2,
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1)

figure_core <- figure_core + theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))
figure_core 
  
ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_core_prevalence_2024.png", 
       figure_core, width = 16, height = 6, dpi = 600)

#Otra forma de determinar el "core"
Po1.core <- core(Po1, detection = .2/100, prevalence = 65/100)
Po1.core

core.taxa.standard <- core_members(Po.rel, detection = 0, prevalence = 65/100)
pseq.core <- core(Po.rel, detection = 0, prevalence = .65)
pseq.core
pseq.core2 <- aggregate_rare(Po.rel, "Family", detection = 0, prevalence = .65)
pseq.core2
core.taxa <- taxa(pseq.core)
core.taxa
core.abundance <- sample_sums(core(Po.rel, detection = .01, prevalence = .95))
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(Po.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")


#-------------------------------------------------------------------------------


