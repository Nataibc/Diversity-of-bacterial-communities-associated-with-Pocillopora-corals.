####https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
#vasc <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/vasc_plants.txt', row.names = 1)
#chem <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/chemistry.txt', row.names = 1)

library(vegan)
require(readxl)

setwd("/Users/Usuario/Desktop/RDAgenero")
Metadata <- read_excel("Datos_ambientales.xlsx", sheet = "Metadata")#datos1
dataAbunAbsS <- read_excel("Datos_ambientales.xlsx", sheet = "Generos")#datos1
dataAbun<-dataAbunAbsS[,-1]
rownames(dataAbun) <- dataAbunAbsS$Samples
summary(dataAbun)

dataVar <- read_excel("Datos_ambientales.xlsx", sheet = "Ambiente")#datos2
Var<-dataVar[,-1]
rownames(Var) <- dataVar$Samples
summary(Var)

# transform data using Hellinger transformation to prepare them for tb-RDA:
# Hellinger transformation converts spp. abundances from absolute to relative values
bact.hell <- decostand (dataAbun, 'hell')

###  Square root transfomartion
#str(dataAbun)
### root transformation
#bac_familias_transf <- sqrt(dataAbun)
#rownames(bac_familias_transf) <- dataAbunAbsS$Samples

#set of (multiple) linear regression analyses, where species abundances (for each species in the species composition matrix 
#separately) are regressed against (one or several) environmental variable(s). The result is that variation in species 
#composition is decomposed into variation related to environmental variables (represented by constrained/canonical axes) and 
#not related to environmental variables (unconstrained axes). While in the case of unconstrained ordination the information 
#we are interested is mostly about the configuration of samples and species in the ordination diagram, the relative importance
#of individual ordination axes (measured by their eigenvalues) and ecological interpretation of ordination axes, in the case 
#of constrained ordinations we are more interested in the effect of environmental variables on species composition, namely 
#in the amount of variation these variables explain and whether this variation is significant or not
#transformation-based RDA ordination (tb-RDA). First, we need to test the global model (including all explanatory 
#variables from which you plan to do selection) to see whether it is significant; if it is not, the forward selection should not be done.
tb_rda.all <- rda (bact.hell ~ ., Var)
anova (tb_rda.all)  #   0.001 *** - it is significant
adjR2.tbrda <- RsquareAdj (tb_rda.all)$adj.r.squared 
adjR2.tbrda # adjusted R2 explained by all variables is 20.27764

##### a) Use of forward.sel function
library (adespatial)
sel.fs <- forward.sel (bact.hell, Var, nperm = 49999)
pval.adj <- p.adjust (sel.fs$pval, method = 'bonferroni', n = ncol (Var))
sel.fs$pval.adj <- pval.adj
sel.fs$pval.adj.stars <- gtools::stars.pval (pval.adj)
sel.fs


#### b) use ordi2step
#the criteria for including the variable is based on both significance of the newly selected variables, and the comparison of
#adjusted variation (R2adj) explained by the selected variables to R2adj explained by the global model (with all variables); 
#if the new variable is not significant or the R2adj of the model including this new variable would exceed the R2adj of the ç
#global model, the selection will be stopped.
tb_rda.symb.0 <- rda (bact.hell ~ 1, data = Var) # model containing only species matrix and intercept
tb_rda.symb.all <- rda (bact.hell  ~ ., data = Var) # model including all variables from matrix chem1 (the dot after tilda (~) means "include all from data")
sel.osR2 <- ordiR2step (tb_rda.symb.0, scope = formula (tb_rda.symb.all), R2scope = adjR2.tbrda, direction = 'forward', permutations = 4999)
sel.osR2
sel.osR2$anova

sel.osR2_adj <- sel.osR2
sel.osR2_adj$anova$`Pr(>F)` <- p.adjust (sel.osR2$anova$`Pr(>F)`, method = 'bonferroni', n = ncol (Var))
sel.osR2_adj$anova

#c) Use of ordistep function
#step-wise selection of environmental variables based two criteria: if their inclusion into the model leads to significant 
#increase of explained variance (the same as in forward.sel and ordiR2step), and if the AIC of the new model is lower than 
#AIC of the more simple mode

rda.symb.0 <- rda (bact.hell ~ 1, data = Var) # model containing only species matrix and intercept
rda.symb.all <- rda (bact.hell ~ ., data = Var) # model including all variables from matrix chem1 (the dot after tilda (~) means ALL!)
sel.os <- ordistep (rda.symb.0, scope = formula (rda.symb.all), direction = 'forward')
sel.os
sel.os$anova
summary(sel.os)

#### Graph
#install.packages("BiodiversityR")
library(BiodiversityR) # also loads vegan
library(ggplot2)
library(readxl)
library(dplyr)
library(ggsci)
library(ggrepel)
library(ggforce)
library(tibble)
library(grid)

plot<-ordiplot (sel.osR2_adj, display = c('Sites', 'bp'), type = 't')
biplot_data<-scores(sel.osR2_adj)$biplot

biplot_data1<- data.frame(
  name = c("Chlo-a", "CDOM", "SST"),
  RDA1 = c(-0.3347648,-0.6762569, 0.5961254),
  RDA2 = c(-0.8483552, -0.1292929, 0.7342318))

biplot_data2<- data.frame(
  name = c("SSSprom", "CDOM", "SST", "Anom inf"),
  RDA1 = c(0.0005803718,-0.6342298603, 0.5726607420,0.6050617934),
  RDA2 = c(-0.5981892, -0.2578197, 0.7081072, 0.5394763))

summary(sel.osR2)$cont$importance[2, 1:2]
 
scores(sel.osR2)$sites %>%
  as_tibble(rownames = "Samples") %>%
  inner_join(., Metadata, by="Samples") %>%
  ggplot(aes(x=RDA1, y=RDA2, color=Localidad)) +
  geom_point(aes(shape=Localidad), size=3)+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm")) +
  stat_ellipse(aes(fill=Localidad, color = Localidad),geom="polygon" ,level=0.95, alpha=0.2, show.legend = F)+
  scale_color_manual(values=c("Golfo Cupica" = "#bcbddc","Gorgona"= "#7fcdbb","Utria"= "#fdae6b"))+
  scale_fill_manual(values=c("Golfo Cupica" = "#bcbddc","Gorgona"= "#7fcdbb","Utria"= "#fdae6b"))+
  geom_segment(data=biplot_data1, aes(x=0, xend=RDA1, y=0, yend=RDA2), inherit.aes = FALSE)+
  #geom_text(data=biplot_data1, aes(x=RDA1, y=RDA2, label=name), inherit.aes = FALSE)+
  labs(x = "RDA1 (19,13%)", y = "RDA2 (2,35%)")+
  geom_text_repel(data = biplot_data1,
                  aes(x = RDA1 * 1.1 , y = RDA2 * 1.2, label = name),
                  inherit.aes = FALSE,
                  size = 4, # tamaño de texto
                  segment.color = "grey50")+
  coord_fixed() 


ggsave("~/Downloads/td-rda-var.sel.svg") 


