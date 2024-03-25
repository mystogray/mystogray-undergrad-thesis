library(vegan)
library(dplyr)
library(readr)
library(ggplot2)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(sysfonts)
library(glue)
library(MASS)
library(reshape)
library(reshape2)
library(showtext)
library(reshape2)
library(directlabels)

#Fonts
showtext_auto()
showtext_opts(dpi = 300)
font_add_google("Montserrat", "montserrat")
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
ann_colors = list(Kategori = c("Group A" = "#DC4379", "Group B" = "#11245E"))
color_kelimpahanspesies = randomcoloR::distinctColorPalette
colors35 = c("#67de7e", "#6838c8", "#80da3a", "#c450d8", "#dbd439", "#6a66d2", "#c5de73", "#622981", "#4f9b4b", "#d146a6", "#959631", "#322559", "#e09133", "#6a8ed9", "#da4e2c", "#65d4c3", "#c84074", "#b5d9a8", "#511b2f", "#d8b374", "#ca90d7", "#5b7439", "#da5159", "#50857b", "#843229", "#9ac0d9", "#a56c36", "#62749e", "#2e4224", "#d88fa1", "#2e2525", "#cbb1a0", "#394a5f", "#785c47", "#8b5076")
pos <- position_jitter(width = 0.3, seed = 2)
Kategori=c("Group A", "Group A", "Group A", "Group B", "Group B", "Group B", "Group B", "Group B", "Group B", "Group B", "Group B")

#Raw Data
allspecies_frac = allspecies %>% 
  dplyr::select(-contains('_num')) %>%
  rename_all(~toupper(sub("_frac", "", .))) %>%
  dplyr::select(-c(1,2))

#Alpha Diversity
alpha_div = diversity(as.data.frame(t(allgenus_frac))) %>%
  cbind.data.frame(., Kategori) %>%
  dplyr::rename("Value"=".")
ggplot(alpha_div, aes(x=Kategori, y=Value, colour=Kategori)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.8,position=pos) +
  ggpubr::stat_pvalue_manual(pvalue, label="{p.signif}", y.position = 4.3) +
  geom_text_repel(aes(label=rownames(alpha_div)),position=pos,size=4) +
  theme(panel.border=element_rect(fill=NA, color='black'),
        legend.text=element_text(face = "italic"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=18, face="bold"),
        plot.title = element_text(size=25, family="montserrat", face="bold")) +
  labs(title="Alpha Diversity Between Two Groups", y="Shannon Index", x="Category", color=NULL) +
  ggeasy::easy_center_title() +
  scale_color_manual(values=c("#DC4379","#11245E"))

#Beta Diversity
allgenus_frac = allgenus %>% 
  dplyr::select(-contains('_num')) %>%
  rename_all(~toupper(sub("_frac", "", .))) %>%
  dplyr::select(-c(1,2))

allgenus_frac_vis = t(allgenus_frac) %>%
  cbind.data.frame(metadata) %>%
  melt() %>%
  filter(value > 0.01)

dummy = group_by(allgenus_frac_vis,Sampel) %>%
  summarize(Sum_Value = sum(value, na.rm = TRUE)) %>%
  mutate("Genus dengan kelimpahan <1%" = 1-Sum_Value) %>%
  `[`(-2) %>%
  `[<-`(1:11, "V3", value = "Genus with less than 1% diversity") %>%
  `$<-`("Kategori", metadata$Kategori) %>%
  `[`(c(1,4,3,2)) %>%
  setNames(c("Sampel", "Kategori", "variable", "value"))

betadiv_vis = rbind.data.frame(allgenus_frac_vis,dummy)

ggplot(betadiv_vis, aes(fill=variable, y=value*100, x=Sampel)) + geom_bar(position="stack", stat="identity", width=0.7) +
  theme(panel.border=element_rect(fill=NA, color='black'),
        legend.text=element_text(face = "italic"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype=2, color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(size=21, family="montserrat", face="bold"),
        strip.text.x = element_text(size=13)) +
  labs(title="Bacterial Genus Diversity \nof Each Sample", y="Relative Abundance (%)", fill="Genus", x="Sample") +
  facet_grid(~Kategori,scales="free_x",space="free") +
  scale_fill_manual(values=unname(randomcoloR::distinctColorPalette(34))) +
  ggeasy::easy_center_title()

#Beta Diversity
groups = factor(c(rep(1,3), rep(2,8)), labels = c("Group A","Group B"))
betadiver = vegdist(t(allgenus_frac),index="Bray") %>%
  betadisper(group=metadata$Kategori)

betadiver_1 = as.data.frame(betadiver[["vectors"]]) %>%
  cbind.data.frame("Kategori"=metadata$Kategori)

ggplot(betadiver_1, aes(x=PCoA1, y=PCoA2)) +
  geom_point(size=4,aes(colour=Kategori)) +
  geom_text_repel(label=rownames(betadiver_1)) +
  ggtitle("Beta Diversity \nof Each Sample") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  coord_fixed(ratio=1) +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(panel.grid=element_blank(),
        axis.title = element_text(size=18, face="bold"),
        plot.title = element_text(size=25, family="montserrat", face="bold"),
        axis.text = element_text(size=14)) +
  labs(color=NULL) +
  scale_color_manual(values=c("#DC4379","#11245E")) +
  ggeasy::easy_center_title()
