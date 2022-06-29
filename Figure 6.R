library(Seurat)
library(BiocNeighbors)
library(circlize)
library(ggridges)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)

#----Figure 6A----
load('all_data_harmony.Rdata')
pca_exp <- Embeddings(all_data_harmony,'pca')

OC <- Cells(subset(all_data_harmony,orign == 'OC'))
AC <- Cells(subset(all_data_harmony,orign == 'AC'))
ML <- Cells(subset(all_data_harmony,orign == 'ML'))
AW <- Cells(subset(all_data_harmony,orign == 'AW'))
NO <- Cells(subset(all_data_harmony,orign == 'NO'))

pca_OC <- pca_exp[OC,]
pca_AC <- pca_exp[AC,]
pca_ML <- pca_exp[ML,]
pca_NO <- pca_exp[NO,]
pca_AW <- pca_exp[AW,]

out_OC <- queryKNN(pca_OC, query=pca_AC, k=1,get.index=F)
out_ML <- queryKNN(pca_ML, query=pca_AC, k=1,get.index=F)
out_NO <- queryKNN(pca_NO, query=pca_AC, k=1,get.index=F)
out_AW <- queryKNN(pca_AW, query=pca_AC, k=1,get.index=F)

AC_orign <- as.data.frame(cbind(out_OC$distance,out_ML$distance,out_NO$distance,out_AW$distance))
rownames(AC_orign) <- rownames(pca_AC)
colnames(AC_orign) <- c('OC','ML','NO','AW')

AC_orign$celltype <- 'o'
AC_orign$fine_celltype <- 'o'

for (i in rownames(AC_orign)) {
  AC_orign[i,]$celltype <- as.character(all_data_harmony@meta.data[i,]$celltype)
  AC_orign[i,]$fine_celltype <- all_data_harmony@meta.data[i,]$fine_celltype
}

myeloid_AC <- AC_orign[AC_orign$celltype== 'Myeloid',]
c4_AC <- myeloid_AC[myeloid_AC$fine_celltype == 'mac_c4_PLIN2',]
c4_AC_orign <- apply(c4_AC, 1, function(x) colnames(c4_AC)[which.min(x)])

sectors <- factor(c('AC','NO','OC','AW','ML'))

circos.par("track.height"= 0.2)
circos.par("circle.margin"=c(1,1),ADD = T)
circos.initialize(sectors, xlim = c(0,200))

bg <- brewer.pal(5,'Set1')
circos.track(sectors, ylim = c(0, 1),bg.col = bg,bg.border = 'white',panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 
              CELL_META$cell.ylim[2] + mm_y(4), 
              CELL_META$sector.index,
              cex = 1.1,col ='black')
})

circos.link('AC',c(170,124.4),'OC',c(50,95.6),col='#CB6651')
circos.link('AC',c(120,79.4),'ML',c(100,140.6),col='#CB6651')
circos.link('AC',c(50,36.2),'AW',c(150,163.8),col='#CB6651')

circos.clear()

all_data_harmony@meta.data <-all_data_harmony@meta.data %>%  mutate(pt = case_when(
  orig.ident == 'OC1' ~ "P1",
  orig.ident == 'OC2' ~ "P2",
  orig.ident == 'OC3' ~ "P3",
  orig.ident %in% c('OC4','TC1') ~ "P4",
  orig.ident %in% c('OC5','AC1') ~ "P5",
  orig.ident %in% c('OC6','TC2','AC2') ~ "P6",
  orig.ident %in% c('OC7','TC3','AC3') ~ "P7",
  orig.ident == 'PF1' ~ "P8",
  orig.ident == 'PF2' ~ "P9",
  orig.ident == 'PF3' ~ "P10",
  orig.ident == 'NO1' ~ "P11",
  orig.ident == 'NO2' ~ "P12",
  orig.ident == 'NO3' ~ "P13"))

all_data_harmony <- subset(all_data_harmony,pt == 'P7')
ggplot(tmp,aes(x = myeloid_AC_orign,y = Freq,fill=myeloid_AC_orign))+
  geom_col()+
  xlab("Position")+
  ylab("Fraction of acites cells")+
  theme(panel.background=element_rect(fill='transparent', color='black',size = 2),panel.grid = element_blank())+
  theme(axis.title = element_text(size = 18,color="black"),axis.text = element_text(size = 14,color="black"))+
  theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF'),breaks = c('NO','OC','TC','PF'))+
  scale_y_continuous(labels = percent,limits = c(0,1))+labs(title = 'P7')

#----Figure 6B----
KEGG_name <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_exact_source,gs_description) %>% unique() %>% as.data.frame() %>% column_to_rownames(.,var = 'gs_exact_source')
KEGG_name <- KEGG_name[c(3,13,14,32,33,34,44,49,57,126,172),,drop = F]

KEGG_df <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_exact_source,gene_symbol)
s.sets <-  split(KEGG_df$gene_symbol,KEGG_df$gs_description)

robj <- get(load('myeloid.Rdata'))
data <- t(as.matrix(robj@assays$RNA@data))
ndata <- data[apply(data,1,function(x) !all(x==0)),]
tdata<-t(ndata)
tdata <- tdata[apply(tdata,1,function(x) !all(x==0)),]

gbm_es <- gsva(tdata,
               s.sets,
               mx.diff=TRUE, 
               verbose=TRUE,
               parallel.sz=1,
               kcdf="Gaussian")
t_gbm_es <- as.data.frame(t(gbm_es))
t_gbm_es$celltype <- 'unkown'
for (i in rownames(t_gbm_es)){
  t_gbm_es[i,]$celltype <- as.character(robj@meta.data[i,]$fine_celltype)
}
t_gbm_es <- t_gbm_es[,c(rownames(KEGG_name),'celltype')]
colnames(t_gbm_es)[1:11] <- KEGG_name$gs_description

score <- score[c(1,4,9,11,5,7,2,3,6,8,10),]

anonotation_row <- data.frame(class = c(rep('Adhesion',4)))  
rownames(anonotation_row) <- rownames(score)
ann_col <- list(class = c(Adhesion = "#E7695DFF"))

pheatmap::pheatmap(score,
         annotation_row = anonotation_row,
         annotation_colors = ann_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_legend = TRUE,
         border_color = 'white',
         scale = 'row',
         angle_col = 45,
        fontsize = 14,
         breaks = seq(-1,1,length.out = 100)
       ,cellwidth = 40,cellheight = 20)

#----Figure 6C----      
obj <- subset(myeloid,idents = levels(mac)[c(4,6)])
df <- obj@meta.data[,c('CAMs','fine_celltype')]
compare_means(CAMs~fine_celltype,df,'t.test')

scores.ridgeplot <- scores.ridgeplot$data %>%
  ggplot(aes(x = CAMs, y = ident, fill = ident)) +
  geom_density_ridges(jittered_points=TRUE, scale = .95, rel_min_height = .01,
                                point_shape = "|", point_size = 3, size = 0.25,
                                position = position_points_jitter(height = 0)) +
  scale_y_discrete(expand = c(.01, 0), name = "Cluster") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c('#CB6651','#FFBF19')) +
  guides(fill = guide_legend(
    override.aes = list(fill = c('#CB6651','#FFBF19'),  color = NA, point_color = NA))) +
  ggtitle('CAMs') +
  theme_ridges(center = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 12),
                 axis.title = element_text(size = 12),
                 panel.grid = element_blank())

#----Figure 6D----
marker <- FindAllMarkers(tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
a <-  marker[marker$cluster == 'Macrophage',c('avg_logFC','gene')]
ge <- a$avg_logFC
b <- select(org.Hs.eg.db, keys=a$gene, keytype="SYMBOL", columns=c("ENTREZID"))
names(ge) <- b$ENTREZID
ge <-  sort(ge,decreasing = T)

Go_gseresult <- gseGO(ge, 'org.Hs.eg.db', keyType = "SYMBOL", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
gseaplot2(Go_gseresult, c(790,787,736,672,643,640,610,603), pvalue_table = TRUE,rel_heights = c(2, 0.5, 0.5))

#----Figure 6F----
tmp <- read_excel('adhesion genes.xlsx')
tmp <- c(tmp$Group1,tmp$Group2,tmp$Group3)
tmp <- na.omit(tmp)
tmp <- toupper(tmp)

data <- GetAssayData(mac,slot = 'data')
data <- data[rownames(data)%in%tmp,]

data <- as.data.frame(t(data))
data$celltype <- 'o'
for (i in rownames(data)) {
  data[i,]$celltype <- as.character(myeloid@meta.data[i,]$fine_celltype)
}

mac_c1_C1QC <- data[which(data$celltype == 'mac_c1_C1QC'),] %>% subset(select = c(-celltype))
mac_c2_NLRP3 <- data[which(data$celltype == 'mac_c2_NLRP3'),]%>% subset(select = c(-celltype))
mac_c3_FN1 <- data[which(data$celltype == 'mac_c3_FN1'),]%>% subset(select = c(-celltype))
mac_c4_PLIN2 <- data[which(data$celltype == 'mac_c4_PLIN2'),]%>% subset(select = c(-celltype))
mac_c5_MALAT1 <- data[which(data$celltype == 'mac_c5_MALAT1'),]%>% subset(select = c(-celltype))
mac_c6_ISG15 <- data[which(data$celltype == 'mac_c6_ISG15'),]%>% subset(select = c(-celltype))
mac_c7_proliferation <- data[which(data$celltype == 'mac_c7_proliferation'),]%>% subset(select = c(-celltype))

mac_c1_C1QC_score <- as.data.frame(apply(mac_c1_C1QC, 2, mean))
colnames(mac_c1_C1QC_score) <- gsub('_score',' ','mac_c1_C1QC_score')
mac_c2_NLRP3_score <- as.data.frame(apply(mac_c2_NLRP3, 2, mean))
colnames(mac_c2_NLRP3_score) <- gsub('_score',' ','mac_c2_NLRP3_score')
mac_c3_FN1_score <- as.data.frame(apply(mac_c3_FN1, 2, mean))
colnames(mac_c3_FN1_score) <- gsub('_score',' ','mac_c3_FN1')
mac_c4_PLIN2_score <- as.data.frame(apply(mac_c4_PLIN2, 2, mean))
colnames(mac_c4_PLIN2_score) <- gsub('_score',' ','mac_c4_PLIN2')
mac_c5_MALAT1_score <- as.data.frame(apply(mac_c5_MALAT1, 2, mean))
colnames(mac_c5_MALAT1_score) <- gsub('_score',' ','mac_c5_MALAT1')
mac_c6_ISG15_score <- as.data.frame(apply(mac_c6_ISG15, 2, mean))
colnames(mac_c6_ISG15_score) <- gsub('_score',' ','mac_c6_ISG15')
mac_c7_proliferation_score <- as.data.frame(apply(mac_c7_proliferation, 2, mean))
colnames(mac_c7_proliferation_score) <- gsub('_score',' ','mac_c7_proliferation')

score <- merge(mac_c1_C1QC_score , mac_c2_NLRP3_score ,by = 'gene')
for (i in paste0(levels(mac),sep = '_score')[-c(1,2)]) {
  score <- merge(score , get(i) ,by = 'gene')
}
score <- remove_rownames(score)
score <- column_to_rownames(score,'gene')
score <- round(score,2) 
score <- score[rowSums(score)>0.05,]
nrow(score)
  
pheatmap::pheatmap(
    mat = score, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 14,border_color = 'white',
    main = 'Adhesion gene',
    angle_col = 45,
    display_numbers = F,
    scale = 'row',
    cellwidth = 115,cellheight = 15
  )
#----Figure 6G----
load('NKT.Rdata')
se <- c('FN1','CD40LG','L1CAM','PLA2G2A','FBN1','PSEN1','RhoGDI1','ICAM1','ICAM2','ICAM3','ICAM4','F11R'
        'ITGA4','ITGB1','ITGA9','ITGAD','ITGB2','ITGB7','PODXL','CD34','SELPLG','ITGAM','ITGAX','AREG','SPN','ITGAL')
score <- AverageExpression(NKT,features = rownames(NKT)[rownames(NKT)%in%se],slot = 'data') %>% .[['RNA']]
score <- score[se[se %in% rownames(score)],]
pheatmap::pheatmap(
  mat = score, 
  cluster_rows = T, 
  cluster_cols = F,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 14,border_color = 'white',
  breaks = seq(-2,2,length.out = 100),
  legend_breaks = c(-2,-1,0,1,2),
  main = 'Adhesion gene',
  angle_col = 45,
  gaps_row = c(4),
  display_numbers = F,
  scale = 'row',
  cellwidth = 40,cellheight = 15
)

#----Figure 6I----
Idents(all_data_harmony) <- all_data_harmony$fine_celltype
cd4 <- grep(pattern = "^CD4", all_data_harmony$fine_celltype, value = T)
cd8 <- grep(pattern = "^CD8", all_data_harmony$fine_celltype, value = T)
mac4 <- grep(pattern = "mac_c4_PLIN2", all_data_harmony$fine_celltype, value = T)

all_data_harmony$Tcell <- 'uk'
for (i in rownames(all_data_harmony@meta.data)) {
  all_data_harmony@meta.data[i,]$Tcell <- ifelse(i %in% names(cd4),'CD4',
                                                 ifelse(i %in% names(cd8), 'CD8',
                                                        ifelse(i %in% names(mac4),'mac_c4_PLIN2', 'nt'))))
}

tmp <- subset(all_data_harmony,Tcell %in% c('CD4','CD8','mac_c4_PLIN2') & orign == 'AC')
count <- as.matrix(tmp@assays$RNA@data)
count <- data.frame(Gene=rownames(count), count)
count <- remove_rownames(count)
write.table(count,'t-m_AC_count.txt', sep='\t', quote=F,row.names = F)

eta <- data.frame(Cell=rownames(tmp@meta.data), cell_type=tmp@meta.data$Tcell)
write.table(meta,'t-m_AC_meta.txt', sep='\t', quote=F, row.names=F)

immu_resp <- c(grep('FASLG$',intr_pairs,value = T),
               grep('^TNF.*TNF',intr_pairs,value = T),
               grep('^IFNG',intr_pairs,value = T),
               grep('^IL.*IL',intr_pairs,value = T))

immu_stim <- c(grep('^CSF.*CSF',intr_pairs,value = T), "CRTAM_CADM1","CD48_CD244","CD28_CD80","CD28_CD86","ICOSLG_ICOS","TNF_ICOS","TNFSF14_TNFRSF14","CD226_NECTIN2","PVR_CD226","CD27_CD70","CD226_NECTIN1","CD226_NECTIN3","CD226_NECTIN4","CD40LG_CD40","TNFRSF4_TNFSF4","TNFRSF8_TNFSF8", "TNFRSF9_TNFSF9","TNFRSF18_TNFSF18")

immu_inhib <-c("CD52_SIGLEC10","CD160_TNFRSF14","BTLA_TNFRSF14","CD80_CD274","SIRPA_CD47","CD96_NECTIN1","PVR_CD96","CTLA4_CD80","TIGIT_NECTIN2","PVR_TIGIT","CTLA4_CD86","CTLA4_CD80","LGALS9_HAVCR2","PDCD1_PDCD1LG2","LAIR1_PTPN6","LAIR1_PTPN11","PDCD1_CD274","TIGIT_NECTIN1","TIGIT_NECTIN3","TIGIT_NECTIN4")

kp <- grep(pattern = "mac_c4_PLIN2", colnames(pvalues), value= T)
ex <- c(grep(pattern = "^B", kp, value= T),grep(pattern = "B$", kp, value= T),
        grep(pattern = "DC", kp, value= T),grep(pattern = "Epithelial", kp, value= T),
        grep(pattern = "Fibroblast", kp, value= T),grep(pattern = "mono", kp, value= T),
        grep(pattern = "plasma B", kp, value= T))
pos <- setdiff(kp,ex) %>% .[c(1:11,15:25)]
pos <- which(colnames(pvalues) %in% pos)
choose_pvalues <- pvalues[,c(c(1,2,5,6,8,9),pos)]
choose_means <- means[,c(c(1,2,5,6,8,9),pos)]

logi <- apply(choose_pvalues[,6:ncol(choose_pvalues)]<0.05, 1, sum) 

choose_pvalues <- choose_pvalues[logi>=1,]

choose_pvalues <- choose_pvalues[choose_pvalues$interacting_pair %in% immu_resp,]

logi1 <- choose_pvalues$gene_a != ""
logi2 <- choose_pvalues$gene_b != ""
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi,]

choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]

meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(interacting_pair = paste0(meansdf$gene_a,"_",meansdf$gene_b),
                      CC = meansdf$variable,
                      means = meansdf$value)
pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(interacting_pair = paste0(pvalsdf$gene_a,"_",pvalsdf$gene_b),
                      CC = pvalsdf$variable,
                      pvals = pvalsdf$value)

pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >0))$means)
pldf%>% filter(means >0) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="#B2182B",mid = "#ffdbba",low = "#4393C3",midpoint = 1  )+ 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8))+
  labs(x = '',y = '')



