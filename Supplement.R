library(clusterProfiler)
library(Seurat)
library(ggplot2)

#----sFigure 1A/B----
df <- as.data.frame(table(all_data_harmony$orig.ident))
df <- df[order(df$Freq),] %>% mutate(id = seq(1,19,1))
df$Var1 <- factor(df$Var1,levels = df$Var1)
df <- mutate(df,label = case_when(
  id <= 3 ~ paste0(df[id,1],'\n',Freq) ,
  id <= 9 ~ paste0(Freq, '\n',df[id,1]),
  T ~ paste0(Freq,' ',df[id,1])
))

df <- as.data.frame(table(all_data_harmony$orign))
df <- rbind(df,data.frame(Var1=c('n','a'),Freq = c(0,0)))
df <- df[order(df$Freq),] %>% mutate(id = seq(1,7,1))
df$Var1 <- factor(df$Var1,levels = df$Var1)

ggplot(df, aes(Var1,Freq, fill =id,label=label)) +
  geom_col(width = 1, color = 'white') +
  geom_col(aes(y = 400), width = 1, alpha = 0.2, fill = 'white') + 
  geom_col(aes(y = 200), width = 1, alpha = 0.2, fill = 'white') +
  coord_polar(theta = 'x',start = 0,direction = 1) + 
  theme_void() + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  xlab(label = '')+ylab(label = '')+
  scale_fill_gradientn(colors = c("#54778f", "#4EB043", "#E69D2A", "#DD4714", "#A61650"),guide = 'none')

ggplot(df, aes(Var1,Freq, fill =Var1)) +
  geom_col(width = 1, color = 'white') +
  coord_polar(theta = 'x',start = 0,direction = -1) + 
  theme(axis.text.x=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.text.y.left = element_text(size = 10,colour = 'black'))+
  xlab(label = '')+ylab(label = '')+
  scale_fill_manual(values = cl)

#----sFigure 1B----
Idents(all_data_harmony) <- all_data_harmony$orign
deg <- FindAllMarkers(all_data_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
deg <-  dplyr::group_by(deg,cluster) %>%  as.data.frame
list <- list(OC = deg[deg$cluster == 'OC','gene'],NO = deg[deg$cluster == 'NO','gene'],
             PF = deg[deg$cluster == 'PF','gene'],AC = deg[deg$cluster == 'AC','gene'],
             TC = deg[deg$cluster == 'TC','gene'])

upset(fromList(list), sets.bar.color = cl[c(4,2,1,5,3)],sets.x.label = NULL,point.size = 7,
      queries = list(list(query = intersects, params = 'NO', color = cl[1],active = T),
                                    list(query = intersects, params = 'OC', color = cl[2],active = T),
                                    list(query = intersects, params = 'TC', color = cl[3],active = T),
                                    list(query = intersects, params = 'PF', color = cl[4],active = T),
                                    list(query = intersects, params = 'AC', color = cl[5],active = T)))

#----sFigure 2G----
M <- table(mye_meta[,c('orign','fine_celltype')])
(Xsq <- chisq.test(M))
mye_meta_ob=as.data.frame(table(mye_meta[,c('orign','fine_celltype')]))
mye_meta_expect=as.data.frame(Xsq$expected)
mye_meta_ob$expect=0
for (i in unique(mye_meta_ob$orign)) {
  for (j in unique(mye_meta_ob$fine_celltype)) {
    mye_meta_ob[mye_meta_ob$orign==i&mye_meta_ob$fine_celltype==j,]$expect=mye_meta_expect[i,j]
  }
}
mye_meta_ob$OE=mye_meta_ob$Freq/mye_meta_ob$expect

result_df=data.frame() 
for (i in unique(mye_meta_ob$orign)) {
  for (j in unique(mye_meta_ob$fine_celltype)) {
    result_df[i,j]=mye_meta_ob[mye_meta_ob$orign==i&mye_meta_ob$fine_celltype==j,]$OE
  }
}

tmp <- result_df[,levels(mac)]

df <- data.frame()
for (i in levels(mac)) {
  df[i,1] <- tmp[1,i]/tmp[4,i]
}

df <- rownames_to_column(df,var = 'cell')
df$cell <- factor(df$cell,levels = sort(levels(mac)))
ggplot(df,aes(cell,V1))+
  geom_col(aes(fill = cell))+
  labs(y = 'The ratio of Ro/e (AC/PF)',x = '')+
  theme(panel.background = element_rect(fill = 'transparent'),axis.line = element_line(size = 0.7),
        axis.title.y = element_text(size =14,colour = 'black'),axis.text = element_text(size = 12,colour = 'black'),legend.position = 'none',axis.text.x = element_text(angle = 45,vjust = 0.62))+
  scale_fill_manual(values = cl)

#----sFigure 3A----
tmp <- features.scores[features.scores$orign == 'PF',]
ggplot(tmp,aes(M1,M2))+geom_point(aes(color = orign),alpha = 0.5)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed",lwd = 2)+
  scale_color_manual(values = cl)+
  theme(panel.background = element_rect(fill = 'transparent',colour = 'black',linetype = 'solid',size = 2),
        panel.grid = element_blank(),legend.key = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),axis.title.y = element_text(size = 18),axis.title.x = element_text(size = 18))

#----sFigure 3B----
M1_signature <- read.csv('M1_signature.csv')
M2_signature <- read.csv('M2_signature.csv')

tmp <- GetAssayData(subset(mac,fine_celltype == 'mac_c4_PLIN2'),slot = 'data')
ge <- c(as.character(M1_signature$M1)[-c(24:37)],as.character(M2_signature$M2))[c(as.character(M1_signature$M1)[-c(24:37)],as.character(M2_signature$M2)) %in% rownames(tmp)]
tmp <- tmp[ge,]
or <- as.data.frame(rowMeans(tmp)) %>% rownames_to_column(var = 'gene')
or$group <- c(rep('M1',20),rep('M2',42))
or <- or %>% arrange(group,desc(rowMeans(tmp))) %>% as.data.frame()
tmp <- tmp[or$gene,]
pheatmap(tmp,border_color = 'white',cluster_rows = F,cluster_cols = F,fontsize = 12,show_colnames = F,
         scale = 'none',color = colorRampPalette(brewer.pal(9,'Reds'))(100),gaps_row = c(20))

#----sFigure 3C----
tmp <- subset(myeloid,idents = c('mac_c4_PLIN2','mac_c6_ISG15'))
df <- data.frame(exp = tmp@assays$RNA@data['gene',],cell = tmp$fine_celltype)
ggplot(df,aes(cell,exp))+
    geom_violin(aes(fill = cell),scale = 'width')+
    eom_boxplot(width = 0.1)+
    scale_fill_manual(values = c("#CB6651", "#FFBF19")) + 
    labs(title = 'gene',y = 'Expression', x = '')+
    theme(panel.background = element_rect(fill = 'transparent'),
                axis.line.x.bottom = element_line(color = 'black'),
                axis.line.y.left = element_line(color = 'black'),
                axis.title = element_text(size = 14,color="black"),
                axis.text = element_text(size = 12,color="black"),
                legend.key =element_blank(),
                legend.text = element_text(size = 12),
                legend.title = element_blank(),
                plot.title = element_text(size = 15,hjust = 0.5))+
          stat_compare_means(label="p.signif",label.y = 4,label.x = 1.35,bracket.size = 0.5,size =10)
  )

#----sFigure 3D----
mac_deg <- FindAllMarkers(mac, min.pct = 0.25, logfc.threshold = 0.25,only.pos = F)#only.pos默认值为F，此时logfc.threshold为绝对值

for (i in cluster) {
  mac_cluster[[i]] <- mac_deg[[i]]
  mac_cluster[[i]]$threshold = factor(ifelse(mac_cluster[[i]]$p_val_adj < 0.01 & abs(mac_cluster[[i]]$avg_logFC) >= 0.5, ifelse(mac_cluster[[i]]$avg_logFC>= 0.5 ,'Up','Down'),'NS'),levels=c('Up','Down','NS'))
}
mac_cluster <- do.call('rbind',mac_cluster)
mac_cluster <- mac_cluster[mac_cluster$threshold == 'Up',]
mac_cluster <- mac_cluster[!duplicated(mac_cluster$gene),]

for (i in cluster) {
  mac_cluster[[i]] <- mac_cluster[[i]][mac_cluster[[i]]$threshold != 'NS',]
  gene_id <- mac_cluster[[i]]$gene
  tmp <- bitr(gene_id,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db,drop = FALSE)
  mac_cluster[[i]]$ENTREZID <- tmp$ENTREZID
}

Up <- list()
result_up <- list()
for (i in cluster) {
  Up[[i]] <- mac_cluster[[i]][mac_cluster[[i]]$threshold == 'Up',]
  gene_entrezid <- na.omit(Up[[i]]$ENTREZID)
  result_up[[i]] <- enrichGO(gene          = gene_entrezid,
                             OrgDb         = 'org.Hs.eg.db',
                             keyType       = 'ENTREZID',
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.01)
}

Dn <- list()
result_dn <- list()
for (i in cluster) {
  Dn[[i]] <- mac_cluster[[i]][mac_cluster[[i]]$threshold == 'Down',]
  gene_entrezid <- na.omit(Dn[[i]]$ENTREZID)
  result_dn[[i]] <- enrichGO(gene          = gene_entrezid,
                             OrgDb         = 'org.Hs.eg.db',
                             keyType       = 'ENTREZID',
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.01)
}

result_all <- list()
for (i in cluster[-5]) {
  result_up[[i]] <- as.data.frame(result_up[[i]])
  result_up[[i]]$group <- 'Up'
  result_up[[i]]$p.adjust <- -log10(result_up[[i]]$p.adjust)
  result_up[[i]] <- result_up[[i]][1:10,] 
  result_dn[[i]] <- as.data.frame(result_dn[[i]])
  result_dn[[i]]$group <- 'Dn'
  result_dn[[i]]$p.adjust <- log10(result_dn[[i]]$p.adjust)
  result_dn[[i]] <- result_dn[[i]][1:10,] 
  result_all[[i]] <- rbind(result_up[[i]], result_dn[[i]] )
}

library(ggpubr)
for (i in cluster[-5]) {
  ggbarplot(result_all[[i]], x = "Description", y = "p.adjust",
            fill = "group",           
            color = "white",            
            palette = "lancet",         
            sort.val = "desc",          
            sort.by.groups = TRUE,      
            x.text.angle = 90,         
            ylab = "-/+ log10(Pvalue)",
            xlab = 'Top 10 GO terms',
            main = i,
            legend.title = " ",
            rotate = TRUE,
            ggtheme = theme_minimal(),
            font.tickslab = c(10,'plain','black'),
            font.main = c(18,'bold','black'),
            font.x = c(14,'bold','black'),
            font.y = c(14,'bold','black'))

#----sFigure4 A-D----
position_deg_MLvAC <- FindMarkers(mac, ident.1 = 'AC', ident.2 = 'ML',logfc.threshold = 0.25)
position_deg_MLvAC <- rownames_to_column(position_deg_MLvAC,var = 'gene')
position_deg_MLvAC$threshold = factor(ifelse(position_deg_MLvAC$p_val_adj < 0.001 & abs(position_deg_MLvAC$avg_logFC) >= 1, 
                                  ifelse(position_deg_MLvAC$avg_logFC>= 1 ,'AC','NL'),'NS'),levels=c('AC','ML','NS'))

position_deg_OCvAC <- FindMarkers(mac, ident.1 = 'AC', ident.2 = 'OC',logfc.threshold = 0.25)
position_deg_OCvAC <- rownames_to_column(position_deg_OCvAC,var = 'gene')
position_deg_OCvAC$threshold = factor(ifelse(position_deg_OCvAC$p_val_adj < 0.001 & abs(position_deg_OCvAC$avg_logFC) >= 1, 
                                             ifelse(position_deg_OCvAC$avg_logFC>= 1 ,'AC','OC'),'NS'),levels=c('AC','OC','NS'))

position_deg_AWvAC <- FindMarkers(mac, ident.1 = 'AC', ident.2 = 'AW',logfc.threshold = 0.25)
position_deg_AWvAC <- rownames_to_column(position_deg_AWvAC,var = 'gene')
position_deg_AWvAC$threshold = factor(ifelse(position_deg_AWvAC$p_val_adj < 0.001 & abs(position_deg_AWvAC$avg_logFC) >= 1, 
                                             ifelse(position_deg_AWvAC$avg_logFC>= 1 ,'AC','AW'),'NS'),levels=c('AC','PF','NS'))

position_deg_NOvAC <- FindMarkers(mac, ident.1 = 'AC', ident.2 = 'NO',logfc.threshold = 0.25)
position_deg_NOvAC <- rownames_to_column(position_deg_NOvAC,var = 'gene')
position_deg_NOvAC$threshold = factor(ifelse(position_deg_NOvAC$p_val_adj < 0.001 & abs(position_deg_NOvAC$avg_logFC) >= 1, 
                                             ifelse(position_deg_NOvAC$avg_logFC>= 1 ,'AC','NO'),'NS'),levels=c('AC','NO','NS'))

ggplot(position_deg_MLvAC,aes(x=avg_logFC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=cl)+
  geom_text_repel(
    data = filter(position_deg_MLvAC,threshold != 'NS'),
    aes(label = gene),
    size = 5,
    segment.color = "black", show.legend = FALSE )+
  theme(panel.background = element_rect(fill = 'white', colour = 'black',size = 1))+
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18,hjust = 0.5),
        axis.title = element_text(size = 18),axis.text = element_text(size = 14),
        panel.grid = element_blank())+
  ylab('-log10 (p-adj)')+
  xlab('log2 (FoldChange)')+
  labs(title = 'AC vs ML')+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=1)+
  geom_hline(yintercept = -log10(0.001),lty=3,col="black",lwd=1)

#----sFigure 4E----
Pt5 <- subset(mac,idents = c('OC5','AC1'))
Pt6 <- subset(mac,idents = c('OC6','AC2','TC2'))
Pt7 <- subset(mac,idents = c('OC7','AC3','TC3'))

df <- list()
for (i in c('Pt5','Pt6','Pt7')) {
df[[i]] <- GetAssayData(get(i),slot = 'data') %>%.[c('PLIN2'),] 
df[[i]] <- as.matrix(df[[i]]) %>% as.data.frame() 
df[[i]]$orign <- get(i)$orign
df[[i]]$pt <- get(i)$pt
names(df[[i]])[1] <- 'PLIN2'
}
df <- do.call('rbind',df)
data.frame(compare_means(PLIN2~orign,df,group.by = 'pt'))
ggplot(df,aes(pt,PLIN2,fill = orign))+
  stat_summary(fun = "mean", geom = "bar", 
              position = position_dodge(0.75),width = 0.75) +
   stat_summary(fun = mean,
                geom = "errorbar",width = 0.5,size = 0.3,
                fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
                fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),position = position_dodge(0.75))+
  labs(x ='',y = 'PLIN2 expression')+
  scale_fill_manual(values = cl[c('AC','TC','OC')])+
  theme(panel.background=element_rect(fill='transparent', color='black',size = 2),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14,color="black"),
        axis.text = element_text(size = 12,color="black"),
        legend.key =element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 15,hjust = 0.5))

#----sFigure 4F----
mac <- GetAssayData(tmp,slot = 'data') %>% .[c('C1QC','C1QB','CD163'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
fibr <- GetAssayData(tmp,slot = 'data') %>% .[c('COL1A2','ACTA2','DCN'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
can <- GetAssayData(tmp,slot = 'data') %>% .[c('EPCAM','KRT17','MMP7'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
B <- GetAssayData(tmp,slot = 'data') %>% .[c('CD19','CD79A'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
Tc <- GetAssayData(tmp,slot = 'data') %>% .[c('CD3D','CD3E'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
DC <- GetAssayData(tmp,slot = 'data') %>% .[c('HLA-DQA2','FCER1A','CD1C'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()

umap_pos=Embeddings(tmp,'tsne') 
a <- cbind(DC,umap_pos) %>% as.data.frame()
p6 <- ggplot(a,aes(x = tSNE_1,y = tSNE_2))+
  geom_point(aes(color = .),size = 0.85)+
  labs(title = 'Dentritic cell',subtitle = '(HLA-DQA2,FCER1A,CD1C)')+
  guides(color = F)+
  scale_color_gradientn(colors = rev(brewer.pal(n=11, name = 'Spectral')))+
  theme(panel.background = element_rect(fill = 'transparent',colour = 'black',linetype = 'solid',size = 2),
        axis.ticks = element_blank(),axis.title = element_blank(),axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 18),plot.subtitle = element_text(hjust = 0.5,size = 14),
        legend.title = element_text(size = 14),panel.grid = element_blank())
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(guides='collect',heights = c(1,1),nrow = 2)#library(patchwork)

#----sFigure6 A----
data <- data.frame()
for (i in c(levels(myeloid)[c(8,10,12,14)])) {
  tmp <- myeloid_AC[myeloid_AC$fine_celltype == i,]
  tmp <- apply(tmp, 1, function(x) colnames(tmp)[which.min(x)])
  tmp <- as.data.frame(prop.table(table(tmp)))
  if(nrow(tmp) < 4){
     library(miscTools)
     tmp <- as.matrix(tmp) ;tmp
     tmp <- insertRow(tmp,1,c('NO',0))
     tmp <- transform(tmp,Freq=as.numeric(as.character(Freq)))
  }
  tmp$cluster <- i
  data <- rbind(data,tmp)
  data$tmp <- factor(data$tmp,levels = c('NO','OC','ML','AW','AC'))
}

col <- structure(c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF','#F39B7FFF'),names = c('NO','OC','ML','AW','AC'))
ggplot(data,aes(x=cluster,y=Freq,fill=tmp))+
  geom_bar(position = 'fill',stat="identity")+
  labs(x="Cluster",y="Fraction of metastasis cells")+
  theme(panel.background=element_rect(fill='transparent', color='black',size = 2),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"),
        axis.text.x =  element_text(size = 14,angle = 45,vjust = 0.8,hjust = 0.8),
        axis.text.y = element_text(size = 14),
        axis.ticks.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 18),
        legend.text = element_text(size = 14))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = ' '))+
  scale_fill_manual(values = col,breaks = c('NO','OC','ML','AW','AC'))

#----sFigure 6C----
NKT <- subset(all_data_harmony,idents=c("T","NK"))
NKT <- FindVariableFeatures(NKT,selection.method = "vst", nfeatures = 2000) 
NKT <-ScaleData(NKT,features=VariableFeatures(NKT)) 
NKT <- RunPCA(NKT,features = VariableFeatures(NKT))
NKT <-RunHarmony(NKT,group.by.vars = 'orig.ident') 
NKT <- RunUMAP(NKT,reduction = 'harmony',dims = 1:20,min.dist = 0.5,seed= 19)
NKT <-FindNeighbors(NKT,reduction = "harmony", dims = 1:20)
NKT <- FindClusters(NKT,resolution = 0.8)
col.cluster <- colorRampPalette(brewer.pal(9, "Set1"))
col <- col.cluster(16)

new.cluster.ids=c("CD8_c1_GZMA","CD4_c1_IL7R","CD8_c2_TTN","NK_GNLY","CD4_c2_FOXP3","CD8_c3_HAVCR2",
                  "CD4_c3_CXCL13","CD8_c4_MKI67","CD8_c5_IFIT3","CD8_c6_FGFBP2","CD4_c4_IKZF2",
                  "CD8_c4_MKI67","CD8_c1_GZMA")
names(new.cluster.ids) <- levels(NKT)
NKT <- RenameIdents(NKT, new.cluster.ids)

levels(NKT) <- sort(levels(NKT))
DimPlot(NKT, reduction = "umap", label = T, pt.size = 0.5,label.size = 5,cols = col)+
  theme(legend.text = element_text(size = 18))

#----sFigure 6D----
tmp <- AverageExpression(NKT,features = unique(top20$gene)) %>% .[['RNA']]

annotation_col <- data.frame(cluster = colnames(tmp))
rownames(annotation_col) <- annotation_col$cluster
colo <- c(paletteer_d("awtools::spalette"),paletteer_d("fishualize::Taeniura_lymma"))
names(colo) <- colnames(tmp)
ann_colors = list(
  cluster = colo
)
pheatmap(tmp,show_colnames =F,show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,  
         cluster_rows = F,
         cluster_cols = F,
         scale = 'row')

#----sFigure 6E----
T_func <- list(c('CD4','CD8A','CD8B'),
               c('TCF7','SELL','LEF1','CCR7'),
               c('IL2RA','FOXP3','IKZF2'),
               c('CXCL13','CD40LG','CXCR4'),
               c('BCL6','GZMK','CCL4','IFNG','ANXA1'),
                c('PDCD1','TIGIT','LAG3','HAVCR2'),
               c('MKI67'),
                c('GZMA','EOMES','KLRG1','NKG7')
               )

T_func_ave <- AverageExpression(NKT,features = unlist(T_func)) %>% .[['RNA']] 
outlier <- boxplot(as.numeric(scale(T_func_ave)))$out
r <-  range(scale(T_func_ave)[! scale(T_func_ave )  %in% outlier])
col.heatmap <- colorRampPalette(c(rev(brewer.pal(n = 5, name ="RdYlBu"))))(1000)
pheatmap(T_func_ave,show_colnames = T,show_rownames = T,
         cluster_rows = F,
         cluster_cols = F,
         color = col.heatmap,
         border_color='white',
         breaks = seq(r[1],r[2],length.out = 1000),
         legend_breaks = c(-1,0,1),
         fontsize = 10,scale='row',angle_col = '45',
         gaps_row = c(3,7,10,13,18,22,23), gaps_col = c(4,10),
         cellwidth = 20,cellheight = 10
         )
