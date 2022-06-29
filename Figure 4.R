library(Seurat)
library(ggplot2)
library(ggvenn)
library(ggsci)
library(scales)

#----Figure 4A----
load('mac.Rdata')
mac$fine_celltype <- factor(mac$fine_celltype,levels = sort(as.character(unique(mac$fine_celltype)),decreasing = F))

myeloid.markers <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
myeloid.markers <-  myeloid.markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05) 
modify_mye_new <-  myeloid.markers %>% group_by(cluster) %>% arrange(.,cluster,desc(avg_logFC)) %>% top_n(n = 10, wt = avg_logFC)
tmp <- modify_mye_new[modify_mye_new$cluster == 'mac_c4_PLIN2', c('gene','cluster','avg_logFC')] %>% as.data.frame()

VlnPlot(mac, features = tmp$gene, pt.size = 1,stack = TRUE,flip = T,fill.by = 'ident',cols = c("#4DAF4A","#727E76","#984EA3","#CB6651","#FF7F00","#FFBF19","olivedrab3"))+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y.right = element_text(face = 'plain'))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = 'none')+
  labs(title = paste0('TOP 10 marker of',' ','mac_c4_PLIN2'))+
  scale_x_discrete(limit = sort(as.character(unique(mac$fine_celltype))))

#----Figure 4B----
position_deg_MLvAC <- FindMarkers(mac, ident.1 = 'AC', ident.2 = 'ML',logfc.threshold = 0.25)
position_deg_MLvAC <- rownames_to_column(position_deg_MLvAC,var = 'gene')
position_deg_MLvAC$threshold = factor(ifelse(position_deg_MLvAC$p_val_adj < 0.001 & abs(position_deg_MLvAC$avg_logFC) >= 1, 
                                             ifelse(position_deg_MLvAC$avg_logFC>= 1 ,'AC','ML'),'NS'),levels=c('AC','ML','NS'))

position_deg_OCvAC <- FindMarkers(mac, ident.1 = 'AC', ident.2 = 'OC',logfc.threshold = 0.25)
position_deg_OCvAC <- rownames_to_column(position_deg_OCvAC,var = 'gene')
position_deg_OCvAC$threshold = factor(ifelse(position_deg_OCvAC$p_val_adj < 0.001 & abs(position_deg_OCvAC$avg_logFC) >= 1, 
                                             ifelse(position_deg_OCvAC$avg_logFC>= 1 ,'AC','OC'),'NS'),levels=c('AC','OC','NS'))

position_deg_AWvAC <- FindMarkers(mac, ident.1 = 'AC', ident.2 = 'AW',logfc.threshold = 0.25)
position_deg_AWvAC <- rownames_to_column(position_deg_AWvAC,var = 'gene')
position_deg_AWvAC$threshold = factor(ifelse(position_deg_AWvAC$p_val_adj < 0.001 & abs(position_deg_AWvAC$avg_logFC) >= 1, 
                                             ifelse(position_deg_AWvAC$avg_logFC>= 1 ,'AC','AW'),'NS'),levels=c('AC','AW','NS'))

position_deg_NOvAC <- FindMarkers(mac, ident.1 = 'AC', ident.2 = 'NO',logfc.threshold = 0.25)
position_deg_NOvAC <- rownames_to_column(position_deg_NOvAC,var = 'gene')
position_deg_NOvAC$threshold = factor(ifelse(position_deg_NOvAC$p_val_adj < 0.001 & abs(position_deg_NOvAC$avg_logFC) >= 1, 
                                             ifelse(position_deg_NOvAC$avg_logFC>= 1 ,'AC','NO'),'NS'),levels=c('AC','NO','NS'))
x <- list(
  NO = position_deg_NOvAC[position_deg_NOvAC$threshold == 'AC','gene'],
  OC = position_deg_OCvAC[position_deg_OCvAC$threshold == 'AC','gene'],
  TC = position_deg_MLvAC[position_deg_MLvAC$threshold == 'AC','gene'],
  PF = position_deg_AWvAC[position_deg_AWvAC$threshold == 'AC','gene']
)

ggvenn(
  x, 
  show_elements = F,
  fill_color = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"),fill_alpha = 0.8,
  stroke_size = 1, set_name_size = 8,text_color = "white",text_size = 5
)

#----Figure 4C----
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = F, scale = "width",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

df <- GetAssayData(mac,slot = 'data') %>%.[c('PLIN2'),] 
df <- as.matrix(df) %>% as.data.frame() 
df$orign <- mac$orign
names(df) <- c('PLIN2','orign')

d <- group_by(df, orign) %>%
  summarize(mean = mean(PLIN2),
            sd = sd(PLIN2))

ggplot(df, aes(orign,PLIN2, fill=orign))  +
  geom_flat_violin(position=position_nudge(x=.2)) +
  geom_jitter(aes(color=orign), width=.15) +
  geom_pointrange(aes(y=mean, ymin=mean-sd, ymax=mean+sd),
                  data=d, size=1, position=position_nudge(x=.2)) +
  coord_flip() + 
  theme(panel.background=element_rect(fill="white",colour="black",size = 1.3),
        panel.grid = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1,size = 14,color = 'black'),
        axis.text.y = element_text(size = 14,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'black'),
        legend.position = 'none')+
  scale_x_discrete(limits = c('NO','OC','ML','AW','AC'))+
  scale_fill_manual(values =  cl)+
  scale_color_manual(values = cl)+
  labs(x=NULL,
       y= 'Expression level')+
  stat_compare_means(comparisons=my_comparisons,label = 'p.signif',label.y = c(7.5,8.5,9.5,10.5))

#----Figure 4D----
Pt5 <- subset(myeloid,idents = c('OC5','AC1'))
Pt6 <- subset(myeloid,idents = c('OC6','AC2','ML2'))
Pt7 <- subset(myeloid,idents = c('OC5','AC3','ML3'))
df <- list()
for (i in c('Pt5','Pt6','Pt7')) {
  tmp <- as.data.frame(prop.table(table(tmp <- get(i)@meta.data[which(get(i)$fine_celltype == 'mac_c4_PLIN2'),'orign'])))
  df[[i]] <- tmp
}

df <- do.call('rbind',df)
df <- rownames_to_column(df,'pt')
df$pt <- substr(df$pt,1,3)
tmp$Var1 <- as.character(tmp$Var1)
tmp[which(tmp$Var1 != 'AC'),]$Var1 <- 'nAC'
compare_means(Freq~Var1, data=tmp, paired = F,method = 't.test')

ggplot(df, aes(x=Var1, y=Freq, group=pt, shape=pt)) + 
  geom_line(aes(color = pt),size=1) + 
  geom_point(aes(color = pt),size=5)+
  theme(panel.background=element_rect(fill='transparent'),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left = element_line(color = 'black'),
        axis.title = element_text(size = 14,color="black"),
        axis.text = element_text(size = 14,color="black"),
        legend.key =element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  labs(x='',y = 'Propotion of mac_c4_PLIN2(%)',shape = 'Patient',color = 'Patient')+
  scale_color_nejm()+
  scale_y_continuous(labels = percent)+
  geom_text(x = 3,y=0.7,label = 'p = 0.012',size=5.5)

#----Figure 4E----
tmp <- read.table('GSE146026_Izar_HGSOC_ascites_10x_log.tsv',header = T,row.names = 1)
mt <- t(tmp[1:7,]) %>% as.data.frame()
tmp <- tmp[-c(1:7),]
tmp <- CreateSeuratObject(counts = tmp)

tmp <- AddMetaData(object = tmp, mt)

tmp<- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
tmp<- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
tmp<- RunPCA(tmp, features = VariableFeatures(object = tmp))
tmp <- FindNeighbors(tmp, dims = 1:20)
tmp <- FindClusters(tmp, resolution = 0.2)
tmp <- RunUMAP(tmp,reduction = "pca", dims = 1:20)
tmp <- RunTSNE(tmp, verbose = T)
marker <- FindAllMarkers(tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top <-  dplyr::group_by(marker,cluster) %>% 
  top_n(a,n = 50, wt = avg_logFC)

new.cluster.ids=c("Macrophage","Macrophage","Fibroblast","Cancer cell","Fibroblast","B","T","Macrophage","Cancer cell","Cancer cell","Macrophage","DC")
names(new.cluster.ids) <- levels(tmp)
tmp <- RenameIdents(tmp, new.cluster.ids)
tmp[['celltype']] <- Idents(tmp)

tmp <- subset(tmp,idents = c('Macrophage','T','B','DC'))
phe <- data.frame(cell = rownames(tmp@meta.data),cluster = tmp@meta.data$celltype)
tsne_pos <- Embeddings(tmp,'tsne') 
dat <- cbind(tsne_pos ,phe)
ggplot(dat,aes(x = tSNE_1,y = tSNE_2,color = cluster))+
  geom_point(size = 0.01)+
  scale_color_manual(values = c('#C71000FF','#008EA0FF','#1A5354FF','#8A4198FF','#FF95A8FF','#FF6348FF'))+
  theme(panel.grid = element_blank()) + 
  theme(panel.border = element_blank(),panel.background = element_blank()) +  
  theme(axis.line = element_line(size=1, colour = "black"),axis.text = element_text(size = 14),axis.title = element_text(size = 14))+
  theme(legend.key = element_blank())+ 
  theme(legend.text = element_text(size = 15),legend.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 5),ncol = 1))

#----Figure 4F----
plin2 <- GetAssayData(tmp,slot = 'data') %>% .['PLIN2',] %>% as.data.frame()

umap_pos <- Embeddings(tmp,'tsne') 
a <- cbind(plin2,umap_pos) %>% as.data.frame()
ggplot(a,aes(x = tSNE_1,y = tSNE_2))+
  geom_point(aes(color = .),size = 0.85)+
  guides(color = guide_colorbar(title = 'Average expression',label = F))+
  scale_color_gradientn(colors = rev(brewer.pal(n=11, name = 'Spectral')))+
  theme(panel.background = element_rect(fill = 'transparent',colour = 'black',linetype = 'solid',size = 2),
        axis.ticks = element_blank(),axis.title = element_blank(),axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 18),plot.subtitle = element_text(hjust = 0.5,size = 14),
        legend.title = element_text(size = 14),panel.grid = element_blank())

VlnPlot(tmp,features = 'PLIN2',
        cols = c('#C71000FF','#008EA0FF','#1A5354FF','#8A4198FF','#FF95A8FF','#FF6348FF'),
        pt.size =  0) + scale_x_discrete(limits = rev(levels(tmp))) + labs(title = '',x = '') + guides(fill = FALSE) +coord_flip()
compare_means(PLIN2~ident,data = p$data,ref.group = 'Macrophage')