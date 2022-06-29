library(Seurat)
library(harmony)
library(dplyr)
library(paletteer)
library(ggplot2)

#----Figure 1B----
clinical.df <- data.frame(
  patient = paste("P",seq(1:13),sep = ""),
  pathology = c(rep('HGSOC',7),rep('hysteromyoma',6)),
  tumor = c(rep("YES",7),rep('NO',6)),
  metastasis = c(rep('NO',3),'YES','NO',rep('YES',2),rep('NO',6)),
  ascites = c(rep('NO',4),rep('YES',3),rep('NO',6)),
  `normal ovary` = c(rep('NO',10),rep('YES',3)),
  `Abdominal flushing fluid` = c(rep('NO',7),rep('YES',3),rep('NO',3)),check.names = F
)
clinical.df2 <- melt(clinical.df,id.vars = "patient")

cols=c(
  "HGSOC"="#66C2A5","hysteromyoma"="#FC8D62",
  "YES"="black","NO"="lightgrey"
)

clinical.df2$patient=factor(clinical.df2$patient,levels = paste("P",seq(1:13),sep = ""))
clinical.df2$variable=factor(clinical.df2$variable,levels = c("pathology","tumor","metastasis","ascites","normal ovary","Abdominal flushing fluid"))

clinical.df2 %>% ggplot(aes(x=patient,y=variable))+
  geom_tile(aes(fill=value),color="white",size=1)+ 
  scale_x_discrete("",expand = c(0,0))+ 
  scale_y_discrete("",expand = c(0,0))+
  scale_fill_manual(values = cols,breaks = c('HGSOC','hysteromyoma','YES','NO'))+ 
  theme(
    axis.text.x.bottom = element_text(size=14,colour = 'black'),axis.text.y.left = element_text(size = 18,colour = 'black'), 
    axis.ticks = element_blank(), 
    legend.title = element_blank(),
    legend.key = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,color="black")
  )+
  guides(fill = guide_legend(override.aes = list(size = 8)))

#----Figure 1C----
for (i in unique(all_data_harmony$orig.ident)) {
  ls[[i]] <- subset(all_data_harmony,orig.ident==i) 
  ls[[i]] <- GetAssayData(ls[[i]],slot='data')  
  ls[[i]] <- apply(ls[[i]],1,mean)  
  ls[[i]] <- t(as.data.frame(ls[[i]])) 
  rownames(ls[[i]]) <- i
}
sc_bulk <- do.call('rbind',ls)
sc_bulk <- as.data.frame(sc_bulk)
name <- c(rep('OC',7),rep('NO',3),rep('AW',2),rep('AC',3),rep('ML',3))
sc_bulk$orig.ident <- name
#pca
res.pca <- prcomp(sc_bulk[,-33539])
par(omi = c(0,0,0,1.5),pin = c(4,4))
plot(res.pca$x,cex = 2,main = "PCA analysis", 
     col = c(rep('#4DBBD5FF',7),rep('#E64B35FF',3),rep('#3C5488FF',3),rep('#F39B7FFF',3),rep('#00A087FF',3)),pch = 19)
abline(h=0,v=0,lty=2,col="gray")
legend('topright', legend=c('NO','OC','ML','AW','AC'), fill=c(as.character(paletteer_d("ggsci::nrc_npg"))[1:5]), border ='white',box.col='white',cex = 1.2,inset = c(-0.23,0),xpd = T)

#----Figure 1D----
all_data_harmony <- NormalizeData(all_data, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features =VariableFeatures(all_data_harmony) ) %>% 
  RunPCA(features = VariableFeatures(all_data_harmony)) %>% 
  RunHarmony(group.by.vars = "orig.ident", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) 

DimPlot(all_data_harmony, reduction = "umap", pt.size = .1,label = T)
all_data.marker_harmony <- FindAllMarkers(all_data_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- group_by(all_data.marker_harmony,cluster) %>% top_n(a,n = 10, wt = avg_logFC)
new.cluster.ids=c("T","T","Myeloid","T","Myeloid","NK","T","Myeloid","T","Myeloid","Myeloid","NK","B","T","Myeloid","Myeloid","Epithelial","plasma B","Myeloid","Fibroblast")
names(new.cluster.ids) <- levels(all_data_harmony)
all_data_harmony <- RenameIdents(all_data_harmony, new.cluster.ids)
all_data_harmony$celltype <- levels(all_data_harmony)

phe <- data.frame(cell=rownames(all_data_harmony@meta.data),cluster = all_data_harmony@meta.data$celltype)
ggplot(dat,aes(x=UMAP_1,y=UMAP_2,color=cluster))+
  geom_point(size=0.1)+
  scale_color_paletteer_d("awtools::mpalette")+
  theme(panel.grid =element_blank()) +  
  theme(panel.border = element_blank(),panel.background = element_blank()) +  
  theme(axis.line = element_line(size=1, colour = "black"),axis.text = element_text(size = 14),axis.title = element_text(size = 14))+
  theme(legend.key = element_blank())+ 
  theme(legend.text = element_text(size = 15),legend.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))

#----Figure 1E----
m <- GetAssayData(all_data_harmony,slot = 'data') %>% .[c('CD68','CD1E','LAMP3','IL3RA'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
TT <- GetAssayData(all_data_harmony,slot = 'data') %>% .[c('CD2','CD3D','CD3E','CD3G','NKG7'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
B <- GetAssayData(all_data_harmony,slot = 'data') %>% .[c('MS4A1','CD79A','MZB1','IGHA1'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
epi <- GetAssayData(all_data_harmony,slot = 'data') %>% .[c('EPCAM','KRT19','KRT18','KRT14'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
fri <- GetAssayData(all_data_harmony,slot = 'data') %>% .[c('COL1A1','COL3A1','DCN','FAP'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()

umap_pos <- Embeddings(all_data_harmony,'umap') 
tmp <- cbind(fri,umap_pos) %>% as.data.frame()
p5 <- ggplot(tmp,aes(x = UMAP_1,y = UMAP_2))+
  geom_point(aes(color = .),size = 0.85)+
  labs(title = 'Fibroblast',subtitle = '(COL1A1,COL3A1,DCN,FAP)')+
  guides(color = guide_colorbar(title = 'Average expression',label = F))+
  scale_color_gradientn(colors = rev(brewer.pal(n=11, name = 'Spectral')))+
  theme(panel.background = element_rect(fill = 'transparent',colour = 'black',linetype = 'solid',size = 2),
        axis.ticks = element_blank(),axis.title = element_blank(),axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 18),plot.subtitle = element_text(hjust = 0.5,size = 14),
        legend.title = element_text(size = 14))
  
#----Figure 1F----
tmp <- table(all_data_harmony$orig.ident,Idents(all_data_harmony))
all_data_harmony@meta.data$orign <- sub('[0-9]$', '', all_data_harmony@meta.data$orig.ident)
tmp <- table(all_data_harmony$orign,Idents(all_data_harmony))

prop_all_data <- as.data.frame(prop.table(tmp, margin = 1))
prop_all_data$Var1 <- factor(prop_all_data$Var1,levels = rev(c(paste0('NO',1:3),paste0('OC',1:7),paste0('ML',1:3),paste0('AW',1:3),paste0('AC',1:3))))
mycolors <- c(brewer.pal(12,'Paired'),brewer.pal(6,'Dark2'))
plotCol(mycolors)
ggplot(prop_all_data,aes(x=prop_all_data[,1],y=prop_all_data[,3],fill=prop_all_data[,2]))+
  geom_bar(position = 'fill',stat="identity")+
  labs(x=" ",y=" ")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"),
        axis.text.x =  element_text(size = 14,angle = 90,vjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(face = 'bold',size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = 'bottom')+
  scale_y_continuous(expand=c(0,0))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = ' '))+
  scale_fill_manual(values = mycolors)+
  guides(fill = guide_legend(lable.position = "bottom",direction = "horizontal",title = 'Main Celltype',nrow = 2))+
  coord_flip()
         
(p1 + p2 + p3) /(p4 + p5) + plot_layout(guides='collect')

#----Figure 1G----
df <- as.data.frame(prop.table(table(all_data_harmony$orig.ident,all_data_harmony$celltype),margin = 1)) %>% .[-c(96:133),]

df <- mutate(df,position = substr(df$Var1,1,2))
df$position <- factor(df$position,levels = c('NO','OC','ML','AW','AC'))

df.summary <- df %>%
  group_by(position,Var2) %>%
  summarise(
    sd = sd(Freq, na.rm = TRUE),
    m = mean(Freq)
  ) %>% as.data.frame()
df.summary$position <- factor(df.summary$position,levels = c('NO','OC','ML','AW','AC'))

ggplot(data = df.summary,aes(Var2,m,fill = position))+
  geom_bar(stat = 'identity',width = 0.75,position = position_dodge(0.75))+
  geom_errorbar(data = df.summary,aes(ymin = m,ymax = m+sd),width = 0.2,position = position_dodge(0.75))+
  geom_jitter(data = df,aes(Var2,Freq),color = 'black',size = 2,position = position_jitterdodge(jitter.width = 0.05))+
  labs(y = 'Propotion of cell type',x = '')+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'transparent',colour = 'black',linetype = 'solid',size = 2),
        axis.title.y = element_text(size = 18),axis.text = element_text(size = 14,color = 'black'), 
        legend.title = element_blank(),legend.key.size = unit(5,'mm'),legend.text = element_text(size = 12))+
  scale_fill_paletteer_d("ggsci::nrc_npg")
compare_means(Freq~position, data=df[20:38,], paired = F,method = 't.test')


