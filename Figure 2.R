library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(GSVA)
library(readr)
library(stringr)
library(msigdbr)
library(pheatmap)

#----Figure 2A----
myeloid <- subset(all_data_harmony,idents = "Myeloid")
myeloid <- NormalizeData(myeloid, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features=VariableFeatures(myeloid)) %>%
  RunHarmony("orig.ident", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.6)

myeloid.markers <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_myeloid <- group_by(myeloid.markers,cluster) %>% top_n(a,n = 10, wt = avg_logFC)

new.cluster.ids=c("mac_c1_C1QC","mono_c1_CD14","mac_c2_NLRP3","DC_c1_cDC2","mac_c3_FN1","mac_c5_SPP1",
                  "mac_c6_MALAT1","mac_c7_ISG15","mac_c3_FN1","monolike","mono_c2_CD16","myeloid_proliferation",
                  "DC_c2_cDC1",'mac_c5_SPP1','mono_c3_CD14CD16','DC_c3_pDC','DC_c4_cDC3')

names(new.cluster.ids) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, new.cluster.ids)

levels(myeloid) <- sort(levels(myeloid))
DimPlot(myeloid, reduction = "umap", label = F, pt.size = 0.5,cols = col)+ 
  theme(legend.text = element_text(size = 18))

#----Figure 2B----
mac <- GetAssayData(myeloid,slot = 'data') %>% .[c('C1QC','CD163','APOE'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
dc <- GetAssayData(myeloid,slot = 'data') %>% .[c('CLEC9A','CD1C','LAMP3','IL3RA'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
mo <- GetAssayData(myeloid,slot = 'data') %>% .[c('VCAN','S100A9','S100A12'),] %>% as.matrix %>% t() %>% apply(., 1, mean) %>% as.data.frame()
umap_pos <- Embeddings(myeloid,'umap') 
tmp <- cbind(dc,umap_pos) %>% as.data.frame()
p2 <- ggplot(tmp,aes(x = UMAP_1,y = UMAP_2))+
  geom_point(aes(color = .),size = 0.85)+
  labs(title = 'DC',subtitle = '(CLEC9A,CD1C,LAMP3,IL3RA)')+
  guides(color = F)+
  scale_color_gradientn(colors = rev(brewer.pal(n=11, name = 'Spectral')))+
  theme(panel.background = element_rect(fill = 'transparent',colour = 'black',linetype = 'solid',size = 2),
        axis.ticks = element_blank(),axis.title = element_blank(),axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,face = 'bold',size = 18),plot.subtitle = element_text(hjust = 0.5,size = 14),
        legend.title = element_text(size = 14))

p1 + p2 + p3 + plot_layout(guides='collect',heights = c(1,1),nrow = 1)

#----Figure 2C----
tmp <- table(myeloid$orign,Idents(myeloid))
prop_all_data <- as.data.frame(prop.table(tmp, margin = 1))
prop_all_data$Var1 <- factor(prop_all_data$Var1,levels = rev(c(paste0('NO',1:3),paste0('OC',1:7),paste0('ML',1:3),paste0('AW',1:3),paste0('AC',1:3))))

ggplot(prop_all_data,aes(x=prop_all_data[,1],y=prop_all_data[,3],fill=prop_all_data[,2]))+
  geom_bar(position = 'fill',stat="identity")+
  labs(x=" ",y="Percentage")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"),
        axis.text.x =  element_text(size = 14,angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(face = 'bold',size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))+
  scale_y_continuous(expand=c(0.001,0.001))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = ' '))+
  scale_x_discrete(limits = c('NO','OC','TC','PF','AC'))+
  scale_fill_manual(values = coll)

#----Figure 2D----
DC_mono <- list(c('CD1C','FCER1A','CLEC10A'),
                c('CLEC9A','FLT3','IDO1'),
                c('LILRA4','GZMB','IL3RA'),
                c('LAMP3','CCR7','FSCN1'),
                c('C1QC','C1QA','APOE'),
                c('NLRP3','EREG','IL1B'),
                c('FN1','MARCO','CFD'),
                c('SPP1','PLIN2','ERO1A'),
                c('MALAT1','NEAT1','JUN'),
                c('ISG15','CXCL10','CXCL11'),
                c('MKI67','TOP2A','H2AFZ'),
                c('FCN1','S100A8','S100A9'),
                c('FCGR3A','LST1','LILRB2'),
                c('VAMP8','ELOB','ALOX5AP'))

cols <- c('dodgerblue4','#50b490','#fff337','#E64E00','firebrick')
pal<-colorRampPalette(cols)
DotPlot(myeloid,features = unlist(DC_mono),col.min = -2,col.max = 2)+
  scale_color_gradientn(colors = pal(20))+
  scale_size_area(breaks=c(0,20,40,60,80))+
  theme(panel.grid = element_line(colour = 'grey',size = 0.1,linetype = 'dashed'),
        panel.background = element_rect(fill = 'white', colour = 'black',size = 1,linetype = 'solid'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 0.9))+
  scale_y_discrete(limits = rev(levels(myeloid)))

#----Figure 2E----
robj <- get(load('myeloid.Rdata'))
data <- t(as.matrix(robj@assays$RNA@data))

cellname<-rownames(data)
genename<-colnames(data)

ndata <- data[apply(data,1,function(x) !all(x==0)),]
tdata<-t(ndata)
tdata <- tdata[apply(tdata,1,function(x) !all(x==0)),]

gmt2list <- function(gmtfile){
  sets <- as.list(read_lines(gmtfile))
    for(i in 1:length(sets)){
      tmp = str_split(sets[[i]], '\t')
      n = length(tmp[[1]])
      names(sets)[i] = tmp[[1]][1]
      sets[[i]] = tmp[[1]][3:n]
      rm(tmp, n)
    }
    return(sets)
  }
s.sets = gmt2list("h.all.v7.1.symbols.gmt")

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

tmp <- split(t_gbm_es,t_gbm_es$celltype)

ls <- lapply(tmp, function(x){colMeans(x[,-51])})
df <- do.call(rbind,ls)

colnames(df) <- c('Myc targets v1','Myc targets v2','Estrogen response early','epithelial-mesenchymal transition',
                  'Inflammatory response','Estrogen response late','Androgens response','KRAS activation dn',
                  'UV response dn','E2f targets','Apical junction','Coagulation','Peroxisome',
                  'Complement','Glycolysis','Fatty acids metabolism','Oxidative phosphorylation',
                  'Xenobiotics metabolism','Apical surface','Mitotic spindle','Bile acids metabolism',
                  'Cholesterol homeostasis','Myogenesis','DNA repair','Heme metabolism','p53 pathways',
                  'Protein secretion','G2/M checkpoint','Apoptosis','TNFa signaling via NF-kB',
                  'pancreatic beta cells','Hedgehog signaling','Notch signaling','PI3K/AKT/mTOR pathway',
                  'WNT beta catenin signaling','IL6 STAT3 signaling','KRAS activation up','ROS pathway',
                  'IL2 STAT5 signaling','Adipogenesis','Angiogenesis','Spermatogenesis',
                  'Allograft rejection','Unfolded protein response','Interferon alpha response',
                  'Interferon gamma response','Hypoxia','TGFB1 signaling','UV response dn','mTORC1 complex activation')

df <- df[,c('Allograft rejection','Coagulation','Complement','Interferon alpha response',
            'Interferon gamma response','IL6 STAT3 signaling','Inflammatory response',
            'Bile acids metabolism','Cholesterol homeostasis','Fatty acids metabolism',
            'Glycolysis','Heme metabolism','Oxidative phosphorylation','Xenobiotics metabolism',
            'Androgens response','Estrogen response early','Estrogen response late','IL2 STAT5 signaling',
            'KRAS activation dn','KRAS activation up','mTORC1 complex activation',
            'Notch signaling','PI3K/AKT/mTOR pathway','Hedgehog signaling','TGFB1 signaling',
            'TNFa signaling via NF-kB','WNT beta catenin signaling','Apoptosis','Hypoxia','Protein secretion',
            'Unfolded protein response','ROS pathway','E2f targets','G2/M checkpoint',
            'Myc targets v1','Myc targets v2','p53 pathways','Mitotic spindle')]

annotaion_col <- data.frame(class = rep(factor(rep(c('Immune','Metabolism','Signaling','Proliferation'),c(7,7,18,6)))))
rownames(annotaion_col) <- colnames(df)
bk <- seq(-2,2,length.out = 100)

pheatmap(df,scale = 'column',cluster_rows = F,cluster_cols = F,
         annotation_col = annotaion_col,
         color = colorRampPalette(rev(c("#9E0142","#F46D43","#FFFFBF","#66C2A5","#5E4FA2")))(100),
         breaks = bk,cellwidth = 15,cellheight = 15,border_color = 'white',
         annotation_colors = list(class = c(Immune = '#2A9D8E',Metabolism = '#4B74B2',Signaling = '#E66F51',Proliferation = '#DC9E13')))
#----Figure 2F (python & R)----
#python
import scvelo as scv
import scanpy as sc
import pandas as pd
scv.settings.verbosity = 3 
scv.settings.set_figure_params(dpi=300)  

adata = sc.read("loom_file.h5ad")

n = ["mac_c1_C1QC","mac_c2_NLRP3","mac_c3_FN1","mac_c4_SPP1" , "mac_c5_MALAT1","mac_c6_ISG15","mac_c7_proliferation", 
   "mono_c1_CD14","mono_c2_CD16","mono_c3_CD14CD16","monolike"];type(n)

myeloid_cell = adata.obs['fine_celltype'].isin(n)
myeloid_cell_adata = adata[myeloid_cell,:]

#R
cells <- names(myeloid$fine_celltype[myeloid$fine_celltype %in% levels(myeloid)[5:15]])
emb <- Embeddings(myeloid,reduction = 'umap')[cells,]
write.csv(emb,file = 'mac_mono.csv')

#python
myeloid_UMAP = pd.read_csv('mac_mono.csv')
myeloid_UMAP = myeloid_UMAP.rename(columns = {"Unnamed: 0":'Cell ID'})
myeloid_UMAP = myeloid_UMAP.set_index("Cell ID")
myeloid_UMAP = myeloid_UMAP.loc[myeloid_UMAP.index.intersection(myeloid_cell_adata.obs.index.values)] 
myeloid_cell_adata.obsm['New_umap']=myeloid_UMAP.values

scv.pl.proportions(myeloid_cell_adata,groupby='fine_celltype',show=False,fontsize=7)
scv.pp.filter_and_normalize(myeloid_cell_adata,min_shared_counts=30, n_top_genes=2000)
scv.pp.neighbors(myeloid_cell_adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(myeloid_cell_adata)
scv.tl.velocity(myeloid_cell_adata, mode='dynamical')
scv.tl.velocity_graph(myeloid_cell_adata)

col = pd.read_csv('col.csv')
del col['Unnamed: 0']
import numpy as np
col = np.array(col)
col_list=col.tolist()
from tkinter import _flatten
color = list(_flatten(col_list))
color=color[:15]

myeloid_cell_adata.obs['fine_celltype']=myeloid_cell_adata.obs['fine_celltype'].astype("category")
clu = pd.read_csv('cluster.csv')
clu= np.array(clu)
clu_list=clu.tolist()
clu_list= list(_flatten(clu_list))[:15]
myeloid_cell_adata.obs['fine_celltype']=myeloid_cell_adata.obs['fine_celltype'].cat.set_categories(clu_list)

scv.pl.velocity_embedding_stream(myeloid_cell_adata,color=['fine_celltype'],basis='New_umap',linewidth=2,dpi = 300,palette =color,legend_loc='none',arrow_size=2,alpha=1,size = 25,
                                 frameon=False,figsize=(12,7),xlabel='UMAP1',ylabel='UMAP2',title='',show=False,
                                 save = 'RNA velocity,pdf')
#----Figure 2G----
mye_meta <- myeloid@meta.data
M <- table(mye_meta[,c('orign','fine_celltype')])
(Xsq <- chisq.test(M))
mye_meta_ob <- as.data.frame(table(mye_meta[,c('orign','fine_celltype')]))
mye_meta_expect <- as.data.frame(Xsq$expected)
mye_meta_ob$expect <- 0
for (i in unique(mye_meta_ob$orign)) {
  for (j in unique(mye_meta_ob$fine_celltype)) {
    mye_meta_ob[mye_meta_ob$orign == i&mye_meta_ob$fine_celltype == j,]$expect = mye_meta_expect[i,j]
  }
}
mye_meta_ob$OE <- mye_meta_ob$Freq/mye_meta_ob$expect

result_df = data.frame() 
for (i in unique(mye_meta_ob$orign)) {
  for (j in unique(mye_meta_ob$fine_celltype)) {
    result_df[i,j]=mye_meta_ob[mye_meta_ob$orign == i&mye_meta_ob$fine_celltype == j,]$OE
  }
}

outlier <- boxplot(as.numeric(as.matrix(result_df)))$out
r <- range(as.numeric(as.matrix(result_df))[! as.numeric(as.matrix(result_df))  %in% outlier])
result_df <- round(as.matrix(result_df),2)

pheatmap(t(result_df),breaks = seq(r[1],r[2],length.out = 100),border_color = 'white',cluster_rows = F,cluster_cols = F,display_numbers = TRUE,number_format = "%.2f")