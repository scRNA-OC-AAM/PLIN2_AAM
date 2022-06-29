library(Seurat)
library(ggplot2)
library(tidyr)
library(limma)

#----Figure 3A----
load('mac.Rdata')
tmp <- data.frame()
for (i in levels(mac)) {
  tmp[1,i] <- mean(features.scores[features.scores$fine_celltype == i,'M1'])
  tmp[2,i] <- mean(features.scores[features.scores$fine_celltype == i,'M2'])
}
rownames(tmp) <- c('M1','M2')

t <- t(tmp) %>% as.data.frame()
t <- rownames_to_column(t,'cluster')
colnames(t)[c(2,3)] <- 'value'
t <- rbind(t[,c(1,2)],t[,c(1,3)])
t[['group']] <- c(rep('M1',11),rep('M2',11))
t$group <- factor(t$group,levels = c('M2','M1'))

ggplot(t,aes(cluster,group))+geom_tile(aes(fill = value),color = 'white',size = 0.5)+
  theme(panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
        axis.text = element_text(size = 14,face = 'bold'),axis.title = element_blank(),axis.text.x = element_text(angle = 45,vjust = 1.07,hjust = 1),
        legend.title = element_blank(),legend.text = element_text(size = 14,vjust = 0),legend.position = 'right',
        axis.ticks = element_blank())+
  coord_fixed(ratio=1)+
  scale_fill_gradientn(colors = brewer.pal(9,'Reds'),label = c('Min',' ',' ','Max',''))+
  guides(fill = guide_colorbar(ticks = F,barwidth = unit(3,'mm'),barheight = unit(3,'cm')))

#----Figure 3B----
M1_M2 <- read_excel('M1_M2.xlsx',col_names = T) %>% as.data.frame() %>% .[-1,]

features.scores <- matrix(data = numeric(length = 1L), nrow = 2, ncol = ncol(x = mac))
for (i in 1:2) {
  features.use <- features[[i]]
  data.use <- assay.data[features.use, , drop = FALSE]
  features.scores[i, ] <- Matrix::colMeans(x = data.use)
}

features.scores <- t(features.scores) %>% as.data.frame()
features.scores <- cbind(features.scores,as.character(mac$fine_celltype))
rownames(features.scores) <- colnames(mac)
colnames(features.scores) <- c(names(features),'fine_celltype')

tmp <- features.scores[features.scores$fine_celltype == 'mac_c4_PLIN2',]
tmp <- gather(tmp,key = feature,value = score,-fine_celltype)

ggplot(tmp,aes(x = feature, y = score, fill = feature))+
  geom_violin()+
  labs(title = 'mac_c4_PLIN2')+
  theme(panel.background = element_rect(fill = 'transparent',colour = 'black',linetype = 'solid',size = 2),
        plot.title = element_text(size = 18,hjust = 0.5),
        axis.text = element_text(size = 14),axis.title.y = element_text(size = 18),axis.title.x = element_blank(),
        legend.title = element_blank(),legend.text = element_text(size = 14))+
  scale_fill_manual(values = c('#006064FF','#96281BFF'))

#----Figure 3C----
tmp <- subset(myeloid,idents = c('mac_c4_PLIN2','mac_c6_ISG15'))
topgene <-  FindAllMarkers(tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

topgene$threshold = factor(ifelse(topgene$p_val_adj < 0.05 & abs(topgene$avg_logFC) >= 1, 
                                  ifelse(topgene$avg_logFC>= 1 ,'Up','Down'),'NS'),levels=c('Up','Down','NS'))

exp <- AverageExpression(tmp,features = topgene$gene)[['RNA']] %>% log1p()
df <- cbind(exp,as.data.frame(topgene[,c('p_val_adj','threshold','cluster')]))
df$group <- ifelse(df$threshold == 'Up' & df$cluster == 'mac_c4_PLIN2','mac_c4_PLIN2',
                   ifelse(df$threshold == 'Up' & df$cluster == 'mac_c6_ISG15','mac_c6_ISG15','NS'))
df <- filter(df,threshold != 'NS')
df <- rownames_to_column(df,'gene')
dff <- AverageExpression(tmp,features = rownames(tmp)[!rownames(tmp) %in% df$gene])[['RNA']] %>% log1p()

ggplot(df, aes(x = mac_c6_ISG15, y = mac_c4_PLIN2))+
  geom_point(aes(color = group), size = 2)+  
  geom_point(data = dff, aes(x = mac_c6_ISG15, y = mac_c4_PLIN2), alpha = 0.5, color = 'black')+
  geom_text_repel(data = df, aes(label = gene,color = group), size = 4.5, segment.color = "black", show.legend = FALSE)+
  scale_color_manual(values = c("#CB6651", "#FFBF19"))+ 
  theme(panel.background = element_rect(fill = 'transparent'),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left = element_line(color = 'black'),
        axis.title = element_text(size = 14,color="black"),
        axis.text = element_text(size = 12,color="black"),
        legend.key =element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())+
  labs(x = 'Average expression in mac_c6_ISG15', y = 'Average expression in mac_c4_PLIN2', color = '') + 
  geom_abline(intercept = 0, slope = 1, col = 'red', linetype = 'dashed', size = 1)

#----Figure 3D----
t_gbm_es <- read.table("myeloid_hall.txt", header = T, check.names = FALSE)
robj <- get(load('myeloid.Rdata'))
t_gbm_es$celltype <- 'unkown'
for (i in rownames(t_gbm_es)){
  t_gbm_es[i,]$celltype <- as.character(robj@meta.data[i,]$fine_celltype)
}

t_gbm_es <- t_gbm_es[,-c(1,2,3,6,7,9,11,12,18,19,23,25,31,40,42,43,44,49)]
colnames(t_gbm_es) <- c('EMT','Inflammatory response','KRAS activation dn','E2F transcription factors',
                        'Peroxisome','Complement system','Glycolysis and gluconeogenesis','Fatty acids metabolism',
                        'Oxidative phosphorylation','Mitotic spindle assembly','Bile acids and salts metabolism',
                        'Cholesterol homeostasis','DNA repair','p53 pathways','Protein secretion pathway',
                        'G2/M checkpoint','Apoptosis','TNFA signaling via NF-kB','Hedgehog signaling',
                        'Notch signaling','PI3K/AKT/mTOR pathway','WNT signaling','IL6 STAT3 signaling','KRAS activation up',
                        'ROS','IL2 STAT5 signaling','Angiogenesis','Interferon alpha response','Interferon gamma response',
                        'Hypoxia','TGFB1 response','mTORC1 complex signaling','celltype')

t_gbm_es <- t_gbm_es[t_gbm_es$celltype == c('mac_c4_PLIN2','mac_c6_ISG15'),]
gbm_es <- t(t_gbm_es) %>% as.data.frame()

gbm_es <- gbm_es[c(1:32),]
gbm_es_n <- apply(gbm_es, 2, as.numeric)
rownames(gbm_es_n) <- rownames(gbm_es)

group <- as.factor(t_gbm_es$celltype)
desigN <- model.matrix(~ group + 0)
colnames(desigN) <- levels(group)
rownames(desigN) <- colnames(gbm_es)

contr.matrix <- makeContrasts(mac_c6_ISG15 - mac_c4_SPP1, levels = colnames(desigN))#差异比较矩阵(谁在前谁是分子)

fiT <- lmFit(gbm_es_n, desigN)
fiT2 <- contrasts.fit(fiT, contr.matrix)
fiT3 <- eBayes(fiT2)
sceDiff <- topTable(fiT3, coef=1, number=Inf, p.value = 0.05)

sceDiff <- sceDiff[order(sceDiff$t, decreasing = TRUE),]

t.value <- sceDiff[, 't', drop = FALSE]

t.value$group <- "unknown"
t.value[which(t.value$t > 0 ),]$group <- "mac_c6_ISG15"
t.value[which(t.value$t < 0 ),]$group <- "mac_c4_PLIN2"

t.value <- t.value[order(t.value$t,decreasing = T),]
t.value$pathway <- factor(rownames(t.value), levels = rownames(t.value)[length(rownames(t.value)):1])

ggplot(t.value,aes(x = pathway, y = t, fill=group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#CB6651", "#FFBF19")) +
  coord_flip() +
  geom_hline(yintercept = 0,lwd = 1)+
  theme(panel.background = element_rect(fill = 'transparent',colour = 'transparent'),
        axis.line.x = element_line(size = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.85,0.8),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12))+
  geom_text(aes(label = pathway, y = ifelse(t >= 0, -0.5, 0.5), hjust = ifelse(t > 0, 1, 0),colour = group), size = 4.5)+
  labs(y="t value of GSVA score",  x=" ", fill = "")+
  scale_y_continuous(limits = c(-30,30),breaks = seq(-30, 30, 15))+
  scale_color_manual(values = c("#CB6651", "#F2990CFF"))+
  guides(color = FALSE)

#----Figure 3E----
load('mac.Rdata')
TAM_marker <- FindAllMarkers(mac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
C4set_m <- c(TAM_marker[TAM_marker$cluster == 'mac_c4_PLIN2',]$gene)
C6set_m <- c(TAM_marker[TAM_marker$cluster == 'mac_c6_ISG15',]$gene)

data <- TCGAexpr[C2set_m[-c(1:2)],] %>% colMeans()
data <- TCGAexpr[C3set_m,] %>% colMeans()
data <- TCGAexpr[C4set_m[-c(7,8,12,20)],] %>% colMeans()
data <- TCGAexpr[C6set_m,] %>% colMeans()

data <- data.frame(exp = data,check.rows = F)
tmp <- merge(TCGAclin, data, by.x="sample", by.y=0)
tmp$group <- ifelse(tmp$exp < quantile(tmp$exp,0.25),'low',
                    ifelse(tmp$exp > quantile(tmp$exp,0.75),'high','NA'))
tmp <- tmp[!tmp$group == 'NA',]

ggsurvplot(survfit(Surv(OS.time, OS)~group, data=tmp), conf.int=F, pval=TRUE,
           surv.median.line = "hv",ggtheme = theme_pubr(),
           legend.title = "",palette = "uchicago",,risk.table = TRUE,
           legend = c(0.8,0.75), 
           legend.labs = c("high", "low"), 
           break.x.by = 1000,
           tables.height = 0.2,font.x = 18,font.y = 18,font.legend = 14,risk.table.fontsize = 5,
           title = 'mac_c4_SPP1',font.title = 18)