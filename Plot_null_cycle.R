#Run a PCA for proteomics data
#Lauren Baugh 12/12/19
#uses http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
#http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining

setwd("C:/Users/Lauren/Google Drive/Lauffenburger Lab/Proteomics/Combined_uterine/Mixed Model")
file="786_830_cycle_null_PCA_v2.csv"
data=read.table(file, header=TRUE, sep=',', stringsAsFactors = FALSE)

data = data[,-1]
data.trans = t(data)

colnames(data.trans)=data.trans[1,]
data.r=data.trans[-1,]
data.ready=data.frame(data.matrix(data.r))
#rownames(data.ready)=c(1:nrow(data.ready))

data.pca=prcomp(data.matrix(data.ready[1:ncol(data.ready)]),center=TRUE, scale.=TRUE)
summary(data.pca)

library(factoextra)
fviz_eig(data.pca)

rn=row.names(data.r)
remove.num=gsub('[[:digit:]]+', '', rn) #remove numbers from strings
colors=remove.num

# colors=c('c','o','c','c','o','c','o','o','c','c')
# colors=c('o','c','c','o','o','o','o','o','o')

fviz_pca_ind(data.pca,
             axes = c(1, 2),
             col.ind=colors,
             #col.ind = "cos2", # Color by the quality of representation
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_cos2(data.pca, choice = "var", axes = 1:2)+ coord_flip()

# Contributions of variables to PC1
fviz_contrib(data.pca, choice = "var", axes = 1, top = 10) #redline is average contribution
# Contributions of variables to PC2
fviz_contrib(data.pca, choice = "var", axes = 2, top = 10)

fviz_pca_ind(data.pca,
             axes = c(1, 2),
             col.ind=colors,
             addEllipses=TRUE, ellipse.level=0.95,
             repel = TRUE     # Avoid text overlapping
)

library("corrplot")
var=get_pca_var(data.pca)
corrplot(var$contrib, is.corr=FALSE) 


#volcano plot
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
library(EnhancedVolcano)
setwd("C:/Users/Lauren/Dropbox (MIT)/Endometriosis_scRNAseq/Mixed_Models/Proteomics/786_830")

LRT=read.table('nc_t.value.pLRT.value.csv', header=TRUE, sep=',', stringsAsFactors = FALSE)
rownames(LRT)=LRT$X
pvals=data.frame(LRT$pvalue)
rownames(pvals)=LRT$X

data.vol=data[,c(-1,-2)] 
rownames(data.vol)=data$Protein.Descriptions

control=data.vol[,c(1,3,5,6)]#control patients
osis=data.vol[,c(2,4,7,8)]#osis patients

c.avg=data.frame(rowMeans(control)) #average of log2 fold change
o.avg=data.frame(rowMeans(osis))#average of log2 fold change

log2fd=data.frame(c.avg-o.avg)
colnames(log2fd)='Log2FoldChange'

lrt.subset=merge(log2fd, pvals, by = 0,all.x=TRUE)

EnhancedVolcano(toptable=lrt.subset,
                lab = lrt.subset$Row.names,
                x = 'Log2FoldChange',
                y = 'LRT.pvalue',
                title= "Control Vs Osis",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labSize = 4.0)
