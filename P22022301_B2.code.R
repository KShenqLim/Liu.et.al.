## Seurat 
library(dplyr)
library(Seurat)
library(harmony)
sample <- c('PT_84','PT_86','PT_88','RT_83','RT_85')
group <- c('PT','PT','PT','RT','RT')
#exp_list was the path of expression matrix corresponding to the sample
for(i in 1:length(exp_list)){
	if(i == 1){
		input <- Read10X(data.dir = exp_list[i])
		data <- CreateSeuratObject(counts = input, project = sample[i], min.cells = 5)
		data <- RenameCells(data,add.cell.id = sample[i])
		data@meta.data$sample <- sample[i]
		data@meta.data$group <- group[i]
		data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
		ngenemax <- unname(quantile(data$nFeature_RNA,0.95))
		umimax <- unname(quantile(data$nCount_RNA,0.95))
		data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax & percent.mt < 50,nCount_RNA > 0 &  nCount_RNA  < umimax)
		data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
	}else{
		input <- Read10X(data.dir = exp_list[i])
		input <- Read10X(data.dir = exp_list[i])
        tmp <- CreateSeuratObject(counts = input, project = sample[i], min.cells = 5)
        tmp <- RenameCells(tmp,add.cell.id = sample[i])
        tmp@meta.data$sample <- sample[i]
        tmp@meta.data$group <- group[i]
        tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
        ngenemax <- unname(quantile(tmp$nFeature_RNA,0.95))
        umimax <- unname(quantile(tmp$nCount_RNA,0.95))
        tmp <- subset(tmp, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax & percent.mt < 50,nCount_RNA > 0 &  nCount_RNA  < umimax)
        tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
		data <- merge(data,tmp)
	}
	data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(data)
	data <- ScaleData(data, features = all.genes)
	data <- RunPCA(data, features = VariableFeatures(object = data))
	DimPlot(data, reduction = "pca")
	data <- JackStraw(data, num.replicate = 100)
	data <- ScoreJackStraw(data, dims = 1:20)
	pdf(paste(outdir,'/',prefix,'.PCJackStrawPlot.pdf',sep=''))
	JackStrawPlot(data, dims = 1:20)
	dev.off()
	plot2<-ElbowPlot(data)
	CombinePlots(plots = list(plot1, plot2))
	data <- RunHarmony(PRO,"sample", plot_convergence = FALSE)
	data <- FindNeighbors(data, reduction = "harmony",dims = 1:20)
	data <- FindClusters(data, resolution = 1.2)
	data <- RunUMAP(data, reduction = "harmony",dims = 1:20)
	data <- RunTSNE(data, reduction = "harmony",dims = 1:20)
	pdf(paste(outdir,'/',prefix,'.umap.cluster.pdf',sep=''))
	DimPlot(data, reduction = "umap",label = F)
	dev.off()
	pdf(paste(outdir,'/',prefix,'.tsne.cluster.pdf',sep=''))
	DimPlot(data, reduction = "tsne",label = F)
	dev.off()
	
	for (l in levels (data)) {
	    cluster1.markers <- FindMarkers(data, ident.1 = l, min.pct = 0.01,logfc.threshold = 0.01)
	    cluster1.markers <- data.frame (gene = rownames (cluster1.markers))
	    write.table (cluster1.markers,file.path(outdir,paste0(prefix,".cluster.",l,".diffgenes.xls")),row.names = F,col.names = T, sep='\t',quote = F)
	}
}

## annotation
new.cluster.ids <- celltypes
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
pdf(paste(outdir,'/',prefix,'.umap.label.pdf',sep=''))
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

## visualization
FeaturePlot(data, features = features)
top10 <- data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(data, features = top10$gene) + NoLegend()
DotPlot(data, features = features) 



## barcodes ratio stat
library(dplyr)
#
outdir = '/path/to/out/directory';dir.create(outdir)
prefix = 'xxx_prefix'

expm = readRDS('expm.rds')
# this rds contains data.frame object which shows sample, group, expM(single cell matrix), barcodefile(adding barcode patth)
#     reportname sample group expM    barcodefile
# 1  Renca84-pre  PT_84    PT Renca84-pre_auto/   Renca84-pre_auto.xls
# 2  Renca86-pre  PT_86    PT Renca86-pre_auto/   Renca86-pre_auto.xls
# 3  Renca88-pre  PT_88    PT Renca88-pre_auto/   Renca88-pre_auto.xls
# 4 Renca83-post  RT_83    RT Renca83-post_auto/  Renca83-post_auto.xls
# 5 Renca85-post  RT_85    RT Renca85-post-CUS_auto/  Renca85-post_auto.xls

b = readRDS('b.rds')
# this rds contains list object, which shows sc_barcode adding_barcode corresponding information for each sample:
#                                              sc_barcode adding_barcode
# AACACACAGAACGCTAGTAACAGGAAC AACACACAGAACGCTAGTAACAGGAAC TTTCCTTGAAATGG
# AACACACAGAACGCTAGTGAGAGGAGT AACACACAGAACGCTAGTGAGAGGAGT CAAGTTAAGGATAA
# AACACACAGAACGCTAGTGAGATAGGA AACACACAGAACGCTAGTGAGATAGGA GAATTATTCGTAAG

gps = unique(expm$group)
sps = unique(expm$sample)
PT_sps = expm$sample[expm$group == 'PT']
RT_sps = expm$sample[expm$group == 'RT']

venn_list = lapply(gps, function(i){
    sps = expm[expm$group == i,'sample']
    barcode_group = lapply(sps, function(x) {
        adding_barcode = b[[x]]$adding_barcode %>% unique()
        out = adding_barcode[adding_barcode != "norecord"]
        return(out)
    }) %>% unlist() %>% unique()

    return(barcode_group)
})
names(venn_list) = gps

intersect_barcode = intersect(venn_list[[gps[1]]], venn_list[[2]])
saveRDS(intersect_barcode, file.path(outdir,'intersect_barcode.seq.rds'))

### all intersect_barcode stat
barcode_num = sapply(sps, function(x){
    b1 = b[[x]]
    b1_record = b1[b1$adding_barcode != 'norecord',]
    b1_intersect_barcode = b1_record[b1_record$adding_barcode %in% intersect_barcode,]
    
    
    all_CellNumbers = nrow(b1)
    withAddingBarcode_CellNumbers = nrow(b1_record)
    intersectAddingBarcode_CellNumbers = nrow(b1_intersect_barcode)
    intersectAddingBarcode_CellFrequency = intersectAddingBarcode_CellNumbers/withAddingBarcode_CellNumbers
    
    all_intersectAddingBarcode_type = length(unique(intersect_barcode))
    intersectAddingBarcode_type = length(unique(b1_intersect_barcode$adding_barcode))
    intersectAddingBarcode_typeFrequency = intersectAddingBarcode_type/all_intersectAddingBarcode_type
    
    
    out = c(all_CellNumbers,
            withAddingBarcode_CellNumbers,
            intersectAddingBarcode_CellNumbers,
            intersectAddingBarcode_CellFrequency,
            all_intersectAddingBarcode_type,
            intersectAddingBarcode_type,
            intersectAddingBarcode_typeFrequency)
    out
}) %>% as.data.frame()

rownames(barcode_num) = c('all_CellNumbers',
                          'withAddingBarcode_CellNumbers',
                          'intersectAddingBarcode_CellNumbers',
                          'intersectAddingBarcode_CellFrequency',
                          'all_intersectAddingBarcode_type',
                          'intersectAddingBarcode_type',
                          'intersectAddingBarcode_typeFrequency')

write.table(data.frame(stat= rownames(barcode_num),barcode_num),
            file.path(outdir,paste0(prefix,'.intersect_barcode.stat.xls')),
            row.names = F,
            col.names = T,
            sep = '\t',
            quote = F)
saveRDS(barcode_num, file.path(outdir,paste0(prefix,'.intersect_barcode.stat.rds')))

### each intersect_barcode rds
intersect_barcode_rds = lapply(sps, function(x){
    b1 = b[[x]]
    b1_record = b1[b1$adding_barcode != 'norecord',]
    b1_record$intersect = sapply(b1_record$adding_barcode,function(y) {ifelse(y %in% intersect_barcode, 'yes', 'no')})
    return(b1_record)
})
names(intersect_barcode_rds)= sps
saveRDS(intersect_barcode_rds, file.path(outdir,'intersect_barcode.list.rds'))

### each intersect_barcode stat
barcode_num = readRDS(file.path(outdir,paste0(prefix,'.intersect_barcode.stat.rds')))
freq_table = sapply(sps, function(x){
    b1 = intersect_barcode_rds[[x]]
    b1_inter = b1[b1$intersect == 'yes',]
    tb = table(b1_inter$adding_barcode) %>% as.matrix %>% as.data.frame()
    vec = rep(0,length(intersect_barcode))
    names(vec) = intersect_barcode
    vec[rownames(tb)] = tb$V1
    return(vec)
}) %>% as.data.frame()

freq_table$PT_Number = sapply(rownames(freq_table), function(x){ sum(freq_table[x,PT_sps])} )
freq_table$RT_Number = sapply(rownames(freq_table), function(x){ sum(freq_table[x,RT_sps])} )

freq_table$PT_CellNumber = sum(barcode_num['withAddingBarcode_CellNumbers',PT_sps])
freq_table$RT_CellNumber = sum(barcode_num['withAddingBarcode_CellNumbers',RT_sps])

freq_table$PT_CellNumber_intersect = sum(barcode_num['intersectAddingBarcode_CellNumbers',PT_sps])
freq_table$RT_CellNumber_intersect = sum(barcode_num['intersectAddingBarcode_CellNumbers',RT_sps])

freq_table$PT_Frequency = sapply(freq_table$PT_Number, function(x){ x/freq_table[x,'PT_CellNumber_intersect']} )
freq_table$RT_Frequency = sapply(freq_table$RT_Number, function(x){ x/freq_table[x,'RT_CellNumber_intersect']} )

freq_table$PT_Frequency_x10000 = freq_table$PT_Frequency*10000
freq_table$RT_Frequency_x10000 = freq_table$RT_Frequency*10000


freq_table$RT_Frequency_div_PT_Frequency = freq_table$RT_Frequency_x10000/freq_table$PT_Frequency_x10000
freq_table = freq_table[order(freq_table$RT_Frequency_div_PT_Frequency,decreasing=T),]

df_out = data.frame(adding_barcode = rownames(freq_table), freq_table)
saveRDS(df_out, file.path(outdir,paste0(prefix,'.intersect_barcode.stat.byBarcode.rds')))
write.table(df_out,
            file.path(outdir,paste0(prefix,'.intersect_barcode.stat.byBarcode.xls')),
            row.names = F,
            col.names = T,
            sep = '\t',
            quote = F)


write.table(df_out[df_out$RT_Frequency_div_PT_Frequency>1,],
            file.path(outdir,paste0(prefix,'.intersect_barcode.stat.byBarcode.greater_than_1.xls')),
            row.names = F,
            col.names = T,
            sep = '\t',
            quote = F)

write.table(df_out[df_out$RT_Frequency_div_PT_Frequency<=1,],
            file.path(outdir,paste0(prefix,'.intersect_barcode.stat.byBarcode.less_equal_than_1.xls')),
            row.names = F,
            col.names = T,
            sep = '\t',
            quote = F)


## mapping and plot
#!/usr/bin/env Rscript
suppressMessages({
library(Seurat)
library(ggplot2)
library(argparser)
library(reshape2)
library(dplyr)
library(Cairo)
library(ggnewscale)
})

PRO <- readRDS('sc.seurat.object.rds')

outdir <- '.'
compare <- 'prefix'

intersect_barcode.stat <- readRDS('intersect_barcode.stat.byBarcode.rds')

intersectAddingBarcode <- rownames(intersect_barcode.stat)
write.table(data.frame(intersectAddingBarcode), file.path(outdir, 'intersectAddingBarcode.txt'), row.names = F ,col.names = F, quote = F, sep = '\t')

PRO[['identUSE']] <- PRO@active.ident

PRO$intersectAddingBarcode <- sapply(PRO$addingbarcode, function(x){ ifelse(x %in% intersectAddingBarcode, "intersectAddingBarcode", "Other")})
# saveRDS(PRO, file.path(outdir, 'withintersectAddingBarcode.rds'))


### ------------------- using Seurat Dimplot function -------------------
PRO$intersectAddingBarcode_raw <- sapply(PRO$addingbarcode, function(x){ ifelse(x %in% intersectAddingBarcode, "intersectAddingBarcode", "Other")})
PRO$intersectAddingBarcode_group <- paste0(PRO$intersectAddingBarcode_raw, '_', PRO$group)

lvs = table(PRO$intersectAddingBarcode_group) %>% names() ### "intersectAddingBarcode_PT" "intersectAddingBarcode_RT" "Other_PT"    "Other_RT"
color_use = c("#0067AA","#FF7F00","#bdbdbd","#bdbdbd")
names(color_use) = lvs

### change ident
Idents(PRO) = factor(PRO$intersectAddingBarcode_group, levels = lvs)

p1 <- DimPlot(object=PRO, reduction = 'umap', cols = color_use, order= lvs)
pdf(file.path(outdir,paste0(compare,'.','UMAP','.intersectAddingBarcode.plot.pdf')))
p1
dev.off()
png(file.path(outdir,paste0(compare,'.','UMAP','.intersectAddingBarcode.plot.png')))
p1
dev.off()
## split
p2 <- DimPlot(object=PRO, reduction = 'umap', cols = color_use, split.by = 'sample', , order= lvs)

pdf(file.path(outdir,paste0(compare,'.','UMAP','.intersectAddingBarcode.plot.splitted.pdf')), width = 4 * length(levels(PRO$sample)))
p2
dev.off()
png(file.path(outdir,paste0(compare,'.','UMAP','.intersectAddingBarcode.plot.splitted.png')), width = 400 * length(levels(PRO$sample)))
p2
dev.off()

p3 <- DimPlot(object=PRO, reduction = 'umap', cols = color_use, split.by = 'group', order= lvs)

pdf(file.path(outdir,paste0(compare,'.','UMAP','.intersectAddingBarcode.plot.splitted.group.pdf')), width = 4 * length(levels(PRO$group)))
p3
dev.off()
png(file.path(outdir,paste0(compare,'.','UMAP','.intersectAddingBarcode.plot.splitted.group.png')), width = 400 * length(levels(PRO$group)))
p3
dev.off()


