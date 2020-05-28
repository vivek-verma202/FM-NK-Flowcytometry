# FM-NK-Flowcytometry
This pipeline describes agnostic analyses of multi-parameter flow-cytometry data from peripheral blood mononuclear cells (PBMCs) of fibromyalfia cases and controls. Specifically, we focus on NK cells here.


## Preprocessing steps for NK fcs files using FlowJo v_10.6.2:
### step 1: export compensated single PBMCs files for Cytonorm
| Variable            | Files         | 
|--|--| 
|input fcs location  | S:/FM_FLOW_CYTOMETRY/raw/|
|output fcs location | S:/FM_FLOW_CYTOMETRY/cleaned_fcs/1_compensated_raw/D_NK/|
|flowjo workspace    | S:/FM-NK-Flowcytometry/step1.wsp|

### step 2: run cytonorm.R to remove batch effects:
```R
setwd("/scratch/vverma3/FM-NK-Flowcytometry")
## libraries:
x <- c("CytoNorm","ggplot2", "CytoML", "FlowSOM", "pheatmap",
"flowCore","readxl","cowplot","SingleCellExperiment")
lapply(x, require, character.only = TRUE)
## get data
data <- read.delim("S:/FM-NK-Flowcytometry/data.txt")
data$Type <- as.character(data$Type)
data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]
train_data <- dplyr::filter(data, Type == "Train")
validation_data <- dplyr::filter(data, Type == "Validation")
ff <- flowCore::read.FCS(data$Path[1])
channels <- flowCore::colnames(ff)[c(1:15)]
transformList <- flowCore::transformList(channels,
                                         cytofTransform)
transformList.reverse <- flowCore::transformList(channels,
                                                 cytofTransform.reverse)
fsom <- prepareFlowSOM(train_data$Path,
                       channels,
                       nCells = 6000,
                       FlowSOM.params = list(xdim = 5,
                                             ydim = 5,
                                             nClus = 10,
                                             scale = FALSE),
                       transformList = transformList,
                       seed = 1)
cvs <- testCV(fsom,
              cluster_values = c(5, 10, 15))

cvs$pctgs$`10`

model <- CytoNorm.train(files = train_data$Path,
                        labels = train_data$Batch,
                        channels = channels,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = 6000,
                                              xdim = 5,
                                              ydim = 5,
                                              nClus = 10,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)
CytoNorm.normalize(model = model,
                   files = validation_data$Path,
                   labels = validation_data$Batch,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "Normalized",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE)
```

### step 3: Gate NK cells and export
| Variable            | Files         | 
|--|--| 
|input fcs location  | S:/FM-NK-Flowcytometry/Normalized/|
|output fcs location | S:/FM-NK-Flowcytometry/NK_fcs/|
|flowjo workspace    | S:/FM-NK-Flowcytometry/step2.wsp|

### step 4: Prepare data
#### step 4a: Prepare fcs data
```R
library("diffcyt", lib.loc = "/home/vverma3/R/x86_64-redhat-linux-gnu-library/3.6")
library("CATALYST", lib.loc = "/home/vverma3/R/x86_64-redhat-linux-gnu-library/3.6")
files <- list.files(path = "./NK_fcs",
                    pattern = "\\.fcs$", full.names = T)
d_flowSet <- read.flowSet(files, transformation = F,
                          truncate_max_range = F)
```
#### step 4b: experiment_info
```R
filenames <- as.character(pData(d_flowSet)$name)
sample_id <- gsub("^[A-D]_", "", gsub("\\.fcs$", "", filenames))
pheno <- read_excel("./FM_pheno.xlsx")
pheno <- pheno[,c(1:4)]
experiment_info <- pheno[pheno$ID %in% sample_id,]
names(experiment_info) <- c("sample_id","group_id","age","sex")
experiment_info$patient_id <- factor(experiment_info$sample_id)
experiment_info$sex <- factor(experiment_info$sex)
experiment_info$group_id <- factor(experiment_info$group_id,
                                   levels = c("Control","Case"))
experiment_info <- experiment_info[,c(2,5,1,3,4)]

```
#### step 4c: marker_info 
```R
channel_name <- gsub("FJComp-","",gsub("-A$","",colnames(d_flowSet)))
marker_name <- c("TIGIT","CD16","CD57","CD226","CD56","CD107a",
                 "CD335","CD159c","CD158e","CD314","CD96","CD8a","CD159a")
marker_class <- c("state","type","type","state","type","state",
                  "type","type","type","type","state","type","type")
marker_info <- data.frame(channel_name, marker_name, marker_class,
                          stringsAsFactors = F)
marker_info
```
#### step 4d: combine all and transform
```R
d_se <- prepareData(d_flowSet, experiment_info, marker_info)
daf <- daFrame(cth_fs,cth_panel, cth_md, cofactor = 150)
```
### regression matrices set up
```R
design <- createDesignMatrix(experiment_info,
                             cols_design = c("group_id","age","sex"))
contrast <- createContrast(c(0,1,0,0))
nrow(contrast) == ncol(design) # TRUE
data.frame(parameters = colnames(design), contrast)
```
### run FlowSOM:
```R
d_se <- generateClusters(d_se, seed_clustering = 123,
                         meta_clustering = T, meta_k = 20)
```



















