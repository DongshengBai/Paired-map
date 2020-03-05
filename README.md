# Paired-map

Paired-seq/tag Multi-omics Analysis Pipeline

In Test!!

runJaccard and normOVE borrowed from SnapATAC.

## Usage

### 1. preprocess matrices

<p>Load DNA and RNA matrices:</p>

<code>data<-loadPairedDataSet(rna="..path to RNA matrix..", dna="..path to RNA matrix..", project="demo project", min.umi=500, min.fragments = 500)</code>

<p>Filter matrices:</p>
<p>Convert dna matrix into binary matrix</p>
<code>data<-filtMatrix(data, input="dna", low.threshold = -2, high.threshold = 2, method="binary")</code>
and LogMean normallize RNA matrix.
<code>data<-filtMatrix(data, input="rna", low.threshold = -5, high.threshold = 2, method="logMean")</code>

### 2. Calculated distance matrices
<p>Calculate Jaccard matrix for DNA, and normOVE</p>
<code>data<-runDistance(data,input="dna", cell.downsample = 1, feature.downsample = 1, method="jaccard")</code>
<code>data<-runNormDistance(h1, input="dna")</code>
<p>Calculate Euclidean matrix for RNA, and normOVE</p>
<code>data<-runDistance(data,input="rna", cell.downsample = 1, feature.downsample = 1, method="euclidean")</code>
<code>data<-runNormDistance(h1, input="rna", method="none")</code>
  
### 3. Matrices fusion
<p> Fuse matrices from multi-modal datasets</p>
<code>data<-runFusionMatrices(data, use.rna=T, use.dna=T, use.tag=F, method="hadamard")</code>

### 4. Dimension Reduction
<p>Dimension reduction using PCA</p>

<p> Run PCA for DNA matrix:<code>data<-runPCA(data, input="dna")</code></p>
<p> Run PCA for RNA matrix: <code>data<-runPCA(data, input="rna")</code></p>
<p> Run PCA for fused matrix: <code>data<-runPCA(data, input="int")</code></p>

### 5. Clustering
<p>Find clusters using louvain</p>
<p>Run lovain for DNA:<code>data<-runCluster(data,input="dna",use.dims=c(1,10))</code></p>
<p>Run lovain for RNA:<code>data<-runCluster(data,input="rna",use.dims=c(1,10))</code></p>
<p>Run lovain for fused matrix:<code>data<-runCluster(data,input="int",use.dims=c(1,10))</code></p>
  
### 6. Run UMAP
<p>Run UMAP for DNA:<code>data<-runUMAP(data, input="dna", k=15, use.dims=c(1:10))</code></p>
<p>Run UMAP for RNA:<code>data<-runUMAP(data, input="rna", k=15, use.dims=c(1:10))</code></p>
<p>Run UMAP for fused matric:<code>data<-runUMAP(data, input="int", k=15, use.dims=c(1:10))</code></p>
  
### 7. Visuallization
<p>Gene plot:<code>plotFeature(obj, feature="Snap25", feature.type=c("rna", "dna", "int"), embedding.use=c("UMAP", "tSNE"), embedding.type=c("dna", "rna", "int"))</code></p>
<p>Cluster plot:<code>plotFeature(obj, feature="cluster", feature.type=c("rna", "dna", "int"), embedding.use=c("UMAP", "tSNE"), embedding.type=c("dna", "rna", "int"))</code></p>
  
To be continued..
