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
<code>data<-runNormDistance(h1, input="rna")</code>
  
### 3. Matrices fusion
<p> Fuse matrices from multi-modal datasets</p>
<code>data<-fusionMatrix(data, input1="RNA", input2="DNA", method="hadamard")</code>

### 4. Dimension Reduction
<p>Dimension reduction using PCA</p>

<p> Run PCA for DNA matrix:</p>
<code>data<-runPCA(data, input="dna")</code>


<p> Run PCA for RNA matrix:</p>
<code>data<-runPCA(data, input="rna")</code>


<p> Run PCA for fused matrix:</p>
<code>data<-runPCA(data, input="int")</code>

  
To be continued..
