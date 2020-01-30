### library
print("Loading dependencies.")
{
	library(Matrix)
	library(GenomicRanges)
	library(uwot)
	library(BiocGenerics)
	library(GenomeInfoDb)
	library(graphics)
	library(IRanges)
	library(irlba)
	library(igraph)
	library(parallel)
	library(RANN)
	library(RColorBrewer)
	library(Rtsne)
	library(S4Vectors)
	library(stats)
	library(stats4)
	library(scales)
	library(leiden)
}

print("Loading SnapATAC modules. Citation: Fang et.al. bioRxiv 615179")
{
	### SnapATAC class
	{
		# kgraph class
		{
			methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
			kgraph <- setClass(
			  Class = "kgraph",
			  slots = list(
			    mat = "MatrixOrmatrix",
				file = "character",
			    k = "numeric",
				snn = "logical",
				snn.prune = "numeric"
			  )
			)

			setMethod(
			  f = 'show',
			  signature = 'kgraph',
			  definition = function(object) {
			    cat(
			      'number of cells:', nrow(object@mat), '\n',
			      'K: ', object@k, '\n',
			      'graph file: ', object@file, '\n',
			      'snn: ', object@snn, '\n',
			      'snn.prune: ', object@snn.prune, '\n'
			    )
			  }
			)

			setMethod(
				f="[", 
			    signature = 'kgraph',
				function(x,i,j, drop="missing"){
					.mat = x@mat;
					# a single row or column
			       if(!missing(i)){
					   if(max(i) > nrow(.mat)){
						   stop("idx exceeds number of cells");
					   }
					   if(nrow(.mat) > 0){.mat <- .mat[i,i,drop=FALSE]}
				   }
				   if(!missing(j)){
					   stop("kgraph does not support subsetting for columns");
				   }
				   x@mat = .mat;
				   return(x);
			})

			newKgraph <- function (mat=NULL, file=NULL, k=NULL, snn=NULL, snn.prune=NULL) {
				if(is.null(mat)){
					mat = Matrix::Matrix(0,0,0, sparse=TRUE);
				}
				if(is.null(file)){
					file = character();
				}

				if(is.null(k)){
					k = numeric();
				}

				if(is.null(snn)){
					snn = FALSE;
				}

				if(is.null(snn.prune)){
					snn.prune = numeric();
				}
				
				res = new("kgraph", 
						  mat=mat,
						  file=file,
						  k=k,
						  snn=snn,
						  snn.prune=snn.prune
						  );
				return(res)	
			}

			isKgraphEmpty <- function(obj){
				if(is.null(obj)){
					stop("obj is empty")
				}else{
					if(!is(obj, "kgraph")){
						stop("obj is not a kgraph object");
					}else{
						if((x = nrow(obj@mat)) > 0L){
							return(FALSE)
						}
						if((x=length(obj@file) > 0L)){
							if(file.exists(obj@file)){
								return(FALSE)
							}
						}
					}
				}
				return(TRUE);
			}

			getGraph <- function(obj){
				if(is.null(obj)){
					stop("obj is empty")
				}else{
					if(!is(obj, "kgraph")){
						stop("obj is not a kgraph object");
					}
					if(isKgraphEmpty(obj)){
						stop("obj is empty");			
					}
				}
				
				if((x=nrow(obj@mat) != 0L)){
					return(obj@mat);
				}
				
				if(file.exists(obj@file)){
					edgeList = read.table(obj@file, header=FALSE);
					if((x=ncol(edgeList)) != 3L){
						stop(paste(obj@file, " does not have 3 columns"));
					}
					num.node = max(edgeList[,c(1,2)]);
					M1 = sparseMatrix(i=edgeList[,1], j=edgeList[,2], x=edgeList[,3], dims=c(num.node,num.node));
					M2 = sparseMatrix(i=edgeList[,2], j=edgeList[,1], x=edgeList[,3], dims=c(num.node,num.node));
					M = M1 + M2;
					rm(M1, M2);
					gc();
					return(M);
				}
			}

			writeEdgeListToFile <- function(edges, file, ...
			){
				if(missing(edges) || missing(file)){
					stop("missing edges or file inputs")
				}
				
				if(!file.create(file)){
					stop("fail to create file")
				}			
				
			    write.table(edges, file = file, append = FALSE, quote = FALSE, sep = "\t",
			                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
			                col.names = FALSE)
			}

			readEdgeListFromFile <- function(file)
			{
				edgeList = read.table(file, header=FALSE);
				return(edgeList);
			}
		}

		# jaccard class
		{
			methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
			jaccard <- setClass(
			  Class = "jaccard",
			  slots = list(
			    jmat = "MatrixOrmatrix",
			    nmat = "MatrixOrmatrix",
			    p1 = "numeric",
					p2 = "numeric",
					norm = "logical",
					input.mat="character",
					method = "character"
			  )
			)

			setMethod(
			  f = 'show',
			  signature = 'jaccard',
			  definition = function(object) {
			    cat(
			      'Number of cells:', nrow(object@jmat), '\n',
			      'Number of dims: ', ncol(object@jmat), '\n',
				  	'Input matrix:', object@input.mat, "\n",
			      'Normalized: ', object@norm, '\n'
			    )
			  }
			)

			#' subsetting for jaccard objects
			#'
			#' This function takes a jaccard object and returns the subset of jaccard object.
			#' @param x A jaccard object
			#' @param i selected rows
			#' @param j selected columns
			#' @param drop drop unused levels
			#' @examples
			#' data(demo.sp);
			#' demo.sub.sp = demo.sp[1:5,]
			#' @export
			setMethod(
				f="[", 
			    signature = 'jaccard',
				function(x,i,j, drop="missing"){
					.jmat = x@jmat;
					.nmat = x@nmat;	
					.p1 = x@p1;		
					.p2 = x@p2;			
					# a single row or column
			       if(!missing(i)){
					   if(max(i) > nrow(.jmat)){
						   stop("idx exceeds number of cells");
					   }
					   if(nrow(.jmat) > 0){.jmat <- .jmat[i,,drop=FALSE]}
					   if(nrow(.nmat) > 0){.nmat <- .nmat[i,,drop=FALSE]}
					   if(length(.p1) > 0){.p1 <- .p1[i,drop=FALSE]}
				   }
				   if(!missing(j)){
					   if(max(j) > ncol(.jmat)){
						   stop("idy exceeds number of dimentions");
					   }
			 	 	   if(ncol(.jmat) > 0){.jmat <- .jmat[,j,drop=FALSE]}
					   if(ncol(.nmat) > 0){.nmat <- .nmat[,j,drop=FALSE]}
					   if(length(.p2) > 0){.p2 <- .p2[j,drop=FALSE]}
				   }
				   x@jmat = .jmat;
				   x@p1 = .p1;
				   x@p2 = .p2;
				   x@nmat = .nmat;
				   return(x);
			})

			#' @importFrom methods new
			newJaccard <- function () {
				res = new("jaccard", 
						  jmat=matrix(0,0,0), 
						  p1=numeric(),
						  p2=numeric(),
						  norm=FALSE,
						  nmat=matrix(0,0,0),
						  input.mat="None"
						  )	
			}

			isJaccardComplete <- function (obj) {
				if(missing(obj)){
					stop("obj is missing")
				}else{
					if(!is(obj, "jaccard")){
						stop("obj is not a jaccard object")
					}
				}
				if(!((nrow(obj@jmat) > 0) && (length(obj@p1) > 0) && (length(obj@p2) > 0))){
					return(FALSE);
				}	
				return(TRUE);
			}


			isJaccardNorm <- function (obj) {
				if(missing(obj)){
					stop("obj is missing")
				}else{
					if(!is(obj, "jaccard")){
						stop("obj is not a jaccard object")
					}
				}
				return(obj@norm)
			}
		}

		# dimReduct class
		{
			dim.reduct <- setClass(
			  Class = "dim.reduct",
			  slots = list(
			    imat = "character",
			    dmat = "matrix",
			    sdev = "numeric",
				iter = "numeric",
			    method = "character"
			  )
			)

			setMethod(
			  f = 'show',
			  signature = 'dim.reduct',
			  definition = function(object) {
			    cat(
			      'Input matrix:', object@imat, '\n',
			      'Number of dimensions:', ncol(x = object@dmat), '\n',
			      'Dimentionality reduction method:', object@method, '\n'
			    )
			  }
			)

			setMethod("[", "dim.reduct",
				function(x,i,j, drop="missing"){
					.dmat = x@dmat;
					.sdev = x@sdev;		
					# a single row or column
			       if(!missing(i)){
					   if(max(i) > nrow(.dmat)){
						   stop("idx exceeds number of cells");
					   }
					   if(nrow(.dmat) > 0){.dmat <- .dmat[i,,drop=FALSE]}
				   }
				   if(!missing(j)){
					   if(max(j) > ncol(.dmat)){
						   stop("idy exceeds number of dimentions");
					   }
					   if(length(.sdev) > 0){.sdev <- .sdev[j];}	   
			 	 	   if(ncol(.dmat) > 0){.dmat <- .dmat[,j,drop=FALSE]}
				   }
				   x@dmat = .dmat;
				   x@sdev = .sdev;
				   return(x);
			})

			newDimReduct <- function () {
				res = new("dim.reduct", 
						  dmat=matrix(0,0,0), 
						  sdev=numeric(),
						  method=character(),
						  iter=numeric()
						  )	
				return(res)
			}

			isDimReductComplete <- function(obj) {
				if(missing(obj)){
					stop("obj is missing")
				}else{
					if(!is(obj, "dim.reduct")){
						stop("obj is not a dim.reduct")
					}
				}
				
				if(!((nrow(obj@dmat) > 0) && (length(obj@sdev) > 0))){
					return(FALSE);
				}	
				
				if(ncol(obj@dmat) != length(obj@sdev)){
					return(FALSE);
				}
				
				return(TRUE);
			}

			weightDimReduct <- function(obj, pca.dims, weight.by.sd=TRUE){
				if(missing(obj)){
					stop("obj is missing")
				}else{
					if(!is(obj, "dim.reduct")){
						stop("obj is not a dim.reduct")
					}
				}
				
				if(weight.by.sd){
					data.use = obj@dmat[,pca.dims] %*% diag(obj@sdev[pca.dims]) ;
				}else{
					data.use = obj@dmat[,pca.dims];
				}
				return(data.use);
			}

			dimReductDim <- function(obj){
				if(missing(obj)){
					stop("obj is missing")
				}else{
					if(!is(obj, "dim.reduct")){
						stop("obj is not a dim.reduct")
					}
				}
				return(length(obj@sdev));
			}
		}

		# snap class
		{
			methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
			setClass("snap",
				slots=list(
				des="character",
				barcode="character",
				file="character",
				sample="character",
				metaData="data.frame",
				feature="GRanges",
				peak = "GRanges",
				bmat = "Matrix",
				pmat = "Matrix",
				gmat = "Matrix",
				mmat = "matrix",
				jmat = "jaccard",
				regModel = "numeric",
				smat = "dim.reduct",
				graph = "kgraph",
				tsne = "MatrixOrmatrix",
				umap = "MatrixOrmatrix",
				cluster = "factor"
				)
			)

			.valid.snap.feature <- function(object)
			{
				if(length(object@feature) != ncol(object@bmat)){
					return("slot 'feature' have different length from 'bmat' column")		
				}
				NULL;
			}

			.valid.snap.peak <- function(object)
			{
				if(length(object@peak) != ncol(object@pmat)){
					return("slot 'peak' have different length from 'pmat' column")		
				}
				NULL;
			}

			.valid.snap.barcode <- function(object)
			{
				if(length(object@barcode) != nrow(object@metaData)){
					return("slot 'barcode' have different length from 'metaData'")		
				}
				NULL
			}

			.valid.snap <- function(object)
			{
			    #c(.valid.snap.barcode(object), .valid.snap.feature(object))
			    c(.valid.snap.barcode(object), .valid.snap.peak(object), .valid.snap.feature(object));
			}
			methods::setValidity("snap", .valid.snap)

			setMethod("show", signature = "snap",
				definition = function(object) {
					if((x=length(object@des)) > 0L){
						cat("description: ", object@des, "\n");
						cat("===================", "\n");			
					}
					cat("number of barcodes: ", ifelse(is.null(length(object@barcode)), 0, length(object@barcode)), "\n", sep="");
					cat("number of bins: ", ncol(object@bmat), "\n", sep="");
					cat("number of genes: ", ncol(object@gmat), "\n", sep="");
					cat("number of peaks: ", ncol(object@pmat), "\n", sep="");
					cat("number of motifs: ", ncol(object@mmat), "\n", sep="");
				}                              
			)
		}
	}
	### SnapATAC methods
	{
		###
		### methods 
		###
		newSnap <- function () {
			metaData=data.frame()
			des = character()
			file = as.character(c())
			sample = as.character(c())
			barcode = as.character(c())
			feature = GRanges()
			peak = GRanges()
			metaData = data.frame()
			bmat=Matrix(nrow=0, ncol=0, sparse=TRUE)
			pmat=Matrix(nrow=0, ncol=0, sparse=TRUE)
			gmat=Matrix(nrow=0, ncol=0, sparse=TRUE)
			mmat=matrix(0,0,0)
			jmat=newJaccard()
			smat=newDimReduct()
			graph=newKgraph()
			regModel=c()
			tsne=matrix(nrow=0, ncol=0)
			umap=matrix(nrow=0, ncol=0)
			cluster=factor()
			res = new("snap", 
					  des=des,
					  file=file,
					  sample=sample,
					  barcode=barcode, 
					  feature=feature, 
					  peak=peak, 
					  metaData=metaData, 
					  bmat=bmat, 
					  pmat=pmat, 
					  gmat=gmat,
					  mmat=mmat, 
					  jmat=jmat, 
					  smat=smat, 
					  graph=graph, 
					  tsne=tsne, 
					  umap=umap, 
					  cluster=cluster
					  )
		}

		setGeneric("summarySnap", function(obj) standardGeneric("summarySnap"))
		setGeneric("is.snap", function(obj) standardGeneric("is.snap"))
		setMethod("is.snap", "snap", function(obj) return(is(obj, "snap")))

		setMethod("summarySnap", "snap", function(obj){
			if(nrow(obj@metaData) == 0){stop("metaData is empty")}
			barcode = obj@metaData;
			message("Total  number of barcodes: ", length(obj@barcode))
			message("Median number of sequencing fragments: ", median(barcode$TN))
			message("Median number of uniquely mapped fragments: ", median(barcode$UQ))
			message("Median number of mappability ratio: ", round(median((barcode$UM+1)/(barcode$TN+1)),2))
			message("Median number of properly paired ratio: ", round(median((barcode$PP+1)/(barcode$UM+1)),2))
			message("Median number of duplicate ratio: ", round(median(1 - (barcode$UQ+1)/(barcode$PP+1)),2))
			message("Median number of chrM ratio: ", round(median((barcode$CM+1) / (barcode$UQ+1)),2))
			message("Median number of unique molecules (UMI): ", median(barcode$UQ));
		})

		setMethod("nrow", "snap", function(x) length(x@barcode))

		setMethod("colSums", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=TRUE){
			mat = match.arg(mat)
			mat.use = methods::slot(x, mat)
			if((x=nrow(mat.use))==0L){
				stop("mat is empty")
			}
			res = Matrix::colSums(mat.use, na.rm)
			return(res)
		})
		setMethod("rowSums", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=TRUE){
			mat = match.arg(mat)
			mat.use = methods::slot(x, mat)
			if((x=nrow(mat.use))==0L){
				stop("mat is empty")
			}
			res = Matrix::rowSums(mat.use, na.rm)
			return(res)
		})
		setMethod("rowMeans", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=TRUE){
			mat = match.arg(mat)
			mat.use = methods::slot(x, mat)
			if((x=nrow(mat.use))==0L){
				stop("mat is empty")
			}
			res = Matrix::rowMeans(mat.use, na.rm)
			return(res)
		})
		setMethod("colMeans", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=TRUE){
			mat = match.arg(mat)
			mat.use = methods::slot(x, mat)
			if((x=nrow(mat.use))==0L){
				stop("mat is empty")
			}
			res = Matrix::colMeans(mat.use, na.rm)
			return(res)
		})

		setMethod("[", "snap",
			function(x,i,j,mat=c("bmat", "pmat", "gmat"), drop="missing"){
				.barcode = x@barcode
				.file = x@file
				.sample = x@sample
				.feature = x@feature
				.peak = x@peak
				.bmat = x@bmat
				.pmat = x@pmat
				.gmat = x@gmat
				.mmat = x@mmat
				.jmat = x@jmat
				.smat = x@smat
				.graph = x@graph
				.cluster = x@cluster
				.tsne = x@tsne
				.umap = x@umap
				.metaData = x@metaData
		       if(!missing(i)){
				   if(max(i) > nrow(x)){
					   stop("idx exceeds number of cells");
				   }
				   if(nrow(.bmat) > 0){.bmat <- .bmat[i,,drop=FALSE]}
				   if(nrow(.pmat) > 0){.pmat <- .pmat[i,,drop=FALSE]}
				   if(nrow(.gmat) > 0){.gmat <- .gmat[i,,drop=FALSE]}	   
				   if(nrow(.mmat) > 0){.mmat <- .mmat[i,,drop=FALSE]}	   
				   if(nrow(.jmat@jmat) > 0){.jmat <- .jmat[i,,drop=FALSE]}
				   if(nrow(.smat@dmat) > 0){.smat <- .smat[i,,drop=FALSE]}
				   if(nrow(.tsne) > 0){.tsne <- .tsne[i,,drop=FALSE]}
				   if(nrow(.umap) > 0){.umap <- .umap[i,,drop=FALSE]}
				   if(nrow(.graph@mat) > 0){.graph <- .graph[i,,drop=FALSE]}
				   if(nrow(.metaData) > 0){.metaData <- .metaData[i,,drop=FALSE]}
				   if(length(.cluster) > 0){.cluster <- .cluster[i,drop=FALSE]}
				   if(length(.barcode) > 0){.barcode <- .barcode[i,drop=FALSE]}
				   if(length(.file) > 0){.file <- .file[i,drop=FALSE]}
				   if(length(.sample) > 0){.sample <- .sample[i,drop=FALSE]}
			   }
			   if(!missing(j)){
		   			mat = match.arg(mat);
			   		if(mat == "bmat"){
			 		   if(ncol(.bmat) > 0){.bmat <- .bmat[,j,drop=FALSE]}
					   if(length(.feature) > 0){.feature <- .feature[j];}	   
			   		}else if(mat == "pmat"){
		 	 		   if(ncol(.pmat) > 0){.pmat <- .pmat[,j,drop=FALSE]}
					   if(length(.peak) > 0){.peak <- .peak[j];}	   
			   		}else if(mat == "gmat"){
		 	 		   if(ncol(.gmat) > 0){.gmat <- .gmat[,j,drop=FALSE]}
			   		}
			   }
			   x@bmat = .bmat
			   x@pmat = .pmat
			   x@gmat = .gmat
			   x@mmat = .mmat
			   x@barcode = .barcode
			   x@file = .file
			   x@sample = .sample
			   x@peak = .peak
			   x@feature = .feature
			   x@metaData = .metaData
			   x@umap = .umap
			   x@feature = .feature
			   x@jmat = .jmat
			   x@smat = .smat
			   x@graph = .graph
			   x@cluster = .cluster
			   x@tsne = .tsne
			   return(x)
		})
		findCentrod <- function(x, y){
			x.ls = split(data.frame(x),y);
			centroid.ls = lapply(split(data.frame(x),y), function(xx) apply(xx, 2, median))
			centroid.df = data.frame(do.call(rbind, centroid.ls))
			centroid.df$Y = names(centroid.ls);
			
			return(centroid.df);		
		}
	}
	### SnapATAC functions
	{
		###
		### functiosn snapATAC
		###
		createSnapFromPmat <- function(mat, barcodes, peaks) {
		  UseMethod("createSnapFromPmat", mat)
		}
		createSnapFromPmat.default <- function(mat, barcodes, peaks){
			if(missing(mat) || missing(barcodes) || missing(peaks)){
				stop("mat or barcodes or peaks is missing")
			}

			if(!(is(mat, "dsCMatrix") || is(mat, "dgCMatrix") || is(mat, "dgTMatrix"))){
				stop("'mat' is not a sparse matrix")
			}

			if(length(barcodes) != nrow(mat)){
				stop("'mat' has different number of rows with number of barcodes")
			}
			
			if(!is(peaks, "GRanges")){
				stop("'peaks' is not a GRanges object")
			}
			if(length(peaks) != ncol(mat)){
				stop("'mat' has different number of columns with number of peaks")
			}
			
			obj = newSnap()
			obj@pmat = mat
			obj@barcode = barcodes
			obj@peak = peaks
			return(obj)
		}

		createSnapFromGmat <- function(mat, barcodes, gene.names) {
		  UseMethod("createSnapFromGmat")
		}
		createSnapFromGmat.default <- function(mat, barcodes, gene.names){
			if(missing(mat) || missing(barcodes) || missing(gene.names)){
				stop("mat or barcodes or gene.names is missing")
			}

			if(!(is(mat, "dsCMatrix") || is(mat, "dgCMatrix") || is(mat, "dgTMatrix"))){
				stop("'mat' is not a sparse matrix")
			}

			if(length(barcodes) != nrow(mat)){
				stop("'mat' has different number of rows with number of barcodes")
			}
			
			if(!is(gene.names, "character")){
				stop("'gene.names' is not a character object")
			}
			if(length(gene.names) != ncol(mat)){
				stop("'mat' has different number of columns with number of gene.names")
			}
			
			obj = newSnap()
			obj@gmat = mat
			obj@barcode = barcodes
			colnames(obj@gmat) = gene.names
			return(obj)
		}

		makeBinary <- function(obj, mat, outlier.filter) {
		  UseMethod("makeBinary", obj);
		}
		makeBinary.default <- function(obj, mat=c("bmat", "pmat", "gmat"), outlier.filter=1e-3){
			if(!is(obj, "snap")){
				stop("obj is not a snap obj")
			}
			
			if(!is.numeric(outlier.filter) || outlier.filter < 0 || outlier.filter > 1){
				stop("incorrect outlier.filter")		
			}
			
			mat = match.arg(mat);
			if(mat == "bmat"){
				if(nrow(obj@bmat) == 0){
					stop("@bmat does not exist")
				}
				
				if(max(obj@bmat) == 1){
					stop("@bmat is already binarized")
				}

				x = obj@bmat;
				count = x@x;
				count_cutoff = max(1, quantile(count, 1 - outlier.filter));
				x@x[x@x > count_cutoff] = 0
				x@x[x@x > 0] = 1
				obj@bmat = x;
			}else if(mat == "pmat"){
				if(nrow(obj@pmat) == 0){
					stop("@pmat does not exist")
				}
				
				if(max(obj@pmat) == 1){
					stop("@pmat is already binarized")
				}

				x = obj@pmat;
				count = x@x;
				count_cutoff = max(1, quantile(count, 1 - outlier.filter));
				x@x[x@x > count_cutoff] = 0
				x@x[x@x > 0] = 1
				obj@pmat = x;
			}else if(mat == "gmat"){
				if(nrow(obj@gmat) == 0){
					stop("@gmat does not exist")
				}
				
				if(max(obj@gmat) == 1){
					stop("@gmat is already binarized")
				}

				x = obj@gmat;
				count = x@x;
				count_cutoff = max(1, quantile(count, 1 - outlier.filter));
				x@x[x@x > count_cutoff] = 0
				x@x[x@x > 0] = 1
				obj@gmat = x;
			}
			return(obj);
		}

		filterBins <- function(obj, low.threshold, high.threshold, mat) {
		  UseMethod("filterBins", obj);
		}
		filterBins.default <- function(obj, low.threshold=-2, high.threshold=2, mat=c("bmat", "pmat")){
			if(missing(obj)){
				stop("obj is missing")
			}else{
				if(!is(obj, "snap")){
					stop("'obj' is not a snap obj")
				};		
			}
			
			mat = match.arg(mat);
			data.use = methods::slot(obj, mat);

			if((x=nrow(data.use)) == 0L){		
				stop("count matrix is empty")
			}

			idy = which(Matrix::colSums(data.use) > 0);
			cov = log(Matrix::colSums(data.use)[idy] + 1, 10);
			zcov = (cov - mean(cov)) / stats::sd(cov);	
			idy2 = which(zcov >= low.threshold & zcov <= high.threshold);
			idy = idy[idy2];
			methods::slot(obj, mat) = data.use[,idy,drop=FALSE];
			if(mat=="bmat"){
				obj@feature = obj@feature[idy];
			}else if(mat=="pmat"){
				obj@peak = obj@peak[idy];		
			}
			return(obj)
		}

		calJaccard <- function(X_i, X_j){
			A = Matrix::tcrossprod(X_i, X_j);
			bi = Matrix::rowSums(X_i);
			bj = Matrix::rowSums(X_j);
			jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));
			return(jmat)				
		}

		runJaccard <- function(obj, bin.downsample, mat, max.var, seed.use) {
		  UseMethod("runJaccard", obj);
		}
		runJaccard.default <- function(
			obj, 
			bin.downsample=1,
			mat = c("bmat", "pmat", "gmat"),
			max.var = 1000, 
			seed.use=10
		){
			if(missing(obj)){
				stop("obj is missing")
			}else{
				if(!is(obj, "snap")){
					stop("obj is not a snap obj")
				}
			}
			mat = match.arg(mat)
			mat.use = methods::slot(obj, mat)
			p1 = Matrix::rowMeans(mat.use)
			if((x=nrow(mat.use)) == 0L){
				stop("input matrix is empty")
			}
			if((x=max(mat.use)) > 1L){
				stop("input matrix is not a binary matrix, run 'makeBinary' first")	
			}
			if(any(Matrix::rowSums(mat.use) == 0)){
				stop("input matrix contains empty rows, remove empty rows first")	
			}
			col.covs = log(Matrix::colSums(mat.use)+1, 10)
			mat.use = mat.use[,which(col.covs > 0)]
			col.covs = col.covs[which(col.covs > 0)]

			# randomly select a subset of cells as reference 
			if(bin.downsample > 1 | bin.downsample <= 0){
				stop("bin.downsample must be between 0 and 1");
			}else{
				if(bin.downsample < 1){
					col.covs.dens <- density(x = col.covs, bw = 'nrd', adjust = 1)
					sampling_prob <- 1 / (approx(x = col.covs.dens$x, y = col.covs.dens$y, xout = col.covs)$y + .Machine$double.eps)
					set.seed(seed.use);
					idy <- sort(sample(x = seq(col.covs), size = bin.downsample*length(col.covs), prob = sampling_prob));
					mat.use = mat.use[,idy];				
				} 
			}
			
			# down-sample for jaccard index matrix
			if(max.var < nrow(obj)){
				max.var = min(max.var, nrow(mat.use));
				row.covs = log(Matrix::rowSums(mat.use)+1,10);		
				row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1)
				sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
				set.seed(seed.use);
				idx <- sort(sample(x = seq(row.covs), size = max.var, prob = sampling_prob));
				mat.ref = mat.use[idx,];				
				p2 = p1[idx];		
			}else{
				mat.ref = mat.use;				
				p2 = p1;				
			}
			jmat = calJaccard(mat.use, mat.ref);

			rm(mat.ref);
			rm(mat.use);
			# remove large objects
			obj@jmat@jmat = jmat;
			obj@jmat@p1 = p1;
			obj@jmat@p2 = p2;
			obj@jmat@norm = FALSE;
			obj@jmat@method = character();
			return(obj);
		}

		runNormJaccard <- function(obj, method, row.center, row.scale, low.threshold, high.threshold, do.par, ncell.chunk, num.cores, seed.use){
		  UseMethod("runNormJaccard");
		}
		runNormJaccard.default <- function(
			obj, 
			method=c("residual", "zscore"), 
			row.center=TRUE,
			row.scale=TRUE, 
			low.threshold=-5, 
			high.threshold=5, 
			do.par=FALSE,
			ncell.chunk=1000, 
			num.cores=1,
			seed.use=10
		){
			tmp.folder="./tmp_norm/"
			if(!dir.exists(tmp.folder)){
				system(paste("mkdir ", tmp.folder, sep=""))
			}

			if(missing(obj)){
				stop("obj is missing")
			}else{
				if(!is(obj, "snap")){
					stop("'obj' is not a snap obj")
				}
			}
			
			if(missing(tmp.folder)){
				stop("tmp.folder is missing")
			}else{
				if(!dir.exists(tmp.folder)){
					stop("tmp.folder does not exist");			
				}
			}

			if(!isJaccardComplete(obj@jmat)){
				stop("jaccard object is not complete, run 'runJaccard' first")
			}else{
				if(isJaccardNorm(obj@jmat)){
					stop("jaccard index matrix has been normalized")
				}		
			}
			
			if(!is.logical(row.center)){
				stop("row.center is not a logical")
			}

			if(!is.logical(row.scale)){
				stop("row.scale is not a logical")
			}
			
			if(low.threshold > high.threshold){
				stop("low.threshold must be smaller than high.threshold");
			}

			if(low.threshold > 0 || high.threshold < 0){
				stop("low.threshold must be smaller than 0 and high.threshold must be greater than 0");
			}
			
			method = match.arg(method);
			jmat = obj@jmat@jmat;
			b1 = obj@jmat@p1;
			b2 = obj@jmat@p2;

			if(do.par){
			    # input checking for parallel options
				if(num.cores > 1){
			        if (num.cores == 1) {
			          num.cores = 1
			        } else if (num.cores > detectCores()) {
			          num.cores <- detectCores() - 1
			          warning(paste0("num.cores set greater than number of available cores(", parallel::detectCores(), "). Setting num.cores to ", num.cores, "."))
			        }
			      } else if (num.cores != 1) {
			        num.cores <- 1
				}
			
				# step 2) slice the orginal obj into list
				id = seq(nrow(obj));
				id.ls = split(id, ceiling(seq(id)/ncell.chunk));
				
				if(length(id.ls) > 1){
					id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
					# remove the last item of the list
					id.ls = id.ls[-length(id.ls)];
				}	
				
				prefix_tmp = tempfile(pattern = "file", tmpdir = tmp.folder);
				backingfile_tmp <- paste(prefix_tmp, ".bin", sep="");
				descriptorfile_tmp <- paste(prefix_tmp, ".desc", sep="");
			
				x <- as.big.matrix(x = obj@jmat@jmat, 
								   type = "double", 
				                   separated = FALSE, 
								   backingpath=tmp.folder,
				                   backingfile = basename(backingfile_tmp), 
				                   descriptorfile = basename(descriptorfile_tmp)
								   );
			
				cl <- makeCluster(num.cores);
				registerDoParallel(cl);	
				
				nmat <- foreach(i=1:length(id.ls), .verbose=FALSE, .packages="bigmemory", .combine = "rbind") %dopar% {
				    t_mat <- attach.big.matrix(descriptorfile_tmp);
					return(normObservedJmat2(jmat=t_mat[id.ls[[i]],], b1=b1[id.ls[[i]]], b2=b2, method=method));
				}
				
				stopCluster(cl);
				closeAllConnections();
				rm(x);
				file.remove(backingfile_tmp);
				file.remove(descriptorfile_tmp);
				gc();
			}else{
				model.init = trainRegressModel(jmat, b1, b2);
				nmat = normObservedJmat(obj@jmat@jmat, model.init, obj@jmat@p1, obj@jmat@p2, method=method);
			}
			
			if(row.center || row.scale){
				nmat = t(scale(t(nmat), center=row.center, scale=row.scale));
			}
			
			nmat[nmat >= high.threshold] = high.threshold;
			nmat[nmat <= low.threshold]  = low.threshold;

			obj@jmat@nmat = nmat;
			obj@jmat@method = method;
			obj@jmat@norm = TRUE;	
			return(obj);
		}
		.normOVE <- function(p1, p2){
		    pp = tcrossprod(p1, p2);
			ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
			ee = pp/(ss - pp)
			return(ee)	
		}
		trainRegressModel <- function(jmat, b1, b2){
			# remove the diag elements in the jmat
			idx.pairwise = which(jmat == 1, arr.ind=TRUE);

			# calculate the expected jaccard index matrix given the read depth
			emat = .normOVE(b1, b2);

			# estimate the global scaling factor
			scale.factor = mean(jmat / emat);

			# fill the missing value for the diagnoal elements
			jmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
			data = data.frame(x=c(emat), y=c(jmat));	
			# 2. polynomial regression
			model <- lm(y ~ x + I(x^2), data);
			return(model);	
		}
		normObservedJmat <- function(jmat, model, b1, b2, method){
			.normOVE2 <- function(p1, p2){
			    pp = tcrossprod(p1, p2);
				ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
				ee = pp/(ss - pp)
				return(ee)	
			}
			# 1. remove the "1" elements in the jaccard matrix
			idx.pairwise = which(jmat == 1, arr.ind=TRUE);
			emat = .normOVE2(b1, b2);
			scale.factor = mean(jmat / emat);
			jmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
			
			# 2. Expansion parameters from subset of cells to all cells
			preds = predict(model, data.frame(x=c(emat)), se.fit = TRUE)
			
			# 3. calculate residuals or zscore
			if(method == "zscore"){
				norm = (c(jmat) - preds$fit) / (preds$se.fit);		
			}else if(method == "residual"){
				norm = c(jmat) -  preds$fit;
			}
			nmat = matrix(norm, nrow(emat), ncol(emat));	
			return(nmat);
		}
		normObservedJmat2 <- function(jmat, b1, b2, method){
			.normOVE2 <- function(p1, p2){
			    pp = tcrossprod(p1, p2);
				ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
				ee = pp/(ss - pp)
				return(ee)	
			}
			
			# remove the diag elements in the jmat
			idx.pairwise = which(jmat == 1, arr.ind=TRUE);

			# calculate the expected jaccard index matrix given the read depth
			emat = .normOVE(b1, b2);

			# estimate the global scaling factor
			scale.factor = mean(jmat / emat);

			# fill the missing value for the diagnoal elements
			jmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
			data = data.frame(x=c(emat), y=c(jmat));	
			
			# 2. polynomial regression
			model <- lm(y ~ x + I(x^2), data);
			
			# 2. Expansion parameters from subset of cells to all cells
			preds = predict(model, data.frame(x=c(emat)), se.fit = TRUE);
			
			# 3. calculate residuals or zscore
			if(method == "zscore"){
				norm = (c(jmat) - preds$fit) / (preds$se.fit);		
			}else if(method == "residual"){
				norm = c(jmat) -  preds$fit;
			}
			nmat = matrix(norm, nrow(emat), ncol(emat));	
			return(nmat);
		}

		runDimReduct <- function(obj, pc.num, input.mat, method, center, scale, seed.use, maxit, ...){
		  UseMethod("runDimReduct", obj)
		}
		runDimReduct.default <- function(
			obj, 
			pc.num=50,
			input.mat = c("jmat", "gmat", "bmat", "gmat", "nmat"), 
			method=c("svd", "pca.whiten", "pca"), 
			center=TRUE, 
			scale=FALSE, 
			seed.use=10,
			maxit=1000,
			...
		){
			if(!is(obj, "snap")){
				stop("obj is not a snap obj")
			}
			input.mat = match.arg(input.mat)	
			if(input.mat == "jmat"){
				x = obj@jmat@jmat;
			}else if(input.mat == "bmat"){
				x = obj@bmat;
			}else if(input.mat == "pmat"){
				x = obj@pmat;
			}else if(input.mat == "gmat"){
				x = obj@gmat;
			}else if(input.mat == "nmat"){
				x = obj@jmat@nmat;
			}else{
				stop("input.mat does not exist in obj")
			}
			if(nrow(x) * ncol(x) == 0){
				stop("input.mat is empty")		
			}
			ncell = nrow(obj);
			if(nrow(x) != ncell){
				stop("input.mat has wrong number of cells");
			}
			method = match.arg(method);
			nvar = ncol(x); # number of variables
			if (is.numeric(pc.num)) {
			    pc.num <- as.integer(pc.num)
			  } else{
		  		stop("pc.num must be an integer")	  	
			}
		  if (pc.num > min(ncell, nvar)) {
		      message("'n.comp' is too large: reset to ", min(ncell, nvar))
		      pc.num <- min(nvar, ncell)
		  }
			if(!is.logical(center)){
				stop("center must be logical variable TRUE or FALSEs")
			}
			if(!is.logical(scale)){
				stop("scale must be logical variable TRUE or FALSEs")
			}
			x <- t(x)
		    a <- names(as.list(match.call()))
		    ans <- list(scale=scale)
		    if (!is.matrix(x)) x <- as.matrix(x)
		    args <- list(A=x, nv=pc.num)
		    if (is.logical(center))
		    {
		      if (center) args$center <- colMeans(x)
		    } else{
		    	stop("center must be logical")
		    }
		    if (is.logical(scale))
		    {
		        if (is.numeric(args$center))
		        {
		          f <- function(i) sqrt(sum((x[, i] - args$center[i]) ^ 2) / (nrow(x) - 1L))
		          scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE)
		          if (ans$scale) ans$totalvar <- ncol(x)
		          else ans$totalvar <- sum(scale. ^ 2)
		        } else
		        {
		          if (ans$scale)
		          {
		            scale. <- apply(x, 2L, function(v) sqrt(sum(v ^ 2) / max(1, length(v) - 1L)))
		            f <- function(i) sqrt(sum((x[, i] / scale.[i]) ^ 2) / (nrow(x) - 1L))
		            ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE) ^ 2)
		          } else
		          {
		            f <- function(i) sum(x[, i] ^ 2) / (nrow(x) - 1L)
		            ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE))
		          }
		        }
		        if (ans$scale) args$scale <- scale.
		    } else{
		    	stop("scale must be logical")
		    }
			if(center){
				x.norm = sweep(args$A, 2, args$center, FUN=`-`)		
			}else{
				x.norm = args$A;		
			}
			if(scale){
				x.norm = sweep(x.norm, 2, args$scale, FUN=`/`)		
			}else{
				x.norm = x.norm
			}
			if(method == "svd"){
				set.seed(seed.use)
				S <- irlba(A = x.norm, nv = pc.num, nu = pc.num, maxit=maxit, scale.=FALSE, center=FALSE, ...);
				obj@smat@dmat = S$v;
				obj@smat@sdev = S$sdev;
				obj@smat@iter = maxit;
				obj@smat@method = "svd";
			}else if(method == "pca.whiten"){
				# calculate covariance matrix
				V <- x %*% t(x) / ncell;
				# calculate SVD against the covariance matrix
				set.seed(seed.use);
				S <- irlba(A = V, nv = pc.num, nu = pc.num);
				# PCA whitening
				D <- diag(c(1/sqrt(S$d)))
				K <- D %*% t(S$u)
				obj@smat@dmat = t(K %*% x);
				obj@smat@sdev = S$sdev;
				obj@smat@iter = maxit;
				obj@smat@method = "pca.whiten";		
			}else{
				set.seed(seed.use);
				S <- prcomp_irlba(x.norm, n=pc.num, scale.=FALSE, center=FALSE, maxit=maxit, ...);
				obj@smat@dmat = S$x;
				obj@smat@sdev = S$sdev;
				obj@smat@iter = maxit;
				obj@smat@method = "pca";
			}
			obj@smat@imat = input.mat;
			return(obj);
		}

		runViz<- function(obj, dims, eigs.dims, method, fast_tsne_path, Y.init, seed.use, num.cores, ...) {
		  UseMethod("runViz", obj);
		}
		runViz.default <- function(
			obj, 
			dims=2,
			eigs.dims=NULL, 
			method=c("Rtsne", "umap", "fast_tsne"), 
			fast_tsne_path = NULL, 
			Y.init=NULL,
			seed.use=131,
			num.cores=8,
			...
			){
			tmp.folder="./tmp_norm/"
			if(!dir.exists(tmp.folder)){
				system(paste("mkdir ", tmp.folder, sep=""))
			}
			# check input
			if(!is(obj, "snap")){
				stop("obj is not a snap object")
			}
			
			if(nrow(obj@smat@dmat) == 0){
				stop("PCA is empty, runDimReduct first")
			}

			if(nrow(obj@smat@dmat) != nrow(obj)){
				stop("PCA has different length with obj, data has been subsetted by mistake")
			}

			
			# check PCA dimentions
			ncell = nrow(obj);
			nvar = ncol(obj@smat@dmat);	
			if(is.null(eigs.dims)){
				eigs.dims=1:nvar;	
			}else{
				if(any(eigs.dims > nvar) ){
					stop("'eigs.dims' exceeds reduced dimentions variables number");
				}		
			}
				
			data.use = obj@smat@dmat[,eigs.dims];
			
			if(missing(tmp.folder)){
				stop("tmp.folder is missing")
			}else{
				if(!dir.exists(tmp.folder)){
					stop("tmp.folder does not exist");			
				}
			}
			
			# check input parameters
			method = match.arg(method);
			# check if fi method exists or not
			if(method=="fast_tsne"){
				if(is.null(fast_tsne_path)){
					stop("fast_tsne_path is missing");
				}

				if(!file.exists(fast_tsne_path)){
					stop("'fast_tsne_path' fast tsne does not exist")
				}		

				fast_tsne_path <- normalizePath(fast_tsne_path);
				if (!file_test('-x', fast_tsne_path)) {
					stop(fast_tsne_path, " is not executable; check your fast_tsne_path parameter")
				}				
				obj@tsne = fftRtsne(
						X = data.use, 
						dims=dims, 
						fast_tsne_path=fast_tsne_path, 
						rand_seed=seed.use,
						nthreads=num.cores,
						initialization=Y.init,
						tmp.folder=tmp.folder,
						...
					);	
				colnames(obj@tsne) = c("tsne-1", "tsne-2")				
			}else if(method=="Rtsne"){
				set.seed(seed.use);
				obj@tsne = Rtsne(
					data.use, 
					dims=dims, 
					verbose = FALSE, 
					pca = FALSE, 
					is_distance = FALSE, 
					check_duplicates=FALSE,
					num_threads=num.cores,
					Y_init=Y.init,
					rand_seed=seed.use,
					...
					)$Y;
				colnames(obj@tsne) = c("tsne-1", "tsne-2")
			}else{
				if (requireNamespace("umap", quietly = FALSE)) {
						set.seed(seed.use);
						obj@umap = uwot::umap(data.use, n_threads=num.cores, verbose=TRUE, ...)
						colnames(obj@umap) = c("umap-1", "umap-2")
				  } else {
				      stop("Please install umap - learn more at https://cran.r-project.org/web/packages/umap/index.html")
				  }
				  
			}
			return(obj);
		}


		fftRtsne <- function(X, 
				     dims=2, perplexity=30, theta=0.5,
				     check_duplicates=TRUE,
				     max_iter=1000,
				     fft_not_bh = TRUE,
				     ann_not_vptree = TRUE,
				     stop_early_exag_iter=250,
				     exaggeration_factor=12.0, no_momentum_during_exag=FALSE,
				     start_late_exag_iter=-1.0,late_exag_coeff=1.0,
		             mom_switch_iter=250, momentum=.5, final_momentum=.8, learning_rate=200,
				     n_trees=50, search_k = -1,rand_seed=-1,
				     nterms=3, intervals_per_integer=1, min_num_intervals=50, 
				     K=-1, sigma=-30, initialization=NULL,
				     data_path=NULL, result_path=NULL,
				     load_affinities=NULL,
				     fast_tsne_path=NULL, nthreads=0, perplexity_list = NULL, 
		             get_costs = FALSE, df = 1.0,
					 tmp.folder,
					 ... ) {
		        version_number = '1.1.0'

			if (is.null(fast_tsne_path)) {
				stop("fast_tsne_path is NULL")
			}
			
			if(missing(tmp.folder)){
				stop("tmp.folder is missing")
			}else{
				if(!dir.exists(tmp.folder)){
					stop("tmp.folder does not exist");			
				}
			}
			
			if (is.null(data_path)) {
				data_path <- tempfile(pattern='fftRtsne_data_', tmpdir = tmp.folder, fileext='.dat')
			}
			if (is.null(result_path)) {
				result_path <- tempfile(pattern='fftRtsne_result_', tmpdir = tmp.folder, fileext='.dat')
			}
			if (is.null(fast_tsne_path)) {
				fast_tsne_path <- system2('which', 'fast_tsne', stdout=TRUE)
			}

			fast_tsne_path <- normalizePath(fast_tsne_path)
			if (!file_test('-x', fast_tsne_path)) {
				stop(fast_tsne_path, " does not exist or is not executable; check your fast_tsne_path parameter")
			}

			is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

			if (!is.numeric(theta) || (theta<0.0) || (theta>1.0) ) { stop("Incorrect theta.")}
			if (nrow(X) - 1 < 3 * perplexity) { stop("Perplexity is too large.")}
			if (!is.matrix(X)) { stop("Input X is not a matrix")}
			if (!(max_iter>0)) { stop("Incorrect number of iterations.")}
			if (!is.wholenumber(stop_early_exag_iter) || stop_early_exag_iter<0) { stop("stop_early_exag_iter should be a positive integer")}
			if (!is.numeric(exaggeration_factor)) { stop("exaggeration_factor should be numeric")}
			if (!is.numeric(df)) { stop("df should be numeric")}
			if (!is.wholenumber(dims) || dims<=0) { stop("Incorrect dimensionality.")}
			if (search_k == -1) {
		       if (perplexity>0) {
		          search_k = n_trees*perplexity*3
		       } else if (perplexity==0) {
		          search_k = n_trees*max(perplexity_list)*3
		       } else { 
		          search_k = n_trees*K
		       }
		    }

			if (fft_not_bh){
			  nbody_algo = 2;
			}else{
			  nbody_algo = 1;
			}

			if (is.null(load_affinities)) {
				load_affinities = 0;
			} else {
				if (load_affinities == 'load') {
					load_affinities = 1;
				} else if (load_affinities == 'save') {
					load_affinities = 2;
				} else {
					load_affinities = 0;
				}
			}
			
			if (ann_not_vptree){
			  knn_algo = 1;
			}else{
			  knn_algo = 2;
			}
			tX = c(t(X))

			f <- file(data_path, "wb")
			n = nrow(X);
			D = ncol(X);
			writeBin(as.integer(n), f,size= 4)
			writeBin( as.integer(D),f,size= 4)
			writeBin( as.numeric(theta), f,size= 8) #theta
			writeBin( as.numeric(perplexity), f,size= 8) #theta

		    if (perplexity == 0) {
		    	writeBin( as.integer(length(perplexity_list)), f, size=4)
			    writeBin( perplexity_list, f) 
		    }

			writeBin( as.integer(dims), f,size=4) #theta
			writeBin( as.integer(max_iter),f,size=4)
			writeBin( as.integer(stop_early_exag_iter),f,size=4)
			writeBin( as.integer(mom_switch_iter),f,size=4)
			writeBin( as.numeric(momentum),f,size=8)
			writeBin( as.numeric(final_momentum),f,size=8)
			writeBin( as.numeric(learning_rate),f,size=8)
			writeBin( as.integer(K),f,size=4) #K
			writeBin( as.numeric(sigma), f,size=8) #sigma
			writeBin( as.integer(nbody_algo), f,size=4)  #not barnes hut
			writeBin( as.integer(knn_algo), f,size=4) 
			writeBin( as.numeric(exaggeration_factor), f,size=8) #compexag
			writeBin( as.integer(no_momentum_during_exag), f,size=4) 
			writeBin( as.integer(n_trees), f,size=4) 
			writeBin( as.integer(search_k), f,size=4) 
			writeBin( as.integer(start_late_exag_iter), f,size=4) 
			writeBin( as.numeric(late_exag_coeff), f,size=8) 
			
			writeBin( as.integer(nterms), f,size=4) 
			writeBin( as.numeric(intervals_per_integer), f,size=8) 
			writeBin( as.integer(min_num_intervals), f,size=4) 
			tX = c(t(X))
			writeBin( tX, f) 
			writeBin( as.integer(rand_seed), f,size=4) 
		        writeBin(as.numeric(df), f, size=8)
			writeBin( as.integer(load_affinities), f,size=4) 
			if (! is.null(initialization)){ writeBin( c(t(initialization)), f) }		
		        print(df)
			close(f) 

			flag= system2(command=fast_tsne_path, args=c(version_number,data_path, result_path, nthreads));
			if (flag != 0) {
				stop('tsne call failed');
			}
			f <- file(result_path, "rb")
			n <- readBin(f, integer(), n=1, size=4);
			d <- readBin(f, integer(), n=1, size=4);
			Y <- readBin(f, numeric(), n=n*d);
		        Y <- t(matrix(Y, nrow=d));
		        if (get_costs ) {
		            tmp <- readBin(f, integer(), n=1, size=4);
		            costs <- readBin(f, numeric(), n=max_iter,size=8);
		            Yout <- list( Y=Y, costs=costs);
		        }else{
		            Yout <- Y;
		        }
		        close(f)
		        file.remove(data_path)
		        file.remove(result_path)
		        return(Yout)
		}

		textHalo <- function(
			x, y=NULL, 
			labels, 
			col='white', 
			bg='black', 
			r=0.1,
			... 
		){

			theta= seq(0, 2*pi, length.out=50);
		    xy <- xy.coords(x,y)
		    xo <- r*strwidth('A')
		    yo <- r*strheight('A')

		    # draw background text with small shift in x and y in background colour
		    for (i in theta) {
		        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
		    }
		    # draw actual text in exact xy position in foreground colour
		    text(xy$x, xy$y, labels, col=col, ... )
		}
		createColorPanel <- function(num.color){
			colPanel = c(
				"grey", "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
				"#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744", 
				"#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
			    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
			    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
			    "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
			    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
			    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
			    "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
			    "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
			    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
			    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
			    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
			   )
			if(num.color > length(colPanel)){
				colPanel = c(colPanel, colVector(num.color - length(colPanel)));
			}else{
				colPanel = colPanel[1:num.color];
			}
			return(colPanel)
		}
		colVector <- function(num.color=60, type=c("qual", "div", "seq")){
			type <- match.arg(type);
			set.seed(10)
			qual_col_pals = brewer.pal.info[brewer.pal.info$category == type,];
			set.seed(10)
			col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))));
			return(col_vector)
		} 

		plotViz <- function(obj, 
			method, 
			point.size, 
			point.shape, 
			point.alpha, 
			point.color, 
			text.add, 
			text.size, 
			text.color, 
			text.halo.add, 
			text.halo.color, 
			text.halo.width, 
			legend.add, 
			legend.pos, 
			legend.text.size,
			legend.text.color,
			down.sample, 
			pdf.file.name, 
			pdf.width, 
			pdf.height, 
			...
		){
		  UseMethod("plotViz", obj);
		}
		plotViz.default <- function(obj, 
				method=c("tsne", "umap"), 
				point.size=1, 
				point.shape=19, 
				point.alpha=0.8, 
				point.color=NULL,
				text.add=TRUE,
				text.size=1, 
				text.color="black",
				text.halo.add=TRUE,
				text.halo.color="white",
				text.halo.width=0.2,
				legend.add=FALSE,
				legend.pos=c("bottomleft", "bottom", "left", "topleft", "top", "topright", "right", "center"),
				legend.text.size=1,
				legend.text.color="black",
				down.sample=10000,
				pdf.file.name=NULL,
				pdf.width=7, 
				pdf.height=7,
				...
		){	
			if(missing(obj)){
				stop("obj is missing");
			}else{
				if(!is(obj, "snap")){
					stop("obj is not a snap object");
				}
				ncell = nrow(obj);
			}
				
			if(is.integer(down.sample)){
				stop("down.sample must be an integer");
			}
			
			method = match.arg(method);
			data.use = as.data.frame(slot(obj, method));
			if(method=="tsne"){
				colnames(data.use) = c("TSNE-1", "TSNE-2")
			}else{
				colnames(data.use) = c("UMAP-1", "UMAP-2")
			}
			if((x=nrow(data.use)) == 0L){
				stop("visulization method does not exist, run runViz first!")
			}
			xlims = c(-max(abs(data.use[,1])) * 1.05, max(abs(data.use[,1])) * 1.2);
			ylims = c(-max(abs(data.use[,2])) * 1.05, max(abs(data.use[,2])) * 1.05);
			
			point.color = obj@cluster
			cluster = point.color;
			if(((x=length(cluster)) == 0L) | (is.null(point.color))){
				warning("cluster does not exist, text.add is ignored")
				text.add = FALSE;
			}
			
			if(length(cluster) != 0L){	
				data.use$col = factor(cluster);
			}else{
				data.use$col = factor(1);
			}
			
			if(!is.null(pdf.file.name)){
				if(file.exists(pdf.file.name)){
					warning("pdf.file already exists");
					file.remove(pdf.file.name);
				}else{
					if(!file.create(pdf.file.name)){
						stop("cannot create pdf.file, not a directory")				
					}
					file.remove(pdf.file.name);
				}	
				pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
			}
			
			legend.pos = match.arg(legend.pos);
			down.sample = min(down.sample, ncell);
			idx.ds = sort(sample(seq(ncell), down.sample));
			data.use = data.use[idx.ds,,drop=FALSE]
												
			colPanel = createColorPanel(length(unique(data.use$col)));
			graphics::plot(
						   data.use[,c(1,2)],
				 		   cex=point.size, 
				 		   pch=point.shape, 
				 		   col=scales::alpha(colPanel[factor(data.use$col)], point.alpha),
						   bty="l",
						   font.lab=2,
						   col.axis = 'darkgrey',
						   xlim=xlims,
						   ylim=ylims,
				 		   ...
						   );
		  	box(bty="l", lwd=2)
			if(text.add){
				xx = findCentrod(data.use[,c(1,2)], data.use$col);
				textHalo(x=xx[,1], y=xx[,2], labels = xx[,3], col=text.color, bg=text.halo.color, r=text.halo.width, cex=text.size);
		  	}
			
			if(legend.add){
				legend.ncol = as.integer(length(levels(data.use$col)) / 25) +1;
				legend(
				  "topright", 
				  legend = levels(factor(data.use$col)), 
				  col = colPanel,
				  pch = point.shape,
				  pt.cex=1, 
				  bty = "n", 
				  cex = legend.text.size, 
				  text.col = legend.text.color,
				  horiz = FALSE,
				  ncol=legend.ncol
				  )		
			}

			if(!is.null(pdf.file.name)){
				dev.off()		
			}
		}

		plotDimReductElbow <- function(obj, point.size, point.shape, point.color, point.alpha, pdf.file.name, pdf.height, pdf.width, ...){
		  UseMethod("plotDimReductElbow", obj);
		}
		plotDimReductElbow.default <- function(
			obj, 
			point.size=1.5,
			point.shape=19,
			point.color="red",
			point.alpha=1,
			pdf.file.name=NULL,
			pdf.height=7,
			pdf.width=7,
			...
		){
			if(missing(obj)){
				stop("obj is missing");
			}else{
				if(!is(obj, "snap")){
					stop("obj is not a snap object");
				}
				ncell = nrow(obj);
				if(!isDimReductComplete(obj@smat)){
					stop("obj does not have valid dim.reduct object, run 'runDimReduct' first");			
				}
			}

			if(!is.null(pdf.file.name)){
				if(file.exists(pdf.file.name)){
					warning("pdf.file already exists");
					file.remove(pdf.file.name);
				}else{
					if(!file.create(pdf.file.name)){
						stop("cannot create pdf.file, not a directory")				
					}
					file.remove(pdf.file.name);
				}	
				pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
			}
			
			data.use = data.frame(PC=1:length(obj@smat@sdev), sd=obj@smat@sdev);	
			plot(x=data.use[,1], y=data.use[,2], cex=point.size, pch=point.shape, col=alpha(point.color, point.alpha), xlab="PCs", ylab="Standard Deviation of PCs");
				
			if(!is.null(pdf.file.name)){
				dev.off()
			}
		}

		plotDimReductPW <- function(obj, eigs.dims, point.size, point.color, point.shape, point.alpha, down.sample, pdf.file.name, pdf.height, pdf.width){
		  UseMethod("plotDimReductPW", obj);
		}
		plotDimReductPW.default <- function(
			obj, 
			eigs.dims=1:50,
			point.size=0.5,
			point.color="grey",
			point.shape=19,
			point.alpha=0.5,
			down.sample=3000,
			pdf.file.name=NULL, 
			pdf.height=7, 
			pdf.width=7
		){
			
			if(missing(obj)){
				stop("obj is missing");
			}else{
				if(!is(obj, "snap")){
					stop("obj is not a snap object");
				}
				ncell = nrow(obj);
				if(!isDimReductComplete(obj@smat)){
					stop("dim.reduct is not complete, run 'runDimReduct' first")
				}
			}
				
			down.sample = min(down.sample, ncell);
			idx.ds = sort(sample(seq(ncell), down.sample));
			obj = obj[idx.ds,,drop=FALSE];

			if(max(eigs.dims) > length(obj@smat@sdev)){
				stop(paste("eigs.dims exceeds PCA dimentions ", length(obj@smat@sdev)));
			}
			
			if((x=length(eigs.dims)) > 50L){
				stop("eigs.dims must be within 1:50")
			}
			
			
			if(!is.null(pdf.file.name)){
				if(file.exists(pdf.file.name)){
					warning("pdf.file already exists");
					file.remove(pdf.file.name);
				}else{
					if(!file.create(pdf.file.name)){
						stop("cannot create pdf.file, not a directory")				
					}
					file.remove(pdf.file.name);
				}	
				pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
			}
			
			op <- par(mfrow = c(5,5), oma = c(3,3,1,1) + 0.2, mar = c(0,0,1,1) + 0.2);
			PCA.plot <- split(sort(eigs.dims), ceiling(seq(eigs.dims)/2));
			if((length(x = eigs.dims)  %% 2) == 1){
				PCA.plot = PCA.plot[1:(length(PCA.plot) - 1)]
			}
			
			for(x in PCA.plot){
				data.use = data.frame(obj@smat@dmat[,c(x[1],x[2])]);	
				colnames(data.use) = c("dim1", "dim2");
				plot(x=data.use[,1], 
					 y=data.use[,2],
					 cex=point.size, 
					 col=scales::alpha(point.color, point.alpha),
					 mtext(paste(paste("eigs", x[1]), x[2], sep=" vs "), side=3),
					 yaxt='n', 
					 xaxt="n",
					 xlab="", 
					 ylab=""
				);
			}

			if(!is.null(pdf.file.name)){
				dev.off()		
			}	
			graphics::par(mfrow=c(1,1));
		}
		runKNN <- function(obj, eigs.dims, weight.by.lambda, k, nn.eps, save.knn, filename, snn, snn.prune) {
		  UseMethod("runKNN", obj);
		}
		runKNN.default <- function(
		  obj,
		  eigs.dims,
		  weight.by.lambda = FALSE,
		  k = 15,
		  nn.eps = 0,
		  save.knn = FALSE,
		  filename = NULL,
		  snn = FALSE,
		  snn.prune = 1/15
		){
			
			cat("Epoch: checking input parameters\n", file = stderr())
			if(missing(obj)){
				stop("obj is missing")
			}else{
				if(!is(obj, "snap")){
					stop("obj is not a snap obj")
				}
			}
			
			if(!(isDimReductComplete(obj@smat))){
				stop("dimentionality reduction is not complete, run 'runDimReduct' first")
			}
			
			ncell = nrow(obj);
			nvar = dimReductDim(obj@smat);
			
			if(missing(eigs.dims)){
				stop("eigs.dims is missing")
			}else{
				if(is.null(eigs.dims)){
					eigs.dims=1:nvar;	
				}else{
					if(any(eigs.dims > nvar) ){
						stop("'eigs.dims' exceeds PCA dimentions number");
					}		
				}
			}
			
			if(save.knn){
				if(is.null(filename)){
					stop("save.knn is TRUE but filename is NULL")
				}else{			
					if(!file.create(filename)){
						stop("fail to create filename")
					}			
				}
			}
			
			if(!is.logical(weight.by.lambda)){
				stop("weight.by.lambda must be a logical variable")
			}
			
			data.use = weightDimReduct(obj@smat, eigs.dims, weight.by.lambda);
			
		    if (ncell < k) {
		      warning("k set larger than number of cells. Setting k to number of cells - 1.")
		      k <- ncell - 1
		    }
			
			if(is.na(as.integer(k))){
				stop("k must be an integer")
			}else{
				if(k < 10 || k > 50){
					warning("too small or too large k, recommend to set k within range [10 - 50]")
				}
			}
			
			cat("Epoch: computing nearest neighbor graph\n", file = stderr())
			
			# exclude self neibours
		    nn.ranked <- nn2(
		        data = data.use,
		        k = k,
		        searchtype = 'standard',
		        eps = nn.eps)$nn.idx;

			j <- as.numeric(x = t(x = nn.ranked))
			i <- ((1:length(x = j)) - 1) %/% k + 1	
			edgeList = data.frame(i, j, 1);
				
			if(snn){
				cat("Epoch: converting knn graph into snn graph\n", file = stderr())	
				g = graph_from_edgelist(as.matrix(edgeList[,c(1,2)]), directed=FALSE);
				adj = as(similarity(g), "sparseMatrix");
				i = adj@i+1;
				j = findInterval(seq(adj@x)-1,adj@p[-1])+1;
				w = adj@x;
				idx = which(w >= snn.prune);
				edgeList = data.frame(i[idx], j[idx], w[idx]);
			}
			
			if(save.knn){
				cat("Epoch: writing resulting graph into a file\n", file = stderr())
				writeEdgeListToFile(edgeList, filename);
				obj@graph = newKgraph(file=filename, k=k, snn=snn, snn.prune=snn.prune);
			}else{
				kmat = Matrix(0, ncell, ncell, sparse=TRUE);
				kmat = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3]);
				obj@graph = newKgraph(mat=kmat, k=k, snn=snn, snn.prune=snn.prune);
			}
			gc();
			return(obj);
		} 

		runCluster <- function(obj, tmp.folder, louvain.lib, resolution, seed.use, ...) {
		  UseMethod("runCluster", obj);
		}
		runCluster.default <- function(
			obj, 
			tmp.folder, 
			louvain.lib=c("R-igraph", "leiden"),
			resolution=1.0,
			seed.use=131,
			...
		){
			cat("Epoch: checking input parameters\n", file = stderr())
			tmp.folder="./tmp_norm/"
			if(!dir.exists(tmp.folder)){
				system(paste("mkdir ", tmp.folder, sep=""))
			}
			if(missing(obj)){
				stop("obj is missing");
			}else{
				if(!is.snap(obj)){
					stop("obj is not a snap object");
				}		
			}
			
			if(missing(tmp.folder)){
				stop("tmp.folder is missing")
			}else{
				if(!dir.exists(tmp.folder)){
					stop("tmp.folder does not exist");			
				}
			}
			
			louvain.lib = match.arg(louvain.lib);
			
			if(louvain.lib == "leiden"){
				if (!requireNamespace("leiden", quietly = TRUE)) {
				      stop("Please install leiden - learn more at https://github.com/TomKellyGenetics/leiden")
				}
			}
			
			if(isKgraphEmpty(obj@graph)){
				stop("knn graph is empty, run 'runKNN' first")		
			}
						
			if(is.na(as.numeric(resolution))){
				stop("resolution must be numeric class!")
			}
			
			if(louvain.lib == "R-igraph"){
				data.use = getGraph(obj@graph);
				data.use = data.use + t(data.use);
				g = graph_from_adjacency_matrix(data.use, weighted=TRUE, mode="undirected");
				cat("Epoch: finding clusters using R-igraph\n", file = stderr())		
				set.seed(seed.use);
				cl = cluster_louvain(g);
				obj@cluster = factor(cl$membership);		
			}else if(louvain.lib == "leiden"){
				cat("Epoch: finding clusters using leiden\n", file = stderr())		
				data.use = getGraph(obj@graph);		
				set.seed(seed.use);
				obj@cluster <- factor(leiden(data.use, resolution_parameter=resolution, ...));	
			}else{
				stop("unrecognized louvain.lib option")
			}
			gc();
			return(obj);
		}
	}
}

print("Loading Paired-seq/tag modules.")
####
#### custom  functions
####

#### Loading matrices
{
	readRawMatrix <- function(path){
		data<-readMM(paste(path,"matrix.mtx",sep=""))
		pid<-read.csv(paste(path,"genes.tsv",sep=""), sep="\t", head=F)
		cid<-read.csv(paste(path, "barcodes.tsv", sep=""), sep="\t", head=F)
		colnames(data)<-cid[,1]
		rownames(data)<-pid[,1]
		return(data)
	}

	filterMatrix <-function(cutoff, data){
		nReads<-colSums(data)
		data_f<-data[,nReads>cutoff]
		nCells<-rowSums(data_f)
		data_f<-data_f[nCells>3,]
		return(data_f)
	}

	writeMatrix <- function(path, data){
		system(paste("mkdir", path, sep=" "))
		writeMM(data, paste(path, "matrix.mtx", sep="/"))
		write.table(colnames(data), col.names=F, row.names=F, file=paste(path, "barcodes.tsv", sep="/"), quote=F)
		write.table(rownames(data), col.names=F, row.names=F, file=paste(path, "genes.tsv", sep="/"), quote=F)
		cmd<-paste("/projects/ps-renlab/chz272/scripts/scripts/tr_genestsv.sh ", path, sep="")
		system(cmd)
	}

	createSnapFromPairedDNA <- function(path){
		dna_matrix<-t(readMM(paste(path,"matrix.mtx", sep="/")))
		bc<-read.csv(paste(path,"barcodes.tsv", sep="/"), head=F)
		rownames(dna_matrix)<-bc[,1]
		peaks=read.table(paste(path,"peaks.bed", sep="/"))
		peaks.gr = GRanges(peaks[,1], IRanges(peaks[,2], peaks[,3]))
		obj = createSnapFromPmat(dna_matrix, barcode=rownames(dna_matrix), peaks=peaks.gr)
		return(obj)
		gc()
	}

}

#### Matrices operation
{
	runJaccardAll <- function(obj, mat) {
	  UseMethod("runJaccardAll", obj);
	}
	runJaccardAll.default <- function(
		obj, 
		mat = c("bmat", "pmat", "gmat")
	){
		if(missing(obj)){
			stop("obj is missing")
		}else{
			if(!is(obj, "snap")){
				stop("obj is not a snap obj")
			}
		}
		mat = match.arg(mat)
		mat.use = methods::slot(obj, mat)
		p1 = Matrix::rowMeans(mat.use)
		if((x=nrow(mat.use)) == 0L){
			stop("input matrix is empty")
		}
		if((x=max(mat.use)) > 1L){
			stop("input matrix is not a binary matrix, run 'makeBinary' first")	
		}
		if(any(Matrix::rowSums(mat.use) == 0)){
			stop("input matrix contains empty rows, remove empty rows first")	
		}

		mat.ref = mat.use;				
		p2 = p1;				

		jmat = calJaccard(mat.use, mat.ref);

		rm(mat.ref);
		rm(mat.use);
		obj@jmat@jmat = jmat;
		obj@jmat@p1 = p1;
		obj@jmat@p2 = p2;
		obj@jmat@norm = FALSE;
		obj@jmat@method = character();
		return(obj);
	}

}



test_Paired_map <- function(){
	print("Finished update.")
	}
					    












print("Finished loading...")
