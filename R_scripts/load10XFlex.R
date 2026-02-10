load10XFlex <- function(dataDir, channelName=NULL) {
  # load raw and filtered data
  toc <- Read10X(data.dir=file.path(dataDir, 'sample_filtered_feature_bc_matrix'))
  tod <- Read10X(data.dir=file.path(dataDir, 'sample_raw_feature_bc_matrix'))
  # make sure the same genes are used
  tod <- tod[row.names(toc),]
  
  # load corse clustering
  tgt <- file.path(dataDir,'analysis','clustering', 'gene_expression_graphclust','clusters.csv')
  clusters <- read.csv(tgt)
  mDat <- data.frame(clusters=clusters$Cluster,row.names=clusters$Barcode)
  
  # load fine clustering
  tgt <- file.path(dataDir, 'analysis', 'clustering', 'gene_expression_kmeans_10_clusters', 'clusters.csv')
  clusters <- read.csv(tgt)
  stopifnot(all(row.names(mDat) == clusters$Barcode))
  mDat$clustersFine <- clusters$Cluster
  
  # load projections
  tgt <- file.path(dataDir, 'analysis', 'tsne', 'gene_expression_2_components', 'projection.csv')
  tsne <- read.csv(tgt)
  stopifnot(all(row.names(mDat) == tsne$Barcode))
  mDat$tSNE1 <- tsne$TSNE.1
  mDat$tSNE2 <- tsne$TSNE.2
  
  sc <- SoupChannel(tod=tod, toc=toc, metaData=mDat, channelName=channelName)
  sc
}
