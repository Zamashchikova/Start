options(digits=10)
options(scipen=10)
library(RODBC)
library(corrplot)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(gridExtra)
library(LICORS)
library( BBmisc )
library(Rtsne)
library(umap)
library( cluster )
library(rgl)
library(mclust)

#' Get data for clustering
#'
#' Get data of Sql Server
#' @param name_object is name of oject
#' @param query is query text
#' @param connectionStringDataMine is a connection string
#' @return list of data and y_name
#' @export
get_data_clust<-function(name_object=nULL,query,connection_string=nULL)
{
	dbhandle <- odbcDriverConnect(connection_string)
	data<-sqlQuery(dbhandle,query)
	on.exit(odbcClose(dbhandle))
	return(list(name_object=name_object,data=data))
}

#' Prepare evidences before clustering
#'
#' #' Prepare evidences before clustering
#' @param data data for clustering
#' @param name_object is name of oject
#' @return plot
#' @export
prepare_evidences<-function(data,name_object)
{
	x_index<-which(!names(data) %in% name_object)
	for (i in x_index)
	{
		if(is.factor(data[,i])==FALSE)
		{
			x_i<-log(data[,i]+abs(min(data[,i]))+1)
			print(names(data)[i])


			w<-2*IQR(x_i)*length(x_i)^(-1/3)
			if(w!=0)
			{
				z2<-cut(x_i, breaks=pretty(x_i,n=(ceiling(max(x_i)/w))), include.lowest=TRUE,labels=F,right = F)
				data[,i]<-z2
			}
			else
			{
				z3<-cut(x_i,breaks=1+log2(length(x_i)), labels = FALSE,right=TRUE,include.lowest=TRUE)
				data[,i]<-z3
			}
		}
		else
		{
			print(past("Warning: evidences type is factor:",names(x)[i]))
		}
		data[,x_index]<-data.frame(scale(data[,x_index]))
	}
	return(data=data)
}


#' Get wss
#'
#' #' Get wss for clustering
#' @param x data for clustering
#' @param cluster is clustering result
#' @param center is center for cluster
#' @return value
#' @export
get_wss<-function(x, cluster, center)
{
	x<-data.frame(cluster=cluster,x=x)
	z<-merge(x,center,by="cluster",all.x = TRUE)
	dist<-(z[,seq(2,(ncol(x)))]-z[,seq((ncol(x)+1),ncol(z))])^2
	sum_dist<-data.frame(cluster=z$cluster,dist=(apply(dist,1,sum)))
	wss<-sum(sum_dist$dist)
	return(wss)
}

#' Get bss
#'
#' #' Get bss for clustering
#' @param x data for clustering
#' @param cluster is clustering result
#' @param center is center for cluster
#' @return value
#' @export
get_bss<-function(x,cluster,center)
{
	center_all<-t(data.frame(apply(x,2,mean)))
	center_all<-center_all[rep(1,nrow(center)),]
	sum_dist<-data.frame(cluster=center$cluster,dist=apply((center[,-1]-center_all)^2,1,sum))
	count_object<-aggregate(data.frame(n_c=cluster), list(cluster=cluster), length)
	sum_dist<-merge(sum_dist,count_object,by="cluster",all.x = TRUE)
	bss<-sum(sum_dist$dist*sum_dist$n_c)
	return(bss)
}

#' Get plot
#'
#' Plot clustering result
#' @param data data for clustering
#' @param name_object is name of oject
#' @param cluster is vector of cluster
#' @param nameTest is name of plot
#' @param nameTest is metod for plot for count evidences > 3
#' @return plot
#' @export
get_plot<-function(data,name_object,cluster,nameTest="test",metod="umap")
{
	x_index<-which(!names(data) %in% name_object)
	x<-data[,x_index]
	if(ncol(x)==2)
	{
		p<-ggplot(data.frame('cluster'=as.factor(cluster),x), aes_string(x=x[,1], y=x[,2], color='cluster')) +
			geom_point(size=1) +
			guides(colour=guide_legend(override.aes=list(size=6))) +
			xlab("") + ylab("") +
			ggtitle(nameTest) +
			theme_light(base_size=20) +
			theme(axis.text.x=element_blank(),
						axis.text.y=element_blank(),
						legend.direction = "horizontal",
						legend.position = "bottom",
						legend.box = "horizontal")
		grid.arrange(p,ncol=1, nrow=1)
	}
	if(ncol(x)==3)
	{
		plot3d(x[,1], x[,2], x[,3], xlab=names(x)[1], ylab=names(x)[2], zlab=names(x)[3],
					 type="p", size=5, lit=F, main = nameTest, col = cluster)
		legend3d("topright",pch = 16, legend =unique(cluster) , col =unique(cluster), cex=2, inset=c(0.06))
	}
	if(ncol(x)>3)
	{
		if(metod=="tsne")
		{
			xy<- as.data.frame(Rtsne(as.matrix(x), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)$Y)
		}
		if(metod=="umap")
		{
			xy<-umap(x[,-1])$layout
		}

		p<-ggplot(data.frame('cluster'=as.factor(cluster),xy), aes_string(x=xy[,1], y=xy[,2], color='cluster')) +
			geom_point(size=1) +
			guides(colour=guide_legend(override.aes=list(size=6))) +
			xlab("") + ylab("") +
			ggtitle(nameTest) +
			theme_light(base_size=20) +
			theme(axis.text.x=element_blank(),
						axis.text.y=element_blank(),
						legend.direction = "horizontal",
						legend.position = "bottom",
						legend.box = "horizontal")
		grid.arrange(p,ncol=1, nrow=1)
	}
}

#' Get cluster info
#'
#' Get cluster info
#' @param data data for clustering
#' @param name_object is name of oject
#' @param cluster is vector of cluster
#' @return data frame
#' @export
get_cluster_info<-function(data,name_object,cluster)
{
	x_index<-which(!names(data) %in% name_object)
	x<-data[,x_index]
	r<-data.frame("cluster"=cluster,x)
	g<-NULL
	g_i<-data.frame()
	for( i in 1:length(unique(r$cluster)))
	{
		g<-rbind(g,apply((r[which(r$cluster==i),-1]),2,median))
		g<-rbind(g,apply((r[which(r$cluster==i),-1]),2,mean))
		g<-rbind(g,apply((r[which(r$cluster==i),-1]),2,max))
		g<-rbind(g,apply((r[which(r$cluster==i),-1]),2,min))
		g_i<-rbind(g_i,data.frame("agr"=c("median","mean","max","min"), "cluster"=i,"count_object"=nrow(r[which(r$cluster==i),-1])))
	}
	return(cbind(g_i,g))
}

#' Plot box plot
#'
#' plot box plot clustering result
#' @param data data for clustering
#' @param name_object is name of oject
#' @param cluster is vector of cluster
#' @return box plot
#' @export
get_box_plot<-function(data,name_object,cluster)
{
	x_index<-which(!names(data) %in% name_object)
	x<-data[,x_index]
	r<-data.frame("cluster"=cluster,x)
	par(mfrow=c(3,1))
	for(i in 2:ncol(r))
	{
		y<-log(r[,i]+abs(min(r[,i]))+1)
		boxplot(y ~ r$cluster,
						xlab = 'cluster',
						ylab = names(r)[i],
						main = 'BoxPot',
						col = "coral", data = r)
	}
}

#' Get box plot
#' @param data data for clustering
#' @param name_object is name of oject
#' @param cluster is vector of clusters
#' @param center is center for clusters
#' @param k is count of clusters
#' @return list
#' @export
get_quality<-function(data,name_object,cluster,center,k)
{
	x_index<-which(!names(data) %in% name_object)
	x<-data[,x_index]
	N<-length(cluster)
	K<-k
	D<-ncol(x)
	WSS<-get_wss(x, cluster, center)
	BSS<-get_bss(x, cluster,center)
	BH<-WSS/K
	Xu<-D*log(sqrt(WSS/(D*N^2)))+log(K)
	H<-log(WSS/BSS)
	WB <-(K*WSS)/BSS
	CH <-(WSS/(N-K))/(BSS/(K-1))
	BIC<-WSS+log(N)*D*K
	return(data.frame(WSS,BH,Xu,H,WB,BIC))
}

#' Get clusters
#'
#' Get clusters
#' @param data data for clustering
#' @param name_object is name of oject
#' @param model is model for clustering
#' @param k is count of clusters
#' @return data frame
#' @export
get_cluster<-function(data, name_object, model='kmeans', k=5)
{
	print(model)
	x_index<-which(!names(data) %in% name_object)
	x<-data[,x_index]
	silhouetteAvg<-NA
	if(model=='kmeans')
	{
		fit<-kmeans(x,k,iter.max = 3000, nstart = 3,algorithm="Forgy")
		cluster<- fit$cluster
		center<-data.frame(cluster=as.numeric(row.names(fit$center)),fit$center)
	}
	if(model=='kmeanspp')
	{
		fit<-kmeanspp(x,k)
		cluster<- fit$cluster
		center<-data.frame(cluster=as.numeric(row.names(fit$center)),fit$center)
	}
	if(model=='clara')
	{
		fit<- clara(x,k=k, metric = "euclidean", samples=500,sampsize =500)
		cluster <-fit$clustering
		center<-data.frame(cluster=unique(fit$clustering),fit$medoids)
	}
	if(model=='diana')
	{
		fit<- diana(x, metric = "euclidean")
		cluster <- cutree(fit, k=k)
		center<-data.frame(cluster=unique(cluster),apply(x,2,function(y){tapply(y,list(cluster),mean)}))
	}
	if(model=='agnes')
	{
		fit<- agnes(dist(x), metric = "manhattan")
		cluster <- cutree(fit, k=k)
		center<-data.frame(cluster=unique(cluster),apply(x,2,function(y){tapply(y,list(cluster),mean)}))
	}
	if(model=="Mclust")
	{
		fit <- Mclust(x,G=k)
		cluster <-fit$classification
		center<-data.frame(cluster=unique(cluster),apply(x,2,function(y){tapply(y,list(cluster),mean)}))
	}
	return(list(cluster=cluster,
							center=center,
							model_info= list(model=model,
															 k=k)))
}

#' Get clusters
#'
#' Start clustering
#' @param data data for clustering
#' @param name_object is name of oject
#' @param model is model for clustering
#' @param k is count of clusters
#' @param prepare is logical value for prepare evidences before clustering
#' @param plot_cluster is logical value for plot cluster
#' @param cluster_info is logical value for get cluster information
#' @param box_plot is logical value for box plot cluster
#' @return list
#' @export
clustering_go<-function(data,name_object,model='Mclust',k=3,prepare=FALSE,
												plot_cluster=TRUE,cluster_info=TRUE,box_plot=FALSE)
{
	data_start<-data
	if(prepare)
	{
		data<-prepare_evidences(data,name_object)
	}
	cluster<-get_cluster(data, name_object, model, k)
	quality<-get_quality(data,name_object,cluster$cluster,cluster$center,k)
	if(plot_cluster)
	{
		get_plot(data,name_object,cluster$cluster,nameTest="test",metod="umap")
	}
	if(cluster_info)
	{
		cluster_info<-get_cluster_info(data_start,name_object,cluster$cluster)
	}
	if(box_plot)
	{
		get_box_plot(data_start,name_object,cluster$cluster)
	}
	result<-data.frame(cluster=cluster$cluster,data_start)
	return(list(result=result,
							quality=quality,
							cluster_info=cluster_info))
}

#' find k
#'
#' Get informationt about index for finding count clusters
#' @param data data for clustering
#' @param name_object is name of oject
#' @param model is model for clustering
#' @param range_k is vector of counts of clusters
#' @param prepare is logical value for prepare evidences before clustering
#' @return list
#' @export
find_k<-function(data,name_object,model='kmeans',range_k=2:30,prepare=TRUE)
{
	if(range_k %in% 1)
	{
		range_k<-range_k[-which(range_k %in% 1)]
	}
	if(prepare)
	{
		data<-prepare_evidences(data,name_object)
	}
	else
	{
		data[,which(!names(data) %in% name_object)]<-scale(data[,which(!names(data) %in% name_object)])
	}
	#get measure by k
	quality<-data.frame()
	for(k in range_k)
	{
		cluster<-get_cluster(data, name_object, model, k=k)
		quality<-rbind(quality,data.frame(k=k,get_quality(data,name_object,cluster$cluster,cluster$center,k=k)))
	}
	par(mfrow=c(2,2))
	#get diff_WSS
	cluster<-get_cluster(data, name_object, model, k=(range_k[1]-1))
	dop_wss<-get_quality(data,name_object,cluster$cluster,cluster$center,k=(range_k[1]-1))$WSS
	quality$diff_WSS<--diff(c(dop_wss,quality$WSS))/quality$WSS[1]
	#find k by min
	best_k<-apply(quality[,-1],2,function(y){quality$k[which.min(y)]})
	#find k by cumsum+plot
	breaks<-t(as.data.frame(combn(range_k[-c(1,length(range_k))], 2)))
	best_break_min<-c()
	best_break_max<-c()
	a<-range_k[1]
	b<-range_k[length(range_k)]
	for(i in 2:ncol(quality))
	{
		plot(range_k,quality[,i],type = "l",main = names(quality)[i],xlab='k',ylab='value')
		points(best_k[i-1],quality[which(quality$k==best_k[i-1]),i],col="red")
		#cumsum<-max(quality[,i])-quality[,i]
		cumsum<-cumsum(quality[,i]+max(abs(quality[,i]))+1)
		cumsum<-cumsum/max(cumsum)
		fa<-cumsum[1]
		fb<-cumsum[length(range_k)]
		square<-apply(breaks,1,function(x)
		{
			s<-(fa+cumsum[which(range_k==x[1])])*(x[1]-a)/2+
				(fb+cumsum[which(range_k==x[2])])*(b-x[2])/2+
				(cumsum[which(range_k==x[1])]+cumsum[which(range_k==x[2])])*(x[2]-x[1])/2
			return(s)
		})
		best_break<-breaks[which.max(square),]
		best_break_min<-c(best_break_min,best_break[1])
		best_break_max<-c(best_break_max,best_break[2])
		plot(range_k,cumsum,type = "l",main = names(quality)[i],xlab='k',ylab='cumulative value')
		lines(c(a,best_break,b),c(fa,cumsum[which(range_k %in% best_break)],fb),col="blue",lty=3)
		points(best_break,cumsum[which(range_k %in% best_break)],col="red")
	}
	best_k<-rbind(best_k,best_break_min)
	best_k<-rbind(best_k,best_break_max)
	return(list(best_k=best_k,quality=quality))
}
