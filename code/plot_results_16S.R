library(ggplot2)
library(argparser)
library(reshape2)
library(vegan) 
options(stringsAsFactors=F)

#This R script will:
#1. Normalize (rel abd) / rarefy user provided abd table
#2. Parse metadata file and generate a long format abd matrix for plotting
#3. Quantify alpha and beta-diversity (within and between groups)
#4. Perform ordination (nmds, prin comp, prin coord)
#5. Perform kw/wilcox or anova/t-test across group labels
#6. Perform posthoc tests for kw/anova 
#7. Optionally perform indicator species analysis
#8. Optionally perform random forests
#9. Optionally calculate network statistics
#10. Optionally graph alpha and beta diversity
#11. Optionally graph rel abd by group for significant hypothesis tests
#12. Optionally plot a barplot of top 10 OTUs/seq vars
#13. Optionally perform 1:6 for phylotypes
#14. Optional heatmap of OTU/seq var abundance


#------------------------------------------------#
#             Function make_dict                 #
#------------------------------------------------#

#This function makes a lookup table which maps IDs
#generated as in adjust_colnames with the sequence variant names  

 
make_dict <- function(df){
	x <- cbind(id = paste0(rep("seq", times = ncol(df)), c(1:ncol(df))),
		sequence = colnames(df)
		)
	return(x)
}

#------------------------------------------------#
#          Function adjust_colnames              #
#------------------------------------------------#

#This function makes substitute colnames for data frames 
#with excessively long colnames (dada2) 


adjust_colnames <- function(df){
	x <- paste0(rep("seq", times = ncol(df)), c(1:ncol(df)))
	colnames(df) <- x
	return(df)
}

#------------------------------------------------#
#              Function filter_df                #
#------------------------------------------------#

filter <- function(df) {
	df <- df[which(rowSums(df) > 0), which(colSums(df) > 0)]
	return(df)
} 

#------------------------------------------------#
#             Function normalize                 #
#------------------------------------------------#

#normalize counts either by relative abundance of rarefying
#depends on vegan
 
normalize <- function(df, method="rel", depth=depth){
	#default method = relative abundance
	if(method == "rare"){
		if( is.null(depth)){
			depth <- min(rowSums(df))
		}
		ndf <- rrarefy(df, depth)
		ndf <- ndf[,which(colSums(ndf) > 0), drop =F] 			
	}else{ 
		ndf <- sweep(df,
            	1,
                rowSums(df),
                `/`)
      #  ndf <- ndf[,which(colSums(ndf) > 0), drop =F] 
	}
return(ndf)
} 

#------------------------------------------------#
#        Function diversity_analysis             #
#------------------------------------------------#

#this function performs alpha and beta diversity analysis on the 
#user provided data frame

#depends on vegan

diversity_analysis <- function(obj){
	obj$shannon <- diversity(obj$data, index = "shannon", MARGIN = 1) 	
	obj$simpson <- diversity(obj$data, index = "simpson", MARGIN = 1)
	obj$invsimpson <- diversity(obj$data, index = "invsimpson", MARGIN = 1)
	return(obj)
}

#------------------------------------------------#
#              Function ordinate                 #
#------------------------------------------------#

#This function creates ordination objects for plotting later


ordinate <- function(obj){
	obj_mmds <- metaMDS(obj$data,
							k =5,
							distance = "bray"
							)
	
	
	obj_prcomp <- prcomp(obj$data,
							scale =T,
							center = T
							)
	
	obj_mmds <- as.data.frame(obj_mmds$points)
	obj_mmds <- obj_mmds[rownames(obj$meta),]
	obj$mds <- obj_mmds
	
	obj_prcomp <- as.data.frame(obj_prcomp$x[,1:5])
	obj_prcomp <- obj_prcomp[rownames(obj$meta),]
	obj$prcomp <- obj_prcomp

	return(obj)	
}


#------------------------------------------------#
#            Function MultiWilcox                #
#------------------------------------------------#

#performs wilcox test across a dataframe and returns pvals

MultiWilcox <- function(df, group_x, group_y) {
  #df: dataframe of OTUS, phylotype etc. with otus as rows and samples as columns
  #group_x: vector of names of the members of group x
  #group_y: vector of names of the members of group y
  #Value: a dataframe of pvalues
  temp.vec <- NULL
  for(i in 1:length(rownames(df))){
    temp1 <- wilcox.test(x= as.numeric(df[i,group_x]), y= as.numeric(df[i,group_y] ))
    pval <- temp1$p.value
    temp.vec <- c(temp.vec, pval)
  }
  temp.df <- data.frame(row.names= rownames(df), pval= temp.vec)
  return(temp.df)
}

#------------------------------------------------#
#            Function MultiKruskal               #
#------------------------------------------------#

#performs kruskal.test and returns pvalues across a vector

MultiKruskal <- function(x, grouping) {
  #x: numeric vector of values
  #grouping: a factor or numeric vector indicating which group each value belongs to
  #Value: a p-value across groups
  temp.k <- kruskal.test(x=x, g=grouping)
  return(temp.k$p.value)
}

#------------------------------------------------#
#           Function KruskalPostHoc              #
#------------------------------------------------#

#This function computes posthoc tests appropriate for determining the groups that 
#significantly differ across a significant Kruskal test.

KruskalPostHoc <- function(x, padj='fdr', grouping) {
  #x: A vector of values for which the posthoc test will be conducted
  #padj: The adjustment method for p-values (see p.adjust) default = 'fdr'
  #grouping: a factor of values indicating grouping of the data.
  temp_test <- pairwise.wilcox.test(x=as.numeric(x), g= grouping, p.adjust.method = padj, paired=F)
  lev <- as.character(levels(grouping))
  pval <- temp_test$p.value
  my_labels <- NULL
  my_pvals <- NULL
  for( i in 1:length(lev)){
    for( j in 1:length(lev)){
      if(j > i){
        my_labels <- c(my_labels, (paste(lev[i], lev[j], sep=':')))
        my_pvals  <- c(my_pvals, pval[lev[j],lev[i]])
      }
    }
  }
  names(my_pvals) <- my_labels
  return(my_pvals)
}

#------------------------------------------------#
#             Function hypo_test                 #
#------------------------------------------------#

#for each metadata column in meta this will conduct wilcox/kw test
#if kw test then it will also conduct posthoc optionally
#parametric tests will be added at a later data

hypo_test <- function(obj){
	cats <- ncol(obj$meta) - 1
	for(i in 1:cats){
		if(length(unique(obj$meta[,i+1])) ==1){
			next
		} else if (length(unique(obj$meta[,i+1])) ==2){
			#wilcox
			w <- MultiWilcox(df = obj$data,
        	group_x = 
        		rownames(obj$meta[which(obj$meta[,i+1] == unique(obj$meta[,i+1])[1]),]),
            group_y = 
            	rownames(obj$meta[which(obj$meta[,i+1] == unique(obj$meta[,i+1])[2]),])
            )
			w <- cbind(pval = w,
				fdr = p.adjust(w$pval, method = "fdr") 
				) 
			name <- paste0("wilcox_",colnames(obj$meta)[i+1])  
			obj[[name]] <- w
			#fdr
		} else {
			#kruskal
			
			k <- apply(obj$data, 2, MultiKruskal,
                        grouping = factor(obj$meta[,i+1]))
			k <- cbind(pval = k,
                        fdr = p.adjust(k, method = 'fdr'))
			name <- paste0("kruskal_",colnames(obj$meta)[i+1]) 
			obj[[name]] <- k
			
			
			kph <- 
  				t(as.data.frame(apply(obj$data,
                      2,
                      KruskalPostHoc,
                      grouping = factor(obj$meta[,i+1]),
                      padj = 'holm')
                      )
                )
            name <- paste0("kruskal_ph_",colnames(obj$meta)[i+1]) 
			obj[[name]] <- kph
		
		}
	}	
	
	return(obj)	
}

#------------------------------------------------#
#             Function plotpoints                #
#------------------------------------------------#

#Function to make pretty ordination plots using ggplot

plotpoints <- function(data,
                       name,
                       dims,
                       fill,
                       shape = NULL,
                       path,
                       grad =F,
                       size = 2
){
  #data: data frame to be used
  #name: basename of the figures to be created
  #dims: number of dimensions to be plotted
  #fill: variable that ggplot will used to fill the points
  #shape: **optional** can be specified to define variable that ggplot will use
  #to determine the shape of the points
  #path: path to the output directory
  #grad: should color point be colored by gradient instead of solid colors
  #(not compatible with shape !=F)
  if(!file.exists(path)){
    dir.create(path, recursive=T)
  }
  lshape <- T
  if(is.null(shape)){
    lshape <- F
  }
  if(lshape){
    for( i in 1:dims ){
      for( j in 1:dims ){
        if(j > i){
          d1 <- colnames(data)[i]
          d2 <- colnames(data)[j]
          ordination.plot <- ggplot(data = data,
                                    aes(x = data[,d1], y = data[,d2],
                                        colour = factor(data[,fill]),
                                        shape = factor(data[,shape])
                                    ))+
            geom_point(size=size) +
            theme(text = element_text(size=20),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key = element_blank())+
            ylab(d2)+
            xlab(d1)+
            guides(colour = guide_legend(fill), shape = guide_legend(shape))
          pdf(paste(path, name ,"_", d1, "_", d2, ".pdf",sep = ""))
          print(ordination.plot)
          dev.off()
        }
      }
    }
  }else if(grad){
    for( i in 1:dims ){
      for( j in 1:dims ){
        if(j > i){
          d1 <- colnames(data)[i]
          d2 <- colnames(data)[j]
          ordination.plot <- ggplot(data = data,
                                    aes(x = data[,d1], y = data[,d2],
                                        colour = log(as.numeric(data[,fill]))
                                    ))+
            geom_point(size=size) +
            theme(text = element_text(size=20),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key = element_blank())+
            scale_colour_gradient(paste0("log(", fill, ")"))+
            ylab(d2)+
            xlab(d1)

          pdf(paste(path, name ,"_", d1, "_", d2, ".pdf",sep = ""))
          print(ordination.plot)
          dev.off()
        }
      }
    }
  }else {
    for( i in 1:dims ){
      for( j in 1:dims ){
        if(j > i){
          d1 <- colnames(data)[i]
          d2 <- colnames(data)[j]
          ordination.plot <- ggplot(data = data,
                                    aes(x = data[,d1], y = data[,d2],
                                        colour = factor(data[,fill])
                                    ))+
            geom_point(size=size) +
            theme(text = element_text(size=20),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key = element_blank())+
            ylab(d2)+
            xlab(d1)+
            guides(colour = guide_legend(fill))
          pdf(paste(path, name ,"_", d1, "_", d2, ".pdf",sep = ""))
          print(ordination.plot)
          dev.off()
        }
      }
    }
  }
}



#------------------------------------------------#
#          Function plot_bdiversity              #
#------------------------------------------------#

#this function plots nmds and pca plots 

plot_bdiversity <- function(obj, outdir){
	outdir <- outdir
	cats <- ncol(obj$meta) - 1
	for(i in 1:cats){		 
		x <- as.data.frame(obj$mds)  
		x$param <- obj$meta[,i+1]
		
		y <- as.data.frame(obj$prcomp)
		y$param <- obj$meta[,i+1]
		
		plotpoints(data = x,
			name = colnames(obj$meta)[i+1],
			dims = 5,
			fill = 'param',
    		path  = paste0(outdir, "/plots/beta_diversity/", "nmds_",colnames(obj$meta)[i+1], "/")	 		   
		)
		
		plotpoints(data = y,
			name = colnames(obj$meta)[i+1],
			dims = 5,
			fill = 'param',
    		path  = paste0(outdir, "/plots/beta_diversity/", "pca_",colnames(obj$meta)[i+1], "/")	 		   
		)
		
	}
}

#------------------------------------------------#
#          Function plot_adiversity              #
#------------------------------------------------#

#plot alpha diversity boxplots
plot_adiversity <- function(obj, outdir){
	outdir <- outdir
	cats <- ncol(obj$meta) - 1
	for(i in 1:cats){		 
		
		if(!file.exists(paste0(outdir, "/plots/alpha_div/"))){
   			 dir.create(paste0(outdir, "/plots/alpha_div/"), recursive=T)
  		}
  		
		shannon <- as.data.frame(obj$shannon)
		shannon <- shannon[rownames(obj$meta),,drop =F] 
		shannon$param <- obj$meta[,i+1]
		colnames(shannon)[1] <- "value"

		
		simpson <- as.data.frame(obj$simpson)
		simpson$param <- obj$meta[,i+1]
		simpson <- simpson[rownames(obj$meta),,drop =F]
		colnames(simpson)[1] <- "value"
		
		invsimp <- as.data.frame(obj$invsimpson)
		invsimp <- invsimp[rownames(obj$meta),,drop =F]
		invsimp$param <- obj$meta[,i+1]
		colnames(invsimp)[1] <- "value"


		pdf(paste0(outdir, "/plots/alpha_div/", "shannon_",colnames(obj$meta)[i+1], ".pdf"))
			plot <- ggplot(shannon, 
							aes(x = param,
                            	y = value,
                                fill = param))

			plot <- plot +
  				geom_boxplot() +
 				ylab("Shannon")+
  				xlab(colnames(obj$meta)[i+1])+
  				theme(text = element_text(size=20, colour = "black"),
    			    panel.grid.major = element_blank(),
    			    panel.grid.minor = element_blank(),
     			    panel.background = element_blank(),
   			        axis.line = element_line(colour = "black"),
        			axis.text = element_text(colour = "black"),
        			axis.text.x = element_blank(),
        			axis.ticks.x = element_blank(),
        			legend.key = element_blank())+
  				scale_fill_hue("",l = 40, c = 90)
  				
  				print(plot)
  		dev.off()

		pdf(paste0(outdir, "/plots/alpha_div/", "simpson_",colnames(obj$meta)[i+1],".pdf"))
			plot <- ggplot(simpson, 
							aes(x = param,
                            	y = value,
                                fill = param))

			plot <- plot +
  				geom_boxplot() +
 				ylab("simpson")+
  				xlab(colnames(obj$meta)[i+1])+
  				theme(text = element_text(size=20, colour = "black"),
    			    panel.grid.major = element_blank(),
    			    panel.grid.minor = element_blank(),
     			    panel.background = element_blank(),
   			        axis.line = element_line(colour = "black"),
        			axis.text = element_text(colour = "black"),
        			axis.text.x = element_blank(),
        			axis.ticks.x = element_blank(),
        			legend.key = element_blank())+
  				scale_fill_hue("",l = 40, c = 90)
  				
  				print(plot)
  		dev.off()

		pdf(paste0(outdir, "/plots/alpha_div/", "invsimpson_",colnames(obj$meta)[i+1], ".pdf"))
			plot <- ggplot(invsimp, 
							aes(x = param,
                            	y = value,
                                fill = param))

			plot <- plot +
  				geom_boxplot() +
 				ylab("inverse simpson")+
  				xlab(colnames(obj$meta)[i+1])+
  				theme(text = element_text(size=20, colour = "black"),
    			    panel.grid.major = element_blank(),
    			    panel.grid.minor = element_blank(),
     			    panel.background = element_blank(),
   			        axis.line = element_line(colour = "black"),
        			axis.text = element_text(colour = "black"),
        			axis.text.x = element_blank(),
        			axis.ticks.x = element_blank(),
        			legend.key = element_blank())+
  				scale_fill_hue("",l = 40, c = 90)
  				
  				print(plot)
  		dev.off()

		
	}
}

#------------------------------------------------#
#             Function write_dfs                 #
#------------------------------------------------#

write_dfs <- function(obj, outdir){
	outdir <- outdir
	if(!file.exists(paste0(outdir, "/flatfiles/"))){
    dir.create(paste0(outdir, "/flatfiles/"), recursive=T)
  	}
	
	dfs <- names(obj)
	dfs <- dfs[-which(dfs %in% c("data", "meta"))] 
	for(i in 1:length(dfs)){
	write.table(obj[[dfs[i]]], 
		paste0(outdir, "/flatfiles/", dfs[i], ".txt"),
		quote = F, sep = "\t"
		)	
	}
}


#------------------------------------------------#
#             Function run_adonis                #
#------------------------------------------------#

run_adonis <- function(obj, outdir){
	outdir <- outdir
	if(!file.exists(paste0(outdir, "/flatfiles/adonis/"))){
    dir.create(paste0(outdir, "/flatfiles/adonis/"), recursive=T)
  	}
	
	obj$data <- obj$data[rownames(obj$meta),,drop =T] 
	
	cats <- ncol(obj$meta) - 1
	for(i in 1:cats){
		x <- adonis(obj$data ~ obj$meta[,i+1], permutations = 5000)
		rsq <- x$aov.tab$R2[1]
		pval <- x$aov.tab$`Pr(>F)`[1]
		y <- cbind( r2 = rsq, 
			pval = pval
			)
		rownames(y) <- colnames(obj$meta)[i+1]
		write.table(y, 
		paste0(outdir, "/flatfiles/adonis/", colnames(obj$meta)[i+1], "_adonis.txt"),
		quote = F, sep = "\t"
		)  
	}		 
}


#------------------------------------------------#
#             Function hypo_alpha                #
#------------------------------------------------#

hypo_alpha <- function(obj,outdir){
	outdir <- outdir
	if(!file.exists(paste0(outdir, "/flatfiles/alpha_sig/"))){
    dir.create(paste0(outdir, "/flatfiles/alpha_sig/"), recursive=T)
  	}
  	apval <- NULL # pval 
  	amet  <- NULL # the diversity stat plus meta data param
  	atest <- NULL # test type (wilcoxon, kw)
  	astat <- NULL # statistic
  	
  	cats <- ncol(obj$meta) - 1
	for(i in 1:cats){
		
		shannon <- as.data.frame(obj$shannon)
		shannon <- shannon[rownames(obj$meta),,drop =F] 
		shannon$param <- obj$meta[,i+1]
		colnames(shannon)[1] <- "value"

		
		simpson <- as.data.frame(obj$simpson)
		simpson$param <- obj$meta[,i+1]
		simpson <- simpson[rownames(obj$meta),,drop =F]
		colnames(simpson)[1] <- "value"
		
		invsimp <- as.data.frame(obj$invsimpson)
		invsimp <- invsimp[rownames(obj$meta),,drop =F]
		invsimp$param <- obj$meta[,i+1]
		colnames(invsimp)[1] <- "value"
		
		alphas <-list(shannon = shannon,
			simpson =simpson, 
			invsimpson = invsimp)
	
		if(length(unique(obj$meta[,i+1])) ==1){
			next
		} else if (length(unique(obj$meta[,i+1])) ==2){
			#wilcox
			names <- unique(obj$meta[,i+1])
			for(j in 1:length(alphas)){
				df <- alphas[[j]]
				xvals  <- unlist(df[which(m[,2] == names[1]),1])
				yvals  <- unlist(df[which(m[,2] == names[2]),1]) 
				w      <- wilcox.test(x = xvals ,y = yvals)
				
				apval  <- c(apval, w$p.value)
				astat  <- c(astat, unname(w$statistic))
				atest  <- c(atest, "wilcoxon")
				amet   <- c(amet, paste0(names(alpha)[j], "_", colnames(obj$meta)[i+1])) 
				}
		} else {
			#kruskal
			names <- unique(obj$meta[,i+1])
			for(j in 1:length(alphas)){
				df     <- alphas[[j]]
				#print(colnames(df))
				#print(df[,2])
				w      <- kruskal.test(df[,1], g = as.factor(df[,2]) )
				
				apval  <- c(apval, w$p.value)
				astat  <- c(astat, unname(w$statistic))
				atest  <- c(atest, "kruskal")
				amet   <- c(amet, paste0(names(alphas)[j], "_", colnames(obj$meta)[i+1]))
			
			kph <- pairwise.wilcox.test(df[,1], df[,2], p.adj = "holm")
			
  				
            name <- paste0("kruskal_ph_",names(alphas)[j],"_", colnames(obj$meta)[i+1]) 
			write.table(kph$p.value,
			paste0(outdir, "/flatfiles/alpha_sig/", name, ".txt"),
			quote = F, sep = "\t"
			)
			
			}
		}
	}	
	#make a dataframe of findings
	alpha_div.df <- cbind(method = amet, 
		test = atest, statistic = astat, pval = apval)
	
	write.table(alpha_div.df, 
		paste0(outdir, "/flatfiles/alpha_sig/alpha_sig.txt"),
		quote = F, sep = "\t", row.names = F
		) 
}

#------------------------------------------------#
#             Function run_beta                  #
#------------------------------------------------#

#calculates beta-div between all samples

run_beta <- function(obj){
	bdiv <- as.matrix(vegdist(obj$data))
	obj$bdiv <- bdiv
	return(obj)
}

#------------------------------------------------#
#             Function beta_intra                #
#------------------------------------------------#

beta_intra <- function(obj, group){
	#group: names of samples in group
	bdiv <- obj$bdiv
	bdiv_in <- bdiv[group, group]
	bdiv_in <- bdiv_in[upper.tri(bdiv_in)]
	return(bdiv_in)
}

#------------------------------------------------#
#             Function beta_inter                #
#------------------------------------------------#

beta_inter <- function(obj,x,y){
	#x: names of samples in group 1
	#y: names of samples in group 2 
	bdiv <- obj$bdiv
	bdiv_in <- as.numeric(bdiv[x, y])
	return(bdiv_in)
}


#------------------------------------------------#
#             Function plot_beta                 #
#------------------------------------------------#

#plot boxplots of beta diversity and prints stats 
plot_beta <- function(obj, outdir){
	outdir <- outdir
	if(!file.exists(paste0(outdir, "/plots/bdiv/"))){
    dir.create(paste0(outdir, "/plots/bdiv/"), recursive=T)
  	}
  	
  	if(!file.exists(paste0(outdir, "/flatfiles/bdiv/"))){
    dir.create(paste0(outdir, "/flatfiles/bdiv/"), recursive=T)
  	}
  	
  	cats <- ncol(obj$meta) - 1
	
	for(i in 1:cats){
  		intra <- NULL
  		# get list of all rownames across all groups
  		ug <- unique(obj$meta[,i+1])
  		groups <- list()
  		
  		for(j in 1:length(ug)){
  			rn <- rownames(obj$meta[which(obj$meta[,i+1] == ug[j]),])
  			groups[[ug[j]]] <- rn
  		}
  		
  		intra <- lapply(groups, beta_intra, obj= obj)
  		inter <- list()
  		for(k in 1:length(ug)){
  			for(l in 1:length(ug)){
  				if(l > k ) {
  					x <- beta_inter(obj, unlist(groups[ug[k]]), unlist(groups[ug[l]]))
  					name <- paste0(ug[k], "_", ug[l])
  					inter[[name]] <- x
  				}
  			}	
  		}
  		
  		bdiv_intra.lens  <- sapply(intra, length)
  		bdiv_intra.names <- names(intra)
  		bdiv_intra.names <- rep(bdiv_intra.names, times = bdiv_intra.lens)
  		
  		bdiv_intra.df <- cbind(names = bdiv_intra.names,
  		 						value = unlist(intra))
  		bdiv_intra.df <- as.data.frame(bdiv_intra.df)
  		bdiv_intra.df$value <- as.numeric(bdiv_intra.df$value)
  		
  		
  		pdf(paste0(outdir, "/plots/bdiv/", "bray_intra_",colnames(obj$meta)[i+1], ".pdf"))
			plot <- ggplot(bdiv_intra.df, 
							aes(x = names,
                            	y = value,
                                fill = names))

			plot <- plot +
  				geom_boxplot() +
 				ylab("bray curtis")+
  				xlab(colnames(obj$meta)[i+1])+
  				theme(text = element_text(size=20, colour = "black"),
    			    panel.grid.major = element_blank(),
    			    panel.grid.minor = element_blank(),
     			    panel.background = element_blank(),
   			        axis.line = element_line(colour = "black"),
        			axis.text = element_text(colour = "black"),
        			axis.text.x = element_blank(),
        			axis.ticks.x = element_blank(),
        			legend.key = element_blank())+
  				scale_fill_hue("",l = 40, c = 90)
  				
  				print(plot)
  		dev.off()
  		
  		
  		bdiv_inter.lens  <- sapply(inter, length)
  		bdiv_inter.names <- names(inter)
  		bdiv_inter.names <- rep(bdiv_inter.names, times = bdiv_inter.lens)
  		
  		bdiv_inter.df <- cbind(names = bdiv_inter.names,
  		 						value = unlist(inter))
  		bdiv_inter.df <- as.data.frame(bdiv_inter.df)
  		bdiv_inter.df$value <- as.numeric(bdiv_inter.df$value)
  		
  		
  		pdf(paste0(outdir, "/plots/bdiv/", "bray_inter_",colnames(obj$meta)[i+1], ".pdf"))
			plot <- ggplot(bdiv_inter.df, 
							aes(x = names,
                            	y = value,
                                fill = names))

			plot <- plot +
  				geom_boxplot() +
 				ylab("bray curtis")+
  				xlab(colnames(obj$meta)[i+1])+
  				theme(text = element_text(size=20, colour = "black"),
    			    panel.grid.major = element_blank(),
    			    panel.grid.minor = element_blank(),
     			    panel.background = element_blank(),
   			        axis.line = element_line(colour = "black"),
        			axis.text = element_text(colour = "black"),
        			axis.text.x = element_blank(),
        			axis.ticks.x = element_blank(),
        			legend.key = element_blank())+
  				scale_fill_hue("",l = 40, c = 90)
  				
  				print(plot)
  		dev.off()
  		
  		#Test of significance 
  		if(length(unique(bdiv_intra.df[,1])) ==1){
			next
		} else if (length(unique(bdiv_intra.df[,1])) ==2){
			#wilcox
			#intra
			bdiv_intra.wilcox <- wilcox.test(
				x = unlist(bdiv_intra.df[
					which(bdiv_intra.df[,1] == unique(bdiv_intra.df[,1])[1] ),2]), 
				y = unlist(bdiv_intra.df[
					which(bdiv_intra.df[,1] == unique(bdiv_intra.df[,1])[2] ),2])
			)
			
			vec1 <- cbind( name = paste0("beta_intra", "_", colnames(obj$meta)[i+1]),
				type      = "wilcoxon",
				statistic = unname(bdiv_intra.wilcox$statistic),
				pvalue    = bdiv_intra.wilcox$p.value
			)
  		
	  		name <- paste0("wilcoxon_","beta_intra_", colnames(obj$meta)[i+1]) 
			
			write.table(vec1,
				paste0(outdir, "/flatfiles/bdiv/", name, ".txt"),
				quote = F, sep = "\t"
			)
				
			#inter 
			#there is only one group so no comparisons are needed
			
			#bdiv_inter.wilcox <- wilcox.test(x = 
			#	unlist(bdiv_inter.df[
			#		which(bdiv_inter.df[,1] == unique(bdiv_inter.df[,1])[1] ),2]), 
			#	y = unlist(bdiv_inter.df[
			#		which(bdiv_inter.df[,1] == unique(bdiv_inter.df[,1])[2] ),2])
			#)
			#
			#vec1 <- cbind( name = paste0("beta_inter", "_", colnames(obj$meta)[i+1]),
			#	type      = "wilcoxon",
			#	statistic = unname(bdiv_inter.wilcox$statistic),
			#	pvalue    = bdiv_inter.wilcox$p.value
			#)
  			#
			#name <- paste0("wilcoxon_","beta_inter_", colnames(obj$meta)[i+1]) 
			#write.table(vec1,
			#	paste0(outdir, "/flatfiles/bdiv/", name, ".txt"),
			#	quote = F, sep = "\t"
			#)
				
		} else {
			#kruskal	
			#intra
  			bdiv_intra.kw <- kruskal.test(bdiv_intra.df$value ~ factor(bdiv_intra.df$names))
  			
  			vec1 <- cbind(name = paste0("beta_intra", "_", colnames(obj$meta)[i+1]), 
  		 		statistic = unname(bdiv_intra.kw$statistic),
  		 		test = "kruskal",  
  		 		pvalue    = bdiv_intra.kw$p.value		 
			)
  			name <- paste0("kruskal_", "beta_intra_", colnames(obj$meta)[i+1]) 
  			
  			write.table(vec1,
				paste0(outdir, "/flatfiles/bdiv/", name, ".txt"),
				quote = F, sep = "\t"
			)
  			
  			#posthoc
  			bdiv_intra.kwph <- pairwise.wilcox.test(bdiv_intra.df$value,
  				factor(bdiv_intra.df$names))		
  		
	  		name <- paste0("kruskal_ph_","beta_intra_", colnames(obj$meta)[i+1]) 
			
			write.table(bdiv_intra.kwph$p.value,
				paste0(outdir, "/flatfiles/bdiv/", name, ".txt"),
				quote = F, sep = "\t"
			)
  		
  			#inter	
  			bdiv_inter.kw <- kruskal.test(bdiv_inter.df$value ~ factor(bdiv_inter.df$names))
  			
  			vec1 <- cbind(name = paste0("beta_inter", "_", colnames(obj$meta)[i+1]), 
  		 		statistic = unname(bdiv_inter.kw$statistic), 
  		 		test = "kruskal", 
  		 		pvalue    = bdiv_inter.kw$p.value		 
			)
  			name <- paste0("kruskal_", "beta_inter_", colnames(obj$meta)[i+1]) 
  			
  			write.table(vec1,
				paste0(outdir, "/flatfiles/bdiv/", name, ".txt"),
				quote = F, sep = "\t"
			)
  			
  			#posthoc
  			bdiv_inter.kwph <- pairwise.wilcox.test(bdiv_inter.df$value,
  				factor(bdiv_inter.df$names))		
  		
	  		name <- paste0("kruskal_ph_","beta_inter_", colnames(obj$meta)[i+1]) 
			
			write.table(bdiv_inter.kwph$p.value,
				paste0(outdir, "/flatfiles/bdiv/", name, ".txt"),
				quote = F, sep = "\t"
			)
		
		}
  	}
}	



#------------------------------------------------#
#            Function group_means                #
#------------------------------------------------#

#This function will calculate the group means across all columns in a
#dataframe
  

group_means <- function(df, mapping, t = F){
  #df : data frame with rows as sample IDs and columns as objects (e.g. OTUs)
  #mapping : a mapping df with rownames df as rownames and group id in col2
  #t : boolean, defaults to FALSE. Use if df has sample ids as colnames
  if(t == T){
    df <- t(df)
  }
  groups <- base::unique(x=mapping[,2])
  my_df <- data.frame(matrix(nrow = length(groups), ncol = ncol(df)))
  for(i in 1:ncol(df)){
    tgvec <- NULL
    for(j in 1:length(groups)){
      s <- base::rownames(mapping[base::which(mapping[,2] %in% groups[j]),
                                  ,drop=F])
      m <- base::mean(df[s,i])
      tgvec <- c(tgvec, m)
    }
    my_df[,i] <- tgvec
  }
  rownames(my_df) <- groups
  colnames(my_df) <- colnames(df)
  return(my_df)
}


#------------------------------------------------#
#        Function plot_seq_barplots              #
#------------------------------------------------#

plot_seq_barplots <- function(obj, outdir){
	outdir <- outdir
	if(!file.exists(paste0(outdir, "/plots/barplots/"))){
    dir.create(paste0(outdir, "/plots/barplots/"), recursive=T)
  	}
  	
  	obj$data <- obj$data[rownames(obj$meta),,drop =T] 
	top <- obj$data[,order(colMeans(obj$data), decreasing =T), drop =F] 
	top <- top[,1:10,drop=F] # get 10 most abundant things
	 
	cats <- ncol(obj$meta) - 1
	xm <- NULL
	for(i in 1:cats){		 
		x <- t(group_means(top, mapping = obj$meta[,c(1,i+1)]))
		x <- as.data.frame(x)
		x$taxa <- rownames(x)
		xm <- melt(x) 
		colnames(xm) <- c("taxa", colnames(obj$meta)[i+1], "value")
		
		#prep the plot
		plot <- ggplot(xm, aes(x = colnames(obj$meta)[i+1], 
								y = value,
								fill = taxa
								)
						)
		plot <- plot +
			geom_bar(stat = "identity") +xlab("")+
 		ylab("Taxa Abundance")+
  		theme(
  		text = element_text(size=16, colour = "black"),
  		panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line        = element_line(colour = "black"),
        legend.key = element_blank()
  		)
  		
  		#plot
  		pdf(paste0(outdir, "/plots/barplots/", "barplot_",colnames(obj$meta)[i+1], ".pdf"))
  		print(plot)
  		dev.off()
		}
}

#------------------------------------------------#
#                Function main                   #
#------------------------------------------------#
		 
main <- function(){
		
	p <- arg_parser("plot_results_16S.R")
	p <- add_argument(p, "infile", help="the file path to the input file. Otus (columns) x samples (rows) matrix.")
	p <- add_argument(p, "meta", help="the file path to the metadata file.")
	p <- add_argument(p, "outdir", help="the file path to an output directory")
	p <- add_argument(p, 
		"--method",
		help="Normalization method (relative abundance (\"rel\"default), or rarefaction (\"rare\").",
		default = "rel")
	p <- add_argument(p, 
		"--depth",
		help="Normalization depth for rarefaction).",
		default = NULL)
	p <- add_argument(p, 
		"--source",
		help="source of data (default = dada2). ",
		default = "dada2")
	p <- add_argument(p,
		"--t",
		help="should the infile be transposed?",
		default = F)
	
	args <- parse_args(p)

	trans <- args$t
	outdir <- args$outdir
	
	print("reading input files...")
	
	df <- read.table(args$infile, sep = "\t", row.names = 1, header = T)
	metadf <- read.table(args$meta, sep = "\t", header = T) 
	rownames(metadf) <- metadf[,1]
	
	#does this need to be transposed? 
	if(trans){
		df <- t(df)
	} 	
	
	
	#dada2 names sequence variants by their actual seqeunce. Not sure why. 
	if(args$source == "dada2"){
		print("making variant dictionary...")
		var_dict <- make_dict(df) #this matches 
		df <- adjust_colnames(df)
		write.table(var_dict, 
		paste0(outdir, "id_to_sequence.txt"), 
		sep ="\t",
		quote = F
		)
	}
	
	print("normalizing data...")
	#filter data 
	df <- filter(df)
	metadf <- metadf[rownames(df),,drop=F] 
	#normalize data
	ndf <- normalize(df, method = args$method, depth = args$depth)
	print("populating object...")
	
	obj <- NULL
	obj$data <- ndf
	obj$meta <- metadf # might want to write a check in for metadata
	
	
	#merge_metadata() # not sure if this is totally necessary until later maybe
	
	print("Running alpha diversity analyses...")
	obj <- diversity_analysis(obj)
	
	print("Running beta diversity analyses...")
	obj <- run_beta(obj)
	
	print("Making ordinations...")
	obj <- ordinate(obj) # need to write checks into this to make sure it gets right data type
	print("Running hypothesis tests...")
	obj <- hypo_test(obj)
	
	print("Plotting...")
	plot_bdiversity(obj, outdir = outdir)
	plot_adiversity(obj, outdir = outdir)
	plot_seq_barplots(obj, outdir = outdir)
	write_dfs(obj, outdir = outdir)
	hypo_alpha(obj, outdir = outdir)
	#hypo_beta(obj, outdir = outdir)# doesnt exist
	run_adonis(obj, outdir = outdir)
	plot_beta(obj, outdir = outdir)
	
	#optional and coming soon... maybe soonish. 
	#indic_species_analysis()
	#random_forests()
	#network_analysis()
	#make_heatmap()	
}

main()