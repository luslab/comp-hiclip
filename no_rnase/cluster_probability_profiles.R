#!/usr/bin/env Rscript

library(rematch)
library(stringr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(tidyverse)
library(broom)
library(ggpubr)
library(grid)
library(RColorBrewer)
library(cluster)
library(optparse)
library(gridExtra)
library(factoextra)
theme_set(theme_bw() +
            theme(legend.position = "top"))

# ==========
# Define functions
# ==========

std <- function(x) sd(x)/sqrt(length(x))

get_metaprofile_mean <- function(filename) {
  
  prob.df <- read.csv(filename, sep="\t")
  colnames(prob.df) <- seq(1:ncol(prob.df)) - (ncol(prob.df)+1)/2
  prob.df <- drop_na(prob.df, 0) # remove peaks sites with NAs at the xl site
  
  prob.mean <- prob.df %>% 
    summarise(across(where(is.numeric), mean))
  prob.sd <- prob.df %>% 
    summarise(across(where(is.numeric), std))
  prob.mean <- as.data.frame(t(prob.mean))
  colnames(prob.mean) <- "mean_prob"
  prob.sd <- as.data.frame(t(prob.sd))
  colnames(prob.sd) <- "std_prob"
  df <- cbind(prob.mean, prob.sd)
  df <- rownames_to_column(df, var = "pos")
  return(df)
  
}


plot_metaprofile <- function(data.df) {
  
  profile.gg <- ggplot(data.df, aes(x=as.numeric(pos), y=mean_prob, group = Sample)) +
    geom_line(aes(linetype = Sample)) +
    #geom_errorbar( aes(ymin = mean_prob-std_prob, ymax = mean_prob+std_prob),width = 0.2) +
    scale_linetype_manual(values=c("longdash", "solid"))+
    geom_vline(xintercept = 0, linetype = "dashed", color ="grey60", size = 0.5) +
    #ylim(c(0,1))+
    xlab("Distance relative to the peak start (nt)")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ylab("Mean probability") +
    labs(linetype = "Sample") +
    geom_text(data = data.df, 
              aes(x=75,y=0.60,label=paste0("n = ", peaks_count)), inherit.aes=FALSE, size=4) 
  
  return(profile.gg)
}

run_kmeans <- function(data, k) {
  
  set.seed(123)
  data.kmeans <- kmeans(data, centers = k)
  data <- augment(data.kmeans, data) #%>% arrange(.cluster)
  data <- data %>% add_count(.cluster, name = "cluster_size") %>%
    column_to_rownames(var = ".rownames")
  return(data)
  
}


plot_cluster_heatmap <- function(cluster.df, plot.title, plot.name) {
  
  cluster.df <- cluster.df %>% arrange(.cluster)
  clust.mx = as.matrix(dplyr::select(cluster.df, -c(.cluster, cluster_size))) # get numeric matrix
  annotation.df <- select(cluster.df, .cluster)
  gaps_row = cumsum(unique(cluster.df$cluster_size))
  mat_colors <- list(group = col_pal)
  names(mat_colors$group) <- unique(annotation.df$.cluster)
  heatmap <- pheatmap(clust.mx,cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE,
                      gaps_row = gaps_row, fontsize_col = 3.5,annotation_row = annotation.df, border_color=FALSE,
                      annotation_colors = list(.cluster=mat_colors$group), angle_col="45",
                      main = plot.title, filename = plot.name)
  return(heatmap)
  
}

plot_heatmap <- function(data.df, plot.title, plot.name) {
  
  data.df <- data.df %>% arrange("0")
  data.mx = as.matrix(data.df) # get numeric matrix
  heatmap <- pheatmap(data.mx,cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE,
                      fontsize_col = 3.5, border_color=FALSE,
                      annotation_colors = list(.cluster=mat_colors$group), angle_col="45",
                      main = plot.title, filename = plot.name)
  return(heatmap)
  
}

# cl = dataframe containing nt position probabilities +  a "cluster" column and a "cluster_size" column
plot_cluster_mean <- function(cl, left_flank) {
  
  clust.m <- melt(cl)
  clust.mean <- dcast(clust.m, cluster ~ variable, mean)
  clust.mean.m <- melt(clust.mean, id=c("cluster","cluster_size")) # arrange the df to long format for facet_wrap plotting
  clust.mean.profiles <- ggplot(data = clust.mean.m, aes(x = as.numeric(variable)-left_flank, y = value, color=cluster)) +
    geom_line() +
    geom_vline(xintercept = 0, linetype = "dashed", color ="grey60", size = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color ="grey84", size = 0.5) +
    facet_wrap(~ cluster, scales = "free_x") +   #scales = "free_y"
    xlab("Distance relative to the peak start (nt)") +
    ylab("Mean probability")+
    scale_color_manual(values = col_pal)+
    theme(text = element_text(size=14),
          strip.text = element_text(size=10, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle=60, hjust=1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # rename legend title and add the number of observations in each cluster:
  clust.mean.profiles <- clust.mean.profiles +guides(color=guide_legend(title='Cluster')) + geom_text(data=clust.mean, 
                                                                                                      aes(x=20,y=0.74,label=paste0("n = ",cluster_size)), inherit.aes=FALSE, size=3)
  return(clust.mean.profiles)
}

# ==========
# Define options and params
# ==========

option_list <- list(make_option(c("-p", "--profile"), action = "store", type = "character", default=NA, help = "tab-separated file containing nucleotide positions as columns and IDs as rows"),
                    make_option(c("-s", "--shuffled"), action = "store", type = "character", default=NULL, help = "tab-separated file containing nucleotide positions as columns and IDs as rows"),
                    make_option(c("-c", "--clusters"), action = "store", type = "integer", default = 5, help = "Number of kmeans clusters to ask for [default: %default]"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

prob.file <- opt$prob

if (opt$shuffled) {
  shuff.file <- opt$shuffled
}

# prob.file <- "stau1_threeutrs.rnaplfold_prob.df.txt"
# shuff.file <- "stau1_threeutrs.rnaplfold.shuffled_prob.df.txt"

prefix <- str_split(prob.file, pattern = ".df")[[1]][1]
prob.name <- str_to_upper(str_split(prob.file, pattern = "_")[[1]][1]) #RBP name


# ==========
# Load  metaprofile dataframes (df output from get_structure_metaprofile.R)
# ==========

prob.df <- read.csv(prob.file, sep="\t")
colnames(prob.df) <- seq(1:ncol(prob.df)) - (ncol(prob.df)+1)/2
prob.df <- drop_na(prob.df, 0)

# draw heatmap before clustering
plot_heatmap(prob.df, plot.title =prob.name, plot.name = paste0(prefix,"_heatmap.pdf"))

# calculate the mean probability and standard error of the mean
prob.mean.df <- get_metaprofile_mean(prob.file)
prob.mean.df$Sample <- prob.name

if (opt$shuffled) {
  shuff.file <- opt$shuffled
  shuff.mean.df <- get_metaprofile_mean(shuff.file)
  shuff.mean.df$Sample <- "Shuffled control"
  data.df <- rbind(prob.mean.df, shuff.mean.df)
  data.df$peaks_count <- nrow(prob.mean.df)
} else {
  data.df <- prob.mean.df
}

# plot the mean probability and standard error of the mean
profile.gg <- plot_metaprofile(data.df)
profile.gg <- profile.gg+geom_ribbon(aes(ymin=(data.df$mean_prob-data.df$std_prob), ymax=(data.df$mean_prob+data.df$std_prob)), linetype=2, alpha=0.3)
ggsave(paste0(prefix,"_metaprofile.pdf"), profile.gg)


# ==========
# K-means clustering
# ==========

# Focus on the -50 to +75 nt relative to peak starts:
prob.df <- prob.df %>% dplyr::select(51:176)
# Focus on the +10 to +75 nt relative to peak starts:
prob_downstream.df <- prob.df %>% dplyr::select(61:126)


# K-means clustering - "euclidean" dist, 5 clusters
set.seed(123)
kmeans.df <- run_kmeans(prob_downstream.df, opt$clusters)

col_pal <- brewer.pal(5, "Dark2")
plot_cluster_heatmap(kmeans.df, plot.title = prob.name, plot.name = paste0(prefix,"_kmeans.pdf"))

# join cluster information to the data containing the -50: +75 nt positions 
stopifnot(rownames(kmeans.df) == rownames(prob.df))
prob.df$cluster <- kmeans.df$.cluster # match cluster assignment to the prob.df
prob.df$cluster_size <- kmeans.df$cluster_size

clusters.gg <- plot_cluster_mean(prob.df, 51)
ggsave(paste0(prefix,"_cluster_profiles.pdf"), clusters.gg)

prob.df <- prob.df %>%
  rename(.cluster = cluster)
plot_cluster_heatmap(prob.df, plot.title = prob.name, plot.name = paste0(prefix,"_kmeans_minus50.pdf"))

# Export clusters data
prob.df <- rownames_to_column(prob.df, var = "id")
write.table(prob.df, paste0(prefix,"_clusters.df.txt"), quote = F, row.names = F, sep = "\t")


