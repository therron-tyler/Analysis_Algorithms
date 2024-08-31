#USAGE: Rscript PCA_script.R data_filename.txt conditions_filename.csv <string value for Experimental Groups - (Legend Title)>
# SHAPES: Shape #'s 21 through 25 can be filled using specified colors within the ggplot function
# CONDITIONS FILE: 4 columns, No headers. 

# Tyler Therron 
# Winter Lab - Macrophage Genomics
# Northwestern University - Feinberg School of Medicine, Department of Rheumatology

library(ggplot2)
library(dplyr)
library(rlang)

## Read command line arguments
args <- commandArgs(TRUE)
file_name <- args[1]
conds_file <- read.csv(args[2], header=FALSE)
sample_groups <- as.character(args[3])

# Read counts data
data <- read.table(file_name, header=TRUE, sep="\t", row.names = 1, check.names = FALSE)
data$Gene_Symbol <- NULL

## Data cleaning and verification with samples in conditions file
namelist <- conds_file$V1
print(namelist)

print(colnames(data))

data <- data[ ,match(namelist, colnames(data))]


head(data)


# Define group conditions, print user made experimental conditions, colors, and shapes to the console
conds <- conds_file$V2
for(i in 1:length(conds)) {
  print(paste0("experimental conditions [", i, "/", nrow(conds_file), "] -> ", conds[i]))
}
color_mapping <- conds_file$V3
for(i in 1:length(color_mapping)) {
  print(paste0("Colors [", i, "/", nrow(conds_file), "] -> ", color_mapping[i]))
}

shape_mapping <- conds_file$V4
for(i in 1:length(color_mapping)) {
  print(paste0("Shapes [", i, "/", nrow(conds_file), "] -> ", shape_mapping[i]))
}
## Filter out zero genes
data <- data.matrix(data)

rs <- rowSums(data)
use <- (rs > 0)
data <- data[use, ]

## Define file prefix
file_prefix = strsplit(file_name,"[.]")


## Transpose the data and convert to matrix
data_tr <- t(data)
rownames(data_tr) <- colnames(data)

### PCA calculation
pca_comp <- prcomp(data_tr, scale. = TRUE, center = TRUE) # 01/12/2023 - TylerT - setting scale to false, and making center True for visual ease
percentVar <- pca_comp$sdev^2/sum(pca_comp$sdev^2)

## Create DataFrame with PCA loadings for plotting
pca_df <- data.frame(PC1= pca_comp$x[,1], PC2= pca_comp$x[,2], conds, sampleLab = rownames(pca_comp$x), shape_mapping, check.names = FALSE)
colnames(pca_df)[3] <- sample_groups
pca_df$shape_mapping <- as.factor(pca_df$shape_mapping)

#pca_df$combined_factor <- paste(pca_df[[sample_groups]], pca_df$shape_mapping, sep="_")

pca_df$combined_factor <- paste0(pca_df[[sample_groups]], "_", pca_df$shape_mapping)
pca_df$combined_factor <- as.factor(pca_df$combined_factor)
# Generate combined_factor
print(pca_df)
# Manually define mappings
#olor_mapping <- c("CD11clo_C8_1" = "blue", "CD11clo_CD11cC8_1" = "green", "CD11clo_IRF5CD11cC8_1" = "red", "CD11chi_C8_2" = "blue", "CD11chi_CD11cC8_2" = "green", "CD11chi_IRF5CD11cC8_2" = "red")
#shape_mapping <- c("CD11clo_C8_1" = 21, "CD11clo_CD11cC8_1" = 21, "CD11clo_IRF5CD11cC8_1" = 21, "CD11chi_C8_2" = 22, "CD11chi_CD11cC8_2" = 22, "CD11chi_IRF5CD11cC8_2" = 22)

color_mapping_plot <- setNames(color_mapping, pca_df$combined_factor)
shape_mapping_plot <- setNames(shape_mapping, pca_df$combined_factor)
print(color_mapping_plot)
print(shape_mapping_plot)

# Generate custom legend labels from experimental group column
#condition_labels <- unique(pca_df$`CReCOM-180918`)
condition_labels <- unique(pca_df[[sample_groups]])


# Mockup data (replace with your actual data)
combined_factors <- unique(pca_df$combined_factor)

# Create a named vector for the labels
label_mapping <- setNames(condition_labels, combined_factors)
label_mapping <- sapply(label_mapping, function(x) gsub("_", " ", x))

for(i in 1:length(label_mapping)) {
  print(paste0("Label Mapping [", i, "/", length(label_mapping), "] -> ", label_mapping[i]))
}

## PCA plot without labels
jpeg(paste0(file_prefix[[1]][1],"_PCAplot.png"), height= 10, width= 12, units = "in", res = 600)

# combined factor version
pca_df %>%
  ggplot(aes(x = PC1,y = PC2, label = sampleLab, fill = combined_factor,shape = combined_factor)) +
  geom_point(size=5, aes(shape=combined_factor)) +
  scale_shape_manual(values = shape_mapping_plot, labels = label_mapping, name = sample_groups) +
  scale_fill_manual(values = color_mapping_plot,labels = label_mapping, name = sample_groups) + # 12/02/2022 - TylerT manually putting in the sample numbers b/c code not work
  labs(x=paste0("PC1: ",round(percentVar[1]*100,1), "% Variance Explained"),
       y=paste0("PC2: ",round(percentVar[2]*100,1), "% Variance Explained"),
       fill=sample_groups, # Setting the custom legend title for fill
       shape=sample_groups) +
  theme_linedraw(base_size=16) + # 1/17/2023 - increase the PCA font size for generated figure
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size=15),legend.text = element_text(size = 12),
        legend.position = "right", axis.text = element_text(size = 12), axis.title = element_text(size = 18)) + # Adjust this value to change the label size in the legend
  guides(color = guide_legend(ncol = 1))
dev.off()

# ## Write CSV file
pcaData <- data.frame(pca_df, check.names = FALSE)
write.csv(pcaData ,paste0(file_prefix[[1]][1],"_PCA_Data.csv"), row.names=FALSE)


