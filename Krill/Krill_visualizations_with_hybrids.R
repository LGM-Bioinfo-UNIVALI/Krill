#######################Creating visualizations from krill output tsv file
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(gridExtra)
library(plotly)
library(igraph)
library(vegan)
library(reshape2)
library(MASS)


args <- commandArgs(trailingOnly = TRUE)

# Read the TSV file
BGCs <- read.delim(args[1], stringsAsFactors = T)

output_dir <- args[2]

# Extract "Bioproject" and "Source" from the "database" column
BGCs <- BGCs %>%
  mutate(Bioproject = sub("_(.*?)$", "", Database),
         Source = sub("^(.*?)_", "", Database))


# Convert "product" column to character
BGCs$product_bigscape <- as.character(BGCs$product_bigscape)

# Select columns
df_BGCs <- BGCs %>%
dplyr::select(Source, product_bigscape, Bioproject, Size)

# Create a histogram for each product
plots <- list()

for (product_bigscape in unique(df_BGCs$product_bigscape)) {
  # Subset data for the current product
  product_data <- df_BGCs[df_BGCs$product_bigscape == product_bigscape, ]
  
  
  
  # Create the histogram for the current product
  plot <- ggplot(product_data, aes(x = Size, fill = product_bigscape)) +
    geom_histogram(binwidth = 1000, color = "black", alpha = 0.7) +
    labs(x = "Size", y = "Frequency", title = paste("Size Distribution for", product_bigscape)) +
    theme_minimal() +
    theme(
      text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    )
  
  
  # Store the plot in the list
  plots[[product_bigscape]] <- plot
}

# Print each plot separately
for (product_bigscape in names(plots)) {
  ggsave(paste(output_dir, "/", product_bigscape, "_histogram_size_plot.png", sep=''), plots[[product_bigscape]], width = 10, height = 6, dpi = 300)
}



# Count how many times a product category occurs in each category
source_product_counts <- df_BGCs %>%
  group_by(Source, product_bigscape, Bioproject) %>%
  tally() %>%
  rename(Product_Counts = n)

color_palette <- c(
  "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
  "#FF8000", "#800080", "#008000", "#800000", "#000080", "#FFA500", "#008080",
  "#FF5733", "#33FF57", "#5733FF", "#FFFF33", "#FF33FF", "#33FFFF",
  "#FF9966", "#9966FF", "#66FF99", "#FFCC66", "#66FFCC", "#CC66FF", "cyan"
)


# Create the bar plot
bar_plot <- ggplot(source_product_counts, aes(x = reorder(product_bigscape, -Product_Counts), y = Product_Counts, fill = Source)) +
  geom_bar(stat = "identity") +
  labs(x = "Product", y = "Count", title = "Product Total Absolute Count") +
  scale_fill_manual(values = color_palette) +  # Define the color palette
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  )

# Display the bar plot and save

ggsave(paste(output_dir, "/", "product_total_absolute_count_bar_plot.png", sep=''), plot = bar_plot, width = 10, height = 6, dpi = 800)



# Calculate the total count in each source
source_total <- source_product_counts %>%
  group_by(Source) %>%
  summarise(Total_Count = sum(Product_Counts))

# Calculate the total count in each source
source_total <- aggregate(Product_Counts ~ Source, data = source_product_counts, sum)

# Merge the total count back to the original dataframe
product_counts <- merge(source_product_counts, source_total, by = "Source", suffixes = c("", "_total"))

# Calculate the percentage for each product count
product_counts$Percentage <- product_counts$Product_Counts / product_counts$Product_Counts_total * 100

# Calculate the sum of percentages for each sample
source_sums <- aggregate(Percentage ~ Source, data = product_counts, sum)

# Merge the sum of percentages back to the original dataframe
product_counts <- merge(product_counts, source_sums, by = "Source", suffixes = c("", "_sum"))

# Calculate the normalized percentage for each product count
product_counts$Normalized_Percentage <- product_counts$Percentage / product_counts$Percentage_sum * 100

# Remove the sum of percentages column if you don't need it anymore
product_counts <- subset(product_counts, select = -Percentage_sum)

# Create the stacked bar graph with facets for Source
custom_colors <- c(
  '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
  '#00ff00', '#e6beff', '#9a6324', '#fffac8', '#ff0000', '#aaffc3', '#ffff00', '#ffd8b1', '#000075', '#00ffcc', 
  '#ff7f00', '#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3', '#D55E00', '#CC79A7', '#56B4E9', '#009E73', 
  '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#9E0142', '#00FFFF', '#ff99cc', '#4E387E', '#008000', '#FFD700', 
  '#4B0082', '#8B4513', '#20B2AA', '#8A2BE2', '#5F9EA0', '#B22222', '#8B0000', '#CD853F', '#32CD32', '#4682B4', 
  '#800080', '#ADFF2F', '#4B0082', '#FF00FF', '#8A2BE2', '#6495ED', '#9932CC', '#87CEEB', '#8B4513', '#B8860B', 
  '#CD5C5C', '#48D1CC', '#20B2AA', '#6A5ACD', '#8B008B', '#F0E68C', '#D2691E', '#BA55D3', '#00FA9A', '#C71585'
)

bar_graph <- ggplot(product_counts, aes(x = Source, y = Normalized_Percentage, fill = product_bigscape)) +
  geom_bar(stat = "identity") +
  labs(x = "Source", y = "Percentage", title = "Product Percentage by Source") +
  scale_fill_manual(values = custom_colors) +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) +
  coord_cartesian(ylim = c(0, 100))   

# Print and save the graph
ggsave(paste(output_dir, "/", "product_percentage_by_source_stackedbar_plot.png", sep=''), plot = bar_graph, width = 10, height = 6, dpi = 800)


# Heatmap of product count by source
heatmap_plot <- ggplot(product_counts, aes(x = Percentage, y = Source, fill = product_bigscape)) +
  geom_tile() +
  labs(x = "Percentage", y = "Source") +
  ggtitle("Product Frequency") +
  scale_fill_manual(values = custom_colors) +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  )

#Print and save
ggsave(paste(output_dir, "/", "product_frequency_by_source_heatmap_plot.png", sep=''), plot = heatmap_plot, width = 10, height = 6, dpi = 800)

# Create a bubble plot with product size
bubble_plot_s <- ggplot(df_BGCs, aes(x = Source, y = Size, size = Size, color = product_bigscape)) +
  geom_point() +
  labs(x = "Source", y = "Size", size = "Size", color = "Product") +
  scale_fill_manual(values = custom_colors) +
  ggtitle("Product Size") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  )
#Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_bubble_plot.png", sep=''), plot = bubble_plot_s, width = 10, height = 6, dpi = 800)

# Boxplot of Product Absolute Count by Source
boxplot_plot_pc <- ggplot(product_counts, aes(x = product_bigscape, y = Product_Counts, fill = Source)) +
  geom_boxplot(width = 0.7) +  # Adjust the width as needed
  labs(x = "Product", y = "Count", title = "Product Absolute Count by Source") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)  ) +
  facet_grid(rows = vars(Source))

# Print and save
ggsave(paste(output_dir, "/", "product_absolute_count_by_source_box_plot.png", sep=''), plot = boxplot_plot_pc, width = 10, height = 6, dpi = 800)

# Create the boxplot with facet for Bioproject
boxplot_plot2 <- ggplot(product_counts, aes(x = product_bigscape, y = Product_Counts, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Count", title = "Product Absolute Count by Source and Bioproject") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  ) +
  facet_grid(rows = vars(Bioproject))  # Facet for Bioproject

# Print and save
ggsave(paste(output_dir, "/", "product_absolute_count_by_source_and_bioproject_box_plot.png", sep=''), plot = boxplot_plot2, width = 10, height = 6, dpi = 300)


# Create a Boxplot for Product Size 
boxplot_plot_ps <- ggplot(df_BGCs, aes(x = product_bigscape, y = Size, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Size", title = "Product Size by Source") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )


# Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_box_plot.png", sep=''), plot = boxplot_plot_ps, width = 10, height = 6, dpi = 800)

# Create the boxplot with facet for Bioproject
boxplot_plot2s <- ggplot(df_BGCs, aes(x = product_bigscape, y = Size, fill = Bioproject)) +
  geom_boxplot() +
  labs(x = "Product", y = "Size", title = "Product Size by Source and Bioproject") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  ) +
  facet_grid(rows = vars(Bioproject))  # Facet for Bioproject

# Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_and_bioproject_box_plot.png", sep=''), plot = boxplot_plot2s, width = 10, height = 6, dpi = 300)

# Create the boxplot for  Product Percentage by Source
boxplot_plot_p <- ggplot(product_counts, aes(x = product_bigscape, y = Percentage, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Percentage", title = "Product Percentage by Source") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white")# Rotate x-axis labels for better readability
  )

# Print and save
ggsave(paste(output_dir, "/", "product_percentage_by_source_box_plot.png", sep=''), plot = boxplot_plot_p, width = 12, height = 6, dpi = 300)

# Create the boxplot with facet for Bioproject
boxplot_plot2p <- ggplot(product_counts, aes(x = product_bigscape, y = Percentage, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Percentage", title = "Product Percentage by Source and Bioproject") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    panel.background = element_rect(fill = "white")# Rotate x-axis labels for better readability
  ) +
  facet_grid(rows = vars(Bioproject))  # Facet for Bioproject

# Print and save
ggsave(paste(output_dir, "/", "product_percentage_by_source_and_bioproject_box_plot.png", sep=''), plot = boxplot_plot2p, width = 12, height = 6, dpi = 800)

# Create the bubble plot with Product size by Source and Bioproject
bubble_plot_c <- ggplot(df_BGCs, aes(x = product_bigscape, y = Source, size = Size, color = product_bigscape)) +
  geom_point(alpha = 0.7) +
  labs(x = "Product", y = "Source", title = "Product Size by Source and Bioproject") +
  scale_size_continuous(range = c(3, 15)) +  # Set the range of bubble sizes
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) +
  facet_grid(Bioproject ~ ., scales = "free_y")  # Create facets for each Bioproject

# Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_and_bioproject_box_plot.png", sep=''), plot = boxplot_plot2s, width = 10, height = 6, dpi = 300)


# Create the bubble plot with Product Absolute Count by Source 
bubble_plot_p <- ggplot(product_counts, aes(x = product_bigscape, y = Source, size = Product_Counts, color = product_bigscape)) +
  geom_point(alpha = 0.7) +
  labs(x = "Product", y = "Source", title = "Product Absolute Count by Source") +
  scale_size_continuous(range = c(3, 15)) +  # Set the range of bubble sizes
  scale_fill_manual(values = custom_colors) +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) 
# Print and save


ggsave(paste(output_dir, "/", "product_absolute_count_by_source_bubble_plot.png", sep=''), plot = bubble_plot_p, width = 12, height = 6, dpi = 800)

# Create the heatmap for each product by Source
heatmap_plot_c <- ggplot(product_counts, aes(x = product_bigscape, y = Source, fill = Product_Counts)) +
  geom_tile(color = "white") +
  labs(x = "Product", y = "Source", title = "Product Absolute Count by Source") +
  scale_fill_viridis_c(name = "Product Count") +  # Add color scale as percentage
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )
# Print and save

ggsave(paste(output_dir, "/", "product_absolute_count_by_source_heatmap_plot.png", sep=''), plot = heatmap_plot_c, width = 12, height = 6, dpi = 800)

# Convert 'bioproject' to a factor
product_counts$Bioproject <- factor(product_counts$Bioproject)

# Create the heatmap for each product by Bioproject
heatmap_plot_biop <- ggplot(product_counts, aes(x = product_bigscape, y = Bioproject, fill = Product_Counts)) +
  geom_tile(color = "white") +
  labs(x = "Product", y = "Bioproject", title = "Product Absolute Count by Bioproject") +
  scale_fill_viridis_c(name = "Count") + ##if percentage, label = scales::percent_format(scale = 1)) 
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print and save 
ggsave(paste(output_dir, "/", "product_absolute_count_by_bioproject_heatmap_plot.png", sep=''), plot = heatmap_plot_biop, width = 12, height = 6, dpi = 300)

# Create the pie chart by Product_Counts
pie_chart_percentage <- ggplot(product_counts, aes(x = "", y = Normalized_Percentage, fill = product_bigscape)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Source) +  # Facet by Source
  labs(x = NULL, y = NULL, title = "Product Percentage by Source") +  # Remove axis labels
  coord_polar(theta = "y") +  # Use polar coordinates
  scale_fill_manual(values = rainbow(length(unique(product_counts$product_bigscape)))) +  # Color each product differently
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"))
# Print and save
ggsave(paste(output_dir, "/", "product_percentage_by_source_pie_plot.png", sep=''), plot = pie_chart_counts_absolute, width = 12, height = 6, dpi = 800)


pie_chart_counts_absolute <- ggplot(product_counts, aes(x = "", y = Product_Counts, fill = product_bigscape)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Source) +  # Facet by Source
  labs(x = NULL, y = NULL, title = "Product Absolute Count by Source") +  # Remove axis labels
  coord_polar(theta = "y") +  # Use polar coordinates
  scale_fill_manual(values = rainbow(length(unique(product_counts$product_bigscape)))) +  # Color each product differently
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )

ggsave(paste(output_dir, "/", "product_absolute_counts_by_source_pie_plot.png", sep=''), plot = pie_chart_counts_absolute, width = 12, height = 6, dpi = 800)

# Create the pie chart by Size
pie_chart_s <- ggplot(df_BGCs, aes(x = "", y = Size, fill = product_bigscape)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Source) +  # Facet by Source
  labs(x = NULL, y = NULL, title = "Product Size by Source") +  # Remove axis labels
  coord_polar(theta = "y") +  # Use polar coordinates
  scale_fill_manual(values = rainbow(length(unique(product_counts$product_bigscape)))) +  # Color each product differently
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"))

# Print and save

ggsave(paste(output_dir, "/", "product_size_by_source_pie_plot.png", sep=''), plot = pie_chart_s, width = 12, height = 6, dpi = 800)

####PCoA plot based on product occurrence by source
# Calculate the pairwise Euclidean distances between factors
dist_matrix <- dist(as.matrix(xtabs(Product_Counts ~ Source, product_counts)))

if (nrow(as.matrix(dist_matrix)) > 1 && ncol(as.matrix(dist_matrix)) > 1) {
  # Perform Principal Coordinate Analysis (PCoA)
  pcoa_result <- cmdscale(dist_matrix, k = 1)

  # Convert the result to a dataframe
  pcoa_data <- data.frame(PC1 = pcoa_result[, 1], PC2 = pcoa_result[, 1], product_bigscape = rownames(pcoa_result))

  # Plot the PCoA graph
  pcoa_plot <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = product_bigscape)) +
    geom_point(size = 3) +
    labs(title = "Principal Coordinate Analysis (PCoA) based on product occurrence between sources",
         x = "PC1", y = "PC2", color = "Source") +
    theme(
      text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
    )
  # Print and save

  ggsave(paste(output_dir, "/", "PCOA_product_count_plot_by_source.png", sep=''), plot = pcoa_plot, width = 12, height = 6, dpi = 800)
}


####PCoA plot based on product occurrence between products
# Calculate the pairwise Euclidean distances between factors
dist_matrix <- dist(as.matrix(xtabs(Percentage ~ product_bigscape, product_counts)))


if (nrow(as.matrix(dist_matrix)) > 1 && ncol(as.matrix(dist_matrix)) > 1) {
  # Perform Principal Coordinate Analysis (PCoA)
  pcoa_result <- cmdscale(dist_matrix, k = 2)

  # Convert the result to a dataframe
  pcoa_data <- data.frame(PC1 = pcoa_result[, 1], PC2 = pcoa_result[, 2], product_bigscape = rownames(pcoa_result))

  # Plot the PCoA graph
  pcoa_plot <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = product_bigscape)) +
    geom_point(size = 3) +
    labs(title = "Principal Coordinate Analysis (PCoA) based on product occurrence between samples",
         x = "PC1", y = "PC2", color = "Product") +
    theme(
      text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
    )+
    stat_ellipse(aes(group = product_bigscape), geom = "path", color = "red", linetype = "dashed", linewidth = 1)


  # Print and save

  ggsave(paste(output_dir, "/", "PCOA_product_count.png", sep=''), plot = pcoa_plot, width = 12, height = 6, dpi = 800)
}