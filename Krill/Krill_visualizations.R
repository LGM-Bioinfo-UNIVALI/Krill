#######################testing frame
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(gridExtra)
library(igraph)

library(plotly)


args <- commandArgs(trailingOnly = TRUE)

# Read the TSV file
metagenome <- read.delim(args[1], stringsAsFactors = T)

output_dir <- args[2]


# Extract "Bioproject" and "Source" from the "database" column
metagenome <- metagenome %>%
  mutate(Bioproject = sub("_(.*?)$", "", Database),
         Source = sub("^(.*?)_", "", Database))


# Convert "product" column to character
metagenome$product <- as.character(metagenome$product)

# Split values in "product" and create separate rows
df_metagenome <- metagenome %>%
  mutate(product_copy = product) %>%
  separate_rows(product, sep = ",") %>%
  select(Bioproject, product, Size, Source)

# Define the internal class separation criteria
  replacement_rules <- list(
  "PKS I" = c("t1pks", "T1PKS"),
  "PKS other" = c("transatpks", "t2pks", "t3pks", "otherks", "hglks", "transAT-PKS", "transAT-PKS-like", "T2PKS", "T3PKS", "PKS-like", "hglE-KS"),
  "NRPS" = c("nrps", "NRPS", "NRPS-like", "thioamide-NRP"),
  "RiPPs" = c("lantipeptide", "thiopeptide", "bacteriocin", "linaridin", "cyanobactin", "glycocin", "LAP", "lassopeptide", "sactipeptide", "bottromycin", "head_to_tail", "microcin", "microviridin", "proteusin", "lanthipeptide", "lipolanthine", "RaS-RiPP", "fungal-RiPP"),
  "Saccharides" = c("amglyccycl", "oligosaccharide", "cf_saccharide", "saccharide"),
  "Terpene" = "terpene",
  "PKS/NRPS Hybrids" = c("t1pks", "T1PKS", "transatpks", "t2pks", "t3pks", "otherks", "hglks", "transAT-PKS", "T2PKS", "T3PKS", "PKS-like", "hglE-KS", "nrps", "NRPS", "NRPS-like", "thioamide-NRP"),
  "Others" = c("acyl_amino_acids", "arylpolyene", "aminocoumarin", "ectoine", "butyrolactone", "nucleoside", "melanin", "phosphoglycolipid", "phenazine", "phosphonate", "other", "cf_putative", "resorcinol", "indole", "ladderane", "PUFA", "furan", "hserlactone", "fused", "cf_fatty_acid", "siderophore", "blactam", "fatty_acid", "PpyS-KS", "CDPS", "betalactone", "PBDE", "tropodithietic-acid", "NAGGN", "halogenated")
)

# Function to replace values based on the internal class separation criteria
replace_values <- function(value) {
  for (class in names(replacement_rules)) {
    patterns <- replacement_rules[[class]]
    if (grepl(paste0("\\b(", paste(patterns, collapse = "|"), ")\\b"), value)) {
      return(class)
    }
  }
  return(value)
}

# Replace values in the "product" column based on the internal class separation criteria
df_metagenome$product <- sapply(df_metagenome$product, replace_values)


# Define the list of product categories
product_categories <- c("NAPAA", "NRPS", "PKS I", "Others", "PKS other", "ranthipeptide", "redox-cofactor", "RiPP-like", "RiPPs", "RRE-containing", "Saccharides", "Terpene", "thioamitides")

# Create a list to store the plots
plots <- list()

# Iterate over each product category
for (product in product_categories) 
  # Subset the data for the current product category
  subset_data <- df_metagenome[df_metagenome$product == product, ]
  
# Count
product_count <- df_metagenome %>%
  group_by(product) %>%
  summarise(count = n_distinct(Bioproject))  

# Ordena os produtos por contagem decrescente
product_count <- product_count %>%
  arrange(desc(count))

color_palette <- c(
  "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
  "#FF8000", "#800080", "#008000", "#800000", "#000080", "#FFA500", "#008080",
  "#FF5733", "#33FF57", "#5733FF", "#FFFF33", "#FF33FF", "#33FFFF",
  "#FF9966", "#9966FF", "#66FF99", "#FFCC66", "#66FFCC", "#CC66FF", "cyan"
)


# Create the bar plot
bar_plot <- ggplot(product_count, aes(x = reorder(product, -count), y = count, fill = product)) +
  geom_bar(stat = "identity") +
  labs(x = "Product", y = "Count", title = "Product Total Absolute Count") +
  scale_fill_manual(values = color_palette) +  # Define the color palette
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  )

# Display the bar plot and save
ggsave(paste(output_dir, "/", "product_total_absolute_count_bar_plot.png", sep=''), plot = bar_plot, width = 10, height = 6, dpi = 300)

# Create a histogram for each product

for (product in unique(df_metagenome$product)) {
  # Subset data for the current product
  product_data <- df_metagenome[df_metagenome$product == product, ]
  
  
  
  # Create the histogram for the current product
  plot <- ggplot(product_data, aes(x = Size, fill = product)) +
    geom_histogram(binwidth = 1000, color = "black", alpha = 0.7) +
    labs(x = "Size", y = "Frequency", title = paste("Size Distribution for", product)) +
    theme_minimal() +
    theme(
      text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    )
  
  
  # Store the plot in the list
  plots[[product]] <- plot
}

# Print each plot separately
for (product in names(plots)) {
  ggsave(paste(output_dir, "/", product, "_histogram_size_plot.png", sep=''), plots[[product]], width = 10, height = 6, dpi = 300)
}



# Count how many times a product category occurs in each category
source_product_counts <- df_metagenome %>%
  group_by(Bioproject, Source, product) %>%
  tally() %>%
  rename(Product_Counts = n)

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
bar_graph <- ggplot(product_counts, aes(x = Source, y = Normalized_Percentage, fill = product)) +
  geom_bar(stat = "identity") +
  labs(x = "Source", y = "Percentage", title = "Product Percentage by Source") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) +
  coord_cartesian(ylim = c(0, 100))   

# Print and save the graph
ggsave(paste(output_dir, "/", "product_percentage_by_source_stackedbar_plot.png", sep=''), plot = bar_graph, width = 10, height = 6, dpi = 300)


# Heatmap by Source 
heatmap_plot <- ggplot(product_counts, aes(x = Source, y = Product_Counts, fill = product)) +
  geom_tile() +
  labs(x = "Source", y = "Count") +
  ggtitle("Product Frequency by  Source") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  )

#Print and save
ggsave(paste(output_dir, "/", "product_frequecy_by_source_heatmap_plot.png", sep=''), plot = heatmap_plot, width = 10, height = 6, dpi = 300)

# Create a bubble plot with product size
bubble_plot_s <- ggplot(df_metagenome, aes(x = Source, y = Size, size = Size, color = product)) +
  geom_point() +
  labs(x = "Source", y = "Size", size = "Size", color = "Product") +
  ggtitle("Product Size") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  )
#Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_bubble_plot.png", sep=''), plot = bubble_plot_s, width = 10, height = 6, dpi = 300)

# Boxplot of Product Absolute Count by Source
boxplot_plot_pc <- ggplot(product_counts, aes(x = product, y = Product_Counts, fill = Source)) +
  geom_boxplot(width = 0.7) +  # Adjust the width as needed
  labs(x = "Product", y = "Product_Counts", title = "Boxplot: Product Absolute Count by Source") +
    theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print and save
ggsave(paste(output_dir, "/", "product_absolute_count_by_source_box_plot.png", sep=''), plot = boxplot_plot_pc, width = 10, height = 6, dpi = 300)

# Create the boxplot with facet for Bioproject
boxplot_plot2 <- ggplot(product_counts, aes(x = product, y = Product_Counts, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Product_Counts", title = "Boxplot: Product Absolute Count by Source and Bioproject") +
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
boxplot_plot_ps <- ggplot(df_metagenome, aes(x = product, y = Size, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Size", title = "Boxplot: Product Size by Source") +
   theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )


# Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_box_plot.png", sep=''), plot = boxplot_plot_ps, width = 10, height = 6, dpi = 300)

# Create the boxplot with facet for Bioproject
boxplot_plot2s <- ggplot(df_metagenome, aes(x = product, y = Size, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Size", title = "Boxplot: Product Size by Condition and Bioproject") +
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
boxplot_plot_p <- ggplot(product_counts, aes(x = product, y = Percentage, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Percentage", title = "Boxplot: Product Percentage by Source") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white")# Rotate x-axis labels for better readability
  )

# Print and save
ggsave(paste(output_dir, "/", "product_percentage_by_source_box_plot.png", sep=''), plot = boxplot_plot_p, width = 12, height = 6, dpi = 300)

# Create the boxplot with facet for Bioproject
boxplot_plot2p <- ggplot(product_counts, aes(x = product, y = Percentage, fill = Source)) +
  geom_boxplot() +
  labs(x = "Product", y = "Percentage", title = "Boxplot: Product Percentage by Source and Bioproject") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    panel.background = element_rect(fill = "white")# Rotate x-axis labels for better readability
  ) +
  facet_grid(rows = vars(Bioproject))  # Facet for Bioproject

# Print and save
ggsave(paste(output_dir, "/", "product_percentage_by_source_and_bioproject_box_plot.png", sep=''), plot = boxplot_plot2p, width = 12, height = 6, dpi = 300)

# Create the bubble plot with Product size by Source and Bioproject
bubble_plot_c <- ggplot(df_metagenome, aes(x = product, y = Source, size = Size, color = product)) +
  geom_point(alpha = 0.7) +
  labs(x = "Product", y = "Source", title = "Bubble Plot: Product Size by Source and Bioproject") +
  scale_size_continuous(range = c(3, 15)) +  # Set the range of bubble sizes
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) +
  facet_grid(Bioproject ~ ., scales = "free_y")  # Create facets for each Bioproject

# Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_and_bioproject_box_plot.png", sep=''), plot = boxplot_plot2s, width = 10, height = 6, dpi = 300)


# Create the bubble plot with Product Absolute Count by Source and Bioproject
bubble_plot_p <- ggplot(product_counts, aes(x = product, y = Source, size = Product_Counts, color = product)) +
  geom_point(alpha = 0.7) +
  labs(x = "Product", y = "Source", title = "Bubble Plot: Product Absolute Count by Source and Bioproject") +
  scale_size_continuous(range = c(3, 15)) +  # Set the range of bubble sizes
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) +
  facet_grid(Bioproject ~ ., scales = "free_y")  # Create facets for each Bioproject

# Print and save
ggsave(paste(output_dir, "/", "product_absolute_count_by_source_and_bioproject_bubble_plot.png", sep=''), plot = bubble_plot_p, width = 12, height = 6, dpi = 300)

# Create the heatmap for each product by Source
heatmap_plot_c <- ggplot(product_counts, aes(x = product, y = Source, fill = Product_Counts)) +
  geom_tile(color = "white") +
  labs(x = "Product", y = "Source", title = "Heatmap: Product Absolute Count by Source") +
  scale_fill_viridis_c(name = "Product Count") +  # Add color scale as percentage
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )
# Print and save
ggsave(paste(output_dir, "/", "product_absolute_count_by_source_heatmap_plot.png", sep=''), plot = heatmap_plot_c, width = 12, height = 6, dpi = 300)

# Convert 'bioproject' to a factor
product_counts$Bioproject <- factor(product_counts$Bioproject)

# Create the heatmap for each product by Bioproject
heatmap_plot_biop <- ggplot(product_counts, aes(x = product, y = Bioproject, fill = Product_Counts)) +
  geom_tile(color = "white") +
  labs(x = "Product", y = "Bioproject", title = "Heatmap: Product Absolute Count by Bioproject") +
  scale_fill_viridis_c(name = "Product Count") + ##if percentage, label = scales::percent_format(scale = 1)) 
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print and save 
ggsave(paste(output_dir, "/", "product_absolute_count_by_bioproject_heatmap_plot.png", sep=''), plot = heatmap_plot_biop, width = 12, height = 6, dpi = 300)


# Parallel coordinates plot for Bioproject by Product Absolute Counts
parallel_plot <- ggplot(product_counts, aes(x = product, y = Product_Counts, group = Bioproject, color = Bioproject)) +
  geom_point() +  # Use geom_point() instead of geom_line()
  labs(x = "Product", y = "Product Count", title = "Product Absolute Counts") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),
    legend.position = "right"
  )

# Print and save
ggsave(paste(output_dir, "/", "product_absolute_counts_by_bioproject_parallel_plot.png", sep=''), plot = parallel_plot, width = 12, height = 6, dpi = 300)


# Parallel coordinates plot for Source
parallel_plot_si <- ggplot(product_counts, aes(x = product, y = Product_Counts, group = Source, color = Source)) +
  geom_point() +  # Use geom_point() instead of geom_line()
  labs(x = "Product", y = "Product Count", title = "Product Absolute Counts by Source") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),
    legend.position = "right"
  )


# Print and save
ggsave(paste(output_dir, "/", "product_absolute_counts_by_source_parallel_plot.png", sep=''), plot = parallel_plot_si, width = 12, height = 6, dpi = 300)


# Parallel coordinates plot for Source and Product Size
parallel_plot_s <- ggplot(df_metagenome, aes(x = product, y = Size, group = Source, color = Source)) +
  geom_point() +  # Use geom_point() instead of geom_line()
  labs(x = "Product", y = "Size", title = "Product Size by Source") +
  theme(
    text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),
    legend.position = "right"
  )
# Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_parallel_plot.png", sep=''), plot = parallel_plot_s, width = 12, height = 6, dpi = 300)

# Create the pie chart by Percentage
pie_chart_percentage <- ggplot(product_counts, aes(x = "", y = Normalized_Percentage, fill = product)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Source) +  # Facet by Source
  labs(x = NULL, y = NULL, title = "Product Percentage by Source") +  # Remove axis labels
  coord_polar(theta = "y") +  # Use polar coordinates
  scale_fill_manual(values = rainbow(length(unique(product_counts$product)))) +  # Color each product differently
  theme_void()  # Remove unnecessary elements

# Print and save
ggsave(paste(output_dir, "/", "product_percentage_by_source_pie_plot.png", sep=''), plot = pie_chart_percentage, width = 12, height = 6, dpi = 300)

# Create the pie chart by Product_Counts
pie_chart_counts_absolute <- ggplot(product_counts, aes(x = "", y = Product_Counts, fill = product)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Source) +  # Facet by Source
  labs(x = NULL, y = NULL, title = "Product Absolute Count by Source") +  # Remove axis labels
  coord_polar(theta = "y") +  # Use polar coordinates
  scale_fill_manual(values = rainbow(length(unique(product_counts$product)))) +  # Color each product differently
  theme_void()  # Remove unnecessary elements

# Print and save
ggsave(paste(output_dir, "/", "product_absolute_counts_by_source_pie_plot.png", sep=''), plot = pie_chart_counts_absolute, width = 12, height = 6, dpi = 300)


# Create the pie chart by Size
pie_chart_s <- ggplot(df_metagenome, aes(x = "", y = Size, fill = product)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Source) +  # Facet by Source
  labs(x = NULL, y = NULL, title = "Product Size by Source") +  # Remove axis labels
  coord_polar(theta = "y") +  # Use polar coordinates
  scale_fill_manual(values = rainbow(length(unique(product_counts$product)))) +  # Color each product differently
  theme_void()  # Remove unnecessary elements

# Print and save
ggsave(paste(output_dir, "/", "product_size_by_source_pie_plot.png", sep=''), plot = pie_chart_s, width = 12, height = 6, dpi = 300)






