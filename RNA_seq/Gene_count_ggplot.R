# Install and load the ggplot2 package
install.packages("ggplot2")
library(ggplot2)

# Create the data frame
data <- data.frame(
  Category = c("TPE genes", "Stress genes", "Others"),
  Upregulated = c(1, 2, 7),
  Downregulated = c(21, 1, 21)
)

# Convert the data to long format
data_long <- reshape2::melt(data, id.vars = "Category", variable.name = "Condition", value.name = "No_of_genes")

# Define the stacking order of conditions
data_long$Condition <- factor(
  data_long$Condition,
  levels = c("Upregulated", "Downregulated")
)

# Define the order of categories on the x-axis
data_long$Category <- factor(
  data_long$Category,
  levels = c("TPE genes", "Stress genes", "Others")
)

# Plot the stacked bar chart
ggplot(data_long, aes(x = Category, y = No_of_genes, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("red", "forestgreen")) +
  labs(y = "No. of genes", x = "", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(
  aes(label = No_of_genes),
  position = position_stack(vjust = 0.5),
  size = 5,
  fontface = "bold",
  colour = "black"
)

ggsave(
  "Gene_count_ggplot.png",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)
