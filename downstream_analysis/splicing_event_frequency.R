# Load the ggplot2 library
library(ggplot2)

# Define custom colors for the categories
custom_colors <- c("CS" = "#fdaaad", "EX" = "#CC99FF", "INT" = "#99FF99", "ALTA" = "#add5ff", "ALTD" = "#2e8cee")

# HG18

# Create a sample dataset
data <- data.frame(
  category = c("ALTA", "ALTD", "INT", "EX", "CS"),
  count = c(3634, 5332, 5150, 24511, 9043)
)

# Calculate the percentage for each category

data$percentage <- data$count / sum(data$count) * 100
# Create the pie chart
ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +  # Remove background, grid, and numeric labels
  theme(legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(hjust = 0.5)) +  # Center the title
  labs(title = "HG18 splicing events") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = custom_colors)  # Apply custom colors




# HG38

# Create a sample dataset
data <- data.frame(
  category = c("ALTA", "ALTD", "INT", "EX", "CS"),
  count = c(11729,19849,34905,85268,128437)
)

# Calculate the percentage for each category

data$percentage <- data$count / sum(data$count) * 100
# Create the pie chart
ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +  # Remove background, grid, and numeric labels
  theme(legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(hjust = 0.5)) +  # Center the title
  labs(title = "HG38 splicing events") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = custom_colors)  # Apply custom colors


# A.thaliana

# Create a sample dataset
data <- data.frame(
  category = c("ALTA", "ALTD", "INT", "EX", "CS"),
  count = c(85049,43643,115163,92764,29886)
)

# Calculate the percentage for each category

data$percentage <- data$count / sum(data$count) * 100
# Create the pie chart
ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +  # Remove background, grid, and numeric labels
  theme(legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(hjust = 0.5)) +  # Center the title
  labs(title = "A. thaliana splicing events") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = custom_colors)  # Apply custom colors

# Maize

# Create a sample dataset
data <- data.frame(
  category = c("ALTA", "ALTD", "INT", "EX", "CS"),
  count = c(24748,17871,46416,14531,4410)
)

# Calculate the percentage for each category

data$percentage <- data$count / sum(data$count) * 100
# Create the pie chart
ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +  # Remove background, grid, and numeric labels
  theme(legend.title = element_blank(),  # Remove legend title
        plot.title = element_text(hjust = 0.5)) +  # Center the title
  labs(title = "Maize splicing events") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = custom_colors)  # Apply custom colors

################################################



