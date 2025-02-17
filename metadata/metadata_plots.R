# Instalar e carregar os pacotes necessários
# install.packages("ggplot2")
library(ggplot2)

# Species: Zea mays; Sorghum bicolor; Hordeum vulgare; Arabidopsis thaliana; Triticum aestivum; Setaria italica

# Criar um dataframe com os valores fornecidos
data <- data.frame(
  Species = c("Z. mays", "S. bicolor", "H. vulgare", "A. thaliana", "T. aestivum", "S. italica", "O. sativa"),
  Especies = c("Z. mays", "S. bicolor", "H. vulgare", "A. thaliana", "T. aestivum", "S. italica", "O. sativa"),
  Values = c(27942, 7154, 8141, 28156, 16118, 1606, 23034)
)

# Dados filtrados
data <- data.frame(
  Species = c("Z. mays", "S. bicolor", "H. vulgare", "A. thaliana", "T. aestivum", "S. italica", "O. sativa"),
  Especies = c("Z. mays", "S. bicolor", "H. vulgare", "A. thaliana", "T. aestivum", "S. italica", "O. sativa"),
  Values = c(26948, 7115, 8006, 26591, 16053, 1577, 22880)
)

# Definir as cores em formato hex
colors <- c("Z. mays" = "#FFEB3B", "S. bicolor" = "#1E88E5", "H. vulgare" = "#e31a1c",
            "A. thaliana" = "#33a02c", "T. aestivum" = "#6a3d9a", "S. italica" = "#E91E63", "O. sativa" = "#f97e43")

# PORTUGUES
ggplot(data, aes(x = Especies, y = Values, fill = Especies)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Values), vjust = -0.5) + 
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Quantidade de datasets SRA/NCBI",
       x = "Espécie",
       y = "Total de SRA") +
  theme(
    plot.title = element_text(hjust = 0.5), 
    axis.text.x = element_text(face = "italic") 
  )

# INGLES
ggplot(data, aes(x = Species, y = Values, fill = Species)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Values), vjust = -0.5) + 
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Total datasets (SRA/NCBI) available",
       x = "Species",
       y = "Total SRA datasets") +
  theme(
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(face = "italic") 
  )

