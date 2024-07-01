# Instalar e carregar os pacotes necessários
# install.packages("ggplot2")
library(ggplot2)

# Criar um dataframe com os valores fornecidos
data <- data.frame(
  Species = c("Z. mays", "S. bicolor", "H. vulgare", "A. thaliana", "T. aestivum", "S. italica"),
  Especies = c("Z. mays", "S. bicolor", "H. vulgare", "A. thaliana", "T. aestivum", "S. italica"),
  Values = c(27924, 6994, 7860, 34946, 15734, 1595)
)

# Definir as cores em formato hex
colors <- c("Z. mays" = "#FFEB3B", "S. bicolor" = "#1E88E5", "H. vulgare" = "#e31a1c",
            "A. thaliana" = "#33a02c", "T. aestivum" = "#6a3d9a", "S. italica" = "#E91E63")

# PORTUGUES 
# Criar o gráfico de barras com cores personalizadas
ggplot(data, aes(x = Especies, y = Values, fill = Especies)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Quantidade de datasets SRA/NCBI",
       x = "Espécie",
       y = "Total de SRA")

# INGLES 
# Criar o gráfico de barras com cores personalizadas
ggplot(data, aes(x = Species, y = Values, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Total datasets (SRA/NCBI) available",
       x = "Species",
       y = "Total SRA datasets")
