library(ggrepel)

# Função para carregar e processar cada arquivo
process_file <- function(filename, splicing_type, species) {
  data <- read.table(filename, header = FALSE, col.names = c("SpliceSite"))
  # Contar a frequência de cada splice site
  data <- data %>%
    group_by(SpliceSite) %>%
    summarise(Count = n()) %>%
    ungroup()
  # Calcular a porcentagem de uso
  data <- data %>%
    mutate(TotalCount = sum(Count),
           Percentage = (Count / TotalCount) * 100,
           SplicingType = splicing_type,
           Species = species)  # Adiciona a espécie
  return(data)
}
# Carregar os dados de Maize
cs_Maize <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/maize_don_cs_ss.fa", "CS", "Maize")
ALTD_Maize <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/maize_don_ALTD_ss.fa", "ALTD", "Maize")
int_Maize <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/maize_don_INT_ss.fa", "INT", "Maize")
ex_Maize <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/maize_don_EX_ss.fa", "EX", "Maize")
# Carregar os dados de Arabidopsis
cs_arabidopsis <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/atha_don_cs_ss.fa", "CS", "Arabidopsis")
ALTD_arabidopsis <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/atha_don_ALTD_ss.fa", "ALTD", "Arabidopsis")
int_arabidopsis <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/atha_don_INT_ss.fa", "INT", "Arabidopsis")
ex_arabidopsis <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/atha_don_EX_ss.fa", "EX", "Arabidopsis")
# Carregar os dados de Humanos
cs_human <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/human_don_cs_ss.fa", "CS", "Human")
ALTD_human <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/human_acc_ALTA_ss.fa", "ALTD", "Human")
int_human <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/human_don_INT_ss.fa", "INT", "Human")
ex_human <- process_file("/home/bia/genetica2024_plots/BEPE/trimmed_cs/human_don_EX_ss.fa", "EX", "Human")
# Combinar todos os dados em um único dataframe
all_data <- bind_rows(cs_Maize, ALTD_Maize, int_Maize, ex_Maize,
                      cs_arabidopsis, ALTD_arabidopsis, int_arabidopsis, ex_arabidopsis,
                      cs_human, ALTD_human, int_human, ex_human)
# Definir a ordem desejada para a legenda
all_data$SplicingType_Species <- factor(interaction(all_data$Species, all_data$SplicingType),
                                        levels = c("Arabidopsis.CS", "Arabidopsis.ALTA", "Arabidopsis.ALTD",
                                                   "Arabidopsis.INT", "Arabidopsis.EX",
                                                   "Human.CS", "Human.ALTA", "Human.ALTD",
                                                   "Human.INT", "Human.EX",
                                                   "Maize.CS", "Maize.ALTA", "Maize.ALTD",
                                                   "Maize.INT", "Maize.EX"))
# Definir cores específicas para cada combinação de espécie e tipo de splicing
species_colors <- c("Arabidopsis.CS" = "darkgreen", "Arabidopsis.ALTA" = "forestgreen",
                    "Arabidopsis.ALTD" = "springgreen", "Arabidopsis.INT" = "limegreen",
                    "Arabidopsis.EX" = "palegreen",
                    "Human.CS" = "#6A0DAD", "Human.ALTA" = "#9370DB", "Human.ALTD" = "#BA55D3",
                    "Human.INT" = "#DA70D6", "Human.EX" = "#EE82EE",
                    "Maize.CS" = "#FFD700", "Maize.ALTA" = "#FFEC8B", "Maize.ALTD" = "#FFFACD",
                    "Maize.INT" = "#FFEFD5", "Maize.EX" = "#FFF8DC")
# Criar o gráfico com cores diferentes para cada espécie e tipo de splicing
ggplot(all_data, aes(x = SpliceSite, y = Percentage, color = SplicingType_Species)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = SplicingType), size = 3, max.overlaps = 10) +
  labs(title = "Splice site usage per species per splice type",
       x = "Splice Site",
       y = "Percentage Occurence (%)") +
  scale_color_manual(values = species_colors) +
  theme_minimal()
