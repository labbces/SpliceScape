from collections import Counter
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

colunas = ["chrom", "start", "end", "event_id", "score", "strand"]
df = pd.read_csv(
    "/home/bia/LandscapeSplicingGrasses/Results/corrected_variants_plus.bed", sep="\t", names=colunas)

bin_size = 100000
df['bin'] = df['start'] // bin_size

# 3. Encontrar o último bin necessário para cada cromossomo
# Obtém o maior valor de `start` em cada cromossomo para calcular o bin máximo
max_bins = df.groupby('chrom')['start'].max() // bin_size

# 4. Contar a densidade de eventos por cromossomo e bin
heatmap_data = df.pivot_table(
    index="chrom", columns="bin", values="event_id", aggfunc="count", fill_value=0)

# 5. Aplicar a máscara para excluir bins além do tamanho de cada cromossomo
for chrom, max_bin in max_bins.items():
    # Define como NaN para esconder no heatmap
    heatmap_data.loc[chrom, max_bin+1:] = np.nan

# 6. Plotar o Heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_data, cmap="viridis", cbar_kws={
            'label': 'Densidade de Eventos'}, mask=heatmap_data.isna())

plt.title(
    "Mapa de Calor da Densidade de Eventos de Splicing por Cromossomo e Posição")
plt.xlabel("Posição ao longo do cromossomo (bins de {} bp)".format(bin_size))
plt.ylabel("Cromossomo")
plt.show()


# 1. Carregar os dados
colunas = ["chrom", "start", "end", "event_id", "score", "strand"]
df = pd.read_csv(
    "/home/bia/LandscapeSplicingGrasses/Results/corrected_variants_plus.bed", sep="\t", names=colunas)

# 2. Extrair e filtrar os tipos de splicing
# Extrai o tipo após o último "_" da coluna 'event_id' e separa por vírgula
df['splicing_types'] = df['event_id'].apply(
    lambda x: x.split('_')[-1].split(','))

# 3. Filtrar apenas os tipos de splicing desejados
tipos_desejados = {"A5SS", "RI", "A3SS", "MXE", "SE"}
df['splicing_types'] = df['splicing_types'].apply(
    lambda x: [t for t in x if t in tipos_desejados])

# 4. Contar a ocorrência de cada tipo de splicing
# Expande a coluna em várias linhas para contar cada tipo separadamente
splicing_types = df['splicing_types'].explode(
).dropna()  # Remove valores vazios
splicing_count = splicing_types.value_counts()

# 5. Plotar gráfico de pizza
plt.figure(figsize=(8, 8))
plt.pie(splicing_count, labels=splicing_count.index,
        autopct='%1.1f%%', startangle=140, colors=plt.cm.tab20.colors)
plt.title("Distribuição dos Tipos de Splicing em Arabidopsis")
plt.show()
