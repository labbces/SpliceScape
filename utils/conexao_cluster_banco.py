imimport sqlite3
import pandas as pd
import os


CAMINHO_DO_ARQUIVO = "/Storage/data1/ellen.camargo/databaseLandscapeSlicing/grasses_with_genome_info.db"

print("--- INICIANDO TESTE ---")


if os.path.exists(CAMINHO_DO_ARQUIVO):
    print(f"✅ Arquivo encontrado: {CAMINHO_DO_ARQUIVO}")
else:
    print(f"ERRO: Arquivo não encontrado no caminho: {CAMINHO_DO_ARQUIVO}")
    exit()


try:
    conexao = sqlite3.connect(CAMINHO_DO_ARQUIVO)
    print("✅ Conexão realizada!")

    print("\n--- TABELAS DO BANCO ---")
    query = "SELECT name FROM sqlite_master WHERE type='table';"
    tabelas = pd.read_sql_query(query, conexao)
    print(tabelas)

except Exception as e:
    print(f"Erro na conexão: {e}")

finally:
    if 'conexao' in locals():
        conexao.close()
        print("\n--- FIM ---")
