import sqlite3
import pandas as pd
import pysam
import os
import sys


DB_PATH = "/Storage/data1/ellen.camargo/databaseLandscapeSlicing/grasses_with_genome_info.db"

FASTA_PATH = "/Storage/data1/ellen.camargo/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/assembly/Athaliana_447_TAIR10.fa"

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def main():
    print(f"--- ANÁLISE FINAL V2 (CORREÇÃO FINA DE SHIFT) ---")

    if not os.path.exists(FASTA_PATH):
        print("❌ ERRO: Arquivo FASTA não encontrado.")
        sys.exit()

    try:
        genome = pysam.FastaFile(FASTA_PATH)
        cromossomos_validos = set(genome.references)
    except Exception as e:
        print(f"❌ Erro no FASTA: {e}")
        sys.exit()

    conn = sqlite3.connect(DB_PATH)
    query = "SELECT seqid, start, end, strand, event_type FROM splicing_events"
    df = pd.read_sql_query(query, conn)
    conn.close()

    print(f"✅ {len(df)} eventos carregados. Aplicando correção dupla (+1 no Start)...")

    donors = []
    acceptors = []

    for index, row in df.iterrows():
        chrom_db = str(row['seqid'])


        if chrom_db in cromossomos_validos:
            target_chrom = chrom_db
        elif ("Chr" + chrom_db) in cromossomos_validos:
            target_chrom = "Chr" + chrom_db
        else:
            donors.append("IGNORE")
            acceptors.append("IGNORE")
            continue

        try:
            start = int(row['start'])
            end = int(row['end'])
            strand = str(row['strand'])


            seq = genome.fetch(target_chrom, start - 1, end).upper()

            if not seq:
                donors.append("NA")
                acceptors.append("NA")
                continue


            if strand == '+' or strand == '1':


                d_site = seq[1:3]


                a_site = seq[-3:-1]

            elif strand == '-' or strand == '-1':



                site_raw_end = seq[-3:-1]
                d_site = reverse_complement(site_raw_end)


                site_raw_start = seq[1:3]
                a_site = reverse_complement(site_raw_start)

            else:
                d_site = "NA"
                a_site = "NA"

            donors.append(d_site)
            acceptors.append(a_site)

        except:
            donors.append("ERRO")
            acceptors.append("ERRO")

    df['donor'] = donors
    df['acceptor'] = acceptors

    df_clean = df[~df['donor'].isin(["IGNORE", "ERRO", "NA"])]


    df_gold = df_clean[ (df_clean['donor'] == 'GT') & (df_clean['acceptor'] == 'AG') ]

    print("\n" + "="*40)
    print("      RESULTADOS FINAIS V2")
    print("="*40)

    print(f"Total Analisado (Arabidopsis): {len(df_clean)}")
    print(f"Total Perfeito (GT-AG): {len(df_gold)}")

    if len(df_clean) > 0:
        print("\n🔹 TOP 5 DOADORES (Geral):")
        print(df_clean['donor'].value_counts().head())

        print("\n🔹 TOP 5 ACEITADORES (Geral):")
        print(df_clean['acceptor'].value_counts().head())


        df_clean.to_csv("splicing_arabidopsis_final_v2.csv", index=False)
        print("\n✅ Salvo em 'splicing_arabidopsis_final_v2.csv'")
    else:
        print("\n❌ Erro: Nenhum dado sobrou.")

if __name__ == "__main__":
    main()
