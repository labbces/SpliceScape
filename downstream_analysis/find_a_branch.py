import sqlite3
from pyfaidx import Fasta
import argparse
import re
import logomaker
import matplotlib.pyplot as plt
import pandas as pd

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement_dict.get(base, 'N') for base in reversed(seq.upper()))

def main():
    parser = argparse.ArgumentParser(
        description="Predicts branch sites, precisely extracts flanking sequences, generates a colored logo, and calculates summary statistics."
    )
    parser.add_argument("--db", required=True, help="Path to the SQLite database file.")
    parser.add_argument("--fasta", required=True, help="Path to the genome FASTA file.")
    parser.add_argument("--species", required=True, help="Name of the species to filter from the database.")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files.")
    parser.add_argument("--upstream_window", type=int, default=50, help="Bases upstream of 3'ss to search for the motif.")
    parser.add_argument("--motif", type=str, default="[CT]T[AG]A[CT]", help="Regex for YURAY motif.")
    parser.add_argument("--motif_flank", type=int, default=10, help="Bases to extract upstream and downstream of the found motif for the logo.")
    args = parser.parse_args()

    # Output file paths
    csv_output_path = f"{args.output_prefix}_predictions.csv"
    logo_output_path = f"{args.output_prefix}_logo.png"

    print(f"Indexing genome: {args.fasta}...")
    genome = Fasta(args.fasta)
    print("Indexing complete.")

    branch_site_regex = re.compile(args.motif, re.IGNORECASE)
    
    flanked_sequences = []
    total_events_processed = 0
    motifs_found_count = 0

    con = sqlite3.connect(args.db)
    cursor = con.cursor()

    with open(csv_output_path, 'w') as outfile:
        outfile.write("event_id,strand,three_prime_ss_coord,branch_site_motif_found,distance_from_3ss\n")
        
        print(f"Searching for branch sites in '{args.species}'...")
        query = "SELECT event_id, seqid, start, end, strand FROM splicing_events WHERE species = ?"
        
        for event_id, seqid, start, end, strand in cursor.execute(query, (args.species,)):
            try:
                total_events_processed += 1

                if strand == '+':
                    three_prime_ss_coord = end
                    search_window_start_genomic = max(1, three_prime_ss_coord - args.upstream_window)
                    search_window_end_genomic = three_prime_ss_coord
                else: # strand == '-'
                    three_prime_ss_coord = start
                    search_window_start_genomic = three_prime_ss_coord
                    search_window_end_genomic = three_prime_ss_coord + args.upstream_window
                
                search_sequence = genome[seqid][search_window_start_genomic-1:search_window_end_genomic].seq
                if strand == '-':
                    search_sequence = reverse_complement(search_sequence)
                
                matches = list(branch_site_regex.finditer(search_sequence))
                
                if matches:
                    motifs_found_count += 1
                    best_match = matches[-1]
                    motif_found = best_match.group(0)
                    
                    distance = len(search_sequence) - best_match.end()
                    outfile.write(f"{event_id},{strand},{three_prime_ss_coord},{motif_found},-{distance}\n")
                    

                    if strand == '+':
                        motif_start_genomic = search_window_start_genomic + best_match.start()
                        motif_end_genomic = search_window_start_genomic + best_match.end() - 1
                    else: 
                        motif_start_genomic = search_window_end_genomic - best_match.end() + 1
                        motif_end_genomic = search_window_end_genomic - best_match.start()

                    final_extract_start = motif_start_genomic - args.motif_flank
                    final_extract_end = motif_end_genomic + args.motif_flank
                    final_sequence = genome[seqid][final_extract_start-1:final_extract_end].seq
                    
                    if strand == '-':
                        final_sequence = reverse_complement(final_sequence)
                    
                    flanked_sequences.append(final_sequence)
                else:
                    outfile.write(f"{event_id},{strand},{three_prime_ss_coord},No_Motif_Found,NA\n")

            except Exception as e:
                print(f"WARNING: Could not process event {event_id}. Error: {e}")
                continue
    
    con.close()
    print(f"\nSearch complete! Results saved to: {csv_output_path}")

    if flanked_sequences:
        print("Generating branch site sequence logo...")
        info_matrix = logomaker.alignment_to_matrix(sequences=flanked_sequences, to_type='information')
        
        branch_logo = logomaker.Logo(info_matrix, font_name='Arial Rounded MT Bold', color_scheme='classic')
        
        branch_logo.style_spines(visible=False)
        branch_logo.ax.set_title(f"Região flanqueante ao Sítio de Ramificação ({args.species})")
        branch_logo.ax.set_xlabel(f"Posição relativa ao Motif {args.motif_flank}bp ")
        branch_logo.ax.set_ylabel("Bits")
        plt.savefig(logo_output_path, dpi=300)
        plt.close()
        print(f"Logo saved to: {logo_output_path}")
    else:
        print("No branch site motifs were found, skipping logo generation.")

    if total_events_processed > 0:
        percentage_found = (motifs_found_count / total_events_processed) * 100
        print("\n--- Summary ---")
        print(f"Total events processed: {total_events_processed}")
        print(f"Events with branch site motif found: {motifs_found_count}")
        print(f"Percentage of events with motif: {percentage_found:.2f}%")
    else:
        print("\nNo events were processed.")

if __name__ == "__main__":
    main()