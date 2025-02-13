import subprocess
import pandas as pd
import os

spades_output_dir = "../data/FastQ/assembly"
quast_output_dir = "../data/Quast"
reference_genome = "../data/FastA/SRR32230015.fasta"

contigs_file = os.path.join(spades_output_dir, "contigs.fasta")
scaffolds_file = os.path.join(spades_output_dir, "scaffolds.fasta")

if not os.path.exists(contigs_file):
    raise FileNotFoundError(f"Fichier contigs.fasta non trouvé dans {spades_output_dir}")

quast_cmd = [
    "quast",
    "-o", quast_output_dir,
    contigs_file,
    scaffolds_file,
]

# Ajouter la référence si disponible
if os.path.exists(reference_genome):
    quast_cmd.extend(["-r", reference_genome])

# Exécuter QUAST
print(f"Exécution de QUAST avec la commande : {' '.join(quast_cmd)}")
subprocess.run(quast_cmd, check=True)

# Lire les résultats de QUAST
report_file = os.path.join(quast_output_dir, "transposed_report.tsv")
if os.path.exists(report_file):
    df_results = pd.read_csv(report_file, sep="\t")
    print("\nRésultats de QUAST :")
    import ace_tools as tools
    tools.display_dataframe_to_user(name="Résultats QUAST", dataframe=df_results)
else:
    print("Échec : fichier de résultats introuvable.")
