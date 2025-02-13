import subprocess
import pandas as pd
import matplotlib.pyplot as plt

# Définition des fichiers
REFERENCE = "data/FastA/refSeq.fasta"  # Séquence de référence
READS_R1 = "data/FastQ/reads_R1.fastq"
READS_R2 = "data/FastQ/reads_R2.fastq"
SAM_FILE = "aligned_reads.sam"
BAM_FILE = "aligned_reads.bam"
SORTED_BAM = "aligned_reads_sorted.bam"
COVERAGE_FILE = "coverage.txt"


# Indexation de la séquence de référence
def index_reference():
    print("Indexation de la séquence de référence...")
    subprocess.run(f"bwa index {REFERENCE}", shell=True, check=True)


# Mapping des reads sur la séquence de référence
def align_reads():
    print("Alignement des reads...")
    subprocess.run(f"bwa mem -t 5 {REFERENCE} {READS_R1} {READS_R2} > {SAM_FILE}", shell=True, check=True)


# Conversion du fichier SAM en BAM et tri des reads
def convert_sort_bam():
    print("Conversion et tri du fichier BAM...")
    subprocess.run(f"samtools view -Sb {SAM_FILE} > {BAM_FILE}", shell=True, check=True)
    subprocess.run(f"samtools sort {BAM_FILE} -o {SORTED_BAM}", shell=True, check=True)
    subprocess.run(f"samtools index {SORTED_BAM}", shell=True, check=True)


# Calcul de la couverture des reads sur le génome de référence
def compute_coverage():
    print(" Calcul de la couverture...")
    subprocess.run(f"samtools depth {SORTED_BAM} > {COVERAGE_FILE}", shell=True, check=True)


# Visualisation de la couverture des reads
def plot_coverage():
    print("Génération du graphique de couverture...")
    df = pd.read_csv(COVERAGE_FILE, sep="\t", header=None, names=["Chromosome", "Position", "Couverture"])

    plt.figure(figsize=(12, 6))
    plt.plot(df["Position"], df["Couverture"], color="blue", alpha=0.7)
    plt.xlabel("Position sur le génome")
    plt.ylabel("Profondeur de couverture")
    plt.title("Couverture des reads sur le génome du SARS-CoV-2")
    plt.grid()
    plt.show()


# Exécution du pipeline
if __name__ == "__main__":
    index_reference()
    align_reads()
    convert_sort_bam()
    compute_coverage()
    plot_coverage()
