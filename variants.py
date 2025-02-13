import os
import io
import gzip
import shutil
import subprocess
import pandas as pd

# Définition des fichiers
REFERENCE_FASTA = "data/FastA/refSeq.fasta"
DICT_FILE = REFERENCE_FASTA.replace(".fasta", ".dict")
FAI_FILE = REFERENCE_FASTA + ".fai"
SORTED_BAM = "aligned_reads_sorted.bam"
VCF_FILE = "output_variants_sarscov2.vcf.idx"
ANNOTATED_VCF = "annotated_variants.vcf"

# Vérifier si un fichier existe
def check_file(filepath):
    """Vérifie si un fichier existe, sinon affiche un message d'erreur et arrête le script."""
    if not os.path.exists(filepath):
        print(f"ERREUR : Fichier {filepath} introuvable. Vérifiez le chemin.")
        exit(1)

# Générer l’index du FASTA (.fai) avec samtools
def generate_fasta_index():
    if os.path.exists(FAI_FILE):
        print(f"Index FASTA déjà généré : {FAI_FILE}")
        return

    print("Génération de l’index FASTA (.fai)...")
    try:
        subprocess.run(f"samtools faidx {REFERENCE_FASTA}", shell=True, check=True)
        print(f"Index FASTA généré : {FAI_FILE}")
    except subprocess.CalledProcessError:
        print("ERREUR : Échec de la génération de l'index FASTA.")
        exit(1)

# Générer le dictionnaire du FASTA (.dict) avec GATK
def generate_fasta_dict():
    if os.path.exists(DICT_FILE):
        print(f"Dictionnaire FASTA déjà existant : {DICT_FILE}")
        return

    print("Génération du dictionnaire FASTA (.dict)...")
    try:
        subprocess.run(f"gatk CreateSequenceDictionary -R {REFERENCE_FASTA} -O {DICT_FILE}", shell=True, check=True)
        print(f"Dictionnaire FASTA généré : {DICT_FILE}")
    except subprocess.CalledProcessError:
        print("ERREUR : Échec de la génération du dictionnaire FASTA.")
        exit(1)

# Détection des variants avec GATK
def detect_variants():
    check_file(REFERENCE_FASTA)
    check_file(SORTED_BAM)

    if os.path.exists(VCF_FILE):
        print(f"Fichier VCF déjà existant : {VCF_FILE}, saut de l'étape de détection.")
        return

    print("Détection des variants avec GATK...")
    try:
        subprocess.run(f"gatk HaplotypeCaller -R {REFERENCE_FASTA} -I {SORTED_BAM} -O {VCF_FILE} --ploidy 2", shell=True,
                       check=True)
        print(f"Fichier VCF généré : {VCF_FILE}")
    except subprocess.CalledProcessError:
        print("ERREUR : Échec de la détection des variants.")
        exit(1)

# Analyse des variants avec VCFtools
def analyze_variants():
    check_file(VCF_FILE)

    print("Analyse des variants avec VCFtools...")
    try:
        subprocess.run(f"vcftools --vcf {VCF_FILE} --out analysis_output --freq --depth --site-mean-depth --hist-indel-len",
                       shell=True, check=True)
        print("Analyse des variants terminée.")
    except subprocess.CalledProcessError:
        print("ERREUR : Échec de l'analyse des variants.")
        exit(1)

# Annotation des variants avec SnpEff
def annotate_variants():
    check_file(VCF_FILE)

    if os.path.exists(ANNOTATED_VCF):
        print(f"Fichier VCF annoté déjà existant : {ANNOTATED_VCF}, saut de l'étape d'annotation.")
        return

    print("Annotation des variants avec SnpEff...")
    try:
        subprocess.run(f"snpEff ann -o vcf NC_045512.2 {VCF_FILE} > {ANNOTATED_VCF}", shell=True, check=True)
        print(f"Fichier VCF annoté généré : {ANNOTATED_VCF}")
    except subprocess.CalledProcessError:
        print("ERREUR : Échec de l'annotation des variants.")
        exit(1)

# Visualisation avec IGV
def launch_igv():
    print("Ouvrir IGV pour la visualisation...")
    print("Accès en ligne : https://igv.org/app/")

    if shutil.which("igv"):
        os.system("igv")
        print("IGV lancé en local.")
    else:
        print("IGV n'est pas installé en local. Utilisez la version en ligne.")


def read_vcf(file_path):
    """
    Charge un fichier VCF en DataFrame en ignorant les lignes de commentaires.
    Gère aussi bien les fichiers `.vcf` que `.vcf.gz`.
    """
    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]  # Ignore les lignes de métadonnées

    # Charger le fichier en DataFrame
    df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')

    # Vérification et correction des colonnes
    print("\n✅ Colonnes détectées :", df.columns.tolist())  # Debugging

    if '#CHROM' in df.columns:
        df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)

    return df

def compare_variants(user_variants_file, clinvar_file):
    """
    Compare les variants d'une séquence avec ceux connus dans ClinVar et affiche les résultats.
    """

    # Charger les variants de l'utilisateur
    user_variants = read_vcf(user_variants_file)

    # Vérification des colonnes avant sélection
    required_columns = {'CHROM', 'POS', 'REF', 'ALT'}
    actual_columns = set(user_variants.columns)

    if not required_columns.issubset(actual_columns):
        print("\n❌ ERREUR : Les colonnes attendues ne sont pas toutes présentes dans user_variants.")
        print("🔍 Colonnes disponibles :", user_variants.columns.tolist())
        exit(1)  # Stopper l'exécution

    user_variants = user_variants[['CHROM', 'POS', 'REF', 'ALT']]

    # Charger les variants ClinVar
    clinvar_variants = read_vcf(clinvar_file)

    if not required_columns.issubset(set(clinvar_variants.columns)):
        print("\n❌ ERREUR : Les colonnes attendues ne sont pas toutes présentes dans clinvar_variants.")
        print("🔍 Colonnes disponibles :", clinvar_variants.columns.tolist())
        exit(1)

    clinvar_variants = clinvar_variants[['CHROM', 'POS', 'REF', 'ALT']]

    # Comparer les variants
    common_variants = pd.merge(user_variants, clinvar_variants, on=['CHROM', 'POS', 'REF', 'ALT'], how='inner')

    # Affichage des résultats
    print(f"\n✅ Nombre de variants reconnus dans ClinVar : {len(common_variants)}")

    if not common_variants.empty:
        print("\n🧬 Variants reconnus :")
        print(common_variants.to_string(index=False))  # Affichage formaté

if __name__ == "__main__":
    # generate_fasta_index()
    # generate_fasta_dict()
    # detect_variants()
    # analyze_variants()
    # annotate_variants()
    # launch_igv()
    user_variants_path = "annotated_variants.vcf"
    clinvar_path = "clinvar.vcf.gz"

    compare_variants(user_variants_path, clinvar_path)

