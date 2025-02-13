import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


class FastqProcessor:
    def __init__(self, fastq_file, output_dir="data/FastQ"):
        self.fastq_file = fastq_file
        self.output_dir = output_dir
        self.cleaned_fastq_r1 = None
        self.cleaned_fastq_r2 = None
        self.assembly_output_dir = os.path.join(output_dir, "assembly")
        self.r1_fastq = os.path.join(output_dir, "reads_R1.fastq")
        self.r2_fastq = os.path.join(output_dir, "reads_R2.fastq")

        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.assembly_output_dir, exist_ok=True)

    def check_fastq_head(self, num_lines=5):
        print(f"Checking first {num_lines} lines of {self.fastq_file}...\n")
        with open(self.fastq_file, "r") as file:
            for i in range(num_lines):
                print(file.readline().strip())

    def run_fastqc(self):
        print("Running FastQC analysis...")
        cmd = ["fastqc", self.fastq_file, "-o", self.output_dir]
        try:
            subprocess.run(cmd, check=True)
            print(f"FastQC completed. Results saved in {self.output_dir}/")
        except subprocess.CalledProcessError as e:
            print(f"FastQC failed: {e}")

    def separate_paired_reads(self):
        """Sépare les reads appariés en deux fichiers distincts"""
        print("Separating Paired-End reads into two files...")
        cmd = [
            "seqtk", "seq", "-1", self.fastq_file, ">", self.r1_fastq
        ]
        subprocess.run(" ".join(cmd), shell=True, check=True)

        cmd = [
            "seqtk", "seq", "-2", self.fastq_file, ">", self.r2_fastq
        ]
        subprocess.run(" ".join(cmd), shell=True, check=True)

        print(f"Reads separated into {self.r1_fastq} and {self.r2_fastq}")

    def clean_fastq(self, min_quality=20, min_length=50):
        print(f"Cleaning FASTQ files with Trimmomatic (Min Quality: {min_quality}, Min Length: {min_length})...")

        self.cleaned_fastq_r1 = os.path.join(self.output_dir, "cleaned_R1.fastq.gz")
        self.cleaned_fastq_r2 = os.path.join(self.output_dir, "cleaned_R2.fastq.gz")

        cmd = [
            "trimmomatic", "PE", "-phred33",
            self.r1_fastq, self.r2_fastq,
            self.cleaned_fastq_r1, "/dev/null",
            self.cleaned_fastq_r2, "/dev/null",
            f"LEADING:{min_quality}", f"TRAILING:{min_quality}",
            f"SLIDINGWINDOW:4:{min_quality}", f"MINLEN:{min_length}"
        ]

        try:
            subprocess.run(cmd, check=True)
            print(
                f"Cleaning completed. Cleaned FASTQ files saved as {self.cleaned_fastq_r1} and {self.cleaned_fastq_r2}")
        except subprocess.CalledProcessError as e:
            print(f"Trimmomatic failed: {e}")

    def extract_statistics(self, max_reads=100000):
        read_lengths = []
        mean_quality_scores = []

        try:
            with open(self.fastq_file, "r") as f:
                i = 0
                for line in f:
                    i += 1
                    if i % 4 == 2:
                        read_lengths.append(len(line.strip()))
                    elif i % 4 == 0:
                        quality_scores = [ord(char) - 33 for char in line.strip()]
                        if quality_scores:
                            mean_quality_scores.append(np.mean(quality_scores))

                    if len(read_lengths) >= max_reads:
                        break
        except FileNotFoundError:
            print(f"Error: File {self.fastq_file} not found!")
            return None

        return np.array(read_lengths), np.array(mean_quality_scores)

    def plot_statistics(self):
        data = self.extract_statistics()
        if data is None:
            return

        read_lengths, mean_quality_scores = data
        if len(read_lengths) == 0 or len(mean_quality_scores) == 0:
            print("No valid reads found for analysis.")
            return

        df_lengths = pd.DataFrame({"Longueurs de reads": read_lengths})
        df_quality = pd.DataFrame({"Scores de qualité moyens": mean_quality_scores})

        print("\nFastQ Read Length Statistics:")
        print(df_lengths.describe())
        print("\nFastQ Quality Score Statistics:")
        print(df_quality.describe())

        fig, axs = plt.subplots(1, 2, figsize=(14, 5))

        sns.histplot(read_lengths, bins=50, kde=True, ax=axs[0])
        axs[0].set_title("Distribution des longueurs de reads")
        axs[0].set_xlabel("Longueur des reads")
        axs[0].set_ylabel("Fréquence")

        sns.histplot(mean_quality_scores, bins=50, kde=True, ax=axs[1])
        axs[1].set_title("Distribution des scores de qualité moyens (Phred33)")
        axs[1].set_xlabel("Score de qualité moyen")
        axs[1].set_ylabel("Fréquence")

        plt.tight_layout()
        plt.show()

    def run_spades_assembly(self):
        print("Running SPAdes de novo assembly...")

        cmd = ["spades.py",
               "--isolate", "-1",
               self.cleaned_fastq_r1, "-2",
               self.cleaned_fastq_r2, "-o",
               self.assembly_output_dir]

        try:
            subprocess.run(cmd, check=True)
            print(f"SPAdes assembly completed. Results in {self.assembly_output_dir}")
        except subprocess.CalledProcessError as e:
            print(f"SPAdes failed: {e}")

    def process(self, clean_reads=False, assemble=False):
        self.check_fastq_head()
        self.run_fastqc()
        self.separate_paired_reads()
        if clean_reads:
            self.clean_fastq()
            self.plot_statistics()
        if assemble:
            self.run_spades_assembly()

if __name__ == "__main__":
    FASTQ_FILE = "data/FastQ/SRR32230015.fastq" # remplacer avec un chemin de fichier fastq
    processor = FastqProcessor(FASTQ_FILE)
    processor.process(clean_reads=True, assemble=True)
