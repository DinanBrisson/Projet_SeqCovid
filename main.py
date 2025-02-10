from FastQ.fastq_proc import FastqProcessor

if __name__ == "__main__":
    FASTQ_FILE = "data/FastQ/SRR32230015.fastq"
    processor = FastqProcessor(FASTQ_FILE)
    processor.process(clean_reads=True, assemble=True)