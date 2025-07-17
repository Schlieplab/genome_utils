
import os
from pathlib import Path
from genome.builder import GenomeBuilder

def create_dummy_files():
    """Creates dummy FASTA and GTF files for demonstration."""
    # Create dummy directory
    dummy_dir = Path("dummy_data")
    dummy_dir.mkdir(exist_ok=True)

    # DNA FASTA
    dna_fasta_path = dummy_dir / "dna.fa"
    with open(dna_fasta_path, "w") as f:
        f.write(">chr1\n")
        f.write("AGCATGATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n")

    # cDNA FASTA
    cdna_fasta_path = dummy_dir / "cdna.fa"
    with open(cdna_fasta_path, "w") as f:
        f.write(">TRANSCRIPT001\n")
        f.write("CATGATGCATGCATGCATGCATGCATGC\n")

    # GTF
    gtf_path = dummy_dir / "annotations.gtf"
    with open(gtf_path, "w") as f:
        f.write('chr1\tdummy\tgene\t5\t35\t.\t+\t.\tgene_id "GENE001"; gene_name "MYGENE"; gene_source "ensembl";\n')
        f.write('chr1\tdummy\ttranscript\t5\t35\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANSCRIPT001"; transcript_biotype "protein_coding";\n')
        f.write('chr1\tdummy\texon\t5\t15\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANSCRIPT001"; exon_id "EXON001.1";\n')
        f.write('chr1\tdummy\texon\t25\t35\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANSCRIPT001"; exon_id "EXON001.2";\n')

    return dna_fasta_path, cdna_fasta_path, gtf_path

def main():
    """Demonstrates how to use the GenomeBuilder."""
    dna_path, cdna_path, gtf_path = create_dummy_files()

    # 1. Initialize the builder
    builder = GenomeBuilder(id="hg38_toy", species="Homo sapiens", name="Toy Genome")

    # 2. Chain the methods to build the genome
    genome = (
        builder.with_dna_fasta(dna_path)
        .with_cdna_fasta(cdna_path)
        .with_gtf_file(gtf_path)
        .build()
    )

    # 3. Access the constructed genome and its features
    print(f"Successfully built genome: {genome.name}")
    print(f"Number of chromosomes: {len(genome.chromosomes)}")

    # Access a gene and its extra attributes
    my_gene = genome.gene_by_id("GENE001")
    print(f"\nFound Gene: {my_gene.name} (ID: {my_gene.id})")
    print(f"Gene's source from GTF: {my_gene.gene_source[0]}")

    # Access a transcript and its extra attributes
    my_transcript = genome.transcript_by_id("TRANSCRIPT001")
    print(f"\nFound Transcript: {my_transcript.id}")
    print(f"Transcript's biotype from GTF: {my_transcript.transcript_biotype[0]}")

    # Clean up dummy files
    for file_path in [dna_path, cdna_path, gtf_path, gtf_path.with_suffix(".db")]:
        if os.path.exists(file_path):
            os.remove(file_path)
    if os.path.exists(dna_path.parent):
        os.rmdir(dna_path.parent)


if __name__ == "__main__":
    main() 