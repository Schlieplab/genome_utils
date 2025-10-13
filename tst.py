import GenomeUtils.Genome as Genome
import GenomeUtils.Downloaders
from pathlib import Path

downloader = GenomeUtils.Downloaders.EnsemblGenomeDownloader(
    assembly_id="GRCh38",
    ensembl_release=109,
    species="homo_sapiens",
    genomes_root_dir=Path("./data/genomes"),
)

files = downloader.download()
print(files)

builder = Genome.GenomeBuilder(
    id="GRCh38",
    species="homo_sapiens",
    name="GRCh38",
    genomes_root_dir=Path("./data/genomes"),
)

builder.with_dna_fasta(files["dna"])
builder.with_cdna_fasta(files["cdna"])
builder.with_gtf_file(files["annotation"])
genome = builder.build()
print(genome)