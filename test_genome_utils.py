import os
from pathlib import Path

from . import GgetEnsemblGenomeDownloader, GenomeBuilder


def main():
    """
    Tests the full pipeline:
    1. Download genome files using GgetDownloader.
    2. Build a Genome object from the downloaded files.
    3. Print a summary.
    """
    # Define parameters for the test
    ensembl_release = 114 
    assembly_id = "GRCh38"
    species = "homo_sapiens"
    
    # 1. Download the files
    print(f"Starting download for {species} (release {ensembl_release})...")
    downloader = GgetEnsemblGenomeDownloader(assembly_id=assembly_id, ensembl_release=ensembl_release, species=species)
    downloaded_files = downloader.download()
    
    dna_path = downloaded_files['dna']
    cdna_path = downloaded_files['cdna']
    gtf_path = downloaded_files['annotation']
    
    # 2. Build the Genome object from the downloaded files
    print("\nBuilding Genome object...")
    builder = GenomeBuilder(id=assembly_id, species=species, name=f"{species} Genome (release {ensembl_release})")
    genome = (
        builder.with_dna_fasta(dna_path)
        .with_cdna_fasta(cdna_path)
        .with_gtf_file(gtf_path)
        .build()
    )
    
    # 3. Print summary
    print("\n--- Genome Build Summary ---")
    print(f"Genome: {genome.name} (ID: {genome.id})")
    print(f"Species: {genome.species}")
    print(f"Number of chromosomes: {len(genome.chromosomes)}")
    print(f"Chromosome names: {', '.join(genome.chromosomes.keys())}")
    
    # Convert generators to lists to allow multiple iterations and len()
    all_genes = list(genome.genes)
    all_transcripts = list(genome.transcripts)

    total_genes = len(all_genes)
    total_transcripts = len(all_transcripts)
    
    print(f"Total genes found: {total_genes:,}")
    print(f"Total transcripts found: {total_transcripts:,}")

    if total_genes > 0:
        transcripts_per_gene = [len(gene.transcripts) for gene in all_genes]
        avg_transcripts = sum(transcripts_per_gene) / total_genes
        print(f"Transcripts per gene: Min={min(transcripts_per_gene)}, Max={max(transcripts_per_gene)}, Avg={avg_transcripts:.2f}")

    # --- Detailed Gene Type Comparison ---
    print("\n--- Gene Type Comparison (Primary Assembly) ---")
    
    # Numbers from Ensembl GRCh38.p14
    ensembl_stats = {
        "coding_genes": 19871,
        "non_coding_genes": 42126,
        "small_non_coding": 4866,
        "long_non_coding": 35044,
        "misc_non_coding": 2216,
        "pseudogenes": 15198,
        "gene_transcripts": 387954,
    }

    counts = {
        "coding_genes": 0, "non_coding_genes": 0, "small_non_coding": 0,
        "long_non_coding": 0, "misc_non_coding": 0, "pseudogenes": 0,
    }

    # Classify genes based on their biotype from the GTF attributes
    for gene in all_genes:
        # The 'gene_biotype' is usually a list with one item
        biotype = getattr(gene, 'gene_biotype', ['other'])[0]

        if biotype == 'protein_coding':
            counts['coding_genes'] += 1
        elif 'pseudogene' in biotype:
            counts['pseudogenes'] += 1
        else:
            counts['non_coding_genes'] += 1
            # Further classify non-coding genes
            if 'lncRNA' in biotype:
                counts['long_non_coding'] += 1
            elif biotype in ['snRNA', 'snoRNA', 'miRNA', 'rRNA', 'scRNA', 'scaRNA', 'sRNA', 'ribozyme']:
                counts['small_non_coding'] += 1
            else:
                counts['misc_non_coding'] += 1
    
    def print_stat(label, local_count, ensembl_count):
        print(f"{label:<25} {local_count:>10,} (Ensembl: {ensembl_count:>10,})")

    print_stat("Coding genes", counts['coding_genes'], ensembl_stats['coding_genes'])
    print_stat("Non-coding genes (total)", counts['non_coding_genes'], ensembl_stats['non_coding_genes'])
    print_stat("  Small non-coding", counts['small_non_coding'], ensembl_stats['small_non_coding'])
    print_stat("  Long non-coding", counts['long_non_coding'], ensembl_stats['long_non_coding'])
    print_stat("  Misc non-coding", counts['misc_non_coding'], ensembl_stats['misc_non_coding'])
    print_stat("Pseudogenes", counts['pseudogenes'], ensembl_stats['pseudogenes'])
    print_stat("Gene transcripts", total_transcripts, ensembl_stats['gene_transcripts'])
    
    # Verify by fetching a well-known gene
    try:
        gene = genome.gene_by_id("ENSG00000139618") # Example: BRCA2 gene
        print(f"\nSuccessfully fetched example gene: {gene.name} (ID: {gene.id})")
        print(f"Location: Chromosome {gene.chromosome_id}, {gene.start}-{gene.end} ({gene.strand})")
        print(f"Number of transcripts: {len(gene.transcripts)}")
        
        if gene.transcripts:
            first_transcript = gene.transcripts[0]
            print(f"  - First transcript ({first_transcript.id}) has {len(first_transcript.exons)} exons.")

    except (ValueError, KeyError) as e:
        print(f"\nCould not fetch example gene. This might be expected if using a minimal test set. Error: {e}")


if __name__ == "__main__":
    main() 