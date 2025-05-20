from dataclasses import dataclass
from typing import List, Optional

from .transcript import Transcript
from .exon import Exon

@dataclass
class Site:
    sequence: Optional[str] = None
    chromosomal_position: Optional[str] = None
    
    def __len__(self):
        return len(self.sequence) if self.sequence else 0

@dataclass
class TargetSite(Site):
    gene_id: Optional[str] = None
    transcripts: Optional[List[Transcript]] = None
    exons: Optional[List[Exon]] = None
    dG: Optional[float] = None
    oligo_dG: Optional[float] = None  # dG of the homodimer of the oligo(reverse complement of the target site)
    pedersen_steady_state: Optional[float] = None  # percentage of target gene in steady state of the Pedersen model

    def __post_init__(self):
        if self.transcripts is None:
            self.transcripts = []
        if self.exons is None:
            self.exons = [] 