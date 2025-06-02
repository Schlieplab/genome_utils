from dataclasses import dataclass, field
from typing import List, Optional, Tuple, TYPE_CHECKING
import re

if TYPE_CHECKING:
    from .transcript import Transcript
    from .exon import Exon

@dataclass
class Site:
    sequence: Optional[str] = None
    chromosomal_position: Optional[str] = None  # Expected format: "chr:start-end:strand"

    def __len__(self) -> int:
        return len(self.sequence) if self.sequence else 0

    @staticmethod
    def _parse_chrom_pos(pos_str: Optional[str]) -> Tuple[Optional[str], Optional[int], Optional[int], Optional[str]]:
        if not pos_str:
            return None, None, None, None
        
        match = re.match(r"([^:]+):(\d+)-(\d+):([+-])", pos_str)
        if match:
            chrom, start_str, end_str, strand_char = match.groups()
            try:
                return chrom, int(start_str), int(end_str), strand_char
            except ValueError:
                return None, None, None, None # Or log error
        return None, None, None, None # Or log error

    def to_gtf_attributes(self, additional_attrs: Optional[dict] = None) -> str:
        """Helper to create a GTF attribute string."""
        attrs = {}
        if self.sequence:
            attrs["sequence"] = self.sequence
        if self.chromosomal_position:
            attrs["original_chrom_pos"] = self.chromosomal_position
        
        if additional_attrs:
            attrs.update(additional_attrs)
        
        return "; ".join([f'{k} "{v}"' for k, v in attrs.items()])

@dataclass
class RepeatedSite(Site):
    parent_target_id: str = ""
    ddg_to_parent: Optional[float] = None

    def to_gtf_record(self, source: str, feature_type: str = "repeated_site") -> Optional[str]:
        chrom, start, end, strand = Site._parse_chrom_pos(self.chromosomal_position)
        if not all([chrom, start, end, strand]):
            return None # Cannot form valid GTF

        attributes_dict = {
            "parent_target_id": self.parent_target_id,
        }
        if self.ddg_to_parent is not None:
            attributes_dict["ddG_to_parent"] = f"{self.ddg_to_parent:.2f}"
        
        attribute_str = self.to_gtf_attributes(attributes_dict)
        
        # GTF fields: seqname, source, feature, start, end, score, strand, frame, attributes
        return f"{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attribute_str}"

@dataclass
class CandidateTarget(Site): 
    id: str = "" 
    gene_id: Optional[str] = None
    transcripts: List['Transcript'] = field(default_factory=list)
    exons: List['Exon'] = field(default_factory=list)
    dG_binding: Optional[float] = None
    oligo_homodimer_dG: Optional[float] = None
    pedersen_steady_state: Optional[float] = None
    repeated_sites: List[RepeatedSite] = field(default_factory=list)

    def add_repeated_site(self, site: RepeatedSite):
        if site.parent_target_id != self.id:
            # Or raise an error, or auto-correct parent_target_id
            print(f"Warning: Adding repeated site with parent_target_id {site.parent_target_id} to CandidateTarget {self.id}")
            site.parent_target_id = self.id # Auto-correct
        self.repeated_sites.append(site)

    def filter_repeated_sites(self, ddg_threshold: float):
        """Filters self.repeated_sites in-place, keeping sites with ddG <= threshold."""
        self.repeated_sites = [
            rs for rs in self.repeated_sites 
            if rs.ddg_to_parent is not None and rs.ddg_to_parent <= ddg_threshold
        ]

    def to_gtf_record(self, source: str, feature_type: str = "candidate_target") -> Optional[str]:
        chrom, start, end, strand = Site._parse_chrom_pos(self.chromosomal_position)
        if not all([chrom, start, end, strand]):
            return None # Cannot form valid GTF

        attributes_dict = {
            "target_id": self.id,
        }
        if self.gene_id:
            attributes_dict["gene_id"] = self.gene_id
        if self.dG_binding is not None:
            attributes_dict["dG_binding"] = f"{self.dG_binding:.2f}"
        if self.oligo_homodimer_dG is not None:
            attributes_dict["oligo_homodimer_dG"] = f"{self.oligo_homodimer_dG:.2f}"
        if self.pedersen_steady_state is not None:
             attributes_dict["pedersen_steady_state"] = f"{self.pedersen_steady_state:.3f}"
        attributes_dict["num_repeated_sites"] = len(self.repeated_sites)
        
        # Could add transcript_ids, exon_ids if desired
        # attributes_dict["transcript_ids"] = ",".join([t.transcript_id for t in self.transcripts])

        attribute_str = self.to_gtf_attributes(attributes_dict)
        
        return f"{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attribute_str}" 