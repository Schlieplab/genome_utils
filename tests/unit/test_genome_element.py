"""Tests for the GenomeElement base class."""

import pytest
from unittest.mock import Mock
from Bio.Seq import Seq
from GenomeUtils.Genome import GenomeElement, Locus

class ConcreteGenomeElement(GenomeElement):
    """Concrete implementation of GenomeElement for testing purposes."""
    
    def __init__(self, *args, sequence="ATCGATCGATCG", **kwargs):
        super().__init__(*args, **kwargs)
        self._sequence = sequence
    
    @property
    def sequence(self) -> Seq:
        """Return the sequence of this element."""
        return Seq(self._sequence)


class TestGenomeElement:
    """Test cases for the GenomeElement base class."""

    def test_genome_element_creation(self):
        """Test basic GenomeElement creation."""
        locus = Locus("chr1", 100, 200, "+")
        genome_mock = Mock()
        parent_mock = Mock()
        
        element = ConcreteGenomeElement(
            id="test_element",
            locus=locus,
            parent=parent_mock,
            genome=genome_mock,
            custom_attr="custom_value"
        )
        
        assert element.id == "test_element"
        assert element.locus == locus
        assert element._parent == parent_mock
        assert element._genome == genome_mock
        assert element.custom_attr == "custom_value"
        assert element._children == []

    def test_genome_element_creation_minimal(self):
        """Test GenomeElement creation with minimal parameters."""
        locus = Locus("chr2", 300, 400)
        
        element = ConcreteGenomeElement(id="minimal", locus=locus)
        
        assert element.id == "minimal"
        assert element.locus == locus
        assert element._parent is None
        assert element._genome is None
        assert element._children == []

    def test_genome_element_properties(self):
        """Test that properties correctly access locus attributes."""
        locus = Locus("chr3", 500, 600, "-")
        element = ConcreteGenomeElement(id="test", locus=locus)
        
        assert element.chr == "chr3"
        assert element.start == 500
        assert element.end == 600
        assert element.strand == "-"

    def test_genome_element_length(self):
        """Test that __len__ returns locus length."""
        locus = Locus("chr1", 100, 200)  # Length should be 101
        element = ConcreteGenomeElement(id="test", locus=locus)
        
        assert len(element) == 101

    def test_genome_element_repr(self):
        """Test string representation."""
        locus = Locus("chr1", 100, 200, "+")
        element = ConcreteGenomeElement(id="test_elem", locus=locus)
        
        expected = "ConcreteGenomeElement(id='test_elem', locus=Locus(chr1:100-200, strand=+))"
        assert repr(element) == expected

    def test_genome_element_equality(self):
        """Test equality comparison based on ID and locus."""
        locus1 = Locus("chr1", 100, 200)
        locus2 = Locus("chr2", 300, 400)  # Different locus
        
        element1 = ConcreteGenomeElement(id="same_id", locus=locus1)
        element2 = ConcreteGenomeElement(id="same_id", locus=locus1)  # Same ID and locus
        element3 = ConcreteGenomeElement(id="same_id", locus=locus2)  # Same ID, different locus
        element4 = ConcreteGenomeElement(id="different_id", locus=locus1)  # Different ID, same locus
        
        assert element1 == element2  # Same ID and locus
        assert element1 != element3  # Same ID, different locus
        assert element1 != element4  # Different ID, same locus
        
        # Test type checking
        assert element1 != "not_a_genome_element"
        assert element1 != None

    def test_genome_element_hash(self):
        """Test that GenomeElement objects are hashable and work in sets/dicts."""
        locus1 = Locus("chr1", 100, 200)
        locus2 = Locus("chr2", 300, 400)
        
        element1 = ConcreteGenomeElement(id="element1", locus=locus1)
        element2 = ConcreteGenomeElement(id="element1", locus=locus1)  # Same as element1
        element3 = ConcreteGenomeElement(id="element2", locus=locus2)  # Different
        
        # Test that equal elements have same hash
        assert hash(element1) == hash(element2)
        
        # Test that they work in sets (duplicates removed)
        element_set = {element1, element2, element3}
        assert len(element_set) == 2  # element1 and element2 are equal
        
        # Test that they work as dict keys
        element_dict = {element1: "value1", element3: "value2"}
        assert len(element_dict) == 2
        assert element_dict[element2] == "value1"  # element2 equals element1

    def test_genome_element_kwargs_attributes(self):
        """Test that additional kwargs become attributes."""
        locus = Locus("chr1", 100, 200)
        
        element = ConcreteGenomeElement(
            id="test",
            locus=locus,
            gene_type="protein_coding",
            score=95.5,
            tags=["tag1", "tag2"],
            metadata={"key": "value"}
        )
        
        assert element.gene_type == "protein_coding"
        assert element.score == 95.5
        assert element.tags == ["tag1", "tag2"]
        assert element.metadata == {"key": "value"}

    def test_genome_element_children_list(self):
        """Test that children list is properly initialized and accessible."""
        locus = Locus("chr1", 100, 200)
        element = ConcreteGenomeElement(id="parent", locus=locus)
        
        # Should start empty
        assert element._children == []
        assert len(element._children) == 0
        
        # Should be able to add children
        child_mock = Mock()
        element._children.append(child_mock)
        assert len(element._children) == 1
        assert element._children[0] == child_mock

    def test_genome_element_parent_reference(self):
        """Test parent reference functionality."""
        locus = Locus("chr1", 100, 200)
        parent_mock = Mock()
        
        element = ConcreteGenomeElement(id="child", locus=locus, parent=parent_mock)
        
        assert element._parent == parent_mock

    def test_genome_element_genome_reference(self):
        """Test genome reference functionality."""
        locus = Locus("chr1", 100, 200)
        genome_mock = Mock()
        
        element = ConcreteGenomeElement(id="element", locus=locus, genome=genome_mock)
        
        assert element._genome == genome_mock

    def test_genome_element_abstract_nature(self):
        """Test that GenomeElement can be instantiated despite being abstract."""
        # GenomeElement is marked as ABC but doesn't have abstract methods
        # so it should be instantiable
        locus = Locus("chr1", 100, 200)
        
        # This should work
        element = ConcreteGenomeElement(id="test", locus=locus)
        assert isinstance(element, GenomeElement)

    def test_genome_element_locus_immutability(self):
        """Test that the locus reference is maintained correctly."""
        locus = Locus("chr1", 100, 200)
        element = ConcreteGenomeElement(id="test", locus=locus)
        
        # The locus should be the same object
        assert element.locus is locus
        
        # And should maintain its properties
        assert element.locus.chr == "chr1"
        assert element.locus.start == 100
        assert element.locus.end == 200

    def test_genome_element_with_none_values(self):
        """Test GenomeElement creation with None values for optional parameters."""
        locus = Locus("chr1", 100, 200)
        
        element = ConcreteGenomeElement(
            id="test",
            locus=locus,
            parent=None,
            genome=None
        )
        
        assert element._parent is None
        assert element._genome is None
        assert element._children == []


