"""Tests for output generation functions."""

import csv
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock

from pangenomeplus.clustering import FamilyStats
from pangenomeplus.compact_ids import CompactIDManager
from pangenomeplus.core import Feature
from pangenomeplus.outputs import (
    generate_all_outputs,
    generate_family_summary,
    generate_presence_absence_matrix,
    generate_roary_compatible_output,
    generate_transformer_format,
    reconstruct_genomic_coordinates,
)


class TestReconstructGenomicCoordinates(unittest.TestCase):
    """Test genomic coordinate reconstruction."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.id_manager = CompactIDManager()

        # Create test features
        self.feature1 = Feature(
            compact_id="P1",
            genome_id="genome_001",
            contig="contig1",
            start=100,
            end=500,
            strand="+",
            sequence="ATGCGT",
            feature_type="P",
            original_id="gene_001",
            metadata={"product": "protein A"},
        )

        self.feature2 = Feature(
            compact_id="P2",
            genome_id="genome_001",
            contig="contig1",
            start=800,
            end=1200,
            strand="+",
            sequence="ATGAAA",
            feature_type="P",
            original_id="gene_002",
            metadata={"product": "protein B"},
        )

        self.feature3 = Feature(
            compact_id="I1",
            genome_id="genome_001",
            contig="contig1",
            start=600,
            end=700,
            strand=".",
            sequence="ATGATG",
            feature_type="I",
            original_id="intergenic_001",
            metadata={"product": "intergenic_region"},
        )

        # Add features to ID manager
        for feature in [self.feature1, self.feature2, self.feature3]:
            self.id_manager.register_feature(feature)

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_coordinate_reconstruction_basic(self):
        """Test basic coordinate reconstruction."""
        # Mock family assignments
        family_assignments = {
            "P": {"P1": "FAM_P1", "P2": "FAM_P2"},
            "I": {"I1": "FAM_I1"},
        }

        family_ids = ["FAM_P1", "FAM_I1", "FAM_P2"]

        coords = reconstruct_genomic_coordinates(
            "genome_001", family_ids, family_assignments, self.id_manager
        )

        # Should be sorted by genomic coordinate
        self.assertEqual(len(coords), 3)
        self.assertEqual(coords[0], ("FAM_P1", 100))  # First by position
        self.assertEqual(coords[1], ("FAM_I1", 600))  # Second by position
        self.assertEqual(coords[2], ("FAM_P2", 800))  # Third by position

    def test_coordinate_reconstruction_missing_family(self):
        """Test coordinate reconstruction with missing family."""
        family_assignments = {"P": {"P1": "FAM_P1"}}
        family_ids = ["FAM_P1", "FAM_MISSING"]

        coords = reconstruct_genomic_coordinates(
            "genome_001", family_ids, family_assignments, self.id_manager
        )

        # Should only find the existing family
        self.assertEqual(len(coords), 1)
        self.assertEqual(coords[0], ("FAM_P1", 100))

    def test_coordinate_reconstruction_wrong_genome(self):
        """Test coordinate reconstruction for wrong genome."""
        family_assignments = {"P": {"P1": "FAM_P1"}}
        family_ids = ["FAM_P1"]

        coords = reconstruct_genomic_coordinates(
            "genome_999", family_ids, family_assignments, self.id_manager
        )

        # Should find no families for wrong genome
        self.assertEqual(len(coords), 0)


class TestTransformerFormat(unittest.TestCase):
    """Test transformer format generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.id_manager = CompactIDManager()

        # Create test features
        features = [
            Feature(
                compact_id="P1",
                genome_id="E_coli_001",
                contig="contig1",
                start=100,
                end=500,
                strand="+",
                sequence="ATGCGT",
                feature_type="P",
                original_id="gene_001",
                metadata={"product": "protein A"},
            ),
            Feature(
                compact_id="P2",
                genome_id="E_coli_001",
                contig="contig1",
                start=800,
                end=1200,
                strand="+",
                sequence="ATGAAA",
                feature_type="P",
                original_id="gene_002",
                metadata={"product": "protein B"},
            ),
            Feature(
                compact_id="P3",
                genome_id="E_coli_002",
                contig="contig1",
                start=200,
                end=600,
                strand="+",
                sequence="ATGTTT",
                feature_type="P",
                original_id="gene_003",
                metadata={"product": "protein C"},
            ),
        ]

        for feature in features:
            self.id_manager.register_feature(feature)

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_transformer_format_generation(self):
        """Test transformer format generation."""
        family_assignments = {
            "P": {
                "P1": "FAM_P1",
                "P2": "FAM_P2",
                "P3": "FAM_P1",
            }  # P3 in same family as P1
        }

        output_file = generate_transformer_format(
            family_assignments, self.id_manager, self.temp_dir
        )

        # Check file was created
        self.assertTrue(os.path.exists(output_file))

        # Check file content
        with open(output_file, "r") as f:
            lines = f.readlines()

        # Should have 2 genomes
        self.assertEqual(len(lines), 2)

        # Check genome ordering (alphabetical)
        lines.sort()  # Sort to ensure consistent order
        line1_parts = lines[0].strip().split()
        line2_parts = lines[1].strip().split()

        self.assertEqual(line1_parts[0], "E_coli_001")
        self.assertEqual(line2_parts[0], "E_coli_002")

        # Check that families are space-separated tokens
        self.assertTrue(len(line1_parts) >= 2)  # Genome + at least 1 family
        self.assertTrue(len(line2_parts) >= 2)  # Genome + at least 1 family

    def test_transformer_format_empty_input(self):
        """Test transformer format with empty input."""
        family_assignments = {}

        output_file = generate_transformer_format(
            family_assignments, self.id_manager, self.temp_dir
        )

        # File should be created but empty
        self.assertTrue(os.path.exists(output_file))

        with open(output_file, "r") as f:
            content = f.read().strip()

        self.assertEqual(content, "")


class TestPresenceAbsenceMatrix(unittest.TestCase):
    """Test presence/absence matrix generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.id_manager = CompactIDManager()

        # Create test features
        features = [
            Feature(
                compact_id="P1",
                genome_id="genome_A",
                contig="contig1",
                start=100,
                end=500,
                strand="+",
                sequence="ATGCGT",
                feature_type="P",
                original_id="gene_001",
                metadata={"product": "protein A"},
            ),
            Feature(
                compact_id="P2",
                genome_id="genome_A",
                contig="contig1",
                start=600,
                end=900,
                strand="+",
                sequence="ATGAAA",
                feature_type="P",
                original_id="gene_002",
                metadata={"product": "protein B"},
            ),
            Feature(
                compact_id="P3",
                genome_id="genome_B",
                contig="contig1",
                start=100,
                end=500,
                strand="+",
                sequence="ATGCGT",
                feature_type="P",
                original_id="gene_003",
                metadata={"product": "protein A"},
            ),
        ]

        for feature in features:
            self.id_manager.register_feature(feature)

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_presence_absence_matrix_generation(self):
        """Test presence/absence matrix generation."""
        family_assignments = {
            "P": {"P1": "FAM_P1", "P2": "FAM_P2", "P3": "FAM_P1"}  # P3 same as P1
        }

        output_file = generate_presence_absence_matrix(
            family_assignments, self.id_manager, self.temp_dir
        )

        # Check file was created
        self.assertTrue(os.path.exists(output_file))

        # Parse CSV and check structure
        with open(output_file, "r", newline="") as f:
            reader = csv.reader(f)
            rows = list(reader)

        # Check header
        header = rows[0]
        self.assertEqual(header[0], "Genome")
        self.assertIn("FAM_P1", header)
        self.assertIn("FAM_P2", header)

        # Check data rows
        self.assertEqual(len(rows), 3)  # Header + 2 genomes

        # Find genome rows
        genome_rows = {row[0]: row[1:] for row in rows[1:]}

        # Check presence/absence values
        self.assertIn("genome_A", genome_rows)
        self.assertIn("genome_B", genome_rows)

        # Each row should have same number of columns as families
        family_count = len(header) - 1
        for genome_id, values in genome_rows.items():
            self.assertEqual(len(values), family_count)
            for value in values:
                self.assertIn(int(value), [0, 1])  # Only 0 or 1


class TestFamilySummary(unittest.TestCase):
    """Test family summary generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_family_summary_generation(self):
        """Test family summary generation."""
        # Mock family stats
        family_stats = {
            "P": {
                "FAM_P1": FamilyStats(
                    family_id="FAM_P1",
                    feature_type="P",
                    member_count=5,
                    genome_count=10,
                    representative_id="P1",
                    classification="core",
                ),
                "FAM_P2": FamilyStats(
                    family_id="FAM_P2",
                    feature_type="P",
                    member_count=2,
                    genome_count=3,
                    representative_id="P2",
                    classification="accessory",
                ),
            }
        }

        output_file = generate_family_summary(
            family_stats, self.temp_dir, total_genomes=10
        )

        # Check file was created
        self.assertTrue(os.path.exists(output_file))

        # Parse TSV and check content
        with open(output_file, "r", newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            rows = list(reader)

        # Check header
        header = rows[0]
        expected_columns = [
            "Family_ID",
            "Feature_Type",
            "Size",
            "Genome_Count",
            "Genome_Percentage",
            "Classification",
            "Representative",
        ]
        self.assertEqual(header, expected_columns)

        # Check data rows
        self.assertEqual(len(rows), 3)  # Header + 2 families

        # Check family data
        for row in rows[1:]:
            family_id = row[0]
            if family_id == "FAM_P1":
                self.assertEqual(row[1], "P")  # Feature type
                self.assertEqual(row[2], "5")  # Size
                self.assertEqual(row[3], "10")  # Genome count
                self.assertEqual(row[4], "100.0%")  # Percentage
                self.assertEqual(row[5], "core")  # Classification
            elif family_id == "FAM_P2":
                self.assertEqual(row[1], "P")
                self.assertEqual(row[2], "2")
                self.assertEqual(row[3], "3")
                self.assertEqual(row[4], "30.0%")
                self.assertEqual(row[5], "accessory")


class TestRoaryCompatibleOutput(unittest.TestCase):
    """Test Roary-compatible output generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.id_manager = CompactIDManager()

        # Create test features
        features = [
            Feature(
                compact_id="P1",
                genome_id="genome_A",
                contig="contig1",
                start=100,
                end=500,
                strand="+",
                sequence="ATGCGT",
                feature_type="P",
                original_id="gene_001",
                metadata={"product": "hypothetical protein"},
            ),
            Feature(
                compact_id="P2",
                genome_id="genome_B",
                contig="contig1",
                start=100,
                end=500,
                strand="+",
                sequence="ATGCGT",
                feature_type="P",
                original_id="gene_002",
                metadata={"product": "DNA polymerase"},
            ),
        ]

        for feature in features:
            self.id_manager.register_feature(feature)

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_roary_compatible_generation(self):
        """Test Roary-compatible output generation."""
        family_assignments = {"P": {"P1": "FAM_P1", "P2": "FAM_P1"}}

        output_file = generate_roary_compatible_output(
            family_assignments, self.id_manager, self.temp_dir
        )

        # Check file was created
        self.assertTrue(os.path.exists(output_file))

        # Parse CSV and check content
        with open(output_file, "r", newline="") as f:
            reader = csv.reader(f)
            rows = list(reader)

        # Check header
        header = rows[0]
        expected_columns = [
            "Gene",
            "Non-unique Gene name",
            "Annotation",
            "No. isolates",
            "No. sequences",
        ]
        self.assertEqual(header, expected_columns)

        # Check data rows
        self.assertEqual(len(rows), 2)  # Header + 1 family

        family_row = rows[1]
        self.assertEqual(family_row[0], "FAM_P1")  # Family ID
        self.assertIn(family_row[1], ["gene_001", "gene_002"])  # Original ID
        self.assertIn("protein", family_row[2].lower())  # Product annotation
        self.assertEqual(family_row[3], "2")  # Number of genomes
        self.assertEqual(family_row[4], "2")  # Number of sequences


class TestGenerateAllOutputs(unittest.TestCase):
    """Test comprehensive output generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.id_manager = CompactIDManager()

        # Create test features
        features = [
            Feature(
                compact_id="P1",
                genome_id="test_genome",
                contig="contig1",
                start=100,
                end=500,
                strand="+",
                sequence="ATGCGT",
                feature_type="P",
                original_id="gene_001",
                metadata={"product": "protein A"},
            )
        ]

        for feature in features:
            self.id_manager.register_feature(feature)

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_generate_all_outputs(self):
        """Test generation of all output formats."""
        family_assignments = {"P": {"P1": "FAM_P1"}}

        family_stats = {
            "P": {
                "FAM_P1": FamilyStats(
                    family_id="FAM_P1",
                    feature_type="P",
                    member_count=1,
                    genome_count=1,
                    representative_id="P1",
                    classification="singleton",
                )
            }
        }

        output_files = generate_all_outputs(
            family_assignments, family_stats, self.id_manager, self.temp_dir, 1
        )

        # Check all expected outputs were generated
        expected_formats = ["transformer", "matrix", "summary", "roary"]
        for format_name in expected_formats:
            self.assertIn(format_name, output_files)
            self.assertTrue(os.path.exists(output_files[format_name]))

        # Check directory structure
        expected_dirs = ["transformer", "matrices", "roary"]
        for dir_name in expected_dirs:
            dir_path = os.path.join(self.temp_dir, dir_name)
            self.assertTrue(os.path.exists(dir_path))


if __name__ == "__main__":
    unittest.main()
