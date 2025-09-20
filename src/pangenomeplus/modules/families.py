"""Gene family assignment and paralog resolution module."""

import logging
from pathlib import Path
from typing import Dict, List, Any, Tuple
import json

from pangenomeplus.core.types import FeatureType, GenomeFeature


logger = logging.getLogger(__name__)


def assign_gene_families(
    cluster_results: Dict[FeatureType, Dict[str, List[str]]],
    id_mappings: Dict[str, str],
    output_dir: Path,
    family_prefix: str = ""
) -> Tuple[Dict[str, str], Dict[str, List[Dict[str, Any]]]]:
    """
    Create gene families from cluster assignments.

    Args:
        cluster_results: Clustering results by feature type
        id_mappings: Mapping from PanGenomePlus IDs to original IDs
        output_dir: Directory for family assignment outputs
        family_prefix: Prefix for family IDs

    Returns:
        Tuple of (gene_to_family_mapping, family_member_details)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_to_family = {}
    family_members = {}
    family_counter = 1

    # Process each feature type
    for feature_type, clusters in cluster_results.items():
        logger.info(f"Assigning families for {len(clusters)} clusters of type {feature_type.value}")

        type_family_counter = 1

        for representative, members in clusters.items():
            # Create family ID
            if family_prefix:
                family_id = f"{family_prefix}_{family_counter:06d}"
            else:
                family_id = f"{family_counter:06d}"

            # Record family members
            family_members[family_id] = []

            for member_id in members:
                # Map to family
                gene_to_family[member_id] = family_id

                # Store member details
                member_info = {
                    "pangenomeplus_id": member_id,
                    "original_id": id_mappings.get(member_id, member_id),
                    "feature_type": feature_type.value,
                    "is_representative": member_id == representative,
                    "family_id": family_id
                }
                family_members[family_id].append(member_info)

            family_counter += 1
            type_family_counter += 1

        logger.info(f"Created {type_family_counter - 1} families for {feature_type.value}")

    # Handle singletons (genes not in any cluster)
    all_clustered_genes = set(gene_to_family.keys())
    all_genes = set(id_mappings.keys())
    singletons = all_genes - all_clustered_genes

    if singletons:
        logger.info(f"Creating families for {len(singletons)} singleton genes")

        for singleton_id in singletons:
            family_id = f"SINGLETON_{singleton_id}"
            gene_to_family[singleton_id] = family_id

            family_members[family_id] = [{
                "pangenomeplus_id": singleton_id,
                "original_id": id_mappings.get(singleton_id, singleton_id),
                "feature_type": "unknown",  # Would need to lookup from features
                "is_representative": True,
                "family_id": family_id
            }]

    # Save family assignments
    save_family_assignments(gene_to_family, family_members, output_dir)

    logger.info(f"Created {len(family_members)} total gene families")
    logger.info(f"Assigned {len(gene_to_family)} genes to families")

    return gene_to_family, family_members



def save_family_assignments(
    gene_to_family: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    output_dir: Path
) -> None:
    """
    Save gene family assignments to files.

    Args:
        gene_to_family: Gene to family mapping
        family_members: Family member details
        output_dir: Output directory
    """
    # Save gene to family mapping
    mapping_file = output_dir / "gene_to_family.tsv"
    with open(mapping_file, 'w') as f:
        f.write("gene_id\tfamily_id\n")
        for gene_id, family_id in gene_to_family.items():
            f.write(f"{gene_id}\t{family_id}\n")

    # Save family member details
    families_file = output_dir / "family_members.json"
    with open(families_file, 'w') as f:
        json.dump(family_members, f, indent=2)

    # Save family summary
    summary_file = output_dir / "family_summary.tsv"
    with open(summary_file, 'w') as f:
        f.write("family_id\tsize\tfeature_types\trepresentative\n")

        for family_id, members in family_members.items():
            size = len(members)
            feature_types = set(m["feature_type"] for m in members)
            representative = next(
                (m["pangenomeplus_id"] for m in members if m["is_representative"]), ""
            )

            f.write(f"{family_id}\t{size}\t{','.join(feature_types)}\t{representative}\n")

    logger.info(f"Saved family assignments to {output_dir}")