"""Configuration management and validation."""

import yaml
import json
from pathlib import Path
from typing import Dict, Any, Union, Optional

from pangenomeplus.core.types import ValidationResult
from pangenomeplus.core.exceptions import ConfigurationError


# JSON Schema for configuration validation
CONFIG_SCHEMA = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "type": "object",
    "required": ["clustering", "output"],
    "properties": {
        "clustering": {
            "type": "object",
            "properties": {
                "protein": {
                    "type": "object",
                    "properties": {
                        "identity": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                        "coverage": {"type": "number", "minimum": 0.0, "maximum": 1.0}
                    }
                },
                "trna_identity": {
                    "type": "number",
                    "minimum": 0.0,
                    "maximum": 1.0,
                    "default": 0.9
                },
                "rrna_identity": {
                    "type": "number",
                    "minimum": 0.0,
                    "maximum": 1.0,
                    "default": 0.97
                },
                "use_gpu": {
                    "type": "boolean",
                    "default": False
                },
                "memory_limit": {
                    "type": ["string", "null"],
                    "pattern": "^[0-9]+[GMK]?$",
                    "default": None
                }
            }
        },
        "annotation": {
            "type": "object",
            "properties": {
                "tools": {
                    "type": "array",
                    "items": {
                        "type": "string",
                        "enum": ["prodigal", "trnascan", "barrnap", "minced"]
                    },
                    "default": ["prodigal", "trnascan", "barrnap"]
                },
                "genetic_code": {
                    "type": "integer",
                    "minimum": 1,
                    "maximum": 31,
                    "default": 11
                }
            }
        },
        "paralogs": {
            "type": "object",
            "properties": {
                "handling": {
                    "type": "string",
                    "enum": ["keep", "split", "both"],
                    "default": "keep"
                },
                "synteny_window": {
                    "type": "integer",
                    "minimum": 1,
                    "maximum": 20,
                    "default": 5
                }
            }
        },
        "output": {
            "type": "object",
            "properties": {
                "formats": {
                    "type": "array",
                    "items": {
                        "type": "string",
                        "enum": ["transformer", "roary", "presence_absence", "fasta"]
                    },
                    "default": ["transformer", "presence_absence"]
                }
            }
        },
        "hpc": {
            "type": "object",
            "properties": {
                "use_slurm": {
                    "type": "boolean",
                    "default": False
                },
                "partition": {
                    "type": "string",
                    "default": "compute"
                },
                "max_runtime": {
                    "type": "string",
                    "pattern": "^[0-9]+:[0-9]{2}:[0-9]{2}$",
                    "default": "24:00:00"
                }
            }
        }
    }
}


def create_default_configuration(
    n_genomes: Optional[int] = None
) -> Dict[str, Any]:
    """
    Create default configuration with optional optimizations.

    Args:
        n_genomes: Number of genomes for parameter optimization

    Returns:
        Default configuration dictionary
    """
    config = {
        "pipeline": {
            "name": "pangenome_analysis",
            "version": "1.0.0",
            "description": "Adaptive scale pangenome analysis"
        },
        "resources": {
            "threads": 4,
            "memory_limit": None,
            "temp_dir": None,
            "keep_temp": False
        },
        "annotation": {
            "tools": ["prodigal", "trnascan", "barrnap"],
            "prodigal": {
                "genetic_code": 11,
                "procedure": "single",
                "closed_ends": True
            },
            "trnascan": {
                "search_mode": "bacterial",
                "relaxed": False
            },
            "barrnap": {
                "kingdom": "bac",
                "evalue": 1e-6,
                "length_cutoff": 0.8
            },
            "intergenic": {
                "extract": True,
                "min_length": 50
            }
        },
        "clustering": {
            "protein": {
                "identity": 0.8,
                "coverage": 0.8,
                "coverage_mode": 0,
                "cluster_mode": 0,
                "sensitivity": 7.5
            },
            "trna": {
                "identity": 0.9,
                "coverage": 0.9,
                "coverage_mode": 0,
                "sensitivity": 7.5
            },
            "rrna": {
                "identity": 0.97,
                "coverage": 0.9,
                "coverage_mode": 0,
                "sensitivity": 7.5
            },
            "intergenic": {
                "identity": 0.8,
                "coverage": 0.8,
                "coverage_mode": 0,
                "sensitivity": 7.5
            },
            "gpu": {
                "use_gpu": False,
                "devices": [],
                "memory_fraction": 0.8
            }
        },
        "paralogs": {
            "handling": "keep",
            "synteny_window": 5,
            "min_synteny_conservation": 0.8
        },
        "classification": {
            "hard_core_threshold": 1.0,
            "soft_core_threshold": 0.95,
            "shell_threshold": 0.15,
            "cloud_threshold": 0.15
        },
        "output": {
            "formats": ["transformer", "presence_absence"],
            "transformer": {
                "include_singletons": True,
                "singleton_prefix": "S",
                "compress": False,
                "feature_type_prefix": True
            },
            "presence_absence": {
                "format": "tsv",
                "binary": True
            }
        },
        "checkpoints": {
            "enabled": True,
            "frequency": "stage",
            "keep_n": 5,
            "compression": True
        },
        "logging": {
            "level": "INFO",
            "file": None,
            "format": "detailed"
        },
        "validation": {
            "input_files": True,
            "system_requirements": True,
            "tool_versions": True,
            "intermediate_results": True
        }
    }

    # Optimize for specific strategies if provided
    # Simple optimization based on dataset size
    if n_genomes and n_genomes > 500:
        config["clustering"]["protein"]["sensitivity"] = 6.0
        config["resources"]["threads"] = min(8, 4)

    return config


def load_configuration(config_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load and validate configuration from file.

    Args:
        config_path: Path to configuration file (YAML or JSON)

    Returns:
        Validated configuration dictionary

    Raises:
        ConfigurationError: Invalid configuration
        FileNotFoundError: Configuration file not found
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    try:
        with open(config_path, 'r') as f:
            if config_path.suffix.lower() in ['.yaml', '.yml']:
                config = yaml.safe_load(f)
            elif config_path.suffix.lower() == '.json':
                config = json.load(f)
            else:
                raise ConfigurationError(f"Unsupported config file format: {config_path.suffix}")

    except (yaml.YAMLError, json.JSONDecodeError) as e:
        raise ConfigurationError(f"Error parsing configuration file: {e}")

    # Validate against schema
    validation_result = validate_configuration_schema(config)
    if not validation_result.is_valid:
        raise ConfigurationError(f"Configuration validation failed: {validation_result.errors}")

    return config


def validate_configuration_schema(config: Dict[str, Any]) -> ValidationResult:
    """
    Simple configuration validation without external dependencies.

    Args:
        config: Configuration dictionary

    Returns:
        ValidationResult with validation status
    """
    errors = []
    warnings = []

    try:
        # Basic validation checks
        if "clustering" not in config:
            warnings.append("Missing 'clustering' configuration - using defaults")

        if "output" not in config:
            warnings.append("Missing 'output' configuration - using defaults")

        # Validate clustering config structure
        if "clustering" in config:
            clustering = config["clustering"]
            if not isinstance(clustering, dict):
                errors.append("'clustering' must be a dictionary")

        # Validate output config structure
        if "output" in config:
            output = config["output"]
            if not isinstance(output, dict):
                errors.append("'output' must be a dictionary")

        is_valid = len(errors) == 0

    except Exception as e:
        errors.append(f"Configuration validation error: {e}")
        is_valid = False

    return ValidationResult(
        is_valid=is_valid,
        errors=errors,
        warnings=warnings,
        details={"validation": "basic_checks_completed"}
    )


def merge_configurations(
    base_config: Dict[str, Any],
    override_config: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Merge configuration dictionaries with validation.

    Args:
        base_config: Base configuration
        override_config: Override parameters

    Returns:
        Merged configuration
    """
    def merge_dicts(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
        """Recursively merge dictionaries."""
        result = base.copy()
        for key, value in override.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = merge_dicts(result[key], value)
            else:
                result[key] = value
        return result

    merged = merge_dicts(base_config, override_config)

    # Validate merged configuration
    validation_result = validate_configuration_schema(merged)
    if not validation_result.is_valid:
        raise ConfigurationError(
            f"Merged configuration validation failed: {validation_result.errors}"
        )

    return merged


def save_configuration(config: Dict[str, Any], output_path: Union[str, Path]) -> None:
    """
    Save configuration to file.

    Args:
        config: Configuration dictionary
        output_path: Output file path

    Raises:
        ConfigurationError: Error saving configuration
    """
    output_path = Path(output_path)

    try:
        with open(output_path, 'w') as f:
            if output_path.suffix.lower() in ['.yaml', '.yml']:
                yaml.dump(config, f, default_flow_style=False, indent=2)
            elif output_path.suffix.lower() == '.json':
                json.dump(config, f, indent=2)
            else:
                raise ConfigurationError(f"Unsupported output format: {output_path.suffix}")

    except (yaml.YAMLError, json.JSONEncodeError, IOError) as e:
        raise ConfigurationError(f"Error saving configuration: {e}")