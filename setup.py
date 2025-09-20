"""Setup script for PanGenomePlus package."""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "README.md"
if readme_file.exists():
    with open(readme_file, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "Adaptive scale pangenome analysis pipeline"

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    with open(requirements_file, "r", encoding="utf-8") as f:
        requirements = [
            line.strip() for line in f
            if line.strip() and not line.startswith("#")
        ]
else:
    requirements = [
        "click>=8.0.0",
        "pyyaml>=6.0",
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "psutil>=5.8.0",
        "scipy>=1.7.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "biopython>=1.79"
    ]

setup(
    name="pangenomeplus",
    version="1.0.0",
    author="PanGenomePlus Team",
    description="Adaptive scale pangenome analysis pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "pangenomeplus=pangenomeplus.main:cli",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="pangenome genomics bioinformatics clustering mmseqs2",
)