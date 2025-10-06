from setuptools import find_packages, setup

# Import version from single source of truth
from pangenomeplus import __version__

setup(
    name="pangenomeplus",
    version=__version__,
    description="Comprehensive pangenome analysis pipeline for all genomic features",
    author="James McInerney",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.79",
        "matplotlib>=3.5.0",
        "seaborn>=0.12.0",
        "pytest>=7.0",
        "pytest-cov>=4.0",
        "black>=22.0",
        "isort>=5.0",
        "flake8>=5.0",
        "mypy>=0.991",
        "pre-commit>=2.20",
        "rich>=13.0.0",
    ],
    entry_points={
        "console_scripts": [
            "pangenomeplus=pangenomeplus.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
