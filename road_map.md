# Feature Roadmap: Sequence Fetching and Design Region Mining

## Overview
This roadmap outlines the implementation of a Python-based system to:
1. **Fetch target and off-target sequences** from public databases
2. **Mine alignments for candidate "good design regions"** for primer/probe design

The system will front-load region picking to make downstream primer/probe design faster and more specific.

## Phase 1: Database Integration and Sequence Fetching

### 1.1 NCBI Integration Module
**Timeline: 2-3 weeks**

**Dependencies to add:**
```python
# Add to pyproject.toml dependencies:
# Note: entrez-direct is a command-line tool, not a Python package
# NCBI access is handled through BioPython's Entrez module (already available)
"requests>=2.25.0",  # Already present
"pandas>=1.5.0",     # Already present, ensuring version compatibility
```

**Features:**
- **NCBI Nucleotide Database Access**
  - Query by taxon, gene/marker, molecule type, length
  - Targeted pulls for COI, 12S/16S, ITS markers
  - Bulk download complete genomes/mitogenomes by taxonomy ID
  
- **Implementation Structure:**
  ```
  src/genetic_testing/sequence_fetching/
  ├── __init__.py
  ├── ncbi_client.py          # NCBI API wrapper
  ├── sequence_filters.py     # Length, quality filters
  └── database_models.py      # Pydantic models for sequences
  ```

**Key Functions:**
- `fetch_target_sequences(taxon, marker, max_results)`
- `fetch_off_target_sequences(genus, exclude_species, marker)`
- `download_mitogenomes(taxonomy_ids)`

### 1.2 Specialized Database Connectors
**Timeline: 3-4 weeks**

**Dependencies to add:**
```python
"requests-cache>=1.0.0",  # Rate limiting and caching
"aiohttp>=3.8.0",         # Async HTTP requests
```

**BOLD Database Integration:**
- COI barcode sequences via BOLD API
- Rate-limited chunked downloads
- Species and specimen metadata

**SILVA Database Integration:**
- 16S/18S/23S/28S curated rRNA sets
- Broad bacterial/archaeal/eukaryotic coverage
- FASTA download and parsing

**UNITE Database Integration:**
- Fungal ITS/ITS1/ITS2 reference sets
- Species-level fungal identification

**Implementation:**
```
src/genetic_testing/sequence_fetching/
├── bold_client.py          # BOLD API integration
├── silva_client.py         # SILVA database access
├── unite_client.py         # UNITE database access
└── async_downloader.py     # Async batch downloading
```

## Phase 2: Sequence Processing and Quality Control

### 2.1 Sequence Cleaning Pipeline
**Timeline: 2 weeks**

**Leveraging existing dependencies:**
- `biopython>=1.81` (already present)
- `cialign` (already present)

**Features:**
- **Quality Filtering:**
  - Remove short/partial sequences
  - Drop sequences with excessive N's
  - Length-based filtering

- **Deduplication:**
  - Remove highly similar sequences per species
  - Prevent alignment bias

- **Masking:**
  - Low-complexity region masking
  - Poly-run detection and masking
  - Repeat sequence identification

**Implementation:**
```
src/genetic_testing/sequence_analysis/
├── quality_control.py      # Sequence QC pipeline
├── deduplication.py        # Similarity-based deduplication
└── masking.py             # Low-complexity masking
```

### 2.2 Chimera Detection (Optional)
**Timeline: 1 week**

**New dependency:**
```python
"vsearch>=2.22.1",  # Chimera detection
```

- Environmental/amplicon data screening
- Integration with existing sequence analysis tools

## Phase 3: Multiple Sequence Alignment and Analysis

### 3.1 Alignment Engine
**Timeline: 2-3 weeks**

**Dependencies to add:**
```python
"mafft>=7.0",              # Multiple sequence alignment
"muscle>=5.0",             # Alternative aligner
"clustalo>=1.2.4",         # Another alignment option
```

**Features:**
- **Multi-algorithm alignment support:**
  - MAFFT (fast, accurate)
  - MUSCLE (traditional choice)
  - Clustal Omega (large datasets)

- **Alignment processing:**
  - Target vs off-target grouping
  - Alignment quality assessment
  - Gap handling and trimming

**Implementation:**
```
src/genetic_testing/sequence_analysis/
├── alignment_engine.py     # Multi-algorithm alignment
├── alignment_processor.py  # Post-alignment processing
└── group_manager.py       # Target/off-target grouping
```

### 3.2 Phylogenetic Context (Enhanced)
**Timeline: 2 weeks**

**New dependencies:**
```python
"ete3>=3.1.2",            # Phylogenetic analysis
"dendropy>=4.5.0",        # Tree manipulation
```

- Phylogenetic tree construction
- Evolutionary distance calculations
- Clade-aware design region selection

## Phase 4: Design Region Discovery Engine

### 4.1 Signature Region Identification
**Timeline: 3-4 weeks**

**Core Algorithm Implementation:**
- **Group-specific signature windows:**
  - Conserved within target group
  - Divergent from off-targets
  - Configurable window sizes (80-200 bp)

- **Scoring metrics:**
  - Target conservation score
  - Off-target divergence score
  - Combined specificity index

**Implementation:**
```
src/genetic_testing/design_regions/
├── __init__.py
├── signature_finder.py     # Core signature detection
├── conservation_scorer.py  # Conservation metrics
├── specificity_analyzer.py # Target vs off-target analysis
└── region_ranker.py       # Candidate region ranking
```

### 4.2 Advanced Region Analysis
**Timeline: 2-3 weeks**

**Features:**
- **Physicochemical properties:**
  - GC content analysis (target: 40-60%)
  - Melting temperature estimation
  - Secondary structure prediction

- **Complexity filtering:**
  - Homopolymer run detection
  - Low-complexity region exclusion
  - Repeat sequence avoidance

**New dependencies:**
```python
"primer3-py>=2.0.0",      # Primer design constraints
"ViennaRNA>=2.5.0",       # RNA secondary structure
```

### 4.3 In-silico Validation
**Timeline: 2 weeks**

**Features:**
- **BLAST integration:**
  - Local BLAST database creation
  - Specificity checking against nt database
  - Custom database BLAST (SILVA/UNITE/BOLD)

- **Coverage analysis:**
  - Target haplotype coverage assessment
  - Strain/variant compatibility

**Dependencies:**
```python
"blast>=2.13.0",          # Local BLAST
"ncbi-blast+>=2.13.0",    # BLAST+ suite
```

## Phase 5: API and Integration Layer

### 5.1 REST API Development
**Timeline: 2 weeks**

**Leveraging existing infrastructure:**
- `streamlit>=1.26.0` (already present for web UI)
- FastAPI integration for REST endpoints

**Endpoints:**
- `/api/v1/fetch-sequences` - Database querying
- `/api/v1/analyze-regions` - Design region analysis
- `/api/v1/validate-regions` - In-silico validation

### 5.2 Web Interface Enhancement
**Timeline: 2-3 weeks**

**Streamlit Dashboard:**
- Sequence fetching interface
- Alignment visualization
- Design region exploration
- Results export functionality

## Phase 6: Output and Export System

### 6.1 Results Management
**Timeline: 1-2 weeks**

**Features:**
- **Ranked candidate regions:**
  - Conservation scores
  - Specificity metrics
  - Physicochemical properties

- **Export formats:**
  - FASTA sequences for downstream tools
  - CSV/TSV for analysis
  - JSON for API integration

**Implementation:**
```
src/genetic_testing/output/
├── result_formatter.py    # Multiple output formats
├── fasta_exporter.py     # BLAST-ready sequences
└── report_generator.py   # Summary reports
```

## Implementation Timeline

| Phase | Duration | Dependencies |
|-------|----------|--------------|
| Phase 1: Database Integration | 5-7 weeks | NCBI APIs, BOLD, SILVA, UNITE |
| Phase 2: Sequence Processing | 3 weeks | BioPython enhancements |
| Phase 3: Alignment & Analysis | 4-5 weeks | MAFFT, phylogenetics |
| Phase 4: Design Region Engine | 7-9 weeks | Core algorithm development |
| Phase 5: API & Integration | 4-5 weeks | FastAPI, Streamlit |
| Phase 6: Output System | 1-2 weeks | Export functionality |

**Total Estimated Timeline: 24-31 weeks**

## Technical Architecture

```
neglected-diagnostics/
├── src/genetic_testing/
│   ├── sequence_fetching/     # Phase 1
│   │   ├── ncbi_client.py
│   │   ├── bold_client.py
│   │   ├── silva_client.py
│   │   └── unite_client.py
│   ├── sequence_analysis/     # Phases 2-3
│   │   ├── quality_control.py
│   │   ├── alignment_engine.py
│   │   └── phylogenetic_analysis.py
│   ├── design_regions/        # Phase 4
│   │   ├── signature_finder.py
│   │   ├── conservation_scorer.py
│   │   └── region_ranker.py
│   ├── validation/           # Phase 4
│   │   ├── blast_validator.py
│   │   └── coverage_analyzer.py
│   ├── api/                  # Phase 5
│   │   ├── endpoints.py
│   │   └── models.py
│   └── output/               # Phase 6
│       ├── result_formatter.py
│       └── report_generator.py
└── tests/                    # Continuous
    ├── test_fetching/
    ├── test_analysis/
    └── test_regions/
```

## Success Metrics

1. **Database Coverage:**
   - Support for 4+ major databases (NCBI, BOLD, SILVA, UNITE)
   - Successful sequence retrieval for major taxonomic groups

2. **Analysis Quality:**
   - Alignment accuracy >95% for test datasets
   - Design region specificity >90% in validation

3. **Performance:**
   - Process 1000+ sequences in <10 minutes
   - API response time <30 seconds for typical queries

4. **Usability:**
   - Web interface for non-technical users
   - Export compatibility with major primer design tools

This roadmap provides a comprehensive Python-based alternative to the R-centric workflow described in the reference, leveraging the existing project infrastructure and adding specialized bioinformatics capabilities for sequence-based diagnostic design.