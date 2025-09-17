# Phase 1 & 2 Implementation Actions

## Phase 1: Database Integration MCP Server
**Timeline: 5-7 weeks**  
**Container**: `ndiag-database-server:latest`

### Week 1: Project Setup & Infrastructure

#### 1.1 Repository Structure Setup
- [ ] Create `mcp_servers/database_server/` directory
- [ ] Initialize `requirements.txt` with core dependencies:
  ```
  gget>=0.28.0
  biopython>=1.81
  requests>=2.31.0
  pandas>=2.0.0
  pysradb>=1.4.0
  google-cloud-bigquery>=3.11.0
  boto3>=1.28.0
  mcp>=1.0.0
  ```
- [ ] Create base `Dockerfile` structure
- [ ] Set up development environment with Docker Compose
- [ ] Initialize Git submodule for MCP server development

#### 1.2 MCP Server Foundation
- [ ] Install MCP Python SDK and dependencies
- [ ] Create base `database_mcp_server.py` with MCP server initialization
- [ ] Implement MCP server configuration and health check endpoints
- [ ] Set up logging and error handling framework
- [ ] Create unit test structure with pytest

### Week 2: Core Database Integrations

#### 2.1 gget Integration Implementation
- [ ] Implement `gget_search` tool:
  - Input validation for searchwords, species, id_type, andor
  - Error handling for API failures
  - Response formatting and standardization
- [ ] Implement `gget_ref` tool:
  - Species name validation and normalization
  - Support for different Ensembl releases
  - File format handling (GTF, FASTA)
- [ ] Implement `gget_info` tool:
  - Ensembl ID validation
  - Expand parameter handling
  - Rich metadata extraction
- [ ] Implement `gget_seq` tool:
  - Sequence type selection (genomic, transcript, protein)
  - Translation parameter handling
  - FASTA output formatting

#### 2.2 NCBI Direct Integration
- [ ] Install and configure EDirect tools in container
- [ ] Implement NCBI taxonomy lookup functions
- [ ] Create sequence retrieval functions for GenBank/RefSeq
- [ ] Add support for bulk sequence downloads
- [ ] Implement rate limiting and retry logic for NCBI API

### Week 3: Specialized Database Sources

#### 3.1 BOLD Database Integration
- [ ] Implement BOLD API client
- [ ] Create `get_sequences` tool with BOLD source support:
  - Taxon name resolution
  - COI sequence retrieval
  - Metadata extraction (collection data, geography)
- [ ] Add BOLD-specific error handling
- [ ] Implement result caching for repeated queries

#### 3.2 SILVA Database Integration
- [ ] Set up SILVA database access
- [ ] Implement 16S/18S rRNA sequence retrieval
- [ ] Add taxonomic classification support
- [ ] Create quality filtering for SILVA sequences

#### 3.3 UNITE Database Integration
- [ ] Implement UNITE API integration
- [ ] Add ITS sequence retrieval functionality
- [ ] Support for fungal taxonomic classifications
- [ ] Implement sequence quality metrics

### Week 4: SRA/BioProject Integration

#### 4.1 SRA Search Implementation
- [ ] Implement `search_sra_studies` tool:
  - Free text search with Entrez integration
  - Advanced filtering (organism, library_strategy, platform)
  - Date range filtering
  - Result pagination and limiting
- [ ] Add `get_sra_runinfo` tool:
  - Study accession validation
  - Sample metadata extraction
  - Multiple output formats (JSON, TSV, CSV)
- [ ] Implement robust error handling for SRA API timeouts

#### 4.2 Cloud SQL Integration
- [ ] Set up Google Cloud BigQuery client
- [ ] Implement `search_sra_cloud` tool:
  - SQL query validation and sanitization
  - BigQuery result processing
  - Row limiting and cost controls
- [ ] Add AWS Athena support as alternative
- [ ] Implement query result caching

### Week 5: Unified Tools & Advanced Features

#### 5.1 Unified `get_sequences` Tool
- [ ] Implement source routing logic (gget, ncbi, bold, silva, unite)
- [ ] Add intelligent source selection based on query type
- [ ] Implement result merging and deduplication
- [ ] Add format standardization across sources
- [ ] Create comprehensive input validation

#### 5.2 Taxonomic Tools Implementation
- [ ] Implement `get_neighbors` tool:
  - Taxonomic rank traversal
  - Distance-based neighbor finding
  - Common misidentification flagging
- [ ] Implement `get_taxonomy` tool:
  - Multi-source taxonomic resolution
  - Synonym handling
  - Taxonomic hierarchy extraction

### Week 6: Testing & Quality Assurance

#### 6.1 Unit Testing
- [ ] Write comprehensive unit tests for all MCP tools
- [ ] Mock external API calls for reliable testing
- [ ] Test error conditions and edge cases
- [ ] Implement test data fixtures for reproducible tests

#### 6.2 Integration Testing
- [ ] Test real API calls with rate limiting
- [ ] Validate output formats across all tools
- [ ] Test large-scale data retrieval scenarios
- [ ] Performance testing with concurrent requests

#### 6.3 Container Testing
- [ ] Build and test Docker container
- [ ] Validate all dependencies are properly installed
- [ ] Test container startup and health checks
- [ ] Memory and CPU usage profiling

### Week 7: Documentation & Deployment

#### 7.1 Documentation
- [ ] Write comprehensive API documentation for all tools
- [ ] Create usage examples for each MCP tool
- [ ] Document error codes and troubleshooting
- [ ] Create developer setup guide

#### 7.2 Deployment Preparation
- [ ] Optimize Docker image size
- [ ] Set up production logging configuration
- [ ] Configure environment variables for API keys
- [ ] Create Kubernetes deployment manifests
- [ ] Set up monitoring and alerting

---

## Phase 2: Sequence Processing MCP Server
**Timeline: 3 weeks**  
**Container**: `ndiag-processing-server:latest`

### Week 1: Foundation & Quality Control

#### 1.1 Project Setup
- [ ] Create `mcp_servers/processing_server/` directory
- [ ] Initialize `requirements.txt`:
  ```
  biopython>=1.81
  seqkit>=2.3.0
  vsearch>=2.22.0
  cialign>=1.0.0
  mcp>=1.0.0
  numpy>=1.24.0
  pandas>=2.0.0
  ```
- [ ] Create base `Dockerfile` with Ubuntu 22.04 and tool installations
- [ ] Set up `processing_mcp_server.py` with MCP framework

#### 1.2 FASTA Quality Control Implementation
- [ ] Implement `fasta_qc` tool:
  - Sequence length filtering (min_length parameter)
  - N-content percentage calculation and filtering
  - Invalid character detection and removal
  - Basic sequence statistics generation
- [ ] Add FASTA format validation and repair
- [ ] Implement sequence header standardization
- [ ] Create quality metrics reporting

#### 1.3 Duplicate Detection & Removal
- [ ] Implement `dereplicate_sequences` tool:
  - Exact sequence duplicate detection
  - Identity threshold-based clustering
  - Per-species dereplication option
  - Representative sequence selection logic
- [ ] Add sequence similarity calculation using vsearch
- [ ] Implement memory-efficient processing for large datasets

### Week 2: Advanced Processing Features

#### 2.1 Low Complexity Masking
- [ ] Implement `mask_low_complexity` tool:
  - Homopolymer run detection and masking
  - Repetitive sequence identification
  - Complexity scoring algorithm implementation
  - Soft vs hard masking options
- [ ] Integrate external masking tools (dustmasker, segmasker)
- [ ] Add custom complexity thresholds

#### 2.2 Chimera Detection
- [ ] Implement `detect_chimeras` tool:
  - Integration with vsearch chimera detection
  - Reference database selection (SILVA, UNITE)
  - Abundance-based chimera filtering
  - De novo chimera detection methods
- [ ] Add chimera removal and flagging options
- [ ] Implement batch processing for large datasets

#### 2.3 Pipeline Integration
- [ ] Implement `process_sequences` tool:
  - Configurable processing pipeline
  - Parameter passing to individual tools
  - Progress tracking and reporting
  - Error recovery and partial results
- [ ] Add pipeline validation and optimization
- [ ] Create processing reports and statistics

### Week 3: Testing, Optimization & Documentation

#### 3.1 Performance Optimization
- [ ] Implement parallel processing for CPU-intensive operations
- [ ] Add memory usage monitoring and optimization
- [ ] Create streaming processing for large files
- [ ] Optimize temporary file handling

#### 3.2 Comprehensive Testing
- [ ] Unit tests for all processing functions
- [ ] Integration tests with real sequence data
- [ ] Performance benchmarking with large datasets
- [ ] Memory leak testing and profiling
- [ ] Container integration testing

#### 3.3 Error Handling & Validation
- [ ] Robust input validation for all tools
- [ ] Comprehensive error messages and codes
- [ ] Graceful handling of malformed FASTA files
- [ ] Recovery mechanisms for partial failures

#### 3.4 Documentation & Deployment
- [ ] Complete API documentation for all tools
- [ ] Usage examples and best practices
- [ ] Performance tuning guidelines
- [ ] Troubleshooting guide
- [ ] Container deployment documentation

---

## Cross-Phase Integration Tasks

### Integration Testing
- [ ] Test Phase 1 → Phase 2 data flow
- [ ] Validate output format compatibility
- [ ] End-to-end workflow testing
- [ ] Performance testing of combined pipeline

### Documentation
- [ ] Create comprehensive workflow examples
- [ ] Document data formats and schemas
- [ ] API versioning and compatibility guide
- [ ] Troubleshooting and FAQ sections

### Deployment Infrastructure
- [ ] Docker Compose configuration for both phases
- [ ] Kubernetes manifests for production deployment
- [ ] Monitoring and logging setup
- [ ] Backup and recovery procedures

## Success Criteria

### Phase 1 Success Metrics
- [ ] Successfully retrieve sequences from all 5 database sources
- [ ] Handle concurrent requests (>10 simultaneous users)
- [ ] Process queries returning >1000 sequences
- [ ] 99.5% uptime during testing period
- [ ] Complete API documentation with examples

### Phase 2 Success Metrics
- [ ] Process >10,000 sequences in <5 minutes
- [ ] Memory usage <2GB for typical workloads
- [ ] All quality control metrics properly calculated
- [ ] Zero data loss during processing pipeline
- [ ] Comprehensive error reporting and recovery

### Combined Integration Metrics
- [ ] End-to-end workflow (database → processing) <10 minutes for 1000 sequences
- [ ] Seamless data format compatibility
- [ ] Proper error propagation between phases
- [ ] Complete provenance tracking through both phases
