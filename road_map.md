# Feature Roadmap: MCP-Backed Sequence Analysis & Design Region Mining

## Overview
This roadmap outlines the implementation of a **Model Context Protocol (MCP) based system** that enables agentic assistants to:
1. **Fetch target and off-target sequences** from public databases via MCP servers
2. **Mine alignments for candidate "good design regions"** through distributed MCP tools
3. **Provide seamless agent integration** with the existing assistant framework

The system transforms the monolithic approach into a **distributed, agent-friendly architecture** where each major capability is exposed as an MCP server that the assistant can orchestrate through natural language interactions.

## MCP Architecture Philosophy

### Core Principles
- **Agent-First Design**: Every capability is exposed as MCP tools that agents can discover and use
- **Distributed Services**: Each major function runs as an independent MCP server
- **Natural Language Interface**: Users interact through the assistant, not direct API calls
- **Containerized Isolation**: Each MCP server runs in isolated containers with resource limits
- **Secure by Default**: Remote MCP with ephemeral credentials and least-privilege access

## Phase 1: Database Integration MCP Server

### 1.1 Unified Database MCP Server
**Timeline: 5-7 weeks**
**Container**: `ndiag-database-server:latest` [[memory:6143671]]

**Consolidates**: NCBI, SRA/BioProject, BOLD, SILVA, UNITE databases + gget integration

**Unified MCP Tools for All Database Sources:**
```typescript
tools: [
  // Database Tools
  {
    name: "get_sequences",
    inputSchema: {
      type: "object",
      properties: {
        taxon: { type: "string", description: "Taxon name or ID" },
        region: { type: "string", enum: ["COI", "16S", "ITS", "mitogenome", "whole"] },
        source: { type: "string", enum: ["gget", "ncbi", "bold", "silva", "unite"], default: "gget" },
        max_results: { type: "integer", default: 100 },
        format: { type: "string", enum: ["fasta", "genbank"], default: "fasta" }
      },
      required: ["taxon"]
    }
  },
  {
    name: "gget_ref",
    inputSchema: {
      type: "object",
      properties: {
        species: { type: "string", description: "Species name (e.g., 'homo_sapiens')" },
        which: { type: "string", enum: ["all", "gtf", "fasta"], default: "all" },
        release: { type: "integer", description: "Ensembl release number" }
      },
      required: ["species"]
    }
  },
  {
    name: "gget_search",
    inputSchema: {
      type: "object",
      properties: {
        searchwords: { type: "array", items: { type: "string" } },
        species: { type: "string" },
        id_type: { type: "string", enum: ["gene", "transcript"], default: "gene" },
        andor: { type: "string", enum: ["and", "or"], default: "or" }
      },
      required: ["searchwords", "species"]
    }
  },
  {
    name: "gget_info",
    inputSchema: {
      type: "object",
      properties: {
        ens_ids: { type: "array", items: { type: "string" } },
        expand: { type: "boolean", default: false }
      },
      required: ["ens_ids"]
    }
  },
  {
    name: "gget_seq",
    inputSchema: {
      type: "object",
      properties: {
        ens_ids: { type: "array", items: { type: "string" } },
        translate: { type: "boolean", default: false },
        seqtype: { type: "string", enum: ["genomic", "transcript", "protein"], default: "transcript" }
      },
      required: ["ens_ids"]
    }
  },
  {
    name: "get_neighbors",
    inputSchema: {
      type: "object", 
      properties: {
        taxon: { type: "string" },
        rank: { type: "string", enum: ["species", "genus", "family"] },
        distance: { type: "integer", default: 1 },
        common_misIDs: { type: "boolean", default: false }
      },
      required: ["taxon", "rank"]
    }
  },
  {
    name: "get_taxonomy",
    inputSchema: {
      type: "object",
      properties: {
        query: { type: "string", description: "Taxon name or accession" }
      },
      required: ["query"]
    }
  },
  // SRA/BioProject Tools
  {
    name: "search_sra_studies",
    inputSchema: {
      type: "object",
      properties: {
        query: { type: "string", description: "Free text search query" },
        filters: { type: "object", properties: {
          organism: { type: "string" },
          library_strategy: { type: "string", enum: ["AMPLICON", "RNA-Seq", "WGS", "METAGENOMIC"] },
          platform: { type: "string", enum: ["ILLUMINA", "PACBIO", "OXFORD_NANOPORE"] },
          submission_date_start: { type: "string", format: "date" },
          submission_date_end: { type: "string", format: "date" },
          min_read_length: { type: "integer" },
          max_results: { type: "integer", default: 100 }
        }},
        search_method: { type: "string", enum: ["entrez", "cloud_sql"], default: "entrez" }
      },
      required: ["query"]
    }
  },
  {
    name: "get_sra_runinfo",
    inputSchema: {
      type: "object",
      properties: {
        study_accession: { type: "string", description: "BioProject/SRA Study ID" },
        include_sample_metadata: { type: "boolean", default: true },
        format: { type: "string", enum: ["json", "tsv", "csv"], default: "json" }
      },
      required: ["study_accession"]
    }
  },
  {
    name: "search_sra_cloud",
    inputSchema: {
      type: "object",
      properties: {
        query_sql: { type: "string", description: "SQL query for BigQuery/Athena SRA tables" },
        platform: { type: "string", enum: ["bigquery", "athena"], default: "bigquery" },
        max_rows: { type: "integer", default: 1000 }
      },
      required: ["query_sql"]
    }
  }
]
```

**Consolidated Container:**
```dockerfile
FROM python:3.11-slim
RUN apt-get update && apt-get install -y curl wget
# Install EDirect tools for NCBI
RUN curl -s https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh | bash
# Install gget and all database dependencies
RUN pip install gget biopython requests pandas pysradb google-cloud-bigquery boto3
COPY database_mcp_server.py /app/
EXPOSE 8000
CMD ["python", "/app/database_mcp_server.py"]
```

**Agent Integration:**
- **User**: "Find COI sequences for salmon and get reference genome info"
- **Assistant**: 
  1. `gget_search(searchwords=["COI", "cytochrome oxidase"], species="salmo_salar")` for gene IDs
  2. `gget_seq(ens_ids=[...], seqtype="transcript")` for sequences
  3. `gget_ref(species="salmo_salar", which="fasta")` for reference genome
  4. `search_sra_studies(query="Salmo salar COI")` for recent studies
- **Result**: Comprehensive sequence collection with standardized gene annotations

## Phase 2: Sequence Processing MCP Server

### 2.1 Unified Processing MCP Server
**Timeline: 3 weeks**
**Container**: `ndiag-processing-server:latest`

**Consolidates**: Quality control, deduplication, masking, chimera detection

**MCP Tools for Complete Sequence Processing:**
```typescript
tools: [{
  name: "fasta_qc",
  inputSchema: {
    type: "object",
    properties: {
      fasta_content: { type: "string" },
      min_length: { type: "integer", default: 100 },
      max_n_percent: { type: "number", default: 5.0 },
      remove_duplicates: { type: "boolean", default: true }
    },
    required: ["fasta_content"]
  }
}, {
  name: "dereplicate_sequences", 
  inputSchema: {
    type: "object",
    properties: {
      fasta_content: { type: "string" },
      identity_threshold: { type: "number", default: 0.97 },
      per_species: { type: "boolean", default: true }
    },
    required: ["fasta_content"]
  }
}, {
  name: "mask_low_complexity",
  inputSchema: {
    type: "object",
    properties: {
      fasta_content: { type: "string" },
      mask_repeats: { type: "boolean", default: true },
      mask_homopolymers: { type: "boolean", default: true },
      min_complexity: { type: "number", default: 1.5 }
    },
    required: ["fasta_content"]
  }
}, {
  name: "detect_chimeras",
  inputSchema: {
    type: "object",
    properties: {
      fasta_content: { type: "string" },
      reference_db: { type: "string", enum: ["auto", "silva", "unite"], default: "auto" },
      abundance_threshold: { type: "number", default: 2.0 }
    },
    required: ["fasta_content"]
  }
}, {
  name: "process_sequences",
  inputSchema: {
    type: "object",
    properties: {
      fasta_content: { type: "string" },
      pipeline: { type: "array", items: { 
        type: "string", 
        enum: ["qc", "dereplicate", "mask", "chimera"] 
      }, default: ["qc", "dereplicate"] },
      qc_params: { type: "object" },
      derep_params: { type: "object" },
      mask_params: { type: "object" },
      chimera_params: { type: "object" }
    },
    required: ["fasta_content"]
  }
}]
```

**Consolidated Container:**
```dockerfile
FROM ubuntu:22.04
RUN apt-get update && apt-get install -y seqkit vsearch python3-pip
RUN pip3 install biopython cialign
COPY processing_mcp_server.py /app/
EXPOSE 8001
CMD ["python3", "/app/processing_mcp_server.py"]
```

**Agent Integration:**
- **User**: "Clean up these sequences - remove duplicates and check for chimeras"
- **Assistant**: `process_sequences(pipeline=["qc", "dereplicate", "chimera"])`
- **Result**: Complete processed sequence set ready for alignment

## Phase 3: Alignment & Phylogenetics MCP Server

### 3.1 Unified Alignment & Analysis MCP Server
**Timeline: 4-5 weeks**
**Container**: `ndiag-alignment-server:latest`

**Consolidates**: MAFFT, MUSCLE, Clustal Omega alignment + phylogenetic analysis + gget integration

**MCP Tools for Complete Alignment & Analysis:**
```typescript
tools: [{
  name: "align_sequences",
  inputSchema: {
    type: "object",
    properties: {
      fasta_content: { type: "string" },
      algorithm: { type: "string", enum: ["gget_muscle", "mafft", "muscle", "clustalo"], default: "gget_muscle" },
      mafft_strategy: { type: "string", enum: ["auto", "linsi", "ginsi", "einsi"], default: "auto" },
      max_iterations: { type: "integer", default: 1000 },
      target_group: { type: "array", items: { type: "string" } }
    },
    required: ["fasta_content"]
  }
}, {
  name: "gget_muscle",
  inputSchema: {
    type: "object",
    properties: {
      sequences: { type: "array", items: { type: "string" } },
      super5: { type: "boolean", default: false },
      out: { type: "string", description: "Output file path" }
    },
    required: ["sequences"]
  }
}, {
  name: "process_alignment",
  inputSchema: {
    type: "object",
    properties: {
      alignment_content: { type: "string" },
      trim_gaps: { type: "boolean", default: true },
      gap_threshold: { type: "number", default: 0.5 },
      assess_quality: { type: "boolean", default: true }
    },
    required: ["alignment_content"]
  }
}, {
  name: "build_phylogeny",
  inputSchema: {
    type: "object",
    properties: {
      alignment_content: { type: "string" },
      method: { type: "string", enum: ["nj", "ml", "mp"], default: "nj" },
      bootstrap: { type: "integer", default: 100 },
      model: { type: "string", enum: ["p-distance", "jukes-cantor", "kimura"], default: "kimura" }
    },
    required: ["alignment_content"]
  }
}, {
  name: "calculate_distances",
  inputSchema: {
    type: "object",
    properties: {
      alignment_content: { type: "string" },
      model: { type: "string", enum: ["p-distance", "jukes-cantor", "kimura"], default: "kimura" }
    },
    required: ["alignment_content"]
  }
}, {
  name: "align_and_analyze",
  inputSchema: {
    type: "object",
    properties: {
      fasta_content: { type: "string" },
      algorithm: { type: "string", enum: ["mafft", "muscle", "clustalo"], default: "mafft" },
      include_phylogeny: { type: "boolean", default: false },
      include_distances: { type: "boolean", default: false },
      target_group: { type: "array", items: { type: "string" } }
    },
    required: ["fasta_content"]
  }
}]
```

**Consolidated Container:**
```dockerfile  
FROM ubuntu:22.04
RUN apt-get update && apt-get install -y mafft muscle clustalo python3-pip
RUN pip3 install gget biopython numpy ete3 dendropy
COPY alignment_mcp_server.py /app/
EXPOSE 8002  
CMD ["python3", "/app/alignment_mcp_server.py"]
```

**Agent Workflow:**
- **User**: "Align sequences and build a phylogenetic tree to understand relationships"
- **Assistant**: `align_and_analyze(include_phylogeny=true, include_distances=true)`
- **Result**: Complete alignment + phylogenetic analysis for clade-aware design

## Phase 4: Design & Primer Development MCP Server

### 4.1 Unified Design & Primer MCP Server
**Timeline: 7-9 weeks**
**Container**: `ndiag-design-server:latest`

**Consolidates**: Signature region discovery, conservation scoring, specificity analysis, Primer3 design, oligo QC

**MCP Tools for Complete Design Pipeline:**
```typescript
tools: [{
  name: "find_signature_regions",
  inputSchema: {
    type: "object",
    properties: {
      alignment_content: { type: "string" },
      target_sequences: { type: "array", items: { type: "string" } },
      window_size: { type: "integer", default: 150 },
      step_size: { type: "integer", default: 10 },
      min_conservation: { type: "number", default: 0.8 },
      min_divergence: { type: "number", default: 0.3 }
    },
    required: ["alignment_content", "target_sequences"]
  }
}, {
  name: "analyze_specificity", 
  inputSchema: {
    type: "object",
    properties: {
      candidate_regions: { type: "array", items: { type: "object" } },
      target_group: { type: "array", items: { type: "string" } },
      offtarget_group: { type: "array", items: { type: "string" } }
    },
    required: ["candidate_regions", "target_group", "offtarget_group"]
  }
}, {
  name: "rank_regions",
  inputSchema: {
    type: "object",
    properties: {
      scored_regions: { type: "array", items: { type: "object" } },
      weighting: { type: "object", properties: {
        conservation: { type: "number", default: 0.4 },
        specificity: { type: "number", default: 0.4 },
        complexity: { type: "number", default: 0.2 }
      }}
    },
    required: ["scored_regions"]
  }
}, {
  name: "primer3_design",
  inputSchema: {
    type: "object",
    properties: {
      template_fasta: { type: "string" },
      target_regions: { type: "array", items: { type: "array", items: { type: "number" } } },
      constraints: { type: "object", properties: {
        primer_size: { type: "array", items: { type: "number" }, default: [18, 22, 27] },
        tm: { type: "array", items: { type: "number" }, default: [57, 60, 63] },
        gc_content: { type: "array", items: { type: "number" }, default: [40, 50, 60] },
        product_size: { type: "array", items: { type: "number" }, default: [80, 150, 300] }
      }}
    },
    required: ["template_fasta"]
  }
}, {
  name: "oligo_qc",
  inputSchema: {
    type: "object",
    properties: {
      sequence: { type: "string" },
      salt_mM: { type: "number", default: 50 },
      mg_mM: { type: "number", default: 2 },
      oligo_conc_nM: { type: "number", default: 250 }
    },
    required: ["sequence"]
  }
}, {
  name: "design_primers_complete",
  inputSchema: {
    type: "object",
    properties: {
      alignment_content: { type: "string" },
      target_sequences: { type: "array", items: { type: "string" } },
      offtarget_sequences: { type: "array", items: { type: "string" } },
      primer_constraints: { type: "object" },
      region_params: { type: "object" }
    },
    required: ["alignment_content", "target_sequences"]
  }
}]
```

**Consolidated Container:**
```dockerfile
FROM ubuntu:22.04
RUN apt-get update && apt-get install -y primer3 python3-pip
RUN pip3 install primer3-py biopython ViennaRNA numpy scipy
COPY design_mcp_server.py /app/
EXPOSE 8003
CMD ["python3", "/app/design_mcp_server.py"]
```

**Agent Conversation Example:**
- **User**: "Design salmon-specific primers from this alignment"
- **Assistant**: `design_primers_complete(target_sequences=["Salmo_salar_*"], offtarget_sequences=["Oncorhynchus_*", "Thunnus_*"])`
- **Result**: Complete primer design pipeline from region discovery to validated oligos

## Phase 5: Validation & Literature MCP Server

### 5.1 Unified Validation & Research MCP Server
**Timeline: 4-5 weeks**
**Container**: `ndiag-validation-server:latest`

**Consolidates**: BLAST validation, in-silico PCR, coverage analysis, PubMed literature search + gget integration

**MCP Tools for Validation & Research:**
```typescript
tools: [{
  name: "gget_blast",
  inputSchema: {
    type: "object",
    properties: {
      sequence: { type: "string" },
      program: { type: "string", enum: ["blastn", "blastp", "blastx", "tblastn", "tblastx"], default: "blastn" },
      database: { type: "string", enum: ["nt", "nr", "refseq_rna", "refseq_protein"], default: "nt" },
      limit: { type: "integer", default: 50 },
      expect: { type: "number", default: 10.0 },
      low_comp_filt: { type: "boolean", default: false }
    },
    required: ["sequence"]
  }
}, {
  name: "gget_blat",
  inputSchema: {
    type: "object",
    properties: {
      sequence: { type: "string" },
      seqtype: { type: "string", enum: ["DNA", "protein", "translated%20RNA", "translated%20DNA"], default: "DNA" },
      assembly: { type: "string", default: "human" }
    },
    required: ["sequence"]
  }
}, {
  name: "blast_nt",
  inputSchema: {
    type: "object",
    properties: {
      query_fasta: { type: "string" },
      perc_identity: { type: "number", default: 90 },
      max_targets: { type: "integer", default: 50 },
      evalue: { type: "number", default: 0.001 }
    },
    required: ["query_fasta"]
  }
}, {
  name: "in_silico_pcr",
  inputSchema: {
    type: "object",
    properties: {
      forward_primer: { type: "string" },
      reverse_primer: { type: "string" },
      template_fasta: { type: "string" },
      max_mismatches: { type: "integer", default: 2 }
    },
    required: ["forward_primer", "reverse_primer", "template_fasta"]
  }
}, {
  name: "assess_coverage",
  inputSchema: {
    type: "object",
    properties: {
      primers: { type: "object", properties: {
        forward: { type: "string" },
        reverse: { type: "string" }
      }},
      target_panel: { type: "string" },
      offtarget_panel: { type: "string" }
    },
    required: ["primers", "target_panel"]
  }
}, {
  name: "search_pubmed",
  inputSchema: {
    type: "object",
    properties: {
      query: { type: "string" },
      max_results: { type: "integer", default: 20 },
      filters: { type: "object", properties: {
        publication_date: { type: "string" },
        journal: { type: "string" },
        species: { type: "string" }
      }}
    },
    required: ["query"]
  }
}, {
  name: "validate_primers_complete",
  inputSchema: {
    type: "object",
    properties: {
      primers: { type: "object" },
      target_panel: { type: "string" },
      offtarget_panel: { type: "string" },
      include_literature: { type: "boolean", default: true },
      literature_query: { type: "string" }
    },
    required: ["primers", "target_panel"]
  }
}]
```

**Consolidated Container:**
```dockerfile
FROM ubuntu:22.04
RUN apt-get update && apt-get install -y blast2 ncbi-blast+ python3-pip
RUN pip3 install gget biopython requests pandas
COPY validation_mcp_server.py /app/
EXPOSE 8004
CMD ["python3", "/app/validation_mcp_server.py"]
```

**Agent Validation Workflow:**
- **User**: "Validate these salmon primers and find related literature"
- **Assistant**: `validate_primers_complete(primers={...}, literature_query="Salmo salar primers COI")`
- **Result**: Complete validation + literature context

## Phase 6: Export & Provenance MCP Server

### 6.1 Results Export MCP Server
**Timeline: 1-2 weeks**
**Container**: `ndiag-export-server:latest`

**MCP Tools for Results Management:**
```typescript
tools: [{
  name: "export_panel",
  inputSchema: {
    type: "object",
    properties: {
      results_data: { type: "object" },
      format: { type: "string", enum: ["fasta", "csv", "json", "tsv"], default: "fasta" },
      include_metadata: { type: "boolean", default: true }
    },
    required: ["results_data"]
  }
}, {
  name: "generate_report",
  inputSchema: {
    type: "object",
    properties: {
      workflow_results: { type: "object" },
      report_type: { type: "string", enum: ["summary", "detailed", "validation"], default: "summary" },
      format: { type: "string", enum: ["pdf", "html", "markdown"], default: "html" }
    },
    required: ["workflow_results"]
  }
}, {
  name: "record_provenance",
  inputSchema: {
    type: "object",
    properties: {
      workflow_step: { type: "string" },
      inputs: { type: "object" },
      parameters: { type: "object" },
      outputs_uri: { type: "string" },
      timestamp: { type: "string" }
    },
    required: ["workflow_step", "inputs", "outputs_uri"]
  }
}]
```

**Container:**
```dockerfile
FROM python:3.11-slim
RUN pip install pandas reportlab jinja2 markdown
COPY export_mcp_server.py /app/
EXPOSE 8005
CMD ["python", "/app/export_mcp_server.py"]
```

**Agent Export Workflow:**
- **User**: "Export the final primer designs and validation report"
- **Assistant**: Complete export pipeline with full provenance tracking

## Implementation Timeline (6-Container MCP Architecture)

| Phase | Duration | Single MCP Container | Consolidates | Key Benefits |
|-------|----------|---------------------|-------------|--------------|
| Phase 1: Database Integration | 5-7 weeks | `ndiag-database-server` | NCBI, SRA, BOLD, SILVA, UNITE | Unified data access |
| Phase 2: Sequence Processing | 3 weeks | `ndiag-processing-server` | QC, dedup, masking, chimera | Complete processing pipeline |
| Phase 3: Alignment & Phylogenetics | 4-5 weeks | `ndiag-alignment-server` | MAFFT, MUSCLE, phylogenetics | Alignment + evolutionary analysis |
| Phase 4: Design & Primers | 7-9 weeks | `ndiag-design-server` | Region discovery + Primer3 | End-to-end primer design |
| Phase 5: Validation & Literature | 4-5 weeks | `ndiag-validation-server` | BLAST, in-silico PCR, PubMed | Comprehensive validation |
| Phase 6: Export & Provenance | 1-2 weeks | `ndiag-export-server` | Results, reports, provenance | Complete deliverables |

**Total Estimated Timeline: 24-31 weeks** (same capability, 6 containers instead of 12+)

## MCP Architecture Benefits

### For Users
- **Natural Language Interface**: "Find COI primers for salmon" → automated workflow
- **Real-time Progress**: See each MCP server's progress in the assistant
- **Error Recovery**: Failed services don't break the entire workflow
- **Reproducibility**: Full provenance tracking across all MCP calls

### For Developers  
- **Modular Development**: Each MCP server can be developed independently
- **Technology Flexibility**: Mix Python, R, CLI tools as appropriate per server
- **Scalability**: Scale individual services based on demand
- **Testing Isolation**: Test each MCP server independently

## MCP-Based Technical Architecture

```
neglected-diagnostics/
├── mcp_servers/                    # 6 Consolidated MCP Servers
│   ├── database_server/           # Phase 1: All Database Sources
│   │   ├── Dockerfile
│   │   ├── database_mcp_server.py
│   │   └── requirements.txt
│   ├── processing_server/         # Phase 2: Complete Sequence Processing
│   │   ├── Dockerfile
│   │   ├── processing_mcp_server.py
│   │   └── requirements.txt
│   ├── alignment_server/          # Phase 3: Alignment + Phylogenetics
│   │   ├── Dockerfile
│   │   ├── alignment_mcp_server.py
│   │   └── requirements.txt
│   ├── design_server/             # Phase 4: Region Discovery + Primer Design
│   │   ├── Dockerfile
│   │   ├── design_mcp_server.py
│   │   └── requirements.txt
│   ├── validation_server/         # Phase 5: BLAST + Literature
│   │   ├── Dockerfile
│   │   ├── validation_mcp_server.py
│   │   └── requirements.txt
│   └── export_server/             # Phase 6: Results + Provenance
│       ├── Dockerfile
│       ├── export_mcp_server.py
│       └── requirements.txt
├── src/app/assistant/             # Enhanced Assistant with MCP
│   ├── mcp_client.py             # MCP client integration
│   ├── workflow_orchestrator.py  # Workflow coordination
│   ├── core.py                   # Enhanced with MCP support
│   └── ...                       # Existing assistant components
├── docker-compose.yml             # All MCP servers + assistant
├── kubernetes/                    # K8s deployment configs
│   ├── mcp-servers.yaml
│   └── assistant-deployment.yaml
└── tests/                         # MCP Server Testing
    ├── test_mcp_servers/
    ├── test_workflows/
    └── integration_tests/
```

## Success Metrics (MCP-Enhanced)

1. **Database Coverage:**
   - Support for 4+ major databases via MCP servers (NCBI, BOLD, SILVA, UNITE)
   - Agent can seamlessly switch between data sources based on user needs

2. **Analysis Quality:**
   - Alignment accuracy >95% for test datasets via distributed MCP services
   - Design region specificity >90% with multi-server validation pipeline

3. **Performance & Scalability:**
   - Process 1000+ sequences across distributed MCP servers in <10 minutes
   - Individual MCP server response time <30 seconds
   - Horizontal scaling of resource-intensive servers (BLAST, alignment)

4. **Agent Usability:**
   - Natural language → complete workflow execution
   - Real-time progress tracking across all MCP servers
   - Intelligent error recovery and service fallbacks
   - Export compatibility with major primer design tools

5. **Developer Experience:**
   - Independent MCP server development and deployment
   - Comprehensive API documentation for each server
   - Integration testing across the full MCP ecosystem

## Minimal Viable Product (MVP) Stack

## Simplified 6-Container MVP

**Complete pipeline with optimal container count:**

1. **Database Server** (`ndiag-database-server`) - All sequence sources
2. **Processing Server** (`ndiag-processing-server`) - Complete sequence QC
3. **Alignment Server** (`ndiag-alignment-server`) - Alignment + phylogenetics
4. **Design Server** (`ndiag-design-server`) - Region discovery + Primer3
5. **Validation Server** (`ndiag-validation-server`) - BLAST + literature
6. **Export Server** (`ndiag-export-server`) - Results + provenance

**Docker Compose for 6-Container Architecture:**
```yaml
version: '3.8'
services:
  database-server:
    build: ./mcp_servers/database_server
    ports: ["8000:8000"]
    environment:
      - GOOGLE_APPLICATION_CREDENTIALS=/app/gcp-key.json  # For BigQuery/SRA access
  processing-server:
    build: ./mcp_servers/processing_server  
    ports: ["8001:8001"]
  alignment-server:
    build: ./mcp_servers/alignment_server
    ports: ["8002:8002"]
  design-server:
    build: ./mcp_servers/design_server
    ports: ["8003:8003"]
  validation-server:
    build: ./mcp_servers/validation_server
    ports: ["8004:8004"]
  export-server:
    build: ./mcp_servers/export_server
    ports: ["8005:8005"]
  assistant:
    build: ./src/app
    ports: ["8501:8501"]
    depends_on: [database-server, processing-server, alignment-server, design-server, validation-server, export-server]
    environment:
      - MCP_SERVERS=database:8000,processing:8001,alignment:8002,design:8003,validation:8004,export:8005
```

## Container Consolidation Benefits

### **Reduced Operational Overhead**
- **6 containers** instead of 12+ individual services
- **Simplified networking** - fewer inter-service communications
- **Easier monitoring** - 6 health checks vs dozens
- **Streamlined deployment** - single docker-compose command

### **Maintained Functionality**
- **All MCP tools preserved** - no feature loss from consolidation
- **Logical grouping** - related functions bundled together
- **Phase-aligned containers** - matches development workflow
- **Scalability retained** - can still scale individual phases as needed

## Enhanced Capabilities with SRA Integration

The addition of **SRA/BioProject MCP Server** brings powerful new capabilities to the agent:

### **Study-Driven Primer Design**
- **User**: "Design primers for salmon identification using data from recent studies"
- **Workflow**: SRA discovery → sequence extraction → alignment → primer design
- **Benefit**: Leverage real-world sequencing data for more robust primer validation

### **Meta-Analysis Capabilities**  
- **User**: "Compare primer performance across all published bacterial 16S amplicon studies"
- **Workflow**: Large-scale SRA mining → metadata analysis → performance comparison
- **Benefit**: Evidence-based primer selection using comprehensive datasets

### **Temporal Analysis**
- **User**: "How has salmon genetic diversity changed in studies over the last 5 years?"
- **Workflow**: Time-filtered SRA search → sequence extraction → diversity analysis
- **Benefit**: Track evolutionary trends and inform primer durability

### **Platform-Specific Optimization**
- **User**: "Optimize primers for Illumina vs PacBio sequencing platforms"
- **Workflow**: Platform-filtered SRA search → technology-specific validation
- **Benefit**: Platform-aware primer design for optimal performance

### **Cloud-Scale Discovery**
Using BigQuery/Athena integration, the agent can:
- Process metadata from millions of SRA runs
- Identify understudied taxonomic groups
- Find optimal reference datasets for primer validation
- Discover novel sequence variants affecting primer binding

## gget Integration Benefits

The integration of [gget](https://github.com/pachterlab/gget) brings significant advantages to the MCP architecture:

### **🎯 Standardized Database Access**
- **Unified API**: Single interface for Ensembl, NCBI, UniProt, PDB queries
- **Consistent Data Format**: Standardized outputs across all genomic databases  
- **Reliable Performance**: Well-tested, actively maintained genomic toolkit
- **Rich Metadata**: Comprehensive gene annotations and cross-references

### **⚡ Enhanced Capabilities**
- **Multi-Program BLAST**: Support for blastn, blastp, blastx, tblastn, tblastx
- **BLAT Integration**: Fast genomic location mapping
- **Protein Analysis**: AlphaFold structure prediction, motif finding (ELM)
- **Expression Data**: ARCHS4 tissue-specific gene expression
- **Enrichment Analysis**: Ontology and pathway analysis via Enrichr

### **🧬 Agent Workflow Improvements**

**Before gget:**
- **User**: "Find COI sequences for salmon"
- **Assistant**: Custom NCBI API calls → manual parsing → inconsistent formats

**With gget:**
- **User**: "Find COI sequences for salmon with protein structure info"
- **Assistant**: 
  1. `gget_search(["COI", "cytochrome oxidase"], "salmo_salar")` → standardized gene IDs
  2. `gget_seq(ens_ids=[...], translate=true)` → protein sequences
  3. `gget_pdb(pdb_ids=[...])` → 3D structures for primer design optimization
- **Result**: Rich, standardized data with protein context for better primer design

### **🔬 Advanced Use Cases Enabled**
- **Structure-Informed Design**: Use protein structures to avoid critical domains
- **Expression-Guided Selection**: Choose markers based on tissue expression patterns
- **Evolutionary Analysis**: Leverage Ensembl's comparative genomics data
- **Cross-Species Validation**: Consistent gene orthology across species

This **integrated MCP architecture** transforms the system from a basic sequence analysis tool into a **comprehensive, standardized genomic research platform** that leverages the full ecosystem of genomic databases through a unified, agent-friendly interface.