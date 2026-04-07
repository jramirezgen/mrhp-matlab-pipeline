# Shewanella xiamenensis LC6 — GEM XML

Place the SBML Level 3 FBC v2 genome-scale metabolic model here.

## Expected file

- **Format:** SBML Level 3 FBC v2 (`.xml`)
- **Organism:** *Shewanella xiamenensis* BC01 / LC6
- **Recommended name:** `iLC6_base_v*.xml`

## How to obtain

The GEM can be reconstructed from the annotated genome (GenBank: GCA_009684755.1)
using CarveMe, ModelSEED, or equivalent tools, then curated with the pipeline in:
`github.com/jramirezgen/mrhp-matlab-pipeline`

## Usage

The pipeline references this GEM for:
- Route map grounding (`configs/route_maps/`)
- Stoichiometric equation tables (output)
- Hypothesis validation against metabolic network topology
