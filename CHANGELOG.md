# Changelog

## [4.1.0] — 2026-04-07

### Added
- `inputs/` directory with all experimental data for full reproducibility
- Shewanella: decolorization kinetics (~4670 rows), RT-qPCR (6 genes × 4 conditions), expression parameter fits, DB71 reaction catalog (2930 reactions)
- E. coli: L-fucose titers (12 strains, Xia 2025)
- Acidithiobacillus: Fe²⁺ oxidation timeseries (CA/EA), differentially expressed genes
- GEM placeholder directories with instructions per organism
- Updated README with inputs documentation

### Fixed
- `run_all.m`: corrected base_dir resolution (was using double fileparts for scripts/ layout)
- `run_all_systems.m`: corrected out_dir path (removed `..` prefix)

## [4.0.0] — 2026-04-07

### Added
- Universal MATLAB pipeline engine with 13-phase orchestration
- Three organism configurations: Shewanella LC6, E. coli BL21 L-fucose, A. ferrooxidans ATCC 23270
- 13 publication-quality figure types per organism
- GEM traceability via route maps (ModelSEED cpd/rxn identifiers)
- Stoichiometric matrix heatmap (Fig 11), metabolic intermediates (Fig 12), flux profiles (Fig 13)
- Phase 10c: automated route map export with stoichiometric equations
- Phase 10d: S matrix numeric export per condition
- Option D Hybrid coupling for E. coli (18 → 34 species, mRNA + protein as synchronized passengers)
- Stochastic tau-leaping expression layer with deterministic cross-validation
- Hypothesis evaluation with 5-component scoring (R², direction, magnitude, genomic support, parsimony)
- Statistical audit framework (`tests/run_statistical_audit.m`)

### Certified Results
- Shewanella xiamenensis LC6: R² = 0.9930, 12 figures, 106 s
- E. coli BL21(DE3) L-Fucose: R² = 0.9157, 13 figures, 197 s
- Acidithiobacillus ferrooxidans: R² = 0.9912, 12 figures, 51 s
- Total: 78 scenarios, 37 figures, 355.7 s
