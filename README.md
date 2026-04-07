# MRHP — Multiscale Route Hypothesis Platform

**v4.0.0** · MATLAB R2025b · MIT License

A universal MATLAB pipeline for multiscale metabolic analysis: from ODE-based metabolic kinetics through phenotypic bridge coupling to absolute gene expression dynamics (deterministic + stochastic).

## Certified Organisms

| Organism | Conditions | Figures | Runtime | Bridge R² |
|----------|-----------|---------|---------|----------|
| *Shewanella xiamenensis* LC6 | MO, RC, AD | 12 | ~106 s | 0.993 |
| *E. coli* (fucose) | FUC1–FUC12NpG | 13 | ~197 s | 0.916 |
| *Acidithiobacillus ferrooxidans* | CA, EA | 12 | ~51 s | 0.991 |

**Total: 3/3 CERTIFIED · 37 figures · ~356 s**

## Core Equation

The phenotypic bridge links metabolic flux signals to observable phenotype:

$$y(t) = y_{\max} \cdot \left(1 - e^{-H(t)^{\beta}}\right), \quad H(t) = \int_0^t \lambda \cdot \frac{d\Phi}{d\tau}\, d\tau$$

Where Φ(t) is the metabolic activity signal extracted from ODE dynamics, and (λ, β, y_max) are organism/condition-specific parameters fitted via multi-start global optimization.

## Pipeline Architecture

| Phase | Description |
|-------|-------------|
| 0 | Config verification |
| 1 | ODE metabolic simulation |
| 2 | Φ signal extraction |
| 3 | Bridge fitting (50-start global) |
| 4 | Regulatory signal computation |
| 5 | Deterministic gene expression |
| 6 | Stochastic gene expression (tau-leaping) |
| 7 | Full parameter grid |
| 8 | Figure generation (13 types) |
| 9 | Statistical validation |
| 10 | Metadata export |
| 10c | Route-map grounding (GEM traceability) |
| 10d | Coupled expression dynamics |
| 11 | Markdown report |
| 12 | Certification |

## Quick Start

```matlab
% Run all organisms (quick mode)
run_all_systems('quick')

% Run all organisms (full mode with all figures)
run_all_systems('full')

% Run a single organism
addpath('engine', 'configs');
cfg = config_shewanella_lc6();
run_pipeline_generic(cfg);
```

## Directory Structure

```
mrhp-matlab-pipeline/
├── engine/                        % Core pipeline functions
│   ├── run_pipeline_generic.m     % 13-phase orchestrator
│   ├── figure_generation_generic.m% 13 figure types
│   ├── solve_ode_generic.m        % ODE metabolic solver
│   ├── extract_phi_generic.m      % Φ signal extraction
│   ├── bridge_from_phi.m          % Bridge evaluator
│   ├── fit_bridge.m               % 50-start global fitting
│   ├── compute_u_generic.m        % Regulatory signals
│   ├── solve_deterministic_expression.m
│   ├── tau_leaping_expression.m   % Stochastic expression
│   ├── hypothesis_evaluator_generic.m
│   └── phenotype_reference_generic.m
├── configs/                       % Organism configurations
│   ├── config_shewanella_lc6.m
│   ├── config_ecoli_fucose.m
│   ├── config_acidithiobacillus_fe.m
│   └── route_maps/                % GEM-derived stoichiometry
│       ├── route_map_{MO,RC,AD}.tsv
│       ├── species_map_{MO,RC,AD}.tsv
│       └── S_matrix_{MO,RC,AD}.tsv
├── tests/
│   └── run_statistical_audit.m    % Comprehensive validation
├── run_all_systems.m              % Master entry point
├── run_all.m                      % Environment-based runner
├── CHANGELOG.md
├── CITATION.cff
└── LICENSE
```

## Figure Types

The pipeline generates up to 13 publication-ready figures per organism:

1. **ODE Dynamics** — Metabolite concentration profiles
2. **Φ Signal** — Metabolic activity extraction
3. **Bridge Fit** — Phenotypic coupling with R²
4. **Bridge Residuals** — Fit quality diagnostics
5. **Regulatory Signals** — u(t) per gene
6. **Deterministic Expression** — mRNA/protein trajectories
7. **Stochastic Expression** — Tau-leaping ensemble
8. **Stochastic Variance** — Noise characterization
9. **Hypothesis Scores** — 5-component radar plot
10. **Coupling Dynamics** — Coupled expression trajectories
11. **Route Map Network** — GEM-traced reaction topology
12. **Condition Comparison** — Cross-condition bridge overlay
13. **Summary Dashboard** — Certification panel

## Route Maps (GEM Traceability)

Each condition includes TSV-format route maps derived from genome-scale models:

- `route_map_*.tsv` — Reaction IDs, names, equations, EC numbers, genes
- `species_map_*.tsv` — Metabolite IDs, names, compartments, initial concentrations
- `S_matrix_*.tsv` — Stoichiometric matrix linking reactions to species

These ensure every ODE species and reaction traces back to a curated GEM reconstruction.

## Requirements

- MATLAB R2021b or later (tested on R2025b)
- No additional toolboxes required
- No external dependencies

## Citation

```bibtex
@software{ramirez2026mrhp,
  author  = {Ramirez-Bautista, Josue},
  title   = {MRHP: Multiscale Route Hypothesis Platform},
  version = {4.0.0},
  year    = {2026},
  url     = {https://github.com/jramirezgen/mrhp-matlab-pipeline}
}
```

## License

MIT — see [LICENSE](LICENSE)
