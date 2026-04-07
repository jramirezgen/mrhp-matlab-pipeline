function cfg = config_acidithiobacillus_fe()
%CONFIG_ACIDITHIOBACILLUS_FE Configuration for A. ferrooxidans ATCC 23270.
%   Chemolithoautotrophic Fe2+ oxidation via downhill electron transfer.
%   2 conditions: CA (chemoautotrophic) and EA (electroautotrophic).
%   Phenotype: Fe2+ oxidation (%) for CA; cathodic current (uA) for EA.
%
%   Minimal mechanistic ODE: Fe2+ → Cyc2 → Rus → Cyc1 → aa3 → O2
%   Plus Calvin cycle CO2 fixation for biomass.

cfg.organism = 'Acidithiobacillus ferrooxidans ATCC 23270';

%% CONDITIONS
cfg.conditions = {'CA', 'EA'};

%% GENES (electron transfer chain + key DEGs)
cfg.genes = {'cyc2', 'rus', 'coxBA', 'pilA', 'tolC'};

%% SOLVER
cfg.solver.t_span   = [0, 72];
cfg.solver.dt_eval  = 0.5;
cfg.solver.t_eval   = (0:0.5:72)';
cfg.solver.rtol     = 1e-8;
cfg.solver.atol     = 1e-10;
cfg.solver.max_step = 0.5;
cfg.solver.substrate_conversion = 1.0;

%% SWEEP (Fe2+ concentration × O2 availability)
cfg.sweep.substrate_values = [6.1, 12.2, 18.3];  % Fe2+ g/L
cfg.sweep.target_values    = [0.21, 0.10];        % O2 partial pressure

%% REFERENCE CONDITION
cfg.ref_substrate = 12.2;
cfg.ref_target    = 0.21;
cfg.t_sample      = 48.0;

%% EXPRESSION LAYER
cfg.expression.beta_m  = 4.0;    % faster turnover in extremophile
cfg.expression.beta_p  = 0.5;
cfg.expression.ktl     = 3.0;
cfg.expression.stochastic_tau = 0.01;
cfg.expression.n_cells = 20;
cfg.expression.nm_to_mol = 6022;

%% TRANSCRIPTION RATES (estimated from TPM ratios)
cfg.ktx_fits.cyc2  = 0.80;   % highly expressed in CA
cfg.ktx_fits.rus   = 0.60;   % rusticyanin, abundant
cfg.ktx_fits.coxBA = 0.40;   % terminal oxidase
cfg.ktx_fits.pilA  = 0.10;   % low in CA, high in EA
cfg.ktx_fits.tolC  = 0.05;   % low in CA, high in EA

%% SPECIES (12 species for CA model)
% Fe2+, Fe3+, Cyc2_red, Cyc2_ox, Rus_red, Rus_ox, Cyc1_red, Cyc1_ox,
% O2, CO2, RuBP, Biomass
species_ca = {'Fe2','Fe3','Cyc2_red','Cyc2_ox','Rus_red','Rus_ox',...
              'Cyc1_red','Cyc1_ox','O2','CO2','RuBP','Biomass'};

%% ═══ CA MODEL: Chemoautotrophic Fe2+ oxidation ═══
% Reactions:
% R1: Fe2+ + Cyc2_ox → Fe3+ + Cyc2_red             (outer membrane oxidation)
% R2: Cyc2_red + Rus_ox → Cyc2_ox + Rus_red         (periplasm ET)
% R3: Rus_red + Cyc1_ox → Rus_ox + Cyc1_red          (periplasm ET)
% R4: 4 Cyc1_red + O2 → 4 Cyc1_ox                   (aa3 terminal oxidase)
% R5: CO2 + RuBP → Biomass (Calvin cycle simplified) 
% R6: Biomass → RuBP (recycling)
% R7: Fe3+ export (precipitation/removal)
% R8: Growth (logistic)

S_ca = [
    -1, 0, 0, 0, 0, 0, 0, 0;   % Fe2+
     1, 0, 0, 0, 0, 0,-1, 0;   % Fe3+
     1,-1, 0, 0, 0, 0, 0, 0;   % Cyc2_red
    -1, 1, 0, 0, 0, 0, 0, 0;   % Cyc2_ox
     0, 1,-1, 0, 0, 0, 0, 0;   % Rus_red
     0,-1, 1, 0, 0, 0, 0, 0;   % Rus_ox
     0, 0, 1,-4, 0, 0, 0, 0;   % Cyc1_red
     0, 0,-1, 4, 0, 0, 0, 0;   % Cyc1_ox
     0, 0, 0,-1, 0, 0, 0, 0;   % O2
     0, 0, 0, 0,-1, 0, 0, 0;   % CO2
     0, 0, 0, 0,-1, 1, 0, 0;   % RuBP
     0, 0, 0, 0, 1,-1, 0, 1;   % Biomass
];

m_ca.species = species_ca;
m_ca.S = S_ca;
m_ca.base_x0 = [12.2, 0.01, 0.05, 1.0, 0.05, 1.0, ...
                0.05, 1.0, 0.21, 0.04, 0.01, 0.01]';
m_ca.substrate_idx = 1;  % Fe2+
m_ca.target_idx    = 9;  % O2 (sweep target dimension)
m_ca.biomass_idx   = 12;

p_ca.k_cyc2     = 2.0;    % Fe2+ oxidation rate (scaled for 99% in 72h)
p_ca.Km_fe2     = 2.0;    % half-saturation for Fe2+
p_ca.k_et_rus   = 10.0;   % Cyc2→Rus electron transfer
p_ca.k_et_cyc1  = 10.0;   % Rus→Cyc1 electron transfer
p_ca.k_aa3      = 5.0;    % terminal oxidase rate
p_ca.Km_o2      = 0.05;   % O2 half-saturation
p_ca.k_calvin   = 0.5;    % Calvin cycle rate
p_ca.Km_co2     = 0.01;   % CO2 half-saturation
p_ca.k_rubp_rec = 0.2;    % RuBP recycling
p_ca.k_fe3_exp  = 0.05;   % Fe3+ export/precipitation
p_ca.mu_max     = 0.155;  % max growth rate (measured)
p_ca.Ks_growth  = 3.0;    % growth substrate half-sat
p_ca.Bmax       = 0.50;   % max biomass (OD scale)

m_ca.params = p_ca;
m_ca.rate_fn = @ca_rates;
cfg.models.CA = m_ca;

%% ═══ EA MODEL: Electroautotrophic (electrode electron uptake) ═══
% Same structure but electrons come from electrode instead of Fe2+
% Species: Electrode_e, Fe3+, Cyc2_red, Cyc2_ox, Rus_red, Rus_ox,
%          Cyc1_red, Cyc1_ox, O2, CO2, RuBP, Biomass, Pili, EPS

species_ea = {'Electrode_e','Fe3','Cyc2_red','Cyc2_ox','Rus_red','Rus_ox',...
              'Cyc1_red','Cyc1_ox','O2','CO2','RuBP','Biomass','Pili','EPS'};

% R1: Electrode + Cyc2_ox → Cyc2_red (EEU)
% R2: Cyc2_red + Rus_ox → Cyc2_ox + Rus_red
% R3: Rus_red + Cyc1_ox → Rus_ox + Cyc1_red
% R4: 4 Cyc1_red + O2 → 4 Cyc1_ox
% R5: CO2 + RuBP → Biomass
% R6: Biomass → RuBP
% R7: Biomass → Pili (attachment)
% R8: Biomass → EPS (biofilm)
% R9: Growth

S_ea = [
    -1, 0, 0, 0, 0, 0, 0, 0, 0;   % Electrode_e
     0, 0, 0, 0, 0, 0, 0, 0, 0;   % Fe3+ (not produced in EA)
     1,-1, 0, 0, 0, 0, 0, 0, 0;   % Cyc2_red
    -1, 1, 0, 0, 0, 0, 0, 0, 0;   % Cyc2_ox
     0, 1,-1, 0, 0, 0, 0, 0, 0;   % Rus_red
     0,-1, 1, 0, 0, 0, 0, 0, 0;   % Rus_ox
     0, 0, 1,-4, 0, 0, 0, 0, 0;   % Cyc1_red
     0, 0,-1, 4, 0, 0, 0, 0, 0;   % Cyc1_ox
     0, 0, 0,-1, 0, 0, 0, 0, 0;   % O2
     0, 0, 0, 0,-1, 0, 0, 0, 0;   % CO2
     0, 0, 0, 0,-1, 1, 0, 0, 0;   % RuBP
     0, 0, 0, 0, 1,-1,-1,-1, 1;   % Biomass
     0, 0, 0, 0, 0, 0, 1, 0, 0;   % Pili
     0, 0, 0, 0, 0, 0, 0, 1, 0;   % EPS
];

m_ea.species = species_ea;
m_ea.S = S_ea;
m_ea.base_x0 = [5.0, 0, 0.05, 1.0, 0.05, 1.0, ...
                0.05, 1.0, 0.21, 0.04, 0.01, 0.01, 0.001, 0.001]';
m_ea.substrate_idx = 1;  % Electrode electrons
m_ea.target_idx    = 9;  % O2 (sweep target dimension)
m_ea.biomass_idx   = 12;

p_ea = p_ca;
p_ea.k_cyc2    = 0.5;    % slower electrode uptake
p_ea.Km_fe2    = 0.5;    % electrode affinity
p_ea.mu_max    = 0.061;  % measured EA growth rate
p_ea.k_pili    = 0.02;   % pili biogenesis
p_ea.k_eps     = 0.01;   % EPS production
p_ea.Bmax      = 0.20;   % lower max biomass in EA
p_ea.k_rubp_rec = 0.02;  % lower recycling to allow net growth

m_ea.params = p_ea;
m_ea.rate_fn = @ea_rates;
cfg.models.EA = m_ea;

%% PHI DEFINITION
cfg.phi_definition.CA = struct('type','depletion', 'target_idx',1);
cfg.phi_definition.EA = struct('type','accumulation', 'target_idx',12, 'p_max', 0.20);

%% PHENOTYPE REFERENCE
cfg.phenotype_label = 'Fe2+ oxidation (%)';
cfg.phenotype_ref.CA = struct('type','sigmoidal_accumulation',...
    'Pmax',100, 'k',0.10, 'lag',24, 'n',1);
cfg.phenotype_ref.EA = struct('type','sigmoidal_accumulation',...
    'Pmax',75, 'k',0.05, 'lag',12, 'n',1);  % current uptake (uA scale)

%% REGULATORY SIGNALS
% cyc2: activation by Fe2+ availability (CA) or electrode (EA)
cfg.regulatory.cyc2.CA = struct('form','act',...
    'signals',{{1}}, 'K',2.0, 'n',2, 'basal',0.05);
cfg.regulatory.cyc2.EA = struct('form','act',...
    'signals',{{1}}, 'K',0.1, 'n',1, 'basal',0.10);
% rus: activation by reduced Cyc2
cfg.regulatory.rus.default = struct('form','act',...
    'signals',{{3}}, 'K',0.005, 'n',2, 'basal',0.10);
% coxBA: activation by Cyc1_red (electron pressure)
cfg.regulatory.coxBA.default = struct('form','act',...
    'signals',{{7}}, 'K',0.005, 'n',2, 'basal',0.05);
% pilA: low in CA, high in EA (constitutive in EA)
cfg.regulatory.pilA.CA = struct('form','custom', 'basal',0.01,...
    'fn', @(y,d,r) r.basal * ones(1, size(y,2)));
cfg.regulatory.pilA.EA = struct('form','custom', 'basal',0.50,...
    'fn', @(y,d,r) r.basal * ones(1, size(y,2)));
% tolC: same pattern as pilA
cfg.regulatory.tolC.CA = struct('form','custom', 'basal',0.01,...
    'fn', @(y,d,r) r.basal * ones(1, size(y,2)));
cfg.regulatory.tolC.EA = struct('form','custom', 'basal',0.30,...
    'fn', @(y,d,r) r.basal * ones(1, size(y,2)));

%% EXPERIMENTAL DATA (from Wang et al. 2024)
cfg.experimental_data.CA.fe2_oxidation_pct = 99.18;  % at 72h
cfg.experimental_data.CA.doubling_time_h   = 4.46;
cfg.experimental_data.CA.cells_mm3_72h     = 51000;
cfg.experimental_data.EA.current_uA_60h    = -75;
cfg.experimental_data.EA.doubling_time_h   = 11.35;

%% HYPOTHESES
cfg.hypotheses.CA = { ...
  struct('id','CA_downhill','label','Downhill ET via Cyc2-Rus-Cyc1-aa3',...
    'status','BEST', 'diagnostic_up',{{'cyc2','rus','coxBA'}}, 'score_total',0.90), ...
  struct('id','CA_reverse','label','Reverse ET for NADH generation',...
    'status','SUPPORTED', 'diagnostic_up',{{'cyc2','rus'}}, 'score_total',0.75) };
cfg.hypotheses.EA = { ...
  struct('id','EA_direct_EEU','label','Direct electrode uptake via Cyc2',...
    'status','BEST', 'diagnostic_up',{{'cyc2','pilA','tolC'}}, 'score_total',0.85), ...
  struct('id','EA_biofilm','label','Biofilm-mediated electron transfer',...
    'status','SUPPORTED', 'diagnostic_up',{{'pilA','tolC'}}, 'score_total',0.80) };

%% VISUALIZATION
cfg.visualization.colors.CA = [0.800, 0.400, 0.100];  % rust orange
cfg.visualization.colors.EA = [0.200, 0.500, 0.700];  % electrode blue
cfg.visualization.labels.CA = 'Chemoautotrophic (Fe2+)';
cfg.visualization.labels.EA = 'Electroautotrophic';

%% OUTPUT DIRECTORY
cfg.output_dir = '';
end

%% ═══ CA RATE FUNCTIONS ═══
function v = ca_rates(x, p)
    x = max(x, 0);
    Fe2=x(1); Cyc2_ox=x(4); Cyc2_red=x(3);
    Rus_ox=x(6); Rus_red=x(5); Cyc1_ox=x(8); Cyc1_red=x(7);
    O2=x(9); CO2=x(10); RuBP=x(11); B=max(x(12), 1e-12);
    Fe3=x(2);
    v = zeros(8,1);
    % R1: Fe2+ oxidation by Cyc2
    v(1) = p.k_cyc2 * Fe2 / (p.Km_fe2 + Fe2) * Cyc2_ox * B;
    % R2: Cyc2 → Rus ET
    v(2) = p.k_et_rus * Cyc2_red * Rus_ox;
    % R3: Rus → Cyc1 ET
    v(3) = p.k_et_cyc1 * Rus_red * Cyc1_ox;
    % R4: aa3 terminal oxidase (4 Cyc1_red + O2)
    v(4) = p.k_aa3 * Cyc1_red * O2 / (p.Km_o2 + O2);
    % R5: Calvin cycle
    v(5) = p.k_calvin * CO2 / (p.Km_co2 + CO2) * RuBP * B;
    % R6: RuBP recycling
    v(6) = p.k_rubp_rec * B;
    % R7: Fe3+ export
    v(7) = p.k_fe3_exp * Fe3;
    % R8: Growth
    v(8) = p.mu_max * B * Fe2 / (p.Ks_growth + Fe2) * max(0, 1 - B / p.Bmax);
end

%% ═══ EA RATE FUNCTIONS ═══
function v = ea_rates(x, p)
    x = max(x, 0);
    Elec=x(1); Cyc2_ox=x(4); Cyc2_red=x(3);
    Rus_ox=x(6); Rus_red=x(5); Cyc1_ox=x(8); Cyc1_red=x(7);
    O2=x(9); CO2=x(10); RuBP=x(11); B=max(x(12), 1e-12);
    v = zeros(9,1);
    % R1: Electrode → Cyc2 (EEU)
    v(1) = p.k_cyc2 * Elec / (p.Km_fe2 + Elec) * Cyc2_ox * B;
    % R2: Cyc2 → Rus ET
    v(2) = p.k_et_rus * Cyc2_red * Rus_ox;
    % R3: Rus → Cyc1 ET
    v(3) = p.k_et_cyc1 * Rus_red * Cyc1_ox;
    % R4: aa3 terminal oxidase
    v(4) = p.k_aa3 * Cyc1_red * O2 / (p.Km_o2 + O2);
    % R5: Calvin cycle
    v(5) = p.k_calvin * CO2 / (p.Km_co2 + CO2) * RuBP * B;
    % R6: RuBP recycling
    v(6) = p.k_rubp_rec * B;
    % R7: Pili biogenesis
    v(7) = p.k_pili * B;
    % R8: EPS production
    v(8) = p.k_eps * B;
    % R9: Growth
    v(9) = p.mu_max * B * Elec / (p.Ks_growth + Elec) * max(0, 1 - B / p.Bmax);
end