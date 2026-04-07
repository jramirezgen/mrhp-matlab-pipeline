function cfg = config_ecoli_fucose()
%CONFIG_ECOLI_FUCOSE Configuration for E. coli BL21(DE3) L-Fucose MRHP.
%   Engineered 2'-FL -> L-fucose pathway (Xia et al. 2025, J. Agric. Food Chem.).
%   5 conditions = progressive strain optimization with distinct genetic backgrounds.
%   Michaelis-Menten kinetics with cofactors (GTP, NADPH), growth model, and
%   condition-specific knockouts/additions via parameter vector.
%
%   Pathway: F6P -> Man6P -> Man1P -> GDP-Man -> GDP-4keto -> GDP-Fuc
%            GDP-Fuc + Lactose -> 2'-FL + GDP   (wbgL)
%            2'-FL -> L-Fucose + Lactose_recycle (afcA)
%
%   Species (18): Metabolites + cofactors + carbon sources + biomass
%   Reactions (20): Michaelis-Menten, Monod growth, transport, cofactor regeneration

cfg.organism = 'Escherichia coli BL21(DE3) L-Fucose';

%% CONDITIONS (5 distinct genetic backgrounds)
cfg.conditions = {'FUC1', 'FUC8', 'FUC12', 'FUC12NSM', 'FUC12NpG'};

%% GENES (pathway + engineering targets)
cfg.genes = {'manA', 'manB', 'manC', 'gmd', 'wcaG', 'wbgL', 'afcA', 'ndk'};

%% SOLVER
cfg.solver.t_span   = [0, 150];
cfg.solver.dt_eval  = 0.5;
cfg.solver.t_eval   = (0:0.5:150)';
cfg.solver.rtol     = 1e-6;
cfg.solver.atol     = 1e-9;
cfg.solver.max_step = 0.5;
cfg.solver.substrate_conversion = 1.0;

%% SWEEP (lactose feed mM x glycerol supply mM)
cfg.sweep.substrate_values = [15, 30, 45];   % lcts_ext (mM), ~5-15 g/L
cfg.sweep.target_values    = [0, 50];         % glycerol_ext (mM)

%% REFERENCE CONDITION
cfg.ref_substrate = 30;     % 30 mM lactose (~10 g/L)
cfg.ref_target    = 0;      % no glycerol (default)
cfg.t_sample      = 48.0;   % shake-flask sampling time

%% COUPLED EXPRESSION MODE
cfg.coupled_expression = true;  % mRNA/protein as ODE state variables

%% EXPRESSION LAYER
cfg.expression.beta_m  = 8.32;
cfg.expression.beta_p  = 0.924;
cfg.expression.ktl     = 6.0;
cfg.expression.stochastic_tau = 0.01;
cfg.expression.n_cells = 20;
cfg.expression.nm_to_mol = 6022;

%% TRANSCRIPTION RATES (constitutive T7; proportional to avg copy number)
cfg.ktx_fits.manA = 0.20;   % native, single copy
cfg.ktx_fits.manB = 0.25;   % plasmid, single copy
cfg.ktx_fits.manC = 0.25;   % plasmid, single copy
cfg.ktx_fits.gmd  = 0.50;   % CBGW cassette, strong expression
cfg.ktx_fits.wcaG = 0.50;   % CBGW cassette
cfg.ktx_fits.wbgL = 0.40;   % heterologous, strong
cfg.ktx_fits.afcA = 0.35;   % heterologous, strong
cfg.ktx_fits.ndk  = 0.15;   % overexpressed only in N+ strains

%% =====================================================================
%%   SPECIES LIST (18 metabolites)
%% =====================================================================
% Idx  Name           Description                              Units (mM)
%  1   F6P_c          Fructose-6-phosphate                     ~1-5
%  2   man6p_c        Mannose-6-phosphate                      ~0.1-1
%  3   man1p_c        Mannose-1-phosphate                      ~0.1-0.5
%  4   gdpman_c       GDP-D-mannose                            ~0.05-0.5
%  5   gdp4keto_c     GDP-4-keto-6-deoxymannose                ~0.01-0.1
%  6   gdpfuc_c       GDP-L-fucose                             ~0.01-0.1
%  7   lcts_c         Lactose intracellular                    ~1-10
%  8   fl2_c          2'-Fucosyllactose intracellular           ~0-40
%  9   fuc_L_c        L-Fucose intracellular                   ~0-20
% 10   fuc_L_e        L-Fucose extracellular (TARGET)          0->14-560
% 11   fl2_e          2'-FL extracellular (byproduct)          0->32-131
% 12   GTP_c          GTP pool                                 ~1-5
% 13   GDP_c          GDP pool                                 ~0.5-3
% 14   NADPH_c        NADPH pool                               ~0.1-0.5
% 15   biomass        Biomass (OD600 proxy)                    0.1->5-205
% 16   glucose_ext    Glucose in medium                        ~50-200
% 17   glycerol_ext   Glycerol in medium                       0 or ~100
% 18   lcts_ext       Lactose in medium (inducer/substrate)    ~15-44

species_list = {'F6P_c','man6p_c','man1p_c','gdpman_c',...
    'gdp4keto_c','gdpfuc_c','lcts_c','fl2_c',...
    'fuc_L_c','fuc_L_e','fl2_e','GTP_c',...
    'GDP_c','NADPH_c','biomass','glucose_ext',...
    'glycerol_ext','lcts_ext'};

%% =====================================================================
%%   STOICHIOMETRIC MATRIX (18 species x 20 reactions)
%% =====================================================================
% Reactions:
%  1  r_manA       F6P -> Man6P                             (manA, native)
%  2  r_manB       Man6P -> Man1P                           (manB)
%  3  r_manC       Man1P + GTP -> GDPMan + GDP              (manC, rate-limiting)
%  4  r_gmd        GDPMan -> GDP4keto                       (gmd, copy-scaled)
%  5  r_wcaG       GDP4keto + NADPH -> GDPFuc + NADP+       (wcaG, copy-scaled)
%  6  r_wbgL       GDPFuc + Lcts -> 2FL + GDP               (wbgL, copy-scaled)
%  7  r_afcA       2FL -> FucL + Lcts_recycle                (afcA, copy-scaled)
%  8  r_fuc_export FucL_c -> FucL_e                         (cellular export)
%  9  r_mdfA       2FL_c -> 2FL_e                           (mdfA transporter; KO in NSM/NpG)
% 10  r_setA       2FL_c -> 2FL_e                           (setA transporter; KO in NSM/NpG)
% 11  r_ndk        GDP -> GTP                               (ndk; active only in N+ strains)
% 12  r_growth     glucose -> biomass                       (Monod x logistic)
% 13  r_pfkA       F6P -> drain                             (glycolysis; KO in NpG)
% 14  r_glpK       glycerol -> F6P                          (glpK*; only in NpG)
% 15  r_glc_to_F6P glucose -> F6P                           (glycolysis upstream, growth-coupled)
% 16  r_lcts_import lcts_ext -> lcts_c                      (lactose permease)
% 17  r_NADPH_regen -> NADPH                                (pentose phosphate, growth-coupled)
% 18  r_GTP_base   -> GTP                                   (basal nucleotide biosynthesis)
% 19  r_glc_feed   -> glucose_ext                           (continuous glucose feeding)
% 20  r_gly_feed   -> glycerol_ext                          (glycerol co-feed; NpG only)

%        r1  r2  r3  r4  r5  r6  r7  r8  r9 r10 r11 r12 r13 r14 r15 r16 r17 r18 r19 r20
S = [
    -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  1,  1,  0,  0,  0,  0,  0;  %  1 F6P_c
     1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  %  2 man6p_c
     0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  %  3 man1p_c
     0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  %  4 gdpman_c
     0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  %  5 gdp4keto_c
     0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  %  6 gdpfuc_c
     0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0;  %  7 lcts_c
     0,  0,  0,  0,  0,  1, -1,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  %  8 fl2_c
     0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  %  9 fuc_L_c
     0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  % 10 fuc_L_e
     0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;  % 11 fl2_e
     0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0;  % 12 GTP_c
     0,  0,  1,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0;  % 13 GDP_c
     0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0;  % 14 NADPH_c
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0;  % 15 biomass
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0, -1,  0,  0,  0,  1,  0;  % 16 glucose_ext
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1;  % 17 glycerol_ext
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0;  % 18 lcts_ext
];

%% =====================================================================
%%   KINETIC PARAMETERS (base values -- overridden per condition)
%% =====================================================================
% Michaelis-Menten Vmax/Km; copy-number scaling via n_gene in params struct
p_base = struct();
% Upstream pathway (native/plasmid, single copy)
p_base.Vm_manA = 0.15;   p_base.Km_manA = 1.5;    % native high activity - NOT rate-limiting
p_base.Vm_manB = 0.80;   p_base.Km_manB = 0.5;    % plasmid expressed
p_base.Vm_manC = 0.40;   p_base.Km_manC_m1p = 0.3; p_base.Km_manC_gtp = 0.05;  % rate-limiting bi-substrate
% Core engineered pathway (copy-number x Vmax)
p_base.Vm_gmd  = 0.35;   p_base.Km_gmd  = 0.2;    % GDPMan -> GDP4keto
p_base.Vm_wcaG = 0.30;   p_base.Km_wcaG_sub = 0.15; p_base.Km_wcaG_nadph = 0.02;  % bi-substrate, NADPH
p_base.Vm_wbgL = 0.25;   p_base.Km_wbgL_gfuc = 0.1; p_base.Km_wbgL_lcts = 2.0;   % bi-substrate
p_base.Vm_afcA = 0.40;   p_base.Km_afcA = 1.0;    % 2FL -> FucL + Lcts (high Vcat)
% Transport
p_base.k_fuc_export = 0.08;  % L-fucose export (first-order)
p_base.k_mdfA  = 0.015;      % 2'-FL export via mdfA
p_base.k_setA  = 0.010;      % 2'-FL export via setA
% GTP/cofactor
p_base.Vm_ndk  = 0.20;   p_base.Km_ndk = 0.5;     % GDP -> GTP (ndk)
% Growth (Monod x logistic)
p_base.mu_max  = 0.45;  % h^-1 (E. coli BL21)
p_base.Ks_glc  = 0.5;   % mM glucose half-saturation
p_base.Bmax    = 15;     % OD carrying capacity (shake-flask default)
% Glycolysis / Carbon routing
p_base.Vm_pfkA = 0.15;   p_base.Km_pfkA = 2.0;    % F6P drain via glycolysis
p_base.Vm_glpK = 0.0;   p_base.Km_glpK = 5.0;    % glycerol -> F6P (inactive by default)
% Glucose -> F6P supply (growth-coupled)
p_base.k_glc_F6P = 0.060;     % rate constant x [glucose] x [biomass]
% Lactose import
p_base.k_lcts_import = 0.3;  % first-order on lcts_ext
% Cofactor regeneration (growth-coupled)
p_base.k_NADPH_regen = 0.30; % NADPH regeneration per biomass
p_base.k_GTP_base    = 0.25; % basal GTP biosynthesis per biomass
p_base.pathway_pull = 1.0;  % scaling factor for F6P supply based on gene dosage
p_base.expr_factor = 1.0;   % metabolic efficiency factor (expression + recycling)
p_base.k_glc_feed = 2.0;     % glucose continuous feed
p_base.k_gly_feed = 0.0;     % glycerol feed (0 by default, >0 for NpG)
% Copy numbers (1x default)
p_base.n_gmd  = 1;
p_base.n_wcaG = 1;
p_base.n_wbgL = 1;
p_base.n_afcA = 1;
p_base.n_manA = 1;
p_base.n_manB = 1;
p_base.n_manC = 1;
%% =====================================================================
m1.species = species_list;
m1.S = S;
m1.base_x0 = [2.0, 0.1, 0.05, 0.02, ...   % F6P, man6p, man1p, gdpman
               0.01, 0.01, 5.0, 0, ...      % gdp4keto, gdpfuc, lcts_c, fl2_c
               0, 0, 0, 2.0, ...             % fuc_L_c, fuc_L_e, fl2_e, GTP
               1.0, 0.3, 0.1, 100, ...       % GDP, NADPH, biomass(OD), glucose_ext
               0, 30]';                      % glycerol_ext=0, lcts_ext=30mM
m1.substrate_idx = 18;  % lcts_ext (sweep dimension)
m1.target_idx    = 17;  % glycerol_ext (sweep dimension)
m1.biomass_idx   = 15;  % biomass
p1 = p_base;
p1.n_gmd = 1; p1.n_wcaG = 1; p1.n_wbgL = 1; p1.n_afcA = 1;
p1.n_manA = 1; p1.n_manB = 1; p1.n_manC = 1;  % single-copy plasmid
p1.pathway_pull = 1.0;    % single-copy: baseline demand
p1.expr_factor = 1.0;
p1.Vm_glpK = 0;        % no glpK
p1.Bmax = 12;           % shake-flask carrying capacity
p1.k_glc_feed = 1.5;   % moderate glucose feed (shake-flask)
m1.params = p1;
m1.rate_fn = @ecoli_mm_rates;
cfg.models.FUC1 = m1;

%% =====================================================================
%%   FUC8: Gene dosage optimized (4:4:3:1, shake-flask, 6.74 g/L at 48h)
%% =====================================================================
m8 = m1;
p8 = p1;
p8.n_gmd = 4; p8.n_wcaG = 4; p8.n_wbgL = 3; p8.n_afcA = 1;
p8.pathway_pull = 3.5;    % optimized for 41 mM target
p8.expr_factor = 1.0;
p8.n_manA = 1; p8.n_manB = 4; p8.n_manC = 4;  % manA native, manB/C on high-copy plasmid
m8.params = p8;
cfg.models.FUC8 = m8;
%% =====================================================================
%%   FUC12: Balanced optimal (4:4:3:2, bioreactor, 63.77 g/L at 138h)
%% =====================================================================
m12 = m1;
p12 = p_base;
p12.n_gmd = 4; p12.n_wcaG = 4; p12.n_wbgL = 3; p12.n_afcA = 2;
p12.pathway_pull = 3.0;   % same gene dosage as FUC8
p12.expr_factor = 1.0;
p12.n_manA = 1; p12.n_manB = 4; p12.n_manC = 4; % same plasmid backbone
p12.Vm_ndk = 0.20;          % ndk overexpressed
p12.Vm_glpK = 0;           % no glpK
p12.Bmax = 35;              % bioreactor carrying capacity (OD~32)
p12.k_glc_feed = 5.0;      % bioreactor fed-batch glucose feed
m12.base_x0 = m1.base_x0;
m12.base_x0(15) = 0.1;     % biomass start
m12.base_x0(16) = 150;     % glucose_ext = 150 mM (bioreactor)
m12.base_x0(18) = 44;      % lcts_ext = 44 mM (~15 g/L bolus)
m12.params = p12;
m12.rate_fn = @ecoli_mm_rates;
cfg.models.FUC12 = m12;

%% =====================================================================
%%   FUC12NSM: Transporter KOs (mdfA+setA KO, bioreactor, 80.14 g/L)
%% =====================================================================
mNSM = m12;
pNSM = p12;
pNSM.k_mdfA = 0;           % mdfA knocked out
pNSM.k_setA = 0;           % setA knocked out
pNSM.expr_factor = 2.8;  % higher metabolic efficiency: lactose recycling + no 2FL export loss
pNSM.k_glc_feed = 5;     % bioreactor fed-batch
mNSM.params = pNSM;
cfg.models.FUC12NSM = mNSM;

%% =====================================================================
%%   FUC12NpG: Carbon routing (pfkA KO + glpK* + glycerol, 91.90 g/L)
%% =====================================================================
mNpG = m12;
pNpG = p12;
pNpG.k_mdfA  = 0;          % mdfA knocked out (inherited from NSM)
pNpG.k_setA  = 0;          % setA knocked out (inherited from NSM)
pNpG.Vm_pfkA = 0;          % pfkA knocked out -> no F6P drain via glycolysis
pNpG.Vm_glpK = 0.25;        % glpK* G913A mutant (feedback-resistant)
pNpG.Km_glpK = 3.0;        % active glycerol kinase
pNpG.k_gly_feed = 1.8;     % glycerol co-feed at 1:5 ratio
pNpG.Bmax = 210;            % OD600=205.40 (massive growth on dual carbon)
pNpG.mu_max = 0.55;        % enhanced growth on glycerol+glucose
pNpG.k_glc_feed = 4.8;     % higher glucose feed (bioreactor, glycerol:glucose 1:5)
mNpG.base_x0 = m12.base_x0;
mNpG.base_x0(17) = 100;    % glycerol_ext = 100 mM (10 g/L initial)
mNpG.params = pNpG;
cfg.models.FUC12NpG = mNpG;

%% =====================================================================
%%   PHI DEFINITION (per-condition accumulation, calibrated to ODE scale)
%% =====================================================================
% fuc_L_e is species 10; p_max calibrated to expected ODE accumulation
cfg.phi_definition.FUC1     = struct('type','accumulation','target_idx',10, 'p_max', 20);
cfg.phi_definition.FUC8     = struct('type','accumulation','target_idx',10, 'p_max', 55);
cfg.phi_definition.FUC12    = struct('type','accumulation','target_idx',10, 'p_max', 500);
cfg.phi_definition.FUC12NSM = struct('type','accumulation','target_idx',10, 'p_max', 620);
cfg.phi_definition.FUC12NpG = struct('type','accumulation','target_idx',10, 'p_max', 700);

%% PHENOTYPE REFERENCE (sigmoidal accumulation, mg/L)
cfg.phenotype_label = 'L-Fucose titer (mg/L)';
cfg.phenotype_ref.FUC1     = struct('type','sigmoidal_accumulation',...
    'Pmax', 2250,  'k', 0.12, 'lag', 4, 'n', 2);
cfg.phenotype_ref.FUC8     = struct('type','sigmoidal_accumulation',...
    'Pmax', 6740,  'k', 0.10, 'lag', 3, 'n', 2);
cfg.phenotype_ref.FUC12    = struct('type','sigmoidal_accumulation',...
    'Pmax', 63770, 'k', 0.04, 'lag', 6, 'n', 2);
cfg.phenotype_ref.FUC12NSM = struct('type','sigmoidal_accumulation',...
    'Pmax', 80140, 'k', 0.05, 'lag', 5, 'n', 2);
cfg.phenotype_ref.FUC12NpG = struct('type','sigmoidal_accumulation',...
    'Pmax', 91900, 'k', 0.07, 'lag', 4, 'n', 2);

%% =====================================================================
%%   REGULATORY SIGNALS (constitutive T7 promoter; scaled by pathway flux)
%% =====================================================================
% manA: native, activated by F6P availability (species 1)
cfg.regulatory.manA.default = struct('form','act',...
    'signals',{{1}}, 'K', 1.0, 'n', 2, 'basal', 0.20);
% manB: constitutive, scales with Man6P (species 2)
cfg.regulatory.manB.default = struct('form','act',...
    'signals',{{2}}, 'K', 0.3, 'n', 2, 'basal', 0.20);
% manC: constitutive, scales with Man1P and GTP (species 3, 12)
cfg.regulatory.manC.default = struct('form','act',...
    'signals',{{3, 12}}, 'K', [0.2, 1.0], 'n', [2, 2], 'basal', 0.15);
% gmd: constitutive, scales with GDP-Man (species 4)
cfg.regulatory.gmd.default = struct('form','act',...
    'signals',{{4}}, 'K', 0.1, 'n', 2, 'basal', 0.15);
% wcaG: constitutive, scales with GDP-4keto and NADPH (species 5, 14)
cfg.regulatory.wcaG.default = struct('form','act',...
    'signals',{{5, 14}}, 'K', [0.05, 0.1], 'n', [2, 2], 'basal', 0.15);
% wbgL: constitutive, dual substrate (GDP-Fuc + lactose) (species 6, 7)
cfg.regulatory.wbgL.default = struct('form','act',...
    'signals',{{6, 7}}, 'K', [0.05, 2.0], 'n', [2, 2], 'basal', 0.10);
% afcA: constitutive, activated by 2'-FL accumulation (species 8)
cfg.regulatory.afcA.default = struct('form','act',...
    'signals',{{8}}, 'K', 1.0, 'n', 2, 'basal', 0.10);
% ndk: constitutive where present, scales with GDP (species 13)
cfg.regulatory.ndk.default = struct('form','act',...
    'signals',{{13}}, 'K', 0.5, 'n', 2, 'basal', 0.05);
% Per-condition ndk override: absent in FUC1/FUC8
cfg.regulatory.ndk.FUC1 = struct('form','act',...
    'signals',{{13}}, 'K', 0.5, 'n', 2, 'basal', 0.0);
cfg.regulatory.ndk.FUC8 = struct('form','act',...
    'signals',{{13}}, 'K', 0.5, 'n', 2, 'basal', 0.0);

%% EXPERIMENTAL DATA (from Xia et al. 2025, mg/L at sampling time)
cfg.experimental_data.FUC1.fucose_titer     = 2250;   % shake-flask, 48h
cfg.experimental_data.FUC8.fucose_titer     = 6740;   % shake-flask, 48h
cfg.experimental_data.FUC12.fucose_titer    = 63770;  % bioreactor, 138h
cfg.experimental_data.FUC12NSM.fucose_titer = 80140;  % bioreactor, 114h
cfg.experimental_data.FUC12NpG.fucose_titer = 91900;  % bioreactor, 78h

%% =====================================================================
%%   HYPOTHESES (15 mechanistic hypotheses across 5 conditions)
%% =====================================================================
cfg.hypotheses.FUC1 = { ...
  struct('id','FUC1_basic','label','Base pathway, low copy, low flux',...
    'status','BASELINE', 'diagnostic_up',{{'gmd','wcaG','wbgL','afcA'}}, 'score_total',0.50), ...
  struct('id','FUC1_manC_bottleneck','label','ManC rate-limiting (GTP scarcity at 1x)',...
    'status','SUPPORTED', 'diagnostic_up',{{'manC','gmd'}}, 'score_total',0.55), ...
  struct('id','FUC1_growth_limited','label','Low biomass limits total capacity',...
    'status','SUPPORTED', 'diagnostic_up',{{'manA','manB'}}, 'score_total',0.45) };

cfg.hypotheses.FUC8 = { ...
  struct('id','FUC8_dosage_relief','label','4x gmd/wcaG relieves upstream bottleneck',...
    'status','BEST', 'diagnostic_up',{{'gmd','wcaG','wbgL'}}, 'score_total',0.75), ...
  struct('id','FUC8_afcA_unmatch','label','AfcA still 1x, 2FL accumulates as bottleneck',...
    'status','SUPPORTED', 'diagnostic_up',{{'afcA','wbgL'}}, 'score_total',0.70), ...
  struct('id','FUC8_wbgL_saturation','label','WbgL at 3x saturates GDP-Fuc channeling',...
    'status','SUPPORTED', 'diagnostic_up',{{'wbgL'}}, 'score_total',0.65) };

cfg.hypotheses.FUC12 = { ...
  struct('id','FUC12_balanced','label','4:3:2 ratio achieves optimal flux balance',...
    'status','BEST', 'diagnostic_up',{{'gmd','wcaG','wbgL','afcA'}}, 'score_total',0.85), ...
  struct('id','FUC12_afcA_critical','label','2x afcA is the key driver (2FL to FucL)',...
    'status','SUPPORTED', 'diagnostic_up',{{'afcA'}}, 'score_total',0.80), ...
  struct('id','FUC12_ndk_marginal','label','Ndk adds GTP recycling but only +4.5 pct',...
    'status','SUPPORTED', 'diagnostic_up',{{'ndk','manC'}}, 'score_total',0.70) };

cfg.hypotheses.FUC12NSM = { ...
  struct('id','FUC12NSM_retention','label','mdfA+setA KO prevents 2FL export, more afcA substrate',...
    'status','BEST', 'diagnostic_up',{{'afcA','wbgL'}}, 'score_total',0.80), ...
  struct('id','FUC12NSM_channeling','label','Intracellular 2FL retention creates afcA substrate pool',...
    'status','SUPPORTED', 'diagnostic_up',{{'afcA'}}, 'score_total',0.75), ...
  struct('id','FUC12NSM_growth_neutral','label','Transport KOs dont affect growth',...
    'status','SUPPORTED', 'diagnostic_up',{{'manA','gmd'}}, 'score_total',0.65) };

cfg.hypotheses.FUC12NpG = { ...
  struct('id','FUC12NpG_dual_carbon','label','Glycerol for growth, glucose for production',...
    'status','BEST', 'diagnostic_up',{{'manA','gmd','wcaG','wbgL','afcA'}}, 'score_total',0.85), ...
  struct('id','FUC12NpG_pfkA_redirect','label','pfkA KO redirects F6P from glycolysis to pathway',...
    'status','SUPPORTED', 'diagnostic_up',{{'manA','manC','gmd'}}, 'score_total',0.80), ...
  struct('id','FUC12NpG_growth_amplification','label','OD=205 massive enzyme pool, massive production',...
    'status','SUPPORTED', 'diagnostic_up',{{'gmd','wcaG','wbgL','afcA','ndk'}}, 'score_total',0.75) };

%% VISUALIZATION
cfg.visualization.colors.FUC1     = [0.267, 0.467, 0.667];  % steel blue
cfg.visualization.colors.FUC8     = [0.933, 0.533, 0.200];  % orange
cfg.visualization.colors.FUC12    = [0.467, 0.733, 0.267];  % green
cfg.visualization.colors.FUC12NSM = [0.667, 0.267, 0.600];  % purple
cfg.visualization.colors.FUC12NpG = [0.867, 0.200, 0.200];  % red
cfg.visualization.labels.FUC1     = 'FUC1 (1:1:1:1 baseline)';
cfg.visualization.labels.FUC8     = 'FUC8 (4:4:3:1 dosage)';
cfg.visualization.labels.FUC12    = 'FUC12 (4:4:3:2 +ndk)';
cfg.visualization.labels.FUC12NSM = 'FUC12NSM (mdfA/setA KO)';
cfg.visualization.labels.FUC12NpG = 'FUC12NpG (pfkA KO + glpK*)';

%% =====================================================================
%%   COUPLED EXPRESSION SETUP (extend 18-species ODE to 34-species)
%% =====================================================================
if cfg.coupled_expression
    expr_genes = cfg.genes;  % 8 genes
    n_g = numel(expr_genes);

    % Extended species names (mRNA_gene, Protein_gene)
    expr_sp = {};
    for gi = 1:n_g, expr_sp{end+1} = ['mRNA_' expr_genes{gi}]; end
    for gi = 1:n_g, expr_sp{end+1} = ['Protein_' expr_genes{gi}]; end

    % Expression parameters struct (packed into each model.params)
    ep = struct();
    ep.beta_m = cfg.expression.beta_m;
    ep.beta_p = cfg.expression.beta_p;
    ep.ktl    = cfg.expression.ktl;
    ep.n_genes = n_g;
    ep.ktx = zeros(n_g, 1);
    for gi = 1:n_g
        ep.ktx(gi) = cfg.ktx_fits.(expr_genes{gi});
    end

    % Process each condition
    cond_names = fieldnames(cfg.models);
    for ci = 1:numel(cond_names)
        cname = cond_names{ci};
        m = cfg.models.(cname);
        pp = m.params;

        % Copy numbers per gene (genes 1-3 native=1; 4-7 from params; 8=ndk)
        nc = ones(n_g, 1);
        if isfield(pp, 'n_gmd'),  nc(4) = pp.n_gmd;  end
        if isfield(pp, 'n_wcaG'), nc(5) = pp.n_wcaG; end
        if isfield(pp, 'n_wbgL'), nc(6) = pp.n_wbgL; end
        if isfield(pp, 'n_afcA'), nc(7) = pp.n_afcA; end
        nc(8) = double(pp.Vm_ndk > 0);  % ndk present only if Vm_ndk > 0

        % Get regulatory rules (condition-specific or default)
        rr = cell(n_g, 1);
        for gi = 1:n_g
            gene = expr_genes{gi};
            if isfield(cfg.regulatory, gene)
                if isfield(cfg.regulatory.(gene), cname)
                    rr{gi} = cfg.regulatory.(gene).(cname);
                else
                    rr{gi} = cfg.regulatory.(gene).default;
                end
            else
                rr{gi} = struct('form','const','basal',0.5);
            end
        end

        % Compute initial regulatory signals from metabolite x0
        x0_met = m.base_x0(:);
        u0 = zeros(n_g, 1);
        for gi = 1:n_g
            u0(gi) = eval_hill_reg(x0_met, rr{gi});
        end

        % Effective ktx = base_ktx * n_copy
        eff_ktx = ep.ktx(:) .* nc;

        % Steady-state initial conditions
        m_ss = eff_ktx .* u0 / ep.beta_m;       % mRNA
        p_ss = ep.ktl * m_ss / ep.beta_p;        % Protein

        % Extend model
        m.species  = [species_list, expr_sp];
        m.base_x0  = [x0_met; m_ss; p_ss];

        % Pack expression params into model params
        pp.expr = ep;
        pp.expr.n_copy = nc;
        pp.expr.reg_rules = rr;
        m.params = pp;
        m.rate_fn = @ecoli_mm_rates_coupled;

        cfg.models.(cname) = m;
    end
end

%% OUTPUT DIRECTORY
cfg.output_dir = '';
end

%% =====================================================================
%%   RATE FUNCTION: Michaelis-Menten with cofactors and growth
%% =====================================================================
function v = ecoli_mm_rates(x, p)
%ECOLI_MM_RATES Compute 20 reaction rates for E. coli fucose pathway.
%   v = ecoli_mm_rates(x, p)
%   x(18x1) = species concentrations (mM)
%   p = parameter struct with Vmax, Km, copy numbers, growth params
%
%   All enzymatic and transport reactions are BIOMASS-SCALED:
%   v = Vmax * B * S/(Km+S)  where B = biomass (OD proxy).
%   This captures the key biological reality that more cells = more enzyme = more flux.
%   Vmax values are per-OD-unit specific activities.

    x = max(x, 0);
    v = zeros(20, 1);

    % Unpack species for readability
    F6P       = x(1);
    man6p     = x(2);
    man1p     = x(3);
    gdpman    = x(4);
    gdp4keto  = x(5);
    gdpfuc    = x(6);
    lcts_c    = x(7);
    fl2_c     = x(8);
    fuc_L_c   = x(9);
    % fuc_L_e = x(10);  % extracellular product (accumulates only)
    % fl2_e   = x(11);  % extracellular 2'-FL (accumulates only)
    GTP       = x(12);
    GDP       = x(13);
    NADPH     = x(14);
    biomass   = x(15);
    glc_ext   = x(16);
    gly_ext   = x(17);
    lcts_ext  = x(18);

    % Biomass factor: all enzymatic rates scale with cell density
    B = max(biomass, 0.01);

    % Expression/efficiency factor: scales all pathway Vmax per condition
    EF = p.expr_factor;

    % Expression/efficiency factor: scales all pathway Vmax per condition

    % R1: manA -- F6P -> Man6P (biomass-scaled MM)
    v(1) = p.Vm_manA * p.n_manA * EF * B * F6P / (p.Km_manA + F6P);

    % R2: manB -- Man6P -> Man1P (biomass-scaled)
    v(2) = p.Vm_manB * p.n_manB * EF * B * man6p / (p.Km_manB + man6p);

    % R3: manC -- Man1P + GTP -> GDPMan + GDP (bi-substrate, biomass-scaled)
    v(3) = p.Vm_manC * p.n_manC * EF * B * man1p * GTP / ((p.Km_manC_m1p + man1p) * (p.Km_manC_gtp + GTP));

    % R4: gmd -- GDPMan -> GDP4keto (copy-number + biomass scaled)
    v(4) = p.Vm_gmd * p.n_gmd * EF * B * gdpman / (p.Km_gmd + gdpman);

    % R5: wcaG -- GDP4keto + NADPH -> GDPFuc (bi-substrate, copy + biomass scaled)
    v(5) = p.Vm_wcaG * p.n_wcaG * EF * B * gdp4keto * NADPH / ((p.Km_wcaG_sub + gdp4keto) * (p.Km_wcaG_nadph + NADPH));

    % R6: wbgL -- GDPFuc + Lcts -> 2FL + GDP (bi-substrate, copy + biomass scaled)
    v(6) = p.Vm_wbgL * p.n_wbgL * EF * B * gdpfuc * lcts_c / ((p.Km_wbgL_gfuc + gdpfuc) * (p.Km_wbgL_lcts + lcts_c));

    % R7: afcA -- 2FL -> FucL + Lcts_recycle (copy + biomass scaled)
    v(7) = p.Vm_afcA * p.n_afcA * EF * B * fl2_c / (p.Km_afcA + fl2_c);

    % R8: fuc_export -- FucL_c -> FucL_e (biomass-scaled transport)
    v(8) = p.k_fuc_export * B * fuc_L_c;

    % R9: mdfA -- 2FL_c -> 2FL_e (biomass-scaled; k=0 when knocked out)
    v(9) = p.k_mdfA * B * fl2_c;

    % R10: setA -- 2FL_c -> 2FL_e (biomass-scaled; k=0 when knocked out)
    v(10) = p.k_setA * B * fl2_c;

    % R11: ndk -- GDP -> GTP (biomass-scaled; Vm=0 when absent)
    v(11) = p.Vm_ndk * B * GDP / (p.Km_ndk + GDP);

    % R12: growth -- Monod x logistic (glucose-dependent)
    if biomass < p.Bmax && glc_ext > 0
        v(12) = p.mu_max * biomass * glc_ext / (p.Ks_glc + glc_ext) * (1 - biomass / p.Bmax);
    else
        v(12) = 0;
    end

    % R13: pfkA -- F6P -> drain via glycolysis (biomass-scaled; Vm=0 when knocked out)
    v(13) = p.Vm_pfkA * B * F6P / (p.Km_pfkA + F6P);

    % R14: glpK -- glycerol -> F6P (biomass-scaled; Vm=0 when absent)
    v(14) = p.Vm_glpK * B * gly_ext / (p.Km_glpK + gly_ext);

    % R15: glc_to_F6P -- glucose -> F6P (growth-coupled supply, scaled by pathway demand)
    v(15) = p.k_glc_F6P * p.pathway_pull * glc_ext * biomass / (p.Ks_glc + glc_ext);

    % R16: lcts_import -- lcts_ext -> lcts_c (first-order)
    v(16) = p.k_lcts_import * lcts_ext;

    % R17: NADPH_regen -- biosynthetic NADPH regeneration (growth-coupled)
    v(17) = p.k_NADPH_regen * biomass;

    % R18: GTP_base -- basal GTP biosynthesis (growth-coupled)
    v(18) = p.k_GTP_base * biomass;

    % R19: glc_feed -- constant glucose feeding
    v(19) = p.k_glc_feed;

    % R20: gly_feed -- constant glycerol feeding (0 for non-NpG)
    v(20) = p.k_gly_feed;
end

%% =====================================================================
%%   COUPLED RATE FUNCTION: metabolic velocities + expression dynamics
%% =====================================================================
function out = ecoli_mm_rates_coupled(x, p)
%ECOLI_MM_RATES_COUPLED Metabolic velocities (identical to legacy) plus
%   mRNA/protein expression dynamics as synchronized passengers.
%   Returns struct with .v (20x1) and .dxdt_expr (16x1).

    x = max(x, 0);

    % -- Metabolic velocities (EXACT replica of ecoli_mm_rates) ----------
    v = zeros(20, 1);
    F6P      = x(1);  man6p    = x(2);  man1p    = x(3);  gdpman   = x(4);
    gdp4keto = x(5);  gdpfuc   = x(6);  lcts_c   = x(7);  fl2_c    = x(8);
    fuc_L_c  = x(9);  GTP      = x(12); GDP      = x(13); NADPH    = x(14);
    biomass  = x(15); glc_ext  = x(16); gly_ext  = x(17); lcts_ext = x(18);

    v(1)  = p.Vm_manA * F6P / (p.Km_manA + F6P);
    v(2)  = p.Vm_manB * man6p / (p.Km_manB + man6p);
    v(3)  = p.Vm_manC * man1p * GTP / ((p.Km_manC_m1p + man1p) * (p.Km_manC_gtp + GTP));
    v(4)  = p.Vm_gmd * p.n_gmd * gdpman / (p.Km_gmd + gdpman);
    v(5)  = p.Vm_wcaG * p.n_wcaG * gdp4keto * NADPH / ((p.Km_wcaG_sub + gdp4keto) * (p.Km_wcaG_nadph + NADPH));
    v(6)  = p.Vm_wbgL * p.n_wbgL * gdpfuc * lcts_c / ((p.Km_wbgL_gfuc + gdpfuc) * (p.Km_wbgL_lcts + lcts_c));
    v(7)  = p.Vm_afcA * p.n_afcA * fl2_c / (p.Km_afcA + fl2_c);
    v(8)  = p.k_fuc_export * fuc_L_c;
    v(9)  = p.k_mdfA * fl2_c;
    v(10) = p.k_setA * fl2_c;
    v(11) = p.Vm_ndk * GDP / (p.Km_ndk + GDP);
    if biomass < p.Bmax && glc_ext > 0
        v(12) = p.mu_max * biomass * glc_ext / (p.Ks_glc + glc_ext) * (1 - biomass / p.Bmax);
    else
        v(12) = 0;
    end
    v(13) = p.Vm_pfkA * F6P / (p.Km_pfkA + F6P);
    v(14) = p.Vm_glpK * gly_ext / (p.Km_glpK + gly_ext);
    v(15) = p.k_glc_F6P * glc_ext * biomass / (p.Ks_glc + glc_ext);
    v(16) = p.k_lcts_import * lcts_ext;
    v(17) = p.k_NADPH_regen * biomass;
    v(18) = p.k_GTP_base * biomass;
    v(19) = p.k_glc_feed;
    v(20) = p.k_gly_feed;

    % -- Expression dynamics (synchronized passengers) -------------------
    n_g  = p.expr.n_genes;        % 8
    m_rna = x(19:18+n_g);         % mRNA(1:8)
    prot  = x(19+n_g:18+2*n_g);   % Protein(1:8)

    % Regulatory signals from metabolites
    x_met = x(1:18);
    u_reg = zeros(n_g, 1);
    for gi = 1:n_g
        u_reg(gi) = eval_hill_reg(x_met, p.expr.reg_rules{gi});
    end

    % Effective transcription rate = base_ktx * copy_number
    eff_ktx = p.expr.ktx(:) .* p.expr.n_copy(:);

    % dm/dt = ktx_eff * u - beta_m * m
    dmdt = eff_ktx .* u_reg - p.expr.beta_m * m_rna;

    % dp/dt = ktl * m - beta_p * p
    dpdt = p.expr.ktl * m_rna - p.expr.beta_p * prot;

    out.v = v;
    out.dxdt_expr = [dmdt; dpdt];
end

%% =====================================================================
%%   HELPER: Evaluate Hill regulatory signal from metabolite state
%% =====================================================================
function u = eval_hill_reg(x_met, rule)
%EVAL_HILL_REG Compute scalar regulatory signal u in [basal, 1].
    if isfield(rule, 'form') && strcmp(rule.form, 'const')
        u = rule.basal;
        return;
    end

    hill_prod = 1.0;
    for si = 1:numel(rule.signals)
        sig_spec = rule.signals{si};
        if iscell(sig_spec)
            s_val = sum(x_met([sig_spec{:}]));
        else
            s_val = x_met(sig_spec);
        end
        Ki = rule.K(min(si, numel(rule.K)));
        ni = rule.n(min(si, numel(rule.n)));
        hill_prod = hill_prod * s_val^ni / (Ki^ni + s_val^ni);
    end

    u = max(rule.basal, hill_prod);
end