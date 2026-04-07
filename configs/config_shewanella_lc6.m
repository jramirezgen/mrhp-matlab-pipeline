function cfg = config_shewanella_lc6()
%CONFIG_SHEWANELLA_LC6 Configuration for Shewanella xiamenensis LC6 MRHP.
%   Reproduces the canonical 3-dye, 6-gene, 36-scenario pipeline.

cfg.organism = 'Shewanella xiamenensis LC6';

%% CONDITIONS
cfg.conditions = {'MO', 'RC', 'AD'};

%% GENES
cfg.genes = {'AzoR','AceA','bepE','mtrF','CymA','YhhW'};

%% SOLVER
cfg.solver.t_span = [0, 72];
cfg.solver.dt_eval = 0.5;
cfg.solver.t_eval = (0:0.5:72)';
cfg.solver.rtol = 1e-8;
cfg.solver.atol = 1e-10;
cfg.solver.max_step = 0.5;
cfg.solver.substrate_conversion = 3.7; % YE g/L to mM

%% SWEEP (YE x Dye concentration grid)
cfg.sweep.substrate_values = [0.5, 1.0, 3.0, 5.0];
cfg.sweep.target_values    = [0.4, 0.6, 1.0];

%% REFERENCE CONDITION
cfg.ref_substrate = 1.0;
cfg.ref_target    = 0.6;
cfg.t_sample      = 24.0;

%% EXPRESSION LAYER
cfg.expression.beta_m  = 8.32;
cfg.expression.beta_p  = 0.924;
cfg.expression.ktl     = 6.0;
cfg.expression.stochastic_tau = 0.01;
cfg.expression.n_cells = 20;
cfg.expression.nm_to_mol = 6022;

%% TRANSCRIPTION RATES (FROZEN)
cfg.ktx_fits.AzoR = 1.1597;
cfg.ktx_fits.AceA = 0.0218;
cfg.ktx_fits.bepE = 0.1553;
cfg.ktx_fits.mtrF = 0.0427;
cfg.ktx_fits.CymA = 0.0063;
cfg.ktx_fits.YhhW = 0.0104;

%% COMMON KINETIC PARAMS
cp.Vmax_ye   = 2.0;  cp.Km_ye   = 5.0;
cp.Vmax_pdh  = 5.0;  cp.Km_pyr  = 0.5;
cp.Vmax_tca  = 3.0;  cp.Km_accoa = 0.1;
cp.k_resp    = 15.0; cp.Km_nad  = 0.3; cp.Km_nadh = 0.15;
cp.k_frag    = 0.4;  cp.k_hub   = 0.3;
cp.k_reint   = 0.5;  cp.k_sucrecycle = 0.3; cp.k_efflux = 0.6;
cp.mu_max    = 0.20; cp.Ks_growth = 3.0; cp.Bmax = 1.5; cp.Y_xs = 8.0;

%% ═══ MO MODEL ═══
m_mo.species = {'YE_eff','Pyr','AcCoA','NADH','NAD','MO','DMPD',...
                'Sulfanilate','Catechol','OxoAdipate','CO2','SucCoA','Biomass'};
m_mo.S = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1;
            2,-1, 0, 0, 0, 0, 0, 0, 0, 1, 0;
            0, 1,-1, 0, 0, 0, 0, 0, 1, 0, 0;
            2, 1, 3,-1,-2, 0, 0, 0, 0, 0, 0;
           -2,-1,-3, 1, 2, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 1, 1,-1, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0;
            0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0, 0, 0, 1,-1, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
m_mo.base_x0 = [0, 0.1, 0.05, 0.3, 2.7, 0, 0, 0, 0, 0, 0, 0.01, 0.10]';
m_mo.substrate_idx = 1;
m_mo.target_idx = 6;
m_mo.biomass_idx = 13;
p_mo = cp;
p_mo.Vmax_azo = 0.10; p_mo.Km_dye = 1.0;
m_mo.params = p_mo;
m_mo.rate_fn = @mo_rates;
cfg.models.MO = m_mo;

%% ═══ RC MODEL ═══
m_rc.species = {'YE_eff','Pyr','AcCoA','NADH','NAD','RC','DANS','Benzidine',...
                'Benzidine_ext','DHN','Salicylate','Catechol','OxoAdipate',...
                'CO2','SucCoA','Biomass'};
m_rc.S = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1;
            2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
            0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
            2, 1, 3,-1,-4, 0, 0, 0, 0, 0, 0, 0, 0;
           -2,-1,-3, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 2, 0,-1, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0;
            0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
m_rc.base_x0 = [0, 0.1, 0.05, 0.3, 2.7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.10]';
m_rc.substrate_idx = 1;
m_rc.target_idx = 6;
m_rc.biomass_idx = 16;
p_rc = cp;
p_rc.Vmax_azo = 0.12; p_rc.Km_dye = 2.21;
m_rc.params = p_rc;
m_rc.rate_fn = @rc_rates;
cfg.models.RC = m_rc;

%% ═══ AD MODEL ═══
m_ad.species = {'YE_eff','Pyr','AcCoA','NADH','NAD','AD','DAN','J_acid',...
                'Bisulfonate','Bisulfonate_ext','DHN','Salicylate','Catechol',...
                'Protocatechuate','OxoAdipate','CO2','SucCoA','Biomass'};
m_ad.S = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1;
            2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
            0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
            2, 1, 3,-1,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           -2,-1,-3, 1, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0;
            0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,-1, 0, 0;
            0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
m_ad.base_x0 = [0, 0.1, 0.05, 0.3, 2.7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.10]';
m_ad.substrate_idx = 1;
m_ad.target_idx = 6;
m_ad.biomass_idx = 18;
p_ad = cp;
p_ad.Vmax_azo = 0.08; p_ad.Km_dye = 1.0;
m_ad.params = p_ad;
m_ad.rate_fn = @ad_rates;
cfg.models.AD = m_ad;

%% PHI DEFINITION (depletion for all dyes)
cfg.phi_definition.default = struct('type', 'depletion', 'target_idx', 6);

%% PHENOTYPE REFERENCE (V4 Hill)
cfg.phenotype_label = 'Decolorization (%)';
V4.K_yc  = 5.23111;
V4.K_dic = 0.141446;
cfg.phenotype_ref.MO = struct('type','hill_decolorization', ...
    'Dmax',95.88, 'Vmax',3.29, 'n',1.81, 'lag',1.65, ...
    'K_yc',V4.K_yc, 'K_dic',V4.K_dic, ...
    'default_substrate', cfg.ref_substrate, 'default_target', cfg.ref_target);
cfg.phenotype_ref.RC = struct('type','hill_decolorization', ...
    'Dmax',92.26, 'Vmax',2.93, 'n',2.83, 'lag',6.40, ...
    'K_yc',V4.K_yc, 'K_dic',V4.K_dic, ...
    'default_substrate', cfg.ref_substrate, 'default_target', cfg.ref_target);
cfg.phenotype_ref.AD = struct('type','hill_decolorization', ...
    'Dmax',94.63, 'Vmax',0.74, 'n',0.84, 'lag',2.30, ...
    'K_yc',V4.K_yc, 'K_dic',V4.K_dic, ...
    'default_substrate', cfg.ref_substrate, 'default_target', cfg.ref_target);

%% REGULATORY SIGNALS
% AzoR: activation by Dye(6) and NADH(4)
cfg.regulatory.AzoR.default = struct('form','act', ...
    'signals',{{6, 4}}, 'K',[0.3, 0.2], 'n',[2, 2], 'basal',0.01);
% AceA: repression by AcCoA(3) and YE(1)
cfg.regulatory.AceA.default = struct('form','rep', ...
    'signals',{{3, 1}}, 'K',[0.5, 2.0], 'n',[2, 2], 'basal',0.05);
% bepE: activation by toxic (catechol + efflux substrate + oxoadipate)
cfg.regulatory.bepE.MO = struct('form','act', ...
    'signals',{{[9, 10]}}, 'K',0.05, 'n',2, 'basal',0.02);
cfg.regulatory.bepE.RC = struct('form','act', ...
    'signals',{{[12, 8, 13]}}, 'K',0.05, 'n',2, 'basal',0.02);
cfg.regulatory.bepE.AD = struct('form','act', ...
    'signals',{{[13, 9, 15]}}, 'K',0.05, 'n',2, 'basal',0.02);
% mtrF: activation by redox NADH/(NADH+NAD)
cfg.regulatory.mtrF.default = struct('form','custom', 'basal',0.10, ...
    'fn', @(y,d,r) max(r.basal, (y(4,:)./(y(4,:)+y(5,:)+1e-12)).^2 ./ (0.15^2 + (y(4,:)./(y(4,:)+y(5,:)+1e-12)).^2)));
% CymA: activation + constitutive by redox
cfg.regulatory.CymA.default = struct('form','custom', 'basal',0.30, ...
    'fn', @(y,d,r) max(r.basal, 0.5 * ((y(4,:)./(y(4,:)+y(5,:)+1e-12)).^2 ./ (0.15^2 + (y(4,:)./(y(4,:)+y(5,:)+1e-12)).^2)) + 0.5));
% YhhW: activation by catechol
cfg.regulatory.YhhW.MO = struct('form','act', ...
    'signals',{{9}}, 'K',0.02, 'n',2, 'basal',0.05);
cfg.regulatory.YhhW.RC = struct('form','act', ...
    'signals',{{12}}, 'K',0.02, 'n',2, 'basal',0.05);
cfg.regulatory.YhhW.AD = struct('form','act', ...
    'signals',{{13}}, 'K',0.02, 'n',2, 'basal',0.05);

%% EXPERIMENTAL DATA (RT-qPCR nM means)
cfg.experimental_data.MO.AzoR = 7.730739e-02;
cfg.experimental_data.MO.CymA = 1.744733e-04;
cfg.experimental_data.MO.YhhW = 2.383982e-04;
cfg.experimental_data.RC.AzoR = 6.942640e-02;
cfg.experimental_data.RC.CymA = 1.584105e-03;
cfg.experimental_data.RC.YhhW = 7.919690e-04;
cfg.experimental_data.AD.AzoR = 2.395186e-01;
cfg.experimental_data.AD.CymA = 1.290547e-03;
cfg.experimental_data.AD.YhhW = 7.609609e-04;

%% HYPOTHESES (FROZEN)
cfg.hypotheses.MO = { ...
  struct('id','MO_meta_dual', 'label','Meta-cleavage dual', 'status','BEST', ...
    'diagnostic_up',{{'AzoR','YhhW','bepE'}}, 'score_total',0.8054), ...
  struct('id','MO_meta_sulfa', 'label','Meta-cleavage sulfanilate', 'status','SUPPORTED', ...
    'diagnostic_up',{{'AzoR','YhhW'}}, 'score_total',0.8221), ...
  struct('id','MO_dual_efflux', 'label','Dual meta-cleavage + efflux', 'status','SUPPORTED', ...
    'diagnostic_up',{{'AzoR','YhhW','bepE','mtrF'}}, 'score_total',0.7971) };
cfg.hypotheses.RC = { ...
  struct('id','RC_Mtr_EET', 'label','Mtr-dependent EET', 'status','BEST', ...
    'diagnostic_up',{{'AzoR','mtrF','CymA'}}, 'score_total',0.8845), ...
  struct('id','RC_meta_efflux', 'label','Meta + efflux benzidine', 'status','SUPPORTED', ...
    'diagnostic_up',{{'AzoR','bepE','YhhW'}}, 'score_total',0.9179), ...
  struct('id','RC_gentisate', 'label','Gentisate (impaired)', 'status','SUPPORTED*', ...
    'diagnostic_up',{{'AzoR'}}, 'score_total',0.7845) };
cfg.hypotheses.AD = { ...
  struct('id','DB_maxefflux', 'label','Max efflux + meta + EET', 'status','BEST', ...
    'diagnostic_up',{{'AzoR','bepE','mtrF','CymA','YhhW'}}, 'score_total',0.9091), ...
  struct('id','DB_Jacid_meta', 'label','J-acid meta-cleavage', 'status','SUPPORTED', ...
    'diagnostic_up',{{'AzoR','YhhW'}}, 'score_total',0.8424), ...
  struct('id','DB_meta_simple', 'label','Simple meta-cleavage', 'status','SUPPORTED', ...
    'diagnostic_up',{{'AzoR','YhhW'}}, 'score_total',0.9091) };

%% VISUALIZATION
cfg.visualization.colors.MO = [0.902, 0.224, 0.275];
cfg.visualization.colors.RC = [0.271, 0.482, 0.616];
cfg.visualization.colors.AD = [0.165, 0.616, 0.561];
cfg.visualization.labels.MO = 'Methyl Orange';
cfg.visualization.labels.RC = 'Congo Red';
cfg.visualization.labels.AD = 'Direct Blue 71';

%% ROUTE MAPS (GEM→ODE traceability)
cfg.route_maps.MO = struct( ...
    'route_map',   'route_maps/route_map_MO.tsv', ...
    'species_map', 'route_maps/species_map_MO.tsv', ...
    'S_matrix',    'route_maps/S_matrix_MO.tsv');
cfg.route_maps.RC = struct( ...
    'route_map',   'route_maps/route_map_RC.tsv', ...
    'species_map', 'route_maps/species_map_RC.tsv', ...
    'S_matrix',    'route_maps/S_matrix_RC.tsv');
cfg.route_maps.AD = struct( ...
    'route_map',   'route_maps/route_map_AD.tsv', ...
    'species_map', 'route_maps/species_map_AD.tsv', ...
    'S_matrix',    'route_maps/S_matrix_AD.tsv');

%% OUTPUT DIRECTORY (will be set by runner)
cfg.output_dir = '';
end

%% ═══ RATE FUNCTIONS ═══
function v = mo_rates(x, p)
    x = max(x, 0);
    YE=x(1); Pyr=x(2); AcCoA=x(3); NADH=x(4); NAD=x(5);
    MO=x(6); DMPD=x(7); Sulfa=x(8); Cat=x(9); OxoAd=x(10);
    SucCoA=x(12); B=max(x(13), 1e-12);
    nad_lim = NAD / (p.Km_nad + NAD);
    v = zeros(11,1);
    v(1)  = p.Vmax_ye  * YE    / (p.Km_ye   + YE)    * nad_lim * B;
    v(2)  = p.Vmax_pdh * Pyr   / (p.Km_pyr  + Pyr)   * nad_lim * B;
    v(3)  = p.Vmax_tca * AcCoA / (p.Km_accoa + AcCoA) * nad_lim * B;
    v(4)  = p.k_resp   * NADH  * B;
    v(5)  = p.Vmax_azo * MO    / (p.Km_dye  + MO)    * NADH / (p.Km_nadh + NADH) * B;
    v(6)  = p.k_frag   * DMPD  * B;
    v(7)  = p.k_frag   * Sulfa * B;
    v(8)  = p.k_hub    * Cat   * B;
    v(9)  = p.k_reint  * OxoAd * B;
    v(10) = p.k_sucrecycle * SucCoA * B;
    v(11) = p.mu_max * B * YE / (p.Ks_growth + YE) * max(0, 1 - B / p.Bmax);
end

function v = rc_rates(x, p)
    x = max(x, 0);
    YE=x(1); Pyr=x(2); AcCoA=x(3); NADH=x(4); NAD=x(5);
    RC=x(6); DANS=x(7); Benz=x(8);
    DHN=x(10); Sal=x(11); Cat=x(12); OxoAd=x(13);
    SucCoA=x(15); B=max(x(16), 1e-12);
    nad_lim = NAD / (p.Km_nad + NAD);
    v = zeros(13,1);
    v(1)  = p.Vmax_ye  * YE    / (p.Km_ye   + YE)    * nad_lim * B;
    v(2)  = p.Vmax_pdh * Pyr   / (p.Km_pyr  + Pyr)   * nad_lim * B;
    v(3)  = p.Vmax_tca * AcCoA / (p.Km_accoa + AcCoA) * nad_lim * B;
    v(4)  = p.k_resp   * NADH  * B;
    v(5)  = p.Vmax_azo * RC    / (p.Km_dye  + RC)    * NADH / (p.Km_nadh + NADH) * B;
    v(6)  = p.k_efflux * Benz  * B;
    v(7)  = p.k_frag   * DANS  * B;
    v(8)  = p.k_frag   * DHN   * B;
    v(9)  = p.k_frag   * Sal   * B;
    v(10) = p.k_hub    * Cat   * B;
    v(11) = p.k_reint  * OxoAd * B;
    v(12) = p.k_sucrecycle * SucCoA * B;
    v(13) = p.mu_max * B * YE / (p.Ks_growth + YE) * max(0, 1 - B / p.Bmax);
end

function v = ad_rates(x, p)
    x = max(x, 0);
    YE=x(1); Pyr=x(2); AcCoA=x(3); NADH=x(4); NAD=x(5);
    AD_c=x(6); DAN=x(7); Jac=x(8); Bisulf=x(9);
    DHN=x(11); Sal=x(12); Cat=x(13); Pcat=x(14); OxoAd=x(15);
    SucCoA=x(17); B=max(x(18), 1e-12);
    nad_lim = NAD / (p.Km_nad + NAD);
    v = zeros(15,1);
    v(1)  = p.Vmax_ye  * YE    / (p.Km_ye   + YE)    * nad_lim * B;
    v(2)  = p.Vmax_pdh * Pyr   / (p.Km_pyr  + Pyr)   * nad_lim * B;
    v(3)  = p.Vmax_tca * AcCoA / (p.Km_accoa + AcCoA) * nad_lim * B;
    v(4)  = p.k_resp   * NADH  * B;
    v(5)  = p.Vmax_azo * AD_c  / (p.Km_dye  + AD_c)  * NADH / (p.Km_nadh + NADH) * B;
    v(6)  = p.k_efflux * Bisulf * B;
    v(7)  = p.k_frag   * DAN   * B;
    v(8)  = p.k_frag   * DHN   * B;
    v(9)  = p.k_frag   * Sal   * B;
    v(10) = p.k_hub    * Cat   * B;
    v(11) = p.k_frag   * Jac   * B;
    v(12) = p.k_hub    * Pcat  * B;
    v(13) = p.k_reint  * OxoAd * B;
    v(14) = p.k_sucrecycle * SucCoA * B;
    v(15) = p.mu_max * B * YE / (p.Ks_growth + YE) * max(0, 1 - B / p.Bmax);
end