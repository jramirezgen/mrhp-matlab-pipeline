% MRHP Universal Platform — Master Runner
% Runs the generic pipeline for one or all organism configurations.
%
% Usage:
%   cd mrhp-matlab-pipeline
%   matlab -batch "run_all"           % runs all 3 organisms
%
% Or select a specific organism:
%   MRHP_ORGANISM=shewanella matlab -batch "run_all"
%   MRHP_ORGANISM=ecoli      matlab -batch "run_all"
%   MRHP_ORGANISM=acido      matlab -batch "run_all"

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf(  '║     MRHP Universal Platform — Master Runner      ║\n');
fprintf(  '╚══════════════════════════════════════════════════╝\n\n');

%% Setup paths
base_dir   = fileparts(mfilename('fullpath'));  % repo root
engine_dir = fullfile(base_dir, 'engine');
config_dir = fullfile(base_dir, 'configs');
addpath(engine_dir);
addpath(config_dir);

%% Determine which organisms to run
organism_env = strtrim(getenv('MRHP_ORGANISM'));
if isempty(organism_env)
    organism_env = 'all';
end

configs = {};
switch lower(organism_env)
    case 'shewanella'
        configs = {{'shewanella', @config_shewanella_lc6}};
    case 'ecoli'
        configs = {{'ecoli_fucose', @config_ecoli_fucose}};
    case 'acido'
        configs = {{'acidithiobacillus', @config_acidithiobacillus_fe}};
    case 'all'
        configs = { ...
            {'shewanella',       @config_shewanella_lc6}, ...
            {'ecoli_fucose',     @config_ecoli_fucose}, ...
            {'acidithiobacillus', @config_acidithiobacillus_fe} };
    otherwise
        error('Unknown organism: %s. Use shewanella, ecoli, acido, or all.', organism_env);
end

%% Quick or full mode
mode_env = strtrim(getenv('MRHP_MODE'));
if isempty(mode_env)
    mode_env = 'full';
end

%% Run each organism
total_t0 = tic;
n_org = numel(configs);
results_all = struct();

for oi = 1:n_org
    org_name   = configs{oi}{1};
    config_fn  = configs{oi}{2};

    fprintf('\n━━━ [%d/%d] %s ━━━\n', oi, n_org, upper(org_name));

    % Load configuration
    cfg = config_fn();

    % Set output directory
    cfg.output_dir = fullfile(base_dir, 'outputs', org_name);
    if ~exist(cfg.output_dir, 'dir')
        mkdir(cfg.output_dir);
    end

    % Run the generic pipeline
    try
        if strcmpi(mode_env, 'quick')
            run_pipeline_generic(cfg, 'quick');
        else
            run_pipeline_generic(cfg);
        end
        results_all.(org_name) = 'PASS';
        fprintf('  ✓ %s — PASS\n', org_name);
    catch ME
        results_all.(org_name) = sprintf('FAIL: %s', ME.message);
        fprintf('  ✗ %s — FAIL: %s\n', org_name, ME.message);
        fprintf('    Stack: %s line %d\n', ME.stack(1).name, ME.stack(1).line);
    end
end

elapsed = toc(total_t0);

%% Summary
fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf(  '║              MRHP UNIVERSAL SUMMARY              ║\n');
fprintf(  '╠══════════════════════════════════════════════════╣\n');
orgs = fieldnames(results_all);
all_pass = true;
for k = 1:numel(orgs)
    status = results_all.(orgs{k});
    if startsWith(status, 'FAIL')
        all_pass = false;
    end
    fprintf('║  %-20s  %s\n', orgs{k}, status);
end
fprintf('╠══════════════════════════════════════════════════╣\n');
fprintf('║  Total elapsed: %.1f s\n', elapsed);
if all_pass
    fprintf('║  VERDICT: ALL ORGANISMS CERTIFIED ✓\n');
else
    fprintf('║  VERDICT: SOME ORGANISMS FAILED ✗\n');
end
fprintf('╚══════════════════════════════════════════════════╝\n');

%% Write master metadata
meta_dir = fullfile(base_dir, 'metadata');
if ~exist(meta_dir, 'dir'), mkdir(meta_dir); end
meta_path = fullfile(meta_dir, 'master_run_info.json');
fid = fopen(meta_path, 'w');
if fid > 0
    fprintf(fid, '{\n');
    fprintf(fid, '  "timestamp": "%s",\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '  "mode": "%s",\n', mode_env);
    fprintf(fid, '  "elapsed_seconds": %.1f,\n', elapsed);
    fprintf(fid, '  "organisms": {\n');
    for k = 1:numel(orgs)
        comma = ',';
        if k == numel(orgs), comma = ''; end
        fprintf(fid, '    "%s": "%s"%s\n', orgs{k}, results_all.(orgs{k}), comma);
    end
    fprintf(fid, '  },\n');
    fprintf(fid, '  "all_pass": %s\n', mat2str(all_pass));
    fprintf(fid, '}\n');
    fclose(fid);
    fprintf('\nMetadata → %s\n', meta_path);
end
