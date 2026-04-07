%RUN_ALL_SYSTEMS Master runner for all organism pipelines.
%   Executes Shewanella LC6, E. coli fucose, and A. ferrooxidans
%   through the MRHP universal pipeline engine.

function run_all_systems(mode)
    if nargin < 1, mode = 'quick'; end

    t0 = tic;
    base = fileparts(mfilename('fullpath'));
    addpath(fullfile(base, 'engine'));
    addpath(fullfile(base, 'configs'));

    organisms = {
        'config_shewanella_lc6',       'Shewanella LC6';
        'config_ecoli_fucose',         'E. coli fucose';
        'config_acidithiobacillus_fe', 'A. ferrooxidans'
    };

    results_summary = {};
    n_org = size(organisms, 1);

    fprintf('\n========================================\n');
    fprintf('  MRHP UNIVERSAL PLATFORM — ALL SYSTEMS\n');
    fprintf('  Mode: %s\n', mode);
    fprintf('  Organisms: %d\n', n_org);
    fprintf('========================================\n\n');

    for oi = 1:n_org
        config_fn = organisms{oi, 1};
        org_label = organisms{oi, 2};

        fprintf('\n────────────────────────────────────────\n');
        fprintf('  [%d/%d] %s\n', oi, n_org, org_label);
        fprintf('────────────────────────────────────────\n\n');

        try
            cfg = feval(config_fn);

            % Set output directory
            out_dir = fullfile(base, 'outputs', config_fn);
            if ~exist(out_dir, 'dir'), mkdir(out_dir); end
            cfg.output_dir = out_dir;

            % Run pipeline
            run_pipeline_generic(cfg, mode);

            % Read certification
            cert_file = fullfile(out_dir, '..', 'metadata', 'certification.json');
            if exist(cert_file, 'file')
                cert = jsondecode(fileread(cert_file));
                status = cert.status;
            else
                status = 'COMPLETED';
            end

            results_summary{end+1} = struct('organism', org_label, ...
                'config', config_fn, 'status', status); %#ok
            fprintf('\n  >>> %s: %s <<<\n', org_label, status);

        catch ME
            fprintf('\n  !!! %s FAILED: %s !!!\n', org_label, ME.message);
            fprintf('  Stack trace:\n');
            for si = 1:numel(ME.stack)
                fprintf('    %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
            end
            results_summary{end+1} = struct('organism', org_label, ...
                'config', config_fn, 'status', 'FAILED'); %#ok
        end
    end

    %% CROSS-SYSTEM SUMMARY
    fprintf('\n========================================\n');
    fprintf('  CROSS-SYSTEM SUMMARY\n');
    fprintf('========================================\n\n');
    for i = 1:numel(results_summary)
        r = results_summary{i};
        fprintf('  %-25s  %s\n', r.organism, r.status);
    end
    fprintf('\n  Total time: %.1fs\n', toc(t0));
    fprintf('========================================\n');
end
