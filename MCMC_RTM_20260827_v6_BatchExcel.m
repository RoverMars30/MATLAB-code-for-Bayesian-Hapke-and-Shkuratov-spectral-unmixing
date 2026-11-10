%% MCMC_RTM_20260827_v6_BatchExcel
% Batch Excel spectral unmixing based on the v6 workflow.
%
% Workflow:
%   1. Select a folder containing Excel mixture datasets.
%   2. Enter one spectral-range string applied to ALL Excel files:
%        All   -> 0.35-2.5 um
%        SWIR  -> 1.0-2.5 um
%        SWIR2 -> 1.8-2.5 um
%        SWIR3 -> 2.0-2.5 um
%   3. For every Excel file:
%        - load mixture spectra / abundance / pure endmembers
%        - subset the ORIGINAL spectral grid only (NO resampling)
%        - fit each fraction independently and average k for mixtures
%        - run Hapke + MCMC and Shkuratov + MCMC
%        - save one combined *_MCMC_RTM_results.mat
%
% In addition to modelResults, the MAT file saves top-level modeled
% abundances and grain sizes for easy downstream Python/MATLAB analysis:
%   hFrac_m, sFrac_m, hSize_m, sSize_m
% and best-fit spectra:
%   h_fit, s_fit

clear; clc;
tic

%% User configuration
cfg = struct();

cfg.projectRoot = fileparts(mfilename('fullpath')); % portable code/output root

% Initial folder shown by uigetdir. The selected folder is processed in batch.
cfg.datasetSearchRoot = 'F:\Dataset20260909\Dataset\Reflectance_v0\';
cfg.datasetPattern = '*.xlsx';

% All results are written here. Range mode is included in each filename.
cfg.outputDir = fullfile(cfg.projectRoot, 'Results', 'MCMC_RTM_v6_BatchExcel_v3');

% Leave empty to infer components independently for EACH Excel file.
% To force the same subset/order for every file, set for example:
% cfg.componentNames = {'MGS-1', 'Alunite', 'Montmorillonite', 'Kaolinite'};
cfg.componentNames = {};

% Optical constants and densities.
cfg.knownEndmemberNames = {'Nontronite', 'Saponite', 'MGS-1', 'Kaolinite', 'Alunite', 'Montmorillonite'};
cfg.knownNvals = [1.57, 1.51, 1.5, 1.56, 1.58, 1.52];
cfg.knownDensity = [2.29, 2.30, 1.29, 2.62, 2.75, 2.50];

% Pure endmember grain-size range used for k inversion.
cfg.endmemberGrainSizeRange = [1, 74];  % positive effective-size bounds, um
% Example for one measured size fraction: [41, 63].
% Multi-fraction calibration is recommended. Run run_independent_k_from_excel first.
cfg.kCalibrationFile = fullfile(fileparts(mfilename('fullpath')), ...
    'Independent_k_Results', 'independent_k_calibration.mat');
% Set cfg.kCalibrationFile = ''; for single-spectrum profile mode.
cfg.kFitMode = 'independent'; % average independently fitted fraction k for mixtures
cfg.kNumStarts = 3;
cfg.kProfilePoints = 7;
cfg.kRMSEEquivalenceTolerance = 1e-7;
cfg.kAnchorIndex = [];
cfg.kAnchorDiameter = [];

% k inversion settings
cfg.kSolverDisplay = 'off';
cfg.kMaxIterations = 1000;
cfg.kFunctionTolerance = 1e-10;
cfg.kStepTolerance = 1e-10;
cfg.kOptimalityTolerance = 1e-10;

% RTM/MCMC settings
cfg.modelsToRun = {'Hapke', 'SK'};
cfg.Nm = 10000;
cfg.Nsteps = 50;
cfg.sizeLow = 1;
cfg.sizeHigh = 500;
cfg.priorAlpha = 1;
cfg.Cov = 5e-3;
cfg.randomSeed = 20260627;
cfg.saveFigures = true;

% Geometry and model constants
cfg.inc = 30;
cfg.emi = 0;
cfg.phs = 30;
cfg.q = 0.7;
cfg.isSpherical = 0;

% Set to a positive integer only if leading bands must be removed AFTER
% spectral-range selection.
cfg.dropLeadingBands = 0;

%% Paths
addpath(cfg.projectRoot);
addpath(fullfile(cfg.projectRoot, 'MCMC_unmixing_routines'));
% Root copies use the intended six-argument SK function and 0=irregular.
addpath(cfg.projectRoot, '-begin');
addpath(fileparts(mfilename('fullpath')), '-begin');
assert(nargin('Shkuratov_forward_model') == 6, ...
    'Wrong Shkuratov_forward_model on MATLAB path. Check which -all.');
cfg.forwardFunctionPaths = struct('hapke',which('Hapke_nk_to_ssa'), ...
    'shkuratov',which('Shkuratov_forward_model'), ...
    'rb_rf',which('Shkuratov_Calc_rb_rf'));
if ~isempty(cfg.kCalibrationFile) && ~isfile(cfg.kCalibrationFile)
    error('Run run_independent_k_from_excel first, or set cfg.kCalibrationFile to empty.');
end

if ~isfolder(cfg.outputDir)
    mkdir(cfg.outputDir);
end

%% Select input folder
inputFolder = uigetdir(cfg.datasetSearchRoot, ...
    'Select folder containing Excel mixture datasets');

if isequal(inputFolder, 0)
    error('No input folder selected.');
end

cfg.inputFolder = inputFolder;

excelFiles = dir(fullfile(inputFolder, cfg.datasetPattern));

% Fallback: if the strict pattern finds nothing, offer all Excel files.
if isempty(excelFiles)
    xlsxFiles = dir(fullfile(inputFolder, '*.xlsx'));
    xlsFiles = dir(fullfile(inputFolder, '*.xls'));
    excelFiles = [xlsxFiles; xlsFiles]; 
end

if isempty(excelFiles)
    error('No Excel files were found in: %s', inputFolder);
end

[~, sortIdx] = sort(lower(string({excelFiles.name})));
excelFiles = excelFiles(sortIdx);

fprintf('\nSelected input folder:\n%s\n', inputFolder);
fprintf('Found %d Excel files.\n', numel(excelFiles));

%% Spectrum range selection -- applied to every Excel file
rangeInput = strtrim(input( ...
    'Spectrum range? Enter All (0.35-2.5), SWIR (1.0-2.5), or SWIR2 (1.8-2.5), or SWIR3 (2.0-2.5): ', ...
    's'));

switch lower(rangeInput)
    case 'all'
        cfg.spectrumRangeMode = 'All';
        cfg.spectrumRangeUm = [0.35, 2.5];

    case 'swir'
        cfg.spectrumRangeMode = 'SWIR';
        cfg.spectrumRangeUm = [1.0, 2.5];

    case 'swir2'
        cfg.spectrumRangeMode = 'SWIR2';
        cfg.spectrumRangeUm = [1.8, 2.5];

    case 'swir3'
        cfg.spectrumRangeMode = 'SWIR3';
        cfg.spectrumRangeUm = [2.0, 2.5];

    otherwise
        error('Invalid range input "%s". Valid inputs are: All, SWIR, SWIR2, SWIR3.', rangeInput);
end

fprintf('\nSpectrum mode: %s (%.1f-%.1f um)\n', ...
    cfg.spectrumRangeMode, cfg.spectrumRangeUm(1), cfg.spectrumRangeUm(2));

% Preserve the user-level component setting. If empty, components are inferred
% independently for each Excel file rather than inherited from the first file.
requestedComponentNames = cfg.componentNames;

%% Batch processing
failedFiles = {};

for fileNo = 1:numel(excelFiles)

    datasetFile = fullfile(excelFiles(fileNo).folder, excelFiles(fileNo).name);
    [~, datasetStem] = fileparts(datasetFile);

    fprintf('\n============================================================\n');
    fprintf('[%d/%d] Processing: %s\n', fileNo, numel(excelFiles), excelFiles(fileNo).name);
    fprintf('============================================================\n');

    try
        % Make a per-file cfg so inferred component names never leak into the
        % next Excel file.
        cfgFile = cfg;
        cfgFile.datasetFile = datasetFile;
        cfgFile.componentNames = requestedComponentNames;

        rng(cfgFile.randomSeed, 'twister');

        %% Load one Excel dataset
        [wlAll, ObsAll, Abu, Abu_wtpct, spectrumNames, abundanceNames, pureEndmembers] = ...
            load_mixture_dataset(datasetFile);

        %% Select wavelength range -- NO spectral resampling
        rangeMask = wlAll >= cfgFile.spectrumRangeUm(1) & ...
                    wlAll <= cfgFile.spectrumRangeUm(2);

        assert(any(rangeMask), ...
            'No wavelength bands found in selected range %.2f-%.2f um.', ...
            cfgFile.spectrumRangeUm(1), cfgFile.spectrumRangeUm(2));

        wl = wlAll(rangeMask);
        Obs = ObsAll(rangeMask, :);
        pureEndmembers.reflectance = pureEndmembers.reflectance(rangeMask, :);

        if cfgFile.dropLeadingBands > 0
            keepIdxBands = (cfgFile.dropLeadingBands + 1):numel(wl);
            wl = wl(keepIdxBands);
            Obs = Obs(keepIdxBands, :);
            pureEndmembers.reflectance = pureEndmembers.reflectance(keepIdxBands, :);
        end

        wl = wl(:);

        %% Components / constants for this Excel file
        if isempty(cfgFile.componentNames)
            cfgFile.componentNames = abundanceNames;
        end

        cfgFile.componentNames = cellfun( ...
            @full_endmember_name, ...
            cfgFile.componentNames, ...
            'UniformOutput', false);

       

        abuIdx = match_names(cfgFile.componentNames, abundanceNames, 'abundance table');
        constIdx = match_names(cfgFile.componentNames, cfgFile.knownEndmemberNames, 'known optical constants');

        Endm = pureEndmembers.reflectance(:, abuIdx);
        Endm_Name = pureEndmembers.names(abuIdx);

        pureEndmembers.names = pureEndmembers.names(abuIdx);
        pureEndmembers.spectrumNames = pureEndmembers.spectrumNames(abuIdx);
        pureEndmembers.reflectance = Endm;

        n_vals = cfgFile.knownNvals(constIdx);
        density = cfgFile.knownDensity(constIdx);

        Abu = Abu(:, abuIdx);
        Abu_wtpct = Abu_wtpct(:, abuIdx);
        abundanceNames = abundanceNames(abuIdx);

        M = numel(cfgFile.componentNames);
        N = size(Obs, 2);
        L = numel(wl);

        assert(size(Endm, 2) == M, ...
            'Endmember count does not match selected component count.');
        assert(size(Abu, 1) == N, ...
            'Abundance rows do not match spectra columns.');
        assert(size(Obs, 1) == L, ...
            'Observed spectra and wavelength lengths differ.');

        grainSizeBounds = repmat(cfgFile.endmemberGrainSizeRange, M, 1);
        GS = mean(grainSizeBounds, 2);

        fprintf('Components: %s\n', strjoin(cfgFile.componentNames, ', '));
        fprintf('Mixture spectra: %d; wavelength bands: %d\n', N, L);
        fprintf('Actual wavelength range: %.4f-%.4f um\n', wl(1), wl(end));

        %% Calculate k directly on the selected ORIGINAL wavelength grid
        fprintf('\nCalculating h_k and s_k on original spectral grid...\n');
        [h_k, s_k, hD, sD, kDiagnostics] = compute_endmember_k( ...
            wl, Endm, n_vals, grainSizeBounds, cfgFile);

        runStem = string(datasetStem) + "_" + string(cfgFile.spectrumRangeMode);

        plot_endmember_reflectance_and_k( ...
            wl, Endm, h_k, s_k, Endm_Name, cfgFile, runStem);

        %% Common data structure for likelihood functions
        Dbase = struct();
        Dbase.n = repmat(n_vals(:).', L, 1);
        Dbase.rho = density(:);
        Dbase.inc = cfgFile.inc;
        Dbase.emi = cfgFile.emi;
        Dbase.phs = cfgFile.phs;
        Dbase.q = cfgFile.q;
        Dbase.isSpherical = cfgFile.isSpherical;
        Dbase.Cinv = 1 / cfgFile.Cov;
        Dbase.Ncomponents = M;
        Dbase.prior_size_low = cfgFile.sizeLow * ones(1, M);
        Dbase.prior_size_high = cfgFile.sizeHigh * ones(1, M);
        Dbase.prior_alpha = cfgFile.priorAlpha * ones(1, M);
        Dbase.Nm = cfgFile.Nm;
        Dbase.Nsteps = cfgFile.Nsteps;
        Dbase.lambda = wl(:);
        Dbase.lam_SPEC = wl(:);

        assert(numel(Dbase.lambda) == size(Obs, 1), ...
            'D.lambda length does not match observed spectrum length.');
        assert(numel(Dbase.lam_SPEC) == size(Obs, 1), ...
            'D.lam_SPEC length does not match observed spectrum length.');

        %% Run selected RTM models
        modelResults = struct();

        for modelNo = 1:numel(cfgFile.modelsToRun)
            modelName = cfgFile.modelsToRun{modelNo};

            switch lower(modelName)
                case 'hapke'
                    D = Dbase;
                    D.k = h_k;
                    modelSpec = struct( ...
                        'name', 'Hapke', ...
                        'likelihood', @Hapke_likelihood, ...
                        'fitVariable', 'h_fit');

                case 'sk'
                    D = Dbase;
                    D.k = s_k;
                    modelSpec = struct( ...
                        'name', 'SK', ...
                        'likelihood', @sk_likelihood, ...
                        'fitVariable', 's_fit');

                otherwise
                    error('Unknown model name: %s', modelName);
            end

            fprintf('\nRunning %s model for %d spectra, %d components.\n', ...
                modelSpec.name, N, M);

            modelResults.(modelSpec.name) = ...
                run_mcmc_model(modelSpec, D, Obs);
        end

        %% Top-level outputs for convenient downstream analysis
        % These duplicate the key values stored inside modelResults by design.
        hFrac_m = [];
        sFrac_m = [];
        hSize_m = [];
        sSize_m = [];
        h_fit = [];
        s_fit = [];

        if isfield(modelResults, 'Hapke')
            hFrac_m = modelResults.Hapke.massFraction;
            hSize_m = modelResults.Hapke.size;
            h_fit = modelResults.Hapke.fit_spec;
        end

        if isfield(modelResults, 'SK')
            sFrac_m = modelResults.SK.massFraction;
            sSize_m = modelResults.SK.size;
            s_fit = modelResults.SK.fit_spec;
        end

        plot_abundance_scatter( ...
            Abu, modelResults, Endm_Name, cfgFile, runStem);

        %% Save one combined result per Excel file
        combinedOutpath = fullfile( ...
            cfgFile.outputDir, ...
            string(datasetStem) + "_" + string(cfgFile.spectrumRangeMode) + ...
            "_MCMC_RTM_results.mat");

        save(combinedOutpath, ...
            'cfgFile', 'modelResults', ...
            'wl', 'Obs', 'Abu', 'Abu_wtpct', ...
            'Endm', 'Endm_Name', 'GS', 'density', 'n_vals', ...
            'h_k', 's_k', 'hD', 'sD', 'kDiagnostics', 'grainSizeBounds', ...
            'hFrac_m', 'sFrac_m', ...
            'hSize_m', 'sSize_m', ...
            'h_fit', 's_fit', ...
            'spectrumNames', 'abundanceNames', 'pureEndmembers', ...
            '-v7.3');

        fprintf('\nCombined result saved to:\n%s\n', combinedOutpath);

    catch ME
        failedFiles{end + 1, 1} = excelFiles(fileNo).name; %#ok<AGROW>
        failedFiles{end, 2} = ME.message;

        fprintf(2, '\nFAILED: %s\n%s\n', excelFiles(fileNo).name, ME.message);
    end
end

%% Summary
fprintf('\n============================================================\n');
fprintf('Batch processing completed.\n');
fprintf('Total Excel files: %d\n', numel(excelFiles));
fprintf('Successful: %d\n', numel(excelFiles) - size(failedFiles, 1));
fprintf('Failed: %d\n', size(failedFiles, 1));
fprintf('Output directory:\n%s\n', cfg.outputDir);

if ~isempty(failedFiles)
    fprintf('\nFailed files:\n');
    for k = 1:size(failedFiles, 1)
        fprintf('  %s\n    %s\n', failedFiles{k, 1}, failedFiles{k, 2});
    end
end

toc

%% Local functions
function [wl, Obs, Abu, Abu_wtpct, spectrumNames, abundanceNames, pureEndmembers] = load_mixture_dataset(filePath)
    [~, ~, ext] = fileparts(filePath);

    switch lower(ext)
        case '.mat'
            error('v3 calculates k from pure endmembers in Excel files. Please select an Excel dataset.');

        case {'.xlsx', '.xls'}
            spectraTable = readtable(filePath, ...
                'Sheet', 1, ...
                'VariableNamingRule', 'preserve');

            abundanceTable = readtable(filePath, ...
                'Sheet', 2, ...
                'VariableNamingRule', 'preserve');

            wl = spectraTable{:, 1};
            Obs = spectraTable{:, 2:end};
            spectrumNames = spectraTable.Properties.VariableNames(2:end);
            [Obs, spectrumNames] = remove_empty_spectra_columns(Obs, spectrumNames);

            allNames = abundanceTable.Properties.VariableNames;
            abuCols = find(startsWith(allNames, 'Abu_', 'IgnoreCase', true));
            assert(~isempty(abuCols), 'No abundance columns found. Expected columns like Abu_MGS-1_wtpct.');

            abundanceNames = regexprep(allNames(abuCols), '^Abu_', '', 'ignorecase');
            abundanceNames = regexprep(abundanceNames, '_wtpct$', '', 'ignorecase');

            Abu_wtpct_all = abundanceTable{:, abuCols};
            [Abu_all, Abu_wtpct_all] = normalize_abundance_scale(Abu_wtpct_all);

            spectrumCol = find(strcmpi(allNames, 'Spectrum'), 1);
            if ~isempty(spectrumCol)
                tableSpectrumNames = string(abundanceTable{:, spectrumCol});
                [~, loc] = match_spectra_to_abundance_rows(spectrumNames, tableSpectrumNames);
                Abu = Abu_all(loc, :);
                Abu_wtpct = Abu_wtpct_all(loc, :);
            else
                assert(size(Abu_wtpct_all, 1) == numel(spectrumNames), ...
                    'Abundance row count must match spectra column count when Spectrum column is absent.');
                Abu = Abu_all;
                Abu_wtpct = Abu_wtpct_all;
            end

            [Obs, Abu, Abu_wtpct, spectrumNames, pureEndmembers] = ...
                split_pure_endmember_spectra(Obs, Abu, Abu_wtpct, spectrumNames, abundanceNames);

        case '.csv'
            error('v3 calculates k from pure endmembers in Excel files. Please select an Excel dataset.');

        otherwise
            error('Unsupported data file type: %s', ext);
    end

    wl = wl(:);
    Obs = double(Obs);
    Abu = double(Abu);
    Abu_wtpct = double(Abu_wtpct);

    if wl(end) > 1000
        wl = wl / 1000;
    end
end

function [Obs, spectrumNames] = remove_empty_spectra_columns(Obs, spectrumNames)
    if isempty(Obs)
        return;
    end

    emptyNameMask = strlength(strtrim(string(spectrumNames))) == 0;
    emptyDataMask = all(isnan(Obs), 1);
    keepMask = ~(emptyNameMask | emptyDataMask);

    if any(~keepMask)
        fprintf('Removed %d empty spectra columns.\n', nnz(~keepMask));
        Obs = Obs(:, keepMask);
        spectrumNames = spectrumNames(keepMask);
    end
end

function [tf, loc] = match_spectra_to_abundance_rows(spectrumNames, tableSpectrumNames)
    spectrumNames = string(spectrumNames(:));
    tableSpectrumNames = string(tableSpectrumNames(:));

    [tf, loc] = ismember(spectrumNames, tableSpectrumNames);
    if all(tf)
        return;
    end

    canonSpectra = arrayfun(@canonical_spectrum_name, spectrumNames, 'UniformOutput', false);
    canonTable = arrayfun(@canonical_spectrum_name, tableSpectrumNames, 'UniformOutput', false);
    [tfCanon, locCanon] = ismember(canonSpectra, canonTable);

    if all(tfCanon)
        tf = tfCanon;
        loc = locCanon;
        warning('Matched spectra to abundance rows after canonicalizing names.');
        return;
    end

    if numel(spectrumNames) == numel(tableSpectrumNames)
        tf = true(size(spectrumNames));
        loc = (1:numel(spectrumNames)).';
        missingNames = strjoin(cellstr(spectrumNames(~tfCanon)), ', ');
        warning(['Some spectra names were not found in the abundance table; ', ...
            'using row/column order instead. Check these names: %s'], missingNames);
        return;
    end

    missingNames = strjoin(cellstr(spectrumNames(~tfCanon)), ', ');
    error('Some spectra columns were not found in the abundance table: %s', missingNames);
end

function out = canonical_spectrum_name(name)
    out = lower(char(string(name)));
    out = regexprep(out, '[^a-z0-9]+', '');
end
function [Abu, Abu_wtpct] = normalize_abundance_scale(rawAbu)
    rawAbu = double(rawAbu);

    finiteVals = rawAbu(isfinite(rawAbu));
    if isempty(finiteVals)
        Abu = rawAbu;
        Abu_wtpct = rawAbu;
        return;
    end

    if max(finiteVals) <= 1 + 1e-8
        Abu = rawAbu;
        Abu_wtpct = 100 * rawAbu;
    else
        Abu_wtpct = rawAbu;
        Abu = rawAbu / 100;
    end
end

function [Obs, Abu, Abu_wtpct, spectrumNames, pureEndmembers] = split_pure_endmember_spectra(Obs, Abu, Abu_wtpct, spectrumNames, abundanceNames)
    pureEndmembers = struct();
    pureEndmembers.names = abundanceNames;
    pureEndmembers.spectrumNames = cell(1, numel(abundanceNames));
    pureEndmembers.reflectance = zeros(size(Obs, 1), numel(abundanceNames));

    if isempty(Abu)
        error('Abundance matrix is empty. Cannot extract pure endmembers.');
    end

    tol = 1e-10;
    pureMask = any(Abu >= 1 - tol, 2);

    assert(numel(pureMask) == size(Obs, 2), ...
        'Pure-spectrum filter cannot be applied because abundance rows and spectra columns differ.');

    for componentNo = 1:numel(abundanceNames)
        hits = find(Abu(:, componentNo) >= 1 - tol);
        if isempty(hits)
            error('No pure endmember spectrum found for component: %s', abundanceNames{componentNo});
        end

        hit = hits(1);
        pureEndmembers.reflectance(:, componentNo) = Obs(:, hit);
        pureEndmembers.spectrumNames{componentNo} = char(string(spectrumNames{hit}));
    end
 
    keepMask = ~pureMask;

    Obs = Obs(:, keepMask);
    Abu = Abu(keepMask, :);
    Abu_wtpct = Abu_wtpct(keepMask, :);
    spectrumNames = spectrumNames(keepMask);

    fprintf('Extracted %d pure endmember spectra and removed them from mixtures.\n', nnz(pureMask));
end

function idx = match_names(requestedNames, availableNames, sourceName)
    % Match endmember names while allowing both full names and abbreviations.
    %
    % Supported equivalences:
    %   Nontronite       <-> NAu
    %   Saponite         <-> Sap
    %   MGS-1            <-> MGS
    %   Kaolinite        <-> Kao
    %   Alunite          <-> Alu
    %   Montmorillonite  <-> Mont

    idx = zeros(1, numel(requestedNames));

    % Convert all available names to standardized canonical names.
    canonAvailable = cellfun( ...
        @canonical_endmember_name, ...
        availableNames, ...
        'UniformOutput', false);

    for k = 1:numel(requestedNames)

        % Convert requested name to same standardized canonical name.
        target = canonical_endmember_name(requestedNames{k});

        hit = find(strcmp(target, canonAvailable), 1);

        if isempty(hit)
            error( ...
                'Component "%s" was not found in %s. Available names: %s', ...
                requestedNames{k}, ...
                sourceName, ...
                strjoin(availableNames, ', '));
        end

        idx(k) = hit;
    end
end


function name = canonical_endmember_name(name)
    % Standardize full names and abbreviations to one canonical form.
    %
    % Canonical outputs:
    %   nontronite
    %   saponite
    %   mgs1
    %   kaolinite
    %   alunite
    %   montmorillonite

    name = lower(char(string(name)));

    % Remove spaces, hyphens, underscores and other punctuation.
    name = regexprep(name, '[^a-z0-9]+', '');

    switch name

        % ====================================================
        % Nontronite
        % ====================================================
        case {'nontronite', 'nau'}
            name = 'nontronite';

        % ====================================================
        % Saponite
        % ====================================================
        case {'saponite', 'sap'}
            name = 'saponite';

        % ====================================================
        % MGS-1
        % ====================================================
        case {'mgs1', 'mgs'}
            name = 'mgs1';

        % ====================================================
        % Kaolinite
        % ====================================================
        case {'kaolinite', 'kao'}
            name = 'kaolinite';

        % ====================================================
        % Alunite
        % ====================================================
        case {'alunite', 'alu'}
            name = 'alunite';

        % ====================================================
        % Montmorillonite
        % ====================================================
        case {'montmorillonite', 'mont'}
            name = 'montmorillonite';

        % ====================================================
        % Unknown name
        % ====================================================
        otherwise
            % Keep normalized text unchanged.
            % This preserves compatibility with possible
            % additional endmembers in the future.
    end
end

function [h_k,s_k,hD,sD,kDiagnostics] = compute_endmember_k(wl,Endm,n_vals,grainSizeBounds,cfg)
    [h_k,s_k,hD,sD,kDiagnostics] = compute_endmember_k_bounded( ...
        wl,Endm,n_vals,grainSizeBounds,cfg);
end

function result = run_mcmc_model(modelSpec, D, Obs)
    [L, N] = size(Obs);
    M = D.Ncomponents;

    llk_prior_handle = @llk_prior;
    for_model = modelSpec.likelihood;

    all_MODS = cell(1, N);
    all_ind = zeros(1, N);
    bestmod = zeros(2 * M, N);
    bestLogLike = zeros(1, N);
    acceptCount = cell(1, N);
    rejectCount = cell(1, N);

    totalTimer = tic;

    for specNo = 1:N
        fprintf('[%s] Spectrum %d/%d\n', modelSpec.name, specNo, N);

        D.R_SPEC = Obs(:, specNo);
        THETA = initialize_theta(D);

        [THETA, LLK, nAccept, nReject] = catmip_dirichlet( ...
            D.Nm, D.Nsteps, THETA, D, llk_prior_handle, for_model);

        finalLogLike = LLK(2, :, end);
        [bestLogLike(specNo), ind] = max(finalLogLike);

        MODS = THETA(:, :, end);
        bestmod(:, specNo) = MODS(:, ind);

        all_MODS{specNo} = MODS;
        all_ind(specNo) = ind;
        acceptCount{specNo} = nAccept;
        rejectCount{specNo} = nReject;
    end

    fit_spec = synthesize_fit(modelSpec.name, D, bestmod, L, N);
    residual = Obs - fit_spec;
    rmse = sqrt(mean(residual .^ 2, 1)).';

    result = struct();
    result.modelName = modelSpec.name;
    result.all_MODS = all_MODS;
    result.all_ind = all_ind;
    result.bestmod = bestmod;
    result.bestLogLike = bestLogLike;
    result.fit_spec = fit_spec;
    result.residual = residual;
    result.rmse = rmse;
    result.size = bestmod(1:M, :).';
    result.massFraction = bestmod(M+1:2*M, :).';
    result.acceptCount = acceptCount;
    result.rejectCount = rejectCount;
    result.elapsedSeconds = toc(totalTimer);
end

function THETA = initialize_theta(D)
    M = D.Ncomponents;

    low = D.prior_size_low(:);
    high = D.prior_size_high(:);

    sizes = bsxfun(@plus, low, bsxfun(@times, high - low, rand(M, D.Nm)));
    abundances = dirichletrnd_exc(D.prior_alpha, D.Nm);

    if size(abundances, 1) ~= M && size(abundances, 2) == M
        abundances = abundances.';
    end

    assert(size(abundances, 1) == M, 'Initial abundance matrix has wrong component count.');

    THETA = [sizes; abundances];
end

function fit_spec = synthesize_fit(modelName, D, bestmod, L, N)
    M = D.Ncomponents;
    fit_spec = zeros(L, N);

    for specNo = 1:N
        sizes = bestmod(1:M, specNo);
        wtperc = bestmod(M+1:2*M, specNo);
        fi = wtperctof(D.rho, sizes, wtperc);

        switch lower(modelName)
            case 'hapke'
                ssa_model = Hapke_forward_model([fi; sizes], D.n, D.k, D.lambda, D.isSpherical);
                fit_spec(:, specNo) = Hapke_ssa2reff(ssa_model, D.inc, D.emi, D.phs);

            case 'sk'
                fit_spec(:, specNo) = Shkuratov_forward_model([fi; sizes], D.n, D.k, D.lambda, D.q, D.isSpherical);

            otherwise
                error('Unknown model name: %s', modelName);
        end
    end
end

function plot_endmember_reflectance_and_k(wl, Endm, h_k, s_k, Endm_Name, cfg, datasetStem)
    M = numel(Endm_Name);
    fig = figure('Name', 'Endmember reflectance and k', 'Color', 'w');
    tiledlayout(2, M, 'TileSpacing', 'compact', 'Padding', 'compact');

    for componentNo = 1:M
        nexttile(componentNo);
        plot(wl, Endm(:, componentNo), 'k', 'LineWidth', 1.1);
        ylabel('Reflectance');
        xlabel('Wavelength (um)');
        xlim([min(wl), max(wl)]);
        title(Endm_Name{componentNo}, 'Interpreter', 'none');
        grid on;
        box on;

        nexttile(M + componentNo);
        semilogy(wl, h_k(:, componentNo), 'b', 'LineWidth', 1.0);
        hold on;
        semilogy(wl, s_k(:, componentNo), 'r', 'LineWidth', 1.0);
        ylabel('k');
        xlabel('Wavelength (um)');
        xlim([min(wl), max(wl)]);
        grid on;
        box on;

        if componentNo == M
            legend({'Hapke', 'SK99'}, 'Location', 'best');
        end
    end

    save_figure_if_requested(fig, cfg, datasetStem, '01_endmember_reflectance_and_k');
end

function plot_input_spectra(wl, Obs, spectrumNames, cfg, datasetStem)
    fig = figure('Name', 'Observed mixture spectra', 'Color', 'w');
    plot(wl, Obs, 'LineWidth', 0.8);
    xlabel('Wavelength (um)');
    ylabel('Reflectance');
    if isfield(cfg, 'currentResampleStepNm')
        title(sprintf('Observed Mixture Spectra (%s %.2f-%.2f um, %d nm sampling)', ...
            cfg.spectrumRangeMode, cfg.spectrumRangeUm(1), cfg.spectrumRangeUm(2), ...
            cfg.currentResampleStepNm), 'Interpreter', 'none');
    else
        title(sprintf('Observed Mixture Spectra (%s %.2f-%.2f um)', ...
            cfg.spectrumRangeMode, cfg.spectrumRangeUm(1), cfg.spectrumRangeUm(2)), ...
            'Interpreter', 'none');
    end
    xlim([min(wl), max(wl)]);
    grid on;
    box on;

    if numel(spectrumNames) <= 12
        legend(spectrumNames, 'Interpreter', 'none', 'Location', 'best');
    end

    save_figure_if_requested(fig, cfg, datasetStem, '02_mixture_spectra');
end

function plot_abundance_scatter(Abu, modelResults, Endm_Name, cfg, datasetStem)
    modelNames = fieldnames(modelResults);
    if isempty(modelNames)
        return;
    end

    M = size(Abu, 2);
    fig = figure('Name', 'Observed vs modeled abundance diagnostics', 'Color', 'w');
    fig.WindowState = 'maximized';
    tiledlayout(2, M, 'TileSpacing', 'compact', 'Padding', 'compact');

    markerSpec = struct();
    markerSpec.Hapke.marker = 'o';
    markerSpec.Hapke.color = [0.12, 0.46, 0.71];
    markerSpec.SK.marker = 's';
    markerSpec.SK.color = [1.00, 0.49, 0.05];

    for componentNo = 1:M
        nexttile(componentNo);
        hold on;

        scatterHandles = gobjects(0);
        scatterLabels = {};
        rmseLines = {};

        for modelNo = 1:numel(modelNames)
            modelName = modelNames{modelNo};
            modeled = modelResults.(modelName).massFraction(:, componentNo);

            if isfield(markerSpec, modelName)
                marker = markerSpec.(modelName).marker;
                faceColor = markerSpec.(modelName).color;
            else
                marker = 'o';
                faceColor = [0.35, 0.35, 0.35];
            end

            h = scatter( ...
                Abu(:, componentNo), modeled, 35, marker, ...
                'MarkerEdgeColor', [0, 0, 0], ...
                'MarkerFaceColor', faceColor);
            scatterHandles(end + 1) = h; %#ok<AGROW>
            scatterLabels{end + 1} = model_display_name(modelName); %#ok<AGROW>

            rmseVal = sqrt(mean((Abu(:, componentNo) - modeled).^2));
            if strcmpi(modelName, 'Hapke')
                rmseLines{end + 1} = sprintf('rmse_h = %.3f', rmseVal); %#ok<AGROW>
            elseif strcmpi(modelName, 'SK')
                rmseLines{end + 1} = sprintf('rmse_s = %.3f', rmseVal); %#ok<AGROW>
            else
                rmseLines{end + 1} = sprintf('rmse_%s = %.3f', modelName, rmseVal); %#ok<AGROW>
            end
        end

        lineHandles = plot_abundance_reference_lines();
        xlabel('Observed Abundance');
        ylabel('Modeled Abundance');
        title(Endm_Name{componentNo}, 'Interpreter', 'none');
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
        box on;
        xlim([0, 1]);
        ylim([0, 1]);
        daspect([0.9, 1, 1]);

        text(0.04, 0.96, strjoin(rmseLines, newline), ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', [0.75, 0.75, 0.75], ...
            'Margin', 4, ...
            'FontSize', 9);

        if componentNo == M
            legend( ...
                [scatterHandles(:); lineHandles(:)], ...
                [scatterLabels, {'1:1 Line', '5% offset', '10% offset', '20% offset'}], ...
                'Location', 'southeast');
        end

        nexttile(M + componentNo);
        hold on;
        for modelNo = 1:numel(modelNames)
            modelName = modelNames{modelNo};
            modeled = modelResults.(modelName).massFraction(:, componentNo);
            residualAbu = Abu(:, componentNo) - modeled;

            if isfield(markerSpec, modelName)
                marker = markerSpec.(modelName).marker;
                faceColor = markerSpec.(modelName).color;
            else
                marker = 'o';
                faceColor = [0.35, 0.35, 0.35];
            end

            scatter( ...
                Abu(:, componentNo), residualAbu, 35, marker, ...
                'MarkerEdgeColor', [0, 0, 0], ...
                'MarkerFaceColor', faceColor);
        end
        yline(0, 'k', 'LineWidth', 1.0);
        yline(0.05, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-.');
        yline(-0.05, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-.');
        yline(0.10, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--');
        yline(-0.10, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--');
        yline(0.20, 'Color', [0.3, 0.3, 0.3], 'LineStyle', ':');
        yline(-0.20, 'Color', [0.3, 0.3, 0.3], 'LineStyle', ':');

        xlabel('Observed Abundance');
        ylabel('Observed - Modeled');
        xlim([0, 1]);
        ylim([-0.5, 0.5]);
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
        grid on;
        box on;
    end

    save_figure_if_requested(fig, cfg, datasetStem, '03_abundance_diagnostics');
end

function lineHandles = plot_abundance_reference_lines()
    x = [0, 1];

    h1 = plot(x, x, 'k', 'LineWidth', 1.0);

    h5a = plot(x, x + 0.05, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-.');
    plot(x, x - 0.05, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-.');

    h10a = plot(x, x + 0.10, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--');
    plot(x, x - 0.10, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--');

    h20a = plot(x, x + 0.20, 'Color', [0.3, 0.3, 0.3], 'LineStyle', ':');
    plot(x, x - 0.20, 'Color', [0.3, 0.3, 0.3], 'LineStyle', ':');

    lineHandles = [h1; h5a; h10a; h20a];
end
function name = model_display_name(modelName)
    if strcmpi(modelName, 'SK')
        name = 'SK99';
    else
        name = modelName;
    end
end

function save_figure_if_requested(fig, cfg, datasetStem, suffix)
    if ~isfield(cfg, 'saveFigures') || ~cfg.saveFigures
        return;
    end

    figDir = fullfile(cfg.outputDir, 'Figures');
    if ~isfolder(figDir)
        mkdir(figDir);
    end

    figPath = fullfile(figDir, datasetStem + "_" + suffix + ".fig");
    pngPath = fullfile(figDir, datasetStem + "_" + suffix + ".png");

    savefig(fig, figPath);

    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, pngPath, 'Resolution', 300);
    else
        saveas(fig, pngPath);
    end
end

function name = full_endmember_name(name)
    % Convert abbreviations or full names to standardized full names.

    canon = canonical_endmember_name(name);

    switch canon

        case 'nontronite'
            name = 'Nontronite';

        case 'saponite'
            name = 'Saponite';

        case 'mgs1'
            name = 'MGS-1';

        case 'kaolinite'
            name = 'Kaolinite';

        case 'alunite'
            name = 'Alunite';

        case 'montmorillonite'
            name = 'Montmorillonite';

        otherwise
            name = char(string(name));
    end
end
