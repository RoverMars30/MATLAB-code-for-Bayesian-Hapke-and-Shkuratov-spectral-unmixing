function [h_k,s_k,hD,sD,diagnostics] = compute_endmember_k_bounded(wl,Endm,n_vals,grainSizeBounds,cfg)
% Batch adapter. Independent fraction fits return their arithmetic-mean k.
% Endm has one column per DIFFERENT mineral.
% cfg.kCalibrationFile optionally supplies multiple fractions per mineral:
% calibration.wl, .names, .reflectance{m} (L-by-J), .bounds{m} (J-by-2).
% Without that file each mineral is treated as a single-spectrum profile.
% hD/sD are cells: one diameter for each calibration fraction of a mineral.
% Never pass these calibration sizes as mixture size estimates.

persistent savedCalls savedResults
if isempty(savedCalls), savedCalls={}; savedResults={}; end
wl=wl(:); [L,M]=size(Endm);
if L~=numel(wl) || numel(n_vals)~=M || ~isequal(size(grainSizeBounds),[M 2])
    error('Endm, n_vals and grainSizeBounds have inconsistent dimensions.');
end
h_k=zeros(L,M); s_k=h_k; hD=cell(1,M); sD=hD;
diagnostics=struct('Hapke',{cell(1,M)},'SK',{cell(1,M)},...
    'calibrationNames',{cell(1,M)});
cal=[]; source='';
if isfield(cfg,'kCalibrationFile') && ~isempty(cfg.kCalibrationFile)
    z=load(cfg.kCalibrationFile,'calibration');
    if ~isfield(z,'calibration'), error('MAT file must contain calibration.'); end
    cal=z.calibration;
    source=cfg.kCalibrationFile;
    if ~all(isfield(cal,{'wl','names','reflectance','bounds'}))
        error('Calibration needs wl, names, reflectance and bounds.');
    end
end
for m=1:M
    name=cfg.componentNames{m}; R=Endm(:,m); B=grainSizeBounds(m,:);
    if ~isempty(cal)
        idx=find(strcmpi(string(cal.names),string(name)));
        if numel(idx)~=1, error('No unique calibration entry for %s.',name); end
        % Match wavelength channels; do not silently interpolate spectra.
        [found,rows]=ismembertol(wl,cal.wl(:),1e-9,'DataScale',1);
        if ~all(found), error('Calibration grid does not contain the requested channels.'); end
        R=cal.reflectance{idx}(rows,:); B=cal.bounds{idx};
    end
    if any(B(:,1)<=0)
        error('Use a positive lower bound (e.g. [1 74], not [0 74]); D=0 is singular.');
    end
    opt=struct(); fields={'inc','emi','phs','q','isSpherical',...
        'kLowerBound','kUpperBound','kNumStarts','kMaxIterations','kSolverDisplay',...
        'kFunctionTolerance','kStepTolerance','kOptimalityTolerance',...
        'kRMSEEquivalenceTolerance','kProfilePoints','kDiameterTolerance',...
        'kAnchorIndex','kAnchorDiameter'};
    for a=1:numel(fields)
        if isfield(cfg,fields{a}), opt.(fields{a})=cfg.(fields{a}); end
    end
    mode='independent';
    if isfield(cfg,'kFitMode'), mode=cfg.kFitMode; end
    if ~any(strcmpi(mode,{'independent','joint'})), error('Unknown kFitMode.'); end
    call=struct('wl',wl,'R',R,'n',n_vals(m),'B',B,'opt',opt,'source',source,'mode',mode);
    hit=[];
    for a=1:numel(savedCalls)
        if isequaln(call,savedCalls{a}), hit=a; break; end
    end
    if isempty(hit)
        fprintf('%s optical calibration: %s, %d fraction(s).\n',mode,name,size(R,2));
        if strcmpi(mode,'independent')
            [h,s]=fit_endmember_k_by_fraction(wl,R,n_vals(m),B,opt);
        else
            h=invert_endmember_k_size(wl,R,n_vals(m),B,'Hapke',opt);
            s=invert_endmember_k_size(wl,R,n_vals(m),B,'SK',opt);
        end
        savedCalls{end+1}=call; savedResults{end+1}={h,s};
    else
        h=savedResults{hit}{1}; s=savedResults{hit}{2};
    end
    if strcmpi(mode,'independent')
        h_k(:,m)=h.meanK; s_k(:,m)=s.meanK;
    else
        h_k(:,m)=h.k; s_k(:,m)=s.k;
    end
    hD{m}=h.diameter; sD{m}=s.diameter;
    diagnostics.Hapke{m}=h; diagnostics.SK{m}=s;
    diagnostics.calibrationNames{m}=name;
end
end
