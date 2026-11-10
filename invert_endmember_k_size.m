function result = invert_endmember_k_size(wl, R, n, bounds, model, cfg)
% Shared k(lambda), one effective diameter per measured size fraction.
% wl and diameters MUST use micrometres. R is L-by-J for ONE material.
% bounds is J-by-2, e.g. [41 63; 63 74]. Both models are independent fits.
% Single-spectrum mode profiles out k and reports a solution FAMILY.
% Multiple spectra constrain relative sizes but can retain scale ambiguity.
% Optional cfg.kAnchorIndex / cfg.kAnchorDiameter fixes one independently
% justified diameter. No smoothing or size prior is added implicitly.
% Requires Optimization Toolbox (MATLAB R2022a compatible).

if nargin < 6, cfg = struct(); end
cfg = defaults(cfg);
wl = wl(:);
if isvector(R), R = R(:); end
[L,J] = size(R);
validateattributes(wl, {'numeric'}, {'real','finite','positive'});
validateattributes(R, {'numeric'}, {'real','finite','positive','<',1});
validateattributes(n, {'numeric'}, {'scalar','real','finite','>',1});
if L ~= numel(wl) || L < 3 || ~isequal(size(bounds),[J 2])
    error('kSize:Dimensions','Need L-by-J reflectances and J-by-2 size bounds.');
end
validateattributes(bounds, {'numeric'}, {'real','finite','positive'});
if any(bounds(:,2) < bounds(:,1)), error('Reversed diameter bounds.'); end
if any(diff(wl)<=0), error('Wavelengths must be strictly increasing.'); end
if ~any(strcmpi(model,{'Hapke','SK'})), error('Use Hapke or SK.'); end
if exist('lsqnonlin','file')~=2, error('Optimization Toolbox is required.'); end
originalBounds = bounds;
if ~isempty(cfg.kAnchorIndex)
    a=cfg.kAnchorIndex; d=cfg.kAnchorDiameter;
    validateattributes(a,{'numeric'},{'scalar','integer','>=',1,'<=',J});
    validateattributes(d,{'numeric'},{'scalar','real','finite','positive'});
    if d<bounds(a,1) || d>bounds(a,2), error('Anchor outside diameter bounds.'); end
    bounds(a,:)=[d d];
end
D0=mean(bounds,2);
free=find(bounds(:,2)>bounds(:,1));
logLo=log(cfg.kLowerBound); logHi=log(cfg.kUpperBound);
z0=log(min(max(1e-5./mean(R,2),cfg.kLowerBound),cfg.kUpperBound));
patternK=repmat(speye(L),J,1);
baseOptions=optimoptions('lsqnonlin','Display',cfg.kSolverDisplay,...
    'MaxIterations',cfg.kMaxIterations,'MaxFunctionEvaluations',100000,...
    'FunctionTolerance',cfg.kFunctionTolerance,...
    'StepTolerance',cfg.kStepTolerance,...
    'OptimalityTolerance',cfg.kOptimalityTolerance);

profile=struct(); starts=struct([]);
if J==1 && ~isempty(free)
    % For each candidate D, fit all wavelength-dependent k values. The outer
    % bounded search is variable projection of the joint least-squares loss.
    grid=unique([linspace(bounds(1),bounds(2),cfg.kProfilePoints),D0]);
    kvals=zeros(L,numel(grid)); losses=zeros(size(grid)); flags=zeros(size(grid));
    for j=1:numel(grid)
        [kvals(:,j),losses(j),flags(j)]=fit_at_size(grid(j));
    end
    outer=optimset('Display','off','TolX',cfg.kDiameterTolerance);
    ds=fminbnd(@profile_loss,bounds(1),bounds(2),outer);
    [ks,ls,ef]=fit_at_size(ds);
    grid=[grid,ds]; kvals=[kvals,ks]; losses=[losses,ls]; flags=[flags,ef];
    rmses=sqrt(losses/L);
    eligible=find(rmses<=min(rmses)+cfg.kRMSEEquivalenceTolerance & flags>0);
    if isempty(eligible), [~,eligible]=min(losses); end
    % When data fits are equivalent, explicitly choose the closest sampled
    % diameter to the supplied midpoint; do not call numerical noise an optimum.
    [~,near]=min(abs(grid(eligible)-D0)); pick=eligible(near);
    k=kvals(:,pick); D=grid(pick); exitflag=flags(pick);
    profile.diameter=grid(:); profile.rmse=rmses(:);
    profile.k=kvals; profile.exitflag=flags(:);
    nearFit=rmses<=min(rmses)+cfg.kRMSEEquivalenceTolerance;
    profile.equivalentSampledRange=[min(grid(nearFit)),max(grid(nearFit))];
    status='single_spectrum_profile_representative';
    warning('kSize:SingleSpectrum', ...
        '%s: one spectrum with free k cannot establish a unique diameter; inspect profile.',model);
else
    % Joint fit of shared log(k) and scaled free diameters. A fixed anchor,
    % if requested, is removed from the optimization variables exactly.
    P=sparse(L*J,numel(free));
    for j=1:numel(free), P((free(j)-1)*L+(1:L),j)=1; end
    options=optimoptions(baseOptions,'JacobPattern',[patternK P]);
    if isempty(free), nstarts=1; else, nstarts=cfg.kNumStarts; end
    best=Inf;
    for s=1:nstarts
        if nstarts==1, t0=0.5; else, t0=0.2+0.6*(s-1)/(nstarts-1); end
        x0=[z0; repmat(t0,numel(free),1)];
        [x,loss,~,flag,out]=lsqnonlin(@joint_residual,x0,...
            [repmat(logLo,L,1);zeros(numel(free),1)],...
            [repmat(logHi,L,1);ones(numel(free),1)],options);
        [kk,dd]=unpack(x);
        starts(s).diameter=dd; starts(s).rmse=sqrt(loss/(L*J));
        starts(s).exitflag=flag; starts(s).output=out;
        if loss<best, best=loss; k=kk; D=dd; exitflag=flag; end
    end
    status='shared_k_bounded_joint_fit';
end

pred=predict(k,D);
result=struct('model',model,'k',k,'diameter',D(:),...
    'inputBounds',originalBounds,'activeBounds',bounds,'reflectance',pred,...
    'residual',pred-R,'rmse',sqrt(mean((pred-R).^2,1)),...
    'totalRMSE',sqrt(mean((pred(:)-R(:)).^2)),...
    'exitflag',exitflag,'status',status,'profile',profile,'starts',starts,...
    'settings',cfg,'n',n,'wavelength',wl,'scaleDiagnostic',struct());
result.kAtBound=(k<=cfg.kLowerBound*1.01 | k>=cfg.kUpperBound/1.01);
result.diameterAtBound=abs(D(:)-bounds(:,1))<=cfg.kDiameterTolerance | ...
    abs(D(:)-bounds(:,2))<=cfg.kDiameterTolerance;
if exitflag<=0, warning('kSize:Convergence','%s solver did not converge.',model); end
if any(result.kAtBound), warning('kSize:KBounds','Some k estimates reached their bounds.'); end
% Scale-family check: k/c and c*D approximately preserve optical depth.
cl=max([bounds(:,1)./D(:); max(k)/cfg.kUpperBound]);
cu=min([bounds(:,2)./D(:); min(k)/cfg.kLowerBound]);
if cu>=cl
    cs=unique([linspace(cl,cu,21),1]); rr=zeros(size(cs));
    for a=1:numel(cs)
        pp=predict(k/cs(a),D(:)*cs(a)); rr(a)=sqrt(mean((pp(:)-R(:)).^2));
    end
    ok=rr<=result.totalRMSE+cfg.kRMSEEquivalenceTolerance;
    result.scaleDiagnostic=struct('factor',cs,'rmse',rr,...
        'equivalentFactorRange',[min(cs(ok)),max(cs(ok))]);
end
result.hasIndependentSizeAnchor=~isempty(cfg.kAnchorIndex) || any(bounds(:,1)==bounds(:,2));
result.uniqueDiameterEstablished=false; % fit convergence is not identifiability

    function [kk,loss,flag]=fit_at_size(dd)
        o=optimoptions(baseOptions,'JacobPattern',patternK);
        [z,loss,~,flag]=lsqnonlin(@(z)fixed_residual(z,dd),z0,...
            repmat(logLo,L,1),repmat(logHi,L,1),o);
        kk=exp(z);
    end
    function loss=profile_loss(dd)
        [~,loss]=fit_at_size(dd);
    end
    function r=fixed_residual(z,dd)
        pp=predict(exp(z),dd); r=pp(:)-R(:);
    end
    function [kk,dd]=unpack(x)
        kk=exp(x(1:L)); dd=D0;
        dd(free)=bounds(free,1)+(bounds(free,2)-bounds(free,1)).*x(L+1:end);
    end
    function r=joint_residual(x)
        [kk,dd]=unpack(x); pp=predict(kk,dd); r=pp(:)-R(:);
    end
    function pp=predict(kk,dd)
        pp=zeros(L,J);
        for b=1:J
            pp(:,b)=local_forward(wl,n,kk,dd(b),model,cfg);
        end
        if ~isreal(pp) || any(~isfinite(pp(:)))
            error('kSize:ForwardDomain','Nonphysical forward spectrum. Check model settings.');
        end
    end
end

function cfg=defaults(cfg)
values=struct('inc',30,'emi',0,'phs',30,'q',0.7,'isSpherical',0,...
    'kLowerBound',1e-9,'kUpperBound',0.5,'kNumStarts',3,...
    'kMaxIterations',1000,'kSolverDisplay','off',...
    'kFunctionTolerance',1e-10,'kStepTolerance',1e-10,...
    'kOptimalityTolerance',1e-10,'kRMSEEquivalenceTolerance',1e-7,...
    'kProfilePoints',7,'kDiameterTolerance',0.01,...
    'kAnchorIndex',[],'kAnchorDiameter',[]);
f=fieldnames(values);
for a=1:numel(f), if ~isfield(cfg,f{a}), cfg.(f{a})=values.(f{a}); end, end
validateattributes(cfg.kLowerBound,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.kUpperBound,{'numeric'},{'scalar','>',cfg.kLowerBound,'finite'});
validateattributes(cfg.kNumStarts,{'numeric'},{'scalar','integer','>=',1});
validateattributes(cfg.kProfilePoints,{'numeric'},{'scalar','integer','>=',3});
validateattributes(cfg.kRMSEEquivalenceTolerance,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.q,{'numeric'},{'scalar','>',0,'<=',1});
if ~ismember(cfg.isSpherical,[0 1]), error('isSpherical must be 0 or 1.'); end
end

function R=local_forward(wl,n,k,D,model,c)
% Matches ROOT Hapke_nk_to_ref and ROOT Shkuratov_forward_model/rb_rf.
% Does not call duplicate functions in MCMC_unmixing_routines.
if c.isSpherical==0, S=0.2*D;
else, S=(2/3)*(n^2-(1/n)*(n^2-1)^(3/2))*D; end
if strcmpi(model,'Hapke')
    se=((n-1)^2+k.^2)./((n+1)^2+k.^2)+0.05;
    si=1.014-4/(n*(n+1)^2);
    theta=exp(-4*pi*n.*k*S./wl);
    w=se+(1-se).*((1-si)*theta)./(1-si*theta);
    mu=cosd(c.emi); mu0=cosd(c.inc);
    P=1-0.4*cosd(c.phs)+0.25*(1.5*cosd(c.phs)^2-0.5);
    r=sqrt(1-w); r0=(1-r)./(1+r);
    H=1./(1-(1-r)*mu.*(r0+(1-0.5*r0-r0*mu)*log((1+mu)/mu)));
    H0=1./(1-(1-r)*mu0.*(r0+(1-0.5*r0-r0*mu0)*log((1+mu0)/mu0)));
    R=w/(4*(mu0+mu)).*(P+H.*H0-1);
else
    Re=((n-1)/(n+1))^2+0.05; Rb=(0.28*n-0.20)*Re;
    Ri=1.04-1/n^2; Te=1-Re; Ti=Te/n^2;
    e=exp(-4*pi*k*S./wl);
    delta=0.5*Te*Ti*Ri*e.^2./(1-Ri*e);
    rb=Rb+delta; rf=(Re-Rb)+Te*Ti*e+delta;
    b=c.q*rb; f=c.q*rf+1-c.q;
    t=(1+b.^2-f.^2)./(2*b);
    if any(t<1-1e-12), error('Invalid Shkuratov albedo domain.'); end
    t=max(t,1); R=1./(t+sqrt(t.^2-1));
end
end
