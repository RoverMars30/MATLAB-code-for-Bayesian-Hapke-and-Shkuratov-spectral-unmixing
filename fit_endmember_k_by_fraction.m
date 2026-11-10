function [h,s] = fit_endmember_k_by_fraction(wl,R,n,bounds,cfg)
% ONE material, each size fraction fitted independently; no shared k.
% Same input interface as invert_endmember_k_size. Each output .k is L-by-J;
% .meanK is arithmetic mean over ALL J included measured spectra at each wl.
% Full per-fraction solver/profile diagnostics are in .fractionResults{j}.
if nargin<5, cfg=struct(); end
if isvector(R), R=R(:); end
[L,J]=size(R);
if ~isequal(size(bounds),[J 2]), error('Need one bounds row per fraction.'); end
h=empty_result(L,J,'Hapke'); s=empty_result(L,J,'SK');
for j=1:J
    fprintf('  Fraction %d/%d, D in [%g,%g] um\n',j,J,bounds(j,1),bounds(j,2));
    hj=invert_endmember_k_size(wl,R(:,j),n,bounds(j,:),'Hapke',cfg);
    sj=invert_endmember_k_size(wl,R(:,j),n,bounds(j,:),'SK',cfg);
    h=store(h,hj,j); s=store(s,sj,j);
end
h.meanK=mean(h.k,2); s.meanK=mean(s.k,2);
h.meanDefinition='Arithmetic mean of all included size-fraction k spectra';
s.meanDefinition=h.meanDefinition;
h.inputBounds=bounds; s.inputBounds=bounds;
h.wavelength=wl(:); s.wavelength=wl(:);
% meanK is a summary of independently fitted inputs. It was not itself fitted.
end

function r=empty_result(L,J,model)
r=struct('model',model,'mode','independent_by_fraction',...
    'k',nan(L,J),'meanK',nan(L,1),'diameter',nan(J,1),...
    'reflectance',nan(L,J),'residual',nan(L,J),'rmse',nan(1,J),...
    'exitflag',nan(J,1),'diameterAtBound',false(J,1),...
    'fractionResults',{cell(1,J)});
end

function r=store(r,x,j)
r.k(:,j)=x.k; r.diameter(j)=x.diameter;
r.reflectance(:,j)=x.reflectance; r.residual(:,j)=x.residual;
r.rmse(j)=x.totalRMSE; r.exitflag(j)=x.exitflag;
r.diameterAtBound(j)=x.diameterAtBound;
r.fractionResults{j}=x;
end
