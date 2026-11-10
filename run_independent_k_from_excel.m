function [calibration,fits,summary] = run_independent_k_from_excel(excelFile,userCfg)
% EACH measured fraction has its OWN fitted k(lambda) and bounded diameter.
% Plot individual k and their unweighted arithmetic mean for each mineral.
% Usage: run_independent_k_from_excel
% Optional: cfg.includeBulkInJointFit=false excludes <74 from fit AND mean.
% That option name is retained for compatibility with the existing loader.
if nargin<1 || isempty(excelFile)
    excelFile='F:\Dataset20260909\Dataset\Raw_Excel\Endmember_Various_Grain_Size.xlsx';
end
if nargin<2, userCfg=struct(); end
here=fileparts(mfilename('fullpath')); addpath(here,'-begin');
cfg=struct('wavelengthUnit','nm','spectrumRangeUm',[1 2.5],...
    'minimumEffectiveDiameter',1,'includeBulkInJointFit',true,...
    'componentNames',{{}},'sheet',1,...
    'knownEndmemberNames',{{'Nontronite','Saponite','MGS-1','Kaolinite','Alunite','Montmorillonite'}},...
    'knownNvals',[1.57 1.51 1.64 1.56 1.58 1.52],...
    'inc',30,'emi',0,'phs',30,'q',0.7,'isSpherical',0,...
    'kNumStarts',3,'kMaxIterations',1000,'kRMSEEquivalenceTolerance',1e-7,...
    'kAnchorIndex',[],'kAnchorDiameter',[],...
    'outputDir',fullfile(here,'Independent_k_Results'),...
    'doFit',true,'makeFigures',true,'figureVisible','on');
fields=fieldnames(userCfg);
for i=1:numel(fields)
    if ~isfield(cfg,fields{i}) && ~startsWith(fields{i},'k')
        error('Unknown configuration field: %s',fields{i});
    end
    cfg.(fields{i})=userCfg.(fields{i});
end
[calibration,mapping]=load_endmember_grain_excel(excelFile,cfg);
calibration.mode='independent_by_fraction';
calibration.meanDefinition='Arithmetic mean over all included measured fractions';
fprintf('\nINDEPENDENT k for EACH size fraction\nSource: %s\n',excelFile);
disp(mapping);
M=numel(calibration.names); L=numel(calibration.wl);
calibration.n=zeros(1,M);
for m=1:M
    ix=find(strcmpi(calibration.names{m},cfg.knownEndmemberNames));
    if numel(ix)~=1, error('No unique refractive index for %s.',calibration.names{m}); end
    calibration.n(m)=cfg.knownNvals(ix);
end
fits=struct('Hapke',{cell(1,M)},'SK',{cell(1,M)}); summary=table();
if ~cfg.doFit, return; end
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
writetable(mapping,fullfile(cfg.outputDir,'Input_column_mapping.csv'));
runFile=fullfile(cfg.outputDir,'independent_k_calibration.mat');
if isfile(runFile)
    copyfile(runFile,[runFile '.before_' datestr(now,'yyyymmdd_HHMMSS') '.bak']);
end
tables=cell(1,M); h_k_all=cell(1,M); s_k_all=cell(1,M);
h_k_mean=nan(L,M); s_k_mean=nan(L,M);
for m=1:M
    name=calibration.names{m}; R=calibration.reflectance{m}; B=calibration.bounds{m};
    fprintf('\n%s: %d independent spectra; n=%g\n',name,size(R,2),calibration.n(m));
    [h,s]=fit_endmember_k_by_fraction(calibration.wl,R,calibration.n(m),B,cfg);
    fits.Hapke{m}=h; fits.SK{m}=s;
    h_k_all{m}=h.k; s_k_all{m}=s.k;
    h_k_mean(:,m)=h.meanK; s_k_mean(:,m)=s.meanK;
    tables{m}=table(repmat(string(name),size(R,2),1),...
        string(calibration.fractionLabels{m}(:)),B(:,1),B(:,2),...
        h.diameter,s.diameter,h.rmse(:),s.rmse(:),...
        h.diameterAtBound,s.diameterAtBound,...
        'VariableNames',{'Mineral','Fraction','Lower_um','Upper_um',...
        'Hapke_D_um','SK_D_um','Hapke_RMSE','SK_RMSE',...
        'Hapke_D_at_bound','SK_D_at_bound'});
    disp(tables{m});
    writetable(tables{m},fullfile(cfg.outputDir,[name '_diameters.csv']));
    labels=matlab.lang.makeValidName(calibration.fractionLabels{m});
    labels=matlab.lang.makeUniqueStrings(labels);
    ht=array2table([calibration.wl,h.k,h.meanK],...
        'VariableNames',[{'Wavelength_um'},strcat('k_',labels(:)'),{'k_mean_all'}]);
    st=array2table([calibration.wl,s.k,s.meanK],...
        'VariableNames',[{'Wavelength_um'},strcat('k_',labels(:)'),{'k_mean_all'}]);
    writetable(ht,fullfile(cfg.outputDir,[name '_Hapke_k_all_and_mean.csv']));
    writetable(st,fullfile(cfg.outputDir,[name '_SK_k_all_and_mean.csv']));
    save(runFile,'calibration','fits','cfg','h_k_all','s_k_all',...
        'h_k_mean','s_k_mean','-v7.3');
    if cfg.makeFigures, plot_material(calibration,m,h,s,cfg); end
end
summary=vertcat(tables{:});
writetable(summary,fullfile(cfg.outputDir,'All_materials_diameters.csv'));
save(runFile,'calibration','fits','cfg','summary','h_k_all','s_k_all',...
    'h_k_mean','s_k_mean','-v7.3');
if cfg.makeFigures, plot_overview(calibration,fits,cfg); end
fprintf('\nCompleted. Independent k and arithmetic means:\n%s\n',runFile);
end

function plot_material(cal,m,h,s,cfg)
wl=cal.wl; R=cal.reflectance{m}; name=cal.names{m}; J=size(R,2);
colors=lines(J);
fig=figure('Name',[name ' independent k'],'Visible',cfg.figureVisible,...
    'Position',[100 100 1500 440]);
subplot(1,3,1); hold on;
for j=1:J, plot(wl,R(:,j),'Color',colors(j,:),'LineWidth',1.1); end
xlabel('Wavelength (um)'); ylabel('Reflectance'); title(name,'Interpreter','none');
legend(cal.fractionLabels{m},'Location','best'); box on;
subplot(1,3,2); plot_k(wl,h,cal.fractionLabels{m},colors); title('Hapke');
subplot(1,3,3); plot_k(wl,s,cal.fractionLabels{m},colors); title('Shkuratov');
exportgraphics(fig,fullfile(cfg.outputDir,[name '_k_by_fraction_and_mean.png']),'Resolution',300);
exportgraphics(fig,fullfile(cfg.outputDir,[name '_k_by_fraction_and_mean.pdf']),'ContentType','vector');
savefig(fig,fullfile(cfg.outputDir,[name '_k_by_fraction_and_mean.fig']));
end

function plot_overview(cal,fits,cfg)
M=numel(cal.names);
fig=figure('Name','All endmembers: independent k and means','Visible',cfg.figureVisible,...
    'Position',[60 60 max(1200,330*M) 720]);
for m=1:M
    colors=lines(size(cal.reflectance{m},2));
    subplot(2,M,m); plot_k(cal.wl,fits.Hapke{m},cal.fractionLabels{m},colors);
    title([cal.names{m} ' - Hapke'],'Interpreter','none');
    subplot(2,M,M+m); plot_k(cal.wl,fits.SK{m},cal.fractionLabels{m},colors);
    title([cal.names{m} ' - Shkuratov'],'Interpreter','none');
end
exportgraphics(fig,fullfile(cfg.outputDir,'All_endmembers_k_and_mean.png'),'Resolution',300);
exportgraphics(fig,fullfile(cfg.outputDir,'All_endmembers_k_and_mean.pdf'),'ContentType','vector');
savefig(fig,fullfile(cfg.outputDir,'All_endmembers_k_and_mean.fig'));
end

function plot_k(wl,f,labels,colors)
hold on; J=size(f.k,2); hh=gobjects(J+1,1);
for j=1:J
    hh(j)=plot(wl,f.k(:,j),'-','Color',colors(j,:),'LineWidth',1.1);
end
hh(end)=plot(wl,f.meanK,'k--','LineWidth',2.2);
set(gca,'YScale','log'); xlabel('Wavelength (um)'); ylabel('k'); box on;
legend(hh,[labels(:);{'Mean of all fractions'}],'Location','best','Interpreter','none');
end
