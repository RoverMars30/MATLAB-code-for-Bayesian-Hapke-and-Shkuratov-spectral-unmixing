function [calibration,mapping] = load_endmember_grain_excel(excelFile,cfg)
% Read the user's raw size-fraction workbook WITHOUT fixed column numbers.
% Header examples: MGS<43, MGS43<61, MGS61<74, MGS-1<74, NAu-2<74.
% Full-width less-than signs are normalized. Wavelength units are explicit.
% No spectral interpolation, reflectance rescaling, or duplicate averaging.
if nargin<2, cfg=struct(); end
if ~isfield(cfg,'wavelengthUnit'), cfg.wavelengthUnit='nm'; end
if ~isfield(cfg,'spectrumRangeUm'), cfg.spectrumRangeUm=[1 2.5]; end
if ~isfield(cfg,'minimumEffectiveDiameter'), cfg.minimumEffectiveDiameter=1; end
if ~isfield(cfg,'includeBulkInJointFit'), cfg.includeBulkInJointFit=true; end
if ~isfield(cfg,'componentNames'), cfg.componentNames={}; end
if ~isfield(cfg,'sheet'), cfg.sheet=1; end
validateattributes(cfg.minimumEffectiveDiameter,{'numeric'},{'scalar','finite','positive'});
if ~isfile(excelFile), error('Input Excel not found: %s',excelFile); end
C=readcell(excelFile,'Sheet',cfg.sheet);
if size(C,1)<4 || size(C,2)<2, error('Excel needs a header and numeric spectra.'); end
headers=string(C(1,:));
% Column 1 in this source is explicitly Wavelength. Do not guess a spectrum
% column as the wavelength axis when the file layout is unexpected.
if ~contains(lower(headers(1)),'wave') && ~contains(headers(1),'波长')
    error('First column must be Wavelength. Found: %s',headers(1));
end
A=nan(size(C,1)-1,size(C,2));
for j=1:size(C,2)
    for i=2:size(C,1)
        v=C{i,j};
        if isnumeric(v) && isscalar(v), A(i-1,j)=double(v);
        elseif ischar(v) || (isstring(v) && isscalar(v))
            A(i-1,j)=str2double(v);
        end
    end
end
numericRows=isfinite(A(:,1));
sourceRows=find(numericRows)+1; A=A(numericRows,:);
switch lower(string(cfg.wavelengthUnit))
    case "nm", wl=A(:,1)/1000;
    case {"um","micron","micrometre"}, wl=A(:,1);
    otherwise, error('Set wavelengthUnit explicitly to nm or um.');
end
if any(diff(wl)<=0) || any(wl<=0), error('Wavelengths must increase without duplicates.'); end
use=true(size(wl));
if ~isempty(cfg.spectrumRangeUm)
    validateattributes(cfg.spectrumRangeUm,{'numeric'},{'numel',2,'finite','positive'});
    if cfg.spectrumRangeUm(1)>=cfg.spectrumRangeUm(2), error('Reversed spectral interval.'); end
    use=wl>=cfg.spectrumRangeUm(1) & wl<=cfg.spectrumRangeUm(2);
end
if nnz(use)<3, error('Fewer than three wavelength channels selected.'); end
R=A(use,:); wl=wl(use); sourceRows=sourceRows(use);
count=size(C,2)-1;
mineral=strings(count,1); labels=mineral;
lo=zeros(count,1); hi=lo; isBulk=false(count,1);
for j=2:size(C,2)
    [mineral(j-1),lo(j-1),hi(j-1),labels(j-1),isBulk(j-1)]=...
        parse_header(headers(j),cfg.minimumEffectiveDiameter);
    bad=find(~isfinite(R(:,j)) | R(:,j)<=0 | R(:,j)>=1,1);
    if ~isempty(bad)
        error('Invalid reflectance in header %s, Excel row %d. Expected 0<R<1.',...
            headers(j),sourceRows(bad));
    end
end
selected=true(count,1);
if ~isempty(cfg.componentNames)
    selected=ismember(lower(mineral),lower(string(cfg.componentNames)));
    missingNames=string(cfg.componentNames);
    missingNames=missingNames(~ismember(lower(missingNames),lower(mineral)));
    if ~isempty(missingNames), error('Missing requested material: %s',strjoin(missingNames,', ')); end
end
if ~cfg.includeBulkInJointFit, selected=selected & ~isBulk; end
mapping=table((2:size(C,2))',headers(2:end)',mineral,labels,lo,hi,isBulk,selected,...
    'VariableNames',{'ExcelColumn','SourceHeader','Mineral','Fraction',...
    'Lower_um','Upper_um','BulkFraction','Selected'});
names=cellstr(unique(mineral(selected),'stable'))';
if isempty(names), error('No fractions selected for calibration.'); end
calibration=struct('wl',wl,'names',{names},'reflectance',{cell(1,numel(names))},...
    'bounds',{cell(1,numel(names))},'fractionLabels',{cell(1,numel(names))},...
    'sourceHeaders',{cell(1,numel(names))},'sourceColumns',{cell(1,numel(names))},...
    'sourceExcel',excelFile,'sourceSheet',cfg.sheet,'sourceRows',sourceRows,...
    'inputWavelengthUnit',cfg.wavelengthUnit,'mapping',mapping);
for m=1:numel(names)
    ix=find(selected & mineral==string(names{m}));
    [~,order]=sortrows([hi(ix),lo(ix)],[1 2]); ix=ix(order);
    calibration.reflectance{m}=R(:,ix+1);
    calibration.bounds{m}=[lo(ix),hi(ix)];
    calibration.fractionLabels{m}=cellstr(labels(ix));
    calibration.sourceHeaders{m}=cellstr(headers(ix+1));
    calibration.sourceColumns{m}=ix+1;
end
end

function [mineral,lb,ub,label,bulk]=parse_header(raw,floorD)
txt=lower(regexprep(char(raw),'\s+',''));
txt=strrep(txt,char(65308),'<'); % full-width less-than, present in source
txt=strrep(txt,char(8211),'-'); txt=strrep(txt,char(8212),'-');
patterns={'^mgs(?:-?1)?','^sap(?:onite)?','^alu(?:nite)?',...
    '^kao(?:linite)?','^mont(?:morillonite)?','^(?:nau-?2|nontronite)'};
names={'MGS-1','Saponite','Alunite','Kaolinite','Montmorillonite','Nontronite'};
tail=''; mineral="";
for a=1:numel(patterns)
    [~,last]=regexp(txt,patterns{a},'once');
    if ~isempty(last), mineral=string(names{a}); tail=txt(last+1:end); break; end
end
if strlength(mineral)==0, error('Unknown mineral in header: %s',raw); end
% Remove units, retaining the actual numeric bounds rather than a fixed map.
tail=regexprep(tail,'(?:microns?|um|μm|µm)$','');
tokens=regexp(tail,'\d+(?:\.\d+)?','match'); nums=str2double(tokens);
if numel(nums)==2
    lb=nums(1); ub=nums(2); label=sprintf('%g-%g um',lb,ub); bulk=false;
elseif numel(nums)==1 && startsWith(tail,'<')
    lb=floorD; ub=nums(1); label=sprintf('<%g um',ub); bulk=(ub==74);
else
    error('Ambiguous size header: %s. Use e.g. MGS<43 or MGS43-61.',raw);
end
if lb<=0 || ub<lb, error('Invalid positive size bounds in header: %s',raw); end
end
