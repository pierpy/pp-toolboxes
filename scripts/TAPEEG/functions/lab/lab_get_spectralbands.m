function spectralbands = lab_get_spectralbands(spectralbands,spectralbandsI,doedit)

global Main_Path

if ~exist('doedit','var')
    doedit = false;
end
if ~exist('spectralbandsI','var') | ~isnumeric(spectralbands)
    spectralbandsI = false;
end
if exist('spectralbands','var') & ~isempty(spectralbands) & isnumeric(spectralbands) & isnan(spectralbands)
    return
end

if ~exist('spectralbands','var') | isempty(spectralbands)
    if spectralbandsI == false
        if exist(fullfile(fullfile(Main_Path,'.ignore'),'spectras.mat'),'file')
            load(fullfile(fullfile(Main_Path,'.ignore'),'spectras.mat'))
            if isempty(spectralbands)
                spectralbands = {'Delta',1,4; ...
                    'Theta',4,8; ...
                    'Alpha1',8,10; ...
                    'Alpha2',10,13; ...
                    'Beta',13,30};
            end
        else
            spectralbands = {'Delta',1,4; ...
                'Theta',4,8; ...
                'Alpha1',8,10; ...
                'Alpha2',10,13; ...
                'Beta',13,30};
        end
    else
        spectralbands = [1,2,4,5,6,8,9];
    end
else
    
end
if doedit == true
    if spectralbandsI == false
        spectralbands = lab_table_dialog(spectralbands,{'Band','lowfreq','highfreq'},'Spectral Bands',1);
        save(fullfile(fullfile(Main_Path,'.ignore'),'spectras.mat'),'spectralbands');
    else
        List = {'Delta','Theta','Alpha1','Alpha1a','Alpha1b','Alpha2','Beta','Beta1','Beta2','PF'};
        spectralbands = listdlg('PromptString','Frequency bands:','SelectionMode','multiple', ...
            'ListString',List,'InitialValue',spectralbands,'CancelString','None','ListSize',[100 170]);
    end
end