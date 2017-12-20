function [ EEG ] = itab_makeAnalysisStr( fullfilepath, fulllocsfile, saveEEGLABstr)
%ITAB_MAKEEEG ECTSTR Make a single EEG ect analysis structure.
%   Inputs args: 
%                'parfullfilepath' : eventually file with parameters for
%                                analysis;
%                'fullfullfolderpath': path of the folder with all
%                                        subjects
%                'locsfile': channels location file, default EEGLAB locs;
%                'saveEEGLABstr': = 1 if you want to crate an EEGLAB
%                                structure;
%                'condition': conditin name for current subject;
%   check in args
    if nargin < 1
        error('ITAB_MAKEEEG ECTSTR: fullfilepath and analysis folder path are required input arg');
    end
    if nargin < 2
        fulllocsfile = uigetfile('*','Provide the fullpath of the locations file');
    end
    if nargin < 3
        saveEEGLABstr = 0;
    end
    [fullfolderpath, name, ~] = fileparts(fullfilepath);
    condition = name(7:8);
    fulparfile = fullfile(fullfolderpath, strcat(name, '.par'));
    fullmhfile = strcat(fullfilepath,'.mhd');
    [h_filtering, l_filtering, freq_notch, dec, ~, ~, ~, skipped_intervals, ~, ~, ~, ~, ~, ~, ~, ~] = itab_readparfile(fulparfile);
    excludedchs = [48 119 125:128];
    header = itab_lcreadheader(fullmhfile); 
    if header.raw_header_type == header.header_type
        warning('warning! master header not loaded...'); 
    end
    gch = itab_setchannels(header,'ELEC','working');
    for c = 1: numel(excludedchs)
        gch(gch == excludedchs(c)) = [];
    end
    bchs = 1:128;
    bchs(gch) = [];
    for c = 1: numel(excludedchs)
        bchs(bchs == excludedchs(c)) = [];
    end
    nch = length(gch);   
    lab_channels = cell(nch,1);
    for i=1:length(gch)
        lab_channels{i} = header.ch(gch(i)).label;
    end

    [data, newsmplf, newntpdata] = itab_readedf(fullfilepath, h_filtering, l_filtering, freq_notch, dec, gch);
    EEG = struct;
    if saveEEGLABstr
        EEGOUT = pop_importdata( 'data', data, ...
                                 'setname', fullfolderpath,...
                                 'srate', newsmplf, ...
                                 'pnts', newntpdata, ...
                                 'nbchan', length(gch), ...
                                 'subject', name(1:6), ...
                                 'session', name(end-3:end), ...
                                 'chanlocs', fulllocsfile, ...
                                 'condition', condition);
        EEG = EEGOUT;
    else
        EEG.data = data;
        EEG.setname = fullfolderpath;
        EEG.filename = fullfilepath;
        EEG.srate = newsmplf;
        EEG.pnts = newntpdata;
        EEG.nbchan = length(gch);
        EEG.subject = name(1:6);
        EEG.session = name(end-3:end);
        EEG.chanlocs = fulllocsfile;
        EEG.parfile = fulparfile;
        EEG.skipped = skipped_intervals;
        EEG.chlabels = lab_channels;
        EEG.badchannels = bchs;
        EEG.condition = condition;
        EEG.decimationfactor = dec;
        EEG.hfiltfreq = h_filtering;
        EEG.lfiltfreq = l_filtering;
    end 
end

