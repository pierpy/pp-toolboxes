function [cfg,Prompt,Formats] = lab_set_postprocessing(cfg,header,data,Prompt,Formats)

global THeader TData

dodiag = 1;

if ~exist('Formats','var')
    Formats = [];
end
if ~exist('Prompt','var')
    Prompt = cell(0,2);
else
    dodiag = 0;
end
if isempty(TData) & exist('data','var')
    TData = data;
end
if isempty(THeader) & exist('header','var')
    THeader = header;
end
if ~exist('cfg','var')
    cfg = [];
end

Prompt(end+1,:) = {'Postprocessing',''};
Formats(end+1,1).type = 'text';
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 4]; % item is 1 field x 4 fields

Prompt(end+1,:) = {'Write result','RES'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_write_result,'@ALL','@ALL'};

Prompt(end+1,:) = {'Export EDF','EDF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_export2edf,'@ALL','@ALL'};

Prompt(end+1,:) = {'Re-reference','REF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_reference_data,'@ALL','@ALL'};

Prompt(end+1,:) = {'Save epochs','EPOCH'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_save_epochs,'@ALL','@ALL'};

Prompt(end+1,:) = {'Average','AVG'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_average,'@ALL','@ALL'};

Prompt(end+1,:) = {'Inverse solution','IS'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_inversesolution,'@ALL','@ALL'};

Prompt(end+1,:) = {'Inverse solution FFT','ISFFT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_inversesolution_fft,'@ALL','@ALL'};

Prompt(end+1,:) = {'Frequency analysis','FFT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_calculate_spectras,'@ALL','@ALL'};

Prompt(end+1,:) = {'Microstates','MICROST'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_microstates,'@ALL','@ALL'};

Prompt(end+1,:) = {'Connectivity','CONNECT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_calculate_connectivity,'@ALL','@ALL'};

Prompt(end+1,:) = {'Distance matrix','DISTANCE'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Forward solution','FWS'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_forwardsolution,'@ALL','@ALL'};

Prompt(end+1,:) = {'Coregistration LOCS/MRI','COREG'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_coreg,'COREG','COREG',[],1,0};

if dodiag == 1
    [cfg,Cancelled] = inputsdlg(Prompt,'Postprocessing',Formats,cfg);
    if isempty(cfg) | Cancelled == 1
        cfg = [];
        Prompt = 1;
        return
    else
        Prompt = 0;
    end

    TData = [];
    THeader = [];
end

end

