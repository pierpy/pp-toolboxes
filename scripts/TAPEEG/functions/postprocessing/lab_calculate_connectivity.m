% Calculate PLI, dPLI, wPLI, PLV, PLT, SL, AEC and CAE in specified frequency ranges
%
% data                = Matrix channels x timepoints
% header.samplingrate = Samplingrate
% cfg.CONNECT.step    = Size of steps in phases of max frequency
% cfg.CONNECT.freqs   = Matrix with frequencybands to analyze (low1 high1;low2 high2;...)
% cfg.CONNECT.window  = Min window for wavelet analysis in seconds
%
% written by F.Hatz 2013 Vumc Neurophysiology Amsterdam

function [Result,cfg] = lab_calculate_connectivity(data,header,cfg)

Result = [];
if ~exist('data','var')
    cfg =[];
    return
end

if ~exist('cfg','var') | ~isfield(cfg,'Output_file')
    cfg.Output_filepath = tempdir;
    cfg.Output_file = [tempname '.sef'];
    UseTemp = true;
else
    UseTemp = false;
end
if ~isempty(cfg.Output_file)
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
end

if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'CONNECT') & cfg.SKIP.CONNECT == true;
    return
end

disp('Calculation of Connectivity')
if ~exist('cfg','var') | ~isfield(cfg,'CONNECT') | isempty(cfg.CONNECT)
    [cfg,skipprocessing] = lab_set_calculate_connectivity(cfg,header,data);
    if skipprocessing == 1
        return
    else
        pause(0.2);
    end
end

if size(data,2) > 5000000
    disp('   Warning: Connectivity analysis may fail, data-file to large (>5''000''000 TF)')
end

if UseTemp == false
    Output_filepathM = cfg.Output_filepath;
    if ~isfield(cfg.CONNECT,'folder')
        cfg.CONNECT.folder = 'Connectivity';
    end
    cfg.Output_filepath = fullfile(cfg.Output_filepath,cfg.CONNECT.folder);
    warning off %#ok<WNOFF>
    mkdir(cfg.Output_filepath);
    warning on %#ok<WNON>
end

% Calculate new reference
if strcmp(cfg.CONNECT.eegsource,'montage') & isfield(cfg.CONNECT,'montage') & ~isempty(cfg.CONNECT.montage)
    [data,header,cfg.CONNECT] = lab_references(data,header,cfg.CONNECT.montage,cfg.CONNECT);
elseif ~strcmp(cfg.CONNECT.eegsource,'input')
    [data,header,cfg.CONNECT] = lab_references(data,header,cfg.CONNECT.eegsource,cfg.CONNECT);
end

% Interpolate bad channels
if isfield(cfg.CONNECT,'interpolate') & cfg.CONNECT.interpolate == true
    disp('   Interpolate bad channels')
    [data,header] = lab_interpolate_bad(data,header);
end

% Reduce to data channels
if isfield(header,'numdatachannels')
    [data,header] = lab_reduce_channels(data,header,1:header.numdatachannels);
end

if isfield(cfg.CONNECT,'PHASE') & isfield(cfg.CONNECT.PHASE,'phaseestimate')
    if cfg.CONNECT.PHASE.window == 0
        cfg.CONNECT.PHASE.freqwindow = size(data,2);
    else
        cfg.CONNECT.PHASE.freqwindow = floor(header.samplingrate*cfg.CONNECT.PHASE.window);
    end
end

if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerinclude_combine') & ...
        cfg.CONNECT.MARKER.markerinclude_combine == false & ~isempty(cfg.CONNECT.MARKER.markerinclude)
    MaxRun = length(cfg.CONNECT.MARKER.markerinclude);
else
    MaxRun = 1;
end

if isfield(cfg.CONNECT,'spectralbandsI') & cfg.CONNECT.spectralbandsI == true
    if ~isfield(header,'IFREQ') | ~isfield(header.IFREQ,'Bands') | isempty(header.IFREQ.Bands)
        [~,header] = lab_indiv_freqbands(data,header);
    end
    if ~isempty(header.IFREQ.Bands)
        cfg.CONNECT.freqs = cell2mat(header.IFREQ.Bands(:,2:5));
        if isnumeric(cfg.CONNECT.spectralbands) & ~isempty(cfg.CONNECT.spectralbands)
            cfg.CONNECT.freqs = cfg.CONNECT.freqs(cfg.CONNECT.spectralbands,:);
        end
    else
        disp('   Abort: Calculation of individual bands not possible')
        return
    end
elseif ~isempty(cfg.CONNECT.spectralbands)
    cfg.CONNECT.freqs = cell2mat(cfg.CONNECT.spectralbands(:,2:3));
    cfg.CONNECT.freqs = [cfg.CONNECT.freqs cfg.CONNECT.freqs];
else
    cfg.CONNECT.freqs = 0;
end

Result = [];
if ~isfield(cfg,'patient')
    cfg.patient = 'Pat';
end
if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerinclude')
    MarkerInclude = cfg.CONNECT.MARKER.markerinclude;
else
    MarkerInclude = {};
end
if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath)
    Output_fileS = cfg.Output_fileS;
    Output_file = cfg.Output_file;
    Output_filepath = cfg.Output_filepath;
else
    Output_fileS = '';
    Output_file = '';
end
patient = cfg.patient;
for i = 1:MaxRun
    if MaxRun > 1
        cfg.CONNECT.MARKER.markerinclude = MarkerInclude(i);
        cfg.Output_fileS = [Output_fileS '_' MarkerInclude{i}];
        cfg.Output_file = [Output_fileS '_' MarkerInclude{i} '.sef'];
        disp(['   Calculate for Marker: ' MarkerInclude{i}])
        cfg.patient = [patient '_' MarkerInclude{i}];
    end
    for j = 1:length(cfg.CONNECT.measure)
        skipprocessing = 0;
        measure = cfg.CONNECT.measure{j};
        switch measure
            case {'PLI','dPLI','wPLI','PLT','PLV','wPLV','PTE','rPTE'}
                if ~exist('phasedone','var')
                    if strcmp(cfg.CONNECT.PHASE.phaseestimate,'wavelet')
                        [Result,cfg] = do_calc(data,header,cfg,@lab_calculate_connect_wavelet);
                    elseif strcmp(cfg.CONNECT.PHASE.phaseestimate,'rihaczek')
                        [Result,cfg] = do_calc(data,header,cfg,@lab_calculate_connect_tfr);
                    elseif strcmp(cfg.CONNECT.PHASE.phaseestimate,'hilbert') | ...
                            strcmp(cfg.CONNECT.PHASE.phaseestimate,'zerocross') | ...
                            strcmp(cfg.CONNECT.PHASE.phaseestimate,'realphase')
                        [Result,cfg,header] = do_calc(data,header,cfg,@lab_calculate_connect_phase);
                    else
                        disp('    Abort: no valid phase estimation selected')
                        return
                    end
                    phasedone = true; %#ok<NASGU>
                else
                    skipprocessing = 1;
                end
            case 'SL'
                [Result,cfg] = do_calc(data,header,cfg,@lab_calculate_connect_SL);
                clearvars Resulttmp
            case 'SLc'
                [Result,cfg] = do_calc(data,header,cfg,@lab_calculate_connect_SLc);
                clearvars Resulttmp
            case 'EC'
                [Result,cfg] = do_calc(data,header,cfg,@lab_calculate_connect_EC);
                clearvars Resulttmp
            case 'DTF'
                [Result,cfg] = do_calc(data,header,cfg,@lab_calculate_connect_dtf);
                clearvars Resulttmp
            case 'TE'
                [Result,cfg] = do_calc(data,header,cfg,@lab_calculate_connect_TE);
                clearvars Resulttmp
        end
        if isfield(cfg.CONNECT,'doaverage') & cfg.CONNECT.doaverage == true & isfield(cfg,'lastfilefolder') & cfg.lastfilefolder == 0
            skipprocessing = 1;
        end
        if skipprocessing == 0
            if UseTemp == true
                delete(fullfile(cfg.Output_filepath,'Conn_*'));
                cfg = rmfield(cfg,'Output_filepath');
                cfg = rmfield(cfg,'Output_file');
                cfg = rmfield(cfg,'Output_fileS');
                if isfield(cfg,'listold')
                    cfg = rmfield(cfg,'listold');
                end
            end
            
            % correct distances
            if isfield(cfg.CONNECT,'correctdistance') & cfg.CONNECT.correctdistance > 0 & ...
                    isfield(header,'locs') & isfield(header.locs,'x') & ...
                    size(header.locs.x,2) == header.numdatachannels
                Result = lab_calculate_connect_correctdistance(Result,header,cfg);
            end
            
            % write matrixes
            if isfield(cfg,'CONNECT') & isfield(cfg.CONNECT,'writematrix') & cfg.CONNECT.writematrix == true
                cfg = lab_calculate_connect_writematrix(Result,cfg);
            end
            
            % Plot matrices
            if isfield(cfg,'CONNECT') & isfield(cfg.CONNECT,'PLOT') & ~isempty(cfg.CONNECT.PLOT)
                cfg = lab_calculate_connect_plotmatrix(Result,cfg);
            end
            
            % Calculate graph measures
            if isfield(cfg.CONNECT,'GRAPH') & isfield(cfg.CONNECT.GRAPH,'measure') & ~isempty(cfg.CONNECT.GRAPH.measure)
                [Result,cfg] = lab_calculate_connect_dograph(Result,header,cfg);
            end
            
            % calculate kmeans
            if isfield(cfg.CONNECT,'clustering') & cfg.CONNECT.clustering > 1
                Result = lab_calculate_connect_kmeans(Result,cfg);
            end
            
            %--------------------------------------------------------------------------
            % Write verbose file (*.vrb)
            %--------------------------------------------------------------------------
            if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath)
                if ~isfield(cfg.CONNECT,'spectralbandsI') | cfg.CONNECT.spectralbandsI == false
                    Freqs = [];
                    for nfreq = 1:size(cfg.CONNECT.spectralbands,1);
                        Freqs = [Freqs 'F' num2str(cfg.CONNECT.spectralbands{nfreq,2}) 'F' num2str(cfg.CONNECT.spectralbands{nfreq,3})]; %#ok<AGROW>
                        if nfreq ~= size(cfg.CONNECT.spectralbands,1)
                            Freqs = [Freqs ' / ']; %#ok<AGROW>
                        end
                    end
                else
                    Freqs = 'Individual Bands';
                end
                fid=fopen(fullfile(cfg.Output_filepath,[cfg.Output_fileS '.vrb']),'w');
                fprintf(fid,'Connectivity analysis\n');
                fprintf(fid,['File: ' cfg.Output_fileS]);
                fprintf(fid,'\n');
                fprintf(fid,['Spectral bands: ' Freqs '\n']);
                fprintf(fid,'\n');
                fprintf(fid,['Filter: ' cfg.CONNECT.filter '\n']);
                fprintf(fid,'\n');
                if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerexclude') & ...
                        ~isempty(cfg.CONNECT.MARKER.markerexclude) & iscell(cfg.CONNECT.MARKER.markerexclude)
                    fprintf(fid,['  Marker excluded: ' sprintf('%s|',cfg.CONNECT.MARKER.markerexclude{:}) '\n']);
                end
                if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerinclude') & ...
                        ~isempty(cfg.CONNECT.MARKER.markerinclude) & iscell(cfg.CONNECT.MARKER.markerinclude)
                    fprintf(fid,['  Marker included: ' sprintf('%s|',cfg.CONNECT.MARKER.markerinclude{:}) '\n']);
                end
                if exist('phasedone','var')
                    Measures = intersect({'PLI','dPLI','wPLI','PLT','PLV','wPLV','PTE','rPTE'},cfg.CONNECT.measure);
                    fprintf(fid,'Measures: ');
                    fprintf(fid,'%s|',Measures{:});
                    fprintf(fid,'\n');
                    if isfield(cfg.CONNECT,'PHASE') & ~isempty(cfg.CONNECT.PHASE)
                        fprintf(fid,['Phase estimation: ' num2str(cfg.CONNECT.PHASE.phaseestimate) '\n']);
                        fprintf(fid,['  Window size: ' num2str(cfg.CONNECT.PHASE.window) '\n']);
                        fprintf(fid,['  Use hanning window: ' num2str(cfg.CONNECT.PHASE.usehann) '\n']);
                        if isfield(cfg.CONNECT.PHASE,'skipstart')
                            fprintf(fid,['  Skip start: ' num2str(cfg.CONNECT.PHASE.skipstart) ' TFs\n']);
                            fprintf(fid,['  Skip end: ' num2str(cfg.CONNECT.PHASE.skipend) ' TFs\n']);
                        end
                        fprintf(fid,['Step: ' num2str(cfg.CONNECT.PHASE.step) ' ' cfg.CONNECT.PHASE.stepunit '\n']);
                        fprintf(fid,'\n');
                    end
                else
                    fprintf(fid,['Measure: ' measure]);
                    fprintf(fid,'\n');
                    if strcmp(measure,'SL') & isfield(cfg.CONNECT,'SL') & ~isempty(cfg.CONNECT.SL)
                        fprintf(fid,['SL pref: ' num2str(cfg.CONNECT.SL.pref) '\n']);
                        fprintf(fid,'\n');
                    end
                    if strcmp(measure,'SLc') & isfield(cfg.CONNECT,'SLc') & ~isempty(cfg.CONNECT.SLc)
                        fprintf(fid,['SLc pref: ' num2str(cfg.CONNECT.SLc.pref) '\n']);
                        tmp = {'uni-directional symmetrical','uni-directional asymmetrical','bi-directional'};
                        fprintf(fid,['Correct for volume conduction: ' tmp{cfg.CONNECT.SLc.Mcorrect} '\n']);
                        clearvars tmp
                        fprintf(fid,'\n');
                    end
                    if strcmp(measure,'EC') & isfield(cfg.CONNECT,'EC') & ~isempty(cfg.CONNECT.EC)
                        fprintf(fid,['Envelope correlation step: ' num2str(cfg.CONNECT.EC.stepec) '\n']);
                        fprintf(fid,['Use hilbert transform: ' num2str(cfg.CONNECT.EC.EChilbert) '\n']);
                        fprintf(fid,'\n');
                    end
                    if strcmp(measure,'DTF') & isfield(cfg.CONNECT,'DTF') & ~isempty(cfg.CONNECT.DTF)
                        fprintf(fid,['  Lag: ' num2str(cfg.CONNECT.DTF.lag) '\n']);
                        fprintf(fid,['  Length: ' num2str(cfg.CONNECT.DTF.length) '\n']);
                        fprintf(fid,['  Order: ' num2str(cfg.CONNECT.DTF.order) '\n']);
                        if cfg.CONNECT.DTF.sigtest == true
                            fprintf(fid,['  Input significance test (Shuffle times: ' num2str(cfg.CONNECT.DTF.shufftimes) ...
                                ', Significance level: ' num2str(cfg.CONNECT.DTF.siglevel) '\n']);
                        end
                        fprintf(fid,'\n');
                    end
                end
                if isfield(cfg.CONNECT,'correctdistance') & ~isempty(cfg.CONNECT.correctdistance)
                    fprintf(fid,['Correct distance: ' num2str(cfg.CONNECT.correctdistance) '\n']);
                end
                if isfield(cfg.CONNECT,'clustering') & ~isempty(cfg.CONNECT.clustering)
                    fprintf(fid,['Phase-Kmeans clustering (number): ' num2str(cfg.CONNECT.clustering) '\n']);
                end
                fprintf(fid,'\n');
                if isfield(cfg.CONNECT,'GRAPH') & isfield(cfg.CONNECT.GRAPH,'measure')
                    fprintf(fid,'\n');
                    fprintf(fid,'Graph analysis:\n  ');
                    fprintf(fid,'%s|',cfg.CONNECT.measure{:});
                    fprintf(fid,'\n');
                    if cfg.CONNECT.GRAPH.rankmatrix == true;
                        fprintf(fid,['  Rank matrix with order ' num2str(cfg.CONNECT.GRAPH.rankorder) '\n']);
                        fprintf(fid,'\n');
                    end
                    fprintf(fid,['  Max number of matrices analysed per patient: ' num2str(cfg.CONNECT.GRAPH.maxmatrix) '\n']);
                    if cfg.CONNECT.GRAPH.avgmatrix == true
                        if isempty(cfg.CONNECT.GRAPH.avgmatrix)
                            fprintf(fid,'  Number of matrices included in average matrix per patient: all\n');
                        else
                            fprintf(fid,['  Number of matrices included in average matrix per patient: ' num2str(cfg.CONNECT.GRAPH.numavg) '\n']);
                        end
                    end
                    fprintf(fid,['  Number of randomisations for random graphs: ' num2str(cfg.CONNECT.GRAPH.randnumber) '\n']);
                    fprintf(fid,'\n');
                end
                fclose(fid);
            end
            if exist('Output_filepath','var')
                cfg.Output_filepath = Output_filepath;
            end
        end
    end
    clearvars phasedone
end
if exist('Output_filepathM','var')
    cfg.Output_filepath = Output_filepathM;
end
if ~isempty(MarkerInclude)
    cfg.CONNECT.MARKER.markerinclude = MarkerInclude;
end
if ~isempty(Output_fileS)
    cfg.Output_fileS = Output_fileS;
    cfg.Output_file = Output_file;
end
cfg.patient = patient;

end

function [Result,cfg,header] = do_calc(data,header,cfg,fc)
    if isfield(cfg,'CONNECT') & isfield(cfg.CONNECT,'randphase') & cfg.CONNECT.randphase == true & ...
            isfield(cfg.CONNECT,'numrandphase') & ~isempty(cfg.CONNECT.numrandphase)
        % calculate result of original data
        [Result,cfg,header] = fc(data,header,cfg);
        if isfield(cfg.CONNECT,'doaverage') & cfg.CONNECT.doaverage == true
            Result = average_matrices(Result);
        end
        
        % calculate results of data with randomized phases
        disp(['   Create data with random phase (' num2str(cfg.CONNECT.numrandphase) ')'])
        cfg2 = cfg;
        cfg2.CONNECT.PHASE.storephase = false;
        if isfield(cfg,'patient') & ~isempty(cfg.patient)
            cfg2.patient = [cfg.patient '_RandPhase'];
        else
            cfg2.patient = 'RandPhase';
        end
        randdata = phaseran(data',cfg.CONNECT.numrandphase);
        randdata = permute(randdata,[2 1 3]);
        if size(randdata,2) < size(data,2);
            randdata(1,size(data,2),1) = 0;
        end
        for i = 1:size(randdata,3)
            [~,cfg2] = fc(randdata(:,:,i),header,cfg2);
        end
        if isfield(cfg,'listold') & isfield(cfg2,'listold')
            cfg.listold = union(cfg.listold,cfg2.listold);
        elseif isfield(cfg2,'listold')
            cfg.listold = cfg2.listold;
        end
        if size(cfg.listold,1) > 1
            cfg.listold = cfg.listold';
        end
    else
        [Result,cfg,header] = fc(data,header,cfg);
    end
end

function Result = average_matrices(Result)
    variables = fieldnames(Result);
    for i = 1:length(variables)
        if ~isstruct(Result.(variables{i}))
            if isnumeric(Result.(variables{i})) & size(Result.(variables{i}),1) == size(Result.(variables{i}),2)
                Result.(variables{i}) = mean(Result.(variables{i}),3);
            end
        else
            variables2 = fieldnames(Result.(variables{i}));
            for j = 1:length(variables2)
                if isnumeric(Result.(variables{i}).(variables2{j})) & ...
                        size(Result.(variables{i}).(variables2{j}),1) == size(Result.(variables{i}).(variables2{j}),2)
                    Result.(variables{i}).(variables2{j}) = mean(Result.(variables{i}).(variables2{j}),3);
                end
            end
        end
    end
end