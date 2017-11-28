% Create eeg-data using AR-model / Freeman neural mass - model / sinuswaves
%
% [Result,cfg] = lab_create_modeldata(cfg)
%
% Written by F. Hatz 2013 Neurology Basel

function [Result,cfg] = lab_create_modeldata(cfg)

Result = [];
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'MODEL')
    [cfg,skipprocessing] = lab_set_create_modeldata(cfg);
    if skipprocessing == 1;
        return
    else
        pause(0.2);
    end
end
chartf = num2str(cfg.MODEL.numtf','%05.0f');
Output_filepath = cfg.MODEL.output_folder;
if isempty(Output_filepath)
    return
end

for j = 1:size(cfg.MODEL.numtf,2)
    numtf = cfg.MODEL.numtf(1,j);
    warning off %#ok<WNOFF>
    mkdir(fullfile(Output_filepath,['TF' chartf(j,:)]));
    warning on %#ok<WNON>
    cfg.Output_filepath = fullfile(Output_filepath,['TF' chartf(j,:)]);
    MatrixAll = [];
    for i = 1:cfg.MODEL.loops
        if isfield(cfg.MODEL,'randmatrix') & cfg.MODEL.randmatrix == true
            matrix = randmio_und(cfg.MODEL.connections,5);
        else
            matrix = cfg.MODEL.connections;
        end
        if isfield(cfg.MODEL,'matrixformat') & cfg.MODEL.matrixformat == 2
            matrix = tril(matrix);
        elseif isfield(cfg.MODEL,'matrixformat') & cfg.MODEL.matrixformat == 3
            matrix = triu(matrix);
        end
        if strcmp(cfg.MODEL.type,'Freeman')
            if exist('lab_eeg_Freeman') %#ok<EXIST>
                disp(['Create random data with Freeman (' num2str(cfg.MODEL.chans) 'x' num2str(numtf) ')'])
                [data,header] = lab_eeg_Freeman(cfg.MODEL.chans,numtf + cfg.MODEL.lag,cfg.MODEL.fs,cfg.MODEL.SET);
                Oname = 'Freeman';
            else
                disp('Abort, script ''lab_eeg_Freeman'' missing')
                Result = [];
                return
            end
        elseif strcmp(cfg.MODEL.type,'AR')
            if exist('lab_eeg_ARmodel') %#ok<EXIST>
                disp(['Create random data with AR (' num2str(cfg.MODEL.chans) 'x' num2str(numtf) ')'])
                [data,header] = lab_eeg_ARmodel(cfg.MODEL.chans,numtf + cfg.MODEL.lag,cfg.MODEL.fs,cfg.MODEL.SET.coeff);
                Oname = 'ARmodel';
            else
                disp('Abort, script ''lab_eeg_ARmodel'' missing')
                Result = [];
                return
            end
        elseif strcmp(cfg.MODEL.type,'PhaseResetting')
            disp(['Create random data with PhaseResetting (' num2str(cfg.MODEL.chans) 'x' num2str(numtf) ')'])
            [data,header] = lab_create_eeg_phaseresetting(cfg.MODEL.chans,numtf + cfg.MODEL.lag,cfg.MODEL.fs,cfg.MODEL.SET);
        elseif strcmp(cfg.MODEL.type,'Sinuswave')
            if ~isfield(cfg.MODEL,'SET') | ~isfield(cfg.MODEL.SET,'freq') | isempty(cfg.MODEL.SET.freq)
                cfg.MODEL.SET.freq = 8;
                cfg.MODEL.SET.randphase = true;
                disp('  invalid frequency, set to 8');
            end
            Oname = 'Sinuswave';
            [data,header] = lab_eeg_sinuswave(cfg.MODEL.chans,numtf + cfg.MODEL.lag,cfg.MODEL.fs,cfg.MODEL.SET.freq, ...
                cfg.MODEL.SET.randphase);
        elseif strcmp(cfg.MODEL.type,'Roessler')
            if exist('rossler') %#ok<EXIST>
                disp(['Create random data with Roessler attractor (' num2str(cfg.MODEL.chans) 'x' num2str(numtf) ')'])
                [data,header] = lab_eeg_Roessler(cfg.MODEL.chans,numtf + cfg.MODEL.lag,cfg.MODEL.fs,cfg.MODEL.SET);
                Oname = 'Roessler';
            else
                disp('Abort, script ''Chaotic Systems Toolbox'' missing')
                Result = [];
                return
            end
        else
            disp('Abort, invalide Model-Type selected')
            Result = [];
            return
        end
        
        if ~isempty(matrix) & ~isempty(cfg.MODEL.lag)
            % calculate connectivity
            [data,header] = lab_create_modeldata_connectivity(data,header,matrix,cfg.MODEL.lag,cfg.MODEL.K);
        end
        if isfield(cfg.MODEL,'NOISE') & ~isempty(cfg.MODEL.NOISE)
            [data,header,cfg.MODEL] = lab_add_noise(data,header,cfg.MODEL);
        end
        cfg.Output_file = [ Oname '_' num2str(cfg.MODEL.chans) 'x' chartf(j,:) '_' num2str(i) '.sef'];
        lab_save_data(data,header,fullfile(cfg.Output_filepath,cfg.Output_file));
        if isfield(cfg.MODEL,'CONNECT') & ~isempty(cfg.MODEL.CONNECT)
            cfgtmp = cfg;
            cfgtmp.CONNECT = cfg.MODEL.CONNECT;
            if isfield(cfg.MODEL,'GRAPH')
                cfgtmp.GRAPH = cfg.MODEL.GRAPH;
            end
            try %#ok<TRYNC>
                [result,cfgtmp] = lab_calculate_connectivity(data,header,cfgtmp);
                cfg.MODEL.CONNECT = cfgtmp.CONNECT;
                if isfield(cfgtmp,'GRAPH')
                    cfg.MODEL.GRAPH = cfgtmp.GRAPH;
                end
                cfg.listold = cfgtmp.listold;
            end
        end
        if ~isempty(matrix)
            lab_write_matrix(fullfile(cfg.Output_filepath,[cfg.Output_file(1:end-4) '_matrix.txt']),matrix);
            MatrixAll(:,:,i) = matrix; %#ok<AGROW>
            matrix = matrix + matrix';
            lab_write_matrix(fullfile(cfg.Output_filepath,[cfg.Output_file(1:end-4) '_undirect_matrix.txt']),matrix);
        else
            MatrixAll(:,:,i) = zeros(size(data,1),size(data,1)); %#ok<AGROW>
        end
    end
    if isfield(cfg.MODEL,'CONNECT') & ~isempty(cfg.MODEL.CONNECT)
        C_pli = [];
        C_dpli = [];
        C_plv = [];
        C_plt = [];
        C_SL = [];
        C_CAE = [];
        C_AEC = [];
        U_pli = [];
        U_dpli = [];
        U_plv = [];
        U_plt = [];
        U_SL = [];
        U_CAE = [];
        U_AEC = [];
        for i = 1:cfg.MODEL.loops
            matrix = MatrixAll(:,:,i);
            matrix(1:size(matrix,1)+1:end) = 0;
            if max(tril(matrix)) == 0
                matrix = triu(matrix) + triu(matrix)';
            elseif max(triu(matrix)) == 0
                matrix = tril(matrix) + tril(matrix)';
            end
            matrix(1:size(matrix,1)+1:end) = NaN;
            Connected = (matrix>0);
            Unconnected = (matrix==0);
            if isfield(result,'pli') & ~isempty(result.pli) & size(result.pli,3) >= i
                C_pli = cat(1,C_pli,sum(result.pli(:,:,i) .* Connected,1) ./ sum(Connected,1));
                U_pli = cat(1,U_pli,sum(result.pli(:,:,i) .* Unconnected,1) ./ sum(Unconnected,1));
            end
            if isfield(result,'dpli') & ~isempty(result.dpli) & size(result.dpli,3) >= i
                C_dpli = cat(1,C_dpli,sum(result.dpli(:,:,i) .* Connected,1) ./ sum(Connected,1));
                U_dpli = cat(1,U_dpli,sum(result.dpli(:,:,i) .* Unconnected,1) ./ sum(Unconnected,1));
            end
            if isfield(result,'plv') & ~isempty(result.plv) & size(result.plv,3) >= i
                C_plv = cat(1,C_plv,sum(result.plv(:,:,i) .* Connected,1) ./ sum(Connected,1));
                U_plv = cat(1,U_plv,sum(result.plv(:,:,i) .* Unconnected,1) ./ sum(Unconnected,1));
            end
            if isfield(result,'plt') & ~isempty(result.plt) & size(result.plt,3) >= i
                C_plt = cat(1,C_plt,sum(result.plt(:,:,i) .* Connected,1) ./ sum(Connected,1));
                U_plt = cat(1,U_plt,sum(result.plt(:,:,i) .* Unconnected,1) ./ sum(Unconnected,1));
            end
            if isfield(result,'SL') & ~isempty(result.SL) & size(result.SL,3) >= i
                C_SL = cat(1,C_SL,sum(result.SL(:,:,i) .* Connected,1) ./ sum(Connected,1));
                U_SL = cat(1,U_SL,sum(result.SL(:,:,i) .* Unconnected,1) ./ sum(Unconnected,1));
            end
            if isfield(result,'AEC') & ~isempty(result.AEC) & size(result.AEC,3) >= i
                C_AEC = cat(1,C_AEC,sum(result.AEC(:,:,i) .* Connected,1) ./ sum(Connected,1));
                U_AEC = cat(1,U_AEC,sum(result.AEC(:,:,i) .* Unconnected,1) ./ sum(Unconnected,1));
            end
            if isfield(result,'CAE') & ~isempty(result.CAE) & size(result.CAE,3) >= i
                C_CAE = cat(1,C_CAE,sum(result.CAE(:,:,i) .* Connected,1) ./ sum(Connected,1));
                U_CAE = cat(1,U_CAE,sum(result.CAE(:,:,i) .* Unconnected,1) ./ sum(Unconnected,1));
            end
        end
        if ~isempty(C_pli)
            eval(['Result.Connected.TF' num2str(numtf) '.pli = C_pli;']);
        end
        if ~isempty(C_dpli)
            eval(['Result.Connected.TF' num2str(numtf) '.dpli = C_dpli;']);
        end
        if ~isempty(C_plv)
            eval(['Result.Connected.TF' num2str(numtf) '.plv = C_plv;']);
        end
        if ~isempty(C_plt)
            eval(['Result.Connected.TF' num2str(numtf) '.plt = C_plt;']);
        end
        if ~isempty(C_SL)
            eval(['Result.Connected.TF' num2str(numtf) '.SL = C_SL;']);
        end
        if ~isempty(C_AEC)
            eval(['Result.Connected.TF' num2str(numtf) '.AEC = C_AEC;']);
        end
        if ~isempty(C_CAE)
            eval(['Result.Connected.TF' num2str(numtf) '.CAE = C_CAE;']);
        end
        if ~isempty(U_pli)
            eval(['Result.Unconnected.TF' num2str(numtf) '.pli = U_pli;']);
        end
        if ~isempty(U_dpli)
            eval(['Result.Unconnected.TF' num2str(numtf) '.dpli = U_dpli;']);
        end
        if ~isempty(U_plv)
            eval(['Result.Unconnected.TF' num2str(numtf) '.plv = U_plv;']);
        end
        if ~isempty(U_plt)
            eval(['Result.Unconnected.TF' num2str(numtf) '.plt = U_plt;']);
        end
        if ~isempty(U_SL)
            eval(['Result.Unconnected.TF' num2str(numtf) '.SL = U_SL;']);
        end
        if ~isempty(U_AEC)
            eval(['Result.Unconnected.TF' num2str(numtf) '.AEC = U_AEC;']);
        end
        if ~isempty(U_CAE)
            eval(['Result.Unconnected.TF' num2str(numtf) '.CAE = U_CAE;']);
        end
    end
end
if isfield(cfg.MODEL,'CONNECT') & ~isempty(cfg.MODEL.CONNECT)
    save(fullfile(Output_filepath,'MODELResult.mat'),'Result','cfg');
end