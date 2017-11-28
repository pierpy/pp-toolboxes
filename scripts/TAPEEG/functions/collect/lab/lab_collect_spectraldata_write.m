% Helper function for lab_collect_spectras
%
% cfg = lab_collect_spectraldata_write(cfg,calc,result,Mappings)
%
% written by F. Hatz 2012

function cfg = lab_collect_spectraldata_write(cfg,calc,result,Mappings)

cd(calc.result_path)

% Prepare BA & Q
if isfield(result,'BA')
    result.BA = result.BA(1,:);
end
if isfield(result,'Q')
    result.Q = result.Q(1,:);
end
if isfield(result,'EPF')
    result.EPF = result.EPF(1,:);
    result.EMF = result.EMF(1,:);
end
if isfield(result,'PF')
    result.PF = result.PF(1,:);
    result.MF = result.MF(1,:);
end
if isfield(result.R,'bchans')
    result.Q = [result.Q cellstr('Bad channels')];
end
if isfield(result.R,'badact')
    result.Q = [result.Q cellstr('Bad activations')];
end
if isfield(result.R,'good')
    result.Q = [result.Q cellstr('Good')];
end
for N = 1:length(result.R)
    R = result.R(1,N);
    patient = result.patient(1,N);
    
    % Collect BA
    BA = [patient num2cell(R.cogfreq) num2cell(R.peakfreq)];
    result.BA = cat(1,result.BA,BA);
    
    % Collect quality
    peakfreq = R.peakfreq;
    peakfreq(isnan(peakfreq)) = 0;
    Qheader = result.Q;
    Q = [patient num2cell(R.cogfreq) num2cell(R.cogareapower) ...
        num2cell(peakfreq) num2cell(R.areapower) num2cell(R.peakamp) ...
        num2cell(R.peak2min) num2cell(R.peakratio) num2cell(R.bplast)];
    if isfield(R,'bchans')
        if ~isempty(R.bchans)
            Q = [Q num2cell(double(R.bchans))];
        else
            Q = [Q num2cell(0)];
        end
    end
    if isfield(R,'badact')
        if ~isempty(R.badact)
            Q = [Q num2cell(double(R.badact(1,1)))];
        else
            Q = [Q num2cell(0)];
        end
    end
    if isfield(R,'good')
        Q = [Q num2cell(double(R.good))];
    end
    result.Q = cat(1,result.Q,Q);
    
    % Collect every channel
    resulttmp = [patient num2cell(R.peakfreqAll')];
    result.EPF = cat(1,result.EPF,resulttmp);
    resulttmp = [patient num2cell(R.cogfreqAll')];
    result.EMF = cat(1,result.EMF,resulttmp);
    
    % Collect mappings
    if isfield(R,'peakfreqMap')
        resulttmp = [patient num2cell(R.peakfreqMap')];
        result.PF = cat(1,result.PF,resulttmp);
        resulttmp = [patient num2cell(R.cogfreqMap')];
        result.MF = cat(1,result.MF,resulttmp);
    end
end

if size(result.BA,1) < 255
    xlstag = '.xlsx';
else
    xlstag = '.xls';
end
lab_write_xls([calc.output_fileBA xlstag],result.BA');
lab_write_xls([calc.output_fileEPF xlstag],result.EPF');
lab_write_xls([calc.output_fileEMF xlstag],result.EMF');
if isfield(result,'PF') & ~isempty(result.PF) & ~ischar(result.PF) & ~strcmp(result.PF,'empty')
    lab_write_xls([calc.output_filePF xlstag],result.PF');
    lab_write_xls([calc.output_fileMF xlstag],result.MF');
end
if isfield(result,'Q') & ~ischar(result.Q) & ~strcmp(result.Q,'empty')
    lab_write_xls([calc.output_fileQ xlstag],result.Q');
end
if isfield(result,'BP') & ~isempty(result.BP) & ~ischar(result.BP) & ~strcmp(result.BP,'empty')
    if size(result.BP,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileBP xlstag], result.BP');
    lab_write_xls([calc.output_fileBPM xlstag], result.BPM');
    lab_write_xls([calc.output_fileBPR xlstag], result.BPR');
    lab_write_xls([calc.output_fileBPMR xlstag], result.BPMR');
    lab_write_xls([calc.output_fileBP 'Log' xlstag], result.BPlog');
    lab_write_xls([calc.output_fileBPM 'Log' xlstag], result.BPMlog');
    lab_write_xls([calc.output_fileBPR 'Logit' xlstag], result.BPRlog');
    lab_write_xls([calc.output_fileBPMR 'Logit' xlstag], result.BPMRlog');
end
if isfield(result,'EBP') & ~isempty(result.EBP) & ~ischar(result.EBP) & ~strcmp(result.EBP,'empty')
    if size(result.EBP,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileEBP xlstag], result.EBP');
    lab_write_xls([calc.output_fileEBPM xlstag], result.EBPM');
    lab_write_xls([calc.output_fileEBPR xlstag], result.EBPR');
    lab_write_xls([calc.output_fileEBPMR xlstag], result.EBPMR');
    lab_write_xls([calc.output_fileEBP 'Log' xlstag], result.EBPlog');
    lab_write_xls([calc.output_fileEBPM 'Log' xlstag], result.EBPMlog');
    lab_write_xls([calc.output_fileEBPR 'Logit' xlstag], result.EBPRlog');
    lab_write_xls([calc.output_fileEBPMR 'Logit' xlstag], result.EBPMRlog');
end
if isfield(result,'EBE') & ~isempty(result.EBE) & ~ischar(result.EBE) & ~strcmp(result.EBE,'empty')
    if size(result.EBE,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileEBE xlstag], result.EBE');
end
if isfield(result,'ETE') & ~isempty(result.ETE) & ~ischar(result.ETE) & ~strcmp(result.ETE,'empty')
    if size(result.ETE,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileETE xlstag], result.ETE');
end
if isfield(result,'EBEM') & ~isempty(result.EBEM) & ~ischar(result.EBEM) & ~strcmp(result.EBEM,'empty')
    if size(result.EBEM,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileEBEM xlstag], result.EBEM');
end
if isfield(result,'ETEM') & ~isempty(result.ETEM) & ~ischar(result.ETEM) & ~strcmp(result.ETEM,'empty')
    if size(result.ETEM,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileETEM xlstag], result.ETEM');
end
if isfield(result,'BE') & ~isempty(result.BE) & ~ischar(result.BE) & ~strcmp(result.BE,'empty')
    if size(result.BE,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileBE xlstag], result.BE');
end
if isfield(result,'TE') & ~isempty(result.TE) & ~ischar(result.TE) & ~strcmp(result.TE,'empty')
    if size(result.TE,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileTE xlstag], result.TE');
end
if isfield(result,'BEM') & ~isempty(result.BEM) & ~ischar(result.BEM) & ~strcmp(result.BEM,'empty')
    if size(result.BEM,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileBEM xlstag], result.BEM');
end
if isfield(result,'TEM') & ~isempty(result.TEM) & ~ischar(result.TEM) & ~strcmp(result.TEM,'empty')
    if size(result.TEM,1) > 255
        xlstag = '.xlsx';
    else
        xlstag = '.xls';
    end
    lab_write_xls([calc.output_fileTEM xlstag], result.TEM');
end
% Write verbose file (*.vrb)
fid=fopen('FrequencyAnalysis.vrb','w');
fprintf(fid,'Frequency Analysis\n');
fprintf(fid,datestr(now,0));
fprintf(fid,'\n\n');
fprintf(fid,'Channels for Backgroundactivity\n');
fprintf(fid,num2str(calc.mappingBA));
fprintf(fid,'\n\n');
fprintf(fid,'Background activity low freq (peak)\n');
fprintf(fid,num2str(cfg.CollectFFT.lowfreqpeak));
fprintf(fid,'\n\n');
fprintf(fid,'Background activity high freq (peak)\n');
fprintf(fid,num2str(cfg.CollectFFT.highfreqpeak));
fprintf(fid,'\n\n');
fprintf(fid,'Background activity low freq (median)\n');
fprintf(fid,num2str(cfg.CollectFFT.lowfreqcog));
fprintf(fid,'\n\n');
fprintf(fid,'Background activity high freq (median)\n');
fprintf(fid,num2str(cfg.CollectFFT.highfreqcog));
fprintf(fid,'\n\n');
fprintf(fid,'Spectral bands\n');
tmp = char(calc.bandtitles);
tmp = [tmp char(32*ones(size(tmp,1), 3))]';
fprintf(fid,tmp);
fclose(fid);

% Plot Spectras
if size(result.FBM,3) > 1
    SpectChannels = mean(result.FBM,3);
else
    SpectChannels = result.FBM;
end
if isfield(cfg.CollectFFT,'plotchans') & cfg.CollectFFT.plotchans == 1
    % Plot for every channel
    for i = 1:size(SpectChannels,1);
        f = figure('Visible','off');
        area(SpectChannels(i,:));
        set(gca,'XTick',1:length(calc.freqlabel));
        set(gca,'XTickLabel',calc.freqlabel,'FontName','Times','fontsize',9);
        lab_print_figure(fullfile(calc.resultplots_path,['Channel-' num2str(i) '.jpg']),f);
        close;
        xlsout = num2cell(SpectChannels);
        if size(xlsout,2) > 255
            fileout = fullfile(calc.resultplots_path,'Channels.xlsx');
        else
            fileout = fullfile(calc.resultplots_path,'Channels.xls');
        end
        lab_write_xls(fileout,xlsout);
        clearvars xlsout
    end
end
% Plot for mappings
if exist('Mappings','var') & ~isempty(Mappings)
    if isfield(cfg.CollectFFT,'plotmappings') & cfg.CollectFFT.plotmappings == 1
        for i = 1:size(Mappings.mappings,2)
            if size(Mappings.mappings{1,i},2) > 1;
                SpectMappings(:,:,i)= permute(mean(result.FBM(Mappings.mappings{1,i},:,:),1),[3 2 1]); %#ok<AGROW>
            else
                SpectMappings(:,:,i)= permute(result.FBM(Mappings.mappings{1,i},:,:),[3 2 1]); %#ok<AGROW>
            end
        end
        for i = 1:size(SpectMappings,1)
            for j = 1:size(SpectMappings,3)
                f = figure('Visible','off');
                area(SpectMappings(i,:,j));
                set(gca,'XTick',1:length(calc.freqlabel));
                set(gca,'XTickLabel',calc.freqlabel,'FontName','Times','fontsize',9);
                if size(Mappings.mappingstitle,1) == size(Mappings.mappings,2)
                    title(Mappings.mappingstitle{j,1});
                end
                lab_print_figure(fullfile(calc.resultplots_path,['Mapping-' num2str(j) '_' result.FBMpat{1,i} '.jpg']),f);
                close(f);
            end
        end
    else
        for i = 1:size(Mappings.mappings,2)
            if size(Mappings.mappings{1,i},2) > 1;
                SpectMappings(i,:)= mean(SpectChannels(Mappings.mappings{1,i},:)); %#ok<AGROW>
            else
                SpectMappings(i,:)= SpectChannels(Mappings.mappings{1,i},:); %#ok<AGROW>
            end
        end
        for i = 1:size(SpectMappings,1);
            f = figure('Visible','off');
            area(SpectMappings(i,:));
            set(gca,'XTick',1:length(calc.freqlabel));
            set(gca,'XTickLabel',calc.freqlabel,'FontName','Times','fontsize',9);
            if size(Mappings.mappingstitle,1) == size(Mappings.mappings,2)
                title(Mappings.mappingstitle{i,1});
            end
            lab_print_figure(fullfile(calc.resultplots_path,['Mapping-' num2str(i) '.jpg']),f);
            close;
        end
        xlsout = num2cell(SpectMappings);
        if size(xlsout,2) > 255
            fileout = fullfile(calc.resultplots_path,'Mappings.xlsx');
        else
            fileout = fullfile(calc.resultplots_path,'Mappings.xls');
        end
        lab_write_xls(fileout,xlsout);
        clearvars xlsout
    end
end

% Plot Spectra for background activity
if isfield(calc,'mappingBA') & ~isempty(calc.mappingBA)
    SpectBA = mean(SpectChannels(calc.mappingBA,:));
    f = figure('Visible','off');
    area(SpectBA);
    set(gca,'XTick',1:length(calc.freqlabel));
    set(gca,'XTickLabel',calc.freqlabel,'FontName','Times','fontsize',9);
    title('Background Activity Mean');
    lab_print_figure(fullfile(calc.resultplots_path,'BackgroundActivityMean.jpg'),f);
    close;
    SpectBA = mean(permute(result.FBM(calc.mappingBA,:,:),[3 2 1]),3);
    xlsout = [result.patient' num2cell(SpectBA)];
    if isfield(calc,'freqlabelall')
        xlsout = cat(1,[cellstr('Patient') calc.freqlabelall],xlsout);
    else
        xlsout = cat(1,[cellstr('Patient') calc.freqlabel],xlsout);
    end
    if size(xlsout,2) > 255
        fileout = fullfile(calc.resultplots_path,'BackgroundActivityMean.xlsx');
    else
        fileout = fullfile(calc.resultplots_path,'BackgroundActivityMean.xls');
    end
    lab_write_xls(fileout,xlsout);
    clearvars xlsout SpectBA
end

return