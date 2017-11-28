function MinNorm = lab_test_smoothing(data,LOCS,distances,settings,Filename,silent)

if ~exist('silent','var')
    silent = false;
end
if ~exist('settings','var') | isempty(settings)
    settings = lab_set_smoothing_auto;
    if isempty(settings)
        MinNorm = [];
        return
    else
        pause(0.2);
    end
end

if settings.dorelativ == true
    for i = 1:size(data,3);
        data(:,:,i) = data(:,:,i) ./ repmat(sum(data(:,:,i),2),1,size(data,2));
    end
end

% Recalculate LOCS.theta
LOCS.theta = LOCS.theta .* pi/180;

% Main calculation
disp('      Find optimal maximal distance')
Nchans = size(data,1);
DiffNorm = [];
DiffCluster = [];
header.locs = LOCS;
header.samplingrate = [];
for P = 1:settings.loops
    fprintf(['        Loop ' num2str(P) ':'])
    excludeR = 1:Nchans;
    excludeR(2,:) = random('unif',1,Nchans,1,Nchans);
    excludeR = sortrows(excludeR',2)';
    excludeR = excludeR(1,1:settings.omitted);
    excludeR = sortrows(excludeR')';
    for O = 0:settings.steps:settings.maxdistance
        % calculate guassian distribution per channel
        sigma = (O^2 / (-2 * log(settings.weightmaxdistance / 100)))^0.5;
        neighbours = zeros(Nchans,Nchans);
        for i = 1:Nchans;
            for j = 1:Nchans;
                if distances(i,j) <= O
                    if O == 0
                        neighbours(i,j) = 1;
                    else
                        neighbours(i,j) = exp(-distances(i,j)^2 / (sigma^2 * 2));
                    end
                else
                    neighbours(i,j) = 0;
                end
            end
        end
        neighbours(excludeR,excludeR) = 0;
        for i = 1:size(neighbours,1)
            neighbours(i,:) = neighbours(i,:) / sum(neighbours(i,:),2);
        end
        neighbours(excludeR,:) = 0;
        
        % calculate smoothed data
        dataT = data;
        for i = 1:size(dataT,3)
            dataT(:,:,i) = (neighbours *  dataT(:,:,i));
        end
        %clearvars j i sigma neighbours

        %Interpolate bad
        header.badchans = excludeR;
        [dataT,header] = lab_interpolate_bad(dataT,header,settings.method,1);
        
        % calculate differences
        Diff = data(excludeR,:,:) - dataT(excludeR,:,:);
        Diff = sqrt(sum(Diff(:).^2));
        DiffNorm = [DiffNorm Diff]; %#ok<AGROW>
        DiffCluster = [DiffCluster O]; %#ok<AGROW>
        clearvars Diff
        fprintf('.')
    end
    disp(':');
    % scatter(DiffCluster,DiffNorm,'.r','DisplayName','Diff vs DiffCluster','XDataSource','DiffCluster','YDataSource','Diff');figure(gcf)
end

% Plot Result
if exist('Filename','var') & ~isempty(Filename)
    x = 0:0.25:settings.maxdistance;
    pNorm = polyfit(DiffCluster,DiffNorm,5);
    Fnorm = polyval(pNorm,x);
    if silent == true
        fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Find optimal maximal distance','Menubar','none','Visible','off');
    else
        fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Find optimal maximal distance','Menubar','none');
        m1 = uimenu(fig1,'Label','File');
        uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
    end
    plot(DiffCluster,DiffNorm,'.',x,Fnorm,'-');
    title('Summed Error');
    lab_print_figure([Filename '_FindDistance.pdf'],fig1);
    if silent == true
        close(fig1);
    end
    
    % Write .xls
    warning off %#ok<WNOFF>
    if exist([Filename '_FindDistance.xlsx'],'file')
        delete([Filename '_FindDistance.xlsx']);
    end
    xlsout = cat(1,{'Distance','Summed Error'},[num2cell(DiffCluster') num2cell(DiffNorm')]);
    lab_write_xls([Filename '_FindDistance.xlsx'],xlsout,'SummedError');
    MinNorm = x(Fnorm==min(Fnorm));
    if MinNorm <= 0
        MinNorm = NaN;
    end
    xlsout = {'Optimal Distance';'Fitting Curve:';' Factor1';' Factor2';' Factor3';' Factor4';' Factor5';' Factor6'};
    xlsout = [xlsout cat(1,num2cell(MinNorm),{[]},num2cell(pNorm'))];
    lab_write_xls([Filename '_FindDistance.xlsx'],xlsout,'OptimalDistance');
    warning on %#ok<WNON>
end

end


