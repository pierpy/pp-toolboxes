% calculate SL with correction for volume conduction
%
% [RESULT1, RESULT2] = lab_SLcorrect(S_in,lag,m,w1,w2,pref,speed,option,dofull)
%    S_in = timeframes x channels
%
% [RESULT1, RESULT2] = lab_SLcorrect(data,header)
%    data = channels x timeframes
%
% other input values: see Stam et al. 2002
% dofull: true = bidrectional correction / false = undirectional correction

function [RESULT1,RESULT2] = lab_SLcorrect(S_in,lag,m,w1,w2,pref,speed,option,StoreHits,Mcorrect)

if nargin == 2 & isfield(lag,'samplingrate')
    header = lag;
    clearvars lag
    S_in = S_in';
    settings = lab_get_SL([],1);
    pause(0.2);
    pref = settings.SLc.pref;
    speed = floor(header.samplingrate / max(settings.SLc.freqs));
    lag = round(header.samplingrate / (3*max(settings.SLc.freqs)));
    m = round(3 * (max(settings.SLc.freqs) / min(settings.SLc.freqs)))+1;
    w1 = 2 * lag * (m - 1);
    w2 = size(S_in,1);
    clearvars header
    option = ' ';
    Mcorrect = settings.SLc.Mcorrect;
    StoreHits = settings.SLc.storehits;
end

if ~exist('Mcorrect','var')
    Mcorrect = 1;
end
if ~exist('StoreHits','var')
    StoreHits = false;
end
if ~exist('option','var')
    option = ' ';
end

if Mcorrect == 3
    disp('   SL corrected for volume conduction (bidirectional)')
elseif Mcorrect == 2
    disp('   SL corrected for volume conduction (unidirectional asymmetrical)')
else
    disp('   SL corrected for volume conduction (unidirectional symmetrical)')
end

% start counter
tic;

% look for bad channels
for i = 1:size(S_in,2)
    if max(S_in(:,i)) == 0 & min(S_in(:,i)) == 0
        BadFlag(i) = false;
    elseif max(isnan(S_in(:,i)))
        BadFlag(i) = false;
    else
        BadFlag(i) = true;
    end
end

% prepare calculation and control values
[num_samples, num_chans] = size(S_in);
S_in = cat(1,S_in,zeros(lag*m,num_chans));
m1 = 1;           % first sample
m2 = num_samples; % last sample, tested using m2 = 4096
num_it = floor( (m2 - lag*(m-1))/speed ) - ceil( m1/speed ) + 1; % number of iterations
if (num_chans <= 1)
    fprintf(2,'   ERROR: must have at least two channels in use\n');
    RESULT1 = [];
    RESULT2 = [];
    return;
end %if
if (m1 <= 0)
    fprintf(2,'   ERROR: m1 must be positive\n');
    RESULT1 = [];
    RESULT2 = [];
    return;
end %if
if (m1 >= m2 - lag*(m-1))
    fprintf(2,'   ERROR: m1 must be less than m2 - l*(m-1)\n');
    RESULT1 = [];
    RESULT2 = [];
    return;
end %if
if (w2 - w1 < 2 | w2 < 2)
    fprintf(2,'   ERROR: must be at least one positive integer b/t w1 and w2\n');
    RESULT1 = [];
    RESULT2 = [];
    return;
end %if
maxw1 = floor((m2-lag*(m-1)-m1-1)/2);
if (w1 > maxw1 )
    if maxw1 >= 0
        fprintf(2,'   WARNING: w1 to large, set to %d\n',floor((m2-lag*(m-1)-m1-1)/2) );
    else
        fprintf(2,'   ERROR: lag and m are too large. Please reduce these values.\n');
        RESULT1 = [];
        RESULT2 = [];
        return;
    end %if
    w1 = floor((m2-lag*(m-1)-m1-1)/2);
    w2 = 10/pref + w1 - 1;
end %if
if (num_it == 0)
    fprintf(2,'   ERROR: speed must be reduced\n');
    RESULT1 = [];
    RESULT2 = [];
    return;
end %if
epsilon = ones(num_chans,1); % epsilon(of channel k) at time i

% set start and end points for calculation (for all channel combinations the same)
Nspeed = 1:speed:(m2 - lag*(m-1));
j = m1 : (m2 - lag*(m-1));
for i = 1:length(Nspeed)
    valid_range = abs(Nspeed(i)-j)>w1 & abs(Nspeed(i)-j)<w2; % vector of valid range positions
    Vstart = find(diff(valid_range) == 1)+1;
    Vend = find(diff(valid_range) == -1);
    if isempty(Vstart)
        Vstart = 1;
    elseif min(Vend) < min(Vstart)
        Vstart = [1 Vstart];
    end
    if isempty(Vend)
        Vend = length(valid_range);
    elseif length(Vstart) > length(Vend)
        Vend = [Vend length(valid_range)];
    end
    Vstart = Vstart + m1 - 1;
    Vend = Vend + m1 - 1;
    Range{i,1} = Vstart;
    Range{i,2} = Vend;
end

RESULT1 = zeros(num_chans,num_chans,length(Nspeed)); %S(k,l,i) matrix
if StoreHits == true
    RESULT2 = zeros(num_chans,num_chans,length(Nspeed));
else
    RESULT2 = [];
end

% do calculation
if Mcorrect == 3
    calc_bidirectional;
elseif Mcorrect == 2
    calc_unidirectional;
else
    calc_unidirectional;
    for i = 1:size(RESULT1,3)
        RESULT1(:,:,i) = (RESULT1(:,:,i) + tril(RESULT1(:,:,i))' + triu(RESULT1(:,:,i))') ./ 2;
    end   
end

if option ~= 's'
    % overall synchronization likelihood
    S = sum(sum(sum(RESULT1,1)-1))/(num_chans-1)/num_chans/length(Nspeed);
    fprintf('   Results: Overall synchronization likelihood: %f\n', S);
    fprintf('            Elapsed time: %.3f seconds\n',toc);
end

if exist('settings','var') & isfield(settings,'SLc') & ...
        isfield(settings.SLc,'storesingle') & settings.SLc.storesingle == true
    RESULT1 = mean(RESULT1,3);
end

    function calc_bidirectional
        fprintf('     SLc-Progress: ');
        S_matrix = zeros(2,2,length(Nspeed)); %S(k,l,i) matrix
        if ~strcmp(mexext,'mexw32')
            hit_matrix = zeros(2,2,length(Nspeed));
        else
            hit_matrix = [];
        end
        startChan = 0;
        numperdot = floor(num_chans/20);
        for Ci = 1:num_chans
            if mod(Ci,numperdot) == 0
                fprintf('.'); % character to display
            end
            startChan = startChan + 1;
            for Cj = startChan:num_chans
                if Ci ~= Cj & BadFlag(Ci) ~= false & BadFlag(Cj) ~= false
                    [~,~,statstmp] = glmfit(S_in(:,Cj),S_in(:,Ci),'normal');
                    S_tmp(:,1) = statstmp.resid;
                    [~,~,statstmp] = glmfit(S_in(:,Ci),S_in(:,Cj),'normal');
                    S_tmp(:,2) = -statstmp.resid;
                    for Bi = 1:length(Nspeed)
                        seed = S_tmp(Nspeed(Bi)+lag*(0:(m-1)),:);
                        seed = repmat(seed,[1 1 lag]);
                        seed = permute(seed,[3 1 4 2]);
                        euclid4_table = [];
                        for k = 1:length(Range{Bi,1})
                            Vend2 = Range{Bi,2}(k) + (m-1)*lag;
                            for Bj = 1:m
                                startS = Range{Bi,1}(k)+(Bj-1)*lag;
                                shifts = floor((Range{Bi,2}(k) + lag*m - startS) / (lag*m));
                                stopS = startS + (shifts * lag * m) - 1;
                                if stopS > Vend2
                                    skip = stopS - Vend2;
                                else
                                    skip = 0;
                                end
                                if startS < stopS
                                    tmp = reshape(S_tmp(startS:stopS,:),lag,m,shifts,2);
                                    tmp = reshape(sum((repmat(seed,[1 1 shifts 1])-tmp).^2,2).^2,(lag*shifts),2);
                                    euclid4_table = cat(1,euclid4_table,tmp(1:end-skip,:));
                                end
                            end
                        end
                        euclid4_table = euclid4_table(~isnan(euclid4_table(:,1)),:);
                        num_validj = size(euclid4_table,1);
                        euclid4_table = euclid4_table';
                        for k = 1:2
                            sorted_table = sort( euclid4_table(k,:) ); % size (1,validj)
                            epsilon(k) = sorted_table( ceil( pref * num_validj ) );
                        end
                        hit_table = ( euclid4_table <= ( epsilon(1:2) * ones(1,num_validj) ) );
                        hit_table = double(hit_table);
                        hit_table2 = hit_table * hit_table';
                        num_hitsperchan = diag(hit_table2);
                        if ~strcmp(mexext,'mexw32')
                            hit_matrix(:,:,Bi) = hit_table2;
                        end
                        S_tmp2 = hit_table2 ./ ( num_hitsperchan * ones(1,2) );
                        S_tmp2(1:size(S_tmp2,1)+1:end) = 0;
                        S_matrix(:,:,Bi) = S_tmp2;
                    end
                    RESULT1(Ci,Cj,:) = S_matrix(1,2,:);
                    RESULT1(Cj,Ci,:) = S_matrix(2,1,:);
                    if ~strcmp(mexext,'mexw32') & StoreHits == true
                        RESULT2(Ci,Cj,:) = hit_matrix(1,2,:);
                        RESULT2(Cj,Ci,:) = hit_matrix(2,1,:);
                    end
                end
            end
        end
        fprintf('done\n');
    end

    function calc_unidirectional
        fprintf('     SLc-Progress: ');
        S_matrix = zeros(num_chans,num_chans,length(Nspeed)); %S(k,l,i) matrix
        if ~strcmp(mexext,'mexw32') & StoreHits == true
            hit_matrix = zeros(num_chans,num_chans,length(Nspeed));
        else
            hit_matrix = [];
        end
        numberdot = floor(num_chans/20);
        for Ci = 1:num_chans
            if mod(Ci,numberdot) == 0
                fprintf('.'); % character to display
            end
            S_tmp = S_in;
            if BadFlag(Ci) ~= false
                for Cj = 1:num_chans
                    if Ci ~= Cj & BadFlag(Cj) ~= false
                        [~,~,statstmp] = glmfit(S_in(:,Ci),S_in(:,Cj),'normal');
                        S_tmp(:,Cj) = statstmp.resid;
                    end
                end
                for Bi = 1:length(Nspeed)
                    seed = S_tmp(Nspeed(Bi)+lag*(0:(m-1)),:);
                    seed = repmat(seed,[1 1 lag]);
                    seed = permute(seed,[3 1 4 2]);
                    euclid4_table = [];
                    for k = 1:length(Range{Bi,1})
                        Vend2 = Range{Bi,2}(k) + (m-1)*lag;
                        for Bj = 1:m
                            startS = Range{Bi,1}(k)+(Bj-1)*lag;
                            shifts = floor((Range{Bi,2}(k) + lag*m - startS) / (lag*m));
                            stopS = startS + (shifts * lag * m) - 1;
                            if stopS > Vend2
                                skip = stopS - Vend2;
                            else
                                skip = 0;
                            end
                            if startS < stopS
                                tmp = reshape(S_tmp(startS:stopS,:),lag,m,shifts,num_chans);
                                tmp = reshape(sum((repmat(seed,[1 1 shifts 1])-tmp).^2,2).^2,(lag*shifts),num_chans);
                                euclid4_table = cat(1,euclid4_table,tmp(1:end-skip,:));
                            end
                        end
                    end
                    euclid4_table = euclid4_table(~isnan(euclid4_table(:,1)),:);
                    num_validj = size(euclid4_table,1);
                    euclid4_table = euclid4_table';
                    for k = 1 : num_chans
                        sorted_table = sort( euclid4_table(k,:) );
                        epsilon(k) = sorted_table( ceil( pref * num_validj ) );
                    end
                    hit_table = ( euclid4_table <= ( epsilon(1:num_chans) * ones(1,num_validj) ) );
                    hit_table = double(hit_table);
                    hit_table2 = hit_table(Ci,:) * hit_table';
                    if ~strcmp(mexext,'mexw32')
                        hit_matrix(Ci,:,Bi) = hit_table2;
                    end
                    S_tmp2 = hit_table2 ./ hit_table2(Ci);
                    S_tmp2(Ci) = 0;
                    S_matrix(Ci,:,Bi) = S_tmp2;
                end
            end
        end
        tmp = find(BadFlag==false);
        if ~isempty(tmp)
            S_matrix(tmp,:,:) = 0;
            S_matrix(:,tmp,:) = 0;
            if ~strcmp(mexext,'mexw32') & StoreHits == true
                hit_matrix(tmp,:,:) = 0;
                hit_matrix(:,tmp,:) = 0;
            end
        end
        RESULT1 = S_matrix;
        RESULT2 = hit_matrix;
        fprintf('done\n');
    end

end
