% USAGE
% [S_matrix, hit_matrix] = synchronization(S_in,lag,m,w1,w2,pref,speed);
%
% DESCRIPTION
% This function takes in data from the variable S_in and computes S_matrix, described below.
%
% X_in: input 2-D matrix (n,m)
%       M columns, with N samples each column, corresponding to the N samples for each of the M channels.
%
% S_matrix: output 3-D matrix (k,l,i)
%           The value in the (k,l,i) position is the synchronization likelihood of
%            *driver system* k and *response system* l, at calculation time i.
%           NOTE: At each calculation time i, the matrix is not symmetrical!
%                 S(k,l,i) != S(l,k,i)
%
% hit_matrix: output 3-D matrix
%           The value in the (k,l,i) position is the number of hits at both channel k & l, at time i (across all j)
%           The value in the (k,k,i) position is the number of hits at channel k, at time i (across all j)
%
% Other parameters:
%  lag        % lag
%  m          % embedding dimension
%  w1         % window (Theiler correction for autocorrelation)
%  w2         % window (used to sharpen the time resolution of synchronization measure)
%  pref
%  speed

% SYNCHRONIZATION LIKELIHOOD of a driver system X and a response system Y
% 
% PROPERTY (*): If X is in the same state at times i & j, then Y is also in the same state at times i & j.
%
% Generalized synchronization exists when property (*) holds.
% If X & Y exhibit property (*), then the synchronization likelihood between X and Y will be high.
%
% State (def'n): Given l (lag) and m (embedding dimension), the state of X at time i is defined
%  as the vector X_i = ( x_{i}, x_{i+l}, x_{i+2l}, ..., x_{i+(m-1)l} )
%
% Similarity between states (def'n): Two states X_i and X_j are similar if the distance between
%  X_i and X_j is small.  We use the euclidean distance measure (the L2 norm of the vector difference).
%
%
% To have a measure of similarity, we construct an epsilon that depends on the system (X) and time (i).
% We will define that similarity exists if |X_i - X_j| < epsilon_{X,i}.
%
% To determine epsilon_{X,i}, we consider the set of X_j for all j satisfying w1<|i-j|<w2.  From this set,
% we then determine the set of distances between X_i and all X_j.  We then choose epsilon_{X,i} as the
% distance so that the fraction of distances |X_i - X_j| less than epsilon_{X,i} is Pref.
%
% We calculate epsilon_{Y,i} in a similar manner, *using the same Pref* (and same time i).
%
% Now, with i fixed, for each of the values j satisfying w1<|i-j|<w2, we consider a hit wrt X occurs
% (at time j) when |X_i - X_j| < epsilon_{X,i}.  A hit wrt Y occurs when |Y_i - Y_j| < epsilon_{Y,i}.
%
% The synchronization likelihood is defined as the probability (over the valid j's) that
%  there is a hit wrt Y, *given* that there is a hit wrt X.
%
% That is, S(x,y,i) = (# of times there is a hit wrt X and a hit wrt Y) / (# of hits wrt X).
%
%
% If there is no synchronization between X & Y, then the synchronization likelihood is close to Pref.
% If there is synchronization between X & Y, then the synchronization likelihood will approach 1.
%
% Once the whole set of S(x,y,i) is constructed, we can construct variations by averaging over one or
%  more of the three axes (driver system, response system, time).
%
% "Averaging over all [axes] gives S, the overall level of synchronization in a multi channel epoch."
% [2] p.67
%

% SOURCES:
% 1. C.J. Stam, B.W. van Dijk.  "Synchronization likelihood: an unbiased measure of generalized
%     synchronization in multivariate data sets."  Physica D 163 (2002), 236-251.
% 2. C.J. Stam, M. Breakspear, A. van Cappellen van Walsum, B. W. van Dijk.  "Nonlinear
%     Synchronization in EEG and Whole-Head MEG Recordings of Healthy Subjects."  Human Brain
%     Mapping, 19:63-78 (2003).
%

% AUTHOR: Jimmy Chui
%         Modifications 2013: F. Hatz
% LAST REVISION: 18.04.2013
%
% REVISION NOTES
% 18.04.13 - Removed some loops to increase speed (~40x)
% 30.07.03 - Added progress meter
% 18.07.03 - Changed input from (channel, sample) to (sample, channel)
% 05.07.03 - Error message improvement
% 02.07.03 - Algorithm improvements
% 23.06.03 - Code changed to be a secondary file, to be called from a main GUI interface
% 23.06.03 - V1.0
% 23.06.03 - Renamed file to synchronization.m
% 20.06.03 - Optimized code, code cleanup
% 20.06.03 - Code cleanup, added pairwise synchronization output
% 19.06.03 - Reworked code
% 17.06.03 - More comments
% 17.06.03 - Added comments
% 15.06.03 - Added ability to read in data from file
% 15.06.03 - Reverted back to 13.06.03 version, plus minor fix
% 15.06.03 - Changed one algorithm
% 13.06.03 - Minor fix
% 11.06.03 - Optimized parts of code
% 09.06.03 - Implemented draft of code
% 09.06.03 - First version (template)
%--------------------------------------------------

function [RESULT1, RESULT2] = lab_SLslow(S_in,lag,m,w1,w2,pref,speed,option)

% initialize time counter
tic;

% check number of input arguments
narginchk( 7, 8 ) ; % ensure there are 7 (+1 optional) arguments

% set option to null if not given
if nargin == 7
    option = ' ';
end %if
if nargin == 8 && isempty(option)
    option = ' ';
end %if

% read in data, determine size of data
% S_in must be 2-D, size #samples x #chans
[num_samples, num_chans] = size(S_in);

% padding with zeros for calculation
S_in = cat(1,S_in,zeros(lag*m,num_chans));

% initialize variables
m1 = 1;           % first sample
m2 = num_samples; % last sample, tested using m2 = 4096

% output some data, if option ~= 'silent'
if option ~= 's'
    % output data characteristics
    fprintf('Synchronization likelihood (Stam et al. 2002)\n');
    fprintf('   Data: %d channels, %d samples\n',num_chans,num_samples);
    
    % output values of some variables
    fprintf('   Parameters: lag = %d, embedding dimension m = %d\n',lag,m);
    fprintf('               window w1 = %d, w2 = %d\n',w1,w2);
    fprintf('               pref = %.2f\n',pref);
    fprintf('               speed factor: %d\n',speed);
end %if

% calculate number of iterations
num_it = floor( (m2 - lag*(m-1))/speed ) - ceil( m1/speed ) + 1; % number of iterations

% check for errors
% assume m1, m2, w1, w2 positive integers (w1 can be 0)
% need at least 2 channels
if (num_chans <= 1)
    fprintf(2,'   ERROR: must have at least two channels in use\n');
    return;
end %if
% m1 > 0
if (m1 <= 0)
    fprintf(2,'   ERROR: m1 must be positive\n');
    return;
end %if
% m1 < m2 - (m-1)*lag
if (m1 >= m2 - lag*(m-1))
    fprintf(2,'   ERROR: m1 must be less than m2 - l*(m-1)\n');
    return;
end %if
% w1 >= 0, w2 >= w1 + 2
if (w2 - w1 < 2 | w2 < 2)
    fprintf(2,'   ERROR: must be at least one positive integer b/t w1 and w2\n');
    return;
end %if

maxw1 = floor((m2-lag*(m-1)-m1-1)/2);
% ( m2 - (m-1)*lag ) - m1 >= 2 * w1 + 1
if (w1 > maxw1 )
    if maxw1 >= 0
        fprintf(2,'   WARNING: w1 to large, set to %d\n',floor((m2-lag*(m-1)-m1-1)/2) );
    else
        fprintf(2,'   ERROR: lag and m are too large. Please reduce these values.\n');
        return;
    end %if
    w1 = floor((m2-lag*(m-1)-m1-1)/2);
    w2 = 10/pref + w1 - 1;
end %if
if (num_it == 0)
    fprintf(2,'   ERROR: speed must be reduced\n');
    return;
end %if

% initialize variables - inner loop
epsilon = ones(num_chans,1); % epsilon(of channel k) at time i
num_validj = 0; % number of valid j's at time

% initialize variables - outer loop
S_matrix = zeros(num_chans,num_chans,num_it); %S(k,l,i) matrix
hit_matrix = zeros(num_chans,num_chans,num_it);
i_count = 0; % iteration count, used to store matrix entries

% on_display percentage meter
pmdots = 20; % number of dots to display
pm = 0;
if option ~= 's'
    fprintf('   SL-Progress: ');
end

for i = 1:speed:(m2 - lag*(m-1)) %was 'for i = speed:speed:(m2 - lag*(m-1))' (wrong!)
    i_count = i_count + 1;
    for j = (pm+1):floor(i_count/num_it*pmdots)
        if option ~= 's'
            fprintf('.'); % character to display
        end
    end
    pm = floor(i_count/num_it*pmdots);
    
    j = m1 : (m2 - lag*(m-1));
    valid_range = abs(i-j)>w1 & abs(i-j)<w2; % vector of valid range positions
    
    % ---- Old code before 2013 (about 2x slower) ----
    num_validj = sum(valid_range); % number of valid range positions
    valid_range = find(valid_range)+m1-1;
    euclid4_table = zeros(num_chans,num_validj);
    n = 0; % counter
    for j = valid_range
        n = n + 1;
        euclid4_table(:,n) = (sum((S_in(i+lag*(0:(m-1)),:)-S_in(j+lag*(0:(m-1)),:)).^2,1)).^2;
    end
    
    % ---- New Code 2013 ----
%     Vstart = find(valid_range == 1,1);
%     Vend = find(valid_range == 1,1,'last');
%     if Vstart == 1 & Vend == size(valid_range,2)
%         Vend(2) = Vend(1);
%         Vend(1) = find(valid_range == 0,1)-1;
%         Vstart(2) = find(valid_range == 0,1,'last')+1;
%     end
%     Vstart = Vstart + m1 - 1;
%     Vend = Vend + m1 - 1;
%     seed = S_in(i+lag*(0:(m-1)),:);
%     seed = repmat(seed,[1 1 lag]);
%     seed = permute(seed,[3 1 4 2]);
%     euclid4_table = [];
%     for k = 1:length(Vstart)
%         Vend2 = Vend(k) + (m-1)*lag;
%         for j = 1:m
%             startS = Vstart(k)+(j-1)*lag;
%             shifts = floor((Vend(k) + lag*m - startS) / (lag*m));
%             stopS = startS + (shifts * lag * m) - 1;
%             if stopS > Vend2
%                 skip = stopS - Vend2;
%             else
%                 skip = 0;
%             end
%             if startS < stopS
%                 tmp = reshape(S_in(startS:stopS,:),lag,m,shifts,num_chans);
%                 tmp = reshape(sum((repmat(seed,[1 1 shifts 1])-tmp).^2,2).^2,(lag*shifts),num_chans);
%                 euclid4_table = cat(1,euclid4_table,tmp(1:end-skip,:));
%             end
%         end
%     end
%     euclid4_table = euclid4_table(~isnan(euclid4_table(:,1)),:);
%     num_validj = size(euclid4_table,1);
%     euclid4_table = euclid4_table';
    % ---- End New Code 2013 ----
    
    for k = 1 : num_chans
        sorted_table = sort( euclid4_table(k,:) ); % size (1,validj)
        epsilon(k) = sorted_table( ceil( pref * num_validj ) );
    end
    
    % construct 'hit' table, i.e. determine if |X_{k,i} - X_{k,j}| <= epsilon_x for each k & j
    % size (num_chans, num_validj)
    hit_table = ( euclid4_table <= ( epsilon(1:num_chans) * ones(1,num_validj) ) );
    hit_table = double(hit_table); %Matlab 6.5
    
    % construct alternate hit table:
    % at position (k,l), determine the number of hits occuring at both channels k & l (across all j)
    % size (num_chans, num_chans)
    hit_table2 = hit_table * hit_table';
    
    % determine number of hits for each channel, across all j
    % NOTE: this is equivalent to diag(hit_table2)
    %num_hitsperchan = sum( hit_table, 2 ); % size (num_chans,1)
    num_hitsperchan = diag(hit_table2);
    
    % store hit_table2 in a 3-D array
    hit_matrix(:,:,i_count) = hit_table2;
    
    % perform conditional probability calculation
    % divide k^th row by number of hits for channel k
    S_matrix(:,:,i_count) = hit_table2 ./ ( num_hitsperchan * ones(1,num_chans) );
end
if option ~= 's'
    fprintf('done\n');
end

RESULT1 = S_matrix;
RESULT2 = hit_matrix;

% results if option ~= 'silent'
if option ~= 's'
    % overall synchronization likelihood
    S = sum(sum(sum(S_matrix,1)-1))/(num_chans-1)/num_chans/num_it;
    fprintf('   Results: Overall synchronization likelihood: %f\n', S);
    fprintf('            Elapsed time: %.3f seconds\n',toc);
end

return
