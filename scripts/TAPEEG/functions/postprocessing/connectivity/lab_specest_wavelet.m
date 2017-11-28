% FT_SPECEST_WAVELET performs time-frequency analysis on any time series trial
% data using the 'wavelet method' based on Morlet wavelets, doing
% convolution in the time domain by multiplaction in the frequency domain
%
% Use as
%   [spectrum,freqoi,timeoi] = specest_wavelet(dat,time...)
% where
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = array of chan*freqoi*timeoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include:
%   pad        = number, total length of data after zero padding (in seconds)
%   freqoi     = vector, containing frequencies of interest
%   timeoi     = vector, containing time points of interest (in seconds)
%   width      = number or vector, width of the wavelet, determines the temporal and spectral resolution
%   gwidth     = number, determines the length of the used wavelets in standard deviations of the implicit Gaussian kernel
%   verbose    = output progress to console (0 or 1, default 1)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMCONVOL, FT_SPECEST_TFR, FT_SPECEST_HILBERT, FT_SPECEST_MTMFFT

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% $Id: ft_specest_wavelet.m 5106 2012-01-10 13:10:28Z jansch $
%
% Code adapted for TAPEEG 2012-09-24 FHatz Vumc Amsterdam

function [spectrum,freqoi,timeoi] = lab_specest_wavelet(dat, time, varargin)

% get the optional input arguments
freqoi    = lab_getopt(varargin, 'freqoi', 'all');
timeoi    = lab_getopt(varargin, 'timeoi', 'all');
width     = lab_getopt(varargin, 'width', 7);
gwidth    = lab_getopt(varargin, 'gwidth', 3);
pad       = lab_getopt(varargin, 'pad');
fbopt     = lab_getopt(varargin, 'feedback');
verbose   = lab_getopt(varargin, 'verbose', false);

if isempty(fbopt),
  fbopt.i = 1;
  fbopt.n = 1;
end


% Set n's
[nchan,ndatsample] = size(dat);

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad = zeros(1,round((pad - dattime) * fsample));
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data

% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all') % if input was 'all'
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
% check for freqoi = 0 and remove it, there is no wavelet for freqoi = 0
if freqoi(1)==0
  freqoi(1)  = [];
  freqboi(1) = [];
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);

% Set timeboi and timeoi
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
  timeoi   = round(timeoi .* fsample) ./ fsample;
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

% Creating wavelets
% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end
wltspctrm = cell(nfreqoi,1);
for ifreqoi = 1:nfreqoi
  dt = 1/fsample;
  sf = freqoi(ifreqoi) / width(ifreqoi);
  st = 1/(2*pi*sf);
  toi2 = -gwidth*st:dt:gwidth*st;
  A = 1/sqrt(st*sqrt(pi));
  tap = (A*exp(-toi2.^2/(2*st^2)))';
  acttapnumsmp = size(tap,1);
  taplen(ifreqoi) = acttapnumsmp;
  ins = ceil(endnsample./2) - floor(acttapnumsmp./2);
  prezer = zeros(ins,1);
  pstzer = zeros(endnsample - ((ins-1) + acttapnumsmp)-1,1);
  
  % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
  ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./fsample) .* freqoi(ifreqoi));
  
  % create wavelet and fft it
  wavelet = complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer));
  wltspctrm{ifreqoi} = complex(zeros(1,endnsample));
  wltspctrm{ifreqoi} = fft(wavelet,[],1)';
  
  
  %%%% debug plotting
%   figure('name',['wavelet @ ' num2str(freqoi(ifreqoi)) 'Hz' ],'NumberTitle','off');
%   subplot(2,1,1);
%   hold on;
%   plot(real(wavelet));
%   plot(imag(wavelet),'color','r');
%   legend('real','imag');
%   tline = length(wavelet)/2;
%   if mod(tline,2)==0
%     line([tline tline],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--')
%   else
%     line([ceil(tline) ceil(tline)],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--');
%     line([floor(tline) floor(tline)],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--');
%   end;
%   subplot(2,1,2);
%   plot(angle(wavelet),'color','g');
%   if mod(tline,2)==0,
%     line([tline tline],[-pi pi],'color','r','linestyle','--')
%   else
%     line([ceil(tline) ceil(tline)],[-pi pi],'color','r','linestyle','--')
%     line([floor(tline) floor(tline)],[-pi pi],'color','r','linestyle','--')
%   end
  %%%% debug plotting
  
end

% Compute fft
spectrum = complex(nan(nchan,nfreqoi,ntimeboi),nan(nchan,nfreqoi,ntimeboi));
datspectrum = transpose(fft(transpose([dat repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
for ifreqoi = 1:nfreqoi
  % compute indices that will be used to extracted the requested fft output
  nsamplefreqoi    = taplen(ifreqoi);
  reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi < (ndatsample - (nsamplefreqoi ./2))));
  reqtimeboi       = timeboi(reqtimeboiind);
  
  % compute datspectrum*wavelet, if there are reqtimeboi's that have data
  if ~isempty(reqtimeboi)
    dum = fftshift(transpose(ifft(transpose(datspectrum .* repmat(wltspctrm{ifreqoi},[nchan 1])))),2); % double explicit transpose to speedup fft
    dum = dum .* sqrt(2 ./ fsample);
    spectrum(:,ifreqoi,reqtimeboiind) = dum(:,reqtimeboi);
  end
end

function val = lab_getopt(opt, key, default, emptymeaningful)

% FT_GETOPT gets the value of a specified option from a configuration structure
% or from a cell-array with key-value pairs.
%
% Use as
%   val = lab_getopt(s, key, default, emptymeaningful)
% where the input values are
%   s               = structure or cell-array
%   key             = string
%   default         = any valid MATLAB data type
%   emptymeaningful = boolean value (optional, default = 0)
%
% If the key is present as field in the structure, or as key-value
% pair in the cell-array, the corresponding value will be returned.
% 
% If the key is not present, lab_getopt will return an empty array.
%
% If the key is present but has an empty value, then the emptymeaningful
% flag specifies whether the empty value or the default value should
% be returned. If emptymeaningful==true, then an empty array will be
% returned. If emptymeaningful==false, then the specified default will
% be returned.
%
% See also FT_SETOPT, FT_CHECKOPT

% Copyright (C) 2011, Robert Oostenveld
%
% $Id: lab_getopt.m 5076 2011-12-22 13:40:40Z borreu $

if nargin<3
  default = [];
end

if nargin < 4
  emptymeaningful = 0;
end

if isa(opt, 'struct') || isa(opt, 'config')
  % get the key-value from the structure
  fn = fieldnames(opt);
  if ~any(strcmp(key, fn))
    val = default;
  else
    val = opt.(key);
  end
  
elseif isa(opt, 'cell')
  % get the key-value from the cell-array
  if mod(length(opt),2)
    error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
  end
  
  % the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
  keys = opt(1:2:end);
  vals = opt(2:2:end);
  
  % the following may be faster than cellfun(@ischar, keys)
  valid = false(size(keys));
  for i=1:numel(keys)
    valid(i) = ischar(keys{i});
  end
  
  if ~all(valid)
    error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
  end
  
  hit = find(strcmpi(key, keys));
  if isempty(hit)
    % the requested key was not found
    val = default;
  elseif length(hit)==1
    % the requested key was found
    val = vals{hit};
  else
    error('multiple input arguments with the same name');
  end

elseif isempty(opt)
  % no options are specified, return default
  val = default;
end % isstruct or iscell or isempty

if isempty(val) && ~isempty(default) && ~emptymeaningful
  % use the default value instead of the empty input that was specified:
  % this applies for example if you do functionname('key', []), where
  % the empty is meant to indicate that the user does not know or care
  % what the value is
  val = default;
end
