% rihaczek1 -- Compute samples of the type I Rihaczek distribution.
%
%  Usage
%    [tfd, t, f] = rihaczek1(x, fs, nfreq)
%
%  Inputs
%    x     signal vector
%    fs    sampling frequency of x (optional, default is 1 sample/second)
%    nfreq number of samples to compute in frequency (optional, default
%          is the length of x)
%
%  Outputs
%    tfd  matrix containing the Rihaczek distribution of signal x.  If x has
%         length N, then tfd will be nfreq by N. (optional)
%    t    vector of sampling times (optional)
%    f    vector of frequency values (optional)
%
% If no output arguments are specified, then the Rihaczek distribution is 
% displayed using ptfd(tfd, t, f).

% Copyright (C) -- see DiscreteTFDs/Copyright

% specify defaults
%
% Modifications to limit memory use: F. Hatz 2012

function [tfd, t, f] = lab_rihaczek(x, fs, nfreq)

x = x(:);
if size(x,2) > 1;
    x = x.';
end
N = length(x);

if ~exist('nfreq','var')
    nfreq = N;
end
if ~exist('fs','var')
    fs = 1;
end

X = fft(x,nfreq);
X = conj(X);
x = x.';
factor = -1i*[0:2*pi/nfreq:2*pi*(1-1/nfreq)]';
M = ceil(nfreq/2);
tfd = (X * x);
for  i = 1:size(tfd,2)
    tmp = tfd(:,i) .* exp(factor * (i-1));
    tfd(1:(end-M),i) = tmp(M+1:end,1);
    tfd((end-M+1):end,i) = tmp(1:M,1);
end
clearvars factor M tmp X

t = 1/fs * (0:N-1);
f = -fs/2:fs/nfreq:fs/2;
f = f(1:nfreq);

if (nargout == 0)
  ptfd(tfd, t, f);
  clear tfd
end