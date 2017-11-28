function [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices(cfg, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that converts a freq structure into Cf, Cr and Pr
% this is used in sourecanalysis
%
% This function returns data matrices with a channel order that is consistent
% with the original channel order in the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2004-2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: prepare_freq_matrices.m 7451 2013-02-07 12:05:28Z jansch $

% set the defaults
if ~isfield(cfg, 'dicsfix'), cfg.dicsfix = 'yes'; end
if ~isfield(cfg, 'quickflag'), cfg.quickflag = 0; end
if ~isfield(cfg, 'refchan'), cfg.refchan = []; end

quickflag = cfg.quickflag==1;

Cf = [];
Cr = [];
Pr = [];

% the old freqanalysis (up to revision 1.9) used sgn and sgncmb
% for backward compatibility rename these
if isfield(freq, 'sgn')
  freq.label = freq.sgn;
  freq = rmfield(freq, 'sgn');
end
if isfield(freq, 'sgncmb')
  freq.labelcmb = freq.sgncmb;
  freq = rmfield(freq, 'sgncmb');
end

% select the latency of interest for time-frequency data
if strcmp(freq.dimord, 'chan_freq_time')
  tbin = nearest(freq.time, cfg.latency);
  fprintf('selecting timeslice %d\n', tbin);
  freq.time = freq.time(tbin);
  % remove all other latencies from the data structure and reduce the number of dimensions
  if isfield(freq, 'powspctrm'),     freq.powspctrm     = squeeze(freq.powspctrm(:,:,tbin));     end;
  if isfield(freq, 'crsspctrm'),     freq.crsspctrm     = squeeze(freq.crsspctrm(:,:,tbin));     end;
  if isfield(freq, 'fourierspctrm'), freq.fourierspctrm = squeeze(freq.fourierspctrm(:,:,tbin)); end;
  freq.dimord = freq.dimord(1:(end-5));  % remove the '_time' part
elseif strcmp(freq.dimord, 'rpt_chan_freq_time') || strcmp(freq.dimord, 'rpttap_chan_freq_time')
  tbin = nearest(freq.time, cfg.latency);
  fprintf('selecting timeslice %d\n', tbin);
  freq.time = freq.time(tbin);
  % remove all other latencies from the data structure and reduce the number of dimensions
  if isfield(freq, 'powspctrm'),    freq.powspctrm     = squeeze(freq.powspctrm(:,:,:,tbin));      end;
  if isfield(freq, 'crsspctrm'),    freq.crsspctrm     = squeeze(freq.crsspctrm(:,:,:,tbin));      end;
  if isfield(freq, 'fourierspctrm') freq.fourierspctrm = squeeze(freq.fourierspctrm(:,:,:,tbin));  end;
  freq.dimord = freq.dimord(1:(end-5));  % remove the '_time' part
else
  tbin = [];
end

% the time-frequency latency has already been squeezed away (see above)
if strcmp(freq.dimord, 'chan_freq')
  Ntrials = 1;
elseif strcmp(freq.dimord, 'rpt_chan_freq')
  Ntrials = size(freq.cumtapcnt,1);
elseif strcmp(freq.dimord, 'rpttap_chan_freq')
  Ntrials = size(freq.cumtapcnt,1);
else
  error('unrecognized dimord for frequency data');
end

% find the frequency of interest
fbin = nearest(freq.freq, cfg.frequency);

if isfield(freq, 'powspctrm') && isfield(freq, 'crsspctrm')
  % use the power and cross spectrum and construct a square matrix

  % find the index of each sensor channel into powspctrm
  % keep the channel order of the cfg
  [dum, powspctrmindx] = match_str(cfg.channel, freq.label);
  Nchans = length(powspctrmindx);

  % find the index of each sensor channel combination into crsspctrm
  % keep the channel order of the cfg
  crsspctrmindx = zeros(Nchans);
  for sgncmblop=1:size(freq.labelcmb,1)
    ch1 = find(strcmp(cfg.channel, freq.labelcmb(sgncmblop,1)));
    ch2 = find(strcmp(cfg.channel, freq.labelcmb(sgncmblop,2)));
    if ~isempty(ch1) && ~isempty(ch2)
      % this square matrix contains the indices into the signal combinations
      crsspctrmindx(ch1,ch2) = sgncmblop;
    end
  end

  % this complex rearrangement of channel indices transforms the CSDs into a square matrix
  if strcmp(freq.dimord, 'chan_freq')
    % FIXME this fails in case dimord=rpt_chan_freq and only 1 trial
    Cf = complex(nan(Nchans,Nchans));
    % first use the complex conjugate for all reversed signal combinations
    Cf(find(crsspctrmindx)) = freq.crsspctrm(crsspctrmindx(find(crsspctrmindx)), fbin);
    Cf = ctranspose(Cf);
    % and then get get the csd for all signal combinations
    Cf(find(crsspctrmindx)) = freq.crsspctrm(crsspctrmindx(find(crsspctrmindx)), fbin);
    % put the power on the diagonal
    Cf(find(eye(Nchans))) = freq.powspctrm(powspctrmindx, fbin);
  else
    Cf  = complex(nan(Ntrials,Nchans,Nchans));
    tmp = complex(nan(Nchans,Nchans));
    for trial=1:Ntrials
      % first use the complex conjugate for all signal combinations reversed
      tmp(find(crsspctrmindx)) = freq.crsspctrm(trial, crsspctrmindx(find(crsspctrmindx)), fbin);
      tmp = ctranspose(tmp);
      % and then get get the csd for all signal combinations
      tmp(find(crsspctrmindx)) = freq.crsspctrm(trial, crsspctrmindx(find(crsspctrmindx)), fbin);
      % put the power on the diagonal
      tmp(find(eye(Nchans))) = freq.powspctrm(trial, powspctrmindx, fbin);
      Cf(trial,:,:) = tmp;
    end
  end

  % do a sanity check on the cross-spectral-density matrix
  if any(isnan(Cf(:)))
    error('The cross-spectral-density matrix is not complete');
  end

  if isfield(cfg, 'refchan') && ~isempty(cfg.refchan)
    % contruct the cross-spectral-density vector of the reference channel with all MEG channels
    tmpindx = match_str(freq.labelcmb(:,1), cfg.refchan);
    refindx = match_str(freq.labelcmb(tmpindx,2), cfg.channel);
    refindx = tmpindx(refindx);
    flipref = 0;
    if isempty(refindx)
      % first look in the second column, then in the first
      tmpindx = match_str(freq.labelcmb(:,2), cfg.refchan);
      refindx = match_str(freq.labelcmb(tmpindx,1), cfg.channel);
      refindx = tmpindx(refindx);
      flipref = 1;
    end
    if length(refindx)~=Nchans
      error('The cross-spectral-density with the reference channel is not complete');
    end
    if Ntrials==1
      Cr = freq.crsspctrm(refindx, fbin);
    else
      for trial=1:Ntrials
        Cr(trial,:) = freq.crsspctrm(trial, refindx, fbin);
      end
    end
    if flipref
      Cr = conj(Cr);
    end
    % obtain the power of the reference channel
    refindx = match_str(freq.label, cfg.refchan);
    if length(refindx)<1
      error('The reference channel was not found in powspctrm');
    elseif length(refindx)>1
      error('Multiple occurences of the reference channel found in powspctrm');
    end
    if Ntrials==1
      Pr = freq.powspctrm(refindx, fbin);
    else
      for trial=1:Ntrials
        Pr(trial) = freq.powspctrm(trial, refindx, fbin);
      end
      Pr = Pr(:);   % ensure that the first dimension contains the trials
    end
  end

  if strcmp(cfg.dicsfix, 'yes')
    Cr = conj(Cr);
  end

else
  fprintf('computing cross-spectrum from fourier\n');
  [dum, powspctrmindx] = match_str(cfg.channel, freq.label);
  % use the fourier spectrum to compute the complete cross spectrum matrix
  Nchans = length(powspctrmindx);
  if strcmp(freq.dimord, 'chan_freq')
    error('incompatible dimord for computing CSD matrix from fourier');
  elseif strcmp(freq.dimord, 'rpt_chan_freq')
    error('incompatible dimord for computing CSD matrix from fourier');
  elseif strcmp(freq.dimord, 'rpttap_chan_freq'),
    if quickflag,
      Ntrials = 1;
    end
    Cf = zeros(Ntrials,Nchans,Nchans);
    refindx = match_str(freq.label, cfg.refchan);
    if ~isempty(refindx)
      Cr = zeros(Ntrials,Nchans,1);
      Pr = zeros(Ntrials,1,1);
    end

    if quickflag,
      ntap = sum(freq.cumtapcnt);
      dat  = transpose(freq.fourierspctrm(:, powspctrmindx, fbin));
      Cf(1,:,:) = (dat * ctranspose(dat)) ./ ntap;
      if ~isempty(refindx)
        ref = transpose(freq.fourierspctrm(:, refindx, fbin));
    Cr(1,:,1) = dat * ctranspose(ref) ./ ntap;
    Pr(1,1,1) = ref * ctranspose(ref) ./ ntap;
      end
    else
      freq.cumtapcnt = freq.cumtapcnt(:)';
      for k=1:Ntrials
        tapbeg = 1 + sum([0 freq.cumtapcnt(1:(k-1))]);
        tapend =     sum([0 freq.cumtapcnt(1:(k  ))]);
        ntap = freq.cumtapcnt(k);
        dat  = transpose(freq.fourierspctrm(tapbeg:tapend, powspctrmindx, fbin));
        Cf(k,:,:) = (dat * ctranspose(dat)) ./ ntap;
        if ~isempty(refindx)
          ref = transpose(freq.fourierspctrm(tapbeg:tapend, refindx, fbin));
          Cr(k,:,1) = dat * ctranspose(ref) ./ ntap;
          Pr(k,1,1) = ref * ctranspose(ref) ./ ntap;
        end
      end
    end
  else
    error('unsupported dimord for fourier cross-spectrum computation');
  end
end

% update the configuration, so that the calling function exactly knows what was selected
if ~isempty(tbin),
  % a single latency was selected in the freq structure
  cfg.latency = freq.time;
else
  cfg.latency = [];
end
cfg.frequency = freq.freq(fbin);
cfg.channel   = freq.label(powspctrmindx);

