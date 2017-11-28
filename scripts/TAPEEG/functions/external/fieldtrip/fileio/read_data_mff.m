function [data] = read_data_mff(filename,hdr)
  chanindx = 1:hdr.nChans;
  endsample = hdr.nSamples*hdr.nTrials;
  begsample = 1;
  binfiles = dir(fullfile(filename, 'signal*.bin'));
  if isempty(binfiles)
     error('FieldTrip:read_mff_header:nobin', ['could not find any signal.bin in ' filename_mff ])
  end
  %determine which channels are in which signal
  for iSig = 1:length(hdr.orig.signal)
    if iSig == 1
      chan2sig_ind(1:hdr.orig.signal(iSig).blockhdr(1).nsignals(1)) = iSig;
    else
      chan2sig_ind(end+1:end+1+hdr.orig.signal(iSig).blockhdr(1).nsignals(1)) = iSig;
    end
  end
  
  for iSig = 1:length(hdr.orig.signal)
      % adjust chanindx to match with current signal
      [dum1, dum2, chanind_sig] = intersect(chanindx, find(chan2sig_ind==iSig));
      if isempty(chanind_sig)
        % no channels requested from current signal
      else
        blockhdr = hdr.orig.signal(iSig).blockhdr;
        signalname = binfiles(iSig).name;
        fullsignalname = fullfile(filename, signalname);

        % the number of samples per block can be different
        % assume that all channels have the same sampling frequency and number of samples per block
        nsamples = zeros(size(blockhdr));
        for i=1:length(blockhdr)
          nsamples(i) = blockhdr(i).nsamples(1);
        end

        cumsamples = cumsum(nsamples);
        begblock = find(begsample<=cumsamples, 1, 'first');
        endblock = find(endsample<=cumsamples, 1, 'first');
        datsig = read_mff_bin(fullsignalname, begblock, endblock, chanind_sig);

        % concatenate in a matrix
        if exist('dat', 'var')
          dat{length(dat)+1} = cell2mat(datsig(:,:));
        else
          dat{1} = cell2mat(datsig(:,:));
        end
        % select the desired samples from the concatenated blocks
        if begblock==1
          prevsamples = 0;
        else
          prevsamples = cumsamples(begblock-1);
        end
        begsel = begsample-prevsamples;
        endsel = endsample-prevsamples;
        dat{end} = dat{end}(:,begsel:endsel);
      end
  end
  if length(hdr.orig.signal)== 2;
    tmp = dat{1,2};
    tmp = tmp((size(tmp,1)-1),:);
    data = cat(1,dat{1,1},tmp);
  else
    data = dat{1,1};
  end