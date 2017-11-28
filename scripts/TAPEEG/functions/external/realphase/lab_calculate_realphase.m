function phase = lab_calculate_realphase(data,header,cfgfilt)

if exist('cfgfilt','var') & ~isempty(cfgfilt)
    data = lab_filter(data,header,cfgfilt,'novrb');
end
    
phase = zeros(header.numdatachannels,size(data,2));

for nsig = 1:header.numdatachannels
    if max(data(nsig,:)) == 0 & min(data(nsig,:)) == 0
        channels(nsig) = 0; %#ok<AGROW>
    else
        channels(nsig) = 1; %#ok<AGROW>
    end
end
channels = find(channels == 1);

for nsig = channels
    theta= hilbphase(data(nsig,:)');
    [theta,~,sigma] = fbtransf(theta,80,0.05,100);
    phase(nsig,:) = gth2phi(theta,sigma);
end