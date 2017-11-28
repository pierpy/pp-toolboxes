function [connections,channels] = lab_split_connections(connections,numchans)
    
if (numchans^2/2 - numchans/2) == size(connections,2)
    channels = zeros(size(connections,1),numchans);
    return
elseif (numchans^2/2 + numchans/2) == size(connections,2)
    [~,ChanIdx] = intersect(find(tril(true(numchans,numchans))),find(eye(numchans)));
    ConnIdx = setdiff(1:size(connections,2),ChanIdx);
    channels = connections(:,ChanIdx);
    connections = connections(:,ConnIdx);
else
    channels = [];
    disp('Extracting channels from connections not possible, wrong input-data')
end