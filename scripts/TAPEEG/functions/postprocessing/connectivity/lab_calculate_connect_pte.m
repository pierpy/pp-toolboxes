% Calculate phase shifts and p-values
%
% phase = lab_calculate_connect_pte(phase,phaseT);
%
% Input = phaseT (channels x timepoints) lagged
%         phase (channels x timepoints)


function [PTE,IC] = lab_calculate_connect_pte(phase,phaseT,Pbar)
    
if ~exist('Pbar','var')
    Pbar = false;
end
NumChans = size(phaseT,1);
NumTF = size(phaseT,2);
IC =zeros(NumChans,NumChans);
PTE =zeros(NumChans,NumChans);

if Pbar == false
    fprintf('   calculate PTE:');
else
    progressbar('Calculate PTE',0);
end

Tidx = floor(NumChans/20);
for Idx = 1:NumChans
    for Idx2 = 1:NumChans
        phasetmp = cat(1,phase(Idx,:),phaseT(Idx,:),phase(Idx2,:),phaseT(Idx2,:));
        phasetmp = lab_calc_phaseshifts(phasetmp);
        phasetmp = transform2latentspace(phasetmp);
        
        % Calculate multiinformational term 1
        C1 = cat(1,phasetmp(1,:),phasetmp(2,:),phasetmp(3,:),phasetmp(4,:));
        C1 = det(C1*C1'*NumTF^-1);
        
        D1 = cat(1,phasetmp(1,:),phasetmp(2,:),phasetmp(4,:));
        D1 = det(D1*D1'*NumTF^-1);
        
        % calculte multiinformational term 2
        C2 = cat(1,phasetmp(2,:),phasetmp(4,:));
        C2 = det(C2*C2'*NumTF^-1);
        
        D2 = det(phasetmp(2,:)*phasetmp(2,:)'*NumTF^-1);
        
        % Calculate multiinformational term 3
        C3 = cat(1,phasetmp(2,:),phasetmp(3,:),phasetmp(4,:));
        C3 = det(C3*C3'*NumTF^-1);
        
        D3 = cat(1,phasetmp(1,:),phasetmp(2,:));
        D3 = det(D3*D3'*NumTF^-1);
        
        % Calculate multiinformational term 4
        C4 = cat(1,phasetmp(1,:),phasetmp(2,:),phasetmp(4,:));
        C4 = det(C4*C4'*NumTF^-1);
        
        D4 = cat(1,phasetmp(2,:),phasetmp(4,:));
        D4 = det(D4*D4'*NumTF^-1);

        % Calc PTE
        PTE(Idx2,Idx) = 0.5*(-D1 - D2 + D3 + D4);
        
        % Calc IC
        IC(Idx2,Idx) = 0.5*(-C1 - C2 + C3 + C4);
        
        if Pbar == true & mod(Idx2,Tidx) == 0
            progressbar([],Idx2/NumChans-0.0001);
        end
    end
    if Pbar == false & mod(Idx,Tidx) == 0
        fprintf('.');
    else
        progressbar(Idx/NumChans,[]);
    end
end
if Pbar == false
    disp(':');
else
    progressbar(1);
end
% IC = (IC + IC') ./ 2;
PTE = PTE ./ (PTE + PTE' + IC);
PTE(1:NumChans+1:end) = 0;
    
end

function phase = transform2latentspace(phase)
    % rank (0 to 1)
    NumChans = size(phase,1);
    for i = 1:NumChans
        [~,phase(i,:)] = sort(phase(i,:));
    end
    phase = phase / (size(phase,2)+1);
    
    % transform to latent space
    phase = norminv(phase,0,1);
end