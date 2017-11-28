% Calculate correct phase shifts
%
% phase = lab_calc_phaseshifts(phase);
%
% Input = phase (channels x timepoints)

function phase = lab_calc_phaseshifts(phase)

NumChans = size(phase,1);
Ca = fcorrelation(phase);
Hold = [];
P1 = zeros(1,NumChans-1);
P2 = P1;
for i = 1:NumChans - 1
    if isempty(Hold)
        Cx = Ca;
    else
        Cx = zeros(NumChans,NumChans);
        Cx(:,Hold) = Ca(:,Hold);
    end
    P1(i) = find(max(Cx,[],1) == max(Cx(:)),1,'first');
    [~,Idx] = sort(Cx(:,P1(i)),'descend');
    P2(i) = Idx(1);
    Hold = union(Hold,[P1(i) P2(i)]);
    Ca(P1(i),P2(i)) = 0;
    Ca(P2(i),P1(i)) = 0;
end

for i = 1:length(P1)
    shift = 1;
    counter = 1;
    while abs(shift) > 0.001 & counter < 20
        [phase(P2(i),:),shift] = estimate_shift(phase(P2(i),:),phase(P1(i),:));
        counter = counter + 1;
    end
end

end

function [phase,shift] = estimate_shift(phase,phaseR)
    % Estimate phase shifts
    phasediff = abs(phase - phaseR);
    phasecorr = angle(sum(phasediff.*exp(1i*phase),2)) - pi;
    phasecorr(phasecorr < -pi) = phasecorr(phasecorr < -pi) + 2*pi;
    phasecorr(phasecorr > pi) = phasecorr(phasecorr  > pi) - 2*pi;
    phasecorr(sum(phasediff,2) == 0) = 0;
    phasecorr2 = angle(sum(phasediff.*exp(1i*phaseR),2)) - pi;
    phasecorr2(phasecorr2 < -pi) = phasecorr2(phasecorr2 < -pi) + 2*pi;
    phasecorr2(phasecorr2 > pi) = phasecorr2(phasecorr2  > pi) - 2*pi;
    phasecorr2(sum(phasediff,2) == 0) = 0;
    shift = (phasecorr - phasecorr2) / 2;
    
    phase = phase - shift;
    phase(phase < -pi) = phase(phase < -pi) + 2*pi;
    phase(phase > pi) = phase(phase  > pi) - 2*pi;
end

function Ca = fcorrelation(phase)
    NumChans = size(phase,1);
    Ca = zeros(NumChans,NumChans);
    for i = 1:NumChans
        for j = i+1:NumChans
             A = sum(cos(phase(i,:)).*cos(phase(j,:)));  
             B = sum(sin(phase(i,:)).*sin(phase(j,:)));
             C = sum(cos(phase(i,:)).*sin(phase(j,:)));
             D = sum(sin(phase(i,:)).*cos(phase(j,:)));
             E = sum(cos(2*phase(i,:)));
             F = sum(sin(2*phase(i,:)));
             G = sum(cos(2*phase(j,:)));
             H = sum(sin(2*phase(j,:)));
             n = length(phase(i,:)) ;
             denom = sqrt(((n*n)-(E*E)-(F*F))*((n*n)-(G*G)-(H*H)));
             Ca(j,i) = 4*((A*B)-(C*D))/denom;
        end
    end
    Ca = Ca + Ca';
end
