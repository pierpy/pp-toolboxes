function dirin = dirin(N1,N2, omeg1,omeg2)
% 
%
% Given the norms of the coupling functions of two coupled 
% oscillators and their autonomous frequencies, this 
% functions returns the directionality index.
% In case of symmetrical coupling dirin=0 holds, while purely 
% uni-directional coupling yields dirin=1 or dirin=-1.  
%
% Form of call: dirin = dirin(N1,N2,omeg1,omeg2)
% Input:        N1,N2 : norms of the coupling functions
%               omeg1,omeg2: frequencies
%
c1=N1/omeg1; % strength of the external contribution to the phase 
c2=N2/omeg2; % dynamics normalized by the natural frequency
% to check for presence of interaction

dirin= (c1-c2) / (c1+c2); %  Directionality index
end
