function phi=gth2phi(theta,sigma)

% This function computes the true phase phi from the protophase theta,
% using the transformation function sigma, which is given on a grid
%
% Form of call: phi=gth2phi(theta,sigma)
%
pi2=pi+pi; ngrid=length(sigma); ng1=ngrid-1;
argsi=(0:ng1)*pi2/ng1;
sigma=sigma(:); sigma=sigma';        
sigext=[sigma(1:end-1) sigma sigma(2:end)];     
argext=[argsi(1:end-1)-pi2  argsi argsi(2:end)+pi2];
sigext=cumtrapz(sigext)*pi2/ng1;
phi=interp1(argext,sigext,theta,'spline');
phi=mod(phi,pi2); 

end
