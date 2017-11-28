function [pnt, tri] = mysphere(N)
% This is a copy of MSPHERE without the confusing output message
% Returns a triangulated sphere with approximately M vertices
% that are nicely distributed over the sphere. The vertices are aligned
% along equally spaced horizontal contours according to an algorithm of
% Dave Russel.
% 
% Use as
%  [pnt, tri] = msphere(M)
%
% See also SPHERE, NSPHERE, ICOSAHEDRON, REFINE
% Copyright (C) 1994, Dave Rusin 

storeM    = [];
storelen  = [];
increaseM = 0;
while (1)

  % put a single vertex at the top
  phi = [0];
  th  = [0];

  M = round((pi/4)*sqrt(N)) + increaseM;
  for k=1:M
    newphi = (k/M)*pi;
    Q = round(2*M*sin(newphi));
    for j=1:Q
      phi(end+1) = newphi;
      th(end+1)  = (j/Q)*2*pi;
      % in case of even number of contours
      if mod(M,2) & k>(M/2)
        th(end) = th(end) + pi/Q;
      end
    end
  end

  % put a single vertex at the bottom
  phi(end+1) = [pi];
  th(end+1)  = [0];

  % store this vertex packing
  storeM(end+1).th  = th;
  storeM(end  ).phi = phi;
  storelen(end+1) = length(phi);
  if storelen(end)>N
    break;
  else
    increaseM = increaseM+1;
    % fprintf('increasing M by %d\n', increaseM);
  end
end

% take the vertex packing that most closely matches the requirement
[m, i] = min(abs(storelen-N));
th  = storeM(i).th;
phi = storeM(i).phi;

% convert from spherical to cartehsian coordinates
[x, y, z] = sph2cart(th, pi/2-phi, 1);
pnt = [x' y' z'];
tri = convhulln(pnt);