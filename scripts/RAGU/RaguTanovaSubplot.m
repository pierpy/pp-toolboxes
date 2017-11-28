function h = RaguTanovaSubplot(nr,nc,r,c)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

r_width  = 0.5 / nc;
r_height = 1.0 / nr;

h = subplot('Position',[(0.05 + (c-1) * r_width) (0.05 + (nr-r) * r_height) (r_width - 0.04) (r_height - 0.1) ]);
cla(h);
set(h,'Tag','subplot');
