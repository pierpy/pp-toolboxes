function [Vectors,Values,Psi] = itab_pcevectors(A, numvecs, verbose)
%PC_EVECTORS Get the top numvecs eigenvectors of the covariance matrix
%            of A, using Turk and Pentland's trick for numrows >> numcols
%            Returns the eigenvectors as the colums of Vectors and a
%            vector of ALL the eigenvectors in Values.
% 
% This is a slightly modified version of pr_evectors by Matthew Dailey
% Released under Creative Commons License, see:
% http://www.cs.ait.ac.th/~mdailey/matlab/pc_evectors.m
% Matthew Dailey 2000

% Check arguments
nexamp = size(A,2);
if nargin < 2
  numvecs = [];
end
if isempty(numvecs)
  numvecs = nexamp-1;
end
if nargin < 3
  verbose = 1;
end
% Now compute the eigenvectors of the covariance matrix, using
% a little trick from Turk and Pentland 1991

% Compute the "average" vector
% mean(A) gives you a row vector containing the mean of each column of A

if verbose
  fprintf(1,'Computing average vector and vector differences from avg...\n');
end
Psi = mean(A')';

% Compute difference with average for each vector
A = A - repmat(Psi, 1, nexamp);

% Get the patternwise (nexamp x nexamp) covariance matrix
if verbose
  fprintf(1,'Calculating L=A''A\n');
end
% L = spm_atranspa(A);
L = A''*A;
% Get the eigenvectors (columns of Vectors) and eigenvalues (diag of Values)
if verbose
  fprintf(1,'Calculating eigenvectors of L...\n');
end
[Vectors,Values] = eig(L);

% Sort the vectors/values according to size of eigenvalue
if verbose
  fprintf(1,'Sorting evectors/values...\n');
end
[Vectors,Values] = itab_sortempp(Vectors,Values);

% Convert the eigenvectors of A'*A into eigenvectors of A*A'
if verbose
  fprintf(1,'Computing eigenvectors of the real covariance matrix..\n');
end
Vectors = A*Vectors;

% Get the eigenvalues out of the diagonal matrix and
% normalize them so the evalues are specifically for cov(A'), not A*A'.

Values = diag(Values);
Values = Values / (nexamp-1);

% Normalize Vectors to unit length, kill vectors corr. to tiny evalues

num_good = 0;
for i = 1:nexamp
  Vectors(:,i) = Vectors(:,i)/norm(Vectors(:,i));
  if Values(i) < 0.00001
    % Set the vector to the 0 vector; set the value to 0.
    Values(i) = 0;
    Vectors(:,i) = zeros(size(Vectors,1),1);
  else
    num_good = num_good + 1;
  end;
end;
if (numvecs > num_good)
  if verbose
    fprintf(1,'Warning: numvecs is %d; only %d exist.\n',numvecs, ...
	    num_good);
  end
  numvecs = num_good;
end;
Vectors = Vectors(:,1:numvecs);
return
