function prs = degreecorrelation(A) 

%calculates Pearson degree correlation of A
%modified 9-14-06
A=sortbyk(A); 
[rows,colms]=size(A); 
k=kvec(A); 
n=length(find(k~=0)); 
won=[ones(n,1); zeros(rows-n,1)]; 
%k=won'*A;
ksum=won'*k'; 
ksqsum=k*k'; 
xbar=ksqsum/ksum; 
num=(k-won'*xbar)*A*(k'-xbar*won); 
kkk=(k'-xbar*won).*(k'.^.5); 
denom=kkk'*kkk; 
prs=num/denom; 