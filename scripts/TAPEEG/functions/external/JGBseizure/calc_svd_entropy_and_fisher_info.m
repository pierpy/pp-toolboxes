function [SVD_entropy,Fisher_information] = calc_svd_entropy_and_fisher_info(data)

SVD_entropy = zeros(1,size(data,1));
Fisher_information = zeros(1,size(data,1));  
for C = 1:size(data,1)
    %delay vectors y(n)
    tau=1;
    dE=20;
    y=zeros(dE,size(data,2)-dE);
    for n=1:size(data,2)-dE
        y(:,n)= data(C,n:(n+(dE-1)*tau));
    end
    
    singular_values=svd(y);
    norm_singular_values=singular_values/sum(singular_values);
    
    SVD_entropy(1,C)=-sum(norm_singular_values.*log(norm_singular_values));
    
    Fisher_information(1,C) = sum(((norm_singular_values(2:end)-...
        norm_singular_values(1:end-1)).^2)./norm_singular_values(1:end-1));
end
