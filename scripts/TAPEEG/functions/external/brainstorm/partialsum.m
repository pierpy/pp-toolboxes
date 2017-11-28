function x = partialsum(x, pct)
    % Sort
    x = sort(x); 
    % Sum lowest values
    x = sum(x(1:ceil(length(x)*pct))); 
end