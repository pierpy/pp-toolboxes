function y=NLEO(x)
  
y=zeros(1,length(x)-3);
    for n=4:length(x)
        y(n-3)=(x(n-1)*x(n-2)) - (x(n)*x(n-3));
    end
      y=abs(y);%???
%     figure;plot(y)