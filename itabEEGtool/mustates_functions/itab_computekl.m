function [kl] = computekl(mm)
    for i=1:size(mm,2)-1
        dd(i)=mm(i)-mm(i+1);
    end
    for k=1:size(dd,2)
        if (k==1 || dd(k-1)<=0 || dd(k-1)<dd(k))
            kl(k)=0;
        else
            kl(k)=(dd(k-1)-dd(k))/mm(k-1);
        end
    end
end