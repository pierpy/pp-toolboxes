function [OutFeature] = MSFeatures(MSClass,n,MSFit,gfp, reject_small_segments, fs)


OutFeature = struct;
on   = nan(1,n);
off  = nan(1,n);
dur  = nan(1,n);
auc  = nan(1,n);
cog  = nan(1,n);
mgfp = nan(1,n);


%MSClassR = reshape(MSClass,size(MSClass,1)*size(MSClass,2),size(MSClass,3));

%ndiff = sum(var(MSClassR)> 0);

%h = hist(MSClassR,1:n) / size(MSClassR,1);

% mscomplexity = h.*log(h);
% mscomplexity(isnan(mscomplexity)) = 0;
% mscomplexity = -sum(mscomplexity);
% ndiff = sum(mscomplexity);

x = 1:size(MSClass,2);
        for m = 1:n
            o = find((MSClass == m),1,'first');
            f = find((MSClass == m),1,'last');
            if ~isempty(gfp)
                mgfp(1,m) = mean(gfp(1,MSClass == m));
            end
            tdur(1,m) = sum(MSClass == m);
            tcov(1,m) = tdur(1, m)/length(MSClass);
            occ(1,m) = length(find(MSClass == m));
            if ~isempty(o)
                on(1,m)  = o;
            end
            if ~isempty(f)
                off(1,m) = f;
            end
            MSGFP = MSFit;
            MSGFP(MSClass ~= m) = 0;
            auc(1,m) = sum(MSGFP);
            cog(1,m) = sum(MSGFP .* x) / auc(1,m);
%            auc(g,c,m) = sum(MSFit(g,c,MSClass(g,c,:) == m));
        end

        %durata media
%         count=0;
%         kk=1;
%         durata=[];
%         for m=1:n
%             for k = 1:size(MSClass,2)-1
%                 
%                 if MSClass(k) == m 
%                     count=count+1; 
%                     if (MSClass(k+1) ~= MSClass(k))||k==size(MSClass,2)-1
%                          durata(kk)=count;
%                          count=0;
%                          kk=kk+1;
%                     end
%                 else
% %                     durata(kk)=count; 
%                 end
%             
%             end
%             
%             if reject_small_segments
%                 iii1 = find(durata == 1);
%                  iii2 = find(durata == 2);
%                  iii3 = find(durata == 3);
%                 iii=[iii1,iii2,iii3];
% 
%                 durata(iii)=[];
%             end
%             if ~isempty(durata)
%             meandur(m)=mean(durata);
%             mediandur(m)=median(durata);
%             geomedur(m)=geomean(durata);
%             maxdur(m)=max(durata);
%             mindur(m)=min(durata);
% %             clearvars durata
%             count=0;
%             kk=1;
%             end
%         end
        
        %durata media altro metodo
        splitted = SplitVec(MSClass);
        for m = 1:n
            index = cell2mat(cellfun(@(x) any(x == m), splitted, 'UniformOutput', 0));
            dur = splitted(index);
            if reject_small_segments
                dursmall = dur(cell2mat(cellfun(@(x) length(x)>=6, dur, 'UniformOutput', 0)));           
            else
                dursmall = dur;
            end
            dur_sum=0;
            for jj = 1:length(dur)
                dur_sum = dur_sum+length(dur{1,jj});
            end
            meandurNew(m) = mean(cell2mat(cellfun(@(x) mean(length(x)), dursmall, 'UniformOutput', 0)));
            freq_of_occ(m)= length(dursmall)/(length(MSClass)/fs) ;
%           freq_of_occ(m)= length(dur) * (length(MSClass)*(1000/125/length(MSClass) )/(length(MSClass)/125));
%           maxdurNew(m) = max(cell2mat(cellfun(@(x) length(x), dursmall, 'UniformOutput', 0)));
%           mindurNew(m) = min(cell2mat(cellfun(@(x) length(x), dursmall, 'UniformOutput', 0)));
%           geomNew(m) = geomean(cell2mat(cellfun(@(x) length(x), dursmall, 'UniformOutput', 0)));
%           medianNew(m) = median(cell2mat(cellfun(@(x) length(x), dursmall, 'UniformOutput', 0)));
        end
% if nargin < 5
%     NoCenter = true;
% end
% 
% if NoCenter == true
%     for m = 1:n
%         on(1,m)  = on(1,m)  - mean(on(1,m));
%         off(1,m) = off(1,m) - mean(off(1,m));
%         dur(1,m) = dur(1,m) - mean(dur(1,m));
%         auc(1,m) = auc(1,m) - mean(auc(1,m));
%         cog(1,m) = cog(1,m) - mean(cog(1,m));
%         mgfp(1,m) = mgfp(1,m) - mean(mgfp(1,m));
%     end
% end

OutFeature.on = on;
OutFeature.of = off;
OutFeature.totDur = tdur;
OutFeature.auc = auc;
OutFeature.cog = cog;
OutFeature.mgfp = mgfp;
% OutFeature.meandur = meandur;
% OutFeature.maxdur = maxdur;
% OutFeature.mindur = mindur;
% OutFeature.mediandur = mediandur;
% OutFeature.geomedur = geomedur;
OutFeature.meandurNew = meandurNew;
%OutFeature.maxdurNew = maxdurNew;
% OutFeature.mindurNew = mindurNew;
% OutFeature.geomNew = geomNew;
% OutFeature.medianNew = medianNew;
OutFeature.occ = occ;
OutFeature.freq_of_occ = freq_of_occ;
OutFeature.tcov=tcov;
