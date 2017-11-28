function Randomizer_ShowHitCountResults(out,flag)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

figure(104);
clf

if (nargin == 1)
    flag = 0;
end

if flag == 0
    OutData = out.TanovaHits;
    OutP    = out.PHitCount;
end

if (flag == 1)
    OutData = out.TanovaMax;
    OutP    = out.PMaxStat;
end

if (flag == 2)
    OutData = out.TanovaAUC;
    OutP    = out.PAUCCount;
end

if (flag == 3)
    OutData = out.GFPHits;
    OutP    = out.PHitCountGFP;
end

if (flag == 4)
    OutData = out.GFPAUC;
    OutP    = out.PAUCCountGFP;
end


DoF1    = (numel(unique(out.Design(:,1)))> 1);
DoF2    = (numel(unique(out.Design(:,2)))> 1);
if out.TwoFactors == 0
    DoF2 = 0;
end

if DoF1 && DoF2
    nc = 4;
elseif DoF1 && ~DoF2
    nc = 2;
else
    nc = 1;
end

if numel(unique(out.IndFeature(~isnan(out.IndFeature))))> 1
    ng = 2;
else
    ng = 1;
end

t = {'' 'Group main effect' out.strF1 [out.strF1 '* Group'] out.strF2 [out.strF2 ' * Group'] [out.strF1 ' * ' out.strF2] [out.strF1 ' * ' out.strF2 ' * Group']};

mx = max(OutData(:));

if ng == 2
    for i = 1:2
        for j = 1:nc
            if (i * j > 1)
                subplot(nc,ng,ng*j+i-2);
                [hst,x] = hist(squeeze(OutData(i,j,:)),20);
                h = bar(x,hst);
                set(h,'BarWidth',1,'EdgeColor',[0 0 0],'FaceColor',[0.3 0.3 0.3]);
                set(gca,'XLim',[0 mx]);
                hold on
                plot(squeeze(OutData(i,j,1)),0,'or','MarkerSize',10,'MarkerFaceColor','r');
                title([char(t{ng*j+i-2}) ': p=' num2str(OutP(i,j))]);
            end
        end
    end
else
    for j = 2:nc
        subplot(nc-1,1,j-1);
        [hst,x] = hist(squeeze(OutData(1,j,:)),20);
        h = bar(x,hst);
        set(h,'BarWidth',1,'EdgeColor',[0 0 0],'FaceColor',[0.3 0.3 0.3]);
        set(gca,'XLim',[0 mx]);
        hold on
        plot(squeeze(OutData(1,j,1)),0,'or','MarkerSize',10,'MarkerFaceColor','r');
        title([char(t{2*j-1}) ': p=' num2str(OutP(1,j))]);
    end
end