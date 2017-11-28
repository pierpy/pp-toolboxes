 function Randomizer_ShowHitDurationResults(out, flag)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License




if nargin < 2
    flag = 0;
end

figure(105+flag);
clf



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

mx = [];
for i = 1:2
    for j = 1:nc
        if (flag == 0)
            idx = find(out.PHitDuration{i,j} == 0);
        else
            idx = find(out.PHitDurationGFP{i,j} == 0);
        end
        mx = [mx,idx(1)];
    end
end

mx = max(mx);

f = out.DeltaX;
if ng == 2
    for i = 1:2
        for j = 1:nc
            if (i * j > 1)
                subplot(nc,ng,ng*j+i-2);
                if flag == 0
                    p = 1-(cumsum(out.PHitDuration{i,j})/sum(out.PHitDuration{i,j}));
                    crit = out.CritDuration(i,j);
                else
                    p = 1-(cumsum(out.PHitDurationGFP{i,j})/sum(out.PHitDurationGFP{i,j}));
                    crit = out.CritDurationGFP(i,j);
                end
                plot([0 mx*f],[out.Threshold out.Threshold],'-r');    
                hold on
                plot([crit*f crit*f],[0 1],'-g');
                plot((1:mx)*f,p(1:mx),'-k');
                axis([0 mx*f 0 1]);
                title(char(t{ng*j+i-2}));
                set(text(crit*f,0.5,[num2str(crit*out.DeltaX) out.txtX]),'HorizontalAlignment','center','BackgroundColor',[0 1 0]);
            end
        end
    end
else
    for j = 2:nc
        subplot(nc-1,1,j-1);
        if flag == 0
            p = 1-(cumsum(out.PHitDuration{1,j})/sum(out.PHitDuration{1,j}));
            crit = out.CritDuration(1,j);
        else
            p = 1-(cumsum(out.PHitDurationGFP{1,j})/sum(out.PHitDurationGFP{1,j}));
            crit = out.CritDurationGFP(1,j);
        end
        plot([0 mx*f],[out.Threshold out.Threshold],'-r');
        hold on
        plot([crit*f crit*f],[0 1],'-g');
        plot((1:mx)*f,p(1:mx),'-k');
        axis([0 mx*f 0 1]);
        title(char(t{2*j-1}));
        set(text(crit*f,0.5,[num2str(crit*out.DeltaX) out.txtX]),'HorizontalAlignment','center','BackgroundColor',[0 1 0]);
    end
end