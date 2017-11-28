function [gm,Labels] = RaguGrandMeans(V,Group,Cond,L1,L2,Gl)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

if nargin < 3
    Cond = (1:size(V,2))';
end


if isempty(Cond)
    Cond = (1:size(V,2))';
end

if isnan(Cond)
    Cond = ones(size(V,2),1);
end

if isempty(Group)
    Group = ones(size(V,1),1);
end

ValidGroup = Group;
ValidCond  = Cond;

ValidGroup(isnan(ValidGroup)) = [];
ValidCond( isnan(ValidCond(:,1)),:) = [];

if size(ValidCond,2) == 1
    ValidCond(:,2) = 1;
    Cond(:,2) = 1;
end

GroupIndex = unique(ValidGroup);
Cond1Index  = unique(ValidCond(:,1));
Cond2Index  = unique(ValidCond(:,2));

nc1 = numel(Cond1Index);
nc2 = numel(Cond2Index);
ngr = numel(GroupIndex);

gm = nan(ngr,nc1*nc2,size(V,3),size(V,4));
for g = 1:ngr
    s_idx = (Group == GroupIndex(g));
    for c1 = 1:nc1
        for c2 = 1:nc2
            c_idx = (     Cond(:,1) == Cond1Index(c1)) & (     Cond(:,2) == Cond2Index(c2));
            gm(g,c2+(c1-1)*nc2,:,:) = mean(mean(V(s_idx,c_idx,:,:),1),2);
        end
    end
end

if nargout > 1
    Labels = cell(ngr,nc1*nc2);
    
    if ngr == 1 && isempty(Gl)
        Gl = {'All'};
    end
    
    % Set all group labels
    for g = 1:ngr
        for i = 1:nc1*nc2
            Labels{g,i} = sprintf('%s:',Gl{GroupIndex(g)});
        end
    end

    if nc1*nc2 == 1 && numel(L1) > 1 && numel(L2) > 1
        for g = 1:ngr
            Labels{g,1} = sprintf('%s All Conditions',Labels{g,1});
        end
    end

    if nc1 > 1 || numel(L1) == 1
        for c1 = 1:nc1
            c1txt = 'Should not happen';
            for i1 = 1:numel(L1)
                if Cond1Index(c1) == L1(i1).Level
                    c1txt = L1(i1).Label;
                end
            end
            
            for g = 1:ngr
                for c2 = 1:nc2
                    Labels{g,c2+(c1-1)*nc2} = sprintf('%s %s',Labels{g,c2+(c1-1)*nc2},c1txt);
                end
            end
        end
    end

    if (nc1 > 1 && nc2 > 1) || ((numel(L1) * numel(L2) > 1) && size(Cond,1) > 1)
        for g = 1:ngr
            for c1 = 1:nc1
                for c2 = 1:nc2
                    Labels{g,c2+(c1-1)*nc2} = sprintf('%s,',Labels{g,c2+(c1-1)*nc2});
                end
            end
        end
    end

    if nc2 > 1
        for c2 = 1:nc2
            c2txt = 'Should not happen';
            for i2 = 1:numel(L2)
                if Cond2Index(c2) == L2(i2).Level
                    c2txt = L2(i2).Label;
                end
            end
            for g = 1:ngr
                for c1 = 1:nc1
                    Labels{g,c2+(c1-1)*nc2} = sprintf('%s ',Labels{g,c2+(c1-1)*nc2},c2txt);
                end
            end
        end
    end
end
