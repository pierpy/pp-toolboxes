function res = TestForSingleGroups(handles)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

d = get(handles.output,'UserData');
res = 0;

if d.ContBetween == 1
    return;
end
Groups = unique(d.IndFeature(~isnan(d.IndFeature)));
for i = 1:numel(Groups)
    if numel(find(d.IndFeature == Groups(i))) < 2
        res = 1;
    end
end

if res == 1 && d.ContBetween == 0
    errordlg('Some groups contain only one subject, please check before continuing','Randomizer ERROR');
end