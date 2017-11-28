function lab_convert_result2table
    
[data,header,~,~,settings] = lab_read_statistics([],-1,0,1,0,1);

if ~isfield(settings,'clustervars')
    settings.clustervars = 1;
end
if ~isfield(settings,'clustervars2')
    settings.clustervars2 = 1;
end
Nmeasures = floor(size(data,2) / (settings.clustervars*settings.clustervars2));

data = data(:,1:Nmeasures*settings.clustervars*settings.clustervars2);
datatmp = reshape(data,[size(data,1) settings.clustervars settings.clustervars2 Nmeasures]);
datatmp = permute(datatmp,[3 4 1 2]);
if ~isempty(strfind(upper(header.subjects{1,1}),upper('median')))
    datatmp = median(datatmp,4);
else
    datatmp = median(datatmp,4);
end

xlsout = cell(size(datatmp,2),size(datatmp,1));
for i = 1:size(datatmp,1)
    for j = 1:size(datatmp,2)
        xlsout{j,i} = [num2str(datatmp(i,j,1),2) ' (' num2str(datatmp(i,j,2),2) '-' num2str(datatmp(i,j,3),2) ')'];
    end
end

if ~isfield(header,'variables2')
    header.variables2 = {'VAR'};
end

xlsout = cat(1,header.variables2(:)',xlsout);
xlsout = cat(2,cat(1,{''},header.measures(:)),xlsout);

lab_write_xls(settings.file,xlsout);