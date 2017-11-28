function Labels= lab_correct_labels(Labels,Labels2)

NamesStd = {'T3','T4','T5','T6'};
NamesMCN = {'T7','T8','P7','P8'};
FlagStd = 0;
FlagMCN = 0;
for i = 1:length(Labels2)
    for j = 1:4
        if ~isempty(strfind(upper(Labels2{i}),upper(NamesStd{j})))
            FlagStd = 1;
        end
        if ~isempty(strfind(upper(Labels2{i}),upper(NamesMCN{j})))
            FlagMCN = 1;
        end
    end
end
if FlagStd == 1
    for i = 1:length(Labels)
        for j = 1:4
            if ~isempty(strfind(upper(Labels{i}),upper(NamesMCN{j})))
                Labels{i} = regexprep(Labels{i},NamesMCN{j},NamesStd{j});
            end
        end
    end
end
if FlagMCN == 1
    for i = 1:length(Labels)
        for j = 1:4
            if ~isempty(strfind(upper(Labels{i}),upper(NamesStd{j})))
                Labels{i} = regexprep(Labels{i},NamesStd{j},NamesMCN{j});
            end
        end
    end
end

end