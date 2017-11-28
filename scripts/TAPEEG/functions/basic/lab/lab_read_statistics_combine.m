% Helper function for lab_read_statistics
%
% written by F. Hatz 2012

function [dataout,groupout,cfg,postfix,skipprocessing] = lab_read_statistics_combine(datainput,datainput2,group,cfg,postfix)
   skipprocessing = 0;
   dataout = datainput;
   groupout = group;
   
   if isempty(datainput2)
       disp('No match of data, exclude last file and proceed')
       skipprocessing = 1;
       return
   end
   vars1 = datainput(2:end,1);
   vars1 = correct_double(vars1);
   vars2 = datainput2(2:end,1);
   vars2 = correct_double(vars2);
   subjects1 = datainput(1,2:end);
   subjects1 = correct_double(subjects1);
   subjects2 = datainput2(1,2:end);
   subjects2 = correct_double(subjects2);
   datainput = datainput(2:end,2:end);
   datainput2 = datainput2(2:end,2:end);
   Multiple.Vars = {};
   Multiple.Subjects = {};
   if isfield(postfix,'enable') & isfield(postfix,'mode')
       if strcmp(postfix.mode,'subjects (post)') | strcmp(postfix.mode,'subjects (pre)')
           for i = 1:size(subjects1,2)
               tmp = strfind(subjects1{1,i},'_');
               if ~isempty(tmp)
                   if strcmp(postfix.mode,'subjects (post)')
                       subjectstmp{1,i} = subjects1{1,i}(1:tmp(end)-1); %#ok<AGROW>
                       if length(subjects1{1,i}) > tmp(end)
                           Multiple.Subjects{1,i} = subjects1{1,i}(tmp(end)+1:end);
                       end
                   else
                       subjectstmp{1,i} = subjects1{1,i}(tmp(1)+1:end); %#ok<AGROW>
                       if tmp(1) > 1
                           Multiple.Subjects{1,i} = subjects1{1,i}(1:tmp(1)-1);
                       end
                   end
               end
           end
           if ~isempty(Multiple.Subjects)
               if verLessThan('matlab','8')
                   Multiple.Subjects = unique(Multiple.Subjects);
                   subjectstmp = unique(subjectstmp);
               else
                   Multiple.Subjects = unique(Multiple.Subjects,'stable');
                   subjectstmp = unique(subjectstmp,'stable');
               end
               Multiple.Subjects = Multiple.Subjects(:);
               Multiple.Nx = length(Multiple.Subjects);
               Multiple.Ny = 1;
               if mod(size(datainput,2),Multiple.Nx) == 0
                   subjects1 = subjectstmp;
                   if strcmp(postfix.mode,'subjects (post)')
                       datainput = reshape(datainput,[size(datainput,1) Multiple.Nx size(datainput,2)/Multiple.Nx]);
                       datainput = permute(datainput,[1 3 2]);
                   else
                       datainput = reshape(datainput,[size(datainput,1) size(datainput,2)/Multiple.Nx Multiple.Nx]);
                   end
               else
                   Multiple.Subjects = {};
               end
           end
           clearvars i subjectstmp
       elseif strcmp(postfix.mode,'measures (post)') | strcmp(postfix.mode,'measures (pre)')
           for i = 1:size(vars1,1)
               tmp = strfind(vars1{i,1},'_');
               if ~isempty(tmp)
                   if strcmp(postfix.mode,'measures (post)')
                       varstmp{i,1} = vars1{i,1}(1:tmp(end)-1); %#ok<AGROW>
                       if length(vars1{i,1}) > tmp(end)
                           Multiple.Vars{1,i} = vars1{i,1}(tmp(end)+1:end);
                       end
                   else
                       varstmp{i,1} = vars1{i,1}(tmp(1)+1:end); %#ok<AGROW>
                       if tmp(1) > 1
                           Multiple.Vars{1,i} = vars1{i,1}(1:tmp(1)-1);
                       end
                   end
               end
           end
           if ~isempty(Multiple.Vars)
               if verLessThan('matlab','8')
                   Multiple.Vars = unique(Multiple.Vars);
                   varstmp = unique(varstmp);
               else
                   Multiple.Vars = unique(Multiple.Vars,'stable');
                   varstmp = unique(varstmp,'stable');
               end
               Multiple.Vars = Multiple.Vars(:);
               Multiple.Ny = length(Multiple.Vars);
               Multiple.Nx = 1;
               if mod(size(datainput,1),Multiple.Ny) == 0
                   vars1 = varstmp;
                   if strcmp(postfix.mode,'measures (post)')
                       datainput = reshape(datainput,[Multiple.Ny size(datainput,1)/Multiple.Ny size(datainput,2)]);
                       datainput = permute(datainput,[2 3 1]);
                   else
                       datainput = reshape(datainput,[size(datainput,1)/Multiple.Ny Multiple.Ny size(datainput,2)]);
                       datainput = permute(datainput,[1 3 2]);
                   end
               else
                   Multiple.Vars = {};
               end
           end
           clearvars i varstmp
       end
   end
   if ~isfield(postfix,'correct_subjects')
       subjects1 = regexprep(subjects1,{'*','/'},'');
       vars1 = regexprep(vars1,{'*','/'},'');
       cfg2 = [];
       for i = 1:length(subjects1)
           [subjects1{i},cfg2,skipprocessing] = lab_subjectname(subjects1{i},cfg2);
           if skipprocessing == 1
               break;
           end
       end
       postfix.correct_subjects = 1;
   end
   subjects2 = regexprep(subjects2,{'*','/'},'');
   vars2 = regexprep(vars2,{'*','/'},'');
   cfg2 = [];
   for i = 1:length(subjects2)
       [subjects2{i},cfg2,skipprocessing] = lab_subjectname(subjects2{i},cfg2);
       if skipprocessing == 1
           break;
       end
   end
   
   if verLessThan('matlab','8')
       vars = union(vars1,vars2);
       [~,VindexN1,Vindex1] = intersect(vars,vars1);
       [~,VindexN2,Vindex2] = intersect(vars,vars2);
       subjects = union(subjects1,subjects2);
       [~,SindexN1,Sindex1] = intersect(subjects,subjects1);
       [~,SindexN2,Sindex2] = intersect(subjects,subjects2);
   else
       vars = union(vars1,vars2,'stable');
       [~,VindexN1,Vindex1] = intersect(vars,vars1,'stable');
       [~,VindexN2,Vindex2] = intersect(vars,vars2,'stable');
       subjects = union(subjects1,subjects2,'stable');
       [~,SindexN1,Sindex1] = intersect(subjects,subjects1,'stable');
       [~,SindexN2,Sindex2] = intersect(subjects,subjects2,'stable');
   end
   if size(vars,1) == 1
       vars = vars';
   end
   if size(subjects,2) == 1
       subjects = subjects';
   end
   
   if length(subjects1) + length(subjects2) == length(subjects) | ...
           length(vars1) + length(vars2) == length(vars)
       grouptmp = [cellstr('group') repmat(num2cell(group{end}+1),1,size(datainput2,2))];
       dataout = num2cell(NaN(length(vars),length(subjects)));
       dataout(VindexN1,SindexN1) = datainput(Vindex1,Sindex1);
       dataout(VindexN2,SindexN2) = datainput2(Vindex2,Sindex2);
       dataout = cat(1,subjects,dataout);
       dataout = cat(2,[{''};vars(:)],dataout);
       groupout(1,SindexN2+1) = grouptmp(2:end);
       groupout(1,SindexN1+1) = group(2:end);
       clearvars datainput2 grouptmp Vindex1 Vindex2 Sindex
   else
       if ~isfield(postfix,'enable')
           Prompt = cell(0,2);
           Formats = [];
           Prompt(end+1,:) = {'Tag File 1', 'file1'};
           Formats(end+1,1).type = 'edit';
           Formats(end,1).format = 'text';
           Formats(end,1).size = 250;
           if ~isfield(postfix,'mode') | isempty(postfix.mode)
               postfix.mode = 'subjects (post)';
               Prompt(end+1,:) = {'for','mode'};
               Formats(end+1,1).type = 'list';
               Formats(end,1).style = 'popupmenu';
               Formats(end,1).format = 'input';
               Formats(end,1).items = {'measures (post)','measures (pre)','subjects (post)','subjects (pre)'};
           end
           [postfix,Cancelled] = inputsdlg(Prompt,'Tag File 1',Formats,postfix,2);
           if Cancelled == 1
               skipprocessing = 1;
               return
           end
           if strcmp(postfix.mode,'measures (post)')
               Multiple.Vars{end+1,1} = postfix.file1;
               Multiple.Ny = 1;
               Multiple.Nx = 1;
               if isfield(cfg,'clustervars') & cfg.clustervars > 1
                   cfg.clustervars2 = cfg.clustervars;
               end
           elseif strcmp(postfix.mode,'measures (pre)')
               Multiple.Vars{end+1,1} = postfix.file1;
               Multiple.Ny = 1;
               Multiple.Nx = 1;
               if ~isfield(cfg,'clustervars2') | cfg.clustervars2 == 1
                   cfg.clustervars2 = cfg.clustervars;
               end
           else
               Multiple.Subjects{end+1,1} = postfix.file1;
               Multiple.Ny = 1;
               Multiple.Nx = 1;
           end
           postfix.enable = 1;
       end
       Prompt = cell(0,2);
       Formats = [];
       Prompt(end+1,:) = {['Tag File ' num2str(size(datainput,3)+1)], 'file2'};
       Formats(end+1,1).type = 'edit';
       Formats(end,1).format = 'text';
       Formats(end,1).size = 250;
       [postfix,Cancelled] = inputsdlg(Prompt,['Tag File ' num2str(size(datainput,3)+1)],Formats,postfix);
       if Cancelled == 1
           skipprocessing = 1;
           return
       else
           if strcmp(postfix.mode,'measures (post)') | strcmp(postfix.mode,'measures (pre)')
               Multiple.Vars{end+1,1} = postfix.file2;
               Multiple.Ny = Multiple.Ny + 1;
           else
               Multiple.Subjects{end+1,1} = postfix.file2;
               Multiple.Nx = Multiple.Nx + 1;
           end
       end
       dataout1 = num2cell(NaN(length(vars),length(subjects),size(datainput,3)));
       dataout1(VindexN1,SindexN1,:) = datainput(Vindex1,Sindex1,:);
       dataout2 = num2cell(NaN(length(vars),length(subjects)));
       dataout2(VindexN2,SindexN2) = datainput2(Vindex2,Sindex2);
       dataout = cat(3,dataout1,dataout2);
       
       if strcmp(postfix.mode,'measures (post)')
           varsout = {};
           for i = 1:length(vars)
               for j = 1:Multiple.Ny
                   varsout{end+1,1} = [vars{i} '_' Multiple.Vars{j}]; %#ok<AGROW>
               end
           end
           dataout = permute(dataout,[3 1 2]);
           dataout = reshape(dataout,[size(dataout,1)*size(dataout,2) size(dataout,3)]);
           subjectsout = subjects;
           cfg.clustervars = Multiple.Ny;
       elseif strcmp(postfix.mode,'measures (pre)')
           varsout = {};
           for j = 1:Multiple.Ny
               for i = 1:length(vars)
                   varsout{end+1,1} = [Multiple.Vars{j} '_' vars{i}]; %#ok<AGROW>
               end
           end
           dataout = permute(dataout,[1 3 2]);
           dataout = reshape(dataout,[size(dataout,1)*size(dataout,2) size(dataout,3)]);
           subjectsout = subjects;
           cfg.clustervars = Multiple.Ny;
       elseif strcmp(postfix.mode,'subjects (pre)')
           subjectsout = {};
           for i = 1:length(subjects)
               for j = 1:Multiple.Nx
                   subjectsout{1,end+1} = [Multiple.Subjects{j} '_' subjects{i}]; %#ok<AGROW>
               end
           end
           dataout = permute(dataout,[1 3 2]);
           dataout = reshape(dataout,[size(dataout,1) size(dataout,2)*size(dataout,3)]);
           varsout = vars;
       else
           subjectsout = {};
           for j = 1:Multiple.Nx
               for i = 1:length(subjects)
                   subjectsout{1,end+1} = [subjects{i} '_' Multiple.Subjects{j}]; %#ok<AGROW>
               end
           end
           dataout = reshape(dataout,[size(dataout,1) size(dataout,2)*size(dataout,3)]);
           varsout = vars;
       end
       dataout = cat(1,subjectsout,dataout);
       dataout = cat(2,[{''};varsout],dataout);
       groupout = [{'group'} repmat({1},1,size(dataout,2)-1)];
   end
end

function list = correct_double(list)
   tmp = {};
   for i = 1:length(list)
       tmp2 = find(strcmp(tmp,list{i}),1);
       if ~isempty(list) & ~isempty(tmp2)
           list{i} = [list{i} num2str(length(tmp2)+1)];
       end
       tmp{end+1,1} = list{i}; %#ok<AGROW>
   end
end