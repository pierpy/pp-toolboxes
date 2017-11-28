function [datainput,cfg] = lab_getstructure(datainput,cfg)
   if ~exist('cfg','var') | ~isfield(cfg,'clustervars')
       cfg.clustervars = [];
   end
   if ~isfield(cfg,'numresults')
       cfg.numresults = 0;
   end
   if ~isempty(datainput)
       tmp = datainput{1,1};
       if ~isempty(tmp) & ischar(tmp)
           tmp = textscan(tmp,'%s');
           tmp = tmp{1,1};
           if strcmp(tmp{1,1}(1),'C')
               cfg.clustervars = str2num(tmp{1,1}(2:end)); %#ok<ST2NM>
           end
           if size(tmp,1) > 1 & strcmp(tmp{2,1}(1),'R')
               cfg.numresults = str2num(tmp{2,1}(2:end)); %#ok<ST2NM>
           end
           if size(tmp,1) > 2 & strcmp(tmp{3,1}(1),'V')
               cfg.clustervars2 = str2num(tmp{3,1}(2:end)); %#ok<ST2NM>
           end
       end
       clearvars tmp
   end
   flag = [];
   for i = 2:size(datainput,1)
       tmp = strfind(datainput{i,1},'_');
       if ~isempty(tmp)
           flag(i-1) = true; %#ok<AGROW>
       else
           flag(i-1) = false; %#ok<AGROW>
       end
   end
   if ~isempty(flag) & min(flag) == false
       flag = [];
       for i = 2:size(datainput,1)
           tmp = strfind(datainput{i,1},' ');
           if ~isempty(tmp)
               flag(i-1) = true; %#ok<AGROW>
           else
               flag(i-1) = false; %#ok<AGROW>
           end
       end
       if ~isempty(flag) & min(flag) == true
           for i = 2:size(datainput,1)
               datainput{i,1} = regexprep(datainput{i,1},' ','_');
           end
       end
   end
   if isempty(cfg.clustervars)
       for i = 2:size(datainput,1)
           tmp = strfind(datainput{i,1},'_');
           if ~isempty(tmp)
               Vars{i-1,1} = datainput{i,1}(1:tmp(end)-1); %#ok<AGROW>
           else
               Vars{i-1,1} = datainput{i,1}; %#ok<AGROW>
           end
       end
       [~,m,~]=unique(Vars,'last');
       m = sort(m);
       n = mod(m,m(1));
       if length(m) == 1 | (m(1) > 1 & (m(2)-m(1)) == m(1))
           cfg.clustervars = m(1);
           cfg.numresults = size(Vars,1) - m(find(n>0,1)) + 1;
           if isempty(cfg.numresults)
               cfg.numresults = 0;
           end
       else
           cfg.clustervars = 1;
       end
       clearvars tmp m n i
   end
   if ~isempty(cfg.clustervars) & cfg.clustervars > 1 & (~isfield(cfg,'clustervars2') | isempty(cfg.clustervars2))
       Vars = datainput(2:cfg.clustervars:end,1);
       for i = 1:size(Vars,1)
           tmp = strfind(Vars{i,1},'_');
           if length(tmp) > 1
               Vars2{i,1} = Vars{i,1}(1:tmp(end-1)-1); %#ok<AGROW>
           else
               Vars2{i,1} = Vars{i,1}; %#ok<AGROW>
           end
       end
       [~,m,~]=unique(Vars2,'last');
       m = sort(m);
       if length(m) == 1 | (m(1) > 1 & (m(2)-m(1)) == m(1))
           cfg.clustervars2 = m(1);
       else
           cfg.clustervars2 = 1;
       end
       clearvars tmp m n i
   end
   clearvars Vars Vars2
   if isfield(cfg,'clustervars2') & ~isempty(cfg.clustervars2) & cfg.clustervars2 > 1
       datainput{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults) ' V' num2str(cfg.clustervars2)];
   else
       datainput{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults)];
   end
   if cfg.clustervars == 1
       for i = 2:size(datainput,1)
           datainput{i,1} = regexprep(datainput{i,1},'_','');
       end
   end
end