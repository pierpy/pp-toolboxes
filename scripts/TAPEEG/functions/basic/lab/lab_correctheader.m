function datainput = lab_correctheader(datainput)
   if isnumeric(datainput)
       datainput = num2cell(datainput);
   end
   if size(datainput,1) == 1 | size(datainput,2) == 1
       if isnumeric(datainput{1,1})
           for i = 1:size(datainput,2)
               header{1,i} = ['Trial' num2str(i)]; %#ok<AGROW>
           end
           for i = 1:size(datainput,1)
               measure{i,1} = ['Measure' num2str(i)]; %#ok<AGROW>
           end
           measure = [{''};measure];
           datainput = [measure cat(1,header,datainput)];
       end
   elseif isnumeric(datainput{1,1}) & isnumeric(datainput{2,1}) & isnumeric(datainput{1,2})
       for i = 1:size(datainput,2)
           header{1,i} = ['Trial' num2str(i)]; %#ok<AGROW>
       end
       for i = 1:size(datainput,1)
           measure{i,1} = ['Measure' num2str(i)]; %#ok<AGROW>
       end
       measure = [{''};measure];
       datainput = [measure cat(1,header,datainput)];
   elseif isnumeric(datainput{2,1}) & ~isnumeric(datainput{1,2})
       for i = 1:size(datainput,1)-1
           measure{i,1} = ['Measure' num2str(i)]; %#ok<AGROW>
       end
       measure = [{''};measure];
       datainput = [measure datainput];
   elseif isnumeric(datainput{1,2}) & ~isnumeric(datainput{2,1})
       for i = 1:size(datainput,2)-1
           header{1,i} = ['Trial' num2str(i)]; %#ok<AGROW>
       end
       header = [{''} header];
       datainput = cat(1,header,datainput);
   end
   tmp = find(cellfun(@isnumeric,datainput(1,2:end)));
   for i = tmp
       datainput{1,i+1} = num2str(datainput{1,i+1});
   end
   clearvars tmp
   tmp = find(cellfun(@isnumeric,datainput(2:end,1)));
   if isempty(tmp)
       return
   end
   for i = tmp'
       if isempty(datainput{i+1,1})
           datainput{i+1,1} = '';
       else
           datainput{i+1,1} = num2str(datainput{i+1,1});
       end
   end
   clearvars tmp
end