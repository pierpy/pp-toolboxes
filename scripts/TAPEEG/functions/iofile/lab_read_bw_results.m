function Result = lab_read_bw_results(Filename)

if ~exist('Filename','var') | isempty(Filename)
   [Filename,Filepath]=uigetfile('*.*','Select file');
   Result_file = fullfile(Filepath,Filename);
else
   Result_file = Filename;
end

Result.results = {};
Result.resultsmean = {};
Filenames = {};

fid=fopen(Result_file);
tline = fgetl(fid);
doepoch = false;
if ~strcmp(tline(1:9),'Filename:')
    Result = [];
    disp('   No valid BrainWave Result-File');
    return
else
    tmp = textscan(regexprep(tline,'Ch. ','Ch.'),'%s');
    numchannels = size(tmp{1},1) - 4;
    tline = fgetl(fid);
end
while ischar(tline)
    if ~isempty(tline)
        line = textscan(tline,'%s');
    else
        line = [];
    end
    if ~isempty(line) & size(line{1},1) >= 3
        line = line{1};
        name = line{1,1};
        epoch = str2num(line{2,1});
        if isempty(epoch)
            epoch = '';
            while ismember(name(end),48:57) & length(name) > 1
                epoch = [name(end) epoch];
                name = name(1:end-1);
            end
            if ~isempty(epoch)
                epoch = str2num(epoch);
            else
                epoch = 1;
            end
            line = line(2:end,1);

        else
            line = line(3:end,1);
        end
        Nfilename = find(strcmp(Filenames,[name '_' num2str(epoch)]),1);
        if isempty(Nfilename)
            Filenames{1,end+1} = [name '_' num2str(epoch)];
            Nfilename = size(Filenames,2);
            subjects1{1,Nfilename} = name;
            subjects2{1,Nfilename} = [name '_' num2str(epoch)];
            if epoch > 1
                doepoch = true;
            end
        end
        Nresults = find(strcmp(Result.resultsmean,line{1,1}),1);
        if isempty(Nresults)
            Result.resultsmean{end+1,1} = line{1};
            Nresults = size(Result.resultsmean,1);
        end
        Result.datamean(Nresults,Nfilename) = str2double(strrep(line{2,1},',','.'));
        clearvars Nresults
        if size(line,1) > 2
            line2 = line(3:end,1);
            if size(line2,1)  == numchannels
                tmp = [];
                for i = 1:size(line2,1)
                    tmp = [tmp;str2double(strrep(line2{i,1},',','.'))];
                end
            else
                disp('    Reading error, try fixed length of 8')
                Result.readingerror = true;
                tmp = [];
                for i = 1:size(line2,1)
                    tmp2 = line2{i,1};
                    while ~isempty(tmp2)
                        if length(tmp2) >= 8
                            tmp = [tmp;str2double(strrep(tmp2(1,1:8),',','.'))];
                            if length(tmp2) > 8
                                tmp2 = tmp2(1,9:end);
                            else
                                tmp2 = '';
                            end
                        else
                            tmp = [tmp;str2double(strrep(tmp2,',','.'))];
                            tmp2 = '';
                        end
                    end
                end
                if length(tmp) ~= numchannels
                    disp('    Reading not possible, replacing bad values by NaN')
                    tmp = NaN(numchannels,1);
                end   
            end
            Nresults = find(strcmp(Result.results,line{1,1}),1);
            if isempty(Nresults)
                Result.results{end+1,1} = line{1};
                Nresults = size(Result.results,1);
            end
            Result.data(:,Nfilename,Nresults) = tmp;
            clearvars tmp tmp2 Nresults
        end
        clearvars Nfilename
    end
    tline = fgetl(fid);
end
fclose(fid);

if doepoch == true & exist('subjects2','var')
    Result.subjects = subjects2;
elseif exist('subjects1','var')
    Result.subjects = subjects1;
end