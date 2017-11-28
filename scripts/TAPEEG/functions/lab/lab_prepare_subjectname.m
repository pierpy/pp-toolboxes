function [Output,Foldername,Filename,Mainfolder] = lab_prepare_subjectname(Filename)

[~,Filepath,~,Filename] = lab_filename(Filename);
Filename = regexprep(Filename,'Conn_F','F');   
Foldername = '';
Mainfolder = '';
if ~isempty(Filepath)
    tmp = strfind(Filepath,filesep);
    if ~isempty(tmp) & length(tmp) > 1
        Mainfolder = Filepath(1:tmp(end-1));
        for i = 1:length(tmp)-1
            Name = Filepath(tmp(i)+1:tmp(i+1)-1);
            Name = regexprep(Name,{' ','_'},'');
            if isempty(Foldername)
                Foldername = Name;
            else
                Foldername = [Foldername '_' Name]; %#ok<AGROW>
            end
        end
    end
end

if ~isempty(Foldername)
    Output{1} = ['( ' create_name(Foldername,1) ' )'];
else
    Output{1} = [];
end
Output{2} = ['( ' create_name(Filename,0) ' )'];

end

function NameOut = create_name(NameIn,Inverse)

tmp = union(strfind(NameIn,'_'),strfind(NameIn,' '));
if ~isempty(tmp)
    if tmp(1) > 1
        NameOut = NameIn(1:tmp(1)-1);
    else
        NameOut = '';
    end
    for i = 1:length(tmp)
        if Inverse == 0
            Nflag = num2str(i-1);
        else
            Nflag = ['-' num2str(length(tmp)+2-i)];
        end
        numi = '';
        for j = 1:length(Nflag)
            numi = [numi '^' Nflag(j)];
        end
        if i ~= length(tmp)
            NameOut = [NameOut numi ' ' NameIn(tmp(i)+1:tmp(i+1)-1)];
        elseif tmp(i) == length(NameIn)
            NameOut = [NameOut numi];
        else
            NameOut = [NameOut numi ' ' NameIn(tmp(i)+1:end)];
        end
    end
    if Inverse == 0
        Nflag = num2str(i);
        numi = '';
        for j = 1:length(Nflag)
            numi = [numi '^' Nflag(j)];
        end
    else
        numi = '^-^1';
    end
    
    NameOut = [NameOut numi];
else
    NameOut = [NameIn '^0'];
end
clearvars tmp Nflag numi

end