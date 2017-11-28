function [] = txt2edf(filename, desired_filename)
% txt2edf(filename, desired_filename);
% Converts a .txt file into an .edf file.
% Because .txt files are not standardized, you will need to create a function
% that outputs a matrix of your data. Header information can be specified
% by replacing the header. variables below.
% Go to http://www.edfplus.info/specs/edf.html for EDF specifications
% As is, this will convert .txt files from MC_DataTool that were originally
% .mcd files.
% desired_filename is an optional argument that will rename the file.

% if a second argument is passed in, the new file will have that name.
% nargin is the number of input arguments
if nargin == 2
    savename = desired_filename;
    if isempty(strfind(savename,'.edf')) %if .edf isn't included in the desired_filename
        savename = [savename,'.edf']; %append it
    end
else  %if no desired_filename, the new file will have the original .txt filename
    savename = strrep(filename,'.txt','.edf'); %with a .edf extension
end

if exist(savename,'file')  % if the new file exists
    prompt = [savename,' already exists. Do you want to overwrite the file?']; % Create a prompt string
    response = questdlg(prompt,'File exists','Yes','No','No');  % Ask to overwrite
    if strcmp(response,'Yes') %strcmp compares the strings 'Yes' and response.
        try
            edffid = fopen(savename,'w+','ieee-le'); %create the file
        catch err
            error('Make sure to close all programs using the file you want to convert.');
        end
    end
else % if the new file does exist
    try
        edffid = fopen(savename,'w+','ieee-le'); %create the file
    catch err
        error('Make sure to close all programs using the file you want to convert.');
    end
end



if exist('edffid','var') % if we successfully opened the new file
    try
        txtfid = fopen(filename); % open the text file
    catch err
        error('The txt file could not be opened. Check your filename and working directory.');
    end
    
    % <!--
    % These lines are specific to .txt files converted from .mcd files (MC_DataTool).
    % Edit the lines using textscan to retrieve data from your .txt file
    % and return a matrix.
    
    textscan(txtfid,'%s',2,'delimiter','\n');  %Skips the first 2 lines of the document
    head = textscan(txtfid,'%s',2,'delimiter','\n');  %Retrieves unit information in 2nd two lines

    % EDIT ncol to match the number of columns in your data set.
    ncol = length(sscanf(head{1}{2},'%s'))/4;  %number of columns. Each unit column has four characters.
    format = '';
    for i=1:ncol
        format = [format, '%f '];   %Creates the format specifier for textscan
    end
    
    matrix = textscan(txtfid,format,600000,'CollectOutput',1); %grabs the first 600000 rows of the .txt document
    matrix = matrix{1};   %so we can determine 'hertz' without looping
    hertz = 1/(matrix(2,1)-matrix(1,1))*1000;
    matrix_length = 0;  
    
    %-->
    
    % This while loop is only necessary for large datasets.
    
    while ~isempty(matrix)
        % Preallocates a data vector of one row and (# of EEGs times # of samples)
        % columns
        data = ones(1,(ncol-1)*length(matrix(:,1)));
        k = 1;  % k is the count
        for i = 2:length(matrix(1,:))    %nested for loops travel down each column
            for j = 1:length(matrix(:,1))  %of the data matrix, top to bottom, left 
                data(1,k) = matrix(j,i); %to right. 'data' becomes one long array
                k = k+1;              %that will be converted into a binary string
            end
        end
        
        try
            fwrite(edffid,data,'int16','l');
        catch err
            error('Make sure to close all programs using the file you want to convert.');
        end
        
        matrix_length = matrix_length + length(matrix(:,1));        
        matrix = textscan(txtfid,format,600000,'CollectOutput',1); %grabs the next 600000 rows of the .txt document
        matrix = matrix{1};   %Data matrix
    end       

    colunits = textscan(head{1}{2},strrep(format,'%f ','%s '),'delimiter',' ');    %Uses 'format' to create a format specifier with strings to pick out units.
    for i=1:ncol
        colunits{i}(cellfun(@isempty,colunits{i})) = [];    %Deletes blank values in the colunits matrix
    end

    % Header information. If you can parse header information from your
    % .txt file, place the variable names here.
    % Everything must be a string. Put numbers in single quotes or use
    % num2str
    header.version = '0';
    header.patientinfo = 'Unknown';
    header.recordid = 'Unknown';
    header.startdate = '01.01.00';
    header.starttime = '01.00.00';
    header.bytes = num2str(256+((ncol-1)*256));
    header.reserved44 = ' ';
    header.numofrecords = '1';
    header.duration = num2str(matrix_length/hertz);
    header.numofsignals = num2str(ncol-1);
    header.transducertype = 'Pinnacle 4000';
    header.physdimension = 'uV';
    header.physmin = '-819';
    header.physmax = '819';
    header.digmin = '-1';
    header.digmax = '1';
    header.prefilter = 'HP:0.1Hz LP:75Hz';
    header.numofsamples = num2str(matrix_length);
    header.reserved32 = ' ';

    duration = header.duration;
    
    % Adds trailing spaces to header variables
    for i=1:8-length(header.version)
        header.version = [header.version,' '];
    end

    for i=1:80-length(header.patientinfo)
        header.patientinfo = [header.patientinfo,' '];
    end

    for i=1:80-length(header.recordid)
        header.recordid = [header.recordid,' '];
    end

    for i=1:8-length(header.bytes)
        header.bytes = [header.bytes,' '];
    end

    for i=1:43  % 43 because header.reserved44 already has a space.
        header.reserved44 = [header.reserved44,' '];
    end

    for i=1:8-length(header.numofrecords)
        header.numofrecords = [header.numofrecords,' '];
    end

    for i=1:8-length(header.duration)
        header.duration = [header.duration,' '];
    end

    for i=1:4-length(header.numofsignals)
        header.numofsignals = [header.numofsignals,' '];
    end

    header.label = ''; 
    for i=1:ncol-1
        label = ['EEG ',num2str(i)]; %Concatenates label with the value of i to number EEGs
        for j=1:16-length(label) % Adds trailing spaces to labels
            label = [label,' '];
        end
        header.label = [header.label,label]; %Concatenates labels into one string as they are constructed
    end

    % Trailing spaces
    for j=1:80-length(header.transducertype)
        header.transducertype = [header.transducertype,' '];
    end

    for i=1:8-length(header.physdimension)
        header.physdimension = [header.physdimension, ' '];
    end

    for j=1:8-length(header.physmin)
        header.physmin = [header.physmin,' '];
    end

    for j=1:8-length(header.physmax)
        header.physmax = [header.physmax,' '];
    end

    for j=1:8-length(header.digmin)
        header.digmin = [header.digmin,' '];
    end

    for j=1:8-length(header.digmax)
        header.digmax = [header.digmax,' '];
    end

    for j=1:80-length(header.prefilter)
        header.prefilter = [header.prefilter,' '];
    end

    for j=1:8-length(header.numofsamples)
        header.numofsamples = [header.numofsamples,' '];
    end

    for j=1:31
        header.reserved32 = [header.reserved32,' '];
    end

    % Placeholder variables, because the header. subvariables are about to
    % grow in length.
    h.tr = header.transducertype;
    h.pdm = header.physdimension;
    h.pmi = header.physmin;
    h.pma = header.physmax;
    h.dmi = header.digmin;
    h.dma = header.digmax;
    h.pf = header.prefilter;
    h.nos = header.numofsamples;
    h.r = header.reserved32;

    % In EDF, these parameters must be repeated for every column of data.
    for i=1:ncol-2
        header.transducertype = [header.transducertype, h.tr];
        header.physdimension = [header.physdimension, h.pdm];
        header.physmin = [header.physmin, h.pmi];
        header.physmax = [header.physmax, h.pma];
        header.digmin = [header.digmin, h.dmi];
        header.digmax = [header.digmax, h.dma];
        header.prefilter = [header.prefilter, h.pf];
        header.numofsamples = [header.numofsamples, h.nos];
        header.reserved32 = [header.reserved32, h.r];
    end

    % concatenates header strings
    edfheader = [header.version, header.patientinfo, header.recordid, header.startdate, header.starttime, header.bytes, header.reserved44, header.numofrecords, header.duration, header.numofsignals, header.label, header.transducertype, header.physdimension, header.physmin, header.physmax, header.digmin, header.digmax, header.prefilter, header.numofsamples, header.reserved32];

    frewind(edffid) % moves the pointer to the beginning of the edf file.
    try
        fprintf(edffid,edfheader);  
    catch err
        error('Make sure the file you want to convert is closed in all other programs.');
    end

    fclose('all');

end