function bactivity = lab_get_bactivity(numchans)

global Main_Path BactivityN Bactivity

if exist(fullfile(Main_Path,'BackgroundActivity.xls'),'file')
    if ispc
        [BactivityN,~,Bactivity] = xlsread(fullfile(Main_Path,'BackgroundActivity.xls'));
    else
        [BactivityN,~,Bactivity] = xlsread(fullfile(Main_Path,'BackgroundActivity.xls'),1,'','basic');
    end
end

if ~isempty(find(BactivityN==numchans,1)) & size(Bactivity,2) == 2
    bactivity = Bactivity{find(BactivityN==numchans,1),2};
    if ischar(bactivity)
        bactivity = str2num(bactivity); %#ok<ST2NM>
    end
else
    bactivity = [];
end