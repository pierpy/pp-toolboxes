function table = lab_montage2table(montage)

if ~isempty(montage)
    if size(montage,2) > 1
        montage = montage(1,1);
    end
    for i = 1:size(montage.chans,1)
        if montage.chans{i,2} == 0
            tmp{i,1} = num2str(montage.chans{i,1});
        else
            tmp{i,1} = ['A' num2str(montage.chans{i,1})];
        end
        if montage.chans{i,4} == 0
            tmp{i,2} = num2str(montage.chans{i,3});
        else
            tmp{i,2} = ['A' num2str(montage.chans{i,3})];
        end
    end
    table = [montage.label tmp(:,1) tmp(:,2)];
    clearvars tmp
else
    table = [];
end