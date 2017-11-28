% Calculate Kappa
%    Cohen's Kappa for 2 raters
%    Fleiss'es Kappa for > 2 raters
%
% R = lab_kappa(data,settings)
%
% data = array (subjects x raters)
% settings (optional)
%
% Written by F. Hatz 2014

function R = lab_kappa(data,settings)

if ~exist('settings','var')
    settings = [];
end

if ~exist('data','var') | isempty(data)
    [data,settings] = lab_read_statistics([],-1,0,1,0,1);
    if isempty(data)
        return
    end
end



if size(data,2) > 2
    disp('Calculate Fleiss''es Kappa')
    classes = unique(data(:));
    datatmp = zeros(size(data,1),length(classes));
    for i = 1:size(data,1)
        for j = 1:size(data,2)
            tmp = find(classes == data(i,j));
            if ~isempty(tmp)
                datatmp(i,tmp) = datatmp(i,tmp) + 1;
            end
        end
    end

    R = fleiss(datatmp);
    
    if isfield(settings,'file')
        if isfield(settings,'vars')
            for i = 1:length(settings.vars)
                settings.file = [settings.file '_' settings.vars{i}];
            end
        end
        fid=fopen(fullfile(settings.path,[settings.file '_Kappa.txt']),'w');
        fprintf(fid,['kj:   ' num2str(R.kj) '\r\n']);
        fprintf(fid,['s.e.: ' num2str(R.sekj) '\r\n']);
        fprintf(fid,['z:    ' num2str(R.zkj) '\r\n']);
        fprintf(fid,['p:    ' num2str(R.pkj) '\r\n']);
        fprintf(fid,[repmat('-',1,60) '\r\n']);
        fprintf(fid,['Fleiss''es (overall) kappa = ' num2str(R.kappa,'%0.4f') '\r\n']);
        fprintf(fid,['kappa error = ' num2str(R.kappa_error,'%0.4f') '\r\n']);
        fprintf(fid,['kappa C.I. (' num2str((1-R.alpha)*100) ') = ' num2str(R.kappa_ci,'%0.4f\t') '\r\n']);
        if R.kappa<0
            fprintf(fid,'Poor agreement\r\n');
        elseif R.kappa>=0 && R.kappa<=0.2
            fprintf(fid,'Slight agreement\r\n');
        elseif R.kappa>0.2 && R.kappa<=0.4
            fprintf(fid,'Fair agreement\r\n');
        elseif R.kappa>0.4 && R.kappa<=0.6
            fprintf(fid,'Moderate agreement\r\n');
        elseif R.kappa>0.6 && R.kappa<=0.8
            fprintf(fid,'Substantial agreement\r\n');
        elseif R.kappa>0.8 && R.kappa<=1
            fprintf(fid,'Perfect agreement\r\n');
        end
        fprintf(fid,['z = ' num2str(R.z,'%0.4f') '\t p = ' num2str(R.p,'%0.4f') '\r\n']);
        if R.p<0.05
            fprintf(fid,'Reject null hypotesis: observed agreement is not accidental');
        else
            fprintf(fid,'Accept null hypotesis: observed agreement is accidental');
        end
        fclose(fid);
    end
elseif size(data,2) == 2
    disp('Calculate Cohen''s Kappa')
    classes = unique(data(:));
    datatmp = zeros(length(classes),length(classes));
    for i = 1:size(data,1)
        x = find(classes == data(i,1));
        y = find(classes == data(i,2));
        if ~isempty(x) & ~isempty(y)
            datatmp(x,y) = datatmp(x,y) + 1;
        end
    end

    R = kappa(datatmp);
    
    if isfield(settings,'file')
        if isfield(settings,'vars')
            for i = 1:length(settings.vars)
                settings.file = [settings.file '_' settings.vars{i}];
            end
        end
        fid=fopen(fullfile(settings.path,[settings.file '_Kappa.txt']),'w');
        fprintf(fid,['Observed agreement (po) = ' num2str(R.po,'%0.4f') '\r\n']);
        fprintf(fid,['Random agreement (pe) = ' num2str(R.pe,'%0.4f') '\r\n']);
        fprintf(fid,['Agreement due to true concordance (po-pe) = ' num2str(R.po-R.pe,'%0.4f') '\r\n']);
        fprintf(fid,['Residual not random agreement (1-pe) ' num2str(1-R.pe,'%0.4f') '\r\n']);
        fprintf(fid,[repmat('-',1,60) '\r\n']);
        fprintf(fid,['Cohen''s kappa = ' num2str(R.kappa,'%0.4f') '\r\n']);
        fprintf(fid,['kappa error = ' num2str(R.kappa_error,'%0.4f') '\r\n']);
        fprintf(fid,['kappa C.I. (' num2str((1-R.alpha)*100) ') = ' num2str(R.kappa_ci,'%0.4f\t') '\r\n']);
        if R.kappa<0
            fprintf(fid,'Poor agreement\r\n');
        elseif R.kappa>=0 && R.kappa<=0.2
            fprintf(fid,'Slight agreement\r\n');
        elseif R.kappa>0.2 && R.kappa<=0.4
            fprintf(fid,'Fair agreement\r\n');
        elseif R.kappa>0.4 && R.kappa<=0.6
            fprintf(fid,'Moderate agreement\r\n');
        elseif R.kappa>0.6 && R.kappa<=0.8
            fprintf(fid,'Substantial agreement\r\n');
        elseif R.kappa>0.8 && R.kappa<=1
            fprintf(fid,'Perfect agreement\r\n');
        end
        fprintf(fid,['z = ' num2str(R.z,'%0.4f') '\t p = ' num2str(R.p,'%0.4f') '\r\n']);
        if R.p<0.05
            fprintf(fid,'Reject null hypotesis: observed agreement is not accidental');
        else
            fprintf(fid,'Accept null hypotesis: observed agreement is accidental');
        end
        fclose(fid);
    end
end