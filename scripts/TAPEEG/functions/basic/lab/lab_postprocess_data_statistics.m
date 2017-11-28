% Helper function for lab_read_statistics
%
% written by F. Hatz 2012

function Rstat = lab_postprocess_data_statistics(Rstat,header,cfg)

if isfield(cfg,'includevars')
    allvars = union(cfg.includevars,cfg.excludevars);
    if size(allvars,1) > 1
        allvars = allvars';
    end
    clustervars2 = size(allvars,2);
    for Nresult = 1:cfg.resultvars
        Rstat2{Nresult}.p = ones(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        if isfield(Rstat{Nresult},'mxp')
            Rstat2{Nresult}.mxp = ones(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        end
        if isfield(Rstat{Nresult},'mxpC')
            Rstat2{Nresult}.mxpC = ones(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        end
        if isfield(Rstat{Nresult},'mxpV')
            Rstat2{Nresult}.mxpV = ones(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        end
        if isfield(Rstat{Nresult},'mxpS')
            Rstat2{Nresult}.mxpS = ones(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        end
        if isfield(Rstat{Nresult},'T')
            Rstat2{Nresult}.T = zeros(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        elseif isfield(Rstat{Nresult},'R')
            Rstat2{Nresult}.R = zeros(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        elseif isfield(Rstat{Nresult},'F')
            Rstat2{Nresult}.F = zeros(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        elseif isfield(Rstat{Nresult},'Z')
            Rstat2{Nresult}.Z = zeros(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        elseif isfield(Rstat{Nresult},'df')
            Rstat2{Nresult}.df = zeros(1,cfg.numclusters*clustervars2); %#ok<AGROW>
        end
        if isfield(Rstat{Nresult},'mean')
            for j = 1:length(Rstat{Nresult}.mean)
                Rstat2{Nresult}.mean{j} = zeros(1,cfg.numclusters*clustervars2); %#ok<AGROW>
                Rstat2{Nresult}.std{j} = zeros(1,cfg.numclusters*clustervars2); %#ok<AGROW>
            end
        end
        for i = 1:cfg.numclusters
            calcstart = (i-1)*cfg.clustervars + 1;
            calcend = i * cfg.clustervars;
            calc2 = cfg.includevars + (i-1)*clustervars2;
            Rstat2{Nresult}.p(1,calc2) = Rstat{Nresult}.p(1,calcstart:calcend); %#ok<AGROW>
            if isfield(Rstat{Nresult},'mxp')
                Rstat2{Nresult}.mxp(1,calc2) = Rstat{Nresult}.mxp(1,calcstart:calcend); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpC')
                Rstat2{Nresult}.mxpC(1,calc2) = Rstat{Nresult}.mxpC(1,calcstart:calcend); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpV')
                Rstat2{Nresult}.mxpV(1,calc2) = Rstat{Nresult}.mxpV(1,calcstart:calcend); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpS')
                Rstat2{Nresult}.mxpS(1,calc2) = Rstat{Nresult}.mxpS(1,calcstart:calcend); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'T')
                Rstat2{Nresult}.T(1,calc2) = Rstat{Nresult}.T(1,calcstart:calcend); %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'R')
                Rstat2{Nresult}.R(1,calc2) = Rstat{Nresult}.R(1,calcstart:calcend); %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'F')
                Rstat2{Nresult}.F(1,calc2) = Rstat{Nresult}.F(1,calcstart:calcend); %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'Z')
                Rstat2{Nresult}.Z(1,calc2) = Rstat{Nresult}.Z(1,calcstart:calcend); %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'df')
                Rstat2{Nresult}.df(1,calc2) = Rstat{Nresult}.df(1,calcstart:calcend); %#ok<AGROW>
            end
            if isfield(Rstat2{Nresult},'mean')
                for j = 1:length(Rstat2{Nresult}.mean)
                    Rstat2{Nresult}.mean{j}(1,calc2) = Rstat{Nresult}.mean{j}(1,calcstart:calcend); %#ok<AGROW>
                    Rstat2{Nresult}.std{j}(1,calc2) = Rstat{Nresult}.std{j}(1,calcstart:calcend); %#ok<AGROW>
                end
            end
        end
    end
    Rstat = Rstat2;
    cfg.clustervarsI = cfg.clustervars;
    cfg.clustervars = clustervars2;
    clearvars Rstat2 clustervars2
    header.vars = header.varsall;
    header = rmfield(header,'varsall');
end

if isfield(cfg,'Latindex') & cfg.Latindex == 1
    if strcmp(Latmode,'LL..RR..')
        for Nresult = 1:cfg.resultvars
            Rstat2{Nresult}.p = []; %#ok<AGROW>
            if isfield(Rstat{Nresult},'mxp')
                Rstat2{Nresult}.mxp = []; %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpC')
                Rstat2{Nresult}.mxpC = []; %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpV')
                Rstat2{Nresult}.mxpV = []; %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpS')
                Rstat2{Nresult}.mxpS = []; %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'T')
                Rstat2{Nresult}.T = []; %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'R')
                Rstat2{Nresult}.R = []; %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'F')
                Rstat2{Nresult}.F = []; %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'Z')
                Rstat2{Nresult}.Z = []; %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'df')
                Rstat2{Nresult}.df = []; %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mean')
                for j = 1:length(Rstat{Nresult}.mean)
                    Rstat2{Nresult}.mean{j} = []; %#ok<AGROW>
                    Rstat2{Nresult}.std{j} = []; %#ok<AGROW>
                end
            end
            for i = 1:cfg.numclusters
                calcstart = (i-1)*cfg.clustervars + 1;
                calcend = i * cfg.clustervars;
                Rstat2{Nresult}.p = cat(2,Rstat2{Nresult}.p,Rstat{Nresult}.p(1,calcstart:calcend),Rstat{Nresult}.p(1,calcstart:calcend)); %#ok<AGROW>
                if isfield(Rstat{Nresult},'mxp')
                    Rstat2{Nresult}.mxp = cat(2,Rstat2{Nresult}.mxp,Rstat{Nresult}.mxp(1,calcstart:calcend),Rstat{Nresult}.mxp(1,calcstart:calcend)); %#ok<AGROW>
                end
                if isfield(Rstat{Nresult},'mxpC')
                    Rstat2{Nresult}.mxpC = cat(2,Rstat2{Nresult}.mxpC,Rstat{Nresult}.mxpC(1,calcstart:calcend),Rstat{Nresult}.mxpC(1,calcstart:calcend)); %#ok<AGROW>
                end
                if isfield(Rstat{Nresult},'mxpV')
                    Rstat2{Nresult}.mxpV = cat(2,Rstat2{Nresult}.mxpV,Rstat{Nresult}.mxpV(1,calcstart:calcend),Rstat{Nresult}.mxpV(1,calcstart:calcend)); %#ok<AGROW>
                end
                if isfield(Rstat{Nresult},'mxpS')
                    Rstat2{Nresult}.mxpS = cat(2,Rstat2{Nresult}.mxpS,Rstat{Nresult}.mxpS(1,calcstart:calcend),Rstat{Nresult}.mxpS(1,calcstart:calcend)); %#ok<AGROW>
                end
                if isfield(Rstat{Nresult},'T')
                    Rstat2{Nresult}.T = cat(2,Rstat2{Nresult}.T,Rstat{Nresult}.T(1,calcstart:calcend),-Rstat{Nresult}.T(1,calcstart:calcend)); %#ok<AGROW>
                elseif isfield(Rstat{Nresult},'R')
                    Rstat2{Nresult}.R = cat(2,Rstat2{Nresult}.R,Rstat{Nresult}.R(1,calcstart:calcend),-Rstat{Nresult}.R(1,calcstart:calcend)); %#ok<AGROW>
                elseif isfield(Rstat{Nresult},'F')
                    Rstat2{Nresult}.F = cat(2,Rstat2{Nresult}.F,Rstat{Nresult}.F(1,calcstart:calcend),-Rstat{Nresult}.F(1,calcstart:calcend)); %#ok<AGROW>
                elseif isfield(Rstat{Nresult},'Z')
                    Rstat2{Nresult}.Z = cat(2,Rstat2{Nresult}.Z,Rstat{Nresult}.Z(1,calcstart:calcend),-Rstat{Nresult}.Z(1,calcstart:calcend)); %#ok<AGROW>
                elseif isfield(Rstat{Nresult},'df')
                    Rstat2{Nresult}.df = cat(2,Rstat2{Nresult}.df,Rstat{Nresult}.df(1,calcstart:calcend),Rstat{Nresult}.df(1,calcstart:calcend)); %#ok<AGROW>
                end
                if isfield(Rstat2{Nresult},'mean')
                    for j = 1:length(Rstat2{Nresult}.mean)
                        Rstat2{Nresult}.mean{j} = cat(2,Rstat2{Nresult}.mean{j},Rstat{Nresult}.mean{j}(1,calcstart:calcend),Rstat{Nresult}.mean{j}(1,calcstart:calcend)); %#ok<AGROW>
                        Rstat2{Nresult}.std{j} = cat(2,Rstat2{Nresult}.std{j},Rstat{Nresult}.std{j}(1,calcstart:calcend),Rstat{Nresult}.std{j}(1,calcstart:calcend)); %#ok<AGROW>
                    end
                end
            end
        end
        header.vars = [header.vars header.vars];
    else
        for Nresult = 1:cfg.resultvars
            Rstat2{Nresult}.p = zeros(1,2*size(Rstat{Nresult}.p,2)); %#ok<AGROW>
            Rstat2{Nresult}.p(1:2:end) = Rstat{Nresult}.p(1:2:end); %#ok<AGROW>
            Rstat2{Nresult}.p(2:2:end) = Rstat{Nresult}.p(1:2:end); %#ok<AGROW>
            if isfield(Rstat{Nresult},'mxp')
                Rstat2{Nresult}.mxp = zeros(1,2*size(Rstat{Nresult}.mxp,2)); %#ok<AGROW>
                Rstat2{Nresult}.mxp(1:2:end) = Rstat{Nresult}.mxp(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.mxp(2:2:end) = Rstat{Nresult}.mxp(1:2:end); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpC')
                Rstat2{Nresult}.mxpC = zeros(1,2*size(Rstat{Nresult}.mxpC,2)); %#ok<AGROW>
                Rstat2{Nresult}.mxpC(1:2:end) = Rstat{Nresult}.mxpC(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.mxpC(2:2:end) = Rstat{Nresult}.mxpC(1:2:end); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpV')
                Rstat2{Nresult}.mxpV = zeros(1,2*size(Rstat{Nresult}.mxpV,2)); %#ok<AGROW>
                Rstat2{Nresult}.mxpV(1:2:end) = Rstat{Nresult}.mxpV(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.mxpV(2:2:end) = Rstat{Nresult}.mxpV(1:2:end); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mxpS')
                Rstat2{Nresult}.mxpS = zeros(1,2*size(Rstat{Nresult}.mxpS,2)); %#ok<AGROW>
                Rstat2{Nresult}.mxpS(1:2:end) = Rstat{Nresult}.mxpS(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.mxpS(2:2:end) = Rstat{Nresult}.mxpS(1:2:end); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'T')
                Rstat2{Nresult}.T = zeros(1,2*size(Rstat{Nresult}.T,2)); %#ok<AGROW>
                Rstat2{Nresult}.T(1:2:end) = Rstat{Nresult}.T(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.T(2:2:end) = Rstat{Nresult}.T(1:2:end); %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'R')
                Rstat2{Nresult}.R = zeros(1,2*size(Rstat{Nresult}.R,2)); %#ok<AGROW>
                Rstat2{Nresult}.R(1:2:end) = Rstat{Nresult}.R(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.R(2:2:end) = Rstat{Nresult}.R(1:2:end); %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'F')
                Rstat2{Nresult}.F = zeros(1,2*size(Rstat{Nresult}.F,2)); %#ok<AGROW>
                Rstat2{Nresult}.F(1:2:end) = Rstat{Nresult}.F(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.F(2:2:end) = Rstat{Nresult}.F(1:2:end); %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'Z')
                Rstat2{Nresult}.Z = zeros(1,2*size(Rstat{Nresult}.Z,2)); %#ok<AGROW>
                Rstat2{Nresult}.Z(1:2:end) = Rstat{Nresult}.Z(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.Z(2:2:end) = Rstat{Nresult}.Z(1:2:end); %#ok<AGROW>
            elseif isfield(Rstat{Nresult},'df')
                Rstat2{Nresult}.df = zeros(1,2*size(Rstat{Nresult}.dt,2)); %#ok<AGROW>
                Rstat2{Nresult}.df(1:2:end) = Rstat{Nresult}.df(1:2:end); %#ok<AGROW>
                Rstat2{Nresult}.df(2:2:end) = Rstat{Nresult}.df(1:2:end); %#ok<AGROW>
            end
            if isfield(Rstat{Nresult},'mean')
                for j = 1:length(Rstat{Nresult}.mean)
                    Rstat2{Nresult}.mean{j} = zeros(1,2*size(Rstat{Nresult}.mean{j},2)); %#ok<AGROW>
                    Rstat2{Nresult}.mean{j}(1:2:end) = Rstat{Nresult}.mean{j}(1:2:end); %#ok<AGROW>
                    Rstat2{Nresult}.mean{j}(2:2:end) = Rstat{Nresult}.mean{j}(1:2:end); %#ok<AGROW>
                    Rstat2{Nresult}.std{j} = zeros(1,2*size(Rstat{Nresult}.std{j},2)); %#ok<AGROW>
                    Rstat2{Nresult}.std{j}(1:2:end) = Rstat{Nresult}.std{j}(1:2:end); %#ok<AGROW>
                    Rstat2{Nresult}.std{j}(2:2:end) = Rstat{Nresult}.std{j}(1:2:end); %#ok<AGROW>
                end
            end
        end
        varstmp = cell(1,2*size(header.vars,2));
        varstmp(1:2:end) = header.vars;
        varstmp(2:2:end) = header.vars;
        header.vars = varstmp;
        clearvars varstmp
    end
    Rstat = Rstat2;
    clearvars Rstat2
    cfg.clustervarsI = cfg.clustervars;
    cfg.clustervars = 2*cfg.clustervars;
end

return