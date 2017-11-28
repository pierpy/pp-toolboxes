function lab_plot_IS_plotall(cfg)

zerocolor = [0.941176474094391,0.941176474094391,0.941176474094391];
[DATA,PLOT,PlotList,Brain] = lab_plot_prepare_data(cfg,zerocolor);
if isempty(DATA)
    disp(['Wrong input data, must be exactly ' num2str(size(Brain.labels,1)) ' values'])
    return
end
[DATA2,PLOT2] = lab_plot_mix_data(DATA,PLOT,PlotList);
[DATA2,PLOT2] = lab_plot_create_title(DATA2,PLOT2);

% define color of surfaces
for plotnr = 1:size(DATA2,1)
    clearvars facecolor facecolorL facecolorR PlotFaces PlotVertices PlotFacesR PlotVerticesR PlotFacesL PlotVerticesL
    Modes = fieldnames(DATA2(plotnr));
    Matrix = [];
    FlagM = false;
    PlotFaces{plotnr}{1} = Brain.faces; %#ok<AGROW>
    PlotVertices{plotnr}{1} = Brain.vertices; %#ok<AGROW>
    PlotFacesL{plotnr}{1} = Brain.facesL; %#ok<AGROW>
    PlotVerticesL{plotnr}{1} = Brain.verticesL; %#ok<AGROW>
    PlotFacesR{plotnr}{1} = Brain.facesR; %#ok<AGROW>
    PlotVerticesR{plotnr}{1} = Brain.verticesR; %#ok<AGROW>
    facecolor{plotnr}{1} = repmat(zerocolor,size(Brain.faces,1),1); %#ok<AGROW>
    facecolorR{plotnr}{1} = facecolor{plotnr}{1}(Brain.locR,:); %#ok<AGROW>
    facecolorL{plotnr}{1} = facecolor{plotnr}{1}(Brain.locL,:); %#ok<AGROW>
    Valpha = 1;
    for i = 1:length(Modes)
        if isfield(DATA2(plotnr).(Modes{i}),'data') & ~isempty(DATA2(plotnr).(Modes{i}).data) & ...
                (~isfield(PLOT2(plotnr).(Modes{i}),'Valid') | PLOT2(plotnr).(Modes{i}).Valid == true)
            if FlagM == false
                PLOT2(plotnr).(Modes{i}).AddPlot = false;
                FlagM = true;
            else
                PLOT2(plotnr).(Modes{i}).AddPlot = true;
            end
            datatmp = (DATA2(plotnr).(Modes{i}).data - PLOT2(plotnr).(Modes{i}).MinValue) / ...
                (PLOT2(plotnr).(Modes{i}).MaxValue - PLOT2(plotnr).(Modes{i}).MinValue);
            datatmp(datatmp < 0) = 0;
            datatmp(datatmp > 1) = 1;
            Color = zeros(length(datatmp),3);
            for C = 1:length(datatmp)
                Cidx = ceil(datatmp(C)*size(PLOT2(plotnr).(Modes{i}).Color,1));
                if isnan(Cidx)
                    Color(C,:) = [NaN NaN NaN];
                else
                    if Cidx == 0
                        Cidx = 1;
                    elseif Cidx > size(PLOT2(plotnr).(Modes{i}).Color,1)
                        Cidx = size(PLOT2(plotnr).(Modes{i}).Color,1);
                    end
                    Color(C,:) = PLOT2(plotnr).(Modes{i}).Color(Cidx,:);
                end
            end
            Color = cat(1,Color,zerocolor);
            if isfield(DATA2(plotnr).(Modes{i}),'Nr') & DATA2(plotnr).(Modes{i}).Nr > 1
                datatmp = mod(DATA2(plotnr).(Modes{i}).data,1);
            end
            if strcmp(Modes{i},'Surface')
                % MapsAll = Brain.mapsall .* repmat([datatmp(:)' 1],size(Brain.faces,1),1);
                MapsAll = Brain.mapsall;
                MapsAll(:,end) = ones(size(MapsAll,1),1)-sum(MapsAll(:,1:end-1),2);
                facecolor{plotnr}{1} = MapsAll * Color; %#ok<AGROW>
                facecolor{plotnr}{1}(facecolor{plotnr}{1} < 0) = 0; %#ok<AGROW>
                facecolorR{plotnr}{1} = facecolor{plotnr}{1}(Brain.locR,:); %#ok<AGROW>
                facecolorL{plotnr}{1} = facecolor{plotnr}{1}(Brain.locL,:); %#ok<AGROW>
                if isfield(PLOT2(plotnr).(Modes{i}),'Alpha')
                    Valpha(1) = min(PLOT2(plotnr).(Modes{i}).Alpha);
                end
            elseif strcmp(Modes{i},'Volume')
                Valpha(1) = 0.5;
                for regionnr = 1:length(datatmp)
                    if datatmp(regionnr) > 0
                        facecolor{plotnr}{end+1} = repmat(Color(regionnr,:),size(Brain.regions(regionnr).faces,1),1); %#ok<AGROW>
                        PlotFaces{plotnr}{end+1} = Brain.regions(regionnr).faces; %#ok<AGROW>
                        PlotVertices{plotnr}{end+1} = Brain.regions(regionnr).vertices; %#ok<AGROW>
                        if Brain.regionsLR{regionnr} == 'L'
                            facecolorL{plotnr}{end+1} = facecolor{plotnr}{end}; %#ok<AGROW>
                            PlotFacesL{plotnr}{end+1} = Brain.regions(regionnr).faces; %#ok<AGROW>
                            PlotVerticesL{plotnr}{end+1} = Brain.regions(regionnr).vertices; %#ok<AGROW>
                        elseif Brain.regionsLR{regionnr} == 'R'
                            facecolorR{plotnr}{end+1} = facecolor{plotnr}{end}; %#ok<AGROW>
                            PlotFacesR{plotnr}{end+1} = Brain.regions(regionnr).faces; %#ok<AGROW>
                            PlotVerticesR{plotnr}{end+1} = Brain.regions(regionnr).vertices; %#ok<AGROW>
                        end
                        if isfield(PLOT2(plotnr).(Modes{i}),'Alpha')
                            Valpha(1,end+1) = PLOT2(plotnr).(Modes{i}).Alpha(regionnr); %#ok<AGROW>
                        else
                            Valpha(1,end+1) = 1; %#ok<AGROW>
                        end
                    end
                end
            elseif strcmp(Modes{i},'Nodes')
                Matrix.nodes = datatmp;
                Matrix.xyz = Brain.xyz;
                Matrix.ncolor = Color;
                if isfield(PLOT2(plotnr).(Modes{i}),'Size')
                    Matrix.nsize = PLOT2(plotnr).(Modes{i}).Size;
                end
                if isfield(PLOT2(plotnr).(Modes{i}),'Alpha')
                    Matrix.alpha = PLOT2(plotnr).(Modes{i}).Alpha;
                end
                Valpha(1) = 0.3;
            elseif strcmp(Modes{i},'Connections')
                Matrix.matrixedges = Brain.matrixedges;
                Matrix.xyz = Brain.xyz;
                Matrix.edges = datatmp;
                Matrix.ecolor = Color;
                if isfield(PLOT2(plotnr).(Modes{i}),'Size')
                    Matrix.esize = PLOT2(plotnr).(Modes{i}).Size;
                end
                if isfield(PLOT2(plotnr).(Modes{i}),'Alpha')
                    Matrix.ealpha = PLOT2(plotnr).(Modes{i}).Alpha;
                end
                Valpha(1) = 0.3;
            end
        end
    end
    lab_plot_IS_plotsingle(PlotFaces{plotnr},PlotVertices{plotnr},facecolor{plotnr},PlotFacesR{plotnr}, ...
        PlotVerticesR{plotnr},facecolorR{plotnr},PlotFacesL{plotnr},PlotVerticesL{plotnr},facecolorL{plotnr},Matrix,Valpha,PLOT2(plotnr),cfg)
end

end