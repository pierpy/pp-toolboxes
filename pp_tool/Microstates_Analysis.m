function [NMicr] = Microstates_Analysis(algorithm, subject, risposta, folder_path)

        % SegmentationResult : structure with segmentation results
        % algorithm : for specifying the alghorihm to be used in the
        %               segmentation  
        %               ---> 'ModKMeans': modified version of k-means algorithm
        %                    like in Pascual-Marqui, 1995
        %               ---> 'AAHC': atomize and agglomerate hierarchical
        %                    clustering algorithm like in Cartool
        %               ---> 'SimpleKMeans': simple k-means algorithm like in
        %                    MatLab
        % subject_to_analyze : number of subject to analyze (content in the selected folder)
        
        % this function allows to load files to be analized (.txt or .ep)
        % and allows to do some preprocessing to data through a inputdlg.
        % the folder containing data must contain one folder for each
        % subject with .txt or .ep files.
        cd(folder_path);
        subjects = dir;
        subjects(1)=[];
        subjects(1)=[];
        
        h_filtering = str2num(risposta{1});
        l_filtering = str2num(risposta{2});
        freq_notch = str2num(risposta{3});
        dec=str2num(risposta{4});
        pmode = str2num(risposta{5});
        minclusters = str2num(risposta{6});
        maxclusters = str2num(risposta{7});
        zscore = str2num(risposta{8});

      

%        for subject = 1:subjects_to_analyze%
            disp(strcat('subject: ', num2str(subject)));
            path = strcat(folder_path,'\', subjects(subject).name);
            cd(path)
                 files=dir('*.interpolate.ep');
            %files = dir;
%             files(1)=[];
%             files(1)=[];
            total_data=[];

            for k=1:size(files,1)
                segment = load(files(k).name);
                total_data=[total_data;segment];
            end

            % any filtering
            if ~isempty(h_filtering) && ~isempty(l_filtering)
                disp('filtering...')
                [total_data] = eegfilt(total_data',500,h_filtering,l_filtering);
                total_data = total_data';
            end
            %  downsample
            if ~isempty(dec)
                disp('-------decimating...')
                total_data = downsample(total_data,dec);
            end
            % z score
            if zscore
                disp('-------zscore...')
                total_data = zscore(total_data);
            end
            
            % GFP maxima
            
            [data, GFP] = globalFieldPower(total_data');
            data=data';
            ClusterNr = minclusters:maxclusters;
            numclusters = length(ClusterNr);
            [Ntr, Nchans]=size(data);


            for q=1:numclusters
                disp(strcat('------cluster: ', num2str(q)));

                switch algorithm
                    case 'ModKMeans' % Nmicrostates alg.
                        disp('modified version of k-means')
                         [mps,ind,assigned_cluster_ampl,exp_var] = NMicrostates(data,ClusterNr(q),pmode,20);
                         NMicr.exp_var(q)=exp_var;
                    case 'AAHC'% AAHC alg.
                        disp('modified version of AAHC')
                         [mps,ind, c, gfp] = AAHC(data,ClusterNr(q),true);
%                           NMicr.exp_var(q) = (sum(c*gfp'))/(sum(gfp.^2));
                    case 'SimpleKMeans'% simple k-means alg
                        disp('simple k-means')
                         [ind, mps]=kmeans(data,q);
                    otherwise
                        disp('other value')
                end

                if numclusters ~= 1
                    % calcolo CV
                    [NMicr.CV(q)] = CV(data,mps,ind,ClusterNr(q));
                    % calcolo W (dispersion measure)
                    [NMicr.W(q)] = W(data,ind,ClusterNr(q));
                    [NMicr.m(q)] = M(W(data,ind,ClusterNr(q)),Nchans,ClusterNr(q));
                end
                 
                NMicr.template{q}=mps;
                NMicr.clusters{q}=ind;
                NMicr.original_data = total_data;
                NMicr.gfp_maxima_data = data;
                    % NMicr.amplitude{q}=assigned_cluster_ampl;
            end
            if numclusters ~= 1
                NMicr.kl  =  KL(NMicr.m);
            end
%             SegmentationResult{subject} = NMicr;
            clearvars files segment
            disp(strcat('-------------------------end subject:  ', num2str(subject), '-------------------'));
            cd(folder_path);
%         end

end