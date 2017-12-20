classdef mustates 
%     MUSTATES Microstates analysis class
%       Detailed explanation goes here
%       fissare una struttura per i dati in ingresso : una possibile
%       soggetto.cond{runs}
%       ripensarla solo come classe che contiene tutte le funzioni per fare l analisi, i cicicli su condizioni, sogg ecc si fanno fuori
    properties
        confstruct
    end   
    methods
        function [ self ] = mustates (confstruct)
            % funzione che legge tutta la cartella dei file e crea una cell
            % con tutti i soggetti e tutte le condizoni
                self.confstruct = confstruct;
                if ~ isfield(self.confstruct, 'substract_column_at_start')
                    self.confstruct.substract_column_at_start = 1;
                end
                if ~ isfield(self.confstruct, 'method_GFPeak')
                    self.confstruct.method_GFPeak = 'GFPL2';
                end
                if ~ isfield(self.confstruct, 'use_gfp_peaks')
                    self.confstruct.use_gfp_peaks = 1;
                end
                if ~ isfield(self.confstruct, 'setall1')
                    self.confstruct.setall1 = 1;
                end
                if ~ isfield(self.confstruct, 'normalize')
                    self.confstruct.normalize = 1;
                end
                if ~ isfield(self.confstruct, 'minclusters')
                    self.confstruct.minclusters = 1;
                end
                if ~ isfield(self.confstruct, 'maxclusters')
                    self.confstruct.maxclusters = 15;
                end
                if ~ isfield(self.confstruct, 'pmode')
                    self.confstruct.pmode = 0;
                end
                if ~ isfield(self.confstruct, 'algorithm')
                    self.confstruct.algorithm = 1;
                end
                if ~ isfield(self.confstruct, 'similarity_measure')
                    self.confstruct.similarity_measure = 'correlation';
                end
                if ~ isfield(self.confstruct, 'maxima_method')
                    self.confstruct.maxima_method = 'simple';
                end
                if ~ isfield(self.confstruct, 'debuG')
                    self.confstruct.debug = 0;
                end
                if ~ isfield(self.confstruct, 'minDistance')
                    self.confstruct.minDistance=0;
                end
                if ~ isfield(self.confstruct, 'minAmplitude')
                    self.confstruct.minAmplitude=0;
                end
                if ~ isfield(self.confstruct, 'showpeaks')
                    self.confstruct.showpeaks=0;
                end
                if ~ isfield(self.confstruct, 'downsample')
                    self.confstruct.downsample = 0;
                end
                if  ~ isfield(self.confstruct, 'segment_min_len')
                    self.confstruct.segment_min_len = 3;
                end
                if ~ isfield(self.confstruct, 'tf')
                    error('tf is a required field')
                end
                if  ~ isfield(self.confstruct, 'fs') 
                    error('fs is a required field')
                end
                
        end
        [eeg, gfp_peak_indices, gfp_curve]  = preprocess(self, eeg, nch);
        [eeg] = avgref(self, eeg, nch);
        [eeg] = detrendcolumnwise(self, eeg);
        [eeg] = setgfpall1(self, eeg, gfp_curve);
        [gfp_curve] = computegfp(self, eeg, method);
        [gfp_peak_indices, gfp_peak_values, gfp_curve] = computegfppeaks(self, gfp_curve);
        [cv] = computecv(self, data, templates, clust_ind, ClusterNr);
        [w] = computew(self, data, clust_ind, ClusterNr);
        [m] = computem(self, w, Nchans, ClusterNr);
        [kl] = computekl(self, mm);
        [Vectors,Values,Psi] = pcevectors(self, A, numvecs, verbose);
        [vectors values] = sortempp(self, vectors, values);
        [b_model, b_ind, b_loading, exp_var] = modifiedkmeans(self, eeg, n_mod, pmode, reruns);
        [b_model, ind, c, gfp] = aahc(self, eeg, n_mod, IgnorePolarity);
        [ clusteringResults ] = segmentation(self, data, gfp_peaks_indices );
        [ expVar, prototypes, statesSequence, cv, kl ] = extractsegmentationstruct( self, clusteringResults, ntemplatesToExtract );  
        [ output ] = computemsparameters( self, eeg, maps, nch);
        [ meanDur, meanOcc, gev, totalGev, meanCov, statesSequence, singleCorrelations ] = extractbackfittingstruct( self, backfittingResults );
    end
   
end

