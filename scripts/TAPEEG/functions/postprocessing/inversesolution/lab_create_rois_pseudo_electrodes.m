% Read Cartool .rois or atlas(.hdr / eg AAL-atlas)
%
% [ROIS] = lab_create_rois_pseudo_electrode(LOCS,SPI)
%
% LOCS         = structures variable with location info of electrode
% locations
% SPI          = structures variable with location info of solutionpoints
%                (=lab_read_spi); only needed for atlas conversion
%
% written by G. Bogaarts 2016
function [ROIS] = lab_create_rois_pseudo_electrodes(LOCS,SPI)
    
ROIS.labels=LOCS.labels;
ROIS.numsolutionpts=length(SPI.x);
ROIS.numrois=length(LOCS.x);
ROIS.solutionptsAll = [];

distances=pdist2([SPI.x; SPI.y; SPI.z]',[LOCS.x; LOCS.y; LOCS.z]');
[~,ind]=sort(distances);
ROIS.solutionpts=num2cell(ind(1,:));
ROIS.solutionptsAll = ind(1,:);
ROIS.solutionptsAll = sort(ROIS.solutionptsAll);
ROIS = lab_create_rois_shortlabel(ROIS);