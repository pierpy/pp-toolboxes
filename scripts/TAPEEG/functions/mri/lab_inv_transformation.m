% Function to calculate inverse of a spm transformation (using spm_prep2sn)
%
% function  Tout = lab_inv_transformation(Tin,Source,Target)

function  [Tout,Tout2] = lab_inv_transformation(T,Source,Target)

if ischar(T)
    Tin = load(T);
else
    Tin = T;
end
if ischar(Source)
    Source = spm_vol(Source);
end
if ischar(Target)
    Target = spm_vol(Target);
end

P.tpm = Target;
P.image = Source;
P.Twarp = Tin.Tr;
P.Affine = Tin.Affine;

% empty data for flags
P.ngaus = [];
P.mg = [];
P.mn = [];
P.vr = [];
P.warpreg = [];
P.warpco = [];
P.biasreg = [];
P.biasfwhm = [];
P.regtype = [];
P.fudge = [];
P.samp = [];
P.msk = [];
P.Tbias = [];
P.thresh = [];

[Tout,Tout2] = spm_prep2sn(P);
Tout.flags = Tin.flags;
Tout2.flags = Tin.flags;

if ischar(T)
    [~,T_filepath,~,T_fileS] = lab_filename(T);
    
    VF = Tout.VF;
    VG = Tout.VG;
    Affine = Tout.Affine;
    Tr = Tout.Tr;
    flags = Tout.flags;
    Tout = fullfile(T_filepath,[T_fileS '_inv.mat']);
    save(Tout,'VG','VF','Affine','Tr','flags');
    
    VF = Tout2.VF;
    VG = Tout2.VG;
    Affine = Tout2.Affine;
    Tr = Tout2.Tr;
    flags = Tout2.flags;
    Tout2 = fullfile(T_filepath,[T_fileS '_inv2.mat']);
    save(Tout2,'VG','VF','Affine','Tr','flags');
end

end

