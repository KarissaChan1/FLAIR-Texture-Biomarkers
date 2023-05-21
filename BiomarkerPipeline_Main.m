%% Biomarker Pipeline Main Function
% Run any biomarker analysis. Includes NABM texture features, wavelet
% analysis, and tract ROI analysis.

% Required: Experiment directory and function paths. Include Clinical
% sheets for corresponding datasets.

% Universal and Texture Function Paths
addpath('G:\My Drive\IAMLAB\MATLAB tools');
addpath('C:\Users\kchan\Documents\GitHub\Karissa.Chan\WM Tract Analysis Pipeline\Functions\');
addpath('C:\Users\kchan\Documents\GitHub\Karissa.Chan\WM Tract Analysis Pipeline\Functions\Wavelets\');
addpath('C:\Users\kchan\Documents\GitHub\Karissa.Chan\WM Tract Analysis Pipeline\Functions\Texture\');
addpath('C:\Users\kchan\Documents\GitHub\Karissa.Chan\WM Tract Analysis Pipeline\Functions\Penumbra\');
addpath('C:\Users\kchan\Documents\GitHub\Karissa.Chan\WM Tract Analysis Pipeline\Functions\Post-processing\');


% Directories

% Volume Directories
volDir = ('G:\.shortcut-targets-by-id\1okBQi9CCZwdTAAEzFuGM1P61Lj8qTM7V\Karissa Data\08-31-2021\CCN2021 v4 Registered\Volumes\');
wmlDir = 'G:\.shortcut-targets-by-id\1okBQi9CCZwdTAAEzFuGM1P61Lj8qTM7V\Karissa Data\08-31-2021\CCN2021 v4 Registered\WML Masks\';

% Results Directories
dataDir = 'C:\Users\kchan\Documents\Biomarker Pipeline\CCNA2021_V4_Experiment2\';
mkdir(dataDir);


%% NABM TEXTURE BIOMARKERS

mkdir([dataDir,'Masks\NABM\']);
mkdir([dataDir,'Masks\NABM NIFTI\']);
mkdir([dataDir,'Texture\MAD\FeaturesMaps\']);
mkdir([dataDir,'Texture\MII\FeaturesMaps\']);
mkdir([dataDir,'Texture\MID\FeaturesMaps\']);
mkdir([dataDir,'Texture\MAD\MeanFeatureMaps\']);
mkdir([dataDir,'Texture\MII\MeanFeatureMaps\']);
mkdir([dataDir,'Texture\MID\MeanFeatureMaps\']);

% Physical Dimension being analyzed
physical_dim = 3*0.8594;

% Files
files = dir(volDir);
files(1:2)=[];
    
% NABM
for k = 1:length(files)
    for j=1
        
    % CCNA
    load(fullfile(volDir,files(k).name));
    load(fullfile(wmlDir,files(k).name));

    % Matching window size to the physical dimension
    window_size = physical_dim;
      
    % Width Calculation FLAIR
    pixel_size = 0.8594;
    minimum_pixel = window_size.*(1./(pixel_size));
    width = ceil(minimum_pixel);
    rem = mod(width,2);
    if rem == 0
        width = width + 1;          
    end

    % NABM Extraction
    %%% already brain extracted %%%
    [brainClass] = getNABM_wpm(regvol,'wml',wmlMask); % for registered volumes

    % Texture and Volume Measures
    [corr, LBP_img] = TextureAnalysisFinal(brainClass, width);
    Final_window(k,j) = width;
    Min_dim_req(k,j) = minimum_pixel;
    Actual_dim(k,j) = width*(pixel_size);

    % FLAIR Masks
    save(fullfile(dataDir, 'Masks\NABM\',files(k).name), 'brainClass');
%     niftiwrite(brainClass,fullfile(dataDir, 'Masks\NABM NIFTI\',files(k).name));

    % FLAIR Texture
    save(fullfile(dataDir, 'Texture\MAD\FeaturesMaps\',files(k).name), 'corr');
    save(fullfile(dataDir, 'Texture\MII\FeaturesMaps\',files(k).name), 'LBP_img');
    
    [med_mvf_MII, volMean_Integrity, medC_mvf_MAD, volMean_Macro, medC_mvf_MID, Micro_FM, volMean_Micro] = Post_processing_Final(corr,LBP_img, width);
    
    % FLAIR
    Median_mvf_MII(k,j) =  med_mvf_MII;
    Median_mvf_MAD(k,j) = medC_mvf_MAD;
    Median_mvf_MID(k,j) = medC_mvf_MID;
    save([dataDir, 'Texture\Texture_Metrics.mat'], 'Median_mvf_MII','Median_mvf_MAD','Median_mvf_MID');

    % Micro FM Vol
    save(fullfile(dataDir, 'Texture\MID\FeaturesMaps\', files(k).name), 'Micro_FM');
    
    % FM Imgs
    save(fullfile(dataDir, 'Texture\MID\MeanFeatureMaps\', files(k).name), 'volMean_Micro');
    save(fullfile(dataDir, 'Texture\MII\MeanFeatureMaps\', files(k).name), 'volMean_Integrity');
    save(fullfile(dataDir, 'Texture\MAD\MeanFeatureMaps\', files(k).name), 'volMean_Macro');    
    
    % compute NABM intensity
    Int_NABM(k,j) = median(brainClass(brainClass(:)>0));

    
    
    k
   
    end
end
Feat_Table = table(Median_mvf_MII,Median_mvf_MAD,Median_mvf_MID, Int_NABM);
idx = find(~Median_mvf_MII);
Feat_Table(idx,:)=[];
Feat_Table.Properties.VariableNames = {'Integrity_LBP' 'Damage_Macro' 'Damage_Micro' 'Intensity'};
save([dataDir, 'NABM_Metrics_Table.mat'], 'Feat_Table');


%% WAVELET BIOMARKERS
addpath('G:\My Drive\IAMLAB\MATLAB tools\');
addpath('C:\Users\karis\Documents\GitHub\Karissa.Chan\WM Tract Analysis Pipeline\Functions\Wavelets\');
addpath('C:\Users\karis\Documents\GitHub\Karissa.Chan\WM Tract Analysis Pipeline\Functions\Texture\');
addpath('C:\Users\karis\Documents\GitHub\Karissa.Chan\WM Tract Analysis Pipeline\Functions\Post-processing\');

nabmdir = fullfile(dataDir,'Masks\NABM\');
mkdir([dataDir,'Wavelets\Level3\']);

mkdir([dataDir,'Wavelets\MID\FeaturesMaps']);
mkdir([dataDir,'Wavelets\MID\MeanFeatureMaps']);


clear avg_vol_energy avg_vol_entropy
files = dir(nabmdir);
files(1:2)=[];

% Get physical dimensions of FLAIR vols
physical_dim = 3*0.8594;
window_size = physical_dim;
pixel_size = 0.8594;
minimum_pixel = window_size.*(1./(pixel_size));
width = ceil(minimum_pixel);
rem = mod(width,2);
if rem == 0
    width = width + 1;
end
File_Name = strings(length(files),1);
for j = 1:length(files)
    load(fullfile(nabmdir,files(j).name));
    
    f = files(j).name;
    [~,File_Name(j,1),~]=fileparts(f);
    
    for a = 1:size(brainClass,3)
        
        [ca,chd,cvd,cdd] = swt2(brainClass(:,:,a),3,'haar');

        A3 = wcodemat(ca(:,:,3),255);
        H3 = wcodemat(chd(:,:,3),255);
        V3 = wcodemat(cvd(:,:,3),255);
        D3 = wcodemat(cdd(:,:,3),255);

        approx_holder3(:,:,a) = A3;
        hor_holder3(:,:,a) = H3;
        vert_holder3(:,:,a) = V3;
        diag_holder3(:,:,a) = D3;
        
        % Compute Wavelet Statistical features
        [energy(a,1),entropy(a,1)]=waveletStatisticalFeatures(A3);
        
    end
    
    % Store Wavelet Statistical Features
    avg_vol_energy(j,1)=mean(energy(:));
    avg_vol_entropy(j,1)=mean(entropy(:));
    Mean_A3(j,1) = mean(approx_holder3(approx_holder3(:)>1));
    
    vol_decomp.A3 = approx_holder3;
    vol_decomp.H3 = hor_holder3;
    vol_decomp.V3 = vert_holder3;
    vol_decomp.D3 = diag_holder3;
    save(fullfile(dataDir,'Wavelets\Level3\',files(j).name),'vol_decomp');
     
    % Compute texture features from A3 volumes
    [corr, LBP_img] = TextureAnalysisFinal(approx_holder3, width);
    [median_MID(j,1),Micro_FM,volMean_Micro]=getWaveletMID(corr);

    save(fullfile(dataDir,'Wavelets\MID\FeaturesMaps',files(j).name),'Micro_FM');
    save(fullfile(dataDir,'Wavelets\MID\MeanFeatureMaps',files(j).name),'volMean_Micro');
    save([dataDir, 'Wavelets\Wavelet_Metrics.mat'], 'avg_vol_entropy','avg_vol_energy','median_MID','Mean_A3');

    j
end


Feat_Table = table(File_Name,avg_vol_entropy,avg_vol_energy,median_MID,Mean_A3);
Feat_Table.Properties.VariableNames = {'File_Name' 'Entropy_A3' 'Energy_A3' 'Damage_A3' 'Mean_A3'};
save([dataDir, 'Wavelet_Metrics_Table.mat'], 'Feat_Table');


%% WHITE MATTER LESION VOLUMES + PENUMBRA ROI
corrDir = fullfile(dataDir,'Texture\MAD\FeaturesMaps\');
lbpDir = fullfile(dataDir,'Texture\MII\FeaturesMaps\');
a3Dir = fullfile(dataDir,'Wavelets\Level3\');
nabmDir = fullfile(dataDir,'Masks\NABM\');

saveMasks = fullfile(dataDir,'Masks\Penumbra\');
mkdir(saveMasks);

% registered volume resolution
SpacingX = 0.8594;
SpacingY = 0.8594;
SliceThickness = 3.5;

clear penum_mad penum_mii penum_a3 penum_int

files = dir(wmlDir);
files(1:2)=[];
N = 5;
for j = 1:length(files)     
    load(fullfile(wmlDir,files(j).name));
    num_vox = sum(wmlMask,'all');
    WML_Volume(j,1) = (num_vox*SpacingX*SpacingY*SliceThickness)/1000;
    
    load(fullfile(corrDir,files(j).name));
    load(fullfile(lbpDir,files(j).name));
    load(fullfile(a3Dir,files(j).name));
    load(fullfile(nabmDir,files(j).name));
    
    
    [penumbra_masks,penumbra_features]=getPenumbraFeatures(wmlMask,N,brainClass,corr,LBP_img,vol_decomp.A3);
    save(fullfile(saveMasks,files(j).name),'penumbra_masks');

    for k = 1:N
        penum_mad(j,k)=penumbra_features{1,k}.corr;
        penum_mii(j,k)=penumbra_features{1,k}.lbp;
        penum_a3(j,k)=penumbra_features{1,k}.a3_mean;
        penum_int(j,k)=penumbra_features{1,k}.intensity;
        
    end
    
    j
end

save([dataDir, 'WML_Volume.mat'], 'WML_Volume');

for k = 1:N
    region_names{1,k} = convertStringsToChars(strcat('P',string(k)));
end
penum_mad = array2table(penum_mad);
penum_mad.Properties.VariableNames = region_names;

penum_mii = array2table(penum_mii);
penum_mii.Properties.VariableNames = region_names;

penum_a3 = array2table(penum_a3);
penum_a3.Properties.VariableNames = region_names;

penum_int = array2table(penum_int);
penum_int.Properties.VariableNames = region_names;

save([dataDir, 'Penumbra_Features.mat'],'penum_mad','penum_mii','penum_a3','penum_int');

%% Blood Supply Territory Features
clear bst_mad bst_mii bst_a3 bst_mid


regionDir = 'G:\Shared drives\_NeuroMRI_ValDB\Atlases\FLAIR\BST\v4_final\regions\';
corrDir = fullfile(dataDir,'Texture\MAD\FeaturesMaps\');
lbpDir = fullfile(dataDir,'Texture\MII\FeaturesMaps\');
a3Dir = fullfile(dataDir,'Wavelets\Level3\');
mida3Dir = fullfile(dataDir,'Wavelets\MID\FeaturesMaps\');

regionlist = {'ACA','MCA','PCA'};

files = dir(corrDir);
files(1:2)=[];
clear bst_mad bst_mii bst_a3 bst_mid

for j = 1:length(files)     

    load(fullfile(corrDir,files(j).name));
    load(fullfile(lbpDir,files(j).name));
    load(fullfile(a3Dir,files(j).name));
    load(fullfile(mida3Dir,files(j).name));
    
    for k = 1:length(regionlist)
        load(fullfile(regionDir,regionlist{k}+".mat"));
        mask = eval(regionlist{k});
        
        corr_region = corr .*mask;
        lbp_region = LBP_img .*mask;
        a3_region = vol_decomp.A3 .*mask;
        mid_region = Micro_FM.*mask;
        
        [bst_mii(j,k), volMean_Integrity, bst_mad(j,k), volMean_Macro,bst_mid(j,k),volMean_Micro] = Post_processing_BST(corr_region,lbp_region,mid_region);
        bst_a3(j,k) = median(a3_region(a3_region(:)>0));

        
    end
    
    j
end

bst_mad = array2table(bst_mad);
bst_mad.Properties.VariableNames = regionlist;

bst_mii = array2table(bst_mii);
bst_mii.Properties.VariableNames = regionlist;

bst_a3 = array2table(bst_a3);
bst_a3.Properties.VariableNames = regionlist;

bst_mid = array2table(bst_mid);
bst_mid.Properties.VariableNames = regionlist;

save([dataDir, 'BST_Features.mat'],'bst_mad','bst_mii','bst_a3','bst_mid');

%% WM TRACT BIOMARKER EXTRACTION
% use LBP (mean), A3 (energy and mean), and MAD (median)
clear TractA3Energy TractA3Median TractA3Mean TractMAD TractMII TractA3MID

texDir = fullfile(dataDir,'Texture');
wmlDir = 'G:\.shortcut-targets-by-id\1okBQi9CCZwdTAAEzFuGM1P61Lj8qTM7V\Karissa Data\08-31-2021\CCN2021 v4 Registered\WML Masks\';

wvDir = fullfile(dataDir,'\Wavelets\Level3\');
wvmidDir = fullfile(dataDir,'\Wavelets\MID\FeaturesMaps\');

tractDir = fullfile(dataDir, 'Masks\WM Tracts_FINAL\');
patientid = dir(tractDir);
patientid(1:2)=[];

textures = {'MAD','MII'};
mkdir(fullfile(dataDir,'Tract Texture Maps\MAD\'));
mkdir(fullfile(dataDir,'Tract Texture Maps\MII\'));
mkdir(fullfile(dataDir,'Tract Texture Maps\A3 MID\'));


% T = readtable('G:\My Drive\IAMLAB\Data\CAIN2\FLAIR_DTI_Matched.xlsx');
T = readtable('G:\My Drive\IAMLAB\Data\CCNA_2021\FLAIR_DTI_Matched.xlsx');

File_Name = strings([height(T),1]);
f = 0;

for k = 1:height(T)
    if ~strcmp(T.FA(k),"Missing") && ~strcmp(T.MD(k),"Missing")
        f = f+1;
        File_Name(f,1)=T.FLAIRNames(k);
        for t = 1:length(textures)
            mapsDir = fullfile(texDir,textures{t},'FeaturesMaps\');

            load(fullfile(mapsDir,T.FLAIRNames{k}));
            load(fullfile(wvDir,T.FLAIRNames{k}));
            load(fullfile(wvmidDir,T.FLAIRNames{k}));
            load(fullfile(wmlDir,T.FLAIRNames{k}));


            tractfolder = fullfile(tractDir,T.PatientID{k}); 
%             tractfolder = fullfile(tractDir,T.DTINames{k}); %T.DTINames for CAIN
            tract = dir(tractfolder);
            tract(1:2)=[];

            for j = 1:length(tract)

                load(fullfile(tractfolder,tract(j).name));

                %%% GET WAVELET FEATURES %%%
                [TractA3Energy(f,j),TractA3Median(f,j),TractA3Mean(f,j)]=waveletTractFeatures(vol_decomp.A3,tractMask);
                a3_mid_seg = tractMask.*Micro_FM;
                txseg = pixel_wise_avg(a3_mid_seg);
                TractA3MID(f,j)=median(txseg(txseg(:)>0));
                volMean_Micro_tr{j} = txseg;
                save(fullfile(dataDir,'Tract Texture Maps\A3 MID\',File_Name(f)),'volMean_Micro_tr');
                          
                %%% GET TEXTURE FEATURES %%%
                if t==1
                    txseg = tractMask.*corr;
                    txseg = pixel_wise_avg(txseg);
                    TractMAD(f,j) = median(txseg(txseg(:)>0));
                    volMean_Macro_tr{j} = txseg;
                    
                    save(fullfile(dataDir,'Tract Texture Maps\MAD\',File_Name(f)),'volMean_Macro_tr');
                elseif t==2
                    txseg = tractMask.*LBP_img;
                    txseg = pixel_wise_avg(txseg);
                    TractMII(f,j)=mean(txseg(txseg(:)>0));
                    volMean_Integrity_tr{j} = txseg;
                    
                    save(fullfile(dataDir,'Tract Texture Maps\MII\',File_Name(f)),'volMean_Integrity_tr');
                end
                
                %%% GET WML VOLUMES %%%
                wml_intersect = tractMask.*wmlMask;
                WML_Tract_Ratio(f,j) = 100*(sum(wml_intersect(:))/sum(tractMask(:)));
                
            
            end
        end
    end
    k
end

save(fullfile(dataDir,'Tract_Metrics_ALL.mat'),'File_Name','TractA3Energy','TractA3Median','TractA3Mean','TractMAD','TractMII','TractA3MID','WML_Tract_Ratio');

%% Add Clinical Data to NABM FEATURES: CCNA or CAIN options
dataset = "CCNA";

wv = load([dataDir,'Wavelet_Metrics_Table.mat']);
tx = load([dataDir,'NABM_Metrics_Table.mat']);
load([dataDir,'WML_Volume.mat']);
p = load([dataDir,'Penumbra_Features.mat']);
tr = load([dataDir,'Tract_Metrics_ALL.mat']);
bst = load([dataDir,'BST_Features.mat']);

addClinicalData(tx.Feat_Table,wv.Feat_Table,WML_Volume,p,tr,bst,dataDir,dataset)

