addpath('/usr/local/MATLAB/R2017a/toolbox/fMRItoolbox/')
addpath('/home/xiaojue/bin/nifti_tools/')
addpath('/data1/2018_ActionDecoding/pilot/analysis_class/')
addpath('/data1/2018_ActionDecoding/analysis_fc')
addpath('/data1/2018_ActionDecoding/analysis_fc/Scripts/Misc/')
DataDir = '/data1/2018_ActionDecoding/data/';
addpath('/data1/2018_ActionDecoding/analysis_fc/Misc/BME295/')
addpath('/data/home/xiaojue/bin/bsmart/')
OutFileDir = '/data1/2018_ActionDecoding/analysis_fc/Misc/BME295/';
saveGCTxtFileName = 'AllSubGCIFGpSTS_RH.txt';
fileID = fopen ('/data1/2018_ActionDecoding/analysis_fc/InputList/Subjects_list_FC_redo_02092020.txt','r');
file = textscan(fileID,'%q');
subList = file{1};
fclose(fileID);
NumSubs = length(subList);
SubPathPrefix = '/bv';
OutIFG2pSTS=zeros(NumSubs,1);
OutpSTS2IFG=zeros(NumSubs,1);

for sub =  1:NumSubs 
 
    subID = char(subList(sub));
    subPath = strcat(DataDir,subID, SubPathPrefix);
    %Granger causality 
    %run it on 
    SeedAtlasfName = strcat(subID,'_rh_cba.mat'); 
    ROITS = load(strcat(DataDir,subID,'/bv/',SeedAtlasfName));
    IFGRHfName = strcat(subID,'_rh_IFG.mat'); 
    IFGTS = load(strcat(DataDir,subID,'/bv/',IFGRHfName));
    ROITS = ROITS.Data;
    fil_TCIFG_RH = IFGTS.Data.patterns{1};
    %get two ROI of interest 
    fil_TCpSTS_RH = ROITS.patterns{3};
    numRun = size(fil_TCIFG_RH,3);
    %use this 
    tempdata = [reshape(mean(fil_TCpSTS_RH,2),[],1)'; reshape(mean(fil_TCIFG_RH,2),[],1)'];

    %after picking the order 5
    [~,Ex] =  armorf(tempdata(1,:),numRun,size(fil_TCIFG_RH,1),5);
    [~,Ey] =  armorf(tempdata(2,:),numRun,size(fil_TCIFG_RH,1),5);
    [~,E] =  armorf(tempdata,numRun,size(fil_TCIFG_RH,1),5);

    %causal estimates 
    OutIFG2pSTS(sub)=log(Ex/E(1,1));
    OutpSTS2IFG(sub)=log(Ey/E(2,2));
end

%save output to txt for used in R
if ~exist(strcat(OutFileDir,saveGCTxtFileName),'file')
    header = cellstr(['IFG2pSTS';'pSTS2IFG'])';
    fid = fopen(fullfile(OutFileDir,saveGCTxtFileName),'wt');
    fprintf(fid,'%s\t',header{1:end-1});
    fprintf(fid,'%s\n',header{end});
    fclose(fid);    
    dlmwrite(fullfile(OutFileDir,saveGCTxtFileName),[OutIFG2pSTS OutpSTS2IFG],'delimiter','\t','-append')
else
    disp(['already saved' strcat(OutFileDir,saveGCTxtFileName)])
end