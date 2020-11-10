clear;clc
addpath(pwd)
addpath('/usr/local/MATLAB/R2017a/toolbox/fMRItoolbox/')
addpath('/home/xiaojue/bin/nifti_tools/')
addpath('/data1/2018_ActionDecoding/pilot/analysis_class/')
%addpath('/data1/2018_ActionDecoding/pilot/analysis_fc/')
addpath('/data1/2018_ActionDecoding/analysis_fc')
addpath '/data1/2018_ActionDecoding/analysis_fc/Scripts/CalculateBeta/'

DataDir = '/data1/2018_ActionDecoding/data/';
OutFileDir = '/data1/2018_ActionDecoding/analysis_fc/DesignMat/GSR/';

DesignMnewLSSDir = '/data1/2018_ActionDecoding/analysis_fc/Despike/DesignM/LSS_noBeh/';
SubPathPrefix = '/bv';
behPathPrefix = '/beh/sess-02/';
GSRBetaTxtSuffix = '_GSRBetaBeh';

designMotionfNameSuffix = '_LSSDespike_GSR.mat';
NumVols = 202;
NumTrials = 24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       make a list of subjects to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen ('/data1/2018_ActionDecoding/analysis_fc/InputList/Subjects_list_FC_redo_02092020.txt','r');
file = textscan(fileID,'%q');
subList = file{1};
fclose(fileID);
NumSubs = length(subList);

for sub =  1:NumSubs 
 
    subID = char(subList(sub));
    %example: sub-07
    subID
    saveGSRBetaTxtFileName = strcat(subID,GSRBetaTxtSuffix,'.txt');
    subPath = strcat(DataDir,subID, SubPathPrefix);
    behPath = strcat(DataDir, subID,behPathPrefix);
    OutPrefix  = subID;
    
    designMotionfName = strcat(subID,'DespikeLSS_GSR.mat');
    load(fullfile(DesignMnewLSSDir,designMotionfName));
    
    %27:28 is the GSR signal 

    eventsList = dir(strcat(behPath,'*_events.tsv'));
    behavList = dir(strcat(behPath,'*_responses.tsv'));
    vtcList = dir(strcat(subPath,'/*actdecode*NATIVE.vtc'));
    NumVTCs = length(vtcList);
    
    OutGSRBetaBeh = zeros(NumTrials*NumVTCs,2+5);
    for scan = 1:NumVTCs
        %%
        %Getting all the files and timing here first 
        %
        indexBetaHatStart = 1+(scan-1)*NumTrials; 
        %trying to find which trial index in 192 the current Beta
        indexBetaHatEnd = NumTrials+(scan-1)*NumTrials;
        
        %events = importdata([behPath '/' eventsList(scan).name]);
        behav = importdata([behPath '/' behavList(scan).name]);
        OutGSRBetaBeh(indexBetaHatStart:indexBetaHatEnd,3:end) = behav.data(:,4:8);
        %column 4: instruction
        %column 5: actor
        %column 6: action
        %column 7: viewpoint
        %column 8: motion direction 
        %%%function1 Get timing for different event types
        GSRTS = OutDesignMatrix.Xlss{scan};
        GSRTS = GSRTS(:,27:28,1);
        
        Xlss = OutDesignMatrix.Xlss{scan}(:,1:2,:);
        Xlss(:,size(Xlss,2)+1,:) = 1;
        c = zeros(1,size(Xlss,2));
        c(1,1) = 1;
        
        
        for ig = 1:size(GSRTS,2)
            gsrin = GSRTS(:,ig);
            BetaHat = zeros(NumTrials,1);
            %calculate beta for GSR  for two sets of GSR
            for i = 1:size(Xlss,3)
                b = regress(gsrin,Xlss(:,:,i));
                BetaHat(i) = c*b;

                %BetaHat(i) = c*((X_lss(:,:,i)'*X_lss(:,:,i))\X_lss(:,:,i)'*Y);
            end 
            OutGSRBetaBeh(indexBetaHatStart:indexBetaHatEnd,ig) = BetaHat;
        end
        
    end
    %save the OutGSRBetaHat with correct header 
    if ~exist(strcat(OutFileDir,saveGSRBetaTxtFileName),'file')
        header = cellstr(['GSRBeta1','GSRBeta2',behav.textdata(1,5:9)]);
        fid = fopen(fullfile(OutFileDir,saveGSRBetaTxtFileName),'wt');
        fprintf(fid,'%s\t',header{1:end-1});
        fprintf(fid,'%s\n',header{end});
        fclose(fid);    
        dlmwrite(fullfile(OutFileDir,saveGSRBetaTxtFileName),OutGSRBetaBeh,'delimiter','\t','-append')
    else
        disp(['already saved' strcat(OutFileDir,saveGSRBetaTxtFileName)])
    end
    
end%sub