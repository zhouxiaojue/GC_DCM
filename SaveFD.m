SublistTxt = '/data1/2018_ActionDecoding/analysis_fc/InputList/Subjects_list_FC_redo_02092020.txt';
addpath('/data1/2018_ActionDecoding/analysis_fc/');
DataDir = '/data1/2018_ActionDecoding/data/';
fileID = fopen (SublistTxt,'r');
file = textscan(fileID,'%q');
subList = file{1};
fclose(fileID);

NumSubs = length(subList);
radius = 50;
OutFileDir = '/data1/2018_ActionDecoding/analysis_fc/MotionAnalysis/IndFD/';

for sub = 1:NumSubs
    subID = char(subList(sub));
    subPath = strcat(DataDir,subID, '/bv/');
    sdmList = dir(strcat(subPath,'/*actdecode*run*3DMC.sdm'));
    outFD = zeros(length(sdmList)*202,2);
    OutFDFileName = strcat(subID,'_FD.txt');
    for scan = 1:length(sdmList)
        sdmFileN = fullfile(subPath,sdmList(scan).name);
        sdm = readBvSDM(sdmFileN,6);     %read in sdm file
        nVols = size(sdm,1);
        %read in all six columns by 202 times points motion parameters 
        sdm(:,4:6,:) = sdm(:,4:6,:)*(pi/180); %First convert deg to rad 
        sdm(:,4:6,:) = sdm(:,4:6,:)*radius;   %... arc length (in mm) 

        fd = getFwd(sdm);
        outFD((1+202*(scan-1)):(202+202*(scan-1)),1) = fd;
        outFD((1+202*(scan-1)):(202+202*(scan-1)),2) = scan;

    end 
    %write out subject cumulative FD here 
    if ~exist(strcat(OutFileDir,OutFDFileName),'file')
        header = cellstr(['FD ';'run'])';
        fid = fopen(fullfile(OutFileDir,OutFDFileName),'wt');
        fprintf(fid,'%s\t',header{1:end-1});
        fprintf(fid,'%s\n',header{end});
        fclose(fid);    
        dlmwrite(fullfile(OutFileDir,OutFDFileName),outFD,'delimiter','\t','-append')
    else
        disp(['already saved' strcat(OutFileDir,OutFDFileName)])
    end
end %sub
