addpath('/usr/local/MATLAB/R2017a/toolbox/fMRItoolbox/')
addpath('/home/xiaojue/bin/nifti_tools/')
addpath('/data1/2018_ActionDecoding/pilot/analysis_class/')
addpath('/data1/2018_ActionDecoding/analysis_fc')
addpath('/data1/2018_ActionDecoding/analysis_fc/Scripts/Misc/')
DataDir = '/data1/2018_ActionDecoding/data/';
addpath('/data1/2018_ActionDecoding/analysis_fc/Misc/BME295/')
addpath('/data/home/xiaojue/bin/bsmart/')
subID = 'sub-08';
%this is the script for BME 295 final project which includes band-pass
%filtering and granger-causality (or other connectivity measure that's
%supposed to be interesting) 

%need field map correction for one functional run, either task run or
%resting state connectivity scan 




%%
% %simulation block 
% %load design matrix 
% designPath = '/data1/2018_ActionDecoding/analysis_fc/Despike/DesignM/NewLSS/';
% load(strcat(designPath,'sub-08DespikeLSS_all.mat'));
% %loaded outDesignMatrix 
% %add up the first 24 as hrf convolved test signal 
% Ytest = zeros(1,202*8);
% for i=1:8
%     Xlss = OutDesignMatrix.Xlss{i};
%     Ytest((1+202*(i-1)):(202+202*(i-1))) = sum(Xlss(:,1,:),3);
% end
% %and concatenate into long by 8 run signal 
% 
%add noise 
% Y_high = Ytest + normrnd(0,0.4,1,length(Ytest));
% Y_med = Ytest + normrnd(0,0.15,1,length(Ytest));
% Y_low = Ytest + normrnd(0,0.05,1,length(Ytest));
% %Ytest signal std = 0.2
% %plot to check simulated fMRI data 
% 
% figure;
% subplot(4,1,1)
% plot(1:length(Ytest),Ytest);
% title('test signal')
% ylabel('Simulated BOLD signal');
% subplot(4,1,2)
% plot(1:length(Ytest),Y_low);
% title('low noise')
% ylabel('Simulated BOLD signal');
% subplot(4,1,3)
% plot(1:length(Ytest),Y_med);
% ylabel('Simulated BOLD signal');
% title('medium noise')
% subplot(4,1,4)
% plot(1:length(Ytest),Y_high);
% title('high noise');
% xlabel('Volumes');
% ylabel('Simulated BOLD signal');
% %find trial index 
% 
% 
% 
% %Simulation for GC analysis 
% 
% %low noise data 
% %shift data from 2-7 data points
% orderall = 2:7;
% OutBIClow = zeros(length(orderall));
% OutGCxylow = zeros(length(orderall));
% OutGCyxlow = zeros(length(orderall));
% 
% OutBICmed = zeros(length(orderall));
% OutGCmed = zeros(length(orderall));
% OutGCxymed = zeros(length(orderall));
% OutGCyxmed = zeros(length(orderall));
% 
% OutBIChigh = zeros(length(orderall));
% OutGChigh = zeros(length(orderall));
% OutGCxyhigh = zeros(length(orderall));
% OutGCyxhigh = zeros(length(orderall));
% 
% for ord = 1:length(orderall)
%     order = orderall(ord);
%     Yxtest_low  = circshift(Ytest,order) + normrnd(0,0.05,1,length(Ytest));
%     Yxtest_med  = circshift(Ytest,order) + normrnd(0,0.15,1,length(Ytest));
%     Yxtest_high  = circshift(Ytest,order) + normrnd(0,0.4,1,length(Ytest));
% 
%     for p = 1:length(orderall)
%         timeord = orderall(p);
%         [~,Elow] = armorf([Y_low;Yxtest_low],8,length(Xlss),timeord);
%         OutBIClow(ord,p) = log(det(Elow)) + (log(length(Y_low))*timeord*2^2)/length(Y_low);
%         [~,Ex] =  armorf(Y_low,8,length(Xlss),5);
%         [~,Ey] =  armorf(Yxtest_low,8,length(Xlss),5);
%         %causal estimates 
%         OutGCyxlow(ord,p)=log(Ex/Elow(1,1));
%         OutGCxylow(ord,p)=log(Ey/Elow(2,2));
%         
%         [~,Emed] = armorf([Y_med;Yxtest_med],8,length(Xlss),timeord);
%         OutBICmed(ord,p) = log(det(Emed)) + (log(length(Y_med))*timeord*2^2)/length(Y_med);
%         [~,Ex] =  armorf(Y_med,8,length(Xlss),5);
%         [~,Ey] =  armorf(Yxtest_med,8,length(Xlss),5);
%         OutGCyxmed(ord,p)=log(Ex/Emed(1,1));
%         OutGCxymed(ord,p)=log(Ey/Emed(2,2));
%         
%         [~,Ehigh] = armorf([Y_high;Yxtest_high],8,length(Xlss),timeord);
%         OutBIChigh(ord,p) = log(det(Ehigh)) + (log(length(Y_high))*timeord*2^2)/length(Y_high);
%         [~,Ex] =  armorf(Y_high,8,length(Xlss),5);
%         [~,Ey] =  armorf(Yxtest_high,8,length(Xlss),5);
%         OutGCyxhigh(ord,p)=log(Ex/Ehigh(1,1));
%         OutGCxyhigh(ord,p)=log(Ey/Ehigh(2,2));
%     end
% end
% 
% figure(1);
% for i=1:6
%     subplot(6,2,1+2*(i-1))
%     plot(orderall,OutBIClow(i,:),'-o')
%     ylabel(['shift ' num2str(i+1) 'BIC'])
%     xlabel('fitted order')
%     subplot(6,2,2+2*(i-1))
%     plot(orderall,OutGCxylow(i,:),'-o')
%     hold on;
%     plot(orderall,OutGCyxlow(i,:),'-o')
%     hold off;
%     legend('Y on shifted Y','shifted Y on Y','Location','northwest')
%     xlabel('fitted order')
%     ylabel('GC index')
% end
% 
% figure(2);
% for i=1:6
%     subplot(6,2,1+2*(i-1))
%     plot(orderall,OutBICmed(i,:),'-o')
%     ylabel(['shift ' num2str(i+1) 'BIC'])
%     xlabel('fitted order')
%     subplot(6,2,2+2*(i-1))
%     plot(orderall,OutGCxymed(i,:),'-o')
%     hold on;
%     plot(orderall,OutGCyxmed(i,:),'-o')
%     hold off;
%     legend('Y on shifted Y','shifted Y on Y','Location','northwest')
%     xlabel('fitted order')
%     ylabel('GC index')
% end
% 
% figure(3);
% for i=1:6
%     subplot(6,2,1+2*(i-1))
%     plot(orderall,OutBIChigh(i,:),'-o')
%     ylabel(['shift ' num2str(i+1) 'BIC'])
%     xlabel('fitted order')
%     subplot(6,2,2+2*(i-1))
%     plot(orderall,OutGCxyhigh(i,:),'-o')
%     hold on;
%     plot(orderall,OutGCyxhigh(i,:),'-o')
%     hold off;
%     legend('Y on shifted Y','shifted Y on Y','Location','northwest')
%     xlabel('fitted order')
%     ylabel('GC index')
% end
%
%%
% %loading the data (fmr, unfiltered)
% vtcNonFiltered = xff('/data1/2018_ActionDecoding/analysis_fc/Misc/BME295/sub-08_actdecode_run-5_undist_NATIVE.vtc');
% vtcNonFiltered = vtcNonFiltered.VTCData;
% 
% %Extract ROI TS
% %pSTS_RH, IFG_RH
% vox = strcat(DataDir,subID,'/bv/','sub-08_rh_cba_VOI2Mat.mat');
% load(vox);
% %ROI.Data.MatlabIn are all the coordinates for voxels corresponding to the Time
% %series 
% %choose pSTS_RH (3rd one)
% 
% TCpSTS_RH = ExtractTS_mat(vtcNonFiltered,ROI(1).Data(3).MatlabIn);
% 
% voxIFG = strcat(DataDir,subID,'/bv/','sub-08_rh_IFG_VOI2Mat.mat');
% load(voxIFG);
% TCIFG_RH = ExtractTS_mat(vtcNonFiltered,ROI(1).Data(1).MatlabIn);
% %plot data 

%%
%band-pass filtering, get scripts from SPM 
%SPM is not using the filtfilt, now I have to write my own high pass
%filtering 
%filterin cuttoff at 0.01 Hz
%still needs to know the width of the filter 
% fs = 1/1.5; %TR is 1.5s, sampling rate is 1/1.5 
% nyquist = fs/2;
% lower_filter_bound = 0.01; % Hz
% upper_filter_bound = 0.01; % Hz
% transition_width   = 0.2;
% filter_order       = round(2*(fs/lower_filter_bound)); %ideally with 3 cycles times the minimal frequency
% %high-pass filtering, don't know what to do with the filter order, it cant
% %be zero, right
% 
% ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
% idealresponse = [ 0 0 1 1 1 0 ];
% HPfilter = firls(filter_order,ffrequencies,idealresponse);
% 
% %test on one voxel 
% filter_testdata = filtfilt(HPfilter,1,TCpSTS_RH(:,1)');
% %the sampling frequency doesn't allow me to use filter with filter order =
% %2
%%
%this is from signal processing tool box, can also try. should filter each
%voxel and then average the signal 
%filter length is usually 100-128 seconds long 
% Fstop = ;
% Fpass = 0.01;
% Astop = 65;
% Apass = 0.5;
% Fs = 1/1.5;
% 
% d = designfilt('highpassiir','StopbandFrequency',Fstop ,...
%   'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
%   'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','butter');
% 
% fvtool(d)
% %there is no nyquist frequency for fMRI (I think, or I don't know, maybe TR*2?) and also based on online
% %discussion. So I am using spm package for high-pass filtering 


%%
%lastly, from SPM tool box 
%spm_filter(K,Y)--> can't find manual for what is the filter k. but from
%Mumford it seems like it's using a gaussian linear model with high
%frequency sine waves as regressors to filter out high frequency (assumed
%to be noise signal) 
%spm cannot find the filter matrix, try to write one myself. 
%--> this seems impossible with limited information 

%Plot power spectrum of data before and after applying high-pass filtering.
%
%
%%
%plotting for high-pass filter results 
% figure;
% subplot(2,1,1);
% time  = linspace(0,size(TCpSTS_RH,1)*1.5,size(TCpSTS_RH,1));
% plot(time,mean(TCpSTS_RH,2));
% xlabel('Time (s)');
% ylabel('BOLD Activity');
% title('Average pSTS RH Time series')
% 
% subplot(2,1,2);
% plot(time,mean(TCpSTS_RH,2));
% 
% figure;
% subplot(2,1,1);
% time  = linspace(0,size(TCIFG_RH,1)*1.5,size(TCIFG_RH,1));
% plot(time,mean(TCIFG_RH,2));
% xlabel('Time (s)');
% ylabel('BOLD Activity');
% title('Average IFG RH Time series')
% 
% subplot(2,1,2);
% plot(time,mean(TCIFG_RH,2));

%ok Not doing high pass filtering. Go straight to high-level connectivity
%measure. 
%%
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


figure;

% time  = linspace(0,size(TCpSTS_RH,1)*1.5,size(TCpSTS_RH,1));
% plot(time,mean(TCpSTS_RH,2));
% hold on;
% plot(time,mean(fil_TCpSTS_RH(:,:,5),2));
% xlabel('Time (s)');
% ylabel('BOLD Activity');
% legend('raw Time series','high pass filtered TS');
% title('Average pSTS RH Time series')
% hold off;
% 
% figure;
% time  = linspace(0,size(TCIFG_RH,1)*1.5,size(TCIFG_RH,1));
% plot(time,mean(TCIFG_RH,2));
% hold on;
% plot(time,mean(fil_TCIFG_RH(:,:,5),2));
% xlabel('Time (s)');
% ylabel('BOLD Activity');
% legend('raw Time series','high pass filtered TS');
% title('Average IFG RH Time series')
% hold off;


%%
%using BSMART to calculate grainger causality index 

%input:
%average pSTS and IFG time series: 
bic = zeros(7,8);
for run = 1:8
    tempdata = [mean(fil_TCpSTS_RH(:,:,run),2)'; mean(fil_TCIFG_RH(:,:,run),2)'];
    
    for orderp=1:7
        % run model
        [~,E] = armorf(tempdata,1,length(tempdata),orderp);
        % compute Bayes Information Criteria
        bic(orderp,run) = log(det(E)) + (log(length(tempdata))*orderp*2^2)/length(tempdata);
    end
end 
%I think I would choose this one because order of 
figure;
plot(1:7,bic,'-o')
legend('run 1','run 2', 'run 3', 'run 4', 'run 5', 'run 6', 'run 7', 'run 8')
xlabel('order number')
ylabel('BIC')

y2x = zeros(1,8);
x2y = zeros(1,8);
for run = 1:8
   tempdata = [mean(fil_TCpSTS_RH(:,:,run),2)'; mean(fil_TCIFG_RH(:,:,run),2)'];
    
    %after picking the order 5
    [~,Ex] =  armorf(tempdata(1,:),1,length(tempdata),5);
    [~,Ey] =  armorf(tempdata(2,:),1,length(tempdata),5);
    [~,E] =  armorf(tempdata,1,length(tempdata),5);

    %causal estimates 
    y2x(run)=log(Ex/E(1,1));
    x2y(run)=log(Ey/E(2,2));
end


%use this 
tempdata = [reshape(mean(fil_TCpSTS_RH,2),[],1)'; reshape(mean(fil_TCIFG_RH,2),[],1)'];
bic = zeros(1,7);
for orderp=1:7
    % run model
    [~,E] = armorf(tempdata,8,size(fil_TCIFG_RH,1),orderp);
    % compute Bayes Information Criteria
    bic(orderp) = log(det(E)) + (log(length(tempdata))*orderp*2^2)/length(tempdata);
end
figure;
plot(1:7,bic,'-o')
xlabel('order number')
ylabel('BIC')
%pick order five from this one 

%after picking the order 5
[~,Ex] =  armorf(tempdata(1,:),8,size(fil_TCIFG_RH,1),5);
[~,Ey] =  armorf(tempdata(2,:),8,size(fil_TCIFG_RH,1),5);
[~,E] =  armorf(tempdata,8,size(fil_TCIFG_RH,1),5);

%causal estimates 
y2x=log(Ex/E(1,1));
x2y=log(Ey/E(2,2));


%%
%permutation testing block 
%shift pSTS region's data from 2 to 201 by 1k time (start) 
%calculate GC index x2y and y2x 
%and then use the distribution as baseline for significance cut off
numi = 10000;
pSTS_RH_all = reshape(mean(fil_TCpSTS_RH,2),[],1)';
IFG_RH_all = reshape(mean(fil_TCIFG_RH,2),[],1)'; 
blocklength = size(fil_TCIFG_RH,1);
OutIFG2pSTS = zeros(1,numi);
OutpSTS2IFG = zeros(1,numi);
for i= 1:numi
    shiftn = randi([2,201],1,1);
    pSTS_in = circshift(pSTS_RH_all,shiftn);

    [~,Ex] =  armorf(pSTS_in,8,blocklength,5);
    [~,Ey] =  armorf(IFG_RH_all,8,blocklength,5);
    [~,E] =  armorf([pSTS_in;IFG_RH_all],8,blocklength,5);

    %causal estimates 
    OutIFG2pSTS(i)=log(Ex/E(1,1));
    OutpSTS2IFG(i)=log(Ey/E(2,2));
end

%sort the results and report the two-tail 5% 
CIIFG2pSTS = sort(OutIFG2pSTS);
CIIFG2pSTS(numi*(1-0.01))
CIpSTS2IFG = sort(OutpSTS2IFG);
CIpSTS2IFG(numi*(1-0.01))


figure;
subplot(2,1,1)
hist(OutIFG2pSTS);
xlabel('GC IFG on pSTS')
subplot(2,1,2)
hist(OutpSTS2IFG);
xlabel('GC pSTS on IFG')
title('Non parametric testing of GC between IFG and pSTS')