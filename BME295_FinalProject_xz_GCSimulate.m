%%
%simulation block 
%load design matrix 
addpath('/data/home/xiaojue/bin/bsmart/')
subID = 'sub-08';
designPath = '/data1/2018_ActionDecoding/analysis_fc/Despike/DesignM/NewLSS/';
load(strcat(designPath,'sub-08DespikeLSS_all.mat'));
%loaded outDesignMatrix 
%add up the first 24 as hrf convolved test signal 
Ytest = zeros(1,202*8);
for i=1:8
    Xlss = OutDesignMatrix.Xlss{i};
    Ytest((1+202*(i-1)):(202+202*(i-1))) = sum(Xlss(:,1,:),3);
end
%and concatenate into long by 8 run signal 

%add noise 
Y_high = Ytest + normrnd(0,0.4,1,length(Ytest));
Y_med = Ytest + normrnd(0,0.15,1,length(Ytest));
Y_low = Ytest + normrnd(0,0.05,1,length(Ytest));
%Ytest signal std = 0.2
%plot to check simulated fMRI data 

figure;
subplot(4,1,1)
plot(1:length(Ytest),Ytest);
title('test signal')
ylabel('Simulated BOLD signal');
subplot(4,1,2)
plot(1:length(Ytest),Y_low);
title('low noise')
ylabel('Simulated BOLD signal');
subplot(4,1,3)
plot(1:length(Ytest),Y_med);
ylabel('Simulated BOLD signal');
title('medium noise')
subplot(4,1,4)
plot(1:length(Ytest),Y_high);
title('high noise');
xlabel('Volumes');
ylabel('Simulated BOLD signal');
%find trial index 



%Simulation for GC analysis 

%low noise data 
%shift data from 2-7 data points
orderall = 2:7;
OutBIClow = zeros(length(orderall));
OutGCxylow = zeros(length(orderall));
OutGCyxlow = zeros(length(orderall));

OutBICmed = zeros(length(orderall));
OutGCmed = zeros(length(orderall));
OutGCxymed = zeros(length(orderall));
OutGCyxmed = zeros(length(orderall));

OutBIChigh = zeros(length(orderall));
OutGChigh = zeros(length(orderall));
OutGCxyhigh = zeros(length(orderall));
OutGCyxhigh = zeros(length(orderall));

for ord = 1:length(orderall)
    order = orderall(ord);
    Yxtest_low  = circshift(Ytest,order) + normrnd(0,0.05,1,length(Ytest));
    Yxtest_med  = circshift(Ytest,order) + normrnd(0,0.15,1,length(Ytest));
    Yxtest_high  = circshift(Ytest,order) + normrnd(0,0.4,1,length(Ytest));

    for p = 1:length(orderall)
        timeord = orderall(p);
        [~,Elow] = armorf([Y_low;Yxtest_low],8,length(Xlss),timeord);
        OutBIClow(ord,p) = log(det(Elow)) + (log(length(Y_low))*timeord*2^2)/length(Y_low);
        [~,Ex] =  armorf(Y_low,8,length(Xlss),5);
        [~,Ey] =  armorf(Yxtest_low,8,length(Xlss),5);
        %causal estimates 
        OutGCyxlow(ord,p)=log(Ex/Elow(1,1));
        OutGCxylow(ord,p)=log(Ey/Elow(2,2));
        
        [~,Emed] = armorf([Y_med;Yxtest_med],8,length(Xlss),timeord);
        OutBICmed(ord,p) = log(det(Emed)) + (log(length(Y_med))*timeord*2^2)/length(Y_med);
        [~,Ex] =  armorf(Y_med,8,length(Xlss),5);
        [~,Ey] =  armorf(Yxtest_med,8,length(Xlss),5);
        OutGCyxmed(ord,p)=log(Ex/Emed(1,1));
        OutGCxymed(ord,p)=log(Ey/Emed(2,2));
        
        [~,Ehigh] = armorf([Y_high;Yxtest_high],8,length(Xlss),timeord);
        OutBIChigh(ord,p) = log(det(Ehigh)) + (log(length(Y_high))*timeord*2^2)/length(Y_high);
        [~,Ex] =  armorf(Y_high,8,length(Xlss),5);
        [~,Ey] =  armorf(Yxtest_high,8,length(Xlss),5);
        OutGCyxhigh(ord,p)=log(Ex/Ehigh(1,1));
        OutGCxyhigh(ord,p)=log(Ey/Ehigh(2,2));
    end
end

figure(1);
for i=1:6
    subplot(6,2,1+2*(i-1))
    plot(orderall,OutBIClow(i,:),'-o')
    ylabel(['shift ' num2str(i+1) 'BIC'])
    xlabel('fitted order')
    subplot(6,2,2+2*(i-1))
    plot(orderall,OutGCxylow(i,:),'-o')
    hold on;
    plot(orderall,OutGCyxlow(i,:),'-o')
    hold off;
    legend('Y on shifted Y','shifted Y on Y','Location','northwest')
    xlabel('fitted order')
    ylabel('GC index')
end

figure(2);
for i=1:6
    subplot(6,2,1+2*(i-1))
    plot(orderall,OutBICmed(i,:),'-o')
    ylabel(['shift ' num2str(i+1) 'BIC'])
    xlabel('fitted order')
    subplot(6,2,2+2*(i-1))
    plot(orderall,OutGCxymed(i,:),'-o')
    hold on;
    plot(orderall,OutGCyxmed(i,:),'-o')
    hold off;
    legend('Y on shifted Y','shifted Y on Y','Location','northwest')
    xlabel('fitted order')
    ylabel('GC index')
end

figure(3);
for i=1:6
    subplot(6,2,1+2*(i-1))
    plot(orderall,OutBIChigh(i,:),'-o')
    ylabel(['shift ' num2str(i+1) 'BIC'])
    xlabel('fitted order')
    subplot(6,2,2+2*(i-1))
    plot(orderall,OutGCxyhigh(i,:),'-o')
    hold on;
    plot(orderall,OutGCyxhigh(i,:),'-o')
    hold off;
    legend('Y on shifted Y','shifted Y on Y','Location','northwest')
    xlabel('fitted order')
    ylabel('GC index')
end