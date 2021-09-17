function [ccgZlagN] = pairwise_reactivation(HSE_raw, HSE_coeff, th_coeff, myPSTH, beltLength, binsize, binshift);

% HSE_raw -> Detection HSE
ccg_max = [];
ccg_maxID = [];
ccg_kurtosis = [];
ccg_raw = {};
ccg_zerolag = [];

[maxPF.val,maxPF.bin] = max(myPSTH,[],2);
PF_peak = maxPF.bin;
% place field distance for treadmill
pf_distance(1,:) = [1:1:beltLength*2];
pf_distance(2,:) = [[1:1:beltLength] [1:1:beltLength]];
pf_dist_win = (beltLength - 1) / 2;

halfLength = floor(beltLength/2);
bin_win = [0:1:halfLength]';
% binsize = 5;
% binshift = 1;
framebin1 = [0:binshift:(bin_win(end)-binsize)]';
framebin2 = [0+binsize:binshift:bin_win(end)]';
binWindow = [framebin1 framebin2];
clear binsize binshift framebin1 framebin2

for tt = 1:length(HSE_raw);
myraw = HSE_raw{tt}(:,:);
react_roi = find(HSE_coeff(:,tt) > th_coeff);
if length(react_roi) == 0;
else
myraw_h = myraw(react_roi,:);

% cross-correlogram
c_cov = {};
for i = 1:length(myraw_h(:,1));
    myref = myraw_h(i,:);
    clear myxcorr
    for ii = 1:length(myraw_h(:,1));
        myagainst = myraw_h(ii,:);
        [c,lags] = xcorr(myref,myagainst);
        myxcorr(ii,:) = c;
    end
    c_cov{i} = myxcorr;
    
end
clear i myxcorr myref myagaist c
zero_lag = find(lags==0);

% centerized position
mypeak_pf = PF_peak(react_roi);
cov_PFdistID = zeros([length(c_cov),length(c_cov)]);
for i = 1:length(c_cov); % num of react roi
    % centerize place window
    if mypeak_pf(i) < pf_dist_win;
        pf_center = pf_distance(1,(205+mypeak_pf(i)));
    elseif mypeak_pf(i) >= pf_dist_win;
        pf_center = pf_distance(1,(mypeak_pf(i)));
    end
    mydist = pf_distance(1,:) - pf_center;
    mywindow = find(mydist>=-102&mydist<=102);
    mydistance(1,:) = mydist(mywindow); 
    mydistance(2:3,:) = pf_distance(:,mywindow);
    
    for ii = 1:length(mydistance(1,:));
        distID = find(mypeak_pf == mydistance(3,ii));
        cov_PFdistID(distID,i) = mydistance(1,ii);
    end
    clear mydistance mywindow mydist
end
clear i ii distID mydist mywindow mydistance pf_center

% all psth
C_cov_distbase = {};
for t = 1:length(binWindow);
    mydata_all = [];
    mydata_all2 = [];
for i = 1:length(c_cov);
    mydata = c_cov{i}(:,:);
    mydata_ref = mydata(find(cov_PFdistID(:,i)>=binWindow(t,1)&cov_PFdistID(:,i)<=binWindow(t,2)),:);
    mydata_all = [mydata_all; mydata_ref];
end
C_cov_distbase{t,1} = mydata_all;
C_cov_distbase{t,2} = binWindow(t,:);
C_cov_distbase{t,3} = mydata_all2;
C_cov_distbase{t,4} = -1*binWindow(t,:);
end
clear mydata_all mydata mydata_ref t

mypsth_cov = [];
for i = 1:length(C_cov_distbase);
    mydata = C_cov_distbase{i,1};
    mydata_mn = mean(mydata,1);    
    mypsth_cov = [mypsth_cov; mydata_mn];
end
clear i mydata mydata_nm mydata_mn

%%%%%%%%%%%%%%%%%%%%%%%%%
[mymaxVal,mymaxID] = max(mypsth_cov,[],2);
ccg_max(:,tt) = mymaxVal;
ccg_maxID(:,tt) = mymaxID;
ccg_zerolag(:,tt) = mypsth_cov(:,zero_lag);
ccg_kurtosis(:,tt) = kurtosis(mypsth_cov,1,2);
ccg_raw{tt} = mypsth_cov;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
tt
end
clear tt mymaxVal mymaxID

noise_ref = sum(ccg_kurtosis,1);
noise_ref(noise_ref>0) = 1; 
noise_ref(noise_ref==0) = 0; 

myccgZlag = ccg_zerolag(2:end,find(noise_ref==1));
ccgZlagN = (myccgZlag - mean(myccgZlag,1)) ./ std(myccgZlag,0,1);
end


