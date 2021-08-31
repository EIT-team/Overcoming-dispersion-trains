%% PLOTTING TRAINS, I Tarotin 2021
% Based on the data processed with "post_processing_trains_clean.m"

tr_freq = 10; cuff_n = 3; % Train frquency cuff number. For correct plot titles

LPF = 10; % LPF frequency

% Fs = 10^5; T_window = 6; Fc = 1500; % uncomment if necessary
BV = 0; % Plot BV

ylimit = 1; % ylim
dy = 0.5; % ytick

exchan = 1; % exclude N="exchan" last channels

f_hpf = 1.0; % High pass filter to remove low-frequency noise
[c,d] = butter(3,f_hpf/(Fs/2),'high');
dZ_hpf = cell(1,4);dZ_p1 = cell(1,length(dZ));
dZ_mean_d1 = cell(1,4);
for i = 1 : length(dZ)
    dZ_mean_d1{i} = filtfilt(c,d,dZ_mean_d{i});
    dZ_p1{i} = filtfilt(c,d,dZ_p0{i}); % dZ_p0
end


[c1,d1] = butter(1,LPF/(Fs/2),'low');
t2 = round(1+0*T_window*Fs:1*T_window*Fs);
ylim1 = ylimit; ylim2 = ylimit; % 3 2
xlim1 = -500; xlim2 = 3000;

dztr = cell(1,2);
figure;
subplot(211);
cnt = 1;
for i =  1 : length(dZ_mean_d)-exchan
    dztr{1}(:,cnt) = detrend(filtfilt(c1,d1,dZ_mean_d1{i}(t2)));
    plot(T(t2),dztr{1}(:,cnt),'linewidth',1.2);
    hold on;
    cnt = cnt + 1;
end
ylim(1.*[-ylim1 ylim2]);
xlim([xlim1 xlim2]);
set(gca,'fontsize',7,'ytick',-5:1:5);
title([num2str(Fc) ' Hz, ' num2str(tr_freq) ' Hz train, Cuff ' num2str(cuff_n) ', HPF']);
xlabel('Time (ms)');ylabel('\muV');
legend('location','southeast');set(gca,'fontsize',7,'ytick',-10:dy:10);
subplot(212);
cnt = 1;
for i =  1 : length(dZ_p1)-exchan
    dztr{2}(:,cnt) = (filtfilt(c1,d1,dZ_p1{i}));
    plot(T,dztr{2}(:,cnt),'linewidth',1.2);
    hold on;
    cnt = cnt + 1;
end
ylim(1.*[-0.001*ylim1 0.001*ylim2]);
xlim([xlim1 xlim2]);
title([num2str(Fc) ' Hz, ' num2str(tr_freq) ' Hz train, Cuff ' num2str(cuff_n) ', percent, HPF']);
xlabel('Time (ms)');ylabel('%');
set(gca,'fontsize',7);set(gca,'fontsize',7,'ytick',(-10:dy:10).*1e-3);

