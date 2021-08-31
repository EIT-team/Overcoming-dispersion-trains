%% Plotting single dZ, I Tarotin 2021
% Based on the data processed with "post_processing_single_clean.m"

[b1,a1] = butter(3,200/(Fs/2)); % low pass filter
xl1 = 0; xl2 = 200; 
exchan = 0; % exclude N channels

figure;
subplot(211); % uV
chan = 2 : length(dZ_mean_d1_filt)-exchan;
for i =  chan
    h(i) = plot(T,smooth(filtfilt(b1,a1,detrend(dZ_mean_d1_filt{i})),201)+0.1,'linewidth',1.2); % smooth span = 401
    hold on;
end
xlim([xl1 xl2]);
xlabel('Time (ms)');ylabel('\muV');
% title([num2str(Fc) ' Hz, absolute, BW = ' num2str(dZ_BW) ' Hz']);
title('Porcine SN dZ');
leg1 = cellstr(num2str(chan', 'ch %d'));
legend (leg1,'location','southeast');
ylim([-1 0.4]);
set(gca,'fontsize',8,'xtick',0:50:200,'ytick',-1:0.2:1); % 'ytick',-4:0.5:4,

subplot(212); % percent
for i = chan
    h1(i) = plot(T,(filtfilt(b1,a1,(dZ_p{i}))),'linewidth',1.2);
    hold on;
end
xlim([xl1 xl2]);
xlabel('Time (ms)');ylabel('%');
title(['dZ, ' num2str(Fc/1000) ' kHz, Cuff 2 (3 cm)']); % percentage
leg1 = cellstr(num2str(chan', 'ch %d'));
legend (leg1,'location','southeast');
ylim([-5 5].*3e-4); 
set(gca,'xtick',0:20:200,'fontsize',10); % 'ytick',(-15:5:5).*1e-4,

%% Plotting boundary voltages
figure;bar(bv./1000);ylim([0 300]);
title('Boundary voltages');
xlabel('Electrode number');ylabel('BV (mV)')