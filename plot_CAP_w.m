%% Plotting Compound action potentials, I Tarotin 2021
% Based on the raw data (Porcine subdiaphragmatic ex-vivo) 

folder = 'G:\16062021\';
name = '20mA_50us_st56_8.eeg';
filename = [folder name];
EEG_cap=pop_fileio(filename);

Fs = EEG_cap.srate;

t=cell2mat({EEG_cap.event.latency})';

jj=zeros(length(t),1);
for i=1:length(t)
    if strcmp([EEG_cap.event(1,i).type], ['R 15']); % 'S  2' for scouseTom, 'R 15' for ard
        jj(i)=1;
    end
end

T_trig=t(jj==1);
T_trig = T_trig(3:end-1);

tau=500;
size_bin=floor(tau*Fs/1000);

Data_ap= double(EEG_cap.data)';

[b,a] = butter(3,3000/(Fs/2),'low');
Data_ap = filtfilt(b,a,Data_ap);


T1=[1:size_bin]*1000/Fs-tau/2;

EP=zeros(length(T_trig)-2,size_bin,size(Data_ap,2));
for i=2:length(T_trig)-1
    EP(i-1,:,:)=Data_ap([T_trig(i)-floor(size_bin/2):T_trig(i)+floor(size_bin/2)-1],:);
end

EP_avg=detrend(squeeze(mean(EP,1)));
EP_avg = EP_avg(:,:) - mean(EP_avg(1000:floor(size(EP_avg,1)/2.1),:),1);

% Figure
yl = 100; xl = 200;%60+160;
figure
chan = 7:size(EP_avg,2);
plot(T1,EP_avg(:,chan),'linewidth',1.2);
leg1 = cellstr(num2str(chan', 'ch %d'));
xlim([0,xl]);ylim([-yl,yl]);
ylabel('\muV');xlabel('Time (ms)')
% title(['Pig nerve, ' num2str(name(1:4)) ', ' num2str(name(6:7)) '\mus, absolute CAPs']);
title('Porcine SN CAPs');
set(gca,'fontsize',8,'ytick',-200:20:200);


