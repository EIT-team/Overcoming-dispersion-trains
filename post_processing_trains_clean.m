%% Post-processing trains, Ilya Tarotin 2021
% Only simple plot is built here. See "plot_traindZ.m" with the improved appearance of the resulting data
clear;

Name = '20ma_50us_cuffs_st56_trains_10Hz_100x.vhdr';
N_spikes = 5; % Number of spikes per train -1
Channels = 5; % cuff 3 - 13:16; cuff4 - 19:21
Fs = 100000; % sampling frequency
T_window = 6; % time window for cutting
dZ_BW = 100; % initial BW. WIll then apply LPF = 10 Hz
Fc = 6000; % AC frequency
N_butter_dZ = 3;
OPEN = 1; % Open the file    
SAVE = 0; % save the file
ScouseTom = 1; % 0 if not using ScouseTom, 1 otherwise
if OPEN
    EEG = pop_loadbv( '',Name,[],[Channels]);
end
Data = double(EEG.data');

% Execute to onfirm the correct frequency
%{
inj_l = 3; % injecting electrode number
V_inj = detrend(Data(:,inj_l),'constant');
NFFT = 2^nextpow2(length(V_inj));           % Next power of 2 from length of y
Y = fft(V_inj,NFFT)/length(V_inj);
f = Fs/2*linspace(0,1,NFFT/2+1);
w_inj=2*abs(Y(1:NFFT/2+1)); % , inj_l
[~,maxw] = max((w_inj(50:end)));
Fc1 = f(maxw); % Fc = ExpSetup.Freq;
fprintf('****** Detected carrier frequency: Fc = %i Hz ******\n',Fc1);
%}

T_trig=cell2mat({EEG.event.latency})';

% Arduino Trig code, I Tarotin 2021
if ScouseTom == 0 % no ScouseTom - triggering Keithley with Arduino
    Tdelta = T_trig(2:end)-T_trig(1:end-1);

    for i = 1 : length(Tdelta)
        if Tdelta(i) <= 205 && Tdelta(i) >= 195
            start = i; % Looking for the start of the train
            break;
        end
    end
    T_trig(1:start-1) = 0;
    for i=start:length(T_trig) % Starting from the first train
        if ~strcmp([EEG.event(1,i).type], ['R  7'])% 'R 15'no! 'R  7' for ArdTrig, 'S  2' for ScT
            T_trig(i)=0;
        end
    end
    T_trig(T_trig==0 | T_trig < 1e5)=[]; 
    T_trig1=zeros(length(T_trig),1);
    for i = 1 : length(T_trig)
        if mod(i,6) == 0
            T_trig1(i) = T_trig(i-5);
        end
    end
    T_trig1(T_trig1==0) = []; % Remove triggers between trains (added by Arduino)
elseif ScouseTom == 1
    % ScouseTom code - if ScouseTom was involved
    T_trig0 = T_trig;
    for i=1:length(T_trig)
        if ~strcmp([EEG.event(1,i).type], ['S  2']); % S2 - stim, S1 - start of train
            T_trig0(i)=0;
        end
    end
    T_trig0(T_trig0==0 | T_trig0 < 1e5)=[]; % Edited by Ilya
    T_trig0=T_trig0(3:end);

    for i=1:length(T_trig)
        if ~strcmp([EEG.event(1,i).type], ['S  1']); % S2 - stim, S1 - start of train
            T_trig(i)=0;
        end
    end
    T_trig(T_trig==0 | T_trig < 1e5)=[]; % Edited by Ilya
    T_trig=T_trig(1:end);

    % clear EEG; % saving RAM

    % Compute number of spikes in each train to be sure they are the same
    cnt1 = 1;
    clear T_trig1 T_trig2 Nbad;
    for i = 1 : length(T_trig)
        cnt = 1;
        for k = cnt1 : length(T_trig0)
            if i ~= length(T_trig)
                if T_trig0(k) < T_trig(i+1)
                    T_trig1(i,cnt) = T_trig0(k);
                    cnt = cnt + 1;cnt1 = cnt1 + 1;
                end
            else
                T_trig1(i,cnt) = T_trig0(k);
                cnt = cnt + 1;
            end
        end
    end

    cnt = 1;
    for i = 1 : size(T_trig1)
        if nnz(T_trig1(i,:)) <= N_spikes-1 % Number of spikes - 1
            Nbad(cnt) = i; % bad rows, after looking at T_trig1
            cnt = cnt + 1;
        end
    end

    if exist('Nbad','var')
        T_trig1(Nbad,:) = [];
    end

    T_trig1(:,[size(T_trig1,2)-1 size(T_trig1,2)]) = []; % first 2 stims
    T_trig2 = reshape(T_trig1',1,numel(T_trig1));
    %  T_trig(Nbad) = [];
end

N_chan = size(Data,2);
N_trig = size(T_trig1,1);
N_bin = round(T_window*Fs);
w = (-round(N_bin/2):round(N_bin/2)); % window around trigger
% w = (1:N_bin);    % window after trigger
T = 1e3*w/Fs;

% Processing as trains
[b,a] = butter(N_butter_dZ,(Fc+dZ_BW*[-1,1])/(Fs/2));
Data_filt = zeros(size(Data,1),N_chan);
Data_hilb = zeros(size(Data,1),N_chan);
Data_hilb_lpf = zeros(size(Data,1),N_chan);
f_lpf = 3; % low-pass
[c,d] = butter(3,f_lpf/(Fs/2),'low');
for i = 1 : N_chan
    Data_filt(:,i) = filtfilt(b,a,Data(:,i));
    Data_hilb(:,i) = abs(hilbert(Data_filt(:,i)));
    Data_hilb_lpf(:,i) = filtfilt(c,d,abs(hilbert(Data_filt(:,i))));
end

clear Data Data_filt; % To save RAM

dZ = cell(1,N_chan); dZ_lpf = cell(1,N_chan); 
for iChan = 1:N_chan
    for jTrig = 1:N_trig-1
        ival = T_trig1(jTrig,1)+w;%-5e5; %+2*st_width*Fs*10^(-3);
        dZ{iChan}(:,jTrig) = Data_hilb(ival,iChan);
        dZ_lpf{iChan}(:,jTrig) = Data_hilb_lpf(ival,iChan);
    end
end

% clear Data_hilb Data_hilb_lpf; % uncomment to save RAM
    
bv = zeros(1,length(dZ));
for iChan = 1:length(dZ)
    bv(iChan) = mean(dZ{iChan}(round(0.2*T_window*Fs:0.8*T_window*Fs)));
end

figure;bar(bv./1000); % boundary voltages - to find out if the data was recorded correctly
title('Boundary voltages - Cuff3');
xlabel('Electrode number');ylabel('BV (mV)')

t0 = round(0.5*T_window*Fs:0.7*T_window*Fs);

dZ_mean = cell(1,length(dZ)); dZ_mean_d = cell(1,length(dZ)); 
dZ_mean_d_lpf = cell(1,length(dZ)); 
for iChan = 1:length(dZ)
    dZ_mean{iChan} = mean(dZ{iChan},2);
    dZ_mean_d{iChan} = detrend(dZ_mean{iChan});
    dZ_mean_d_lpf{iChan} = detrend(mean(dZ_lpf{iChan},2));
end

% Efficient algorithm for percents
t1 = round(1+0.7*T_window*Fs:0.9*T_window*Fs);
dZ_p0 = cell(1,length(dZ));
for iChan = 1:length(dZ)
    dZ_p0{iChan} = 100*(-1+dZ_mean{iChan}/mean(dZ_mean{iChan}(t1)));
end


[c1,d1] = butter(3,50/(Fs/2),'low');
t2 = round(1+0*T_window*Fs:1*T_window*Fs);
figure;
subplot(211);
for i =  1 : length(dZ_mean_d)
    h(i) = plot(T(t2),detrend(filtfilt(c1,d1,dZ_mean_d{i}(t2))),'linewidth',1.2);
    hold on;
end
ylim([-2 2]);set(gca,'ytick',-2:1:2);
xlim([-30 3000]);
% set(gca,'ytick',-50:10:50);
title([num2str(Fc) ' Hz, 0.5 s trains, 20 cm']);
xlabel('Time (ms)');ylabel('\muV');
subplot(212);
for i =  1 : length(dZ_mean_d)
    h1(i) = plot(T,detrend(dZ_p0{i}),'linewidth',1.2);
    hold on;
end
ylim([-0.002 0.002]);set(gca,'ytick',-0.005:0.001:0.005);
xlim([-30 3000]);
title([num2str(Fc) ' Hz, trains, percent']);
xlabel('Time (ms)');ylabel('%');
legend('location','southeast');

% Save to file
if SAVE
    if exist('art_chan','var')
        save([Name(1:end-5) '.mat'],'dZ','dZ_mean','dZ_mean_d','dZ_p0','dZ_badrec','T','Fs','Fc','T_window','bv','-v7.3');
    else
        save([Name(1:end-5) '_BW' num2str(dZ_BW) '.mat'],'dZ','dZ_mean','dZ_mean_d','dZ_p0','T','Fs','Fc','T_window','bv','-v7.3');
    end
end