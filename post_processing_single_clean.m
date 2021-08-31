%% Post-processing single spikes, Ilya Tarotin 2021
% Only simple plot is built here. See "plot_singledZ.m" with the improved appearance of the resulting data
clear;
freq = [225 425 625];
CS = [10 20];
inj_l = 2;
f_lpf = 300;

name = '20ma_50us_2hz_cuffs_ref13_eit1025_120ua';
Channels = 5:6; % Channels to open, 7:9
% EIT I frequency
Fc = 1025; % can be computed, see below
Fc_randphase = 200; % Fc above this - randphase. below - sumsub
SAVE = 0; % save the file
OPEN = 1;
event = 'S  2';% 'R 15' for ArdTrig, 'S  2' for ScouseTom
T_window = 0.8; % time to take around trigger (cutting window)

% Bandwidth for demodulation
if Fc == 225
    dZ_BW = 100;
elseif Fc == 300
    dZ_BW = 200;
elseif Fc == 400
    dZ_BW = 200;
elseif Fc == 625
    dZ_BW = 400;
elseif Fc == 800
    dZ_BW = 200;
elseif Fc == 1025
    dZ_BW = 500;
elseif Fc > 1025
    dZ_BW = 800;
elseif Fc < 225
    dZ_BW = 50;
else
    disp('No BW entered!');
    return;
end

files=dir([name '.vhdr']);

files={files.name};

for ffil = 1:length(files)
    
    EIT_fname = files{ffil};


    % Parameters for post-processing
    EP_cutoff = 100;  % cutoff frequency for EPs (low-pass freq.)
    N_butter_EP = 3; % Butterworth order for filtering EPs
    N_butter_dZ = 3;  % Butterworth order for filtering dZs
    N_butter_notch = 3; % Butterworth order for notch filtering
    Noise_thres = 1000;  % Noise threshold (in uV)
    Filter_50Hz = 0; % 50 Hz notch filter on EPs
    st_width = 0; % milliseconds

    Fs = 100000;


    if OPEN
        EEG = pop_loadbv( '',EIT_fname,[],[Channels]);
        Data = double(EEG.data');
    end

    %{ 
    % Execute to confirm if the frequency is correct
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

    for i=1:length(T_trig)
        if ~strcmp([EEG.event(1,i).type], [event]); % 'R 15' for ArdTrig, 'S  2' for ScT
            T_trig(i)=0;
        end
    end

    T_trig(T_trig==0 | T_trig < 1e5)=[]; % Edited by Ilya
    T_trig=T_trig(3:end-2);


    if Filter_50Hz
        [b_n,a_n] = iirnotch(50/(Fs/2),(50/(Fs/2))/(35));
        Data = filtfilt(b_n,a_n,Data);
    end


    % Number of channels, triggers, bins
    N_chan = size(Data,2);
    N_trig = length(T_trig);
    N_bin = round(T_window*Fs);
    w = (1:N_bin)-round(N_bin*0.5);    % window around trigger
    T = 1e3*w/Fs;

    % Filtering before cutting. No need to sumsub
    [b,a] = butter(N_butter_dZ,(Fc+dZ_BW*[-1,1])/(Fs/2));
    Data_filt = zeros(size(Data,1),N_chan);
    Data_hilb = zeros(size(Data,1),N_chan);
    Data_hilb_filt = zeros(size(Data,1),N_chan);
    for i = 1 : N_chan
        Data_hilb(:,i) = abs(hilbert(filtfilt(b,a,Data(:,i))));
    end

    clear Data; % save RAM

    dZ = cell(1,N_chan);dZ_lpf = cell(1,N_chan); % EP0_hilb = cell(1,N_chan); 
    EP_lpf = cell(1,N_chan);
    for iChan = 1:N_chan
        dZ{iChan} = zeros(N_bin,N_trig); 
        EP_lpf{iChan} = zeros(N_bin,N_trig); 
    for jTrig = 1:N_trig
        ival = T_trig(jTrig)+w; %+2*st_width*Fs*10^(-3);
        dZ{iChan}(:,jTrig) = Data_hilb(ival,iChan);
    end
    end
    dZ_mean_high = cell(1,N_chan); dZ_mean_d1_filt = cell(1,N_chan); 
    dZ_mean_d_filt = cell(1,N_chan); 
    for iChan = 1:N_chan
        t0 = round(0.2*T_window*Fs:0.3*T_window*Fs);
        dZ_mean_high{iChan} = mean(dZ{iChan},2);
        dZ_mean_d_filt{iChan} = detrend(dZ_mean_high{iChan});
        dZ_mean_d1_filt{iChan} = dZ_mean_high{iChan} - mean(dZ_mean_high{iChan}(t0));
    end

    % Efficient algorithm for percents
    t1 = round(1+0.51*T_window*Fs:0.6*T_window*Fs);
    dZ_p = cell(1,length(dZ));
    for iChan = 1:length(dZ)
        dZ_p{iChan} = 100*(-1+dZ_mean_high{iChan}/mean(dZ_mean_high{iChan}(t1)));
    end

    % Boundary voltages
    bv = zeros(1,N_chan);
    for iChan = 1:N_chan
        bv(iChan) = mean(dZ{iChan}(round(0.2*T_window*Fs:0.8*T_window*Fs)));
    end


    %%% NEW 09032021
    dZ_mean_high = cell(1,N_chan); dZ_mean_d1_filt = cell(1,N_chan); 
    dZ_mean_d_filt = cell(1,N_chan); 
    for iChan = 1:N_chan
        t0 = round(0.2*T_window*Fs:0.3*T_window*Fs);
        dZ_mean_high{iChan} = mean(dZ{iChan},2);
        dZ_mean_d_filt{iChan} = detrend(dZ_mean_high{iChan});
        dZ_mean_d1_filt{iChan} = dZ_mean_high{iChan} - mean(dZ_mean_high{iChan}(t0));
    end

    % Efficient algorithm for percents
    t1 = round(1+0.51*T_window*Fs:0.6*T_window*Fs);
    dZ_p = cell(1,length(dZ));
    for iChan = 1:length(dZ)
        dZ_p{iChan} = 100*(-1+dZ_mean_high{iChan}/mean(dZ_mean_high{iChan}(t1)));
    end


    [b1,a1] = butter(1,500/(Fs/2),'low');
    xl1 = -10; xl2 = 100; 
    figure;
    subplot(211);
    chan = 1 : length(dZ_mean_d1_filt)-0;
    for i =  chan
        h(i) = plot(T,filtfilt(b1,a1,detrend(dZ_mean_d1_filt{i})),'linewidth',1.2);
        hold on;
    end
    ylim([-2 1]);
    xlim([xl1 xl2]);
    xlabel('Time (ms)');ylabel('\muV');
    title([num2str(Fc) ' Hz, absolute, BW = ' num2str(dZ_BW) ' Hz, 10 mins']);
    leg1 = cellstr(num2str(chan', 'ch %d'));
    legend (leg1,'location','southeast');
    subplot(212);
    for i = chan
        h1(i) = plot(T,filtfilt(b1,a1,(dZ_p{i})),'linewidth',1.2);
        hold on;
    end
    ylim([-0.002 0.001]);
    xlim([xl1 xl2]);
    xlabel('Time (ms)');ylabel('%');
    title([num2str(Fc) ' Hz, percentage, 10 mins']);
    leg1 = cellstr(num2str(chan', 'ch %d'));
    legend (leg1,'location','southeast');

    figure;bar(bv) % plot boundary voltages t ocheck id data was correct

    if SAVE
        save([EIT_fname(1:end-5) '_BW' num2str(dZ_BW) '_processed.mat'],'bv','dZ','dZ_mean_d1_filt','dZ_p','T','Fs','Fc','dZ_BW','T_trig','-v7.3');
    end
    
end
