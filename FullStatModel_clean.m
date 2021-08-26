%% Repetitive stimulation - clean version
% // Computation of SNR. MAIN code. 
% Ilya Tarotin, 2021

clear;

dt = 0.1; 
Fs = 1000/dt;
modN = 50; % Number of models (50). Decrease for speed.

% For detailed explanation and derivation of the used coefficients, see the
% manuscript: "Method for overcoming temporal dispersion in unmyelinated
% nerves to image them with Electrical Impedance Tomography (EIT)"
% For data on coefdxi and coefdxr, see IEEE TBME Tarotin et al., 2019
% Note that coef_el, coefdxr and coefdxi are reciprocals to the ones used
% in the manuscript
coef_el = ((1-1/sqrt(2))*1600/10)^-4; % weighted correctly
Icoef = 43.5/4; % 4mA/cm^2 -> 43.5mA/cm^2, dZ+BV increase 10 times
dximod = 0.02; dxiexp = 3.3; % 3.3mm experiment /0.02mm model, 
dz_mod = 0.99; % dz at 0.02 mm, ref tbme 2019
dz_exp = 0.6; % dz at 0.33 mm to find dZ at 3.3 mm, ref tbme 2019
dzdxi = 0.9/0.2; % dz change (slope) in the model, dz1/dz2
coefdxi = dz_exp/(dzdxi*dz_mod); 
dxrmod = 0.11; dxrexp = 1.1; 
dzdxr_uv = 0.01/1; % dz change in percent in the model, equals to uV 
dxr0 = 0.1/0.01; % dxr in the model, tbme
coefdxr = (dzdxr_uv/dxr0)*(dxrexp/dxrmod); % dz falls 10 times faster than dxR, ref tbme
coefVdZ = (coef_el*coefdxr*coefdxi)*(Icoef*coef_cond^2); 

% BV = dxi*coef_cond*0.8.*(Icoef*((1-1/sqrt(2))*1600/10)^-3); % boundary voltage
dC = 1;
vC = 0.8; devC = 0.3; latC = 2; latAd = 1;
numC = 40000;  % Number of fibres, was 5000
vC0_arr = zeros(numC,modN); % Number of models (50) DECREASE FOR SPEED
for i = 1 : size(vC0_arr,2)
    vC0 = []; % Distribution of velocities
    vC0 = [vC0;devC*randn(5*numC,1)+vC];
    vC0(vC0<0.1)=[];
    if length(vC0)> numC % to have exact number of fibres
        vC0 = vC0(1:numC);
    end
    vC0_arr(:,i) = vC0; 
end

% Myelinated fibres (B(A-delta))
dAd = 4; % Comment this variable only to exclde Ad fibres
vAd = 8; devAd = 3;  % Based on CAPs from ex-vivo pig subdiaphragmatic
numAd = 4000;
vAd0_arr = zeros(numAd,modN); % Number of models (50) DECREASE FOR SPEED
for i = 1 : size(vAd0_arr,2)
    vAd0 = [];
    vAd0 = [vAd0;devAd*randn(5*numAd,1)+vAd];
    vAd0(vAd0<0.1)=[];
    if length(vAd0)> numAd % to have exact number of fibres
        vAd0 = vAd0(1:numAd);
    end
    vAd0_arr(:,i) = vAd0; 
end

dist = [30 150 200 500];
cnt = 1;
tC = cell(1,size(vC0_arr,2)); 
tAd = cell(1,size(vAd0_arr,2));
for k = 1 : size(vC0_arr,2)
    for i = dist
        tC{k}(:,cnt) = dist(cnt)./vC0_arr(:,k);
        if exist('dAd','var') % if myelinated fibres are added. TO leave just C fibres - remove dAd
            tAd{k}(:,cnt) = dist(cnt)./vAd0_arr(:,k);
        end
        cnt = cnt + 1;
    end
    cnt = 1;
end

% Taking only 6 peaks at 10 and 20 Hz
y = flip([0.125 0.3 0.6 1.2 10 10]); % 50-20-10-5Hz - 6 pulses, 5 Hz WAS 10 pulses
fr = [1 2 5 10 20 50]; % frequencies


t0 = 1 : 100001;
dz1 = cell(length(fr),1); 
load('dZ_trains_FEM.mat'); % dZ trains obtained with the FEM model (provided on GitHub)
for i = 1 : 3
    dz1{i} = C_RepStim_freqsw_521Hz(t0+100001*(4-i-1),7); %
end
dz1{4} = C_RepStim_freqsw(26003:end,7); %
dz1{5} = C_RepStim_freqsw(6002:26002,7); %
dz1{6} = C_RepStim_freqsw(1:6001,7); %


% CHANGING LATENCY OF DZ WITHOUT CHANGING FREQUENCY - if latC != latAd
% CUTTING -> SCALING -> ADDING ZEROS -> PUTTING BACK TO ARRAY
if exist('dAd','var') % if myelinated fibres are added. TO leave just C fibres - remove dAd
    Np = [10, 20, 40, 24, 14, 10]; % Computed manually ([1 5 10 20 50]Hz)
    % cut windows of signals
    w = 1./fr;
    w = w + 0.0006; % correct for stim.width + AP initiation time
    start = 0.01;
    dz0 = cell(length(fr),1);
    for i = 1 : length(fr)
        dz0{i} = [detrend(dz1{i}); zeros(1000,1)];
    end
    dz00 = cell(length(fr),1);
    for k = 1 : length(fr)
    for i = 1 : Np(k)
        dz00{k}(:,i) = dz0{k}(1+round(1e4*(start+(i-1)*w(k))):round(1e4*(start+i*w(k))));
    end
    end
    dz001 = cell(length(fr),1);
    for k = 1 : length(fr)
        dz001{k} = zeros(length(dz00{k}),Np(k));
        for i = 1 : Np(k)
            buf = imresize(((dAd/dC)^2).*dz00{k}(:,i),latAd/latC);
            dz001{k}(1+1e4*0.0015:1e4*0.0015+length(buf),i) = buf; % +1e4*0.0015 - for temporal adjustment
        end
    end

    % Timnings of new dZ are fine for first 6 pulses - so we don't consider ADS
    % (activity dependent slowing) in this case
    dz_back{k} = cell(length(fr),1);
    for k = 1 : length(fr)
        dz_back{k} = reshape(dz001{k},[1,length(dz001{k})*Np(k)]);
    end
    for i = 1 : length(fr)
        dz_back{i} = [zeros(1,1e4*start) dz_back{i}];
    end
end
% figure;plot(((1/dAd)^2).*dz_back{6});hold on;plot(detrend(dz1{6}))
if exist('dAd','var')
    dz2 = cell(length(y),1);
    for i = 1 : length(y)
        if length(dz_back{i}) > length(dz1{i})
            dz2{i} = dz_back{i}(1:length(dz1{i}))';
        else
            dz2{i} = [dz_back{i} zeros(1,length(dz1{i})-length(dz_back{i}))];
            dz2{i} = dz2{i}';
        end
    end
end


%{
% Use when latC == latAd, instead of the above
if exist('dAd','var')
    dz2 = cell(length(y),1);
    for i = 1 : length(y)
        dz2{i} = imresize(((dAd/dC)^2).*dz1{i},latAd/latC); % Correct for dC
    end
end
%}


dz_cut = cell(length(dz1),2); 
for i = 1 : length(dz1)
    if exist('dAd','var')
        for k = 1 : 2
        if k == 1
            dz_cut{i,k} = dz1{i}(1:y(i)*1000/dt)';
            dz_cut{i,k} = detrend(dz_cut{i,k}); % DETRENDED (offset), -0.793
        else
            dz_cut{i,k} = dz2{i}(1:y(i)*1000/dt)';
            dz_cut{i,k} = detrend(dz_cut{i,k}); % DETRENDED (offset), -0.793
        end
        end
    else
        dz_cut{i,1} = dz1{i}(1:y(i)*1000/dt)';
        dz_cut{i,1} = detrend(dz_cut{i,1}); % DETRENDED (offset), -0.793
    end
end

% Figure - not diapersed dZ
figure;
for i = 1 : 6
    subplot(1,6,i);
    plot((1:length(dz_cut{i,1}))./Fs,dz_cut{i,1});
    ylim([-0.04 0.04]);
end

% Duration of simulation
t_end = max(max(cell2mat(tC)))+max(y)*1000/dt/10; % max(tC{1}(:,end))+5*latC/dt;
t_sim = 0:dt:t_end;

% MAIN COMPUTATION - takes several hours!
sigC = cell(length(dz_cut),size(vC0_arr,2));
for f = 1 : length(dz_cut)
for l = 1 : size(vC0_arr,2)
    sigC{f,l} = zeros(length(t_sim),length(dist)); % {f}
    for k = 1 : length(dist)
        % C fibres
        for i = 1 : size(tC{l},1)
            sigC{f,l}(round(tC{l}(i,k)/dt):round(tC{l}(i,k)/dt)+length(dz_cut{f,1})-1,k) = sigC{f,l}(round(tC{l}(i,k)/dt):round(tC{l}(i,k)/dt)+length(dz_cut{f,1})-1,k) + dz_cut{f,1}';
        end
        if exist('dAd','var')
            % B(Adelta)
            for i = 1 : size(tAd{l},1)
                sigC{f,l}(round(tAd{l}(i,k)/dt):round(tAd{l}(i,k)/dt)+length(dz_cut{f,2})-1,k) = sigC{f,l}(round(tAd{l}(i,k)/dt):round(tAd{l}(i,k)/dt)+length(dz_cut{f,2})-1,k) + dz_cut{f,2}';
            end
        end
    end
end
end
% save('SigC_TRAINS_C_09_03_1um_B_8_3_4um_23072021.mat','sigC','-v7.3');
%{
% If no B fibres added (fully unmyelinated nerve)
sigC = cell(length(dz_cut),size(vC0_arr,2));
% t_sim = cell(1,length(dz_cut));
for f = 1 : length(dz_cut)
for l = 1 : size(vC0_arr,2)
    sigC{f,l} = zeros(length(t_sim),length(dist)); % {f}
    for k = 1 : length(dist)
        % C fibres
        for i = 1 : size(tC{l},1)
            sigC{f,l}(round(tC{l}(i,k)/dt):round(tC{l}(i,k)/dt)+length(dz_cut{f,1})-1,k) = sigC{f,l}(round(tC{l}(i,k)/dt):round(tC{l}(i,k)/dt)+length(dz_cut{f,1})-1,k) + dz_cut{f,1}';
        end
    end
end
end
save('SigC_TRAINS_C_09_03_1um_NO_B_23072021.mat','sigC','-v7.3');
%}


% Figure - dispersed dZ
cnt = 1;
figure;
for i = 1 : length(fr) % 3:4
for d = 1 : 4 % 1:3
subplot(length(fr),4,cnt); % 2,3,cnt
plot((1:length(sigC{i,1}(:,d)))./1e4,1e3.*sigC{i,1}(:,d).*coef); % Scaled
xlim([0 15]); %ylim([-0.3 0.1]);%xlim([0 length(sigC{i,1}(:,d))]./1e4);
title([num2str(fr(i)) ' Hz, ' num2str(dist(d)./10) ' cm']);
if cnt > 20
    xlabel('Time (s)');
end
if mod(cnt,4) == 1
    ylabel('\muV');
end
cnt = cnt + 1;
box off;
end
end

%% SINGLE PULSES - COHERENT SPIKE AVERAGING
% Np = [10, 40, 24, 14, 10]; % Computed manually ([1 5 10 20 50]Hz)
% load('SigC_TRAINS_C_08_03_1um_NO_B_22072021.mat'); % Saved computed database with no B fibres (for simplicity) 
Np = [10, 20, 6, 6, 6, 6]; % amended with 6 pulses/train for fr >= 5Hz
% cut windows of signals
w = 1./fr; % duration of cut segment, s 
w = w + 0.0005; % correct for stim.width + AP initiation time

% Conduction velocity correction - shift in cutting
if exist('dAd','var')
    shift = dist/(1000*max(max(vAd0_arr))); % NEED TO ACCOUNT FOR PROPAGATION VELOCITY, [s]
else
    shift = dist/(1000*max(max(vC0_arr)));
end

% Cutting windows of signals at various distances (dist) corrected by
% conduction velocity (shift)
sigC_single = cell(length(dz_cut),size(vC0_arr,2),length(Np)); 
start = 0.01;
for k = 1 : length(Np) % frequencies
for m = 1 : size(sigC,2) % models (stats)
for d = 1 : size(sigC{1,1},2) % distances
    for i = 1 : Np(k) % peaks, starting from the 2nd
        sigC_single{k,m,i}(:,d) = sigC{k,m}(1+round(Fs*(start+shift(d)+(i-1)*w(k))):round(Fs*(start+shift(d)+i*w(k))),d).*coef;
    end
end
end
end

% Averaging across peaks in a train; 
% If high stim freq - ADS OCCURS (activity dependent slowing)
sigC_single_sum = cell(length(dz_cut),size(vC0_arr,2));
for k = 1 : size(sigC,1) % frequencies
for m = 1 : size(sigC,2) % models (stats)
for d = 1 : size(sigC{1,1},2) % distances
    for i = 1 : Np(k) % peaks, starting from the 2nd
        if i == 1
            sigC_single_sum{k,m}(:,d) = sigC_single{k,m,i}(:,d)./(Np(k)); % pure signal - exclude 1st peak due to large LF stuff (not dZ?) in the model
        else
            sigC_single_sum{k,m}(:,d) = sigC_single_sum{k,m}(:,d) + sigC_single{k,m,i}(:,d)./(Np(k));
        end
    end
end
end
end


% All BW = 200 Hz (lowest possible). 
sigC_single_bw_opt_av = 200.*ones(length(fr),1); % 200 Hz everywhere, or 100?

% ADDING NOISE
% Number of trains
Dtr = y; % train duration, s
int = 5; % 5 s interval
Trec = 30*60; % 30 minutes
Ntr = zeros(length(int),length(Dtr)); % Number of trains per Trec (100-300 s)
for i = 1 : length(int)
    Ntr(i,:) = floor(Trec./(int(i) + Dtr)); 
end
Ntr(:,1:2) = Trec/10; % No intervals in 1 and 2 Hz trains
An0 = 3.5e-3; % 3.5 uV before averaging
for i = 1 : length(Np)
    An(i) = An0./(sqrt(Np(i))*sqrt(Ntr(i))); % noise after averaging
end
sigC_single_sum_noise = cell(length(dz_cut),size(vC0_arr,2));
for k = 1 : size(sigC,1) % frequencies
for m = 1 : size(sigC,2) % models (stats)
for d = 1 : size(sigC{1,1},2) % distances
    sigC_single_sum_noise{k,m}(:,d) = sigC_single_sum{k,m}(:,d) + An(k).*randn(length(sigC_single_sum{k,m}(:,d)),1);
end
end
end

% Filtering the resultant signal. Zeros added.
sigC_single_filt = cell(length(dz_cut),size(vC0_arr,2));
for k = 1 : size(sigC,1) % frequencies
    [b, a] = butter(3, sigC_single_bw_opt_av(k)/(Fs/2)); % Fs = 10e3, 5th order butter
for m = 1 : size(sigC,2) % models (stats)
for d = 1 : size(sigC{1,1},2) % distances
    buf = filtfilt(b,a,[zeros(1e4,1); sigC_single_sum_noise{k,m}(:,d); zeros(1e4,1)]); 
    sigC_single_filt{k,m}(:,d) = buf(1e4+1:length(sigC_single_sum_noise{k,m}(:,d))+1e4); 
end
end
end

%{
% For debugging. Comparison with the filtered signal. Noise was added
cnt = 1;
figure;
for d = 1 : 4
for k = 1 : length(fr)
subplot(4,length(fr),cnt);
plot(sigC_single_filt{k,1}(:,d));
hold on;
plot(sigC_single_sum{k,1}(:,d));
xlim([0 length(sigC_single_sum{k,1}(:,d))]);
cnt = cnt + 1;
end
end
%}

% Averaging + STD across N models
sigC_single_forav = cell(length(dz_cut),size(sigC{1,1},2));
for k = 1 : size(sigC,1) % frequencies
for d = 1 : size(sigC{1,1},2) % distances
    for m = 1 : size(sigC,2) % models (stats)
        sigC_single_forav{k,d}(m,:) = sigC_single_filt{k,m}(:,d);
    end
end
end
% Av + SD
sigC_single_av = cell(length(dz_cut),size(sigC{1,1},2));
sigC_single_SD = cell(length(dz_cut),size(sigC{1,1},2));
for k = 1 : size(sigC,1) % frequencies
for d = 1 : size(sigC{1,1},2) % distances
    sigC_single_av{k,d} = mean(sigC_single_forav{k,d});
    sigC_single_SD{k,d} = std(sigC_single_forav{k,d});
end
end


% Find minimum + indices
sigC_single_min = zeros(length(dz_cut),size(sigC{1,1},2));
sigC_single_minind = zeros(length(dz_cut),size(sigC{1,1},2));
for k = 1 : size(sigC,1) % frequencies
for d = 1 : size(sigC{1,1},2) % distances
    [sigC_single_min(k,d), sigC_single_minind(k,d)] = max(abs(sigC_single_av{k,d})); % or min???
end
end

% Take SD in the index of minimum
sigC_single_minstd = zeros(length(dz_cut),size(sigC{1,1},2));
for k = 1 : size(sigC,1) % frequencies
for d = 1 : size(sigC{1,1},2) % distances
    sigC_single_minstd(k,d) = sigC_single_SD{k,d}(sigC_single_minind(k,d));
end
end

% SNR - coefficient (add amplitude of noise)
% SNR - NEW. HOW LARGE IS NOISE AFTER FILTER
sigC_single_noise0 = cell(length(dz_cut),1);
sigC_single_noise1 = cell(length(dz_cut),1);
for k = 1 : size(sigC,1) % frequencies
    sigC_single_noise0{k} = zeros(2*length(sigC{k,1}(:,1)),200);
    sigC_single_noise1{k} = zeros(2*length(sigC{k,1}(:,1)),200);
for m = 1 : 200 % for stats - need to average noise for consistency
    [b, a] = butter(3, sigC_single_bw_opt_av(k)/(Fs/2)); % Fs = 10e3, 5th order butter
    sigC_single_noise0{k}(:,m) = An(k).*randn(2*length(sigC{k,1}(:,1)),1); % Original noise
    sigC_single_noise1{k}(:,m) = filtfilt(b,a,sigC_single_noise0{k}(:,m));
end
end
sigC_single_noisestd = cell(length(dz_cut),1);
for k = 1 : size(sigC,1) % frequencies
    for m = 1 : 200
        sigC_single_noisestd{k}(:,m) = std(sigC_single_noise1{k}(round(0.2*end):round(0.8*end),m));
    end
end
clear sigC_single_noise0 sigC_single_noise1

% SNR
sigC_single_snr = cell(1,length(int));
for t = 1 : length(int)
    sigC_single_snr{t} = zeros(length(dz_cut),size(sigC{1,1},2));
    for k = 1 : size(sigC,1) % frequencies
    for d = 1 : size(sigC{1,1},2) % distances
        if sigC_single_min(k,d) == 0
            sigC_single_snr{t}(k,d) = 0;
        else
            sigC_single_snr{t}(k,d) = abs(sigC_single_min(k,d))./mean(sigC_single_noisestd{k}); %std(sigC_single_noise1{k});% sigC_single_SD_noise{t}(k,d);...(:,d)
        end
    end
    end
end
for k = 1 : length(fr)
    for d = 1 : size(sigC{1,1},2)
        std_single(k,d) = sigC_single_minstd(k,d)./mean(sigC_single_noisestd{k}); %std(sigC_single_noise1{k}); %(:,d)
    end
end

% FIGURE - SINGLE SPIKES (COHERENT SPIKES AVERAGING)
figure;
for t = 1 : length(int)
for k = 1 : length(fr)
    errorbar(dist./10+0.1*k,sigC_single_snr{t}(k,:),std_single(k,:));
    hold on;
end
end
legend(cellstr(num2str(fr')));
xlabel('Distance (cm)');ylabel('SNR');
title('Coherent spike averaging: C fibres only');
%% TRAINS OF PULSES
% load('SigC_TRAINS_C_08_03_1um_B_8_3_4um_22072021.mat'); % Load sigC with B fibres (for accuracy and agreement with experiment)
% Dtr = y; % train duration, s
% int = 5; % already defined in single pulses
sigC_train_bw_opt = 2./Dtr'; % Hz
% Determine minimal BW
bw_min = 10; % Minimal BW = 10 Hz. See paper for explanation
for i = 1 : length(sigC_train_bw_opt)
    if sigC_train_bw_opt(i) < bw_min
        sigC_train_bw_opt(i) = bw_min;
    end
end

An0 = 3.5e-3; % 3.5 uV before averaging
An = An0./sqrt(Ntr); % after averaging

% Adding noise
sigC_train_noise = cell(length(dz_cut),size(vC0_arr,2));
for k = 1 : size(sigC,1) % frequencies
for m = 1 : size(sigC,2) % models (stats)
for d = 1 : size(sigC{1,1},2) % distances
    sigC_train_noise{k,m}(:,d) = sigC{k,m}(:,d).*coef + An(k).*randn(length(sigC{k,m}(:,d)),1);
end
end
end

% Filtering
sigC_train_filt = cell(length(dz_cut),size(vC0_arr,2));
for k = 1 : size(sigC,1) % frequencies
    [b1, a1] = butter(3, sigC_train_bw_opt(k)/(Fs/2)); % Fs = 10e3, 5th order butter
for m = 1 : size(sigC,2) % models (stats)
for d = 1 : size(sigC{1,1},2) % distances
    buf = filtfilt(b1,a1,[zeros(1e4,1); sigC_train_noise{k,m}(:,d); zeros(1e4,1)]); 
    sigC_train_filt{k,m}(:,d) = buf(1e4+1:length(sigC_train_noise{1,1}(:,1))+1e4); 
end
end
end

Nav = 10; % for representativeness of randon noise
% Noise after filtering
sigC_train_noise0 = cell(length(dz_cut),1);
sigC_train_noise1 = cell(length(dz_cut),1);
for k = 1 : size(sigC,1) % frequencies
    sigC_train_noise0{k} = zeros(2*length(sigC{k,1}(:,1)),Nav);
    sigC_train_noise1{k} = zeros(2*length(sigC{k,1}(:,1)),Nav);
for m = 1 : Nav % for stats - need to average noise for consistency
    [b1, a1] = butter(3, sigC_train_bw_opt(k)/(Fs/2)); % Fs = 10e3, 5th order butter
    sigC_train_noise0{k}(:,m) = An(k).*randn(2*length(sigC{k,1}(:,1)),1); % Original noise
    sigC_train_noise1{k}(:,m) = filtfilt(b1,a1,sigC_train_noise0{k}(:,m));
end
end
sigC_train_noisestd = cell(length(dz_cut),1);
for k = 1 : size(sigC,1) % frequencies
    for m = 1 : Nav
        sigC_train_noisestd{k}(:,m) = std(sigC_train_noise1{k}(round(0.2*end):round(0.8*end),m));
    end
end
clear sigC_train_noise0 sigC_train_noise1 % freeing RAM

% Averaging
sigC_train_forav = cell(length(dz_cut),size(sigC{1,1},2));
sigC_train_forav_0 = cell(length(dz_cut),size(sigC{1,1},2));
for k = 1 : size(sigC,1) % frequencies
for d = 1 : size(sigC{1,1},2) % distances
    for m = 1 : size(sigC,2) % models (stats)
        sigC_train_forav{k,d}(m,:) = sigC_train_filt{k,m}(:,d); 
        sigC_train_forav_0{k,d}(m,:) = sigC{k,m}(:,d).*coef;
    end
end
end
% Av + SD 
sigC_train_av = cell(length(dz_cut),size(sigC{1,1},2));
sigC_train_av_0 = cell(length(dz_cut),size(sigC{1,1},2));
sigC_train_SD = cell(length(dz_cut),size(sigC{1,1},2));
for k = 1 : size(sigC,1) % frequencies
for d = 1 : size(sigC{1,1},2) % distances
    sigC_train_av{k,d} = mean(sigC_train_forav{k,d});
    sigC_train_av_0{k,d} = mean(sigC_train_forav_0{k,d});
    sigC_train_SD{k,d} = std(sigC_train_forav{k,d});
end
end


% Signal before and after filtering
cnt = 1;
figure;
for d = 1 : 4
for k = 1 : 6
subplot(4,6,cnt);
plot(sigC{k,1}(:,d)./coef);
% hold on;
% plot(sigC_train_filt{k,1}(:,d));
hold on;
% plot(sigC_train_filt{k,1}(:,d));
plot(sigC_train_av{k,d});
xlim([0 length(sigC_train_noise{k,1}(:,d))]);
cnt = cnt + 1;
end
end


%% Find amplitude of dZ
nn = -0.7; % if stim artefact -> nn=-0.5, if no stim artefact nn = 3
sigC_train_min = zeros(length(dz_cut),size(sigC{1,1},2));
sigC_train_minind = zeros(length(dz_cut),size(sigC{1,1},2));
for k = 1 : size(sigC,1) % frequencies
for d = 1 : size(sigC{1,1},2) % distances
    [sigC_train_min(k,d), sigC_train_minind(k,d)] = (min(sigC_train_av{k,d}(round(length(dz_cut{k,1})-nn*Fs/fr(k)+Fs*shift(d)):round(length(dz_cut{k,1})-nn*Fs/fr(k)+Fs*shift(d)+0/fr(k)+10000)))); % look only at the expected time point 
end
end
% Take SD in the index of minimum
sigC_train_minstd = zeros(length(dz_cut),size(sigC{1,1},2));
for k = 1 : size(sigC,1) % frequencies
for d = 1 : size(sigC{1,1},2) % distances
    sigC_train_minstd(k,d) = sigC_train_SD{k,d}(sigC_train_minind(k,d));
end
end

% SNR
sigC_train_snr = cell(1,length(int));
for t = 1 : length(int)
    sigC_train_snr{t} = zeros(length(dz_cut),size(sigC{1,1},2));
    for k = 1 : size(sigC,1) % frequencies
    for d = 1 : size(sigC{1,1},2) % distances
        if sigC_train_min(k,d) == 0
            sigC_train_snr{t}(k,d) = 0;
        else
            sigC_train_snr{t}(k,d) = abs(sigC_train_min(k,d))./mean(sigC_train_noisestd{k}); % std(sigC_train_noise1{k}(round(0.3*end):round(0.8*end))); % round(0.2*end):round(0.8*end)
        end
    end
    end
end

for k = 1 : length(fr)
    for d = 1 : size(sigC{1,1},2)
        std_train(k,d) = sigC_train_minstd(k,d)./mean(sigC_train_noisestd{k}); %std(sigC_train_noise1{k}(round(0.3*end):round(0.8*end))); % round(0.2*end):round(0.8*end)
    end
end

% Figure - trains
figure;
for t = 1 : 1
for k = 1 : length(fr)
    errorbar(dist./10+0.1*k,sigC_train_snr{t}(k,:),std_train(k,:));% sigC_train_minstd(k,:)./std(sigC_train_noise1{k}));
    hold on;
end
end
ylim([0 10]);
xlabel('Distance (cm)');ylabel('SNR');
title('Trains: B+C fibres');
%% COMBINED FIGURE
figure;
for t = 1 : 1
for k = 1 : length(fr)
    p(k) = errorbar(dist./10+0.1*k,sigC_single_snr{t}(k,:),std_single(k,:),'marker','.','markersize',8,'linewidth',0.5);
    hold on;
end
Color = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940, 0.1840, 0.5560;0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];
for k = 1 : length(fr)
    p1(k) = errorbar((dist)./10+0.1*k,sigC_train_snr{t}(k,:),std_train(k,:),'--','marker','.','markersize',8,'Color',Color(k,:),'linewidth',0.5);
    hold on;
end
end
xlabel('Distance (cm)');ylabel('SNR');
title('SNR, final model'); 
xlim([0 51]); ylim([0 10])
set(gca,'FontSize',9,'XTick',0:5:50,'YTick',0:1:15); % font 9
leg = cell(length(fr),1);
for i = 1 : length(fr)
    leg{i,1} = [num2str(fr(i)) ' Hz'];
end
leg1=legend(p,leg,'Position',[0.69 0.64 0.15 0.2]); % [left bottom width height] [0.69 0.68 0.15 0.2]
set(leg1,'FontSize',7); % 7
ah1=axes('position',get(gca,'position'),'visible','off');
leg2=legend(ah1,p1,leg,'Position',[0.43 0.64 0.15 0.2]); % [0.43 0.68 0.15 0.2]
set(leg2,'FontSize',7); % 7
