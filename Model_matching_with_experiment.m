% Code (1) for manuscript: "Method for overcoming temporal dispersion in unmyelinated
% nerves to image them with Electrical Impedance Tomography (EIT)"

% Computation of dispersed CAPs and dZs to match the model with experimental data
% - by Ilya Tarotin, 2021

clear;

coef_cond = (0.97*1.211+0.03*47.8)/0.1;
coef_elV = ((1-1/sqrt(2))*1600/10)^-3;
coefVAP = coef_elV*coef_cond;

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

dt = 0.01;

% C
vC = 0.8; devC = 0.3; % Based on CAPs from ex-vivo pig subdiaphragmatic
dC = 1; % um, needed for dZ amplitudes, Soltanpour - 0.8
numC = 40000;
latC = 2; % ms, Ghitani
vC0 = []; % Distribution of velocities
vC0 = [vC0;devC*randn(5*numC,1)+vC];
vC0(vC0<0.1)=[];
if length(vC0)> numC % to have exact number of fibres
    vC0 = vC0(1:numC);
end

% Ad 
vAd = 8; devAd = 3; dAd = 4; % Based on CAPs from ex-vivo pig subdiaphragmatic, dAd = 4
numAd = 4000;
latAd = 1; % ms, Ghitani
vAd0 = [];
vAd0 = [vAd0;devAd*randn(5*numAd,1)+vAd];
vAd0(vAd0<0.1)=[];
if length(vAd0)> numAd % to have exact number of fibres
    vAd0 = vAd0(1:numAd);
end

% Distances to build the signal
dist = [30 150 200 500]; % mm 
clear tVC tVAd
cnt = 1;
for i = dist

    tAd(:,cnt) = dist(cnt)./vAd0;
    tC(:,cnt) = dist(cnt)./vC0;

    cnt = cnt + 1;
end

% Signals with FEM model with CUFF
load('ap_shape_C_cuff.mat'); % AP shape obtained with FEM model. Example uploaded to GitHub.
ap_shape_C0 = ap_shape_C_arr(1200:1800,3)';
dz_shape_C0 = ap_shape_C_arr(1200:1800,5)'-0.79878; % BV = 0.79878

%{
% Simple plots of CAP and dZ 
Fs = 1e5;t = 1000:2200;
figure; % single CAP
plot(1000.*(1:length(ap_shape_C_arr(t,5)))./Fs,ap_shape_C_arr(t,3),'linewidth',1.5); 
figure; % single dZ
plot(1000.*(1:length(ap_shape_C_arr(t,5)))./Fs,ap_shape_C_arr(t,5)-0.79878,'linewidth',1.5); 

% Old model
figure;plot(1000.*(1:length(ap_c_noI(t)))./Fs,ap_c_noI(t),'linewidth',1.5);
figure;plot(1000.*(1:length(c_dz_abs(t)))./Fs,c_dz_abs(t),'linewidth',1.5);
%}

dz_shape_C = imresize(((dC/1)^2).*dz_shape_C0,1); % ()^2 as dZ~cross-sectional area, latency?
ap_shape_C = imresize(((dC/1)^1).*ap_shape_C0,1); %()^1 as AP~circumference
dz_shape_Ad = imresize(((dAd/1)^2).*dz_shape_C0,(latAd/latC)^1); % ()^2 as dZ~cross-sectional area. For myelinated - use diameter of whole fibre with myelin
ap_shape_Ad = imresize(ap_shape_C0,latAd/latC); % ap Ad same size as C. Remaining area is assumed to be myelin


T_fin = max(tC(:,end));
t_sim = 0:dt:T_fin+5*latC/dt; % time of simulation, 

% Sum all dZ squares which came at the specified times
sigV = zeros(length(t_sim),length(dist));
for k = 1 : length(dist)
    for i = 1 : size(tC,1)
        sigV(round(tC(i,k)/dt):round(tC(i,k)/dt)+length(dz_shape_C)-1,k) = sigV(round(tC(i,k)/dt):round(tC(i,k)/dt)+length(dz_shape_C)-1,k) + dz_shape_C';
    end
    for i = 1 : size(tAd,1)
        sigV(round(tAd(i,k)/dt):round(tAd(i,k)/dt)+length(dz_shape_Ad(1,:))-1,k) = sigV(round(tAd(i,k)/dt):round(tAd(i,k)/dt)+length(dz_shape_Ad(1,:))-1,k) + dz_shape_Ad';
    end
end
sigVAP = zeros(length(t_sim),length(dist));
for k = 1 : length(dist)
    for i = 1 : size(tC,1)
        sigVAP(round(tC(i,k)/dt):round(tC(i,k)/dt)+length(ap_shape_C)-1,k) = sigVAP(round(tC(i,k)/dt):round(tC(i,k)/dt)+length(ap_shape_C)-1,k) + ap_shape_C';
    end
    for i = 1 : size(tAd,1)
        sigVAP(round(tAd(i,k)/dt):round(tAd(i,k)/dt)+length(ap_shape_Ad)-1,k) = sigVAP(round(tAd(i,k)/dt):round(tAd(i,k)/dt)+length(ap_shape_Ad)-1,k) + ap_shape_Ad';
    end
end

% Noise
% noise = 0*5e-4.*randn(length(dist),length(t_sim)); % K Aristovich et al, 2018
noise = 10e-3.*randn(length(dist),1*length(t_sim)); % NOISE for CAPs, as in experiment before filter/averaging (5-10 uV)
noisedZ = noise/sqrt(1200); % Noise for dZ (averaged CAPs = 10 mins, 2Hz = 1200 averages)


ndz = 0; %2e-4;
sigV = sigV.*coefVdZ+ndz.*noise(1:end,:)';
sigVAP = sigVAP.*coefVAP+noise(1:end,:)';

% Filtering
Fs = 1000/dt;
[b, a] = butter(1, 200/(Fs/2),'low'); % Fs = 20e3, 1st order butter
sigV_filt = filtfilt(b,a,sigVAP(:,:));
sigVdZ_filt = filtfilt(b,a,sigV(:,:));

%% Figure dZs
t_end = 250;
dz_end_min = -6; % uV
dz_end_max = 2;
dist_cell = sprintfc('%d',dist./10); 
dt_pl = 50; dz_pl = 0.5;
k_uv = 1000; % mV -> uV
figure;
leg = cell(length(dist),1);
Color = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940, 0.1840, 0.5560;0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];
for i = 1 % length(dist):-1:1
% subplot(1,4,i);
if i == 1
    plot(t_sim,k_uv*sigV(:,i),'linewidth',2,'Color',Color(i,:)); % ,'Color',[0.8500 0.3250 0.0980]
else
    plot(t_sim,k_uv*sigV(:,i),'linewidth',0.5,'Color',Color(i,:));
end
hold on;
xlim([0 200]); ylim([-1 0.4]);
end
for i = length(dist):-1:1
    leg{i,1} = [num2str(dist(i)./10) ' cm'];
end
xlabel('Time (ms)'); ylabel('\muV');
title('dZ, modified model');
set(gca,'fontsize',8,'ytick',-1:0.2:1);
% legend('3 cm','15 cm','20 cm','50 cm','location','southeast');
%% Plot APs
figure;
plot(t_sim,k_uv*sigV_filt(:,1),'Color',[0.8500 0.3250 0.0980],'linewidth',2); %,'Color',[0.3010 0.7450 0.9330]
xlim([0 200]);
ylim([-100 100]);
title('CAP, modified model');
xlabel('Time (ms)'); ylabel('\muV');
set(gca,'fontsize',8,'ytick',-200:20:200);
% legend('3 cm','15 cm','20 cm','50 cm','location','southeast');
