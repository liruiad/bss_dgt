%reference:
%[1]  R. R. Bond, et al., Vibration-based condition monitoring: industrial, aerospace and automotive applications, 2011.
%[2]  https://ww2.mathworks.cn/help/signal/ug/vibration-analysis-of-rotating-machinery.html
fs = 20E3;          % Sample Rate (Hz)
Np = 13;            % Number of teeth on pinion
Ng = 35;            % Number of teeth on gear
fPin = 22.5;        % Pinion (Input) shaft frequency (Hz)
fGear = fPin*Np/Ng; % Gear (Output) shaft frequency (Hz)
fMesh = fPin*Np;    % Gear Mesh frequency (Hz)

t = 0:1/fs:1-1/fs;
vfIn = 0.4*sin(2*pi*fPin*t);    % Pinion waveform     
vfOut = 0.2*sin(2*pi*fGear*t);  % Gear waveform
vMesh = sin(2*pi*fMesh*t);      % Gear-mesh waveform

% plot(t, vfIn + vfOut + vMesh)
% xlim([0 0.25])
% xlabel('Time (s)')
% ylabel('Acceleration')


ipf = fGear;
fImpact = 2000;         
tImpact = 0:1/fs:2.5e-4-1/fs; 
xImpact = sin(2*pi*fImpact*tImpact)/3;

xComb = zeros(size(t));
Ind = (0.25*fs/fMesh):(fs/ipf):length(t);
Ind = round(Ind);
xComb(Ind) = 1;
xPer = 2*conv(xComb,xImpact,'same') + 0.02*randn(size(t));
% plot(t,xPer);
% xlim([0 0.05]);
% ylabel('fault');

vNoFault = vfIn + vfOut + vMesh;
vFault = vNoFault + xPer;                              

% vNoFaultNoisy = vNoFault + randn(size(t))/25;
% vFaultNoisy = vFault + randn(size(t))/25;

% subplot(2,1,1)
% plot(t,vNoFaultNoisy)
% xlabel('Time (s)')
% ylabel('Acceleration')
% xlim([0.0 0.3])
% ylim([-2.5 2.5])
% title('Noisy Signal for Healthy Gear')
% 
% 
% subplot(2,1,2)
% plot(t,vFaultNoisy)
% xlabel('Time (s)')
% ylabel('Acceleration')
% xlim([0.0 0.3])
% ylim([-2.5 2.5])
% title('Noisy Signal for Faulty Gear')
% hold on
% MarkX = t(Ind(1:3));
% MarkY = 2.5;
% plot(MarkX,MarkY,'rv','MarkerFaceColor','red')
% hold off
n = 8;         % Number of rolling element bearings
d = 0.02;      % Diameter of rolling elements 
p = 0.12;      % Pitch diameter of bearing
thetaDeg = 15; % Contact angle in degrees

bpfi = n*fPin/2*(1+d/p*cosd(thetaDeg));% n*f_pin/2*(1+d/p*cosd(thetaDeg));104.4889
bpfo = n*fPin/2*(1-d/p*cosd(thetaDeg));
bpfr=  p/2/d*(1-d^2/p^2*cosd(thetaDeg)^2);

fImpact = 3000;
tImpact = 0:1/fs:2e-3-1/fs;
xImpact = sin(2*pi*fImpact*tImpact).*kaiser(length(tImpact),40)';


xComb = zeros(size(t));
xComb(1:round(fs/bpfi):end) = 1;
xBper_i = 0.33*conv(xComb,xImpact,'same') +0.02*randn(size(t));

figure;
plot(t,xBper_i);
xlim([0 0.05]);
ylabel('inner fault');

fImpact = 2000;
tImpact = 0:1/fs:3e-3-1/fs;
xImpact = cos(2*pi*fImpact*tImpact).*kaiser(length(tImpact),40)';

xComb = zeros(size(t));
xComb(1:round(fs/bpfo):end) = 1;
xBper_o = 0.4*conv(xComb,xImpact,'same') +0.01*randn(size(t));

figure;
plot(t,xBper_o);
xlim([0 0.05]);
ylabel('outer fault');

save xBper_i xBper_i;
save xBper_o  xBper_o;

% vNoBFaultNoisy = vNoFault + randn(size(t))/5;
% vBFaultNoisy = xBper + vNoFault + randn(size(t))/5;
% figure;
% plot(t,vBFaultNoisy);
% xlim([0 0.5]);
figure;
subplot(211)
plot(xBper_i(1:1024))
subplot(212)
plot(xBper_o(1:1024))