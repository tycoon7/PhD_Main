% Tim Coon
% Qualifying Exam Question #1
% Covering Material from MECH 719, Vibration Damping and Control
% Advisor: Dr. Richard Cobb

clear all; close all; clc;
scrsz = get(0,'ScreenSize');
scrwidth = scrsz(3);
scrheight = scrsz(4);

%% load the experimental data

load('MECH719_QualDATA.mat');
F = data(:,1);
x = data(:,2:3);

% test parameters
Fs = 5;             % (1/s) sampling frequency
dt = 1/Fs;          % (s) sample time
N = length(F);      % (-) number of samples
T = N*dt;           % (s) time duration of sample
t = (0:N-1)*dt;     % (s) time vector


% frequency vector
fnyq = 1/(2*dt);                % (1/s) nyquist frequency
% NFFT = 2^nextpow2(N);           % number of frequency lines (harmonics?)
% f = Fs/2*linspace(0,T,NFFT/2+1);     % frequency vector corresponding to PSD values
f = 0:1/T:fnyq;

%% process the amplitude spectrum

% FF = fft(F,NFFT)/N;
% for c = 1:2
%     X(:,c) = fft(x(:,c),NFFT)/T;
% end

%% process and plot the Frequency Response Function and Coherence Function

% calculate the auto-PSD of the input
Sff = FFT_PSD(F,F,T);

% initialize matrices
% Sxx = zeros(size(x));
% Sfx = Sxx;
% H = Sxx;
% C = Sxx;

for out = 1:2
    % calculate the auto-PSDs of the outputs
    Sxx(:,out) = FFT_PSD(x(:,out),x(:,out),T);
    % calculate the cross-PSDs of in/outs
    Sfx(:,out) = FFT_PSD(F,x(:,out),T);
    % calculate the FRF data
    H(:,out) = Sfx(:,out)./Sff;
    % calculate the Coherence
    C(:,out) = abs(Sfx(:,out)).^2./(Sff.*Sxx(:,out));
end

Phase = rad2deg(angle(H));
H_mag = 20*log10(abs(H));

%% Plot the resulting FRFs

for fig = 1:2
    figure(fig)
    suptitle(titles1(fig))
    % time plots
%     subplot(411)
%     plot(t,F,t,x(:,fig))
%     legend('Input','Output')
    % frequency plots
%     subplot(412)
%     semilogy(f,2*abs(F(1:NFFT/2+1)),f,abs(X(1:NFFT/2+1,fig)))
%     legend('Input','Output')
    % magnitude plots
    subplot(311)
    semilogx(f(1:length(H_mag)),H_mag(:,fig))
    ylabel('Magnitude (dB)')
    % phase plots
    subplot(312)
    plot(f(1:length(Phase)),Phase(:,fig))
    ylabel('Phase (deg)')
    % coherence plots
    subplot(313)
    plot(f(1:length(C)),C(:,fig))
    xlabel('Frequency (Hz)'); ylabel('Coherence');
end


%% Plot Simulations

m = 1;      % (kg)
c = 0.1;      % (N-s/m)
k = 1;      % (N/m)

% first-order form
A = [   0   1   0    0;
     -3*k/m 0  k/m   0;
        0   0   0    1;
      k/2*m 0 -k/m -c/2*m];
  
B = [0 1 0 0]';

C1 = eye(4);

D = zeros(4,1);

s = tf('s');
TF1 = C1*inv(s*eye(4)-A)*B;

figure(3)
bode(TF1)
    
    
    
    
    