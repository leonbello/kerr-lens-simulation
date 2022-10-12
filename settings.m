n = 2^nextpow2(2000); % number of simulated time-bins, 
                      % we take this to be a power of 2 for efficiency of the FFT
bw = n; % simulated bandwidth
n = n + 1; % to make the space between frequencies 1
w = linspace(-bw/2, bw/2, n); % frequency is in units of reprate, time is in units of round-trip time
dw = bw/(n - 1);
t = linspace(-1/(2*dw), 1/(2*dw), n);
dt = 1/(n*dw);


 % Kerr lens parameters    
 
n2 = 3e-20; % n2 of sapphire m^2/W
L = 3e-3; % crystal length
kerr_par = 4*L*n2;%
N = 5; % number of NL lenses in the crystal
Ikl = kerr_par/N/50;%kerr_par
Is = 2.6*n^2*500; % saturation power
Wp = 30e-6;


mirror_loss = 0.95; % loss of the OC
spec_G_par = 200; %gaussian linewidth parameter
SNR = 0e-3;
lambda = 780e-9;
delta = 0.001; % how far we go into the stability gap 

deltaPlane = -0.75e-3; % position of crystal - distance from the "plane" lens focal
disp_par = 0*1e-3*2*pi/spec_G_par; % net dispersion 
epsilon = 0.2;