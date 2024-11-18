%% Simulation of electromagnetic field in a cavity
clear;
run('settings.m')

% Units:
% time - 1 round-trip
% freq - 1 reprate
% length - 1 cavity-length
% geometrical length - units of focal length
% width - pump-mode waist
% speed of light (c) - 1
% intensity (energy) - in units of saturation energy

% initializations
num_rounds = 2000; % number of simulated round-trips

% definition

A = @(w) rand(1,n);
theta = (-1+2*rand(1,n))*pi;
noise = SNR*A(w).*(cos(theta(1,:))+1i*sin(theta(1,:)));
phiKerr = @(It,W) exp((1i*Ikl.*It)/(lambda*W.^2)); % non-linear instantenous phase accumulated due to Kerr effect

% zeros matrices

Ew = zeros(num_rounds, n); % field in frequency
Et = zeros(num_rounds, n); % field in time
It = zeros(num_rounds, n); % instantaneous intensity in time
Imean = zeros(num_rounds, 1);
R = zeros(num_rounds, n); % instantaneous waist size
waist = zeros(num_rounds, n); % instantaneous waist size
q = zeros(num_rounds, n); % instantaneous waist size
F = zeros(num_rounds, n); % kerr lens focus


% figures
% % % % % % % % % % % % % % % % % % % % %
f1 = figure(1);
ax1 = subplot(1,2,1, 'Parent', f1);
ax2 = subplot(1,2,2, 'Parent', f1);

% % % % % % % % % % % % % % % % % % % % %


% starting terms

Ew(1,:) = 1e2*(-1+2*rand(1,n)+ 2i*rand(1,n) -1i)./2.;% initalize field to noise
waist(1, :) = ones(1,n)*3.2e-5;% initial waist size probably
R(1, :) =-ones(1,n)*3.0e-2;% initial waist size probably
q(1, :) = (1./R(1, :)-1i.*(lambda./(pi.*waist(1, :).^2))).^-1;

g0 = 1/mirror_loss + epsilon ; % linear gain


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for m = 2:num_rounds % for each round-trip

       
    % Initialize fields based on past round-trip
               
    Et(m - 1, :) = ifft(ifftshift(Ew(m - 1, :)));
    Et(m - 1, :) = fftshift(Et(m - 1, :));  
    It(m - 1, :) = abs(Et(m - 1, :)).^2;
  

    % Nonlinear effects calculated in time
  
    [q(m,:)] = MLSpatial( delta,Et(m - 1, :),It(m - 1, :), q(m-1,:), waist(m-1, :), Ikl, L, deltaPlane); % New waist calculation MLSpatial

    waist(m, :) = (-imag(1./q(m, :))*pi./lambda).^(-1/2);
    Et(m, :) = phiKerr(It(m - 1, :),waist(m, :)).*NLloss(waist(m,:),Wp).*Et(m - 1,:);

    % Linear effects calculated in frequency
    Ew(m, :) = fft(ifftshift(Et(m, :))); 
    Ew(m, :) = fftshift(Ew(m, :)); 
 
        
  
    g1 = SatGain(Ew(m-1,:),waist(m-1,:),g0,Is,Wp);


    W = @(w) 1./(1 + (w/spec_G_par).^2);% spectral gain function
    D = @(w) exp(-1i*disp_par*w.^2); % Dispersion
    G1 = @(w) g1.*W(w).*D(w); % Overall gain
    T1 = @(w)  1/2*(1 + mirror_loss*G1(w).*exp(-1i.*2*pi*w));%mirror_loss*G(w); %
 
  
    Ew(m, :) = T1(w).*Ew(m, :);
 
    %  creating noise
    A = @(w) rand(1,n).*(1./(1 + (w/spec_G_par).^2));
    theta = (-1+2*rand(1,n))*pi;
    noise = SNR*A(w).*(cos(theta(1,:))+1i*sin(theta(1,:))); 
%   Ew(m, :) = Ew(m, :) + sqrt(mean(abs(Ew(m,:)).^2)).*noise;
 

%  plots
   plot(ax1, t,  abs(Et(m,:)).^2);%, w, ones(1, n)*Iss);
   title('Power');

  plot(ax2, t, waist(m,:));
%   hold on
%    plot(ax2, t, waist2(m,:));
%    hold off
   title('Waist');

       
  if mod(m,1)== 0 
    [ m ]
  end 
  
   pause(0.000000000001); 



end

     
%     end 
    


