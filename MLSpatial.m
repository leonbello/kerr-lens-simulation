function [qt] = MLSpatial(delta, Et,Pt, q1, W1, Ikl ,L,deltaPlane)


deltaPoint = delta - deltaPlane;

lambda = 780e-9;
RM = 150e-3;
FM = 75e-3;
L1 = 0.5;
L2 = 0.9;
V = 1/(2/RM-1/L2);
N = 5; % number of NL lenses

Mcur = @(RM)[ 1 0 ; -2/RM 1];
distance = @(d)[ 1 d ; 0 1];
lens = @(fL)[ 1 0 ; -1./fL 1];
lensL = @(fL,fl)[ 1 0 ; -1./fL-1i./fl 1];

lens_aperture = 56e-6;

f = ((2*pi*lens_aperture^2)/lambda);


phiKerr = @(Pt,W) exp((1i*Ikl.*Pt)./(lambda*W.^2)); % non-linear instantenous phase accumulated due to Kerr effect


for i=1:length(Et)
    
    
    
    Feff11(i) = (W1(i).^4)./(Ikl.*Pt(i));  
    M(:,:,i) =  distance(L/(N))*lensL(Feff11(i),f);
    
    %       % Lens #2
    
    q2(i) =  (M(1,1,i).*q1(i)+M(1,2,i))./(M(2,1,i).*q1(i)+M(2,2,i));
    W2(i) = (-imag(1./q2(i))*pi./lambda).^(-1/2);
    Et2(i) = phiKerr(abs(Et(i)).^2,W1(i)).*Et(i);    
    Pt2(i) = abs(Et2(i)).^2;   
    Feff21(i) = (W2(i).^4)./(Ikl.*Pt2(i));  
    
    M(:,:,i) = distance(L/(N))*lensL(Feff21(i),f);
    
    %     % Lens #3
    
    q3(i) =  (M(1,1,i).*q2(i)+M(1,2,i))./(M(2,1,i).*q2(i)+M(2,2,i));    
    W3(i) = (-imag(1./q3(i))*pi./lambda).^(-1/2);    
    Et3(i) = phiKerr(abs(Et2(i)).^2,W3(i)).*Et2(i);    
    Pt3(i) = abs(Et3(i)).^2;    
    Feff31(i) = (W3(i).^4)./(Ikl.*Pt3(i));
    
    M(:,:,i) = distance(L/(N))*lensL(Feff31(i),f);
    
    %     % Lens #4
    
    q4(i) =  (M(1,1,i).*q3(i)+M(1,2,i))./(M(2,1,i).*q3(i)+M(2,2,i));
    W4(i) = (-imag(1./q4(i))*pi./lambda).^(-1/2);
    Et4(i) = phiKerr(abs(Et2(i)).^2,W4(i)).*Et3(i);
    Pt4(i) = abs(Et4(i)).^2;
    Feff41(i) = (W4(i).^4)./(Ikl.*Pt4(i));   
    
    M(:,:,i) = distance(L/(N))*lensL(Feff41(i),f);
    
    %     % Lens #5
    
    q5(i) =  (M(1,1,i).*q4(i)+M(1,2,i))./(M(2,1,i).*q4(i)+M(2,2,i));    
    W5(i) = (-imag(1./q5(i))*pi./lambda).^(-1/2);    
    Et5(i) = phiKerr(abs(Et4(i)).^2,W5(i)).*Et4(i);
    Pt5(i) = abs(Et5(i)).^2;
    Feff51(i) = (W5(i).^4)./(Ikl.*Pt5(i));
  
    M(:,:,i) = distance(L/(N * 2))*distance(RM/2+deltaPlane-eps-L/2)*Mcur(RM)*distance(L1)*[1 0; 0 1]*distance(L1)*Mcur(RM)*distance(RM/2+deltaPlane-eps-L/2)*distance(L/(N * 2))*lensL(Feff51(i),f);
    
    % Lens #52
    
    q5(i) =  (M(1,1,i).*q5(i)+M(1,2,i))./(M(2,1,i).*q5(i)+M(2,2,i));    
    W5(i) = (-imag(1./q5(i))*pi./lambda).^(-1/2);
    Et5(i) = phiKerr(abs(Et5(i)).^2,W5(i)).*Et5(i);    
    Pt5(i) = abs(Et5(i)).^2;
    Feff52(i) = (W5(i).^4)./(Ikl.*Pt5(i));

    M(:,:,i) = distance(L/(N))*lensL(Feff52(i),f);
    
    % Lens #42
    
    q4(i) =  (M(1,1,i).*q5(i)+M(1,2,i))./(M(2,1,i).*q5(i)+M(2,2,i));
    W4(i) = (-imag(1./q4(i))*pi./lambda).^(-1/2);
    Et4(i) = phiKerr(abs(Et5(i)).^2,W4(i)).*Et5(i);
    Pt4(i) = abs(Et4(i)).^2;
    Feff42(i) = (W4(i).^4)./(Ikl.*Pt4(i));
    
    M(:,:,i) = distance(L/(N))*lensL(Feff42(i),f);
    
    % Lens #32
    
    q3(i) =  (M(1,1,i).*q4(i)+M(1,2,i))./(M(2,1,i).*q4(i)+M(2,2,i));    
    W3(i) = (-imag(1./q3(i))*pi./lambda).^(-1/2);
    Et3(i) = phiKerr(abs(Et4(i)).^2,W3(i)).*Et4(i);
    Pt3(i) = abs(Et3(i)).^2;
    Feff32(i) = (W3(i).^4)./(Ikl.*Pt3(i));
    
    M(:,:,i) = distance(L/(N))*lensL(Feff32(i),f);
    
    % Lens #22
    
    q2(i) =  (M(1,1,i).*q3(i)+M(1,2,i))./(M(2,1,i).*q3(i)+M(2,2,i));
    W2(i) = (-imag(1./q2(i))*pi./lambda).^(-1/2);
    Et2(i) = phiKerr(abs(Et3(i)).^2,W2(i)).*Et3(i);    
    Pt2(i) = abs(Et2(i)).^2;
    Feff22(i) = (W2(i).^4)./(Ikl.*Pt2(i));
    
    M(:,:,i) = distance(L/(N))*lensL(Feff22(i),f);
    
    % Lens #12
    
    q1(i) =  (M(1,1,i).*q2(i)+M(1,2,i))./(M(2,1,i).*q2(i)+M(2,2,i));    
    W1(i) = (-imag(1./q1(i))*pi./lambda).^(-1/2);    
    Et(i) = phiKerr(abs(Et2(i)).^2,W1(i)).*Et2(i);
    Pt(i) = abs(Et(i)).^2;    
    Feff12(i) = (W1(i).^4)./(Ikl.*Pt(i));
    
    M(:,:,i) = distance(L/(N * 2))*distance(V+deltaPoint-L/2)*lens(FM)*distance(L2)*[1 0; 0 1]*distance(L2)*lens(FM)*distance(V+deltaPoint-L/2)*distance(L/(N * 2))*lensL(Feff12(i),f);
     
    qt(i) =  (M(1,1,i).*q1(i)+M(1,2,i))./(M(2,1,i).*q1(i)+M(2,2,i));
    
   
    
    
end
end


