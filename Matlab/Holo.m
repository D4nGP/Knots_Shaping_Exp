%Hologram generation for SLM:

%% Parameters [um]
pix = 20; %Pixel pitch
N1 = 792; N2 = 600; N = min(N1,N2); %Sampling number
n = -N/2:N/2-1; %Pixel indices for complex filter
X = n*pix; L = pix*N; [x,y] = meshgrid(X); %Cartesian mesh
[theta,r] = cart2pol(x,y); %Polar mesh

wo = 1e3; %Beam size
Ld = 532e-3; %Incident wavelength of light
ko = 2*pi/Ld; %Wavenumber
zR = (pi*wo^2)/Ld; %Rayleigh length
u0 = 1e-3; v0 = 1*u0; carrier = @(u,v) 2*pi*(u*x+v*y); %Grating

%LG function
LG = @(p,m,z) sqrt(factorial(p)/(pi*factorial(abs(m)+p))).*r.^abs(m)./(wo^(abs(m)+1)).*exp(1i*m*theta) ...
    .*(1-1i*z/zR)^p/(1+1i*z/zR)^(p+abs(m)+1).*exp(-r.^2./(2*wo^2*(1+1i*z/zR))) ...
    .*LaguerrePoly([p,abs(m)],r.^2./(wo^2*(1+z^2/zR^2)));

%% Trefoil (Dennis and Opt6)
Uo = 1.51.*LG(0,0,0) -5.06*LG(1,0,0) +7.23*LG(2,0,0) -2.04*LG(3,0,0) -3.97*LG(0,3,0); %Trefoil Dennis Opt
% Uo = 1.2.*LG(0,0,0) -2.85*LG(1,0,0) +7.48*LG(2,0,0) -3.83*LG(3,0,0) -4.38*LG(0,3,0) -0.82*LG(1,-3,0); %Trefoil Opt6

%% Trefoil (w)
% w = 1.3;
% a00 = 1 - w^2 - 2*w^4 + 6*w^6;
% a01 = w^2*(1 + 4*w^2 - 18*w^4);
% a02 = -2*w^4*(1 - 9*w^2);
% a03 = -6*w^6;
% a30 = -8*sqrt(6)*w^3;
% atot2 = a00^2 + a01^2 + a02^2 + a03^2 + a30^2; 
% an_00 = a00/sqrt(atot2)*10;
% an_01 = a01/sqrt(atot2)*10;
% an_02 = a02/sqrt(atot2)*10;
% an_03 = a03/sqrt(atot2)*10;
% an_30 = a30/sqrt(atot2)*10;
% Uo = an_00*LG(0,0,0) +an_01*LG(1,0,0) +an_02*LG(2,0,0) +an_03*LG(3,0,0) +an_30*LG(0,3,0); % Trefoil (w)

%% Hopf (w)
% w = 1.6;
% a00 = 1 - 2*w^2 + 2*w^4;
% a01 = 2*w^2 - 4*w^4;
% a02 = 2*w^4;
% a20 = 4*sqrt(2)*w^2;
% atot2 = a00^2 + a01^2 + a02^2 + a20^2;
% an_00 = a00/sqrt(atot2)*10;
% an_01 = a01/sqrt(atot2)*10;
% an_02 = a02/sqrt(atot2)*10;
% an_20 = a20/sqrt(atot2)*10;
% Uo = an_00*LG(0,0,0) +an_01*LG(1,0,0) +an_02*LG(2,0,0)+an_20*LG(0,2,0); % Hopf(w)

%% load weights for mod trefoil and hopf
% load weights_trefoil_rotated_shifted_new.mat
% load weights_hopf_rotated_shifted_new.mat
% l = double(l); p = double(p); weight = double(weight);
% Uo = zeros(size(r));
% for ii = 1:size(weight,2)
%     Uo = Uo + weight(ii)*LG(p(ii),l(ii),0);
% end

%% Hologram generation using inverse-sinc method
ph = 0; %phase delay for phase measurements
Uz = Uo.*exp(1i*ph*pi/2);
Uz = Uz/max(max(abs(Uz)));
amp = abs(Uz); phi = angle(Uz);

C = phi -pi.*amp;
PSH1 = exp(1i.*amp.*mod(C + carrier(u0,v0),2*pi));

n0 = abs(N1-N2)/2;
holo = (angle(PSH1)+pi);
z = (holo/(2*pi))*255;
z = padarray(z,[0 n0]);

figure(1); clf(1);
imagesc(z); colormap gray; axis image;

%% LG Functions
function y = LaguerrePoly(params,x)

n=params(1)-1;
k=params(2);

L0 = 1;
L1 = 1 + k - x;
L2 = zeros(size(x));
switch (n+1)
    case 0
        y = L0;
    case 1
        y = L1;
    otherwise
        for p = 1:n
            %             disp(p)
            L2 = ((2*p + 1 + k - x).*L1 - (p+k)*L0)./(p+1);
            L0 = L1;
            L1 = L2;
        end
        y = L2;
end
end
