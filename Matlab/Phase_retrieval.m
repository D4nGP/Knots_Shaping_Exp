%Phase retrieval algorithm

ft2 = @(x) fftshift(fft2(fftshift(x)));
ift2 = @(x) fftshift(ifft2(fftshift(x)));

%Load interference measurements
for kk = 1:4
    load(sprintf('trefoilmod_ph_%d_n_%d.mat',kk,nn));
    Im(:,:,kk+1) = im2double(frame(:,:,1));
end

y1 = 1; y2 = 701; x1 = 100; x2 = 700; I = Im(y1:y2,x1:x2,:); %Crop raw snapshots, if needed
phase = atan2(I(:,:,4)-I(:,:,2),I(:,:,1)-I(:,:,3)); %Unfiltered phase
signal = (I(:,:,1)-I(:,:,3)).^2 + (I(:,:,4)-I(:,:,2)).^2; %Unfiltered intensity
U = sqrt(signal).*exp(1i*phase); %Unfiltered complex field

%Filtering the grating contribution
N = 512; pad = 2*N;
U = padarray(U,[pad,pad],'both');
U = ft2(U);
[ymax, xmax] = find(U==max(max(U)));
crop = 15; %Adjustable crop region in the Fourier space

% U = U(1129-crop:1120+crop,1238-crop:1238+crop); % To check other regions
U = U(ymax-crop:ymax+crop,xmax-crop:xmax+crop); % Automated crop

U = ift2(padarray(U,[pad/2,pad/2],'both'));
U = U(405:675,415:685); %Crop for visualization
U = U./max(max(U));

figure(1); clf(1);
imagesc(angle(U));
axis square; axis on; colormap gray; colorbar;

figure(2); clf(2);
imagesc(abs(U).^2);
axis square; axis on; colormap hot; colorbar;
