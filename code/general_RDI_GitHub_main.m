
clear; close all; clc;

%% define a sample field

pixel_num = 200;

amplitude0 = double(imread('autumn.tif'));
amplitude0 = sum(amplitude0, 3);
amplitude0 = crop_center(amplitude0, pixel_num*[1, 1]);
amplitude0 = amplitude0 ./ max(abs(amplitude0(:)));

phase0 = imread('cameraman.tif');
phase0 = double(sum(phase0, 3));
phase0 = crop_center(phase0, pixel_num*[1, 1]);
phase0 = phase0 / max(abs(phase0(:)))*6 - 3;

sp0 = amplitude0 .* exp(1i .* phase0);
imagesc(abs(sp0), [0, max(abs(sp0(:)))]), colormap(gray), axis image
imagesc(angle(sp0)), colormap(turbo), axis image

%% define a support

pad = 2;
thre = 2;
thick = 10;

NA1 = 67;
NA2 = 97;
NA3 = 100;

rr = ~make_ellipse(NA1 + 1, NA1 + 1, pixel_num*pad, pixel_num*pad);
rr = rr * thre;

black1 = make_ellipse(NA1, NA1, size(rr, 1), size(rr, 2));
black2 = ~make_ellipse(NA2, NA2, size(rr, 1), size(rr, 2));
black3 = ~make_ellipse(NA3, NA3, size(rr, 1), size(rr, 2));
black4 = ~(black1 .* black2) .* black3;
black4(end/2 + NA1 + 10:end, end/2 - thick + 1:end/2 + thick) = 0;

%% modulated field in the detector plane and NA-limited ground truth

temp = fftshift(fft2(ifftshift(sp0)));
temp = pad_center(temp, size(temp)*pad);
ground0 = temp .* black4;

ground = ground0;
ground(rr > 0) = (1 ./ 10.^(rr(rr > 0))) .* ground(rr > 0);
sp_ground = fftshift(ifft2(ifftshift(ground))); % modulated field in the detector

temp = (black2 > 0) .* ground0;
temp = crop_center(temp, pixel_num*[1, 1]);
sp_ground0 = fftshift(ifft2(ifftshift(temp))); % NA-limited ground truth
imagesc(abs(sp_ground0)), colormap(gray), axis image
imagesc(angle(sp_ground0)), colormap(turbo), axis image

pgon2 = abs(ground) > 0; % Fourier mask
imagesc(pgon2), colormap(parula), axis image

%% HIO algorithm

SNR = inf;

temp = awgn(sp_ground(:), SNR, 'measured', 'linear');
temp = reshape(temp, pixel_num*pad*[1, 1]);
amp = abs(temp);
intensity = amp.^2;
initial = amp .* exp(1i*2*pi*rand(size(amp)));
initial = fftshift(fft2(ifftshift(initial)));

beta = 0.6;

for ii = 1:1000

    if ii < 500
        lambda = 0.005;
    else
        lambda = 0.1;
    end
    
    if ii == 1 
        temp0 = initial;
    else
        temp0 = temp2;
    end
    
    til_psi = fftshift(ifft2(ifftshift(temp0)));

    amp_prime = (1 - lambda)*amp + lambda*abs(til_psi);
    temp = amp_prime .* (exp(1i .* angle(til_psi)));
    g_k_prime = fftshift(fft2(ifftshift(temp)));

    g_k = zeros(size(temp0));
    g_k(pgon2 == 1) = g_k_prime(pgon2 == 1);
    g_k(pgon2 == 0) = temp0(pgon2 == 0) - beta*g_k_prime(pgon2 == 0);

    temp2 = g_k;

    if corr_c(temp2, temp0) > 0.9999999
        break;
    end
end

for jj = 1:100
    
    temp0 = temp2;
    temp0 = fftshift(ifft2(ifftshift(temp0)));
    temp0 = amp .* exp(1i .* angle(temp0));
    temp2 = pgon2 .* fftshift(fft2(ifftshift(temp0)));
end

temp = crop_center(temp2 .* (rr > 0), pixel_num*[1, 1]);
imagesc(log10(abs(temp))), colormap(jet), axis image
sp_recon = fftshift(ifft2(ifftshift(temp))); % reconstructed field

corr_c(sp_recon, sp_ground0)

imagesc(abs(sp_recon)), colormap(gray), axis image
imagesc(angle(sp_recon)), colormap(turbo), axis image
