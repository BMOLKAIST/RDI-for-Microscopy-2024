
clear; close all; clc;

%% define a support

pixel_num = 1024;

pad = 2;
thre = 2;
thick = 25;

NA1 = 235;
NA2 = 343;
NA3 = 353;

rr = ~make_ellipse(NA1 + 1, NA1 + 1, pixel_num*pad, pixel_num*pad);
rr = rr * thre;

black1 = make_ellipse(NA1, NA1, size(rr, 1), size(rr, 2));
black2 = ~make_ellipse(NA2, NA2, size(rr, 1), size(rr, 2));
black3 = ~make_ellipse(NA3, NA3, size(rr, 1), size(rr, 2));
black4 = ~(black1 .* black2) .* black3;
black4(end/2 + NA1 + 10:end, end/2 - thick + 1:end/2 + thick) = 0;

pgon2 = flipud(black4); % Fourier mask
imagesc(pgon2), colormap(parula), axis image

%% HIO algorithm

load('Figure4a_intensity');

amp = sqrt(intensity);
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

imagesc(abs(sp_recon)), colormap(gray), axis image
imagesc(angle(sp_recon)), colormap(turbo), axis image
