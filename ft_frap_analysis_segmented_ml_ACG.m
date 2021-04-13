% title: ft_frap_analysis_segmented_ml
% author: Andreas C Geiger
% date last updated: 2020/09/01
% purpose: this program analyzes ft-frap data
% instructions:
%   1) change user-defined inputs
%   2) Hit "Run"

%% user-defined inputs
% path = 'Z:\AndreasGeiger\FRAP in sharedlab\03192021\NR 1hour 1PF 532nm 120 bleach\fov2\'; % folder name
% file = 'ChannelB Averaging RT-19-Mar-2021 17_06_14.bin'; % file name
% prebleach = 10;
% t0 = 131; % first frame of fluorescence recovery
% t_end = 1000;
% ft_peak_1 = 15;  % number of lines - 1
% ft_peak_2 = 15 * 2;
% ft_peak_3 = 15 * 3;
% vpix = 512; % number of pixels in vertical axis
% hpix = 512; % number of pixels in horizontal axis
% microns = 512 * 1;% number of microns in one axis

path = '\\10.164.16.234\Data\DustinHarmon\03312021\Mazin-A1\'; % folder name
file = '31-Mar-2021 11_55_09 GalvoRes ChanA Counting Pol1 10000Images.bin'; % file name
prebleach = 14;
t0 = 39; % first frame of fluorescence recovery
t_end = 1050;
ft_peak_1 = 31;  % number of lines - 1
ft_peak_2 = 31 * 2;
ft_peak_3 = 31 * 3;
vpix = 512; % number of pixels in vertical axis
hpix = 512; % number of pixels in horizontal axis
microns = 512 * 0.625;% number of microns in one axis

%% import data
FID = fopen(strcat(path, file));
image = fread(FID, 'float32');
image = reshape(image, vpix, hpix, []);
% image_norm = image ./ mean(image(:, :, 1:prebleach), 3);

%%
image_mask = mean(image(:, :, 1:prebleach), 3);
imageSegmenter;
imagesc(BW) ; mask_segment(:, :, 1) = BW;  mask_segment(:,:,2) = -1 * (BW- 1);

%%
image_1 = image .* mask_segment(:, :, 1);
image_1 = image_1 ./ mean(image_1(:, :, 1:prebleach), 3);
image_1(isnan(image_1)) = 0;
%% extract fluorescence recovery with fourier transform
clear ft
% compute 2d fourier transform
for i = 1:t_end %numel(image(1, 1, :))
    ft(:, :, i) = fft2(image_1(:, :, i));
%     ft_norm(:, :, i) = fft2(image_norm(:, :, i));
    %ft(:, :, i) = fftshift(ft(:, :, i));
end
ft_norm = ft;

%%
clear ft
power = sqrt(ft_norm .* conj(ft_norm)); % compute power from fourier transform
peak_x = 1:3;
peak_y_1 = ft_peak_1:ft_peak_1 + 2;
peak_y_2 = ft_peak_2:ft_peak_2 + 2;
peak_y_3 = ft_peak_3:ft_peak_3 + 2;

%% heterogeneity analysis
ft_crop_1 = ft_norm(:, ft_peak_1 + 1:end - ft_peak_1, :);
ft_crop_1(:, round(ft_peak_1/4):end - round(ft_peak_1/4), :) = 0;
ft_crop_1(round(ft_peak_1/4):end - round(ft_peak_1/4), :, :) = 0;

ft_crop_2 = ft_norm(:, ft_peak_2 + 1:end - ft_peak_2, :);
ft_crop_2(:, round(ft_peak_1/4):end - round(ft_peak_1/4), :) = 0;
ft_crop_2(round(ft_peak_1/4):end - round(ft_peak_1/4), :, :) = 0;

ft_crop_3 = ft_norm(:, ft_peak_3 + 1:end - ft_peak_3, :);
ft_crop_3(:, round(ft_peak_1/4):end - round(ft_peak_1/4), :) = 0;
ft_crop_3(round(ft_peak_1/4):end - round(ft_peak_1/4), :, :) = 0;

%%
tic
for i = 1:numel(ft_norm(1, 1, :))
    ift_crop_1(:, :, i) = ifft2(ft_crop_1(:, :, i) , 512, 512);
    diff_im_1(:, :, i) = sqrt(ift_crop_1(:, :, i) .* conj(ift_crop_1(:, :, i)));
%     diff_im_gauss_1(:, :, i) = imgaussfilt(diff_im_1(:, :, i), 10);
    
    ift_crop_2(:, :, i) = ifft2(ft_crop_2(:, :, i) , 512, 512);
    diff_im_2(:, :, i) = sqrt(ift_crop_2(:, :, i) .* conj(ift_crop_2(:, :, i)));
%     diff_im_gauss_2(:, :, i) = imgaussfilt(diff_im_2(:, :, i), 10);
    
    ift_crop_3(:, :, i) = ifft2(ft_crop_3(:, :, i) , 512, 512);
    diff_im_3(:, :, i) = sqrt(ift_crop_3(:, :, i) .* conj(ift_crop_3(:, :, i)));
%     diff_im_gauss_3(:, :, i) = imgaussfilt(diff_im_3(:, :, i), 10);
end
toc
%%
for i = 1:numel(mask_segment(1, 1, :))
    diff_region_1(i, :) = mean(mean(diff_im_1 .* mask_segment(:, :, i), 1), 2);
    diff_region_2(i, :) = mean(mean(diff_im_2 .* mask_segment(:, :, i), 1), 2);
    diff_region_3(i, :) = mean(mean(diff_im_3 .* mask_segment(:, :, i), 1), 2);
end

%%
A_1 = diff_region_1(:, t0);
a1_1 = mean(diff_region_1(:, 1:prebleach), 2);
diff_1 = (diff_region_1(:, t0:end) - a1_1) ./ (A_1 - a1_1);

A_2 = diff_region_2(:, t0);
a1_2 = mean(diff_region_2(:, 1:prebleach), 2);
diff_2 = (diff_region_2(:, t0:end) - a1_2) ./ (A_2 - a1_2);

A_3 = diff_region_3(:, t0);
a1_3 = mean(diff_region_3(:, 1:prebleach), 2);
diff_3 = (diff_region_3(:, t0:end) - a1_3) ./ (A_3 - a1_3);

%%
lb = [0, -0.5, -0.5, -0.5, 0, 0];
ub = [inf, 1, 1, 1, 2, 4];

diff_map = zeros(512, 512);
s_diff_map = zeros(512, 512);
rec_map = zeros(512, 512);
alpha_map = zeros(512, 512);
s_alpha_map = zeros(512, 512);
mu_map = zeros(512, 512);
s_mu_map = zeros(512, 512);

frap_1 = diff_1;
frap_2 = diff_2;
frap_3 = diff_3;

t = 4/17 * (0:(length(frap_1) - 1));% time in seconds
%%
for i = 1:numel(mask_segment(1, 1, :))
    
    fun = @(r) cat(2, (1 - r(2)) * ml(-t.^(2 * r(5)/r(6)) * (2 * pi * ft_peak_1 / 512)^(r(6)) * r(1), r(5)) + r(2) - frap_1(i, :),...
        (1 - r(3)) * ml(-t.^(2 * r(5)/r(6)) * (2 * pi * ft_peak_2 / 512)^(r(6)) * r(1), r(5)) + r(3) - frap_2(i, :),...
        (1 - r(4)) * ml(-t.^(2 * r(5)/r(6)) * (2 * pi * ft_peak_3 / 512)^(r(6)) * r(1), r(5)) + r(4) - frap_3(i, :));
    
    r0 = rand(1, 6);
    [r, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(fun, r0, lb, ub);
    % r(1) = Diffusion Coefficient, r(5) = alpha parameter, r(6) = mu parameter
    
    residual = reshape(residual, 1, []);
    var_res = var(residual, [], 2);
    J = full(jacobian);
    J = J(:, 1:numel(r));
    covar = inv(J' * J) * var_res;
    
    if sum(diag(covar) >= 0) == numel(r0)
        [sigma, ~] = cov2corr(covar);
        
        diff_map = diff_map + r(1) * (microns)^2 / (vpix)^2 * mask_segment(:, :, i);
        s_diff_map = s_diff_map + sigma(1) * (microns)^2 / (vpix)^2 * mask_segment(:, :, i);
        rec_map = rec_map + r(2) * mask_segment(:, :, i);
        alpha_map = alpha_map + r(5) * mask_segment(:, :, i);
        s_alpha_map = s_alpha_map + sigma(5) * mask_segment(:, :, i);
        mu_map = mu_map + r(6) * mask_segment(:, :, i);
        s_mu_map = s_mu_map + sigma(6) * mask_segment(:, :, i);
        
        fit(:, 1, i) = (1 - r(2)) * ml(-t.^(2 * r(5)/r(6)) * (2 * pi * ft_peak_1 / 512)^(r(6)) * r(1), r(5)) + r(2);
        fit(:, 2, i) = (1 - r(3)) * ml(-t.^(2 * r(5)/r(6)) * (2 * pi * ft_peak_2 / 512)^(r(6)) * r(1), r(5)) + r(3);
        fit(:, 3, i) = (1 - r(4)) * ml(-t.^(2 * r(5)/r(6)) * (2 * pi * ft_peak_3 / 512)^(r(6)) * r(1), r(5)) + r(4);

    end
end

% ratio_map = diff_map ./ s_diff_map;
% diff_map_filt = diff_map;
% diff_map_filt(ratio_map < 3) = NaN;