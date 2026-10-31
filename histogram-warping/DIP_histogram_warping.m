% DIP project work, Histogram Warping
% AI has been used in composing the idea for the code

% example images for testing the histogram warping
imgwindmill = imread('windmill.png');
imgswan = imread('swan.png');
imgankenes = imread("ankenes.png");
imgleinikki = imread('leinikki.png');

% histogram warping function in use for the example images
[windmill_warp, f_windmill, fwarp_windmill, v_wm] = histwarp(imgwindmill);
[swan_warp, f_swan, fwarp_swan, v_s] = histwarp(imgswan);
[ankenes_warp, f_ankenes, fwarp_ankenes, v_ank] = histwarp(imgankenes);
[leinikki_warp, f_leinikki, fwarp_leinikki, v_lei] = histwarp(imgleinikki);

% plot style for histograms
bins = linspace(1,256,256);
f_windmill = 600*800*f_windmill;

figure;
tiledlayout(2,2)

nexttile
bar(bins,f_windmill)
hold on;
% positions of histogram valleys plotted as red dashed lines
v_wm = v_wm * 256;
xline(v_wm, 'r--', 'LineWidth', 0.5);
title('Original image')

nexttile
bar(bins,fwarp_windmill)
title('Histogram warping')

nexttile
imagesc(imgwindmill); axis image; colormap('gray')


nexttile
imagesc(windmill_warp); axis image; colormap('gray')

figure;
tiledlayout(2,2)

nexttile
bar(bins,f_swan)
hold on
v_s = v_s * 256;
xline(v_s, 'r--', 'LineWidth',0.5)
title('Original image')

nexttile
bar(bins,fwarp_swan)
title('Histogram warping')

nexttile
imagesc(imgswan); axis image; colormap('gray')

nexttile
imagesc(swan_warp); axis image; colormap('gray')

figure;
tiledlayout(2,2)

nexttile
bar(bins,f_ankenes)
hold on
v_ank = v_ank * 256;
xline(v_ank, 'r--', 'LineWidth',0.5)
title('Original image')

nexttile
bar(bins,fwarp_ankenes)
title('Histogram warping')

nexttile
imagesc(imgankenes); axis image; colormap('gray')

nexttile
imagesc(ankenes_warp); axis image; colormap('gray')

figure;
tiledlayout(2,2)

nexttile
bar(bins,f_leinikki)
hold on
v_lei = v_lei * 256;
xline(v_lei, 'r--', 'LineWidth',0.5)
title('Original image')

nexttile
bar(bins,fwarp_leinikki)
title('Histogram warping')

nexttile
imagesc(imgleinikki); axis image; colormap('gray')

nexttile
imagesc(leinikki_warp); axis image; colormap('gray')


function [I_out, f, fwarp, v] = histwarp(I_in)



% ---- 0. Pre-processing ----
I = double(I_in);
I = I(:);

% map image to [0,1]
I = (I - min(I)) / (max(I) - min(I));
n = numel(I);

% parameters from the paper
delta  = 0.02;  % analysis resolution
tau    = 0.01;  % outlier threshold
lambda = 5;     % distortion limit (for d_k clamping later)

% ==================
% 1. Calculate smooth estimate for image histogram
% ==================

m = 256;
edges = linspace(0,1,m+1);
centers = 0.5*(edges(1:end-1)+edges(2:end));
x_grid  = centers(:);

% empirical histogram (normalized)
[counts,~] = histcounts(I, edges, 'Normalization','probability');
h = counts(:);

% smooth histogram with smoothness priors
T = length(h);
lambda_smooth = 0.5;                 % smoothing parameter
E  = speye(T);
D2 = spdiags(ones(T-2,1)*[1 -2 1],0:2,T-2,T);
f  = (E + lambda_smooth^2 * (D2' * D2)) \ h;

% enforce non-negativity
if min(f) < 0
    f = f - min(f);
end

% CDF on grid
dx     = x_grid(2) - x_grid(1);
F_grid = cumsum(f) * dx;
F_grid = F_grid / max(F_grid);   % normalize to [0,1]

% =====================
% 2. Find histogram valleys
% =====================

% derivative df/dx (finite differences)
df = gradient(f, dx);

% ignore endpoints [tau, 1-tau]
mask_inner = (x_grid >= tau) & (x_grid <= 1 - tau);
x_in  = x_grid(mask_inner);
df_in = df(mask_inner);
F_in  = F_grid(mask_inner);

v = [];

% partition CDF into equal-probability intervals of width delta
edgesF = 0:delta:1;
if edgesF(end) < 1
    edgesF(end+1) = 1;
end

for k = 1:numel(edgesF)-1
    uF = edgesF(k);
    wF = edgesF(k+1);

    % indices in the inner region
    u_idx = find(F_in >= uF, 1, 'first');
    w_idx = find(F_in <= wF, 1, 'last');
    if isempty(u_idx) || isempty(w_idx) || u_idx >= w_idx
        continue;
    end

    u = x_in(u_idx);
    w = x_in(w_idx);

    % valley condition: df<0 on left, df>0 on right
    if df_in(u_idx) < 0 && df_in(w_idx) > 0
        left  = u;
        right = w;
        % bisection on sign of derivative
        for it = 1:20
            mid = 0.5*(left + right);
            mid_idx = find(x_in >= mid, 1, 'first');
            if df_in(mid_idx) < 0
                left = mid;
            else
                right = mid;
            end
        end
        v(end+1,1) = 0.5*(left + right);
    end
end

% if no valleys detected, use global mode as pseudo-valley
if isempty(v)
    [~,idx_mode] = max(f);
    v = x_grid(idx_mode);
end

% ==========================
% 3. Histogram mapping: a_k, b_k and d_k
% ==========================

% extend valleys to include 0 and 1
v_ext = [0; v(:); 1];
K     = numel(v_ext) - 1;   % number of regions

% The paper defines a_k as midpoints of each region [v_{k-1}, v_k],
% and then defines b_k via equalization of those midpoints
a = zeros(K+1,1);   % a_0,...,a_K  (we will use indices 1..K+1)
b = zeros(K+1,1);

for k = 1:K+1
    if k == 1
        % left endpoint (corresponding to 0)
        a(k) = v_ext(1);
    elseif k == K+1
        % right endpoint (corresponding to 1)
        a(k) = v_ext(end);
    else
        % midpoint of region [v_{k-1}, v_k]
        vL   = v_ext(k);
        vR   = v_ext(k+1);
        a(k) = 0.5 * (vL + vR);
    end
end

% compute b_k using equation (mapping midpoints via equalization).[file:1]
% For endpoints, we want T(0)=0 and T(1)=1.
b(1)     = 0;
b(K+1)   = 1;
for k = 2:K
    % region k with boundaries v_ext(k) and v_ext(k+1)
    vL  = v_ext(k);
    vR  = v_ext(k+1);

    % F of endpoints and midpoint
    Fa   = interp1(x_grid, F_grid, a(k), 'linear','extrap');
    FvL  = interp1(x_grid, F_grid, vL,   'linear','extrap');
    FvR  = interp1(x_grid, F_grid, vR,   'linear','extrap');

    % equation in paper: equalization of midpoint within its local region.
    % b_k is the gray level whose CDF in the stretched domain matches the
    % normalized position of Fa within [F(vL), F(vR)].
    % Map [F(vL),F(vR)] -> [vL,vR] linearly:
    t_mid = (Fa - FvL) / max(FvR - FvL, eps);
    b(k)  = vL + t_mid * (vR - vL);
end

d = zeros(K+1,1);

% slopes at midpoints from local median mapping
% compute local "median" positions in F-space as in paper
for k = 1:K+1
    if k == 1
        % left boundary slope using right side only
        a_plus  = 0.5 * (a(k) + a(k+1));
        b_plus  = 0.5 * (b(k) + b(k+1));
        s_plus  = (b_plus - b(k)) / max(a_plus - a(k), eps);
        d(k)    = s_plus;
    elseif k == K+1
        % right boundary slope using left side only
        a_minus = 0.5 * (a(k-1) + a(k));
        b_minus = 0.5 * (b(k-1) + b(k));
        s_minus = (b(k) - b_minus) / max(a(k) - a_minus, eps);
        d(k)    = s_minus;
    else
        % interior: geometric mean of left and right slopes, weighted by F
        a_minus = 0.5 * (a(k-1) + a(k));
        a_plus  = 0.5 * (a(k)   + a(k+1));
        b_minus = 0.5 * (b(k-1) + b(k));
        b_plus  = 0.5 * (b(k)   + b(k+1));

        % slopes
        s_minus = (b(k) - b_minus) / max(a(k) - a_minus, eps);
        s_plus  = (b_plus - b(k))  / max(a_plus - a(k),  eps);

        % CDF values
        Fa    = interp1(x_grid, F_grid, a(k),     'linear','extrap');
        Fa_m  = interp1(x_grid, F_grid, a_minus,  'linear','extrap');
        Fa_p  = interp1(x_grid, F_grid, a_plus,   'linear','extrap');

        w1 = (Fa  - Fa_m) / max(Fa_p - Fa_m, eps);
        w2 = (Fa_p - Fa)  / max(Fa_p - Fa_m, eps);

        s_minus = max(s_minus, 0);
        s_plus  = max(s_plus,  0);

        d(k) = s_minus.^w1 .* s_plus.^w2; % weighted geometric mean
    end

    % clamp distortions
    d(k) = max(min(d(k), lambda), 1/lambda);
end

% ===========================
% 4. Build monotonic spline T(x) with piecewise rational quadratic
% ===========================

T_grid = zeros(size(x_grid));

for k = 1:K
    % interval [a_k, a_{k+1}] -> indices in MATLAB are [k, k+1] in 1-based
    aL = a(k);
    aR = a(k+1);
    bL = b(k);
    bR = b(k+1);
    dL = d(k);
    dR = d(k+1);

    idx = (x_grid >= aL) & (x_grid <= aR);
    x   = x_grid(idx);

    if isempty(x) || aR == aL
        continue;
    end

    % normalized t in [0,1]
    t = (x - aL) / (aR - aL);

    % r_k (secant slope between knots)
    r = (bR - bL) / (aR - aL);

    % piecewise rational quadratic spline (eq. (1))
    num = r * t.^2 + dL * (1 - t) .* t;
    den = r + (dL + dR - 2*r) * (1 - t) .* t;
    T_grid(idx) = bL + (num ./ max(den, eps)) * (bR - bL);
end

% guarantee T_grid is monotone and in [0,1]
T_grid = max(min(T_grid, 1), 0);
T_grid = cummax(T_grid);   % enforce monotonicity numerically


% original pixel values are on grid centers; interpolate T_grid
I_enh = interp1(x_grid, T_grid, I, 'linear','extrap');
I_enh = max(min(I_enh, 1), 0);

I_out = uint8(255 * I_enh);

% reshape image from vector to matrix
I_out = reshape(I_out,size(I_in));

% histogram of contrast enhanced image
[fwarp,~] = imhist(I_out);

end
