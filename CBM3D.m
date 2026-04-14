function [PSNR, yRGB_est] = CBM3D(yRGB, zRGB, sigma, profile, print_to_screen, colorspace)

image_name = [

   'input.png'

   ];

if (exist('profile') ~= 1)
    profile         = 'np';
end

if (exist('sigma') ~= 1),
   sigma                = 50;
end

if (exist('colorspace') ~= 1),
    colorspace              = 'opp';
end

transform_2D_HT_name     = 'bior1.5';
transform_2D_Wiener_name = 'dct';
transform_3rd_dim_name   = 'haar';

N1                  = 8;
Nstep               = 2;
N2                  = 16;
Ns                  = 39;
tau_match           = 3000;
lambda_thr2D        = 0;
lambda_thr3D        = 2.6;
beta                = 2.0;

N1_wiener           = 8;
Nstep_wiener        = 2;
N2_wiener           = 32;
Ns_wiener           = 39;
tau_match_wiener    = 400;
beta_wiener         = 2.0;

stepFS              = 1;
smallLN             = 'not used in np';
stepFSW             = 1;
smallLNW            = 'not used in np';
thrToIncStep        = 8;

if strcmp(profile, 'lc') == 1,

    Nstep               = 6;
    Ns                  = 25;
    Nstep_wiener        = 5;
    N2_wiener           = 16;
    Ns_wiener           = 25;

    thrToIncStep        = 3;
    smallLN             = 3;
    stepFS              = 6*Nstep;
    smallLNW            = 2;
    stepFSW             = 5*Nstep_wiener;

end

if (strcmp(profile, 'vn') == 1) | (sigma > 40),

    N2                  = 32;
    Nstep               = 4;
 
    N1_wiener           = 11;
    Nstep_wiener        = 6;

    lambda_thr3D        = 2.8;
    thrToIncStep        = 3;
    tau_match_wiener    = 3500;
    tau_match           = 25000;
    
    Ns_wiener           = 39;
    
end

if (strcmp(profile, 'vn_old') == 1) & (sigma > 40),

    transform_2D_HT_name = 'dct'; 
    
    N1                  = 12;
    Nstep               = 4;
 
    N1_wiener           = 11;
    Nstep_wiener        = 6;

    lambda_thr3D        = 2.8;
    lambda_thr2D        = 2.0;
    thrToIncStep        = 3;
    tau_match_wiener    = 3500;
    tau_match           = 5000;
    
    Ns_wiener           = 39;
    
end

decLevel = 0;        
thr_mask = ones(N1); 

if strcmp(profile, 'high') == 1,
    
    decLevel     = 1; 
    Nstep        = 2;
    Nstep_wiener = 2;
    lambda_thr3D = 2.5;
    vMask = ones(N1,1); vMask((end/4+1):end/2)= 1.01; vMask((end/2+1):end) = 1.07; 
    thr_mask = vMask * vMask'; 
    beta         = 2.5;
    beta_wiener  = 1.5;
    
end

dump_output_information = 1;
if (exist('print_to_screen') == 1) & (print_to_screen == 0),
    dump_output_information = 0;
end

[Tfor, Tinv]   = getTransfMatrix(N1, transform_2D_HT_name, decLevel);  
[TforW, TinvW] = getTransfMatrix(N1_wiener, transform_2D_Wiener_name);

if (strcmp(transform_3rd_dim_name, 'haar') == 1) | (strcmp(transform_3rd_dim_name(end-2:end), '1.1') == 1),
    
    hadper_trans_single_den         = {};
    inverse_hadper_trans_single_den = {};
else
    
    for hpow = 0:ceil(log2(max(N2,N2_wiener))),
        h = 2^hpow;
        [Tfor3rd, Tinv3rd]   = getTransfMatrix(h, transform_3rd_dim_name, 0);
        hadper_trans_single_den{h}         = single(Tfor3rd);
        inverse_hadper_trans_single_den{h} = single(Tinv3rd');
    end
end

if beta_wiener==2 & beta==2 & N1_wiener==8 & N1==8
    Wwin2D = [ 0.1924    0.2989    0.3846    0.4325    0.4325    0.3846    0.2989    0.1924;
        0.2989    0.4642    0.5974    0.6717    0.6717    0.5974    0.4642    0.2989;
        0.3846    0.5974    0.7688    0.8644    0.8644    0.7688    0.5974    0.3846;
        0.4325    0.6717    0.8644    0.9718    0.9718    0.8644    0.6717    0.4325;
        0.4325    0.6717    0.8644    0.9718    0.9718    0.8644    0.6717    0.4325;
        0.3846    0.5974    0.7688    0.8644    0.8644    0.7688    0.5974    0.3846;
        0.2989    0.4642    0.5974    0.6717    0.6717    0.5974    0.4642    0.2989;
        0.1924    0.2989    0.3846    0.4325    0.4325    0.3846    0.2989    0.1924];
    Wwin2D_wiener = Wwin2D;
else
    Wwin2D           = kaiser(N1, beta) * kaiser(N1, beta)';
    Wwin2D_wiener    = kaiser(N1_wiener, beta_wiener) * kaiser(N1_wiener, beta_wiener)';
end

if (exist('yRGB') ~= 1) | (exist('zRGB') ~= 1)
    yRGB        = im2double(imread(image_name)); 
    randn('seed', 0);                         
    zRGB        = yRGB + (sigma/255)*randn(size(yRGB)); 
else 
    image_name = 'External image';
    
    zRGB = double(zRGB);

    yRGB = double(yRGB);
    
    if (max(zRGB(:)) > 10),
        zRGB = zRGB / 255;
    end
    
    if (max(yRGB(:)) > 10), % a naive check for intensity range
        yRGB = yRGB / 255;
    end    
end


if (size(zRGB,3) ~= 3) | (size(zRGB,4) ~= 1),
    error('CBM3D accepts only input RGB images (i.e. matrices of size M x N x 3).');
end

yRGB_is_invalid_image = (length(size(zRGB)) ~= length(size(yRGB))) | (size(zRGB,1) ~= size(yRGB,1)) | (size(zRGB,2) ~= size(yRGB,2)) | (size(zRGB,3) ~= size(yRGB,3));
if (yRGB_is_invalid_image),
    dump_output_information = 0;
end


[Xv, Xh, numSlices] = size(zRGB);             

if numSlices ~= 3
    fprintf('Error, an RGB color image is required!\n');
    return;
end

[zColSpace l2normLumChrom] = function_rgb2LumChrom(zRGB, colorspace);

if dump_output_information == 1,
    fprintf(sprintf('Image: %s (%dx%dx%d), sigma: %.1f\n', image_name, Xv, Xh, numSlices, sigma));
end


tic;
y_hat = bm3d_thr_color(zColSpace, hadper_trans_single_den, Nstep, N1, N2, lambda_thr2D,...
    lambda_thr3D, tau_match*N1*N1/(255*255), (Ns-1)/2, sigma/255, thrToIncStep, single(Tfor), single(Tinv)', inverse_hadper_trans_single_den, single(thr_mask), 'unused arg', 'unused arg', l2normLumChrom, Wwin2D, smallLN, stepFS );
estimate_elapsed_time = toc;


tic;
yRGB_est = bm3d_wiener_color(zColSpace, y_hat, hadper_trans_single_den, Nstep_wiener, N1_wiener, N2_wiener, ...
    'unused_arg', tau_match_wiener*N1_wiener*N1_wiener/(255*255), (Ns_wiener-1)/2, sigma/255, 'unused arg', single(TforW), single(TinvW)', inverse_hadper_trans_single_den, 'unused arg', 'unused arg', l2normLumChrom, Wwin2D_wiener, smallLNW, stepFSW );
wiener_elapsed_time = toc;

yRGB_est = double(yRGB_est);

yRGB_est = function_LumChrom2rgb(yRGB_est, colorspace);

PSNR = 0;
if (~yRGB_is_invalid_image),
    PSNR = 10*log10(1/mean((yRGB(:)-yRGB_est(:)).^2));
end

if dump_output_information == 1,
    fprintf(sprintf('FINAL ESTIMATE (total time: %.1f sec), PSNR: %.2f dB\n', ...
        wiener_elapsed_time + estimate_elapsed_time, PSNR));

    figure, imshow(min(max(zRGB,0),1)); title(sprintf('Noisy %s, PSNR: %.3f dB (sigma: %d)', ...
        image_name(1:end-4), 10*log10(1/mean((yRGB(:)-zRGB(:)).^2)), sigma));

    figure, imshow(min(max(yRGB_est,0),1)); title(sprintf('Denoised %s, PSNR: %.3f dB', ...
        image_name(1:end-4), PSNR));
end

return;

function [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)

if exist('dec_levels') ~= 1,
    dec_levels = 0;
end

if N == 1,
    Tforward = 1;
elseif strcmp(transform_type, 'hadamard') == 1,
    Tforward    = hadamard(N);
elseif (N == 8) & strcmp(transform_type, 'bior1.5')==1
    Tforward =  [ 0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274;
       0.219417649252501   0.449283757993216   0.449283757993216   0.219417649252501  -0.219417649252501  -0.449283757993216  -0.449283757993216  -0.219417649252501;
       0.569359398342846   0.402347308162278  -0.402347308162278  -0.569359398342846  -0.083506045090284   0.083506045090284  -0.083506045090284   0.083506045090284;
      -0.083506045090284   0.083506045090284  -0.083506045090284   0.083506045090284   0.569359398342846   0.402347308162278  -0.402347308162278  -0.569359398342846;
       0.707106781186547  -0.707106781186547                   0                   0                   0                   0                   0                   0;
                       0                   0   0.707106781186547  -0.707106781186547                   0                   0                   0                   0;
                       0                   0                   0                   0   0.707106781186547  -0.707106781186547                   0                   0;
                       0                   0                   0                   0                   0                   0   0.707106781186547  -0.707106781186547];   
elseif (N == 8) & strcmp(transform_type, 'dct')==1
    Tforward = [ 0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274;
       0.490392640201615   0.415734806151273   0.277785116509801   0.097545161008064  -0.097545161008064  -0.277785116509801  -0.415734806151273  -0.490392640201615;
       0.461939766255643   0.191341716182545  -0.191341716182545  -0.461939766255643  -0.461939766255643  -0.191341716182545   0.191341716182545   0.461939766255643;
       0.415734806151273  -0.097545161008064  -0.490392640201615  -0.277785116509801   0.277785116509801   0.490392640201615   0.097545161008064  -0.415734806151273;
       0.353553390593274  -0.353553390593274  -0.353553390593274   0.353553390593274   0.353553390593274  -0.353553390593274  -0.353553390593274   0.353553390593274;
       0.277785116509801  -0.490392640201615   0.097545161008064   0.415734806151273  -0.415734806151273  -0.097545161008064   0.490392640201615  -0.277785116509801;
       0.191341716182545  -0.461939766255643   0.461939766255643  -0.191341716182545  -0.191341716182545   0.461939766255643  -0.461939766255643   0.191341716182545;
       0.097545161008064  -0.277785116509801   0.415734806151273  -0.490392640201615   0.490392640201615  -0.415734806151273   0.277785116509801  -0.097545161008064];
elseif (N == 8) & strcmp(transform_type, 'dst')==1
    Tforward = [ 0.161229841765317   0.303012985114696   0.408248290463863   0.464242826880013   0.464242826880013   0.408248290463863   0.303012985114696   0.161229841765317;
       0.303012985114696   0.464242826880013   0.408248290463863   0.161229841765317  -0.161229841765317  -0.408248290463863  -0.464242826880013  -0.303012985114696;
       0.408248290463863   0.408248290463863                   0  -0.408248290463863  -0.408248290463863                   0   0.408248290463863   0.408248290463863;
       0.464242826880013   0.161229841765317  -0.408248290463863  -0.303012985114696   0.303012985114696   0.408248290463863  -0.161229841765317  -0.464242826880013;
       0.464242826880013  -0.161229841765317  -0.408248290463863   0.303012985114696   0.303012985114696  -0.408248290463863  -0.161229841765317   0.464242826880013;
       0.408248290463863  -0.408248290463863                   0   0.408248290463863  -0.408248290463863                   0   0.408248290463863  -0.408248290463863;
       0.303012985114696  -0.464242826880013   0.408248290463863  -0.161229841765317  -0.161229841765317   0.408248290463863  -0.464242826880013   0.303012985114696;
       0.161229841765317  -0.303012985114696   0.408248290463863  -0.464242826880013   0.464242826880013  -0.408248290463863   0.303012985114696  -0.161229841765317];
elseif strcmp(transform_type, 'dct') == 1,
    Tforward    = dct(eye(N));
elseif strcmp(transform_type, 'dst') == 1,
    Tforward    = dst(eye(N));
elseif strcmp(transform_type, 'DCrand') == 1,
    x = randn(N); x(1:end,1) = 1; [Q,R] = qr(x); 
    if (Q(1) < 0), 
        Q = -Q; 
    end;
    Tforward = Q';
else 
    
    dwtmode('per','nodisp');  
    
    Tforward = zeros(N,N);
    for i = 1:N
        Tforward(:,i)=wavedec(circshift([1 zeros(1,N-1)],[dec_levels i-1]), log2(N), transform_type);  %% construct transform matrix
    end
end

Tforward = (Tforward' * diag(sqrt(1./sum(Tforward.^2,2))))'; 

Tinverse = inv(Tforward);

return;

function [y, A, l2normLumChrom]=function_rgb2LumChrom(xRGB, colormode)

if nargin==1
    colormode='opp';
end
change_output=0;
if size(colormode)==[3 3]
    A=colormode;
    l2normLumChrom=sqrt(sum(A.^2,2));
else
    if strcmp(colormode,'opp')
        A=[1/3 1/3 1/3; 0.5  0  -0.5; 0.25  -0.5  0.25];
    end
    if strcmp(colormode,'yCbCr')
        A=[0.299   0.587   0.114;   -0.16873660714285  -0.33126339285715   0.5;   0.5  -0.4186875  -0.0813125];
    end
    if strcmp(colormode,'pca')
        A=princomp(reshape(xRGB,[size(xRGB,1)*size(xRGB,2) 3]))';
        A=A./repmat(sum(A.*(A>0),2)-sum(A.*(A<0),2),[1 3]);
    else
        if nargout==2
            change_output=1;
        end
    end
end

maxV = sum(A.*(A>0),2);
minV = sum(A.*(A<0),2);
yNormal = (reshape(xRGB,[size(xRGB,1)*size(xRGB,2) 3]) * A' - repmat(minV, [1 size(xRGB,1)*size(xRGB,2)])') * diag(1./(maxV-minV)); % put in range [0,1]
y = reshape(yNormal, [size(xRGB,1) size(xRGB,2) 3]);

l2normLumChrom = diag(1./(maxV-minV))*sqrt(sum(A.^2,2));

if change_output
    A=l2normLumChrom;
end

return;

function yRGB=function_LumChrom2rgb(x,colormode)

if nargin==1
    colormode='opp';
end
if size(colormode)==[3 3]
    A=colormode;
    B=inv(A);
else
    if strcmp(colormode,'opp')
        A =[1/3 1/3 1/3; 0.5  0  -0.5; 0.25  -0.5  0.25];
        B =[1 1 2/3;1 0 -4/3;1 -1 2/3];
    end
    if strcmp(colormode,'yCbCr')
        A=[0.299   0.587   0.114;   -0.16873660714285  -0.33126339285715   0.5;   0.5  -0.4186875  -0.0813125];
        B=inv(A);
    end
end

maxV = sum(A.*(A>0),2);
minV = sum(A.*(A<0),2);
xNormal = reshape(x,[size(x,1)*size(x,2) 3]) * diag(maxV-minV) +  repmat(minV, [1 size(x,1)*size(x,2)])'; % put in range [0,1]
yRGB = reshape(xNormal * B', [ size(x,1) size(x,2) 3]);

return;

