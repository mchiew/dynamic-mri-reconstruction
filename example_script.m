% Example script of some low-rank k-t MR image reconstruction with
% sensitivity encoding and non-Cartesian sampling
% Mark Chiew
% mark.chiew@ndcn.ox.ac.uk
%
% Compiled 24/07/19
%
% These examples illustrate a few different methods for reconstructing
% under-sampled k-t data, including:
%   - un-regularised reconstruction (linear)
%   - temporal finite difference L2 regularisation (linear)
%     see Chiew & Miller, NeuroImage 2019 (https://doi.org/10.1016/j.neuroimage.2019.116165)
%   - spatial finite difference L2 regularisation (linear)
%   - spatio-temporal finite difference L2 regularisation (linear)
%   - low-rank recon using soft thresholding (nonlinear)
%   - low-rank recon using hard thresholding (nonlinear)
%     see Chiew et al., MRM 2015 (https://doi.org/10.1002/mrm.25395) 
%     and Chiew et al., MRM 2016 (https://doi.org/10.1002/mrm.26079)
%
% NB: This example relies on MR encoding operators (for the sensitivity and non-Cartesian 
% k-space encoding) that can be downloaded at https://users.fmrib.ox.ac.uk/~mchiew/tools.html

%% xfm_NUFFT with random trajectory, coil sensitivities and multiple time points
k = randn(1024,64,2);   % Gaussian density 2D random with 64 time-points

load('sens');   % load coil sensitivities

% Strongly recommend using NUFFT-based operator (relies on J. Fessler's IRT
% toolbox (http://web.eecs.umich.edu/~fessler/code/index.html)
E = xfm_NUFFT([64,64,1,64],sens,[],k,'wi',1); % no density compensation

% If not available, the following DFT-based transform has no dependencies,
% but is slow and takes ~ 20 GB RAM to initialise
% E = xfm_DFT([64,64,1,64],sens,[],k); % no density compensation


%   Dynamic phantom with some different time-courses in ROIs
[~,P] = phantom(64);
for t = 1:64
    P(5,1)  =   0.1*sin(2*pi*t/10 + 0*pi/4);
    P(6,1)  =   0.1*sin(2*pi*t/20 + 1*pi/4);
    P(7,1)  =   0.1*sin(2*pi*t/30 + 2*pi/4);
    P(8,1)  =   0.1*sin(2*pi*t/40 + 3*pi/4);
    P(9,1)  =   0.1*sin(2*pi*t/50 + 4*pi/4);
    P(10,1) =   0.1*sin(2*pi*t/60 + 5*pi/4);
    x(:,:,1,t)  =   phantom(P,64);
end

y = E*x;            % Forward Transform

%% Linear Recon, no regularisation
z1 = reshape(E.iter(y,@pcg,1E-4,100),[E.Nd E.Nt]);
fprintf(1, 'Linear Recon RMSE: %f\n', norm(z1(:)-x(:))/norm(x(:)));

%% Linear Recon, with temporal finite difference regularisation
z2 = reshape(E.iter(y,@pcg,1E-4,100,[0 0 0 1E4]),[E.Nd E.Nt]);
fprintf(1, 'Linear Recon w/ temporal regularisation RMSE: %f\n', norm(z2(:)-x(:))/norm(x(:)));

%% Linear Recon, with spatial finite difference regularisation
z3 = reshape(E.iter(y,@pcg,1E-4,100,[1E0 1E0 0 0]),[E.Nd E.Nt]);
fprintf(1, 'Linear Recon w/ spatial regularisation RMSE: %f\n', norm(z3(:)-x(:))/norm(x(:)));

%% Linear Recon, with spatio-temporal finite difference regularisation
z4 = reshape(E.iter(y,@pcg,1E-4,100,[1E0 1E0 0 1E4]),[E.Nd E.Nt]);
fprintf(1, 'Linear Recon w/ spatial-temporal regularisation RMSE: %f\n', norm(z4(:)-x(:))/norm(x(:)));

%% Non-linear low-rank recon using singular value thresholding
z5 = reshape(svt(E, y, 1E3, 'step', 5E-5, 'tol', 1E-3),[E.Nd E.Nt]);
fprintf(1, 'Non-Linear Low Rank (SVT) RMSE: %f\n', norm(z5(:)-x(:))/norm(x(:)));

%% Non-linear low-rank recon using iterative hard thresholding with shrinkage
z6 = iht_ms(E, y, 'rank', 8, 'step', 1E-4, 'shrink', 0.5, 'tol', 1E-3);
z6 = reshape(z6.u*z6.v',[E.Nd E.Nt]);
fprintf(1, 'Non-Linear Low Rank (IHT+MS) RMSE: %f\n', norm(z6(:)-x(:))/norm(x(:)));

%% Example reconstructed image
figure();
subplot(2,3,1); imshow(abs(z1(:,:,1,32)),[0 1]); title('Linear');
subplot(2,3,2); imshow(abs(z2(:,:,1,32)),[0 1]); title('T-reg');
subplot(2,3,3); imshow(abs(z3(:,:,1,32)),[0 1]); title('X-reg');
subplot(2,3,4); imshow(abs(z4(:,:,1,32)),[0 1]); title('XT-reg');
subplot(2,3,5); imshow(abs(z5(:,:,1,32)),[0 1]); title('SVT');
subplot(2,3,6); imshow(abs(z6(:,:,1,32)),[0 1]); title('IHTMS');

%% Example reconstructed time-course
figure();
subplot(2,3,1); plot(squeeze(abs(z1(35,33,1,:)))); title('Linear');
subplot(2,3,2); plot(squeeze(abs(z2(35,33,1,:)))); title('T-reg');
subplot(2,3,3); plot(squeeze(abs(z3(35,33,1,:)))); title('X-reg');
subplot(2,3,4); plot(squeeze(abs(z4(35,33,1,:)))); title('XT-reg');
subplot(2,3,5); plot(squeeze(abs(z5(35,33,1,:)))); title('SVT');
subplot(2,3,6); plot(squeeze(abs(z6(35,33,1,:)))); title('IHTMS');