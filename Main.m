clear;

N   = 20;    % Size of matrix K
C   = 0.01;  % threshold for determinant
n_t = 100;   % number of theta points
z0  = -4;    % Spectral shift
M   = 400;   % Precision in complex plane

t  = linspace(0,2*pi,n_t); % theta array
Id = speye(2*N+1);

% Define lattice in C:
L = build_lattice(0-6i, 25+6i, M);

%  Define potential matrix:
n = 1e+6;
a = fourier_coefficient((-2*N-1):(2*N+1), n);
V = compute_potential_matrix(a, N);

%  Define the vector of k's containing integer multiples of 2π:
k = -N:N;
k = 2*pi*k;

% Inverse square root of H0:
sqrt_H0_inv = spdiags(1./sqrt(z0 + k.^2).', 0, 2*N+1, 2*N+1);

Spectrum = zeros(size(L));
%% Main loop:
disp('Main Loop...')
tic
for m=1:length(t)
    theta=t(m);
    Determinant = zeros(size(L));
%     Constant term |θ|^2-1:
    qm = (theta^2 - z0) * Id;
%     First part of potential term:  H_0^(-1/2)(|θ|^2−1+V):
    A = sqrt_H0_inv * (qm+V);
    parfor n=1:length(L(:))
        z = L(n);
%         Gradient term -2iθ∇:
        diag_grad = 2*theta*k./(z - z0 - k.^2);
        grad = spdiags(diag_grad.', 0, 2*N+1, 2*N+1);
%         Resolvent:
        diag_resolvent = sqrt(z0 + k.^2)./(z - z0 - k.^2);
        sqrt_H0_resolvent = spdiags(diag_resolvent.', 0, 2*N+1, 2*N+1);
%         Compact operator K(z,θ):
        K = grad + A * sqrt_H0_resolvent;
        Determinant(n) = det(Id - K); 
    end
    Spectrum_theta = abs(Determinant)<C;
    Spectrum = or(Spectrum, Spectrum_theta);
end
disp('Done!')
disp([num2str(toc/60),' ',' minutes'])

%% ========================================================================
% Plot results:

x = linspace(0,1,1000);
figure('Position',[100 100 1100 600])

subplot(2,1,1)
plot(x, [real(potential(x)); imag(potential(x))]);

subplot(2,1,2)
plot(L(Spectrum), '.')
xlim([min(real(L(:))) max(real(L(:)))]);
ylim([min(imag(L(:))) max(imag(L(:)))]);











