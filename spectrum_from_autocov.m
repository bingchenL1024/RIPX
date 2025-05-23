function [f, S, R, lags] = spectrum_from_autocov(t, x)
% spectrum_from_autocov  Compute PSD via FFT of the autocovariance
%
%   [f, S, R, lags] = spectrum_from_autocov(t, x)
%
%   Inputs:
%     t    = time vector (Nx1 or 1xN), assumed uniformly spaced
%     x    = signal vector (same size as t)
%
%   Outputs:
%     f    = frequency vector (Hz), length = length(R)
%     S    = two‐sided PSD estimate, same length as f
%     R    = biased autocovariance, length = 2*N-1
%     lags = lag‐index vector (in samples), same length as R
%
%   The normalization is such that
%     mean(x.^2) ≈ trapz(f, S)
%
% -----------------------------------------------------------------------

% ensure column vectors
t = t(:);
x = x(:);
N = numel(x);

% sampling parameters
dt = mean(diff(t));      % assume uniform sampling
Fs = 1/dt;               % sampling frequency
T  = N*dt;               % total record length

% 1) remove mean → zero‐mean process
x0 = x - mean(x);

% 2) compute biased autocovariance, lags from -(N-1):(N-1)
[R, lags] = xcov(x0, 'biased');  

% 3) FFT of autocovariance → two‐sided PSD
%    - multiply by dt to approximate the continuous FT integral
%    - no 1/N scaling here, because xcov('biased') already divides by N
S_raw = fft(R);          
S     = fftshift(S_raw) * dt;

% 4) build frequency axis (Hz) centered at zero
L = numel(R);
% frequencies: –Fs/2 : +Fs/2 in steps of 1/(L*dt)
f = ((-floor(L/2)):(ceil(L/2)-1))' / (L*dt);

% 5) (Optional) verify Parseval: average power in time vs. freq
P_time = mean(x0.^2);
P_freq = trapz(f, real(S));
fprintf('Time‐domain mean power = %.6g,  Freq‐domain integral = %.6g\n', ...
        P_time, P_freq);

end