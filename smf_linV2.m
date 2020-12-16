% linear model of single-mode fiber
% the ouput signal is calc. in the f-domain
% NOTE: if x(t) exhbits a DC-offset or steps, x(t) should be filtered first (mux)
function [yt, Ndelay] = smf_linV2(alpha_dB, L, D, S, t0, xt, Bsignal);

% note: lamba is assumed to be 1.55 mum
% alpha_dB;     % loss in dB/km
% L: length in km
% D: 1st order dispersion coefficient in ps/(km*nm); typically 17 @1.55mu
% S: 2nd order dispersion coefficient ps/(km*nm^2); z.b. 0.056 @1.55mu
% t0: sampling time interval
% xt: input signal (field strength)
% Bsignal: signal bandwidth (or 1/Tsym); required to estimate the zero padding interval
% yt: output signal
% Ndelay: signal delay in samples

% first and second order dispersion coefficients beta2 and beta3 derived from D and S
lambda = 1.55;             % in mum
c = 3e5;                   % Lichtgeschwindigkeit in km/s
beta2 = -(D*1e-6)/(2*pi*c)*(lambda*1e-6)^2;
beta3 = (S*1e3 + 2*D*1e-6/(lambda*1e-6)) ...
         * (lambda*1e-6)^4/(2*pi*c*1000)^2 * 1000;


% estimation of group delay difference within Bsignal
% used to estimate pulse spreading and thus N_zp
fp = 1/t0;
DeltaT = abs(beta2)*L*(2*pi)^2*Bsignal;
N_zp = 2*ceil(DeltaT/t0); % zero padding interval

xt = [zeros(N_zp, 1); xt; zeros(N_zp, 1)];

% fft
N = length(xt);
if mod(N, 2)
  N = N+1;
  xt = [xt; 0];
end

tp = N*t0;
f0 = 1/tp; % frequ. resolution
f = f0*[-N/2:N/2-1]';


% tranfer function
phi_f = -beta2*L*2*(pi*f).^2 - beta3*L*(2*pi*f).^3/6;
Gf = fftshift( exp( j*phi_f ) ); % only dispersion part
yt = ifft( fft(xt).*Gf ) * 10^(-alpha_dB/20*L);
Ndelay = N_zp;
