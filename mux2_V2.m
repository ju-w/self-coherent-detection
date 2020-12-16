% version 2 of optical multiplexer / demultiplex /
% super Gauss filter (even symmetry)
% the center wavelength (equiv. lowpass range) is 0

function [yt, Ndelay] = mux2_V2(B3, fc, N, t0, xt);

% input paramters %%%%%%%%%%%%%%%%%%%%%%%%%%
% B3: optical 3 dB bandwidth
% N: filter order 1, 2 or 3
% t0: sampling interval
% xt:input sigal (field strength)
% fc: center wavelength

% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ndelay: delay in number of samples

% estimation of impulse response length (required, since fft is used)

if N == 1
  N_zp = 2*ceil(1/(t0*B3)); % zero padding (at beginning and end)
else
  if N == 2
    N_zp = 5*ceil(1/(t0*B3)); % zero padding (at beginning and end)
  else
    if N == 3
      N_zp = 10*ceil(1/(t0*B3)); % zero padding (at beginning and end)
    end
  end
end

xt = [zeros(N_zp, 1); xt; zeros(N_zp, 1)];

% fft
Nfft = length(xt);
if mod(Nfft, 2)
  Nfft = Nfft+1;
  xt = [xt; 0];
end

tp = Nfft*t0;
f0 = 1/tp; % frequ. resolution
f = f0*[-Nfft/2:Nfft/2-1]';

% tranfer function
Gf = fftshift( exp( -log(2) * (2*(f-fc)/B3).^(2*N) / 2) );
yt = ifft( fft(xt).*Gf );

Ndelay = N_zp;
