function yt = bessel_filt(fg, f_sample, xt);
% 5th order Bessel low-pass filter
% PARAMETERS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fg: 3 dB cut-off
% fsample: sample frequency
% xt: input signal
% yt: filtered signal

t0 = 1/f_sample;
[alpha, beta] = besself(5, (2*pi)/0.6157); %fg=1 Hz
[tNorm, gtNorm] = Gp2gt(alpha, beta, 1/t0/fg, 5);
gtBessel = gtNorm'*fg;
yt = conv(xt, gtBessel)*t0;
