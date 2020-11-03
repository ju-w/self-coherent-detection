% for comparision with self coherent detection
% DC-Biased DMT (as for VLC-transmission)
% suffers from fading of the equivalent electrical channel
% and some signal-sigal-beat interference
% power loading is used; -> inefficient, if path loss to high!
% here channel estimation based on Golay sequence
pkg load communications

clear all
SPEC = 1;

Nsim = 200; % simulated OFDM-symbols

M = 16;      % QAM modulation order
L = 0; % length in km -> kleiner Wert für Back-To-Back
Rb=100*1.07e9;
kClip=2.5;  % determindes the Bias
D=17;
BW_MUX = 85e9; % for WDM and noise rejection at Rx

Nover=10; % determines 'analog' time redulation t0: t0=1/fp/Nover

% Sync Sequence Parameters
Ngolay=256; % length of Golay sequnce
Nhp=100;     % must be longer if a highpass is used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tb=1/Rb;
% OFDM Parameter
N = 512; % FFT-Size
Lcp = 32;  % length of CP
Ncarrier = N/2-30; % number of indep. modulated subcarriers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tofdm= log2(M)/Rb*Ncarrier; % duration of an OFDM symbol in seconds
Tfft = N / (N+Lcp) * Tofdm;   % FFT-window size in seconds
f0 = 1/Tfft;                % spacing between subcarriers in Hz
fp = N*f0;                  % sampling frequency (at Tx DAC)

f0DataMax = Ncarrier*f0/1e9
fgBesselTx=fp/2;
fgBesselRx=fp/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = (1/fp)/Nover;          % time resolution on analog level


% transmitted data symbols with values 0..M-1
S2B = de2bi(0:M-1, log2(M) ); % mapping table: Bits bn to symbols zk
zk = randi(M, Nsim*Ncarrier, 1)-1; % symbols to be transmitted
bk_tx = S2B(zk + 1,: );            % bits to be transmitted
sv_tx = reshape(qammod(zk, M), Ncarrier, Nsim); % each column contains the OFDM-data
Nbits = Nsim*Ncarrier*log2(M);

% electrical Signal in f-domain, each columns contains an OFDM-symbol
Xmu = zeros(N, Nsim);

mu=2:Ncarrier+1;
Xmu(mu,:) = sv_tx; % you may add power loading at this point ...

% for power loading of sub-carriers
lambda=1.55;             % in mum
c=3e5;                   % Lichtgeschwindigkeit in km/s
beta2 = -(D*1e-6)/(2*pi*c)*(lambda*1e-6)^2;
Gest = cos( (2*pi*(mu-1)'*f0).^2*L*beta2/2);
if L>0
  Xmu(mu,:) = Xmu(mu,:)./repmat( abs(Gest), 1, Nsim);
end


Xmu(end:-1:end-Ncarrier+1,:) = conj( Xmu(2:Ncarrier+1,:) ); % symmetry ensures that xk is real

xk = real( ifft(Xmu) )*N;
xk = [xk(end-Lcp+1:end,:) ;xk];   % Tx vector with CP (still: each column is an OFDM-symbol)
xk = reshape(xk, Nsim*(N+Lcp),1); % serial, digital transmit signal (electrical)

sigma_xk = sqrt( sum(xk.^2)/length(xk) );

% training sequence for sync and channel estimation
train = 2*preamble_gen(Ngolay, Lcp, Nhp)-1; % symbol spaced

% assumed analog transmit signal (here: rectangular interpolation)
gt_tx = rect([-Nover/2:Nover/2]/Nover);
xt = conv(gt_tx, upsample([sigma_xk*train; xk], Nover) );


% Clipping / required DC-bias and optical signal generation
Bias = kClip*sigma_xk;
xt_bias = xt + Bias; % 'analog' signal with bias (still bipolar)


xt_bias = conv(xt_bias, gt_bessel(fgBesselTx, 1/t0))*t0;
ind = find( xt_bias < 0 );
xt_bias(ind) = 0; % electrical unipolar 'analog' signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optical transmission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tx-signal
xt_opt = sqrt(xt_bias); % optical field strength; direct modulation of laser current

% mean optical power:
Popt = sum(abs(xt_opt).^2)/length(xt_opt);
Eb = Popt*Tb;

t=[0:length(xt_opt)-1]'*t0;
[yt_opt, dummy] = muxV2(2*BW_MUX, 2,  t0, xt_opt);    % reduce opt. BW in order to avoid steps at fiber inp.

if L > 0
  [yt_opt, dummy] = smf_linV2(0.0, L, D, 0.045, t0, yt_opt, 1/Tb );
end

% normalized optical noise (noise power is 1)
Nyt_opt=length(yt_opt);
nt_cp = (randn(Nyt_opt, 1) + j*randn(Nyt_opt, 1)); % co-polarized
nt_op = (randn(Nyt_opt, 1) + j*randn(Nyt_opt, 1)); % orthogonal-polarization

% optical filtering
[yt_opt_Filt, dummy] = muxV2(BW_MUX, 2,  t0, yt_opt);    % desired signal after optical Rx-filter

[nt_cpFilt] = muxV2(BW_MUX, 2,  t0, nt_cp); % coploarized noise (same polarization as signal)
[nt_opFilt, dummy] = muxV2(BW_MUX, 2,  t0, nt_op); % orthogonal noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electrical receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for now only desired signal for channel estimation and sync
% later channel estimation can also be based on the noise signal

% Bessel filtering (analog domain)
it_Filt = bessel_filt(fgBesselRx, 1/t0, abs(yt_opt_Filt).^2 );

% sampling with fp and random phase, yk_tmp includes the training sequence
delta =0; % random phase between -Nover/2 and Nover/2
ind1 = find( it_Filt > Bias/10 );  % just find 'signal' (like carrier detect in RF')
start1 = ind1(1)+delta;

yk_tmp = it_Filt(start1:Nover:end)-Bias;  % romove Bias, otherwise sync doesn't work well
                                          % in practice, a highpass is used to remove the Bias

% find start of data part and estimate impulse response
% gn: estimated impulse response of length L
[gk, start_idx] = sync_gen(yk_tmp(1:(Nhp+(Ngolay+Lcp))*2), Ngolay, Lcp);

yk = reshape( yk_tmp(start_idx:start_idx+Nsim*(N+Lcp)-1), N+Lcp, Nsim );  % data part only
yk = yk(Lcp+1:end,:); % remove cyclic prefix

gk = gk/sigma_xk; % since training part was weighted
Gmu = fft( [gk; zeros(N-length(gk), 1)] ); % estimated transfer function
Emu = 1./Gmu(2:Ncarrier+1);
if L>0
  Emu = Emu .* abs(Gest(1:Ncarrier)); % power loading
end

lf=0;
for EbN0_dB = 10:1:35
   lf=lf+1;
   EbN0=10^(EbN0_dB/10);
   N0 = Eb/EbN0; % noise power spectral density per polarization

   % filtered optical signal at PD input including noise
   yt_optRx = yt_opt_Filt + nt_cpFilt * sqrt(N0/(2*t0)); % signal plus copolarized noise
   ntRx_op = nt_opFilt * sqrt(N0/(2*t0)); % orthognal noise
   % photodiode output signal including noise
   it_Rx = bessel_filt(fgBesselRx, 1/t0, abs(yt_optRx).^2 + abs(ntRx_op).^2);

   yk_tmp = it_Rx(start1:Nover:end);  %xxx
   yk = reshape( yk_tmp(start_idx:start_idx+Nsim*(N+Lcp)-1), N+Lcp, Nsim );  % data part only
   yk = yk(Lcp+1:end,:); % remove cyclic prefix

   % Decision in f-domain
   Ymu = fft(yk)/N;
   sv_rx = Ymu(2:Ncarrier+1,:) .* repmat(Emu, 1, Nsim) ;
   zk_rx =  qamdemod(sv_rx, M); % each column contains an OFDM-symbol
   zk_rx = reshape(zk_rx, Nsim*Ncarrier, 1);

   bk_rx = S2B(zk_rx + 1,: );            % bits received

   Pb(lf) = sum(sum(abs(bk_tx -bk_rx))) / Nbits;
   OSNR(lf) = 10*log10( Popt / (2*N0*12.5e9) );

end

%figure(2);
%semilogy([10:0.1:22], Pb,'r','linewidth',2);
%xlabel('E_b/N_0 in dB');
%ylabel('BER');
%axis([5 15 1e-3 1]);

figure(1);
semilogy(OSNR, Pb,'linewidth',2);
axis([5 40 1e-3 1]);
xlabel('OSNR in dB');
ylabel('BER');


% Spektrum
Nwin = Nover*N; % bestimmt spektrale Auflösung diese ist:
f0_2 = 1/t0/Nwin;
w = rectwin( Nwin );
% optisches Spektrum am Tx
if SPEC
  figure(2)
  [Phi_xx, f] = pwelch( xt  , w, 0, Nwin, 1/t0, 'twosided' );
  %Phi_xx(1)=Phi_xx(2);
  f=[-Nwin/2+1:Nwin/2]*f0_2;
  plot(f/1e9, mag2db(fftshift(Phi_xx/Phi_xx(2))), 'linewidth', 2);
  axis([-75 75 -30 10])
  title('');
  xlabel('x');
  ylabel('y');

  grid on
end
