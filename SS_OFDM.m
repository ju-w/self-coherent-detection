% selfhomodyne detection of a single-sideband OFDM-signal
% the carrier (unmodulated) is at f=0;
% suffers from signal-signal-beat interference, but it is bandwidth efficient
% here we have to find the best BIAS (which minimized the required OSNR!)
pkg load communications

% TBD: 1.) Noise 2.) Bias  3.) dual drive MZMs %%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all;


Nsim = 50; % simulated OFDM-symbols


Rel_bias = 7; % relative Bias power (the power of the carrier is Rel_bias * mean optical power of modulated component)
NOISE = 1;

M=4;      % QAM modulation order
L=100; % length in km -> kleiner Wert für Back-To-Back
Rb=40*1.07e9;
D=17;
BW_MUX = 43e9;

SPEC = 1;

Nover=10; % determines 'analog' time redulation t0: t0=1/fp/Nover

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tb=1/Rb;
% OFDM Parameter
N   = 4*256; % FFT-Size
Lcp = 32;  % length of CP
Ncarrier = N/2-20; % number of indep. modulated subcarriers
Ntrain = 100;    % Ntrain zeros as preamble

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tofdm= log2(M)/Rb*Ncarrier;   % duration of an OFDM symbol in seconds
Tfft = N / (N+Lcp) * Tofdm;   % FFT-window size in seconds
f0 = 1/Tfft;                % spacing between subcarriers in Hz
fp = N*f0;                  % sampling frequency (at Tx DAC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = (1/fp)/Nover;          % time resolution on simulation level

fgBesselTx = fp*0.75;            % electrical Bessel filter 3 dB cut-off
fgBesselRx = fgBesselTx;

% center frequency of MUX /DEMUX
f_center = Ncarrier/2*f0;
if f_center > BW_MUX/2
  fprintf("Bandwidth too large");
end

% transmitted data symbols with values 0..M-1
S2B = de2bi(0:M-1, log2(M) ); % mapping table: Bits bn to symbols zk
Nsim= Nsim+1;                 % first symbol for channel estimation
zk = randi(M, Nsim*Ncarrier, 1)-1; % symbols to be transmitted
bk_tx = S2B(zk(Ncarrier+1:end)+1,: ); % bits to be transmitted; first OFDM-symbol
                                      % acts for channel estimation
sv_tx =  reshape(qammod(zk, M), Ncarrier, Nsim); % each column contains the OFDM-data
sv_tx(2:2:end,1)=0; % for channel estimation only odd carriers
Nbits = (Nsim-1)*Ncarrier*log2(M);

% electrical Signal in f-domain, each columns contains an OFDM-symbol
Xmu = zeros(N, Nsim);

% Single Sideband signal
mu=2:Ncarrier+1;
Xmu(mu,:) = sv_tx; % all the other subcarrier are set to zero -> single sideband

xk = ifft(Xmu)*N;  % complex time domain signal
xk = [xk(end-Lcp+1:end,:) ;xk];   % Tx vector with CP (still:
                                  % each column is an OFDM-symbol)
xk = reshape(xk, Nsim*(N+Lcp),1); % serial, digital transmit signal (electrical)

% assumed analog transmit signal (here: rectangular interpolation)
gt_tx = rect([-Nover/2:Nover/2]/Nover);
xt = conv(gt_tx, upsample([zeros(Ntrain,1); xk], Nover) );
xt = conv(xt, gt_bessel(fgBesselTx, 1/t0))*t0;

Bias = Rel_bias * sum(abs(xt).^2)/length(xt) ;  % referenz-carrier (LO)

% optical signal (perfect modulator)
xt_opt = xt + sqrt(Bias);

% mean optical power at sender
Popt = sum(abs(xt_opt).^2)/length(xt_opt);
Eb = Popt*Tb;

% only for later analysis
i_ac = abs(xt).^2;
i_misch = 2*real(sqrt(Bias)*conj(xt) );
P_ac = sum(i_ac.^2) / length(xt);
P_misch = sum(i_misch.^2) / length(xt);
P_misch/P_ac

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optical transmission
[yt_opt1, Ndelay1] = mux2_V2(BW_MUX, f_center, 2,  t0, xt_opt); % MUX direkt am Sender

% Glasfaser, lineares Modell: hier: ohne Daempfung, da nur SNR interessiert
if L>0
  [yt_opt2, Ndelay2] = smf_linV2(0.0, L, D, 0.0, t0, yt_opt1, 1/Tb);
else
  Ndelay2=0;
  yt_opt2 = yt_opt1;
end

% after demux
[yt_opt3, Ndelay3] = mux2_V2(BW_MUX, f_center, 2,  t0, yt_opt2);  % Signal nach Demux ohne Rauschen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RX
% current after Rx-filter
it_Rx = abs(yt_opt3).^2; % elektrischer Strom
it_Filt = bessel_filt(fgBesselRx, 1/t0, it_Rx ); % elektrischer Strom
Ndelay = Ndelay1 + Ndelay2 + Ndelay3 + Ntrain*Nover-Lcp/2*Nover + 2*Nover;

yk_tmp = it_Filt(1+Ndelay:Nover:end); % sampled Rx-signal
yk = reshape( yk_tmp(1:Nsim*(N+Lcp)), N+Lcp, Nsim );

yk = yk(Lcp+1:end,:); % remove cyclic prefix interval
Ymu = fft(yk)/N;

% channel estimation on odd subcarriers
Emu = zeros(Ncarrier,1);
mu = 1:2:Ncarrier;
Emu(mu) = sv_tx(mu,1) ./ Ymu(mu+1,1);
Emu(mu+1) = interp1(mu, Emu(mu), mu+1, 'linear','extrap' ); % linear interpol.

sv_rx = Ymu(2:Ncarrier+1,:) .* repmat(Emu, 1, Nsim) ;
zk_rx =  qamdemod(sv_rx, M); % each column contains an OFDM-symbol
zk_rx = reshape(zk_rx, Nsim*Ncarrier, 1);
bk_rx = S2B(zk_rx(Ncarrier+1:end) + 1,: );            % bits received

% BER
Pb = sum(sum(abs(bk_tx -bk_rx))) / Nbits

%stem(real(sv_tx(1:2:end,3:end))-real(sv_rx(:,3:end)))
scatterplot(sv_rx(1:end,2));

if NOISE;
printf("NOISE on!\n")

% normalized optical noise (noise power is 1)
Nyt_opt = length(yt_opt2);
nt_cp = (randn(Nyt_opt, 1) + j*randn(Nyt_opt, 1)); % co-polarized
nt_op = (randn(Nyt_opt, 1) + j*randn(Nyt_opt, 1)); % orthogonal-polarization

[nt_cpFilt] = muxV2(BW_MUX, 2,  t0, nt_cp); % coploarized noise (same polarization as signal)
[nt_opFilt, dummy] = muxV2(BW_MUX, 2,  t0, nt_op); % orthogonal noise

lf = 0;
for EbN0_dB = 10:1:35
   lf += 1;
   t1 = time();

   EbN0=10^(EbN0_dB/10);
   N0 = Eb/EbN0; % noise power spectral density per polarization

  % RX
  % current after Rx-filter
  it_Rx = abs(yt_opt3 + nt_cpFilt * sqrt(N0/(2*t0))).^2 ...
        + abs(          nt_opFilt * sqrt(N0/(2*t0))).^2; % elektrischer Strom

  it_Filt = bessel_filt(fgBesselRx, 1/t0, it_Rx ); % elektrischer Strom
  Ndelay = Ndelay1 + Ndelay2 + Ndelay3 + Ntrain*Nover-Lcp/2*Nover + 2*Nover;

  yk_tmp = it_Filt(1+Ndelay:Nover:end); % down-sampled Rx-signal
  yk = reshape( yk_tmp(1:Nsim*(N+Lcp)), N+Lcp, Nsim );

  yk = yk(Lcp+1:end,:); % remove cyclic prefix interval
  Ymu = fft(yk)/N;

  sv_rx = Ymu(2:Ncarrier+1,:) .* repmat(Emu, 1, Nsim);
  zk_rx = qamdemod(sv_rx, M); % each column contains an OFDM-symbol
  zk_rx = reshape(zk_rx, Nsim*Ncarrier, 1);
  bk_rx = S2B(zk_rx(Ncarrier+1:end) + 1,: ); % bits received

  % BER
  Pb(lf) = sum(sum(abs(bk_tx -bk_rx))) / Nbits;
  OSNR(lf) = 10*log10( Popt / (2*N0*12.5e9) );

  t2 = time();
  printf("Loop %d, t: %d, OSNR: %d, Pb: %d\n", lf, t2-t1, OSNR(lf), Pb(lf));
end

figure(1);
semilogy(OSNR, Pb,'linewidth',2);
axis([5 40 1e-3 1]);
xlabel('OSNR in dB');
ylabel('BER');

endif

% optisches Spektrum
Nwin = (N+L)*Nover; % bestimmt spektrale Auflösung diese ist:
w = rectwin( Nwin );
f0 = 1/t0/Nwin;
if SPEC
  figure(2)
  [Phi_xx, f] = pwelch( yt_opt3  , w, 0, Nwin, 1/t0, 'twosided' );
  f=[-Nwin/2+1:Nwin/2]*f0;
  plot(f/1e9, mag2db(fftshift(Phi_xx/Phi_xx(2))), 'linewidth', 2);
  axis([-80 80 -30 10])
  SetFontSize
  title('');
  xlabel('x');
  ylabel('y');

  grid on
end
