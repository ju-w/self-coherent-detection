% Single Sideband
% Single carrier Frequency domain 
clear all
close all

SPEC = 1;
Nblocks = 2; % round( 1/PbRef / N );

M = 16;
L=100; % length in km -> kleiner Wert für Back-To-Back
Rb=40*1.07e9;
D=17;

BW_MUX = 2*43e9;

Rel_bias = 1; % relative Bias power (the power of the carrier is Rel_bias * mean optical power of modulated component)

Rcos=0.4; % cos-roll off factor; chose not smaller than 0.2 or change delay of filter gtTx
NdelayRRC=10; % don't change ..

% FFT parameters
N = 4*256; % block size in symbols incl. CP
Lcp = 4*32;  % length of CP 
PbRef=1e-3;

% channel estimation based on an m-sequence
N1010 = 200;
xk_train = preamble_gen2(N, Lcp, N1010); % length 2*(L+N)+N1010

% fgBesselTx = fp*0.75;            % electrical Bessel filter 3 dB cut-off
% fgBesselRx = fgBesselTx;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = BW_MUX*4; % simulation bandwidth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = log2(M)/Rb*(N-Lcp)/N; % QAM/ PSK symbol interval
Tb=1/Rb;
DeltaF=3/Ts;               % Diff. between center freq. of QAM-signal and carrier
fc_mux = (DeltaF + (1+Rcos)/(2*Ts))/2;  % center freq. MUX

Nover = ceil(Ts*fs); % samples per symbol interval
if mod(Nover,2);
  Nover=Nover+1;
end
t0 = Ts/Nover;
fs=1/t0; % correction to have integer number of symples per symbol interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transmitted data symbols with values 0..M-1
S2B = de2bi(0:M-1, log2(M) ); % mapping table: Bits bn to symbols zk
zk = randi(M, Nblocks*N, 1)-1; % symbols to be transmitted 
bk_tx = S2B( zk+1,: ); % bits to be transmitted


xk = reshape(qammod(zk, M, 'gray'), N, Nblocks); % each column contains a data block

xk_cp = [xk(end-Lcp+1:end,:) ;xk];   % Tx vector with CP (still: % each column is an OFDM-symbol)
xk_cp = reshape(xk_cp, Nblocks*(N+Lcp),1); % serial, digital transmit signal (electrical) before pulse shaping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PULSE SHPING and analog signal gerneration
%Tx-Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gtTx = rcosflt([1], 1, Nover, 'sqrt', Rcos, NdelayRRC); % NdelayRRC symbols delay
xt = conv(gtTx, upsample( [xk_train; xk_cp], Nover ) ); % complex time domain signal before upconversion for single sideband

t=[0:length(xt)-1]'*t0;
xt = xt.*exp(j*2*pi*t*DeltaF); % frequeny shifting to get a signle sideband signal

Bias = Rel_bias * mean( abs(xt).^2 );  % referenz-carrier (LO)

% optical signal (perfect modulator)
xt_opt = xt + sqrt(Bias);
Popt = mean( abs(xt_opt).^2 );
Eb = Popt*Tb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optical transmission
[yt_opt1, Ndelay1] = mux2_V2(BW_MUX, fc_mux, 2,  t0, xt_opt); % MUX direkt am Sender

% Glasfaser, lineares Modell: hier: ohne Daempfung, da nur SNR interessiert
if L>0
  [yt_opt2, Ndelay2] = smf_linV2(0.0, L, D, 0.0, t0, yt_opt1, 1/Tb);
else
  Ndelay2=0;
  yt_opt2 = yt_opt1;
end

% mittlere optische Leistung vor optischem Empfangsfilter (Rx-Eingang)

% after demux
[yt_opt3, Ndelay3] = mux2_V2(BW_MUX, fc_mux, 2,  t0, yt_opt2);  % Signal nach Demux ohne Rauschen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electrical domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% current after Rx-filter
it_Rx = abs(yt_opt3).^2; % elektrischer Strom
Ndelay = Ndelay1 + Ndelay2 + Ndelay3 + (N1010 - Lcp/2 + 2*NdelayRRC)*Nover;


%  here an anlog down-conversion is assumed
t=[0:length(it_Rx)-1]'*t0;
it_complex = it_Rx.*exp(-j*2*pi*t*DeltaF); % frequeny shifting to get a complex baseband signal

% Rx-filtering
gtRx = gtTx;
it_complex = conv( it_complex, gtRx );


% fractional sampling 
yk_tmp = it_complex(1+Ndelay+2*(Lcp+N)*Nover:Nover/2:end); % 2*(Lcp+N)*Nover due to Golay
yk = reshape( yk_tmp(1:(Nblocks)*(N+Lcp)*2), (N+Lcp)*2, Nblocks ); % only data part
yk = yk(2*Lcp+1:end,:); % remove cyclic prefix interval

Ymu = fft(yk); % spectrum

preamble =  it_complex(1+Ndelay:Nover/2:Ndelay+2*(Lcp+N)*Nover );
G_EST  = channel_est2(preamble, N, Lcp, 2);



% Emu for fraction sampling (equalizer coefficients)
idx=[1:N];
Norm=abs( G_EST(idx) ).^2 + abs( G_EST(idx+N) ).^2 ;
Emu = conj( G_EST(idx) )./ Norm;
idx=[N+1:2*N];
Emu(idx) = conj( G_EST(idx) ) ./ Norm;
Emu = 2*Emu;

%Emu=abs(Emu).*exp(-j*phase(Emu));


%%%%%%%%% end of channel  estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normiertes Signal nach Equalizer
Ymu = Ymu .* repmat(Emu, 1, Nblocks) ; % 
idx=[1:N];
YmuE = 0.5*(Ymu(idx,:) + Ymu(idx+N,:)); % Downsampling in f
sv_rx = ifft( YmuE );


scatterplot(sv_rx(1:end,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optisches Spektrum
Nwin = (N+L)*Nover; % bestimmt spektrale Auflösung diese ist:
w = rectwin( Nwin );
f0 = 1/t0/Nwin;
if SPEC
  figure(3)
  [Phi_xx, f] = pwelch( xt_opt  , w, 0, Nwin, 1/t0, 'twosided' );
  f=[-Nwin/2+1:Nwin/2]*f0;
  plot(f/1e9, db10(fftshift(Phi_xx/Phi_xx(2))), 'linewidth', 2);
  axis([-50 80 -30 70])
  SetFontSize
  title('');
  xlabel('x');
  ylabel('y');

  grid on
end
