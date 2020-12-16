function G_EST = channel_est2(preamble, N, L, Nup)

% channel estimation based on Golay sequences
% Nup: upsampling factor ; 2 for factor 2 fractional sampling
% Xmu = Rx-Spektren of length 2*(N+L)*Nup (received preambel incl. CP-intervals

N_CS = N;          % lengths of complementary sequences, needs to be power of 2
x1 = preamble(Nup*L+1:(L+N)*Nup);
x2 = preamble(Nup*(N+2*L)+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterative generation of Golay sequences
P = [1 1 1 -1].';
Q = [1 1 -1 1].';

for k = 3:log2(N_CS);
  P_old = P;
  P = [P_old; Q];
  Q = [P_old; -Q];
end

% these are the Galoy sequences of length N_CS
%G_A = (-P+1)/2;
%G_B = (Q+1)/2;
G_A = P;
G_B = Q;

X_A = conj( fft( upsample(G_A, Nup) ) );
X_B = conj( fft( upsample(G_B, Nup) ) );

% size(x2)
% size(X_B)

G_EST = (fft(x1).*X_A + fft(x2).*X_B)/(2*N);
