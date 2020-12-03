function [gn, CP_start] = sync_gen(yn, N, L)

% yn: sampled Rx-signal (part which contains the preamble ...)
% N: length of Golay sequence, power of 2
% L: length of CP

% return values:
% gn: estimated impulse response
% CP_start: start of CP-intervall (data part)

N_CS = N;    % lengths of complementary sequences, needs to be power of 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterative generation of Golay sequences
P = [1 1 1 -1].';
Q = [1 1 -1 1].';

for k = 3:log2(N_CS);
  P_old = P;
  P = [P_old; Q];
  Q = [P_old; -Q];
end

N_yn = length(yn);

% 2 matched filters

gn_A = zeros(N_CS,1);
gn_B = zeros(N_CS,1);

gn_A(1:end) = P(end:-1:1);
gn_B(1:end) = Q(end:-1:1);

y_corr1 = conv(yn, gn_A);
y_corr2 = conv(yn, gn_B);

% cross correlation for timing
% sequence A needs to be delayes by N_CS bits, e.g., 2*N_CS samples
psi = -y_corr1(1:end-(N_CS+2*L)) + y_corr2(1+(N_CS+2*L):end);

% now we use a sliding window to find a CP-position with max. energy
gn_slid = ones(L, 1)/(L);
psi2 = conv(gn_slid, (psi).^2);

[tmp, indMax] = max(psi2);
start_optEst = indMax-L+1;
gn = psi(start_optEst:start_optEst+L-1) / (2*N_CS);
CP_start = start_optEst + L + (N_CS+2*L) + 1;
