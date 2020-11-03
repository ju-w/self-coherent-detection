function y = preamble_gen(N, L, N_HP)

% N: length of Golay sequence, power of 2
% L: length of CP
% N_HP: number of 101010.. bits for AGC and highpass
% y : binary preamble

N_CS = N;           % lengths of complementary sequences, needs to be power of 2
G_HP = repmat([1;-1], N_HP/2, 1);  % 101010...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterative generation of Golay sequences
P=[1 1 1 -1].';
Q=[1 1 -1 1].';

for k=3:log2(N_CS);
  P_old = P;
  P = [P_old; Q];
  Q = [P_old; -Q];
end

% these are the Galoy sequences of length N_CS
%G_A = (-P+1)/2;
%G_B = (Q+1)/2; 
G_A = -P;
G_B = Q; 

% preamble
y = [G_HP; -G_A(end-L+1:end); G_A; -G_A(1:L); -G_B(end-L+1:end); G_B; -G_B(1:L)];
y = (y+1)/2; % unipolar

