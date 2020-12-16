% Root-raised-cosine filter impulse response
% taken from https://en.wikipedia.org/wiki/Root-raised-cosine_filter
function y = rrc(t, beta, T_s)

ttt = T_s/4/beta;
y = 1/T_s * ( sin(pi*t/T_s*(1-beta)) + 4*beta*t/T_s .* cos(pi*t/T_s*(1+beta)) ) ...
         ./ ( pi * t/T_s .* (1 - (4*beta*t/T_s).^2) );
y(abs(t)<eps) = 1/T_s * (1 + beta*(4/pi-1));
ii = abs(t-ttt)<eps | abs(t+ttt)<eps;
y(ii) = beta / (T_s * sqrt(2)) * ( (1+2/pi)*sin(pi/4/beta) + (1-2/pi)*cos(pi/4/beta) );
