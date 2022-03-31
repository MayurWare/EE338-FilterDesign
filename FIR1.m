f_s = 540e3;
delta = 0.15;   %Tolerance
%Filter Band Edge speifications
fs1 = 94e3;
fp1 = 97e3;
fp2 = 142e3;
fs2 = 145e3;
Wc1 = fp1*2*pi/f_s;
Wc2  = fp2*2*pi/f_s;
%FIR Bandpass Filter Kaiser Window paramters
A = -20*log10(delta);
beta = 0;  %As A is less than 21 
N_min = ceil((A-7.95) / (2.285*0.013*pi));
%FIR Bandpass Filter 
n=N_min+43;
%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(0.531*pi,n) - ideal_lp(0.3535*pi,n);
%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';
FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response
%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, f_s);
plot(f,abs(H))
grid