f_samp = 400e3;
delta = 0.15;   %Tolerance
%Filter Band Edge speifications
fs1 = 67e3;
fp1 = 70e3;
fp2 = 95e3;
fs2 = 98e3;
%FIR Bandstop Kaiser Window paramters
A = -20*log10(delta);
beta = 0;  %As A is less than 21 
Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-7.95)/(2.285*0.015*pi));       %empirical formula for N_min

%Window length for FIR Bandstop Kaiser Window
n=N_min + 25;
%Ideal bandstop impulse response of length "n"
bs_ideal =  ideal_lp(pi,n) -ideal_lp(0.4825*pi,n) + ideal_lp(0.3425*pi,n);
%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';
FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response
%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
plot(f,abs(H))
grid
