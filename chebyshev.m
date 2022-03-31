%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 5;
%5 LHP Poles Chebyshev Polynomials
p1 = -0.25388;
p2 = -0.078453 + 0.98123i;
p3 = -0.078453 - 0.98123i;
p4 = -0.20539 + 0.60643i;
p5 = -0.20539 - 0.60643i;
%Transfer function evaluating the Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
n3=  [1 -p5];
den = conv(n3,conv(n1,n2));          %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(6)];        % even order, DC Gain set as 1/(1+ epsilon^2)^0.5
%Bandstop IIR Band Edge speifications
fp1 = 67;
fs1 = 70;
fs2 = 95;
fp2 = 98;
%Bilinear Transformation Specifications
f_samp = 400;
wp1 = tan(fp1/f_samp*pi);          
ws1 = tan(fs1/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
%Bandpass Transformation Parameters
W0 = sqrt(wp1*wp2);
B = wp2-wp1;
%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bsf(s) = analog_lpf((B*s)/(s*s +W0*W0));     %bandpass transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));          %bilinear transformation
%Analog BSF Parameters
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;
%Discrete BSF Parameters
[nz, dz] = numden(discrete_bsf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       %frequency response in dB
%Magnitude Plot 
[H,f] = freqz(nz,dz,1024*1024, 400e3);
plot(f,abs(H))
grid
%Plor-Zero Plot
z = roots(nz);
p = roots(dz);
zplane(z,p);
grid
title('Pole Zero Plot')