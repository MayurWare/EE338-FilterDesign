%Butterworth Analog LPF parameters
Wc = 1.03;              %cut-off frequency
N = 20;                  %order 
%20 LHP Poles of Butterworth Polynomial  
p1 = -1.02682 - 0.0808129*i;
p2 = -1.02682 + 0.0808129*i;
p3 = -1.00154 - 0.240449*i;
p4 = -1.00154 + 0.240449*i;
p5 = -0.951596 - 0.394164*i;
p6 = -0.951596 + 0.394164*i;
p7 = -0.878219 - 0.538174*i;
p8 = -0.878219 + 0.538174*i;
p9 = -0.783218 - 0.668931*i;
p10 = -0.783218 + 0.668931*i;
p11 = -0.668931 - 0.783218*i;
p12 = -0.668931 + 0.783218*i;
p13 = -0.538174 - 0.878219*i;
p14 = -0.538174 + 0.878219*i;
p15 = -0.394164 - 0.951596*i;
p16 = -0.394164 + 0.951596*i;
p17 = -0.240449 + 1.00154*i;
p18 = -0.240449 - 1.00154*i;
p19 = -0.0808129 - 1.02682*i;
p20 = -0.0808129 + 1.02682*i;
%plot([p1,p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20], *)
%Bandpass IIR Band Edge speifications
fp1 = 94;
fs1 = 97;
fs2 = 142;
fp2 = 145;
%Bilinear Transformation specifications
f_samp = 540;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
%Bandstop Transformation parameters
W0 = sqrt(wp1*wp2);
B = wp2-wp1;
[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20],Wc^N);   
%TF with poles p1-p20 and numerator Wc^N and no zeroes numerator chosen to make the DC Gain = 1
%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));        %bandstop transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));              %bilinear transformation
%Analog BPF Parameters
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification to collect coeffs
ns = ns/10^100;
ds = ds/10^100;
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;
%Discrete BPF Parameters
[nz, dz] = numden(discrete_bpf(z));                     %numerical simplification to collect coeffs   
nz = nz/10^50;
dz = dz/10^50;
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
%discrete_bpf(z) = vpa(simplify(vpa(expand(discrete_bpf(z)), 20)), 20);
fvtool(nz,dz)                                        %frequency response
%Magnitude Plot 
[H,f] = freqz(nz,dz,1024*1024, 540e3);
plot(f,abs(H))
grid
%Pole-Zero Plot
z = roots(nz);
p = roots(dz);
zplane(z,p);
grid
title('Pole Zero Plot')