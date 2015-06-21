clc
close all
clear all

% Basic Carrier Aggregation for simple OFDM transmit chain
%OFDM 2048 point 
%The available bandwidth is 8 MHz

%Input Parameters-1 
Tu=224e-6; %useful OFDM symbol period
T=Tu/2048; %baseband elementary period
G=0; %choice of 1/4, 1/8, 1/16, and 1/32
delta=G*Tu; %guard band duration
Ts=delta+Tu; %total OFDM symbol period
Kmax=1705; %number of subcarriers
Kmin=0;
FS=4096; %IFFT/FFT length
q=10; %carrier period to elementary period ratio
fc=q*1/T; %carrier frequency
Rs=4*fc; %simulation period
t=0:1/Rs:Tu;

%Data generator-1 (A)
M=Kmax+1;
rand('state',0);
a=-1+2*round(rand(M,1)).'+i*(-1+2*round(rand(M,1))).';
A=length(a);
info=zeros(FS,1);
info(1:(A/2)) = [ a(1:(A/2)).']; %Zero padding
info((FS-((A/2)-1)):FS) = [ a(((A/2)+1):A).'];
%Subcarriers generation (B)
carriers=FS.*ifft(info,FS);
tt=0:T/2:Tu;
figure(1);
subplot(211);
stem(tt(1:20),real(carriers(1:20)));

subplot(212);
stem(tt(1:20),imag(carriers(1:20)));
figure(2);
f=(2/T)*(1:(FS))/(FS);
subplot(211);
plot(f,abs(fft(carriers,FS))/FS);
subplot(212);
pwelch(carriers,[],[],[],2/T);


% D/A simulation-1(A)
L = length(carriers);
chips = [ carriers.';zeros((2*q)-1,L)];
p=1/Rs:1/Rs:T/2;
g=ones(length(p),1); %pulse shape
figure(3);
stem(p,g);
dummy=conv(g,chips(:));
u=[dummy(1:length(t))]; % (C)
figure(4);
subplot(211);
plot(t(1:400),real(u(1:400)));
subplot(212);
plot(t(1:400),imag(u(1:400)));
figure(5);
ff=(Rs)*(1:(q*FS))/(q*FS);
subplot(211);
plot(ff,abs(fft(u,q*FS))/FS);
subplot(212);
pwelch(u,[],[],[],Rs);
[b,a] = butter(13,1/20); %reconstruction filter
[H,F] = freqz(b,a,FS,Rs);
%[H,F] = FREQZ(b,a,FS,Rs);
figure(6);
plot(F,20*log10(abs(H)));
uoft = filter(b,a,u); %baseband signal (D)
figure(7);
subplot(211);
plot(t(80:480),real(uoft(80:480)));
subplot(212);
plot(t(80:480),imag(uoft(80:480)));
figure(8);
subplot(211);
plot(ff,abs(fft(uoft,q*FS))/FS);
subplot(212);
pwelch(uoft,[],[],[],Rs);

%Upconverter-1(A) to 91MHz
s_tilde=(uoft.').*exp(1i*2*pi*fc*t);
s=real(s_tilde); %passband signal (E)
figure(9);
plot(t(80:480),s(80:480));
figure(10);
subplot(211);

%plot(ff,abs(fft(((real(uoft).').*cos(2*pi*fc*t)),q*FS))/FS);
%plot(ff,abs(fft(((imag(uoft).').*sin(2*pi*fc*t)),q*FS))/FS);

plot(ff,abs(fft(s,q*FS))/FS);
subplot(212);
%hold on;
%pwelch(((real(uoft).').*cos(2*pi*fc*t)),[],[],[],Rs);
%pwelch(((imag(uoft).').*sin(2*pi*fc*t)),[],[],[],Rs);
pwelch(s,[],[],[],Rs);


% Input parameters-2(B)
Tu1=200e-6; %useful OFDM symbol period
T1=Tu1/2048; %baseband elementary period
G1=0; %choice of 1/4, 1/8, 1/16, and 1/32
delta1=G1*Tu1; %guard band duration
Ts1=delta1+Tu1; %total OFDM symbol period
Kmax1=1705; %number of subcarriers
Kmin1=0;
FS1=4096; %IFFT/FFT length
q1=10; %carrier period to elementary period ratio
fc1=q1*1/T1; %carrier frequency
Rs1=4*fc1; %simulation period
t1=0:1/Rs1:Tu1;

%Data generator-2 (B)
M1=Kmax1+1;
rand('state',0);
a1=-1+2*round(rand(M1,1)).'+i*(-1+2*round(rand(M1,1))).';
A1=length(a1);
info=zeros(FS1,1);
info(1:(A1/2)) = [ a1(1:(A1/2)).']; %Zero padding
info((FS1-((A1/2)-1)):FS1) = [ a1(((A1/2)+1):A1).'];
%Subcarriers generation (B)
carriers1=FS1.*ifft(info,FS1);
tt1=0:T1/2:Tu1;
figure(11);
subplot(211);
stem(tt1(1:20),real(carriers1(1:20)));

subplot(212);
stem(tt1(1:20),imag(carriers1(1:20)));
figure(12);
f1=(2/T1)*(1:(FS1))/(FS1);
subplot(211);
plot(f1,abs(fft(carriers1,FS1))/FS1);
subplot(212);
pwelch(carriers1,[],[],[],2/T1);


% D/A simulation-2 (B)
L1 = length(carriers1);
chips1 = [ carriers1.';zeros((2*q1)-1,L1)];
p1=1/Rs1:1/Rs1:T1/2;
g1=ones(length(p1),1); %pulse shape
figure(13);
stem(p1,g1);
dummy1=conv(g1,chips(:));
u1=[dummy1(1:length(t1))]; % (C)
figure(14);
subplot(211);
plot(t1(1:400),real(u1(1:400)));
subplot(212);
plot(t1(1:400),imag(u1(1:400)));
figure(15);
ff1=(Rs1)*(1:(q1*FS1))/(q*FS1);
subplot(211);
plot(ff1,abs(fft(u1,q1*FS1))/FS1);
subplot(212);
pwelch(u1,[],[],[],Rs1);
[b1,a1] = butter(13,1/20); %reconstruction filter
[H1,F1] = freqz(b1,a1,FS1,Rs1);
%[H,F] = FREQZ(b,a,FS,Rs);
figure(16);
plot(F1,20*log10(abs(H1)));
uoft1 = filter(b1,a1,u1); %baseband signal (D)
figure(17);
subplot(211);
plot(t1(80:480),real(uoft1(80:480)));
subplot(212);
plot(t1(80:480),imag(uoft1(80:480)));
figure(18);
subplot(211);
plot(ff1,abs(fft(uoft1,q1*FS1))/FS1);
subplot(212);
pwelch(uoft1,[],[],[],Rs1);

%Upconverter-2 (B) to about 110 MHz 
s_tilde1=(uoft1.').*exp(1i*2*pi*fc1*t1);
s1=real(s_tilde1); %passband signal (E)
figure(19);
plot(t1(80:480),s1(80:480));
figure(20);
subplot(211);

%plot(ff,abs(fft(((real(uoft).').*cos(2*pi*fc*t)),q*FS))/FS);
%plot(ff,abs(fft(((imag(uoft).').*sin(2*pi*fc*t)),q*FS))/FS);

%hold on 
plot(ff1,abs(fft(s1,q1*FS1))/FS1);
subplot(212);
%pwelch(((real(uoft).').*cos(2*pi*fc*t)),[],[],[],Rs);
%pwelch(((imag(uoft).').*sin(2*pi*fc*t)),[],[],[],Rs);
pwelch(s1,[],[],[],Rs1);

figure(21);
pwelch(s,[],[],[],Rs)+pwelch(s1,[],[],[],Rs1);
