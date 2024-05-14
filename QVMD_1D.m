%This code is for algorithm QVMD, for extracting the signal components
%sequentially from the non-stationary signal mixture. The paper to describe 
%the algorithm is named "A Queued Variational Mode Decomposition Method". 
%It is finished by Chen Wei and ZhangYong in 2021 at JXUFE, Nanchang, China.
clear all
close all
ts = 0.0002;
fs = 1/ts;
N = 1/ts;

k = 1:N;
t = k*ts;

%y1=1./(1.2+cos(2*pi*t));
%y2=cos(32*pi*t+0.2*cos(64*pi*t))./(1.5+sin(2*pi*t)); 
%y1=sawtooth(t*2*pi*20);
%y2=0.2*sin(2*pi*50*t); no
y1=2*t;
y2=1*sin(2*pi*50*t-10*pi*t.^2);
y3=0.5*(exp(-5*(t-0.5).^2)).*sin(200*pi*t);
y4=randn(1,N)*0.1;
%y2=2*t;


%y3=0;y4=0;
y=y1+y2+y3+y4;
%y1=4*t+1;
%y2=(2+5*(t-0.5).^2).*sin(60*pi*t)*0.6;
%y3=1*sin(120*pi*t-2*t.^2);
%y4=randn(1,N)*0.2;

figure(1), plot(t,y,'r.-',t,y1,'b.',t,y2,'g.',t,y3,'m.',t,y4,'k.'),grid on;

 

%% elongation 1   %%%%%%%%
sp=polyfit(t,y,2);%%trend of the complete region
ytrend=polyval(sp,t);
y=y-ytrend;

LeL=round(N/10);
Np=3;
%elongate in left
xL=t(1:LeL);
yL=y(1:LeL);
%yLtrend=polyval(sp,xL);
%yL=yL-yLtrend;

fL=fft(yL,LeL);
f=(0:LeL-1)*fs/LeL;
[magL,pos]=sort(abs(fL(1:end/2))*2/LeL,'descend');

fre=f(pos(1:Np));
A=magL(1:Np);
ph=angle(fL(pos(1:Np)));

flag=0;
for i=1:Np
    if pos(i)==1
        A(i)=A(i)/2;
        flag=1;
        break;
    end
end

if flag==0
    A=[A,abs(fL(1))*2/LeL/2];
    fre=[fre,f(1)];
    ph=[ph,angle(fL(1))];
end


tL=[t(1)-LeL*ts:ts:t(1)-ts];
ytrendL=polyval(sp,tL);
%figure(2), plot(tL,ytrendL,'b.'),grid on;hold on;
yLe=ytrendL;
for i=1:length(A)
yLe=A(i)*cos(fre(i)*2*pi*tL+ph(i))+yLe;
end


clear  pos A fre ph f
xR=t(end-LeL+1:end);
yR=y(end-LeL+1:end);

%yRtrend=polyval(sp,xR);
%yR=yR-yRtrend;
fR=fft(yR,LeL);
f=(0:LeL-1)*fs/LeL;
[magR,pos]=sort(abs(fR(1:end/2))*2/LeL,'descend');

fre=f(pos(1:Np));
A=magR(1:Np);
ph=angle(fR(pos(1:Np)));

flag=0;
for i=1:Np
    if pos(i)==1
        A(i)=A(i)/2;
        flag=1;
        break;
    end
end

if flag==0
    A=[A,abs(fR(1))*2/LeL/2];
    fre=[fre,f(1)];
    ph=[ph,angle(fR(1))];
end


tR=[t(end)+ts:ts:t(end)+LeL*ts];
ytrendR=polyval(sp,tR);
yRe=ytrendR;
%figure(2), plot(tR,ytrendR,'r.-'),hold on
for i=1:length(A)
yRe=A(i)*cos(fre(i)*2*pi*tR+ph(i))+yRe;
end
y=y+ytrend;
ye=[yLe y yRe];
te=[tL t tR];
figure(2), plot(te,ye,'b.-',t,y,'r.-'),hold off,grid on;xlabel('Time/s'),ylabel('Amplitude');

yh = hilbert(ye);    % 
Nn=length(ye);
fyf=fft(yh,Nn);
fxx=[0:1/(Nn*ts):1/ts];
fxf=fxx(1:end-1);
fy=fyf(1:Nn/2);
fx=fxf(1:Nn/2);
% figure(3),plot(fx,abs(fy),'r.-');grid on
%% decomposition 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[v,p]=max(abs(fy));
u10=[zeros(1,p-1) v zeros(1,length(fy)-p)];
m=100;
alpha=1;
beta=1e-2;  %%the higher the noise level,the smaller the beta
    
k=0;
while k<m    
    wc0=sum(fx.*(abs(u10).^2))/sum(abs(u10).^2);
    wcr0=sum(fx.*(abs(fy-u10).^2))/sum(abs(fy-u10).^2);       
    u11=(fy.*(1+beta*(fx-wcr0).^2))./(1+alpha.*(fx-wc0).^2+beta*(fx-wcr0).^2);
    wc1=sum(fx.*(abs(u11).^2))/sum(abs(u11).^2);
    wcr1=sum(fx.*(abs(fy-u11).^2))/sum(abs(fy-u11).^2);

    u10=u11;
    if abs(wc0-wc1)<1e-3
        break;
    end
    wc0;
    wc1;
    wc0=wc1;
    wcr0=wcr1;
    k=k+1;
    
end
y11=ifft(u11,Nn);
yr11=real(y11);
yr1=yr11(LeL:LeL+N-1);
u1=yr1;
figure(3),plot(t,y1,'b-',t,u1,'r--');
%% elongation 2%%%%%

y=y-yr1;
sp=polyfit(t,y,2);
ytrend=polyval(sp,t);
y=y-ytrend;

xL=t(1:LeL);
yL=y(1:LeL);
fL=fft(yL,LeL);
f=(0:LeL-1)*fs/LeL;

[v1,p1]=max(abs(fL(1:end/2))*2/LeL);
fre1=f(p1);
A1=v1;
ph1=angle(fL(p1));
fLRemaining=[zeros(1,p1-1),v1,zeros(1,LeL/2-p1)];
[v2,p2]=max(abs(fL(1:end/2))*2/LeL-fLRemaining);
fre2=f(p2);
A2=v2;
ph2=angle(fL(p2));

if p1==1
    A1=A1/2;
end
if p2==1
    A2=A2/2;
end
tL=[t(1)-3*LeL*ts:ts:t(1)-ts];
ytrendL=polyval(sp,tL);
yLe=A1*cos(fre1*2*pi*tL+ph1)+A2*cos(fre2*2*pi*tL+ph2)+ytrendL;

%elongate in right
xR=t(end-LeL+1:end);
yR=y(end-LeL+1:end);
fR=fft(yR,LeL);
f=(0:LeL-1)*fs/LeL;

[v1,p1]=max(abs(fR(1:end/2))*2/LeL);
fre1=f(p1);
A1=v1;
ph1=angle(fR(p1));
fRRemaining=[zeros(1,p1-1),v1,zeros(1,LeL/2-p1)];

[v2,p2]=max(abs(fR(1:end/2))*2/LeL-fRRemaining);
fre2=f(p2);
A2=v2;
ph2=angle(fR(p2));

if p1==1
    A1=A1/2;
end
if p2==1
    A2=A2/2;
end

tR=[t(end)+ts:ts:t(end)+3*LeL*ts];
ytrendR=polyval(sp,tR);
yRe=A1*cos(fre1*2*pi*tR+ph1)+A2*cos(fre2*2*pi*tR+ph2)+ytrendR;
y=y+ytrend;
ye=[yLe y yRe];
te=[tL t tR];
yh = hilbert(ye); 
Nn=length(ye);
fyf=fft(yh,Nn);
fxx=[0:1/(Nn*ts):1/ts];
fxf=fxx(1:end-1);

fy=fyf(1:Nn/2);
fx=fxf(1:Nn/2);
% figure(8),plot(fx,abs(fy),'r-');grid on

%% decomposition2 %%%
[v,p]=max(abs(fy));
u10=[zeros(1,p-1) v zeros(1,length(fy)-p)];
m=100;
alpha=1;
beta=1e-2;
    
k=0;
while k<m    
    wc0=sum(fx.*(abs(u10).^2))/sum(abs(u10).^2);
    wcr0=sum(fx.*(abs(fy-u10).^2))/sum(abs(fy-u10).^2);    
    u11=(fy.*(1+beta*(fx-wcr0).^2))./(1+alpha.*(fx-wc0).^2+beta*(fx-wcr0).^2);
    wc1=sum(fx.*(abs(u11).^2))/sum(abs(u11).^2);
    wcr1=sum(fx.*(abs(fy-u11).^2))/sum(abs(fy-u11).^2);
    u10=u11;
    if abs(wc0-wc1)<1e-3
        break;
    end
    wc0=wc1;
    wcr0=wcr1;
    k=k+1;    
end

y11=ifft(u11,Nn);
yr11=real(y11);
yr1=yr11(3*LeL:3*LeL+N-1);
u2=yr1;
figure(4),plot(t,y2,'b-',t,u2,'r--');grid on
%%
%% elongation  3 %%%%
y=y-yr1;
sp=polyfit(t,y,2);
ytrend=polyval(sp,t);

y=y-ytrend;

LeL=500;
%elongate in left
xL=t(1:LeL);
yL=y(1:LeL);
fL=fft(yL,LeL);
f=(0:LeL-1)*fs/LeL;
[v1,p1]=max(abs(fL(1:end/2))*2/LeL);
fre1=f(p1);
A1=v1;
ph1=angle(fL(p1));

fLRemaining=[zeros(1,p1-1),v1,zeros(1,LeL/2-p1)];
[v2,p2]=max(abs(fL(1:end/2))*2/LeL-fLRemaining);
fre2=f(p2);
A2=v2;
ph2=angle(fL(p2));
if p1==1
    A1=A1/2;
end
if p2==1
    A2=A2/2;
end

tL=[t(1)-3*LeL*ts:ts:t(1)-ts];
ytrendL=polyval(sp,tL);
yLe=A1*cos(fre1*2*pi*tL+ph1)+A2*cos(fre2*2*pi*tL+ph2)+ytrendL;

%elongate in right
xR=t(end-LeL+1:end);
yR=y(end-LeL+1:end);
fR=fft(yR,LeL);
f=(0:LeL-1)*fs/LeL;
[v1,p1]=max(abs(fR(1:end/2))*2/LeL);
fre1=f(p1);
A1=v1;
ph1=angle(fR(p1));

fRRemaining=[zeros(1,p1-1),v1,zeros(1,LeL/2-p1)];
[v2,p2]=max(abs(fR(1:end/2))*2/LeL-fRRemaining);
fre2=f(p2);
A2=v2;
ph2=angle(fR(p2));

if p1==1
    A1=A1/2;
end
if p2==1
    A2=A2/2;
end


tR=[t(end)+ts:ts:t(end)+3*LeL*ts];
ytrendR=polyval(sp,tR);
yRe=A1*cos(fre1*2*pi*tR+ph1)+A2*cos(fre2*2*pi*tR+ph2)+ytrendR;

y=y+ytrend;
ye=[yLe y yRe];
te=[tL t tR];
yh = hilbert(ye); 
Nn=length(ye);
fyf=fft(yh,Nn);
fxx=[0:1/(Nn*ts):1/ts];
fxf=fxx(1:end-1);
fy=fyf(1:Nn/2);
fx=fxf(1:Nn/2);

%% decomposition 3%%
[v,p]=max(abs(fy));
u10=[zeros(1,p-1) v zeros(1,length(fy)-p)];
%u10=fy/2;
m=100;
alpha=1;
beta=1e-3;    
k=0;
while k<m    
    wc0=sum(fx.*(abs(u10).^2))/sum(abs(u10).^2);
    wcr0=sum(fx.*(abs(fy-u10).^2))/sum(abs(fy-u10).^2);    
    u11=(fy.*(1+beta*(fx-wcr0).^2))./(1+alpha.*(fx-wc0).^2+beta*(fx-wcr0).^2);
    wc1=sum(fx.*(abs(u11).^2))/sum(abs(u11).^2);
    wcr1=sum(fx.*(abs(fy-u11).^2))/sum(abs(fy-u11).^2);
    u10=u11;
    if abs(wc0-wc1)<1e-3
        break;
    end
    wc0=wc1;
    wcr0=wcr1;
    k=k+1;
    
end

y11=ifft(u11,Nn);
yr11=real(y11);
yr1=yr11(3*LeL:3*LeL+N-1);
u3=yr1;

figure(5),plot(t,y3,'b-',t,u3,'r--');grid on



%% elongation  4 %%%%
y=y-yr1;
sp=polyfit(t,y,2);
ytrend=polyval(sp,t);
y=y-ytrend;

LeL=500;
%elongate in left
xL=t(1:LeL);
yL=y(1:LeL);
fL=fft(yL,LeL);
f=(0:LeL-1)*fs/LeL;
[v1,p1]=max(abs(fL(1:end/2))*2/LeL);
fre1=f(p1);
A1=v1;
ph1=angle(fL(p1));
fLRemaining=[zeros(1,p1-1),v1,zeros(1,LeL/2-p1)];
[v2,p2]=max(abs(fL(1:end/2))*2/LeL-fLRemaining);
fre2=f(p2);
A2=v2;
ph2=angle(fL(p2));

if p1==1
    A1=A1/2;
end
if p2==1
    A2=A2/2;
end

tL=[t(1)-3*LeL*ts:ts:t(1)-ts];
ytrendL=polyval(sp,tL);
yLe=A1*cos(fre1*2*pi*tL+ph1)+A2*cos(fre2*2*pi*tL+ph2)+ytrendL;

%elongate in right
xR=t(end-LeL+1:end);
yR=y(end-LeL+1:end);

fR=fft(yR,LeL);
f=(0:LeL-1)*fs/LeL;

[v1,p1]=max(abs(fR(1:end/2))*2/LeL);
fre1=f(p1);
A1=v1;
ph1=angle(fR(p1));
fRRemaining=[zeros(1,p1-1),v1,zeros(1,LeL/2-p1)];
[v2,p2]=max(abs(fR(1:end/2))*2/LeL-fRRemaining);
fre2=f(p2);
A2=v2;
ph2=angle(fR(p2));

if p1==1
    A1=A1/2;
end
if p2==1
    A2=A2/2;
end

tR=[t(end)+ts:ts:t(end)+3*LeL*ts];
ytrendR=polyval(sp,tR);
yRe=A1*cos(fre1*2*pi*tR+ph1)+A2*cos(fre2*2*pi*tR+ph2)+ytrendR;
y=y+ytrend;
ye=[yLe y yRe];
te=[tL t tR];
yh = hilbert(ye);    
Nn=length(ye);
fyf=fft(yh,Nn);
fxx=[0:1/(Nn*ts):1/ts];
fxf=fxx(1:end-1);
fy=fyf(1:Nn/2);
fx=fxf(1:Nn/2);


%% %% decomposition 4%%
[v,p]=max(abs(fy));
u10=[zeros(1,p-1) v zeros(1,length(fy)-p)];
m=100;
alpha=1;
beta=1e-5;
    
k=0;
while k<m    
    wc0=sum(fx.*(abs(u10).^2))/sum(abs(u10).^2);
    wcr0=sum(fx.*(abs(fy-u10).^2))/sum(abs(fy-u10).^2);
      
    u11=(fy.*(1+beta*(fx-wcr0).^2))./(1+alpha.*(fx-wc0).^2+beta*(fx-wcr0).^2);
    wc1=sum(fx.*(abs(u11).^2))/sum(abs(u11).^2);
    wcr1=sum(fx.*(abs(fy-u11).^2))/sum(abs(fy-u11).^2);
    u10=u11;
    if abs(wc0-wc1)<1e-3
        break;
    end
    wc0=wc1;
    wcr0=wcr1;
    k=k+1;
    
end

y11=ifft(u11,Nn);
yr11=real(y11);
yr1=yr11(3*LeL:3*LeL+N-1);
u32=yr1;

figure(6),plot(t,y3,'b-',t,u32,'r--');grid on
figure(7),plot(t,y3,'b-',t,u32+u3,'r--',t,u3,'k--');grid on

