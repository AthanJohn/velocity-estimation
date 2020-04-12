% Project - Velocity/Distance Estimation
% ATHANASIOU IOANNIS






clear all; clc; clf;

%Receiver vehicle velocity
v1=input('v1= '); 

%Transmitter vehicle velocity
x=(200+v1).*rand(1)-v1;
x=x+v1;

%Loops
N=10000;


SNRdB=0:2:20;

%SNR to linear
SNR=10*log10(SNRdB);



%random channel
r=10000;
k=0;
h1=rand(1,floor(N/r)+1);
h=upsample(h1,r);
for i=1:r:N
    k=k+1;
    for j=i:i+r-1
        h(j)=h1(k);
    end
end
if size(h,2)>N
    for i=size(h,2):-1:N+1
        h(:,i)=[];
    end
end


hx=h*x;
hx_length=length(hx);


y=zeros(length(SNR),N);
Ey=zeros(length(SNR),1);
vary=Ey;
XLMMSE=y;
MSE=Ey;
mse=MSE;
covxy=Ey;
mean_MSE=zeros(1,length(SNR));

Es=sum(abs(hx).^2)/(hx_length); %signal power

%Mse for each SNR value
for i=1:length(SNR)
    for j=1:N
        %noise
        sigp = 10*log10(norm(hx,2)^2/numel(hx));
        snr=sigp-SNRdB(i);
        noisep=10^(snr/10);
        w=sqrt(noisep)*randn(size(hx));
        
        
        y= hx + w;
        
       
        Ex=sum(hx)/N;
        Ey=sum(y)/N;
        vary=sum(abs(y-Ey).^2)/(N-1);
        covxy=sum(x*y')/N - Ex*Ey;
        XLMMSE=Ex+(covxy/vary)*(y-Ey);
        MSE(j,i)=mean((x-XLMMSE).^2);
        mse(j,i)=mean((y-x).^2);
        
    end
end

for i=1:length(SNR)
    mean_MSE(length(SNR)-i+1)=sum(MSE(:,i))/N;
    mean_mse(i)=sum(mse(:,i))/N;
end


semilogy(SNRdB,mean_MSE, 'r--*');
ylabel('MSE');
xlabel('SNR');
title('MSE vs SNR');
