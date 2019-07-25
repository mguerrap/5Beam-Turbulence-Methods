% Spectral Estimates Function
% Maricarmen Guerra

% Inputs
% mts = time series
% dt = seconds
% WOSA Lengths of segments

% Ns1 
% Ns2 


% Returns:
% Periodogram with Hanning Data Taper
% WOSA 1 and WOSA 2, CI and Band Width

function [fh,sx_h,fks1,Swosa1,fks2,Swosa2,BH1,BH2,CIs1_i,CIs1_f,CIs2_i,CIs2_f]=spectralestimate(mts,dt,Ns1,Ns2)

% Detrend
mts=mts-mean(mts);
samplingfreq=1/dt;

N=length(mts);

t=0:dt:dt*(N-1); % Real time in sec
fn=1/(2*dt);

kk=-N:1:N;

dfk=1/(N*dt);

fk=kk*dfk; %-N*dfk:dfk:N*dfk;
fk_1=1./fk;

ini=find(fk==0); %-fn);
fin=find(fk==fn);


Hwdw=hann(N);
Hwdw=Hwdw/sqrt(sum(Hwdw.^2)); % Normalizing the data taper

[sx_h,fh]=periodogram(mts,Hwdw,fk,1/dt);

   
% figure(1)
% clf
% plot(fh(ini:fin),10*log10(sx_h(ini:fin)),'r')
% hold off
% %xlim([0 fx(end)])
% xlabel('f'); ylabel('dB');
% title(['Periodogram and SDF with Hanning data taper']);
%legend('S^{(p)}','S^{(d)}')


%% WOSA Estimate

OL=0.5;

% Solo para las f positivas
ks1=1:1:Ns1;
ks2=1:1:Ns2;

dfs1=1/(Ns1*dt);
dfs2=1/(Ns2*dt);

fks1=ks1*dfs1; %-N*dfk:dfk:N*dfk;
fks2=ks2*dfs2; %-N*dfk:dfk:N*dfk;

taus1=0:1:Ns1-1;
taus2=0:1:Ns2-1;

Sini1=[1:Ns1*OL:N-2*Ns1*OL];
Sini2=[1:Ns2*OL:N-2*Ns2*OL];

Sfin1=Sini1+Ns1-1;
Sfin2=Sini2+Ns2-1;

% Create shorter time series and create hanning data taper for each segment
HwdwS1=hann(Ns1);
HwdwS1=HwdwS1./sqrt(sum(HwdwS1.^2)); % Normalizing the data taper

HwdwS2=hann(Ns2);
HwdwS2=HwdwS2./sqrt(sum(HwdwS2.^2)); % Normalizing the data taper

clear i
for j=1:length(Sini1)
    Xseg1(:,j)=mts(Sini1(j):Sfin1(j),1);
    
    % Data Taper
    XHseg1(:,j)=Xseg1(:,j).*HwdwS1;
    
    % S(d)_j
    
    % Fourier Transform
    XHs1_ft(:,j)=fft(XHseg1(:,j));
    
    
    %     for p=1:length(fks1)
    %         XHs1_ft(p,j)=dt^0.5*sum(XHseg1(:,j).*exp(-1i*2*pi*dt*fks1(p).*taus1'));
    %
    %     end
    
    % DSE
    S_l1(:,j)=dt*(abs(XHs1_ft(:,j))).^2;

end

Swosa1=mean(S_l1,2);

clear i
for j=1:length(Sini2)
    Xseg2(:,j)=mts(Sini2(j):Sfin2(j),1);
    
    % Data Taper
    XHseg2(:,j)=Xseg2(:,j).*HwdwS2;
    
    % S(d)_j
    XHs2_ft(:,j)=fft(XHseg2(:,j));

%     for p=1:length(fks2)
%         XHs2_ft(p,j)=dt^0.5*sum(XHseg2(:,j).*exp(-1i*2*pi*dt*fks2(p).*taus2'));
%         
%     end

    S_l2(:,j)=dt*(abs(XHs2_ft(:,j))).^2;
    %S_l2(:,j)=(abs(XHs2_ft(:,j))).^2;
end

Swosa2=mean(S_l2,2);

iniW1=find(fks1==0);
finW1=find(fks1==fn);

iniW2=find(fks2==0);
finW2=find(fks2==fn);


% Bandwidth
jj=1;
pp=1;
for j=-(Ns1-1):1:(Ns1-1)
    pp=1;
    for p=0:1:(Ns1-abs(j)-1)
        hsh1(pp)=dt*HwdwS1(pp+abs(j))*HwdwS1(pp);
        pp=pp+1;
    end
    hsht1(jj)=sum(hsh1); 
    jj=jj+1;
end

BH1=dt/sum(hsht1.^2);
jj=1;
pp=1;
for j=-(Ns2-1):1:(Ns2-1)
    pp=1;
    for p=0:1:(Ns2-abs(j)-1)
        hsh2(pp)=dt*HwdwS2(pp+abs(j))*HwdwS2(pp);
        pp=pp+1;
    end
    hsht2(jj)=sum(hsh2); 
    jj=jj+1;
end
BH2=dt/sum(hsht2.^2);

% CI
nus1=3.79*N/Ns1;
nus2=3.79*N/Ns2;

% CI real
CIs1_i=nus1/chi2inv(1-0.05/2,nus1).*Swosa1;
CIs1_f=nus1/chi2inv(0.05/2,nus1).*Swosa1;

CIs2_i=nus2/chi2inv(1-0.05/2,nus2).*Swosa2;
CIs2_f=nus2/chi2inv(0.05/2,nus2).*Swosa2;

% CI in dB
%CIs1=(chi2inv(1-0.05/2,nus1)/chi2inv(0.05/2,nus1));
%CIs2=(chi2inv(1-0.05/2,nus2)/chi2inv(0.05/2,nus2));

