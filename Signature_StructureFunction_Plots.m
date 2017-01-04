% Example Plot of Spectral Density Functions

clc
clear
close all

% Select Bin
bin = 10;

% Sampling frequency
fs=8;

% Ensembles
Ens1=10;
EnsN=34;
NEns=EnsN-Ens1+1;

% Maximum velocity for colormap
Umax=2;
Umin=0.5;
dU=0.5;

% Structure Function D(z,r) path and file

fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/Signature_Dissipation'];
fname = ['Dzr_Beam5.mat']; %This file contains the structure function for el depths

load([ fpath '/' fname ])

% Mean Flow
fpath2=['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
fname2=['Signature_Usigned_Bin' int2str(bin) '.mat'];
load([fpath2 '/' fname2])

% Range
fpath3=['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefix3 = ['SignatureData_QC_Bin'];
fname3 = [prefix3 int2str(bin) '.mat'];        
load([ fpath3 '/' fname3 ])

z=SigData.range;

%% --------------------------------------------------------------------------
% D(z,r) Plot

Usigned(1:NEns)=nanmean(Ens_Usigned(:,Ens1:EnsN),1)';
Usigned(abs(Usigned(:))<0.5)=NaN;
maxr = z(bin);

for i=1:NEns
    D(i,:)=Dzr(i).D5(bin,:);
    r(i,:)=Dzr(i).r5(bin,:);
end

figure(1)
clf
colormap jet
cmap=colormap;

for i=1:NEns
 
    goodpts = ~isnan( D(i,:) )  &  r(i,:)< maxr  & r(i,:) >= 0;
    
    if ~isnan(Usigned(i))
    cindex=min([ceil(abs(Usigned(i))/Umax*64); 64]);
    semilogy(r(i,goodpts),D(i,goodpts),'o-','MarkerSize',5,'color',cmap(cindex,:),'MarkerFaceColor',cmap(cindex,:),'LineWidth',2)
    hold on
    end
end
text(4,4.2e-2,'~r ^{2/3}','FontSize',20)
semilogy([1:0.1:10],10^-2*[1:0.1:10].^(2/3),'--k')
hold off
xlabel('r (m)')
ylabel('D(z,r)')
cb = colorbar;
caxis([Umin Umax])
set(cb,'YTickLabel',[Umin:dU:Umax]);
ylabel(cb,'<U> (ms^{-1})','FontSize',16)
set(gca,'FontSize',16)
axis([0 10 10^-5 10^-1])

