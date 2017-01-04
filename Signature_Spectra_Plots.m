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
EnsN=43;
NEns=EnsN-Ens1+1;
% Maximum velocity for colormap
Umax=2;
Umin=0.5;
dU=0.5;
% Spectra path and file
prefix=['Sig_SDF_bin'];
fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/SDF_Signature'];
fname = [prefix int2str(bin) '.mat'];
load([ fpath '/' fname ])

% Mean Flow
fpath2=['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
fname2=['Signature_Usigned_Bin' int2str(bin) '.mat'];
load([fpath2 '/' fname2])


%% --------------------------------------------------------------------------
% WOSA 1 Plot
fks1=SigSDF.fks1(:,Ens1);
fks2=SigSDF.fks2(:,Ens1);

fn1=find(fks1==fs/2);
fn2=find(fks2==fs/2);

Usigned=nanmean(Ens_Usigned,1)';
Usigned(abs(Usigned(:))<Umin)=NaN;

colormap jet
cmap=colormap;

figure(1)
clf
colormap jet
cmap=colormap;

for i=Ens1:EnsN
    
    if ~isnan(Usigned(i))
        cindex=min([ceil(abs(Usigned(i))/Umax*64); 64]);
        
        loglog(SigSDF.fks1(1:fn1,i),SigSDF.DSP_W1(1:fn1,i),'color',cmap(cindex,:),'LineWidth',2)
        hold on
    end
    
    
end
plot([.2 1.5],1e-2*[.2 1.5].^(-5/3),'k--')
text(0.8,2.5e-2,'~f ^{-5/3}','FontSize',20)
hold off
axis([0 10 0 10^2])
set(gca,'FontSize',16)
cb = colorbar;
caxis([Umin Umax])
set(cb,'YTickLabel',[Umin:dU:Umax]);
ylabel(cb,'<U> (ms^{-1})','FontSize',14)
xlabel('Frequency (Hz)','FontSize',18)
ylabel(['S (m^{2}/s^{2}Hz^{-1})'],'FontSize',18)

%saveas(gcf,['Sww_Wosa1_N300_SignatureBin' int2str(bin) ''],'epsc')
