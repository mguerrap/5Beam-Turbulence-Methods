clc
clear all
close all

% 5 beam E-P profiles

% Use this example to plot mean dissipation and production profiles


% Bins
Nbin=20;

% Ensembles
Ens1=10;
EnsN=43;
NEns=EnsN-Ens1+1;

% Dissipation Rate files:
fpathD = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/Signature_Dissipation'];

% Spectra
prefixSp = ['Dissip_Spectra_bin'];

% Structure Function
fnameSF = ['SF_DissipationRate.mat']; %This file contains the structure function for el depths

% Production File
fpathP = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/Signature_Production'];
fnameP=['Production.mat'];

% Mean Flow
fpathMF=['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefixMF=['Signature_Usigned_Bin'];


% Range
bin=1;
fpath3=['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefix3 = ['SignatureData_QC_Bin'];
fname3 = [prefix3 int2str(bin) '.mat'];        
load([ fpath3 '/' fname3 ])

z=SigData.range;


%% ------------------------------------------------------------------------ 
% Load Data for Plots

GoodEnsAI=43;
NensAI=34;

% Data from dissipation rate from the spectra
for bin=1:Nbin
    
    eratefile=[prefixSp int2str(bin) '.mat'];
    load([ fpathD '/' eratefile])
    
    le=length(erate);
    erateall(1:le,bin)=erate'; %mean(Sww_flat) with best fit to 0 slope of Sww_flat
    erate_error(1:le,bin)=erate_error';

    % U Signed
    archivo=[prefixMF int2str(bin) '.mat'];
    load([fpathMF '/' archivo])
    Usigned(:,bin)=nanmean(Ens_Usigned,1)';
    
    %Usigned(abs(Usigned(:,bin))<0.5,bin)=NaN;
    
    clear Ens_Usigned
end

Usigned=Usigned(Ens1:EnsN,:);


esnan=find(erateall==0);
erateall2(esnan)=NaN;
clear esnan

% Saves only the useful ones
erateSp=erateall(Ens1:EnsN,:);
erateSp=erateSp';

% Structure Function Results


load([fpathD '/' fnameSF]);

esnan=find(abs(Usigned)<0.5);

SF.epsilon5(esnan)=NaN;
SF.eflat5(esnan)=NaN;


% Production
load([fpathP '/' fnameP]);

% Ebb and Flood
flood=find(nanmean(Usigned,2)>0);
ebb=find(nanmean(Usigned,2)<=0);


%%
figure(1)
clf
%set(gcf,'Position',[257 85 550 350])

% Ebb
%s1=subplot(2,2,1)
s1=axes('Position',[0.13    0.22    0.3259    0.7])
% SF
h2_ebb=semilogx(nanmean(SF.eflat5(:,ebb),2),z,'k','LineWidth',2)
hold on
% Spectra
h5_ebb=semilogx(nanmean(erateSp(:,ebb),2),z,'b','LineWidth',2)
% Production
h7_ebb=semilogx(nanmean(abs(Prod.P_total(:,ebb)),2),z(1:end-1),'r','LineWidth',2)
hold off
set(gca,'FontSize',16)
xlabel('P, \epsilon (m^{2}s^{-3})','FontSize',16)
ylabel('z (m)','FontSize',16)
title('Ebb')

% Flood
% SF
%s2=subplot(2,2,2)
s2=axes('Position',[0.57    0.22    0.3259    0.7])
h2_flood=semilogx(nanmean(SF.eflat5(:,flood),2),z,'k','LineWidth',2)
hold on
% Spectra
h5_flood=semilogx(nanmean(erateSp(:,flood),2),z,'b','LineWidth',2)
% Production
h7_flood=semilogx(nanmean(abs(Prod.P_total(:,flood)),2),z(1:end-1),'r','LineWidth',2)
hold off
set(gca,'FontSize',16)
xlabel('P, \epsilon (m^{2}s^{-3})','FontSize',16)
ylabel('z (m)','FontSize',16)
title('Flood')

legend([h2_ebb,h5_ebb,h7_ebb],'\epsilon Structure Function','\epsilon TKE Spectra','Production',...
'Position',[280 25 1 0.3],'Orientation','Horizontal','FontSize',16)

