% Estimation of Dissipation Rate from spectra
% Using Sww WOSA 2 and regression to estimate slope and dissipation for beam5
% Define which spectra will use in line 48
% Using compensated spectra and finding slope closer to 0 in fit

% This is only for the 5th Beam spectra

% Uses files created with Signature_Spectrum_w.m

clc
clear all
close all

%--------------------------------------------------------------------------
% Data Info
fs=8; %Sampling frequency

prefix=['Sig_SDF_bin'];
fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/SDF_Signature'];
savepath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/Signature_Dissipation'];  
% Bins
Nbin=20;
% Ensembles to use, which have not NaN spectra
ensembles = [12:1:25, 29:43];
nnoise=4; % Where is the noise estimate for the Wosa2 in the noise matrix.

% Freq ranges to try for inertial subrange
fstart=[0.1 0.2 0.3 0.4 0.5];
fend=[0.6 0.6 0.6 0.6 0.6];
    
% Constants
a=0.69; %Sww constant
vk=0.41; % Von Karman
kolm=0.55; % Kolmogorov constant
rhow=1025; % Water density

%--------------------------------------------------------------------------
for bin=1:Nbin
    
    
    fname = [prefix int2str(bin) '.mat'];        
    load([ fpath '/' fname ])
    
    

    % Using WOSA 2, can change to other spectral estimates for comparison
    % Define here which spectra estimate will use
    f=SigSDF.fks2(:,10);
    Sww=SigSDF.DSP_W2-repmat(SigSDF.noises(:,nnoise)',length(f),1); %WOSA2 Ns=400
    
    Ubar=SigSDF.uH_mean; %Mean horizontal velocity
    
    Ubar=SigSDF.uH_mean'; %Mean horizontal velocity
    
    Nyq=find(f==fs/2);
    
  
    % Find the freq and Sww to calculate E from
    for i=1:length(fstart)
        df=f(:,1)-fstart(i);
        istartaux=find(abs(df)==min(abs(df)));
        istart(i)=istartaux(1);
        clear istartaux
    end
    
    for i=1:length(fend)
        df=f(:,1)-fend(i);
        iendaux=find(abs(df)==min(abs(df)));
        iend(i)=iendaux(end);
        clear iendaux
    end
    
    bflat=NaN(43,1);
    
    for i=ensembles % For the actual ensembles when there are spectrums i=[14:25 31:41 43]
        
        for j=1:length(fstart)
            
            ss=istart(j);
            ee=iend(j);
            
            Sww_flat=Sww(ss:ee,i).*f(ss:ee,1).^(5/3);
            mSflat(j)=mean(Sww_flat);
            stdSflat(j)=std(Sww_flat);
            
            % Regression
            brob = robustfit(f(ss:ee,1),Sww_flat);%robustfit(log(f2(ss:ee,1).^(-5/3)),log(Sww2(ss:ee,i)));
            mS(j)=brob(2);
            cS(j)=brob(1);
            Sfit_flat=cS(j)+mS(j)*f(ss:ee,1);
            %Sfit=exp(Sfit);
            RMSEfit(j)=sqrt(mean((Sww_flat-Sfit_flat).^2));
            SlopeDiff(j)=(abs(mS(j)));
            slopes(j,i)=mS(j);
            
        end
        

        fitmin(i)=find(SlopeDiff==min(SlopeDiff));
        s2=istart(fitmin(i));
        e2=iend(fitmin(i));
        SS_flat2=Sww(s2:e2,i).*f(s2:e2,1).^(5/3);
        Sww_f_mean2=mean(SS_flat2);
        
        % Dissipation Rate
        erate(i)=(1/a*Sww_f_mean2).^(3/2)*2*pi/abs(Ubar(i));
        
        % Flat Slopes
        bflat(i)=mS(fitmin(i));
        Sbestfit= cS(fitmin(i))+mS(fitmin(i))*f(s2:e2,1);

        %confidence interval for the slope
        nn=length(SS_flat2);
        SSE=sum((SS_flat2-Sbestfit).^2); % sum of squared errors
        sigma2=SSE./(nn-2);
        sigma=sqrt(sigma2);
        sB1=sqrt(sigma2/sum(f(s2:e2,1)-mean(f(s2:e2,1))));
        alpha=0.05;
        tstat=-tinv(alpha/2,nn-2);
        bflat_ci(i,1)=tstat*sigma;
        bflat_ci(i,2)=tstat*sigma;
        
        % Errors in Dissipation Rate from flat spectrum
        % Using error propagation
        erate_error(i)=3/2*1/a*(2*pi/abs(Ubar(i))).^(2/3)*std(SS_flat2)*mean(SS_flat2)^0.5;
        
        % U mean effect
        erate_umeanPlus(i)=(1/a*Sww_f_mean2).^(3/2)*2*pi/(1.1*abs(Ubar(i)));
        erate_umeanLow(i)=(1/a*Sww_f_mean2).^(3/2)*2*pi/(0.9*abs(Ubar(i)));
 
    end
    
    eratefile=[savepath '/Dissip_Spectra_bin' int2str(bin) '.mat'];
    save(eratefile,'erate','Ubar','erate_error','erate_umeanPlus','erate_umeanLow');
    
    
    clearvars -except bin slopes fitmin a vk kolm rhow fs nnoise ...
        prefix savepath fpath Nbin ensembles nnoise fstart fend fs
    bin
end