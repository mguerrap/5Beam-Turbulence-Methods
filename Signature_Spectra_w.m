clear all, close all, clc

% Spectral Analysis of beam velocities
% It is set for vertical beam (beam5), but can be changed to any beam or to
% ENU velocities (see line 62)

% For Nortek Signature 1000 measurements
% Uses Signature Files created with Signatue_QC.m

tic


% File Paths and prefix

% Where is the QC Data:
fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefix = ['SignatureData_QC_Bin'];

%Where to save the spectra
savepath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/SDF_Signature'];

% Bins
Nbin=20;
% Ensembles to use
Ens1=10; %First
EnsN=43; %Last

% Segment Length for WOSA
Ns1=300;
Ns2=400;

% Maximum and Minimum U for spectral average m/s
umax=2;
umin=0.5; 
dumean=0.25;

%Frequency at which noise appears
fnoise=2; 

%--------------------------------------------------------------------------
for bin=1:Nbin
    
    fname = [prefix int2str(bin) '.mat'];        
    load([ fpath '/' fname ])
    
    
    [Nt Nens]=size(SigData.vbeam1);
    fs=SigData.fs;
    dt=1/fs;
    
    Ens_size=SigData.Ens_size;
    
    % Spectrums
    
    for fi=Ens1:EnsN 
        
        % Spectra at every burst
        
        fi
        
        
        %%%%% Choose which beam velocity are you going to use
        %%%%% Change here which time series you want to use:
        
        w=SigData.vbeam5(:,fi);
        
        
                
        % All estimated from this function are two-sided, thus power
        % density is multiplied by two to get all variance in one-side
        % Variance is checked bellow from f=0 to f=Nyquist Frequency
        
        [fh,sx_h,fks1,Swosa1,fks2,Swosa2,BH1,BH2,CIs1_i,CIs1_f,CIs2_i,CIs2_f]=spectralestimate(w,dt,Ns1,Ns2);
        SigSDF.vf(:,fi)=fh;
        SigSDF.DSP(:,fi)=sx_h*2; % estimates correspond to double side spectra, multiply by two to get total energy for each frequency
        dfh=fh(2)-fh(1);
        
        SigSDF.fks1(:,fi)=fks1;
        SigSDF.DSP_W1(:,fi)=Swosa1*2;
        dfks1=fks1(2)-fks1(2);
        
        SigSDF.fks2(:,fi)=fks2;
        SigSDF.DSP_W2(:,fi)=Swosa2*2;
        dfks2=fks2(2)-fks2(2);
        
        % BandWidth and CIs for WOSA
        SigSDF.BH1(fi)=BH1;
        SigSDF.BH2(fi)=BH2;
        SigSDF.CIs1_i(:,fi)=CIs1_i;
        SigSDF.CIs2_i(:,fi)=CIs2_i;
        SigSDF.CIs1_f(:,fi)=CIs1_f;
        SigSDF.CIs2_f(:,fi)=CIs2_f;
        

        % Alternative: Use pwelch
        
        [upsd f] = pwelch(  detrend(w) , Ns1, [], [], fs,'onesided');
        %upsd(1) = []; f(1) = []; % remove mean
        SigSDF.PW(:,fi) = real(upsd);  % power spectral density of velocity mag, m^2 s^-2 Hz^-1
        SigSDF.PW_f(:,fi)=f;
        dfpw=f(3)-f(2);
        
        % Horizontal Velocity
        
        % Mean horizontal velocity
        SigSDF.uH(:,fi)=sqrt(SigData.u_east(:,fi).^2+SigData.v_north(:,fi).^2);
        SigSDF.uH_mean(fi)=nanmean(SigSDF.uH(:,fi));
        SigSDF.uH_var(fi)=nanvar(SigSDF.uH(:,fi));
        
        % If uH<0.5 do not use spectrum
        if SigSDF.uH_mean(fi)<=0.5
            
            SigSDF.DSP(:,fi)=NaN;
            SigSDF.PW(:,fi)=NaN;
            SigSDF.DSP_W1(:,fi)=NaN;
            SigSDF.DSP_W2(:,fi)=NaN;
        end
        
        % Check variances
        % Total energy from 0 to Nyquist frequencies==Variance
        f0h=find(fh==0);
        fnh=find(fh==fs/2);       
        fnpw=find(f==fs/2);
        fn1=find(fks1==fs/2);
        fn2=find(fks2==fs/2);
        
        % Variance Check
        var_w(fi,1)=var(detrend(w));
        var_w(fi,2)=trapz(SigSDF.vf(f0h:fnh,fi),real(SigSDF.DSP(f0h:fnh,fi)));
        var_w(fi,3)=trapz(SigSDF.PW_f(1:fnpw,fi),real(SigSDF.PW(1:fnpw,fi)));
        var_w(fi,4)=trapz(SigSDF.fks1(1:fn1,fi),real(SigSDF.DSP_W1(1:fn1,fi)));
        var_w(fi,5)=trapz(SigSDF.fks2(1:fn2,fi),real(SigSDF.DSP_W2(1:fn2,fi)));
        
        % Instrument Noise estimation

        fnnh=find(fh==fnoise);
        fnnpw=find(f==fnoise);
        fnn1=find(fks1==fnoise);
        fnn2=find(fks2==fnoise);
        
        noise_w(fi,1)=nanmean(real(SigSDF.DSP(fnnh:fnh,fi)));
        noise_w(fi,2)=nanmean(real(SigSDF.PW(fnnpw:fnpw,fi)));
        noise_w(fi,3)=nanmean(real(SigSDF.DSP_W1(fnn1:fn1,fi)));
        noise_w(fi,4)=nanmean(real(SigSDF.DSP_W2(fnn1:fn2,fi)));
        
    end

    %% Average spectrums with similar umean:
    umeanNotNaN=SigSDF.uH_mean(find(~isnan(SigSDF.uH_mean)));   
    xumean=umin:dumean:umax;
    [bincounts,ind]= histc(abs(SigSDF.uH_mean),xumean); %histc(abs(umeanNotNaN),xumean);
    [nelements,centers]=hist(abs(SigSDF.uH_mean),xumean);
    SigSDF.uH_mean_mean=centers;
    
    
    for i=1:length(bincounts)
        promediar=find(ind==i);
        
        SigSDF.DSP_mean(:,i)=nanmean(SigSDF.DSP(:,promediar),2);
        SigSDF.PW_mean(:,i)=nanmean(SigSDF.PW(:,promediar),2);
        SigSDF.DSP_W1_mean(:,i)=nanmean(SigSDF.DSP_W1(:,promediar),2);
        SigSDF.DSP_W2_mean(:,i)=nanmean(SigSDF.DSP_W2(:,promediar),2);

        num_ens(i)=length(promediar);
    end
    SigSDF.variances=var_w;
    SigSDF.noises=noise_w;
    
    archivo=['Sig_SDF_bin' int2str(bin) '.mat'];
    

    save([savepath '/' archivo],'SigSDF');
    
    for i=1:5
        var_w_wo(:,i)=var_w(:,i)./var_w(:,1);
    end
    var_w_wo_mean(bin,:)=nanmean(var_w_wo,1);
    noise_w_mean(bin,:)=nanmean(noise_w);
    
    clearvars -except bin fnoise var_w var_w_wo var_w_wo_mean noise_w_mean fpath prefix savepath ...
        Ens1 EnsN Ns1 Ns2 umin umax dumean Nbin
    
    
    
end
toc