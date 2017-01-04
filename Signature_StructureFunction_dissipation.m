% AI May 2015 Dissipation from Structure Function

% Uses dissipation.m and structureFunction.m

clc
clear
close all

% Data: Need along beam velocities for all depths in form of v1=v1(z,time);

% For spectra we have data saved for each time burst for each bin, we need
% to open each file and form v(depth,time);

sigma_w=0.009; %m/s From Spectra

% File Paths and prefix

% Where is the QC Data:
fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefix = ['SignatureData_QC_Bin'];

%Where to save the Structure Function and the Dissipation
savepath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/Signature_Dissipation'];


% Bins
Nbin=20;
% Ensembles to use
Ens1=10; %First
EnsN=43; %Last
% Beam slanted angle
theta=25; %degrees


for bin=1:Nbin
    
    fname = [prefix int2str(bin) '.mat'];        
    load([ fpath '/' fname ])
    
    fs=SigData.fs;
    dt=1/fs;
    
    Ens_time(:,bin)=SigData.Ens_time; %Beggining of each ensemble
    
    [Nt Nens]=size(SigData.vbeam1);
    
    time0=[0:dt:Nt*dt-dt]; % Seconds
    
    fik=1;
    for fi=Ens1:EnsN %Good ensembles
    
        v1(bin,:,fik)=SigData.vbeam1(:,fi);
        v2(bin,:,fik)=SigData.vbeam2(:,fi);
        v3(bin,:,fik)=SigData.vbeam3(:,fi);
        v4(bin,:,fik)=SigData.vbeam4(:,fi);
        v5(bin,:,fik)=SigData.vbeam5(:,fi);
        
        timeEns=Ens_time(fi,bin)+time0/3600/24; %days
        timeV(bin,:,fik)=timeEns;
        fik=fik+1;
    end
    z=SigData.range;
end

%% Structure function in time
[Nz Nt Nens]=size(timeV);

% estimate dissipation [W/Kg] from 8 Hz beam data
disspts = Nt; % number of pts (in time) to use for each structure function

zbeam=z./cos(theta*pi/180);
zbeam5=z;

dr=zeros(length(zbeam));
close all

for fi=1:Nens
    b1=v1(:,:,fi);
    b2=v2(:,:,fi);
    b3=v3(:,:,fi);
    b4=v4(:,:,fi);
    b5=v5(:,:,fi);
    
    % 2/3 power law
    % For All beams: is very slow, so we are only estimating from the 5th
    % beam:
    
    %%beam=1
    %     [tke1(:,fi) epsilon1(:,fi) residual1(:,fi) A1(:,fi) Aerror1(:,fi) N1(:,fi) Nerror1(:,fi) rmax1(:,fi)  eflat1(:,fi) rmaxflat1(:,fi) eflat_E1(:,fi)] = ...
    %         dissipationMG_SF( b1, zbeam, disspts, 0,dr,1,sigma_w);  
    %%beam=2
    %     [tke2(:,fi) epsilon2(:,fi) residual2(:,fi) A2(:,fi) Aerror2(:,fi) N2(:,fi) Nerror2(:,fi) rmax2(:,fi)  eflat2(:,fi) rmaxflat2(:,fi) eflat_E2(:,fi)] = ...
    %         dissipationMG_SF( b2, zbeam, disspts, 0,dr,2,sigma_w);
    %%beam=3
    %     [tke3(:,fi) epsilon3(:,fi) residual3(:,fi) A3(:,fi) Aerror3(:,fi) N3(:,fi) Nerror3(:,fi) rmax3(:,fi)  eflat3(:,fi) rmaxflat3(:,fi) eflat_E3(:,fi)] = ...
    %         dissipationMG_SF( b3, zbeam, disspts, 0,dr,3,sigma_w);
    %%beam=4
    %     [tke4(:,fi) epsilon4(:,fi) residual4(:,fi) A4(:,fi) Aerror4(:,fi) N4(:,fi) Nerror4(:,fi) rmax4(:,fi)  eflat4(:,fi) rmaxflat4(:,fi) eflat_E4(:,fi)] = ...
    %         dissipationMG_SF( b4, zbeam, disspts, 0,dr,4,sigma_w);
    %beam=5
    [tke5(:,fi) epsilon5(:,fi) residual5(:,fi) A5(:,fi) Aerror5(:,fi) N5(:,fi) Nerror5(:,fi) rmax5(:,fi) eflat5(:,fi) rmaxflat5(:,fi) eflat_E5(:,fi) D5 r5] = ...
         dissipationMG_SF( b5, zbeam5, disspts, 0,dr,5,sigma_w);
     
    % Save also the structure function
     Dzr(fi).D5=D5;
     Dzr(fi).r5=r5;

    fi
end



%% Plots
figure(5)
clf
semilogx(nanmean(epsilon5,2),zbeam5,'r','LineWidth',2)
hold on
semilogx(nanmean(eflat5,2),zbeam5,'k','LineWidth',2)
% for fi=1:Nens
%     semilogx(epsilon5(:,fi),zbeam,'b')
% end
hold off
xlabel('Epsilon')
ylabel('Z (m)')
axis([10^-6 10^-4 0 20])

%% Save

% SF.epsilon1=epsilon1;
% SF.epsilon2=epsilon2;
% SF.epsilon3=epsilon3;
% SF.epsilon4=epsilon4;
SF.epsilon5=epsilon5;

% SF.eflat1=eflat1;
% SF.eflat2=eflat2;
% SF.eflat3=eflat3;
% SF.eflat4=eflat4;
SF.eflat5=eflat5;

% SF.eflat_E1=eflat_E1;
% SF.eflat_E2=eflat_E2;
% SF.eflat_E3=eflat_E3;
% SF.eflat_E4=eflat_E4;
SF.eflat_E5=eflat_E5;

SF.timeEns=Ens_time(Ens1:EnsN,1);


savefile=[savepath '/SF_DissipationRate.mat'];

save(savefile, 'SF')

savefile2=[savepath '/Dzr_Beam5.mat'];
save(savefile2, 'Dzr')


