% TKE Production from 5 beam RS


% Uses RS_5beam.m 

clc
clear
close all

% We need to have the Reynolds Stresses first as 'ReynoldsStress.mat'

% RS from 5 beam AD2CP
% RS are in Earth Coordinates
% Uses RS from 5 beams and shear from mean flow

% Where is the QC Data:
fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefix = ['SignatureData_QC_Bin'];

% Where the RS are:
fpath2 = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/Signature_Production'];

%Where to save the produtcion
savepath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/Signature_Production'];

% Bins
Nbin=20;
% Ensembles to use
Ens1=10; %First
EnsN=43; %Last
Nens=EnsN-Ens1+1;

% Mean flow direction
alpha=50; % degrees from east to flood (clockwise)

alpha=alpha*pi/180; %radians

% Treshold for mean flow
Uth=0.5; %m/s

% Noise in m/s, check for every deployment
n_sig=0.0265; %m/s

%% Mean flow Data

%
for bin=1:Nbin

    fname = [prefix int2str(bin) '.mat'];        
    load([ fpath '/' fname ])
    
    fs=SigData.fs;
    dt=1/fs;
    
    Ens_time(:,bin)=SigData.Ens_time; %Beggining of each ensemble
    
    [Nt NensT]=size(SigData.vbeam1);
    
    time0=[0:dt:Nt*dt-dt]; % Seconds
    
    fik=1;
    for fi=Ens1:EnsN %Good ensembles
       
        u_east(bin,:,fik)=SigData.u_east(:,fi);
        v_north(bin,:,fik)=SigData.v_north(:,fi);
        w(bin,:,fik)=SigData.vbeam5(:,fi);
        w2(bin,:,fik)=SigData.w_up(:,fi);


        
        timeEns=Ens_time(fi,bin)+time0/3600/24; %days
        timeV(bin,:,fik)=timeEns;
        fik=fik+1;
    end
    
    z=SigData.range;
end

%% Along and Across Channel Data
u=u_east*cos(alpha)-v_north*sin(alpha);
v=u_east*sin(alpha)+v_north*cos(alpha);

for i=1:Nens %Good Ensembles AI May 2015
    
umean(:,i)=nanmean(u(:,:,i),2);
vmean(:,i)=nanmean(v(:,:,i),2);
wmean(:,i)=nanmean(w(:,:,i),2);
wmean2(:,i)=nanmean(w2(:,:,i),2);

umeanaux=repmat(umean(:,i),[1,Nt]);
wmeanaux=repmat(wmean(:,i),[1,Nt]);

% RS approximation from E-N velocities
uprime=u(:,:,i)-umeanaux;
wprime=w(:,:,i)-wmeanaux;

uw=uprime.*wprime;

uwmean(:,i)=nanmean(uw,2);

end

UHmean=sqrt(umean.^2+vmean.^2);

u_depthAv=nanmean(UHmean,1);

u_depthAv=nanmean(umean,1); %%%%Check for every site!!!!!


i_NaN=find(abs(u_depthAv)<Uth);
i_use=find(abs(u_depthAv)>=Uth);
umean(:,i_NaN)=NaN;
uwmean(:,i_NaN)=NaN;

%save UHmean_SigMay2015.mat UHmean

%% Shear and Production

% U_mean

for i=1:Nens
    
    % Keeping Sign of umean
    % For flood dU/dz is positive
    % For ebb dU/dz is negative
    
    dUdz(:,i)=diff(umean(:,i));
    dVdz(:,i)=diff(vmean(:,i));
    dWdz(:,i)=diff(wmean2(:,i));

end

% RS from 5 beam AD2CP

load([ fpath2 '/ReynoldsStress.mat'])

% RS from 5 beam are now in along channel coordinates

%----Don't need this anymore
% %Transform to Along channel velocities
% u_aw_mean=RS.uw*cos(alpha)-RS.vw*sin(alpha);
% v_aw_mean=RS.uw*sin(alpha)+RS.vw*cos(alpha);
% w_aw_mean=RS.ww;
%-------

% Mean production
P_uw_true=-dUdz.*RS.uw_lc_nr(1:end-1,:);
P_vw_true=-dVdz.*RS.vw_lc_nr(1:end-1,:);
P_ww_true=-dWdz.*RS.ww_lc_nr(1:end-1,:);
P_total=P_uw_true+P_vw_true+P_ww_true;

% VT Production
P_uw_vt=-dUdz.*RS.uw_vt_lc_nr(1:end-1,:);
P_vw_vt=-dVdz.*RS.vw_vt_lc_nr(1:end-1,:);
P_total_vt=P_uw_vt+P_vw_vt+P_ww_true;

figure(1)
clf
for i=1:Nens
    semilogx(abs(P_uw_true(:,i)),z(1:end-1),'b')
    hold on
    semilogx(abs(P_uw_vt(:,i)),z(1:end-1),'g')
end
semilogx(nanmean(abs(P_uw_vt),2),z(1:end-1),'g','LineWidth',3)
semilogx(nanmean(abs(P_uw_true),2),z(1:end-1),'b','LineWidth',3)

hold off
xlabel('P')
ylabel('Z (m)')

%% Noise estimates
% Simple noise RS
M_sig=Nt; %Burst length

var_r=n_sig^4/(M_sig*(sin(2*theta*pi/180))^2); %Reynolds stress variance
sigma_r=sqrt(var_r);

% Simple noise Shear
dz=1; %change dz to the one of your measurements 
var_s=n_sig^2/(M_sig*dz*(sin(2*theta))^2);
sigma_s=sqrt(var_s);

Puw_n=RS.uw_lc_nr(1:end-1,:).^2*sigma_s^2+dUdz.^2*sigma_r^2+sigma_s^2*sigma_r^2;
Pvw_n=RS.vw_lc_nr(1:end-1,:).^2*sigma_s^2+dVdz.^2*sigma_r^2+sigma_s^2*sigma_r^2;
Pww_n=RS.ww_lc_nr(1:end-1,:).^2*sigma_s^2+dWdz.^2*sigma_r^2+sigma_s^2*sigma_r^2;

P_n_var=Puw_n+Pvw_n+Pww_n; %Production variance estimate
%% Save Results
Prod.P_uw_true=P_uw_true;
Prod.P_vw_true=P_vw_true;
Prod.P_ww_true=P_ww_true;
Prod.P_total=P_total;

% Variance technique estimates
Prod.P_uw_vt=P_uw_vt;
Prod.P_vw_vt=P_vw_vt;
Prod.P_total_vt=P_total_vt;

% Noise Variance in (m^2s^-3)^2
Prod.P_noise=P_n;

% Mean along-channel velocity
Prod.ua_mean=umean;

savefile=[savepath '/Production.mat'];

save(savefile, 'Prod')



