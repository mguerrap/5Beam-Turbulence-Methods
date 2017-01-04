% Signature 5-Beam Reynolds Stresses

% Uses RS_5beam.m to estimate the RS
% RS are in Earth coordinates using instrument heading

clc
clear
close all

% Data: 
% Needs along beam velocities for all depths as bi=bi(z,time) for each burst
% Needs pitch and roll for each burst 
% Needs u and v also in u(z,time) form for each burst
% We have data saved for each time burst for each bin, we need
% to open each file and form b(depth,time);

% Where is the QC Data:
fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefix = ['SignatureData_QC_Bin'];

%Where to safe the RS
savepath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/Signature_Production'];

% Bins
Nbin=20;
% Ensembles to use
Ens1=10; %First
EnsN=43; %Last

% Theta: beam slanted angle
theta=25; %Degrees

% Noise in m/s
n_sig=0.0265; % Check noise for every deployment



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
        
        u_x(bin,:,fik)=SigData.u_x(:,fi);
        v_y(bin,:,fik)=SigData.v_y(:,fi);
        
         % Need Pitch and Roll and Heading
         pitch(:,fik)=SigData.Pitch(:,fi);
         roll(:,fik)=SigData.Roll(:,fi);
         heading(:,fik)=SigData.Heading(:,fi);
    
        timeEns=Ens_time(fi,bin)+time0/3600/24; %days
        timeV(bin,:,fik)=timeEns;
        
        fik=fik+1;
        
    end
    
   
    z=SigData.range;
end

%% Reynolds Stresses in time

[Nz Nt Nens]=size(timeV);

% estimate RS [m2/s2] from beam data

for fi=1:Nens
    b1=v1(:,:,fi);
    b2=v2(:,:,fi);
    b3=v3(:,:,fi);
    b4=v4(:,:,fi);
    b5=v5(:,:,fi);
    
    u=u_x(:,:,fi);
    v=v_y(:,:,fi);
    
    Uh(:,fi)=nanmean(sqrt(u.^2+v.^2),2);
    
    r=roll(:,fi);
    p=pitch(:,fi);
    h=heading(:,fi);
    
    % In here we call the Reynolds Stress function from Dewey & Stringer
    % 2007, we need to be consistent with the beams and the tilt angles:
    
    [uu(:,fi) vv(:,fi) ww(:,fi) uw(:,fi) vw(:,fi) uv(:,fi) anisotropy(:,fi) q2(:,fi)]=...
        RS_5beam(b1,b3,b4,b2,b5,theta,h,-r,p,u,v);
     
    %pause
    fi
end

badRS = find( uu<=0 | vv<=0 | ww<=0);
uu(badRS)=NaN;
vv(badRS)=NaN;
ww(badRS)=NaN;
uw(badRS)=NaN;
uv(badRS)=NaN;
vw(badRS)=NaN;

%% Reynolds Stress Tensor
for i=1:Nens
    for j=1:Nbin
        Tij=[uu(j,i) uv(j,i) uw(j,i); uv(j,i) vv(j,i) vw(j,i) ; uw(j,i) vw(j,i) ww(j,i)];
        
        % IF the tensor has no NaN values
        if ~isnan(Tij(1,1)) & ~isnan(Tij(2,1)) & ~isnan(Tij(3,1)) & ~isnan(Tij(2,2)) & ~isnan(Tij(3,2))
            % Estimade Eigenvalues
            EV=eig(Tij);
            
            % If any eigenvalue is negative, then tensor is not positive
            % definite
            if (EV(1)<0 | EV(2)<0 | EV(3)<0)
                % Tensor is not positive definite
                badTS(j,i)=NaN;
                
            else
                badTS(j,i)=0;
            end
        else
            badTS(j,i)=NaN;
        end
    end
end

% Make Reynolds Stresses NaN for bad tensor:
badT=find(isnan(badTS)==1);
uu(badT)=NaN;
vv(badT)=NaN;
ww(badT)=NaN;
uw(badT)=NaN;
uv(badT)=NaN;
vw(badT)=NaN;

%% Variance Technique
for fi=1:Nens
    b1=v1(:,:,fi);
    b2=v2(:,:,fi);
    b3=v3(:,:,fi);
    b4=v4(:,:,fi);
    b5=v5(:,:,fi);
    h=heading(:,fi);
    [uw_vt(:,fi) vw_vt(:,fi)]= RS_VT(b1,b3,b4,b2,25,h);
     
    %pause
    fi
end



%% Save results into RS structure
RS.uu=uu;
RS.vv=vv;
RS.ww=ww;
RS.uw=uw;
RS.vw=vw;
RS.uv=uv;
RS.uw_vt=uw_vt;
RS.vw_vt=vw_vt;

%% Local coordinates Rynolds Stresses
alpha=50; % degrees from east to flood clockwise; check this for every channel!
alpha=alpha*pi/180; %radians

RS.uu_lc=RS.uu*cos(alpha)^2+RS.vv*sin(alpha)^2-2*RS.uv*sin(alpha)*cos(alpha);
RS.vv_lc=RS.uu*sin(alpha)^2+RS.vv*cos(alpha)^2+2*RS.uv*sin(alpha)*cos(alpha);
RS.uw_lc=RS.uw*cos(alpha)-RS.vw*sin(alpha);
RS.vw_lc=RS.uw*sin(alpha)+RS.vw*cos(alpha);
RS.ww_lc=RS.ww;
RS.uv_lc=(RS.uu-RS.vv)*sin(alpha)*cos(alpha)+RS.uv*(cos(alpha)^2-sin(alpha)^2);

RS.uw_vt_lc=RS.uw_vt*cos(alpha)-RS.vw_vt*sin(alpha);
RS.vw_vt_lc=RS.uw_vt*sin(alpha)+RS.vw_vt*cos(alpha);

% TKE
RS.k=q2; %Total Kinetic energy
RS.k_sum=0.5*(RS.uu+RS.vv+RS.ww);
RS.k_sum_lc=0.5*(RS.uu_lc+RS.vv_lc+RS.ww_lc);
RS.anisotropy=anisotropy;

%% Noise removal
M_sig=Nt; %burst length
var_sig=n_sig^4/(M_sig*(sin(2*theta*pi/180))^2);
sigma_sig=sqrt(var_sig);

% Remove noise:

RS.uu_lc_nr=RS.uu_lc-sigma_sig;
RS.vv_lc_nr=RS.vv_lc-sigma_sig;
RS.ww_lc_nr=RS.ww_lc-sigma_sig;
RS.uw_lc_nr=RS.uw_lc-sigma_sig;
RS.vw_lc_nr=RS.vw_lc-sigma_sig;
RS.uv_lc_nr=RS.uv_lc-sigma_sig;

RS.uw_vt_lc_nr=RS.uw_vt_lc-sigma_sig;
RS.vw_vt_lc_nr=RS.vw_vt_lc-sigma_sig;

%% Save File
savefile=[savepath '/ReynoldsStress.mat'];

save(savefile, 'RS')


%% Plots
cmap=colormap;
figure(1)
clf
plot(nanmean(uw,2),z,'r','LineWidth',2)
hold on
%semilogx(nanmean(-uw_simp,2),z,'g','LineWidth',2)
for i=1:Nens
cindex=(floor(nanmean(Uh(:,i))/2*64));
plot(uw(:,i),z,'color',cmap(cindex,:))
end
hold off
xlabel('uw') 
ylabel('Z (m)')

figure(2)
clf
subplot(1,2,1)
semilogx(nanmean(abs(RS.uw_lc_nr),2),z,'r','LineWidth',2)
hold on
semilogy(nanmean(abs(RS.uw_vt_lc_nr),2),z,'b','LineWidth',2)
hold off
xlabel('uw') 
ylabel('Z (m)')
legend('5 beam','VT')
axis([10^-3 10^-2 0 20])
subplot(1,2,2)
semilogx(nanmean(abs(RS.vw_lc_nr),2),z,'r','LineWidth',2)
hold on
semilogx(nanmean(abs(RS.vw_vt_lc_nr),2),z,'b','LineWidth',2)
hold off
xlabel('vw') 
ylabel('Z (m)')
legend('5 beam','VT')
axis([10^-3 10^-2 0 20])
