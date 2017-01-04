% Signature File created with Signature_QC

% This script takes ENU velocities and estimates horizontal velocities and
% give them a sign
% It uses the function sign_speed.m taken from Polagye & Thomson (2013) Matlab codes


% Plot Raw Data 
clc
clear all
close all

tic

theta=50; %flood heading from east clockwise
Nbin=20;

% File Paths and prefix
% Where is the QC Data:
fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefix = ['SignatureData_QC_Bin'];

%Where to safe the spectra
savepath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
    
%------------------------------


flood_heading=theta+90; %FROM NORTH

flood_heading=flood_heading + [-90, +90];
theta=theta*pi/180; %FROM EAST

for bin=1:Nbin

    fname = [prefix int2str(bin) '.mat'];        
    load([ fpath '/' fname ])
       
    [Nt Nens]=size(SigData.vbeam1);
    fs=SigData.fs;
    
    dt=1/fs;
    Ens_size=SigData.Ens_size;
    
    Ens_time=SigData.Ens_time;
    time=SigData.time;
    
    u=reshape(SigData.u_east,Nt*Nens,1);
    v=reshape(SigData.v_north,Nt*Nens,1);
    w=reshape(SigData.w_up,Nt*Nens,1);

    
    U=sqrt(u.^2+v.^2);
    Udir=atan2(u,v)*180/pi;
    
    Usigned=sign_speed(u, v, U, Udir, flood_heading);
    
    
    Ens_Usigned=reshape(Usigned,Nt,Nens);
    
    savefile=['Signature_Usigned_Bin' int2str(bin) '.mat'];
    
    save([savepath '/' savefile],'Ens_Usigned');

end
