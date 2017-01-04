% Maricarmen Guerra Paris

% Open .mat files created with MIDAS Software from Nortek Signature 5 beam
% Transforms to ENU and XYZ velocities
% Save new files containing beam, ENU and XYZ velocities.
% % Save as: '/SS04_Sig_May2015_ENU_00000_' int2str(fnumber) '.mat'];
% Repeat for all MIDAS files
% Plots raw data


clear all, close all, clc

% Raw data files:
fpath = '/Users/Maru/Documents/PhD_UW/AdmiraltyInlet/SignatureMay2015'; % File location
fnumber=2; %Number of files 
fname = ['SS04_Sig_May2015.ad2cp.133.15.AD2CPRaw.00000_' int2str(fnumber) '.mat']; %Sequential name

% Save Files
prefix = ['SS04_Sig_May2015_ENU_00000_']; %where to save the files
savepath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/RawData']; %Sequential name of files

%% Raw Data
load([ fpath '/' fname ])
figure(1), clf

for i = 1:4,
    
    ax(i) = subplot(5,1,i); 
    set(gca,'FontSize',16)
    pcolor(Data.Burst_MatlabTimeStamp-datenum(2015,0,0), double(Data.Burst_Range), double(eval(['Data.Burst_VelBeam' num2str(i)]))' ), 
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Vel Beam ' num2str(i)])
    caxis([-2 2])
    colorbar,
   
    
end

ax(5) = subplot(5,1,5);
 set(gca,'FontSize',16)
pcolor(Data.IBurst_MatlabTimeStamp-datenum(2015,0,0), double(Data.IBurst_Range), double(Data.IBurst_VelBeam5)' ), 
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Vel Beam 5'])
    caxis([-2 2])
    colorbar,

linkaxes(ax,'x')


%% Convert to ENU coordinates

[ Data2, Config2, T_beam2xyz ]=signatureAD2CP_beam2xyz_enu(Data,Config,'burst')

% Save data with ENU coordinates
savefile=[savepath '/' prefix int2str(fnumber) '.mat'];
save(savefile, '-mat', 'Data2');


%% ENU Data
figure(2), clf


    
    ax(1) = subplot(3,1,1); 
    
    pcolor(Data2.Burst_MatlabTimeStamp-datenum(2014,0,0), double(Data2.Burst_Range), double(Data2.Burst_VelX)' ), 
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Vel X'])
    caxis([-3 3])
    colorbar,
      ax(2) = subplot(3,1,2); 
    
    pcolor(Data2.Burst_MatlabTimeStamp-datenum(2014,0,0), double(Data2.Burst_Range), double(Data2.Burst_VelY)' ), 
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Vel Y'])
    caxis([-3 3])
    colorbar,

ax(3) = subplot(3,1,3);
pcolor(Data2.Burst_MatlabTimeStamp-datenum(2014,0,0), double(Data2.Burst_Range), double(Data2.Burst_VelZ)' ), 
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Vel Z'])
    caxis([-2 2])
    colorbar,

linkaxes(ax,'x')

%% Battery life
figure(3)
plot(Data.Burst_MatlabTimeStamp-datenum(2014,0,0),Data.Burst_Battery)


