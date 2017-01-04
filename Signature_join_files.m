% Signature Join .mat files
% This file joins Nortek Signature .mat output from the
% Signature_RawData_ENU.m script

% Since the MIDAS software puts data in separate files, the ENU files are
% also separate.

% This code: 
% 1. Puts all data from each bin in a single file to work with each bin
% separately
% Saves into a .mat file

% This script only works because Nortek output .mat files do not cut
% ensembles. Meaning that each Nortek file has complete ensembles.

clear all, close all, clc

% Introduce number of bins:
Nbin=20;

% Introduce number of files from MIDAS software
fn=3;

% Introduce single path for were data is located:
% We need the data files previously processed to get ENU velocity:
% The file must have the structs Data1 and Data2

fpath=['/Users/Maru/Documents/PhD_UW/5BeamCodes/RawData'];
% Introduce prefix of files
prefix=['SS04_Sig_May2015_ENU_00000_'];

% Path to save files with bin data

savepath=['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];

% Sampling frequency
fs=8; 
    
for bin=1:Nbin
    
    % Ensemble number
    Ens_N=[];
    
    for j=1:fn
                
        fnumber=j;
        
        fname = [prefix int2str(fnumber) '.mat'];
        
        
        load([ fpath '/' fname ])
        
        % Time and Ensemble Size
        time=Data2.Burst_MatlabTimeStamp;
        date=datevec(time);
        Burst_length=1; %10 minutes
        
        % Ensemble treshold
        
        % Modify depending on the sampling frequency and in 
        % Delta time in days use to separate burst
        % In our case it was 10 mins between every burst, so a dt greater
        % than 2/fs should cut out the burst in the right place, but CHECK!
        % Also needs to modify if instrument is recording continuously
        
        dt_days=2*1/fs/3600/24; %Modified because sometimes time step is slightly larger than 1/dt sec
        dt_tstamps=diff(time);
        k=1;
        
        for i=1:length(dt_tstamps)
            if dt_tstamps(i)>dt_days; % Any dt greater than dt_days marks a new ensemble
                % New Ensemble begins
                ens_in(k)=i;
                ens_end(k)=i-1;
                k=k+1;
            else
            end
        end
        
        ens_length=diff(ens_in);
        
        % Good data to start:
        % Here we start with the first burst, even if some burst are out of
        % the water (we just do not use them later)
        
        ens_start=1;
        tstart=1;
        for i=1:ens_start-1;
            tstart=ens_length(i)+tstart;
        end
        
        Nens0=length(ens_length);
        tstart(ens_start)=tstart;
        
        % Everything useful is stored in the structure Sig.
        
        for i=ens_start:Nens0
            tend(i)=tstart(i)+ens_length(i)-1;
            Sig(j).vbeam1(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelBeam1(tstart(i):tend(i),bin);
            Sig(j).vbeam2(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelBeam2(tstart(i):tend(i),bin);
            Sig(j).vbeam3(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelBeam3(tstart(i):tend(i),bin);
            Sig(j).vbeam4(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelBeam4(tstart(i):tend(i),bin);
            Sig(j).vbeam5(1:ens_length(i),i-ens_start+1)=Data2.IBurst_VelBeam5(tstart(i):tend(i),bin);
            Sig(j).u_east(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelEast(tstart(i):tend(i),bin);
            Sig(j).v_north(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelNorth(tstart(i):tend(i),bin);
            Sig(j).w_up(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelUp(tstart(i):tend(i),bin);
            
            Sig(j).u_x(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelX(tstart(i):tend(i),bin);
            Sig(j).v_y(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelY(tstart(i):tend(i),bin);
            Sig(j).w_z(1:ens_length(i),i-ens_start+1)=Data2.Burst_VelZ(tstart(i):tend(i),bin);
            
            Sig(j).heading(1:ens_length(i),i-ens_start+1)=Data2.Burst_Heading(tstart(i):tend(i));
            Sig(j).pitch(1:ens_length(i),i-ens_start+1)=Data2.Burst_Pitch(tstart(i):tend(i));
            Sig(j).roll(1:ens_length(i),i-ens_start+1)=Data2.Burst_Roll(tstart(i):tend(i));
            
            % Correlation and Amplitude for QC
            
            Sig(j).Corbeam1(1:ens_length(i),i-ens_start+1)=Data2.Burst_CorBeam1(tstart(i):tend(i),bin);
            Sig(j).Corbeam2(1:ens_length(i),i-ens_start+1)=Data2.Burst_CorBeam2(tstart(i):tend(i),bin);
            Sig(j).Corbeam3(1:ens_length(i),i-ens_start+1)=Data2.Burst_CorBeam3(tstart(i):tend(i),bin);
            Sig(j).Corbeam4(1:ens_length(i),i-ens_start+1)=Data2.Burst_CorBeam4(tstart(i):tend(i),bin);
            Sig(j).Corbeam5(1:ens_length(i),i-ens_start+1)=Data2.IBurst_CorBeam5(tstart(i):tend(i),bin);
            
            Sig(j).Ampbeam1(1:ens_length(i),i-ens_start+1)=Data2.Burst_AmpBeam1(tstart(i):tend(i),bin);
            Sig(j).Ampbeam2(1:ens_length(i),i-ens_start+1)=Data2.Burst_AmpBeam2(tstart(i):tend(i),bin);
            Sig(j).Ampbeam3(1:ens_length(i),i-ens_start+1)=Data2.Burst_AmpBeam3(tstart(i):tend(i),bin);
            Sig(j).Ampbeam4(1:ens_length(i),i-ens_start+1)=Data2.Burst_AmpBeam4(tstart(i):tend(i),bin);
            Sig(j).Ampbeam5(1:ens_length(i),i-ens_start+1)=Data2.IBurst_AmpBeam5(tstart(i):tend(i),bin);
            
        
            tstart(i+1)=tend(i)+1;
            
        end

        Sig(j).pressure=Data2.Burst_Pressure;
        Sig(j).range=Data2.Burst_Range;
        Sig(j).time=Data2.Burst_MatlabTimeStamp;
        
        % Time ensemble start
        Sig(j).Ens_time=Data2.Burst_MatlabTimeStamp(tstart(1:end-1));       
        Sig(j).Iheading=Data2.IBurst_Heading;
        Sig(j).Ipitch=Data2.IBurst_Pitch;
        Sig(j).Iroll=Data2.IBurst_Roll;
        Sig(j).Ipressure=Data2.IBurst_Pressure;
        Sig(j).Irange=Data2.IBurst_Range;
        Sig(j).Itime=Data2.IBurst_MatlabTimeStamp;
        
        Sig(j).IEns_time=Data2.IBurst_MatlabTimeStamp(tstart(1:end-1));
        
        Nens(j)=Nens0;
        Ens_N=[Ens_N ens_length];
        
        
    end
    
    % Use this figure to check the size of the burst (or ensembles)
    % When battery is low, the Nortek Signature decreases ensemble lenght
    
%     figure(1)
%     clf
%     set(gca,'FontSize',14)
%     plot([1:1:length(Ens_N)],Ens_N)
%     xlabel('ensemble')
%     ylabel('Ensemble Length')
    
    
    % Then we save in only one mat_file for each bin:
    % Choose a different name if you want, here we chose SigData
    
    SigData.vbeam1=[];
    SigData.vbeam2=[];
    SigData.vbeam3=[];
    SigData.vbeam4=[];
    SigData.vbeam5=[];
    
    SigData.u_east=[];
    SigData.v_north=[];
    SigData.w_up=[];
    
    SigData.u_x=[];
    SigData.v_y=[];
    SigData.w_z=[];
    
    SigData.time=[];
    SigData.Itime=[];
    SigData.Ens_time=[];
    SigData.IEns_time=[];
    
    SigData.Heading=[];
    SigData.Roll=[];
    SigData.Pitch=[];
    SigData.Pressure=[];
    
    SigData.IHeading=[];
    SigData.IRoll=[];
    SigData.IPitch=[];
    SigData.IPressure=[];
    
    SigData.Corbeam1=[];
    SigData.Corbeam2=[];
    SigData.Corbeam3=[];
    SigData.Corbeam4=[];
    SigData.Corbeam5=[];
    
    SigData.Ampbeam1=[];
    SigData.Ampbeam2=[];
    SigData.Ampbeam3=[];
    SigData.Ampbeam4=[];
    SigData.Ampbeam5=[];
    
    
    for j=1:fn
        SigData.vbeam1=[SigData.vbeam1 Sig(j).vbeam1];
        SigData.vbeam2=[SigData.vbeam2 Sig(j).vbeam2];
        SigData.vbeam3=[SigData.vbeam3 Sig(j).vbeam3];
        SigData.vbeam4=[SigData.vbeam4 Sig(j).vbeam4];
        SigData.vbeam5=[SigData.vbeam5 Sig(j).vbeam5];
        
        SigData.Corbeam1=[SigData.Corbeam1 Sig(j).Corbeam1];
        SigData.Corbeam2=[SigData.Corbeam2 Sig(j).Corbeam2];
        SigData.Corbeam3=[SigData.Corbeam3 Sig(j).Corbeam3];
        SigData.Corbeam4=[SigData.Corbeam4 Sig(j).Corbeam4];
        SigData.Corbeam5=[SigData.Corbeam5 Sig(j).Corbeam5];
        
        SigData.Ampbeam1=[SigData.Ampbeam1 Sig(j).Ampbeam1];
        SigData.Ampbeam2=[SigData.Ampbeam2 Sig(j).Ampbeam2];
        SigData.Ampbeam3=[SigData.Ampbeam3 Sig(j).Ampbeam3];
        SigData.Ampbeam4=[SigData.Ampbeam4 Sig(j).Ampbeam4];
        SigData.Ampbeam5=[SigData.Ampbeam5 Sig(j).Ampbeam5];
        
        SigData.u_east=[SigData.u_east Sig(j).u_east];
        SigData.v_north=[SigData.v_north Sig(j).v_north];
        SigData.w_up=[SigData.w_up Sig(j).w_up];
        
        SigData.u_x=[SigData.u_x Sig(j).u_x];
        SigData.v_y=[SigData.v_y Sig(j).v_y];
        SigData.w_z=[SigData.w_z Sig(j).w_z];
        
        
        SigData.time=[SigData.time; Sig(j).time];
        SigData.Itime=[SigData.Itime; Sig(j).Itime];
        SigData.Ens_time=[SigData.Ens_time; Sig(j).Ens_time];
        SigData.IEns_time=[SigData.IEns_time; Sig(j).IEns_time];
        
        SigData.Heading=[SigData.Heading Sig(j).heading];
        SigData.Roll=[SigData.Roll Sig(j).roll];
        SigData.Pitch=[SigData.Pitch Sig(j).pitch];
        
        SigData.Pressure=[SigData.Pressure; Sig(j).pressure];
        
        SigData.IHeading=[SigData.IHeading; Sig(j).Iheading];
        SigData.IRoll=[SigData.IRoll; Sig(j).Iroll];
        SigData.IPitch=[SigData.IPitch; Sig(j).Ipitch];
        SigData.IPressure=[SigData.IPressure; Sig(j).Ipressure];
        
        
    end
    
    SigData.fs=8;
    SigData.Irange=[Sig(1).Irange];
    SigData.range=[Sig(1).range];
    SigData.Nens=Nens;
    SigData.Ens_size=Ens_N;
    
    
    save([savepath '/' 'SignatureData_Bin' int2str(bin) '.mat'], 'SigData')
    
    
    bin
    length(Ens_N)
    clearvars -except bin fn fpath prefix savepath fs Nbin
    
end


