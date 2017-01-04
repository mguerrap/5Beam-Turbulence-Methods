% Signature Low Correlation and Low Amp QC

clc
clear all
close all

% Define Min correlation and min amplitude

mincor=50;
minamp=30;

% Number of bins
Nbin=20;

% Ensembles to QC
Ens1=10;
EnsN=43;
NEns=EnsN-Ens1+1;

% File path: where are the bin files
fpath = ['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];
prefix='SignatureData_Bin';

% Save path: where to save the QC data
savepath=['/Users/Maru/Documents/PhD_UW/5BeamCodes/BinDataSignature'];

% Raw Data
for bin=1:Nbin

    
    bin
    
    fname = [prefix int2str(bin) '.mat'];
    
    load([ fpath '/' fname ])
    

    
    % Only for good ensembles: 10 to 43
    
    for fi=Ens1:EnsN
        
        b1=SigData.vbeam1(:,fi);
        b2=SigData.vbeam2(:,fi);
        b3=SigData.vbeam3(:,fi);
        b4=SigData.vbeam4(:,fi);
        b5=SigData.vbeam5(:,fi);
        
        NBurst=length(b1);
        
        c1=SigData.Corbeam1(:,fi);
        c2=SigData.Corbeam2(:,fi);
        c3=SigData.Corbeam3(:,fi);
        c4=SigData.Corbeam4(:,fi);
        c5=SigData.Corbeam5(:,fi);
        
        a1=SigData.Ampbeam1(:,fi);
        a2=SigData.Ampbeam2(:,fi);
        a3=SigData.Ampbeam3(:,fi);
        a4=SigData.Ampbeam4(:,fi);
        a5=SigData.Ampbeam5(:,fi);

        
        bad = find( c1 < mincor | c2 < mincor | c3 < mincor | c4 < mincor | c5 < mincor |...
            a1 < minamp | a2 < minamp | a3 < minamp | a4 < minamp | a5 < minamp); 
        %bad = find(isnan(burstu)==1); 

        all_data=find(c1); %get indices of all data
        good_data=find(ismember(all_data,bad)==0); %finds the indices that are not in bad, i.e. the good ones

        
        largobad(bin,fi)=length(bad);
        largogood(bin,fi)=length(good_data);
    
        if length(bad) < 0.3 * length(b1)
            
        %burstu(bad) = nanmean(burstu); burstv(bad) = nanmean(burstv); burstw(bad) = nanmean(burstw);
        
        SigData.vbeam1(bad,fi)= nanmean(b1(good_data)); 
        SigData.vbeam2(bad,fi)= nanmean(b2(good_data)); 
        SigData.vbeam3(bad,fi)= nanmean(b3(good_data)); 
        SigData.vbeam4(bad,fi)= nanmean(b4(good_data)); 
        SigData.vbeam5(bad,fi)= nanmean(b5(good_data)); 
        
        SigData.u_x(bad,fi)=nanmean(SigData.u_x((good_data),fi));
        SigData.v_y(bad,fi)=nanmean(SigData.v_y((good_data),fi));
        SigData.w_z(bad,fi)=nanmean(SigData.w_z((good_data),fi));
        
        SigData.u_east(bad,fi)=nanmean(SigData.u_east((good_data),fi));
        SigData.v_north(bad,fi)=nanmean(SigData.v_north((good_data),fi));
        SigData.w_up(bad,fi)=nanmean(SigData.w_up((good_data),fi));
       
       % Replace NaN data with an interpolation
       
        else
        %burstu = NaN(size(burstu)); burstv = NaN(size(burstv)); burstw = NaN(size(burstw));
        SigData.vbeam1(bad,fi)= NaN(size(b1)); 
        SigData.vbeam2(bad,fi)= NaN(size(b1)); 
        SigData.vbeam3(bad,fi)= NaN(size(b1)); 
        SigData.vbeam4(bad,fi)= NaN(size(b1));  
        SigData.vbeam5(bad,fi)= NaN(size(b1));  
        
        SigData.u_x(bad,fi)=NaN(size(b1)); 
        SigData.v_y(bad,fi)=NaN(size(b1)); 
        SigData.w_z(bad,fi)=NaN(size(b1)); 
        
        SigData.u_east(bad,fi)=NaN(size(b1)); 
        SigData.v_north(bad,fi)=NaN(size(b1)); 
        SigData.w_up(bad,fi)=NaN(size(b1)); 
        end
    end
        

    save([savepath '/' 'SignatureData_QC_Bin' int2str(bin) '.mat'], 'SigData')
    
    
    clear bad
end


%% Bad data percentage

BadTotal=sum(sum(largobad));
DataTotal=Nbin*NBurst*NEns;

PercentBad=BadTotal/DataTotal*100;
        

