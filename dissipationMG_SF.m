function [tke epsilon residual A Aerror N Nerror rmax_j epsilon_flat rmax_j_flat eps_flat_error D r] = dissipation_SF(v, z, disspts, plots, deltar,beam,sigma_u); 

 %[tke5(:,fi) epsilon5(:,fi) residual5(:,fi) A5(:,fi) Aerror5(:,fi) N5(:,fi) Nerror5(:,fi) fitstats5(:,fi) rmax5(:,fi)]
% function to estimate dissipation rate [m^2/s^3] from velocity profiles
% using subroutine structure function (Wiles et al 2006)
% and a specificied window length (in points)
%
%   [tke epsilon residual A Aerror N Nerror rmax_j epsilon_flat rmax_j_flat eps_flat_error D r]
%   = dissipation_SF(v, z, disspts, plots, deltar,beam,sigma_u);
%
% where velocity is [bins x time] and z is the [1xbins] array of alongbeam locations 
%   disspts is the number of points to using in sub-window processing (with fixed 50% overlap)
%   plots is binary flag *** NO LONGER IMPLIMENTED *** for plotting the structure function and fits
%   and delta r is [1 x bins] offset for SWIFT motion (tilting and bobbing)
%
%
% J. Thomson, 6/2010, 
%   revs: 9/2010 (allow data gaps as NaNs ) 
%         5/2011 (robust fit, noise tracking)
%         7/2011 (limit length scale)
%         9/2011 (use only forward differencing, r>0)
%         4/2012 (use standard m^2/s^3 units)
%         4/2012 (back to double-sided differencing, |r| )
%         4/2012 (include input arguement to correct r for SWIFT tilting)

% Maricarmen
%11/2015 remove the overlap
%03/2016 add flat D(z,r), try different r and choose the one with the best fit to 
% a flat D(z,r) (with slope closer to zero)
% report two estimates: the one with a fit to a max r directly and one with
% the flat one
% Check for N<<A (Bad N)


% constants
Cvsq = 2.1;


% bins and points
[bins, pts ] = size(v);

lengthresults = 1; %floor(pts/(disspts/overlap));

% intialize results
tke = NaN ( [length(z) lengthresults] );
epsilon = NaN( [length(z) lengthresults]  ); 
A = NaN( [length(z) lengthresults]  ); 
Aerror = NaN( [length(z) lengthresults]  ); 
N = NaN( [length(z) lengthresults]  ); 
Nerror = NaN( [length(z) lengthresults]  ); 
residual = NaN( [length(z) lengthresults] );
rmax_j = NaN([length(z) lengthresults] );
epsilon_flat= NaN( [length(z) lengthresults]  ); 
rmax_j_flat = NaN([length(z) lengthresults] );
eps_flat_error = NaN([length(z) lengthresults] );

ept=1;
% estimate tke
tke(:,ept) = 0.5 * (nanstd( v' ).^2)';

% estimate structure function using entire v
[D r] = structureFunction(v, z);

% figure
% clf
% cmap = colormap;
% clf


for j=1:length(z),
    
    % need to limit r values used to be within inertial subrange
    % ad hoc limit is proportion of depth
    % maxr = max(z)./ 2; %How we define this?
    % alternate limit is distance to boundary:
    % Can set a fixed number or choose best fit
    maxr = z(j); % min( abs(z(j) - 50), z(j) ); %4; %min( abs(z(j) - 50), z(j) ); %Approx Depth.
    % option for double sided fit using abs(r)
    %r = abs(r);

    
    % identify points for fitting, option to exclude r == 0
    goodpts = ~isnan( D(j,:) )  &  r(j,:)< maxr  & r(j,:) > 0;
    
    usepts=find(goodpts==1);
    ruse=r(j,usepts);
    Duse=D(j,usepts);
    Ur = unique(ruse);
    
    % fit structure to r^2/3
    if sum(goodpts)>3,
        warning('off','stats:statrobustfit:IterationLimit')
        
        % Here I try different r max and choose the one with the
        % best fit to 2/3
        % Will also delete bad fits
        
        residualAux=NaN(length(Ur),1);
        mDflat=NaN(length(Ur),1);
        stdDflat=NaN(length(Ur),1);
        
        % Set first 3 values to NaN because I am starting from 4
        % Or use a different counter
        residualAux(1:3)=NaN;
        RMSEfit(1:3)=NaN;
        SlopeDiff(1:3)=NaN;
        slopes(1:3,j)=NaN;
        
        for k=4:length(Ur)
            
            maxr=max(abs(Ur(1:k)));
            rtry=find(abs(ruse)<=maxr);           
            
            %length(rtry)
            %[fit stats] = robustfit(r(j,goodpts).^(2/3), D(j,goodpts));
            [fitAux statsAux] = robustfit(ruse(rtry).^(2/3), Duse(rtry));
            Aaux(k)=fitAux(2);
            Naux(k)=fitAux(1);
            residualAux(k) = statsAux.s;
            
            % Flat D
            Dflat= Duse(rtry).*ruse(rtry).^(-2/3);            
            mDflat(k)=nanmean(Dflat);
            stdDflat(k)=nanstd(Dflat);
            brob = robustfit(ruse(rtry),Dflat);
            mS(k)=brob(2);
            cS(k)=brob(1);
            Dfit_flat=cS(k)+mS(k)*ruse(rtry);
            RMSEfit(k)=sqrt(mean((Dflat-Dfit_flat).^2));
            SlopeDiff(k)=(abs(mS(k))); %Difference from zero slope
            slopes(k,j)=mS(k);          
            
        end
        

        kmin_aux=find(residualAux==min(residualAux));
        kmin=kmin_aux(1);
        
        kmin_flat_aux=find(SlopeDiff==min(abs(SlopeDiff)));
        
        kmin_flat=kmin_flat_aux(1);
        
        maxr_best=max(abs(Ur(1:kmin)));
        
        r_best=find(abs(ruse)<=maxr_best);    
        
        maxr_best_flat=max(abs(Ur(1:kmin_flat)));
        
        r_best_flat=find(abs(ruse)<=maxr_best_flat); 
        
        rmax_j(j)=maxr_best; 
        rmax_j_flat(j)=maxr_best_flat; 
        
        % Flat
        DD_flat=Duse(r_best_flat).*ruse(r_best_flat).^(-2/3);
        brob = robustfit(ruse(r_best_flat),DD_flat);
        m_flat(j)=brob(2);
        c_flat(j)=brob(1);
            
        epsilon_flat(j)=(nanmean(DD_flat)./Cvsq).^(3/2);
        eps_flat_error(j)=3/2*nanmean(DD_flat)^0.5/Cvsq.^(3/2)*nanstd(DD_flat);
        
        % Clear for slope to far from zero
        if SlopeDiff(kmin_flat)>0.1
            epsilon_flat(j)=NaN;
            eps_flat_error(j)=NaN;
        end
        
        % For use when not using the r_best finder, just a fixed rbest, 
        
        [fit stats] = robustfit(ruse(r_best).^(2/3), Duse(r_best));
        A(j,ept) = fit(2);
        N(j,ept) = fit(1);
        fitstats(j,ept) = stats;
        residual(j,ept) = stats.s;  %rmse
        Aerror(j,ept) = stats.se(2);
        Nerror(j,ept) = stats.se(1);
        Ar=A(j,ept)*ruse(r_best).^(2/3);
        
        % Check for bad N
        badN=find(Ar<=N(j,ept));
        if ~isempty(badN)
            A(j,ept)=NaN;
            epsilon_flat(j)=NaN;
        end
        
        if N>(2*sigma_u.^2)
            A(j,ept)=NaN;
            epsilon_flat(j)=NaN;
        end
        
            
        
        clear badN

%         figure(j)
%         clf 
%         cmap=colormap;
%         
%         subplot(2,1,1)
%         %plot(real(r(j,goodpts)),D(j,goodpts),'*','MarkerSize',10,'color',cmap(floor(j/20*64),:))
%         plot(ruse(r_best).^(2/3),Duse(r_best),'*','MarkerSize',20,'color',cmap(floor(j/20*64),:))
%         hold on
%         plot(ruse(r_best_flat).^(2/3),Duse(r_best_flat),'o','MarkerSize',20,'color',cmap(floor(j/20*64),:))
%         plot(r(j,r(j,:)>0).^(2/3),D(j,r(j,:)>0),'k*','MarkerSize',5)
%         %hold on
%         %plot(real(ruse(r_best)),A(j,ept)*real(ruse(r_best)).^(2/3)+ N(j,ept),'color',cmap(floor(j/20*64),:),'Linewidth',2)
%         %plot(real(r(j,goodpts)),A(j,ept)*real(r(j,goodpts)).^(2/3)+ N(j,ept),'color',cmap(floor(j/20*64),:),'Linewidth',2)
%         plot(abs(r(j,:)).^(2/3),A(j,ept)*abs(r(j,:)).^(2/3)+ N(j,ept),'color',cmap(floor(j/20*64),:),'Linewidth',2)
%         hold off
%         xlabel('r^{(2/3)} (m)')
%         ylabel('D(z,r)')
%         title(['z= ' int2str(z(j)) ''])
%         set(gca,'FontSize',16)
%         
%         subplot(2,1,2) %Flats
%         plot(ruse(r_best),Duse(r_best).*ruse(r_best).^(-2/3),'*','MarkerSize',20,'color',cmap(floor(j/20*64),:))
%         hold on
%         plot(ruse(r_best_flat),Duse(r_best_flat).*ruse(r_best_flat).^(-2/3),'o','MarkerSize',20,'color',cmap(floor(j/20*64),:))
%         plot(r(j,r(j,:)>0),D(j,r(j,:)>0).*r(j,r(j,:)>0).^(-2/3),'k*','MarkerSize',5)
%         %hold on
%         %plot(real(ruse(r_best)),A(j,ept)*real(ruse(r_best)).^(2/3)+ N(j,ept),'color',cmap(floor(j/20*64),:),'Linewidth',2)
%         %plot(real(r(j,goodpts)),A(j,ept)*real(r(j,goodpts)).^(2/3)+ N(j,ept),'color',cmap(floor(j/20*64),:),'Linewidth',2)
%         plot(abs(r(j,:)),(A(j,ept)*abs(r(j,:)).^(2/3)+ N(j,ept)).*abs(r(j,:)).^(-2/3),'color',cmap(floor(j/20*64),:),'Linewidth',2)
%         hold off
%         xlabel('r (m)')
%         ylabel('D(z,r)*r^{-2/3}')
%         title(['z= ' int2str(z(j)) ''])
%         set(gca,'FontSize',16)
%         
%         pause
        
       

    else
    end
    

    
    clear ruse Duse Ur kmin_flat kmin_flat_aux kmin kmin_aux residualAux SlopeDiff
    
end

% hold off
% xlabel('r (m)')
% ylabel('D (m^2/s^2)')
% saveas(gca,['StructureFunction_Beam' int2str(beam) ''],'jpg')

%end


% dissipation rate [m^2/s^3], switched from W/m^3 on 4/6/2012
posA = find(A>0);
negA=find(A<=0);

%rho = 1024;
epsilon(posA) = ( A(posA) ./ Cvsq ).^(3/2) ; % m^2/s^3
epsilon_flat(negA)=NaN;


