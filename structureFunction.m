function [ D r ] = structureFunction( v , z );
%
% Function to calculate the spatial structure function of an ADCP time series (burst w/stationarity)
%
%   [ D r ] = structureFunction( v , z );
%
% Given a velocity profile that is [bins, time] and a corresponding vector of bin heights, 
% the function returns structure and range arrays that are [bins, range]   
% 
%
% J. Thomson, 7/2009, 
%       rev. 7/2010 (efficiency) 
%       rev. 8/2010 (explicitly remove mean [should be done already], allow nans in profiles) 
%       rev. 9/2011  (return signed r value, to allow upward or downward preference in D(z,r) fit)
% M. Guerra Dec 8 2015, change r deffinition to uplooking

[bins time] = size(v);

v = v - nanmean(v,2)*ones(1,time); % Remove mean

for i = [1:bins],
  
    r(i,:) = z(:) - z(i); % For uplooking z1 is the smallest v(z+r)-v(z), then r=(z+r)-z;
    %For downlooking r(i,:) = z(i) - z(:); 
    
        for t = 1:time, 
            % Here the sign of the dzrt is not important since it is
            % squared: v(z+r)-v(z)
            
        dzrt(:,t) = ( v(:,t) - v(i,t) ).^2; % do the same for velocities for each time
        
        end
        
    D(i,:) = nanmean(dzrt,2);       % Ojo que da vuelta el vector nanmean(dzrt,2) para tener D(z,r)
    
end
