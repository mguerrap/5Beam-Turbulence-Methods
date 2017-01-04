function [uwt vwt]=RS_VT(b1,b2,b3,b4,th,phi1)

% Reynolds Stress for a 4 beam ADCP using the variance technique

% bi= beams 1 to 4
% phi1=Heading
% phi2=Pitch
% phi3=Roll
% theta=beam angle
% [Nz Nt]=[N Bins, Length Time]
% Length in Time must correspond to time for average
% Give burst matrix (5 or 10 minutes for the Signature)
% All angles in degrees

% For the Nortek Signature 1000
% Dewey/Nortek
% b1=b1
% b2=b3
% b3=b4
% b4=b2
% b5=b5
% theta=beam angle, 25º
% Call function like this:
% RS_VT(b1,b3,b4,b2,25)

[Nz,Nt]=size(b1);

%% Transform to rads
th=th*pi/180;
phi1=phi1*pi/180;

phi1=nanmean(phi1);

%% Get beam velocities fluctuations
b1=b1-repmat(nanmean(b1,2),[1,Nt]);
b2=b2-repmat(nanmean(b2,2),[1,Nt]);
b3=b3-repmat(nanmean(b3,2),[1,Nt]);
b4=b4-repmat(nanmean(b4,2),[1,Nt]);

%% Square and take time average
b1_av=nanmean(b1.^2,2);
b2_av=nanmean(b2.^2,2);
b3_av=nanmean(b3.^2,2);
b4_av=nanmean(b4.^2,2);

%% u'w'
uw=-1/(2*sin(2*th)).*(b2_av-b1_av);

%% v'w'
vw=-1/(2*sin(2*th)).*(b4_av-b3_av);

%% Convert to Earth Coordinates using heading
uwt=uw*cos(phi1)+vw*sin(phi1);
vwt=vw*cos(phi1)-uw*sin(phi1);

end

