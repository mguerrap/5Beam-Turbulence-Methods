function [u2t v2t w2t uwt vwt uvt alpha q2]=RS_5beam(b1,b2,b3,b4,b5,th,phi1,phi2,phi3,u,v)

% -----------------------
% Modified May 25th 2017
% Sign of pitch and roll for the Signature/Sentinel in description bellow
% ------------------------
% Modified Jan 3rd 2017
% Reynolds Stress for a 5 beam ADCP for small angel approx for Pitch and
% Roll

% From Dewey & Stringer (2007)
% According to Dewey (2007) 5 beam Janus ADCP configuration:
% bi= beams 1 to 5
% phi1=Heading
% phi2=Pitch
% phi3=Roll
% theta=beam angle
% u=u_x (instrument coordinates)
% v=v_y (instrument coordinates)
% u and v are used to estimate u'v'
% [Nz Nt]=[N Bins, Length Time]
% Lenth in Time must correspond to time for average
% Give burst matrix (5 or 10 minutes for the Signature)
% All angles in degrees


% New May 25th 2017
% For the Nortek Signature 1000
% Dewey/Nortek
% b1=b1
% b2=b3
% b3=b4
% b4=b2
% b5=b5
% phi1=Heading
% phi2=Roll
% phi3=-Pitch


% New May 25th 2017
% For the RDI Sentinel
% Dewey/Sentinel
% b1=b1
% b2=b2
% b3=b3
% b4=b4
% b5=b5
% phi1=Heading
% phi2=+Pitch
% phi3=+Roll


% theta=beam angle, 25º
% Call function like this:
% RS_5beam(b1,b3,b4,b2,b5,25,Heading,-Roll,Pitch,u_east,v_north)


[Nz,Nt]=size(b1);

%% Transform to rads
th=th*pi/180;
phi1=phi1*pi/180;
phi2=phi2*pi/180;
phi3=phi3*pi/180;

%% phi2 and phi3 are assume to vary 
phi1=nanmean(phi1);
phi2=nanmean(phi2);
phi3=nanmean(phi3);

%% Get beam velocities fluctuations
b1=b1-repmat(nanmean(b1,2),[1,Nt]);
b2=b2-repmat(nanmean(b2,2),[1,Nt]);
b3=b3-repmat(nanmean(b3,2),[1,Nt]);
b4=b4-repmat(nanmean(b4,2),[1,Nt]);
b5=b5-repmat(nanmean(b5,2),[1,Nt]);

%% Square and take time average
b1_av=nanmean(b1.^2,2);
b2_av=nanmean(b2.^2,2);
b3_av=nanmean(b3.^2,2);
b4_av=nanmean(b4.^2,2);
b5_av=nanmean(b5.^2,2);


%% u'v'
up=u-repmat(nanmean(u,2),[1,Nt]);
vp=v-repmat(nanmean(v,2),[1,Nt]);
uv_av=nanmean(up.*vp,2);
for i=1:Nz  
    covM=cov(up(i,:),vp(i,:));
    uv_av(i)=covM(1,2);
end


%% u'u'
u2=-1/(4*(sin(th))^6*(cos(th))^2)*(-2*(sin(th))^4*(cos(th))^2*(b2_av+b1_av-2*(cos(th))^2*b5_av)...
    +2*(sin(th))^5*cos(th)*phi3*(b2_av-b1_av));

%% v'v'
v2=-1/(4*(sin(th))^6*(cos(th))^2)*(-2*(sin(th))^4*(cos(th))^2*(b4_av+b3_av-2*(cos(th))^2*b5_av)...
    -2*(sin(th))^4*(cos(th))^2*phi3*(b2_av-b1_av)+2*(sin(th))^3*(cos(th))^3*phi3*(b2_av-b1_av)...
    -2*(sin(th))^5*cos(th)*phi2*(b4_av-b3_av));

%% w2
w2=-1/(4*(sin(th))^6*(cos(th))^2)*(-2*(sin(th))^5*cos(th)*phi3*(b2_av-b1_av)...
    +2*(sin(th))^5*cos(th)*phi2*(b4_av-b3_av)-4*(sin(th))^6*(cos(th))^2*b5_av);

%% u'w'
uw=-1/(4*(sin(th))^6*(cos(th))^2)*((sin(th))^5*cos(th)*(b2_av-b1_av)...
    +2*(sin(th))^4*(cos(th))^2*phi2*(b2_av+b1_av)-4*(sin(th))^4*(cos(th))^2*phi3*b5_av...
    -4*(sin(th))^6*(cos(th))^2*phi2*uv_av);

uw_simp=-1/(2*sin(2*th))*(b2_av-b1_av)+phi3/(sin(th))^2*(1/2*(b2_av+b1_av)-b5_av)-phi2*uv_av;


%% v'w'
vw=-1/(4*(sin(th))^6*(cos(th))^2)*((sin(th))^5*cos(th)*(b4_av-b3_av)...
    -2*(sin(th))^4*(cos(th))^2*phi2*(b4_av+b3_av)+4*(sin(th))^4*(cos(th))^2*phi3*b5_av...
    +4*(sin(th))^6*(cos(th))^2*phi3*uv_av);

%% Convert to Earth Coordinates using heading
u2t=u2*(cos(phi1))^2+2*uv_av*cos(phi1)*sin(phi1)+v2*(sin(phi1))^2;
v2t=v2*(cos(phi1))^2-2*uv_av*cos(phi1)*sin(phi1)+u2*(sin(phi1))^2;
w2t=w2;
uwt=uw*cos(phi1)+vw*sin(phi1);
vwt=vw*cos(phi1)-uw*sin(phi1);
uvt=uv_av*(cos(phi1)^2-sin(phi1)^2)+sin(phi1)*cos(phi1)*(v2-u2);
uw_simp_t=uw_simp*cos(phi1)+vw*sin(phi1);


% TKE and Anisotropy
q2=1/(4*sin(th)^2)*(b1_av+b2_av+b3_av+b4_av-2*(2*cos(th)^2-sin(th)^2)*b5_av...
    -(cot(th)-1)*phi3*(b2_av-b1_av));

F_phi=cot(th)*phi2*(b4_av-b3_av)+(1-2*csc(2*th))*phi3*(b2_av-b1_av);

alpha=(2*sin(th)^2*b5_av+cot(th)*phi3*(b2_av-b1_av)-cot(th)*phi2*(b4_av-b3_av))./...
    b1_av+b2_av+b3_av+b4_av-4*cos(th)^2*b5_av-F_phi;

end

