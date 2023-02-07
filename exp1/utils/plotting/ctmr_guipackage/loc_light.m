function loc_light(th, phi)
%function loc_view(theta, phi) 
%this function orients the brain and always puts the lighting behind you
%theta and phi are in degrees, not radians
%make sure the brain plot is your current axis
%this function rotates the brain and lighting to the spherical angle 
%inputted.   it is non-standard b/c of matlab.  so, phi is really "elevation" and not
%phi from standard physics coordinates.  (0,0) is at the back of the brain.  for example: 
%loc_view(180,90) views from the top with the front downward, and
%loc_view(180,90) has the front upward, loc_view(90,0) is from the right,
%maybe i should alter it later so there is the option
%of inputting cartesian, but that can be ambiguous, for example [0 0 1] has ambiguous orientation.
%kjm 08/06


% view(th,phi),
%%if you want to go the other way
% [th,phi,r]=cart2sph(view_pt(1),view_pt(1),view_pt(1)); %in radians, but "view" uses degrees with different origin
% th=360*(th+pi/2)/(2*pi);
% phi=360*phi/(pi)

view_pt=[cosd(th-90)*cosd(phi) sind(th-90)*cosd(phi) sind(phi)];
%in order to change the direction of the light:
a=get(gca,'Children');
for i=1:length(a)
    b=get(a(i));
    if strcmp(b.Type,'light') %find the correct child (the one that is the light)
        %object for light is the 2nd element, then use a 
        set(a(i),'Position',view_pt) 
        %or something
    end
end