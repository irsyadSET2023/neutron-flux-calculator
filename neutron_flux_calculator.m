na=6.023e23; er=3.2e-11; 



fprintf('Press 1 for Rectangular Parallel-Piped\n')
fprintf('Press 2 for Sphere\n')
fprintf('Press 3 for Cylinder\n')
s=input('Choose your reactor shape:');
clc



if ((s==1))
x=input('Enter the length of x-direction (cm):');
y=input('Enter the length of y-direction (cm):');
zd=input('Enter the length of z-direction (cm):');
elseif((s==2))
radius=input('Enter the length of radius (cm):');
elseif((s==3))
radius=input('Enter the length of radius (cm):');
height=input('Enter the length of height (cm):');
else
clc
fprintf('\n')
fprintf('Error please run the program again')
return
end

fprintf('Press 1 for Uranium-235\n')
fprintf('Press 2 for Uranium-233\n')
fprintf('Press 3 for Plutonium-239\n')
f=input('Choose your fuel (fissile material):');
clc

if ((f==1))
micfiss=582e-24; micabsfuel=681e-24; massfuel=235; fissenergy=200e6;
elseif ((f==2))
micfiss=529e-24; micabsfuel=575e-24; massfuel=233; fissnergy=197.9e6;
elseif ((f==3))
micfiss=749e-24; micabsfuel=1020e-24; massfuel=239; fissenergy=207.1e6;
else
clc
fprintf('\n')
fprintf('Error please run the program again')
return
end
fprintf('Press 1 for Light Water\n')
fprintf('Press 2 for Graphite\n')
fprintf('Press 3 for Heavy Water\n')
fprintf('Press 4 for Berylium\n')

m=input('Choose the moderator:');
clc

if ((m==1))
difflengthmod=8.1; reproductionfact=2.065; neutage=27; 
micabsmod=0.664e-24; massmod=18.0153; density=1;

elseif ((m==2))
difflengthmod=3500; reproductionfact=2.065; neutage=368; 
micabsmod=0.0034e-24; massmod=12; density=1.6;

elseif(m==3)
difflengthmod=8.1; reproductionfact=2.065; neutage=131; 
micabsmod=0.00133e-24; massmod=20.0276; density=1.105;

elseif(m==4)
difflengthmod=480; reproductionfact=2.065; neutage=102; 
micabsmod=0.0092e-24; massmod=9.0122; density=1.85;

else
clc
fprintf('\n')
fprintf('Error please run the program again')
return
end




if ((s==1))
volume=x*y*zd;
buckling=((pi/x)^2)+((pi/y)^2)+((pi/zd)^2);

elseif((s==2))
volume=4*pi*(radius^3);
buckling=((pi/radius)^2);

elseif((s==3))
volume=pi*radius^2*height;
buckling=((2.405/radius)^2)+(pi/height)^2;
end
z=(1+(buckling*(difflengthmod+neutage)))/(reproductionfact-1-(buckling*neutage));

criticalmassmod=(volume*density)/1000;
criticalmassfuel=z*((micabsmod*massfuel)/(micabsfuel*massmod))*criticalmassmod;
thermalutil=z/(1+z);
kinf=reproductionfact*thermalutil;
difflengthfuel=(1-thermalutil)*difflengthmod;
macfiss=na*(criticalmassfuel*1000)/(volume*massfuel)*micfiss;

powerproduct=(criticalmassfuel/massfuel)*na*(fissenergy*1.6e-19)*(1/3600*24);

clc

fprintf('Critical mass is (kg):%3.2e\n',criticalmassfuel)
fprintf('Maximum Power Production is (Watt-thermal-day):%9.2e\n',powerproduct)

p=input('Insert your operating power (watt-thermal):');

if (s==1)
A=(3.87*p)/(volume*er*macfiss);

xx=[-x/2:1:x/2];
yy=[-y/2:1:y/2];
zz=[-zd/2:1:zd/2];

[x1,y1]=meshgrid(xx,yy);
flux=A.*cos(pi.*x1/x).*cos(pi.*y1/y);
flux1=A.*cos(pi.*xx/x);
flux2=A.*cos(pi.*yy/y);
flux3=A.*cos(pi.*zz/zd);



figure(1)
subplot(3,1,1)
plot(xx,flux1)
xlabel('X-Distance (cm)')
ylabel('Flux')
grid
subplot(3,1,2)
plot(yy,flux2)
xlabel('Y-Distance (cm)')
grid
ylabel('Flux')
subplot(3,1,3)
plot(zz,flux3)
xlabel('Z-Distance (cm)')
ylabel('Flux')
grid
figure(2)
surf(x1,y1,flux)

resultx=[xx;flux1];
resulty=[yy;flux2];
resultz=[zz;flux3];

fprintf('      Neutron Flux in X-Direction\n        ')
fprintf('%9s   %9s \n','X-Direction','Neutron flux')
fprintf('%14.2f     %14.4e \n',resultx)


fprintf('      Neutron Flux in Y-Direction\n        ')
fprintf('%9s   %9s \n','X-Direction','Neutron flux')
fprintf('%14.2f     %14.4e \n',resulty)



fprintf('      Neutron Flux in Z-Direction\n        ')
fprintf('%9s   %9s \n','X-Direction','Neutron flux')
fprintf('%14.2f     %14.4e \n',resultz)


elseif(s==2)
A=(p)/((4*radius^2)*er*macfiss);
r=[-radius:1:radius];
rr=meshgrid(r);
flux=A.*(1./rr).*sin((pi.*rr)/radius);
flux2=A.*(1./r).*sin((pi.*r)/radius);

plot(r,flux2)
xlabel('Radius (cm)')
ylabel('Neutron Flux')
title('Neutron Flux Distribution in Spherical Reactor')

results=[r;flux2];

fprintf('      Neutron Flux Distribution in Sphere Reactor\n        ')
fprintf('%9s   %9s \n','Radius ','Neutron flux')
fprintf('%14.2f     %14.4e \n',results)



elseif (s==3)
A=(3.63*p)/(volume*er*macfiss);

r=[-radius:1:radius];
h=[-height/2:1:height/2];

[rr,hh]=meshgrid(r,h);
flux=A.*(real(besselj(0,((2.405.*rr)/radius),0))).*cos((pi.*hh)/height);

flux2=A.*(real(besselj(0,((2.405.*r)/radius),0)));
flux3=A.*cos((pi.*h)/height);




figure (1)
subplot(2,1,1)
plot(r,flux2)
xlabel('Radius(cm)')
ylabel('Neutron Flux')
title('Radial Distrubition')
grid

subplot(2,1,2)
plot(flux3,h)
ylabel('Height(cm)')
xlabel('Neutron Flux')
title('Axial Distrubition')
grid

figure(2)
surf(rr,hh,flux)
xlabel('Radius(cm)')
ylabel('Height (cm)')
zlabel('Neutron Flux')
title('3D-plot for Cylinder Reactor')
grid

resultr=[r;flux2];
resulta=[h;flux3];

fprintf('      Neutron Flux in Radial-Direction\n        ')
fprintf('%9s   %9s \n','X-Direction','Neutron flux')
fprintf('%14.2f     %14.4e \n',resultr)

fprintf('      Neutron Flux in Axial-Direction\n        ')
fprintf('%9s   %9s \n','X-Direction','Neutron flux')
fprintf('%14.2f     %14.4e \n',resulta)



else
fprintf('Error please run the program again')
return
end

fissionrate=(micabsfuel/micfiss)*1.05e-9*p;
fuelrate=criticalmassfuel/fissionrate;

fprintf('Total mass of fuel consumed for a day is (gram):%3.2f\n',fissionrate*1000)
fprintf('Life-span of the fuel is (day):%3.2f\n',fuelrate)
fprintf('Life-span of the fuel is (year):%3.2f\n',fuelrate/365)