function [ FluxDensity ] = FluxDensityLoop(p,loopRadius)

%//////////////////////////////////////////////////////////////////////////
%Author: Julien Leclerc
%Institution: University of Houston
%Initial Verion: April 12, 2020
%Current version: April 12, 2020
%Last modified by: Julien Leclerc
%/////////////////////////////////////////////////////////////////////////

%/////////////////////////////////////////////////////////////////////////
%This function computes the magnetic field produced by a single current
%loop in a cylindrical system (r,theta,z). The loop center is at r=0 and 
%z=0. The loop revolution axis is oriented alog the z axis. The calculation
%uses equations presented in: "Simple Analytic Expressions for the Magnetic 
%Field of a Circular Current Loop" by Simpson James C., Lane, John E., 
%Immer Christopher D. and Youngquist Robert C ; NASA 2001.
%/////////////////////////////////////////////////////////////////////////

%This function accepts matrices as input to compute n points at once

%/////////////////////////////////////////////////////////////////////////

%Inputs:

%p: 
%Size= nx3; 
%Unit= meter; 
%Description= Location of the calculation point [R;theta;Z]. %The value of 
%theta has no effect on the output since the system is axisymetic;

%loopRadius:
%Size= 1x1; 
%Unit= meter; 
%Description= Radius of the current loop;

%Output:

%FluxDensity: 
%Size= nx3; 
%Unit= Telsa; 
%Description= Magnetic flux density vector (Br,Btheta,Bz) in the 
%cylindrical coordinate system.
%loopRadius: radius of the current loop, vector 1xn [m]

%////////////////////////////////////////////////////////////////////////

%obtain the number of points to calculate
[m,n]=size(p);

%Define constants:
mu0=4*pi*1e-7; %vacuum permeability 

%//////////////////////////////////////////////////////////
%////// Here starts the calculation of the magnetic field using the
%equations of the paper.
Rc=p(:,1);
Zc=p(:,3);
k=(4.*loopRadius.*Rc)./((loopRadius+Rc).*(loopRadius+Rc)+Zc.*Zc);
[I1,I2] = ellipke(k,10000000);
Brc=((mu0)./(2.*pi.*Rc)).*(Zc./(sqrt((loopRadius+Rc).*(loopRadius+Rc)+Zc.*Zc))).*(-I1+((loopRadius.*loopRadius+Rc.*Rc+Zc.*Zc)./((loopRadius-Rc).*(loopRadius-Rc)+Zc.*Zc)).*I2);
Bzc=((mu0)./(2.*pi)).*(1./(sqrt((loopRadius+Rc).*(loopRadius+Rc)+Zc.*Zc))).*(I1+((loopRadius.*loopRadius-Rc.*Rc-Zc.*Zc)./((loopRadius-Rc).*(loopRadius-Rc)+Zc.*Zc)).*I2);
%End of the calculation of the magnetic field using the equations of the paper.
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%There is cases where the result is NAN
%    at a=r and z=0: This point has no physical sens since the considered
%    wire is of infinitly small diameter
%    at r=0: The equations for Br returns NAN but the result is 0

%The following fix this problem:
for i=1:m
    if Rc(i)==0
        Brc(i)=0;
    end

    if abs(loopRadius-Rc(i))<(loopRadius./100) && abs(Zc(i))<(loopRadius./100)
        Brc(i)=0;
        Bzc(i)=0;
    end
end

%return result
FluxDensity=[Brc,zeros(m,1),Bzc];
end