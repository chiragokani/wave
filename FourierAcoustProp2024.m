%%  Angular Spectrum Computation of Acoustic Field Radiated by a Velocity Source
%%
%%  ME/EE 384N-8 Wave Phenomena, Spring 2024
%%  Mark F. Hamilton
%%  University of Texas at Austin
%%

clear all

% Dimensionless variables: X=x/a, Y=y/a, Z=z/z0, z0=0.5*k*a^2, K=ka

ka = 50;
Xmax = 20;  % Ymax = Xmax; square domain assumed (-Xmax < X < Xmax)
N = 512;  % must be power of 2, NxN = number of FFT points
KxyMax = 10; % angular spectrum plotted for |k_x*a|,|k_y*a|<KxyMax, KxyMax<(N/2)*(pi/Xmax)
Zxy = 0.3;  % distance where |p(x,y)| is plotted
XoutXY = 2;  % max value of X displayed in plot in XY plane at Z=Zxy
Zmax = 1.5; % end point of propagation curve
Zpoints = 100; % number of output points in Z direction
XoutXZ = 4;  % max value of X displayed in plot of XZ plane (same for YZ plane)

x = (2*Xmax/N)*(-N/2:N/2-1); 
y = (2*Xmax/N)*(-N/2:N/2-1); 
[X,Y] = meshgrid(x,y);
R=sqrt(X.^2+Y.^2); 
KxyPlot = (pi/Xmax)*(-N/2:N/2-1);

% Velocity source functions
%
Source = 1.*(R<=1.0);  % uniform circular piston
% Source = exp(-R.^30);  % circular super-Gaussian
% Source = (abs(X)<=1.5).*(abs(Y)<=0.5);  % uniform rectangular piston
% Source = exp(-(X/1.5).^30-(Y/0.5).^30);  % rectangular super-Gaussian
% Source = sin(pi*X).*(abs(X)<=1).*sin(pi*Y).*(abs(Y)<=1);  % membrane mode
% Source = exp(0.5*X).*sin(pi*X).*(abs(X)<=1).*(abs(Y)<=0.5);  % shaded dipole
% Source = cos(3*pi*X/2).*(abs(X)<=1);  % 1D membrane mode
% Source = cos(2*pi*X);  % infinite 1D standing wave

Spectrum = fft2(Source);  % angular spectrum
Smax = max(max(abs(Spectrum)));  % maximum magnitude
Sdisplay = abs(fftshift(Spectrum))/Smax;  % displayed spectrum
[KxIndex,KyIndex] = meshgrid(-N/2:N/2-1,-N/2:N/2-1);

% Match wavenumbers to fft output

DeltaKxy = pi/Xmax;
for n=1:N/2+1
    Kxy(n)=(n-1)*DeltaKxy;
end
for n=N/2+2:N
    Kxy(n)=-(N+1-n)*DeltaKxy;
end

K = ka;
Kx = repmat(Kxy,N,1);
Ky = repmat(Kxy',1,N);
Kz = sqrt(K^2-Kx.^2-Ky.^2);
Z = Zxy;

P = ifft2(Spectrum.*(K./Kz).*exp(i*K*Kz*Z/2));

% Calculate particle velocity and intensity components

% Ux = ifft2(Spectrum.*(Kx./Kz).*exp(i*K*Kz*Z/2));
% Uy = ifft2(Spectrum.*(Ky./Kz).*exp(i*K*Kz*Z/2));
% Uz = ifft2(Spectrum.*(Kz./Kz).*exp(i*K*Kz*Z/2));

% Ix = 0.5*real(P.*conj(Ux));
% Iy = 0.5*real(P.*conj(Uy));
% Iz = 0.5*real(P.*conj(Uz));

% Plot source function

figure(1)  
xMaxIndex = floor(XoutXY/Xmax*N/2); % compute the index of this value
xPlotIndex = (N/2-xMaxIndex):(N/2+xMaxIndex+2); % indices of the plotted values
surfl(X(xPlotIndex,xPlotIndex),Y(xPlotIndex,xPlotIndex),real(Source(xPlotIndex,xPlotIndex)))
xlim([-2 2]), ylim([-2 2]),
shading interp, colormap copper,
xlabel('x/a'); ylabel('y/a'); zlabel('u_0(x,y)');
title(['Velocity Source:  ' ' ka = ' num2str(ka) ...
    ',  x_{max}/a = y_{max}/a = ' num2str(Xmax) ',  N_{pts} = ' num2str(N) ...
    ' x ' num2str(N)])

% Plot angular spectrum of source velocity

figure(2)  
kMaxIndex = floor(KxyMax/(pi/Xmax))-1; % compute the index of this value
kPlotIndex = (N/2-kMaxIndex):(N/2+kMaxIndex+2); % indices of the plotted values
surfl(KxyPlot(kPlotIndex),KxyPlot(kPlotIndex),abs(Sdisplay(kPlotIndex,kPlotIndex)));
xlim([-KxyMax KxyMax]), ylim([-KxyMax KxyMax]),
shading interp, colormap copper,
xlabel('k_xa'); ylabel('k_ya'); 
zlabel('|U_0(k_x,k_y)|');
title(['Angular Spectrum:  ' ' (k_xa)_{max} = (k_ya)_{max} = ' num2str((N/2)*(pi/Xmax))])

% Plot pressure in XY plane at Z=Zxy

figure(3)  
surfl(X(xPlotIndex,xPlotIndex),Y(xPlotIndex,xPlotIndex),abs(P(xPlotIndex,xPlotIndex)))
xlim([-XoutXY XoutXY]), ylim([-XoutXY XoutXY]),
shading interp, colormap copper,
xlabel('x/a'); ylabel('y/a'); zlabel('|p(x,y)|/\rho_0c_0u_0');
title(['Pressure Field in Plane' ' z/z_0 = ' num2str(Z)...
    '  (z/\lambda = ' num2str(ka^2*Z/(4*pi)) ')'])

% Compute whole pressure field out to Zmax

m1=1+round((N/2)*(1-XoutXZ/Xmax)); % lower index corresponding to -XoutXZ
m2=N-m1+2;  % upper index corresponding to XoutXZ (=YoutYZ)

figure(4)
for n=1:Zpoints
    Z=(n/Zpoints)*Zmax;
    Zaxial(n)=Z;
    P=ifft2(Spectrum.*(K./Kz).*exp(i*K*Kz*Z/2));
    Paxial(n)=abs(P(N/2+1,N/2+1));
    Pxz(:,n)=P(N/2+1,:); % preparing field for plot in xz plane
    Pxzout(1:(m2-m1+1),n)=Pxz(m1:m2,n);
    Pyz(:,n)=P(:,N/2+1); % preparing field for plot in yz plane
    Pyzout(1:(m2-m1+1),n)=Pyz(m1:m2,n);
end
plot(Zaxial,Paxial);
xlabel('z/z_0'),ylabel('|p(0,0)|/\rho_0c_0u_0');
title(['Axial Propagation Curve'])

% Plot pressure in XZ plane

figure(5)
for m=1:(m2-m1+1)
    Xout(m)=-XoutXZ+((m-1)/(m2-m1))*2*XoutXZ; % (=Yout)
end
imagesc(Zaxial,Xout,abs(Pxzout));
xlabel('z/z_0'),ylabel('x/a')
title(['Pressure Field in x-z Plane'])

% Plot pressure in YZ plane

figure(6)
imagesc(Zaxial,Xout,abs(Pyzout));
xlabel('z/z_0'),ylabel('y/a')
title(['Pressure Field in y-z Plane'])
