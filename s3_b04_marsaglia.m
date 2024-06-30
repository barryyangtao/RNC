clear,clc,clf
n = 1e5;    
sigma = 15; % Gaussian FWHM

u1 = rand(floor(n./2),1);
u2 = rand(floor(n./2),1);
z1 = sqrt(-2d0.*log(u1)).*cos(2.*pi.*u2);
z2 = sqrt(-2d0.*log(u1)).*sin(2.*pi.*u2);
z = [z1; z2];
z = sigma./(2d0.*sqrt(2d0.*log(2))).*z;

figure(1)
subplot(2,1,1)
polarplot(atan2(z2,z1),sqrt(z1.^2 + z2.^2),'b.');
title('2D Polar Plot of Rand Numbers');
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
set(gca,'linewidth',1.6);
set(gca,'fontsize',12);

subplot(2,1,2);
histogram(z,100);
title('Histogram of Gaussian Rand Numbers');
xlabel('x');
ylabel('Counts');
legend('z');
box on
set(gca,'linewidth',1.6);
set(gca,'fontsize',12);