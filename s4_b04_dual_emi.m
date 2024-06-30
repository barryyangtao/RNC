clear,clc,clf
%% Simple retrospective MZI Fourier spectroscopy simulation
%% Barry Y. Li and Tim Duong (2024)
%% Parameter Box ##########################################################
delta = -5e-5:2.5e-8:5e-5;                         % stage position (meter)
emitter_r = [1.25 0.40 0.16 -6];            % red {E0,Gamma_g,Gamma_l,beta}
emitter_b = [1.35 0.18 0.08 -3];           % blue {E0,Gamma_g,Gamma_l,beta}
fi = 1.0d0;                                    % energy (E) RNG start point       
ff = 1.5d0;                                      % energy (E) RNG end point
%% ########################################################################


%% ############ you do not have to change anything below :-) ##############
tic

avg_time = 1e5;                                    % average life time (ps)
pulse_space = avg_time.*exp(1);                % pulse in-between time (ps)
photon_per_pulse = 0.05;                      % average # photons per pulse
n_photon_goal = 1e4;                           % desired # photon to detect
n_pulse = 1./photon_per_pulse.*n_photon_goal;

for i = 1:n_pulse
    pulse_t(i) = (i-1).*pulse_space;
end
measure_time_in_s = pulse_t(end)./1e12;

for i = 1:length(delta)
    intere(i,:,:) = spec(delta(i),avg_time,pulse_space,photon_per_pulse,...
        n_photon_goal,emitter_r,emitter_b,fi,ff,n_pulse,pulse_t);
end

figure(1)
for i = 1:length(delta)
    cum_count(i) = sum(intere(i,2,:));
end

[earray,P2] = fft_spec(delta,cum_count);


figure(1)
plot(delta'.*1e6,cum_count'./max(max(cum_count')),'Color',...
    [0 0.4470 0.7410],'LineWidth',1.6)
xlim([-50 50])
ylim([-1.1 1.1])
xlabel('Stage Position (um)')
ylabel('Signal S_A-S_B')
title('Dual-emitter Interferogram')
box on
set(gca,'fontsize',16);
set(gca,'linewidth',1.6);

figure(2)
plot(earray,P2/max(P2),'k:','LineWidth',2.6)
xlim([fi ff])
ylim([0 1.1])
xlabel('Energy (eV)')
ylabel('Norm. Intensity')
title('Dual-emitter FFT Spectrum')
box on
set(gca,'fontsize',16);
set(gca,'linewidth',1.6);
toc


%% Core simulation functions ##############################################
%% ------------------------------------------------------------------------
%% 1. Photon (pulsed) stream simulation - Sub-poisson (antibunching)
function ps = spec(delta,avg_time,pulse_space,photon_per_pulse,...
        n_photon_goal,emitter_r,emitter_b,fi,ff,n_pulse,pulse_t)
c = 299792458;
hbar = 1.054571817e-34;             

detc = 1d0./2d0.*erfc(1000d0.*(rand(n_pulse,1)-photon_per_pulse));
pulse_with_photon = pulse_t(find(detc == 1));
tau(:) = -avg_time.*log(rand(length(pulse_with_photon),1));

T = pulse_with_photon + tau;
T(end) = [];

np = length(T);

%% ------------------------------------------------------------------------
%% 2. Now call out all the frequency-domain info
fd = (ff - fi)./(np - 1d0);
x = 1.0:fd:1.5;   
fr = simu_spec(emitter_r,x);
fb = simu_spec(emitter_b,x);
cdf_r = cdf_gen(x,fr);
cdf_b = cdf_gen(x,fb);
randyr = rand(np,1);
randyb = rand(np,1);
for i = 1:np
    [vr(i) cr(i)] = closest_value(cdf_r, randyr(i));
    [vb(i) cb(i)] = closest_value(cdf_b, randyb(i));
end
new_vecr = x(cr);
new_vecb = x(cb);
new_vec = [new_vecr new_vecb]';

%% ------------------------------------------------------------------------
%% 3. The result will be (t,omega)
omega = new_vec.*1.602176565e-19./hbar;
prob_a = 1d0./2d0.*(1+cos(delta.*omega./c)); 
prob_b = 1d0./2d0.*(1-cos(delta.*omega./c));

for i = 1:length(prob_a)
    detc_choice(i) = detector(prob_a(i));
end
count_array = round(detc_choice(:));

deteca = sum(count_array);
detecb = length(prob_a) - deteca;


[counts_a,centers_a] = histcounts(tau,200);
counts_a(numel(centers_a)) = 0d0;
ps = [centers_a; (deteca-detecb).*counts_a];
end


%% dependency functions ###################################################
%% 1. ---------------------------------------------------------------------
function build_cdf = cdf_gen(x,f)
% this is before normalization: 
for i = 1:(length(x)-1d0)
    intf(i) = 1d0./2d0.*(f(i+1)+f(i)).*(x(i+1)-x(i));         % trapezoidal 
end
cdf_f = 0d0;
for i = 1:(length(x)-1d0)
    cdf_f(i+1) = cdf_f(i)+intf(i);                % get the CDF numerically
end

% this is after normalization:
f = f./cdf_f(end);
for i = 1:(length(x)-1d0)
    intf(i) = 1d0./2d0.*(f(i+1)+f(i)).*(x(i+1)-x(i));% trapezoidal integral
end
cdf_f = 0d0;
for i = 1:(length(x)-1d0)
    cdf_f(i+1) = cdf_f(i)+intf(i);                % get the CDF numerically
end
build_cdf = cdf_f;
end

%% 2. ---------------------------------------------------------------------
function spec = simu_spec(data,xa)
    mu  = [data(1)];
    gba = [data(2)];
    lba = [data(3)];
    gsk = [data(4)];
    
    f1 = g(xa,mu(1),gba(1)./(2.*sqrt(2.*log(2))),gsk(1),mu(1),lba(1)./2d0);

    f1n = f1./max(f1);
    spec = f1n;
end

%% 3. ---------------------------------------------------------------------
function spdf = g(x,mean_g,std_g,skewness,mean_l,w_l)
    % gaussian normalization constant
    gauss_norm = 1d0./(std_g.*sqrt(2d0.*pi));

    % skewed Gaussian
    gauss_pdf = gauss_norm.*exp(-1d0./2d0.*((x-mean_g)./std_g).^2d0);
    gauss_pdf = gauss_pdf.*(1+erf(skewness.*(x-mean_g)./(std_g.*sqrt(2))));

    % lorentzian PDF
    lorentz_pdf = 2d0./pi.*(w_l./(2d0.*((x-mean_l).^2d0+(w_l./2).^2)));

    % compute the product of the two distributions
    spdf = gauss_pdf.*lorentz_pdf;
end

%% 4. ---------------------------------------------------------------------
function [v, inf] = closest_value(arr, val)
len = length(arr);
inf = 1;
sup = len;
while sup - inf > 1
    med = floor((sup + inf)/2);
    if arr(med) >= val 
        sup = med;
    else
        inf = med;
    end
end
if sup - inf == 1 && abs(arr(sup) - val) < abs(arr(inf) - val)
    inf = sup;
end  
v = arr(inf);
end

%% 5. ---------------------------------------------------------------------
function detc = detector(a)
    x = rand(1);
    detc = 1d0./2d0.*erfc(1000d0.*(x-a));
end

%% 6. ---------------------------------------------------------------------
function [earray,P2] = fft_spec(delta,cum_count)
L = length(cum_count);            
X = cum_count;
Y = fft(X);
P2 = abs(Y/L);

ddelta = mean(diff(delta))*1e2;                                     % in cm
ndelta = length(delta);

k_max = 1d0./(ddelta);                                            % in cm-1
emax = k_max.*1.2398e-4;                                        % now in eV
k_min = 1d0./(ndelta.*ddelta);                                    % in cm-1
emin = k_min.*1.2398e-4;                                        % now in eV

earray = emin:((emax-emin)./(ndelta-1)):emax;
earray = earray - emin;                                    % spectral shift
end