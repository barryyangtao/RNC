clear, clc, clf
%% simple numerical random number generator from inv(CDF)
%% Barry Y. Li and Tim Duong (2024)
%% in the "Workspace", the output new_vec is the array of rand #'s you want

lower = 0;                                     % sampling range lower bound
upper = 16;                                    % sampling range upper bound
%upper = pi.^2./4;                             
n = 1e7;                                       % number of draw random #'s
randy = rand(n,1);
x = lower:(1e-4):upper;

% please input f as a known pdf -------------------------------------------          
%f = exp(-x);
%f = 1d0./sqrt(pi).*exp(-(x-3).^2d0);
%f = 2d0./pi.*sin(2d0.*sqrt(x));
f = x.^2.*exp(-x./2).*(sin(x)).^2;


%% ############ you do not have to change anything below :-) ##############

% this is before normalization --------------------------------------------
for i = 1:(length(x)-1d0)
    intf(i) = 1d0./2d0.*(f(i+1)+f(i)).*(x(i+1)-x(i));% trapezoidal integral
end
cdf_f = 0d0;
for i = 1:(length(x)-1d0)
    cdf_f(i+1) = cdf_f(i)+intf(i);                % get the CDF numerically
end
% -------------------------------------------------------------------------

% this is after normalization ---------------------------------------------
f = f./cdf_f(end);
for i = 1:(length(x)-1d0)
    intf(i) = 1d0./2d0.*(f(i+1)+f(i)).*(x(i+1)-x(i));% trapezoidal integral
end
cdf_f = 0d0;
for i = 1:(length(x)-1d0)
    cdf_f(i+1) = cdf_f(i)+intf(i);                % get the CDF numerically
end
% -------------------------------------------------------------------------

figure(1)
yyaxis left
plot(x,f,'LineWidth',2)
ylabel('PDF(x)')
yyaxis right
plot(x,cdf_f,'LineWidth',2)
xlabel('x')
ylabel('CDF(x)')
xlim([lower upper])
ylim([0 1.1])
box on
set(gca,'linewidth',2);
set(gca,'fontsize',16);

for i = 1:n
    c(i) = closest_value(cdf_f, randy(i));
end
new_vec = x(c);

figure(2)
histogram(new_vec,100,'LineWidth',1.36,'EdgeColor','b','FaceAlpha',0)
xlim([lower upper])
xlabel('x')
ylabel('Counts')
box on
set(gca,'linewidth',2);
set(gca,'fontsize',16);


%% function 1: trapezoidal integral ---------------------------------------
function nig = numint(x,f,i)
    nig = trapz([x(1):x(i)],[f(1):f(i)]);
end

%% function 2: binary search subroutine -----------------------------------
function v = closest_value(y,x)
findind = 1d0;
endind = length(y);

% binary search for index
while ((endind - findind) > 1d0)
    midind = floor((endind+findind)./2d0);
    if (y(midind) >= x) 
        endind = midind;
    else
        findind = midind;
    end
end
if ((endind-findind) == 1d0) && (abs(y(endind)-x) < abs(y(findind)-x))
    findind = endind;
end  
v = findind;

end