clear,clc,clf

x = 0d0:0.001:10d0;
pdf = exp(-x);
cdf = 1-exp(-x);
n = 1e5;
v = -log(rand(n,1));

figure(1)
yyaxis left
plot(x,pdf,'b','LineWidth',2.6)
hold on
plot(x,cdf,'b','LineWidth',2.6)
hold off
yyaxis right
histogram(v)
legend('PDF','CDF','Rand #')
box on
set(gca,'linewidth',2);
set(gca,'fontsize',13.6);