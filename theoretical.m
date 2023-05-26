clear
Tc = 2/log(1+sqrt(2));
N = 200; K = 1000; J = 1;
T = (0.5:0.1:2);

figure;
hold on;
Tf = (0.5:0.00000001:3);
plot(Tf,real((1 - sinh(2*J*Tf.^-1).^-4).^(1/8)));
set(gca,'FontSize',20);
xlabel('T/Tc');
ylabel('Magnetization');
