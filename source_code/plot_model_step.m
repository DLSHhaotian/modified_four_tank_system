function plot_model_step(t,y)
figure
subplot(4,1,1)
plot(t,y(1,:),'LineWidth',2.0)
title('y1')
xlabel('t(s)');
ylabel('h[cm]');
subplot(4,1,2)
plot(t,y(2,:),'LineWidth',2.0)
title('y2')
xlabel('t(s)');
ylabel('h[cm]');
subplot(4,1,3)
plot(t,y(3,:),'LineWidth',2.0)
title('y3')
xlabel('t(s)');
ylabel('h[cm]');
subplot(4,1,4)
plot(t,y(4,:),'LineWidth',2.0)
title('y4')
xlabel('t(s)');
ylabel('h[cm]');
%sgtitle("The step responses of 10% step for deterministic model")
end