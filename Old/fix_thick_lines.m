close all

paths_direct=figure(2);
plot(0:domain.dt:domain.endtime,squeeze(output.paths(:,6,:)))
hold on
plot(0:domain.dt:domain.endtime,squeeze(output.paths(:,6,10000:10000:end)),'k-','LineWidth',2)
title('Sampled Paths','FontSize',16)
axis([0,1,-2,2])

saveas(paths_direct,'paths_direct.png')