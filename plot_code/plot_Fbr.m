% Bingchen Liu June 14, 2024
% This code explore the cross-shore structure of Fbr mag and its x and y
% component
load /data1/bliu/data/SS_raw.mat % now S(1--120) struct is in memory

for i = 1:120
    figure(1)
    plot(S(i).X,S(i).Fbr_mag,'Linewidth',3)
    hold on 
    plot(S(i).X,S(i).Fbr_x,'Linewidth',3)
    hold on 
    plot(S(i).X,S(i).Fbr_y,'Linewidth',3)
    legend('Norm','x-component','y-component','Location','northwest')
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$Fbr (s^{-2})$','Interpreter','latex')
    niceplot_nobold_nomintick(20)
    grid on
    hold off 
    pause(0.5)
end 