function [] = plotting(T, H, Qout, name)
% Plotting 
tl = strcat('Water Level in 4-Tank-System - ', name);
figure('Name', tl)
title(tl)
subplot(2,2,1)
plot(T,H(:, 1))
title('Tank 1','interpreter','latex')
xlabel('Time [sec]','interpreter','latex')
ylabel('Water Level [cm]','interpreter','latex')

subplot(2,2,2)
plot(T,H(:, 2))
title('Tank 2','interpreter','latex')
xlabel('Time [sec]','interpreter','latex')
ylabel('Water Level [cm]','interpreter','latex')

if (size(H,2) == 4)
    subplot(2,2,3)
    plot(T,H(:, 3))
    title('Tank 3','interpreter','latex')
    xlabel('Time [sec]','interpreter','latex')
    ylabel('Water Level [cm]','interpreter','latex')
    
    subplot(2,2,4)
    plot(T,H(:, 4))
    title('Tank 4','interpreter','latex')
    xlabel('Time [sec]','interpreter','latex')
    ylabel('Water Level [cm]','interpreter','latex')
end

tl2 = strcat('Outflow Rate in 4-Tank-System - ', name);
figure('Name', tl2)
title(tl2)
subplot(2,2,1)
plot(T,Qout(:, 1))
title('Tank 1','interpreter','latex')
xlabel('Time [sec]','interpreter','latex')
ylabel('Flow rate [cm/s]','interpreter','latex')

subplot(2,2,2)
plot(T,Qout(:, 2))
title('Tank 2','interpreter','latex')
xlabel('Time [sec]','interpreter','latex')
ylabel('Flow rate [cm/s]','interpreter','latex')

subplot(2,2,3)
plot(T,Qout(:, 3))
title('Tank 3','interpreter','latex')
xlabel('Time [sec]','interpreter','latex')
ylabel('Flow rate [cm/s]','interpreter','latex')

subplot(2,2,4)
plot(T,Qout(:, 4))
title('Tank 4','interpreter','latex')
xlabel('Time [sec]','interpreter','latex')
ylabel('Flow rate [cm/s]','interpreter','latex')
end

