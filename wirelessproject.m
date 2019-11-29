%random coordinate generation

    bs = [0,0,0]; %base station
    t1 = -50 + rand(3,1)*100;
    t2 = -50 + rand(3,1)*100;
    i1 = -50 + rand(3,1)*100;
    i2 = -50 + rand(3,1)*100;
    while(1)
        t1 = t1-5+10*rand(3,1);
        t2 = t2-5+10*rand(3,1);
        %i1 = i1-5+10*rand(3,1);
        %i2 = i2-5+10*rand(3,1);
        v=[t1,t2,i1,i2];
        plot3(bs(1),bs(2),bs(3),'vr','MarkerSize',9,'LineWidth',1.5,'Color',[0.6350 0.0780 0.1840]);
        hold on;
        plot3(v(1,(1:2)),v(2,(1:2)),v(3,(1:2)),'ob','LineWidth',1.5); %tracked terminals
        plot3(v(1,(3:4)),v(2,(3:4)),v(3,(3:4)),'x', 'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); %interferers
        xlim([-100,100])
        ylim([-100,100])
        zlim([-100,100])
        grid on;
        legend('Base Station','Terminals','Interferers');
        hold off;
        pause(10)
    end