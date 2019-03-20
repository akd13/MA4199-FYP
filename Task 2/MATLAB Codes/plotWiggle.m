%% Percentage of covered areas on known genes

% countOnKnown=sum(countPos(coveredAreaTrue(1:length(countPos))));
% count = sum(countPos);
% percentage = countOnKnown /count



%% plot

RPM=countPos/lib_size*1000000;

countPos_exist = RPM;
countPos_exist(~coveredAreaTrue(1:length(RPM),:))=0;
countPos_nonexist =RPM;
countPos_nonexist(coveredAreaTrue(1:length(RPM),:))=0;

if doPlot
    figure;hold on;
    plotx=[0:length(RPM)-1,0:length(RPM)-1];
    plot(plotx,[countPos_exist(:,1);countPos_exist(:,2)]);
    plot(plotx,[countPos_nonexist(:,1);countPos_nonexist(:,2)]);
    legend('maps to known region','unknown region');
    title(['chr',int2str(chr),'wiggle plot at T',int2str(T-1)]);
    ylabel('RPM');
    xlabel('Coordinates');
end
%saveas(gcf,['wiggle/chr',int2str(chr),'.png']);
%ylim([0,5]),
% saveas(gcf,['wiggle/zoomedIn chr',int2str(chr),'.png']);
% close all;

%end

cover_percentage = sum(countPos_exist) ./ ...
    (sum(countPos_exist)+sum(countPos_nonexist))