%% parameters

barcodes = ['ATCACG-s_7_1';'TTAGGC-s_7_1';'CAGATC-s_7_1';'GCCAAT-s_7_1';'CTTGTA-s_7_1'];
chr=22;
chr_str=['chr',int2str(22)];
total_length22 = 51304566;
doPlot=0;

window_size=1000;
count_window = zeros(floor(total_length22/window_size),2);
count_windows = zeros(floor(total_length22/window_size),2,3);
RPMs = zeros(total_length22,2,3);
non_exists = zeros(total_length22,2,3);
exists=zeros(total_length22,2,3);


%% collate
for T = 1:5
    % filename
    barcode = barcodes(T,:);
    bl_prefix = [barcode,'_genome.'];
    bl_file= [barcode,'/',bl_prefix];
    if 1
        readBL;
        CountingBL;
        refFlatInfo;
        plotWiggle;
    end
       
    for i=1:window_size:length(countPos_nonexist)-window_size           % 1,1001,2001 and so on
        windowIdx = ceil(i/1000);                                       % 1,2,3,4
        count_window(windowIdx,:) = sum(countPos_nonexist(i:i+999,:));  % 1-1000,1001-2000
    end
    count_windows(:,:,T) = count_window;
    RPMs(:,:,T) = RPM;          %store RPM based on timepoints 
    non_exists(:,:,T) = countPos_nonexist;
    exists(:,:,T) = countPos_exist;
end

%% plot
if doPlot
figure;hold on;
plot(count_windows(:,:,1));
plot(count_windows(:,:,2));
plot(count_windows(:,:,3));
plot(count_windows(:,:,4));
plot(count_windows(:,:,5));

end

%% compare

% i.e. -200 / -50 = 4 >= 2
idx = ((count_windows(:,:,1) ./ count_windows(:,:,2)) >= 2 ...      % double
    | (count_windows(:,:,1) ./ count_windows(:,:,2)) <= 0.5) ...    % half
    & ((count_windows(:,:,2) ./ count_windows(:,:,3)) >= 2 ...      
    | (count_windows(:,:,2) ./ count_windows(:,:,3)) <= 0.5) ...
    & ((count_windows(:,:,3) ./ count_windows(:,:,4)) >= 2 ...      
    | (count_windows(:,:,3) ./ count_windows(:,:,4)) <= 0.5) ...
    & ((count_windows(:,:,4) ./ count_windows(:,:,5)) >= 2 ...      
    | (count_windows(:,:,4) ./ count_windows(:,:,5)) <= 0.5) ...
    & (abs((count_windows(:,:,1)) >= 50 | ...
    abs(count_windows(:,:,2)) >= 50) | abs(count_windows(:,:,3))>=50);   % filter small ones

regions=zeros(total_length22,2);
[r,c]= find(idx);
for i=1:length(r)
    regions(r(i)*1000-999:r(i)*1000,c(i)) = 1;
end
[r,c]=find(regions);

if 0
figure;hold on;
plot([r(c==1);r(c==2)],[RPMs(r(c==1),1,1);RPMs(r(c==2),2,1)]);
plot([r(c==1);r(c==2)],[RPMs(r(c==1),1,2);RPMs(r(c==2),2,2)]);
plot([r(c==1);r(c==2)],[RPMs(r(c==1),1,3);RPMs(r(c==2),2,3)]);
plot([r(c==1);r(c==2)],[RPMs(r(c==1),1,4);RPMs(r(c==2),2,4)]);
plot([r(c==1);r(c==2)],[RPMs(r(c==1),1,5);RPMs(r(c==2),2,5)]);

legend('T0','T1','T2','T3','T4');
title('unknown regions with fold change');
ylabel('RPM');
end

if 0
figure;hold on;
plot(plotx,[non_exists(:,1,1);non_exists(:,2,1)]);
plot(plotx,[non_exists(:,1,2);non_exists(:,2,2)]);
plot(plotx,[non_exists(:,1,3);non_exists(:,2,3)]);
plot(plotx,[non_exists(:,1,4);non_exists(:,2,4)]);
plot(plotx,[non_exists(:,1,5);non_exists(:,2,5)]);

legend('T0','T1','T2','T3','T4');
title('unknown regions');
ylabel('RPM');xlabel('Coordinates');
end

if 1
    [rs,cs]= find(idx);
    for i = 1:length(rs)
        xs = rs(i)*1000-999:rs(i)*1000;
        c = cs(i);
        figure;hold on;
        if 0
            plot(xs,non_exists(xs,c,1));
            plot(xs,non_exists(xs,c,2));
            plot(xs,non_exists(xs,c,3));
            plot(xs,non_exists(xs,c,4));
            plot(xs,non_exists(xs,c,5));

            plot(xs,exists(xs,c,1));
            plot(xs,exists(xs,c,2));
            plot(xs,exists(xs,c,3));
            plot(xs,exists(xs,c,4));
            plot(xs,exists(xs,c,5));
            
            legend('T0\_unknown','T1\_unknown'...
                ,'T2\_unknown','T3\_unknown'...
                ,'T4\_unknown','T0\_known','T1\_known','T2\_known','T3\_known','T4\_known');
            title('window with FC');
            xlabel('Coordinates');ylabel('RPM');
        else
            plot(xs,RPMs(xs,c,1));       
            plot(xs,RPMs(xs,c,2));
            plot(xs,RPMs(xs,c,3));
            plot(xs,RPMs(xs,c,4));
            plot(xs,RPMs(xs,c,5));

            if c==2
                area(xs,0.1*coveredAreaTrue(xs,2));
            else
                area(xs,-0.1*coveredAreaTrue(xs,1));
            end
            legend('T0','T1','T2','T3','T4'...
                ,'known region');
            title('window with FC');
            xlabel('Coordinates');ylabel('RPM');
        end
        path='wiggle_windows/';
        saveas(gcf,[path,'region ',int2str(i),'.png']);
        close;
    end

    
    
end