%% refFlat info
% This file use all 0 based



refFlat = readtable('refFlat_240118.txt');
refFlat.Var3 = categorical(refFlat.Var3);       % chr number
idx = (refFlat.Var3 == chr_str);
refFlat = refFlat(idx,:); 

refFlat.Var4 = categorical(refFlat.Var4);       % sense
coverArea = zeros(total_length22,2);              % (:,1) for +ve sense
                                                % (:,2) for -ve sense,       
for i = 1:height(refFlat)
    num= refFlat.Var9(i);                       % num of pieces of exons
    starts = refFlat.Var10{i};   
    ends = refFlat.Var11{i}; 
    sense = refFlat.Var4(i); 
    startStr = split(starts,',');
    endStr = split(ends,',');
    for j = 1:num
        st = str2double(startStr{j});
        ed = str2double(endStr{j})-1;
        for k = st:ed
            if sense == '+'
                coverArea(k,1) = coverArea(k,1) + 1;
            else
                coverArea(k,2) = coverArea(k,2) + 1;
            end
        end
    end
end

coveredAreaTrue = coverArea~=0;
