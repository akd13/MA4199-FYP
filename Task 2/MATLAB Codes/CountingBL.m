%% Counting
% this file uses all 0 based. 



% left = min([combT.Var9;combT.Var10]);
% right = max([combT.Var9;combT.Var10]);

countPos = zeros(total_length22,2);  % (:,1) for +ve sense
                            % (:,2) for -ve sense, 

for i = 1:height(combT)
    small = min(combT.Var9(i) , combT.Var10(i));
    large = max(combT.Var9(i) , combT.Var10(i));
    if combT.Var9(i) < combT.Var10(i)   % +ve sense
        for read = small:large
           countPos(read,1) = countPos(read,1)+1; 
        end
    else 
        for read = small:large
           countPos(read,2) = countPos(read,2)-1; 
        end
    end

end