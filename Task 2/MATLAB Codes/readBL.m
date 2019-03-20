

combT = table();
lib_size=0; % used for RPM conversion
for len = [18:34 36]
    t = readtable([bl_file,int2str(len),'.txt']);
    t.Var2 = categorical(t.Var2);
    t.Var11 = categorical(t.Var11);
    t.Var12 = categorical(t.Var12);
    % filter   chr number            tRNA                   rRNA
    idx = (t.Var2 == chr_str) & (t.Var11 == 'NA') & (t.Var12 == 'NA');
    selT = t(idx,:); 
    combT = [combT;selT];
    lib_mRNAs = (t.Var11 == 'NA') & (t.Var12 == 'NA');
    lib_size=lib_size+sum(lib_mRNAs);
end