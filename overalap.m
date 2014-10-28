for i=1:9
    for j=1:10
     
        F=ForwardsSelectedAssets{i, j}
        
        mp=MPSelectedAssets{i, j};
        omp=OMPSelectedAssets{i, j};
        ls=LSOMPSelectedAssets{i, j};
        th=thrSelectedAssets{i,j};
        fmp(i, j)=nnz(ismember(F, mp))
        fomp(i, j)=nnz(ismember(F, omp))
        fls(i, j)=nnz(ismember(F, ls))
        fth(i, j)=nnz(ismember(F, th))
    end
    
end