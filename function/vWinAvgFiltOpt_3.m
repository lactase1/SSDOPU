function [outFrameWfilt] = vWinAvgFiltOpt_3(inFrame,inWeight,kRL,kRU,Nsec)
    %input prams
    % inFrame ==> 输入数据
    % inWeight ==> 置信�?
    if nargin == 4, Nsec = kRU - kRL +1; end
    %%
    [nZ,nX] = size(inFrame);
    krs = (kRU-1)/2;kcs = (kRU-1)/2;
    if krs == 0 && kcs ==0, outFrameWfilt = inFrame.*inWeight; return;end
    % 02_inFrame expansion : padding zeros or edge elements
    % padding zero
    inframe = zeros(nZ+krs*2,nX+kcs*2); inweight = inframe;
    inframe(krs+1:end-krs,kcs+1:end-kcs) = inFrame;
    inweight(krs+1:end-krs,kcs+1:end-kcs) = inWeight;
    %% from weigh to kernels
    wRg = linspace(0,1,Nsec);
      kRg = round(linspace(kRU,kRL,Nsec));
    %% Loop to compute weighting averaged results
    for iz = 1:nZ
        for ix = 1:nX
            wInds = wRg >= inWeight(iz,ix);
            pInds = find(wInds);
            fInd = pInds(1); fkrg = kRg(fInd);
            if mod(fkrg,2)
                ks1 = (fkrg-1)/2; ks2 = ks1;
            else
                ks1 = fkrg/2-1; ks2 = ks1+1;
            end
            cols = ix+kcs-ks1:ix+kcs+ks2;
            rows = iz+krs-ks1:iz+krs+ks2;
            inA = inframe(rows,cols);
            inB = inweight(rows,cols);
            outFrameWfilt(iz,ix) = sum(inA.*inB,'all');
        end
    end
end