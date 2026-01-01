 function outFrameWfilt = vWinAvgFiltOpt(inFrame,inWeight,kRL,kRU,Nsec)
    % input params
    % inFrame  ==> ����
    % inWeight ==> ���Ŷ� (Ĭ��0~1)
    % kRL/kRU  ==> ��˹�˳ߴ緶Χ
    if nargin == 4, Nsec = kRU - kRL + 1; end

    [nZ,nX] = size(inFrame);
    inWeight = max(0,min(1,inWeight));
    inWeight(isnan(inWeight)) = 0;

    if kRU <= 1
        outFrameWfilt = inFrame .* inWeight;
        return;
    end

    % Ԥ���ɺ˳ߴ� (ǿ�ȵ� -> ��ˣ�ǿ�ȸ� -> С��)
    kRg = round(linspace(kRU,kRL,Nsec));
    kRg = kRg + mod(kRg+1,2);              % ����Ϊ�������ڶԳƴ���
    maxKer = max(kRg);
    maxRad = (maxKer - 1) / 2;

    % �߽���� replicate padding������ڱ�
    inframe = padarray(inFrame,[maxRad maxRad],'replicate','both');

    gaussHs = cell(maxKer,1);
    uniqSizes = unique(kRg);
    for idx = 1:numel(uniqSizes)
        sz = uniqSizes(idx);
        gaussHs{sz} = fspecial('gaussian',[sz sz],max(1,round(sz/2)));
    end

    outFrameWfilt = zeros(nZ,nX,'like',inFrame);
    for iz = 1:nZ
        for ix = 1:nX
            wVal = inWeight(iz,ix);
            fInd = min(Nsec,max(1,ceil(wVal * (Nsec - 1) + 1)));
            kerSize = kRg(fInd);
            rad = (kerSize - 1) / 2;

            rows = iz + maxRad - rad : iz + maxRad + rad;
            cols = ix + maxRad - rad : ix + maxRad + rad;
            roi = inframe(rows,cols);
            outFrameWfilt(iz,ix) = sum(roi .* gaussHs{kerSize},'all');
        end
    end
end