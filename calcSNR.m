function [SNR, POI] = calcSNR(Imat,y,numClassP,numSamples)

    E_Tc = zeros(numClassP,numSamples);
    var_Tc = zeros(numClassP,numSamples);

    for i = 0:numClassP-1
        x = find(y == i);            % find in y when i appers
        iii = size(x,1);
        if iii == 0
            E_Tc(i+1,:) = 0;
            var_Tc(i+1,:) = 0;
            continue
        end
        countt = Imat(x,:);          % all the rows corresponding to the indexs in x
        E_Tc(i+1,:)= mean(countt);   % signal_c = E(Tc) when Tc is all countt matrix
        var_Tc(i+1,:) = var(countt(:) - E_Tc(i+1,:));   % noise_c = var(Tc - E(Tc))
    end
    signal = var(E_Tc);
    noise = mean(var_Tc);
    SNR = signal ./ noise;
    [~,POI] = max(SNR,[],2);  % getting the POI (the maximum)

end