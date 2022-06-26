function [corrVsTraces] = calcCPAvsTraces(Imat,tracesVector,hypo,keySpace)
% Correlation Coefficient Vs. Number of traces

corrVsTraces = zeros(keySpace,size(tracesVector,2));

for i = 1:size(tracesVector,2)
    R_i = corr(hypo(1:tracesVector(i),:),Imat(1:tracesVector(i),:));

    % max corr for each key with i traces
    corrVsTraces(:,i) = max(abs(R_i),[],2);
end

end