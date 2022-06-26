%% Step 0 -> Create Leakage Traces
clear all
clc

load('AES_Sbox.mat');

leakageModel = "Direct"; % Set the power model for the leakage traces. Implement: HD.
profilingModel = "Direct"; % Set the profiling model of the templates. Implement: HD.
addedNoiseSNR = -30; % SNR for the added noise
keySpace = 256; % We attack an 8bit SBOX so we have a total of 256 keys.
numTraces = 5000; % Number of traces for profiling.
numAttTraces = 1000; % Number of traces for the attack phase.
kC = 123; % Real key of the attack;

iv = randi([0,keySpace-1],numTraces,1); % Draw a set of intermediate values.

if (leakageModel == "HW")
    h = sum(de2bi(iv),2);
    numClass = 9;
elseif (leakageModel == "Direct")
    h = iv;
    numClass = keySpace;
elseif (leakageModel == "HD")
    h = sum(de2bi(iv),2);
    v = zeros(size(h));
    v(1,:) = h(1,:);            % coping the first row
    for i = 2:size(v, 1)         % passing the rows
        for j = 1:size(v, 2)     % passing the coloums
            xor_h = bitxor(dec2bin(h(i,j),8)=='1',dec2bin(h(i-1,j),8)=='1');
            xor_d = double(xor_h+48);
            v(i,j) = sum( char(xor_d).' == '1');
        end
    end
    h = v;
    numClass = 9;
end

% Traces for profiling

Imat = repmat(h,1,20); % Add gaussian noise to the leakage model
Imat = [zeros(numTraces,30) Imat zeros(numTraces,30)]; % Pad the information with zeros
Imat = awgn(Imat,addedNoiseSNR); % Add gaussian noise
numSamples = size(Imat,2);

figure
plot(1:numSamples,Imat(1:100,:));
xlabel("Sample");
ylabel("Information Leakage");
title("Profiling Traces");

figure
plot(1:numSamples,mean(Imat));
xlabel("Sample");
ylabel("Mean Information Leakage");
title("Mean Profiling Traces");

% Traces for attack
plain = randi([0,keySpace-1],numAttTraces,1);
kC = repmat(kC,numAttTraces,1);
ivC = zeros(numAttTraces,1);
XORed = bitxor(plain,kC);
ivC = sbox(1+XORed);

if (leakageModel == "HW")
    yC = sum(de2bi(ivC),2);
    yC = yC';
elseif (leakageModel == "Direct")
    yC = ivC;
elseif (leakageModel == "HD")
    hC = sum(de2bi(ivC),2);
    hC = hC';
    vC = zeros(size(hC));
    vC(1,:) = hC(1,:);            % coping the first row
    for i = 2:size(vC, 1)         % passing the rows
        for j = 1:size(vC, 2)     % passing the coloums
            xor_h = bitxor(dec2bin(hC(i,j),8)=='1',dec2bin(hC(i-1,j),8)=='1');
            xor_d = double(xor_h+48);
            vC(i,j) = sum( char(xor_d).' == '1');
        end
    end
    yC = vC;
end

ImatAtt = repmat(yC',1,20); % Add gaussian noise to the leakage model
ImatAtt = [zeros(numAttTraces,30) ImatAtt zeros(numAttTraces,30)]; % Pad the information with zeros
ImatAtt = awgn(ImatAtt,addedNoiseSNR); % Add gaussian noise

figure
plot(1:numSamples,Imat);
xlabel("Sample");
ylabel("Information Leakage");
title("Attack Traces");

figure
plot(1:numSamples,mean(Imat));
xlabel("Sample");
ylabel("Mean Information Leakage");
title("Mean Attack Traces");

%% Step 1 -> CPA
% The CPA should be based on the leakage model.

hypo = [h randi([0,numClass-1],numTraces,keySpace-1)];
rho = corr(hypo,Imat);

figure
plot(1:numSamples,abs(rho));
xlabel("Sample");
ylabel("Correlation Coefficient");
title("CPA Attack on the Profiling Traces")

%% Step 1.1 -> CPA vs #Traces
tracesVector = 5:5:numTraces;
hypo = [h randi([0,numClass-1],numTraces,keySpace-1)];

corrVsTraces = calcCPAvsTraces(Imat,tracesVector,hypo,keySpace);

figure
plot(tracesVector,corrVsTraces,'Color',[0.7 0.7 0.7]);
hold on
plot(tracesVector,corrVsTraces(1,:),'Color',[0 0 0],'LineWidth',2);
xlabel("Number of Traces");
ylabel("Maximum Correlation");
title("Maximum Correlation vs Number of Traces for all Keys")

%% Step 2 -> Calculate SNR
% For this part we need to decide what profiling model to use
if (profilingModel == "HW")
    y = sum(de2bi(iv),2);
    numClassP = 9;
elseif (profilingModel == "Direct")
    y = iv;
    numClassP = keySpace;
elseif (profilingModel == "HD")
    h_p = sum(de2bi(iv),2);
    v_p = zeros(size(h_p));
    v_p(1,:) = h_p(1,:);            % coping the first row
    for i = 2:size(v_p, 1)         % passing the rows
        for j = 1:size(v_p, 2)     % passing the coloums
            xor_h = bitxor(dec2bin(h_p(i,j),8)=='1',dec2bin(h_p(i-1,j),8)=='1');
            xor_d = double(xor_h+48);
            v_p(i,j) = sum( char(xor_d).' == '1');
        end
    end
    y = v_p;
    numClassP = 9;
end

[SNR, POI] = calcSNR(Imat,y,numClassP,numSamples);

figure
plot(1:numSamples,SNR)
xlabel("Sample");
ylabel("Signal to Noise Ratio");
title(sprintf("SNR vs Samples. POI=%d",POI))

%% Step 2.1 -> Plot both Correlation and SNR
figure
yyaxis left
plot(1:numSamples,abs(rho));
xlabel("Sample");
ylabel("Correlation Coefficient");
yyaxis right
plot(1:numSamples,SNR)
xlabel("Sample");
ylabel("Signal to Noise Ratio");
title("SNR and CPA Combined")

%% Step 3 -> Template - Profiling
profile = zeros(numClassP,numClassP); % main matrix of POIs vs y
ycount = ones(numClassP,1); % For each y, we need to count how much POIs we have

for ii=1:numTraces
    profile(y(ii)+1,ycount(y(ii)+1)) = Imat(ii,POI); % for each input, classify the POI
    ycount(y(ii)+1) = ycount(y(ii)+1) + 1;
end

mu = zeros(numClassP,1); % mean of each y
for ii = 1:numClassP
    mu(ii) = sum(profile(ii,:))/(ycount(ii)-1);
end
profileNaN = profile;
profileNaN(profileNaN==0) = NaN; % Replace all zeros ("empty" cells) with NaN
sigma = std(profileNaN, 0, 2,'omitnan'); % Calculate the std of each y, while omitting the NaN
variance = sigma.^2;

figure
plot(1:numClassP,mu);
xlabel('y value');
ylabel('\mu');
title("\mu(Imat(POI)) vs y")


%% Step 4 -> Template - Attack
ImatAttPOI = ImatAtt(:,POI); % Take only the POI

prob_of_key=zeros(numAttTraces,keySpace); % probability for each input and key
for k = 0:255
    for q = 1:numAttTraces
       i = sbox(1+bitxor(plain(q),k)); % calc hypo
       if (profilingModel == "HW") % select the desired profile
           ya = sum(de2bi(i),2);
       elseif (profilingModel == "Direct")
           ya = i;
       elseif (profilingModel == "HD")
           h_i = sum(de2bi(i),2);
           v_i = zeros(size(h_i));
           v_i(1,:) = h_i(1,:);            % coping the first row
           for i = 2:size(v_i, 1)         % passing the rows
               for j = 1:size(v_i, 2)     % passing the coloums
                   xor_h = bitxor(dec2bin(h_i(i,j),8)=='1',dec2bin(h_i(i-1,j),8)=='1');
                   xor_d = double(xor_h+48);
                   v_i(i,j) = sum( char(xor_d).' == '1');
               end
           end
           ya = v_i;
        end
       prob_of_key(q,k+1) = normpdf(ImatAttPOI(q),mu(ya+1),sigma(ya+1)); % Calc probability
    end
end

prob_of_key(isnan(prob_of_key)) = 1; % take care of problems
prodMat = prod(prob_of_key,1);
[argVal, argIdx] = max(prod(prob_of_key,1));
kE = argIdx - 1;

figure
plot(1:keySpace,prodMat)
xlabel('Key Guess')
ylabel('Product of Probabilities')
title(sprintf('Product of Probabilities vs Key Guess. Recovered key is %d',kE))


% calculating the TR  

TR = zeros(1,numAttTraces);
for t = 1:numAttTraces
    logMat = sum(log2(prob_of_key(1:t,:)),1);
    [argVal, argIdx] = maxk(logMat,2);    % returns the two max values & indexs
    kE = argIdx(1) - 1;                   % best guess that is not the key
    if (kE == kC(1))                      % checking that is not the key
        kE = argIdx(2)-1;
    end

    TR(t) = logMat(kC(1)+1)./logMat(kE+1);   % TR
end

figure
plot(1:numAttTraces,TR)
xlabel('Number of Attack Traces')
ylabel('TR')
title(sprintf('TR vs Number of Attack Traces. Key used = %d',kC(1)))

