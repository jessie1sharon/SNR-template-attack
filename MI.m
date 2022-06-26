%%%%%%%%%%%%%%%% MI %%%%%%%%%%%%%%%%%%%

graph = 'a';
PM = 'HD';
guess_t = zeros(100,16,50);
key_guess_m = zeros(16,50);
success_rate = zeros(16,100);
JointP = zeros(16,16,100);             % P(Key,Key_Guess)

for num_t = 1:100
        for key = 0:1:15
            for exper = 1:1:50
                % Load the measurement file
                key_guess_m(key+1,exper) = CPA_UNITED(key,exper,PM,graph,num_t);
            end
        end
    guess_t(num_t,:,:) = key_guess_m(:,:);
    fprintf("traces = %d",num_t);


    % joint prob
    
    for i = 1:16
        for j = 0:15
            JointP(i,j+1,num_t) = (1/16)*(sum(key_guess_m(i,:)==j)/50);
        end
    end
    
    %calc CPA success rate Vs. Num of traces per key 
    sum_row = sum(JointP(:,:,num_t),2);
    success_rate(:,num_t) = diag(JointP(:,:,num_t))./(sum_row);
end

% figure
% plot(1:100,success_rate,'+');
% xlabel("Number of Traces");
% ylabel("Success Rate");
% title("Success Rate Vs. Number of Traces for each key");

% mutual information

MI_t = zeros(1,100);
for num_t = 1:100

    mi = 0;
    p_k = sum(JointP(:,:,num_t),2);
    p_k_1 = sum(JointP(:,:,num_t),1);
    for k1 = 1:16
        for k2 = 1:16
            pr = JointP(k1,k2,num_t);
            if (pr == 0)
                continue;
            end
            pr1 = p_k(k1)*p_k_1(k2);
            div = pr * abs(log2(pr/pr1));
            mi = mi + div;
        end
    end
    MI_t(1,num_t) = mi;
end

figure
plot(1:100,MI_t,'+');
xlabel("Number of Traces");
ylabel("MI");
title("MI Vs. Number of Traces");



