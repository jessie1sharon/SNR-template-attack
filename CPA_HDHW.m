%%%%%%%%%%%%%%%%%%% CPA %%%%%%%%%%%%%%%%%%%%
% First we load the data
key = 9;
exper = 9;
V = zeros(16,100);   % the V matrix
H = zeros(100,16);   % the H matrix
a = 1;               %counter for the matrix
load(strcat("Attack_Course_New\cmos_key"+ key + "_current_100traces_20files_exper" + exper + ".mat")); % Load the measurement file
load("Attack_Course_New\inputs_cmos_100traces_4bits.mat"); % Load the inputs vector


% Show the power traces
% x_axis = 1:size(Imat,2);
% figure
% plot(x_axis,Imat,'Color',[.8 .8 .8])
% hold on
% plot(x_axis,mean(Imat,1),'LineWidth',1,'color','black')
% xlabel("Sample")
% ylabel("Power")

% Key hypotheses
k = 0:15; % 4 bit key


% Cryptographic algorithm
% Here you need to pass each plaintext and key hypothesis through the
% cryptographic function to create the V matrix. Try to make it efficient.
% vij = SBOX4(pi XOR kij) for i=1..100, j=1..16
% Power model : Implement HD*HW

for i=1:100
    p = de2bi(plain(i),'left-msb');             % converting the plain to binary for xor
    while length(p)<4
            p = [false,p];
    end

    for j=1:16
        k_bin = de2bi(k(j),'left-msb');                % converting the key to binary for xor
        while length(k_bin)<4
            k_bin = [false,k_bin];
        end

        p_k_bin = bitxor(p,k_bin);            % the xor between the select key & plain
        p_k = bi2de(p_k_bin,'left-msb');               % back to dec
        V(a) = SBOX4(p_k);
        a = a +1 ;

    end
end

V = V';

for i=1:16           %passing the coloums
    n = 0;
    while n < size(V,1)-1   % rows
        if n == 0
            H(n+1,i) = sum( dec2bin(V(n+1,i)).' == '1' );  %if its the first in the coloum
        end
        n = n+1;
        v_pre = ~de2bi(V(n,i),'left-msb',4);       % not on the first word
        xor_h = bitand(dec2bin(V(n+1,i),4)=='1',v_pre);  % and with the next word
        xor_d = double(xor_h+48);
        H(n+1,i) = sum( char(xor_d).' == '1');

    end
end

% Correlation Coefficient
% The highest value in matrix R will represent the estimated key and the
% computation time.
R = corr(H,Imat);

%% graph a
figure(1)
max_row = max(abs(R),[],2);      %plot max for each row
figure(1)
plot(k,max_row);                 
xlabel("The Guess key");
ylabel("Maximum Correlation Coefficient");


%% graph b
figure(2)
for i=1:16                        % plot corr coef vs time
    x = (R(i,:));
    plot(x,'DisplayName',sprintf('key num %d',i-1));
    hold on
    legend()                     
end
xlabel("Sample(time)");
ylabel("Correlation Coefficient");
title("Correlation Coefficient Vs. Sample");


%% graph c
figure (3)
% corr coef vs num of traces
max_r_i = zeros(16,100);
max_r_k = zeros(1,100);
f_s_max = zeros(2,100);
indx_max = zeros(2,100);
graph_ = 'd';    % for the next graph
for i = 1:100                  % calc max corr for each key for each num of trace
    R_i = corr(H(1:i,:),Imat(1:i,:));
    
    % max corr for each key with i traces
    max_key_corr_i = max(abs(R_i),[],2);
    max_r_i(:,i) = max_key_corr_i;
    [dummy,max_r_k(i)] = max(max_key_corr_i,[],1);
    if(graph_ == 'd')
    % largest corr k != key
        if ((max_r_k(i)-1) == key)
            [f_s_max(:,i),indx_max(:,i)]= maxk(max_r_i(:,i),2);
            f_s_max(1,i) = f_s_max(2,i);
            indx_max(1,i) = indx_max(2,i);
        end
        if((max_r_k(i)-1) ~= key)
            f_s_max(1,i) = dummy;
            f_s_max(2,i) = 0;
            indx_max(1,i) = max_r_k(i);
            indx_max(2,i) = 0;
        end
    end
    if (i==100)
        curr_max = max_r_k(1);
        curr_int = 0;
        count_max = 1;
        count_int = 0;
        % find the longest repeated occurances of a key
        for j = 2:1:100
            if (curr_max == max_r_k(j))
                count_max = count_max +1;   
                continue;
            end
            if (curr_max ~= max_r_k(j) && (count_int == 0))
                curr_int = max_r_k(j);
                count_int = count_int + 1;
                continue;
            end
            if(curr_int == max_r_k(j) && (count_int ~= 0))
                count_int = count_int + 1;
                if(count_int > count_max)
                    count_max = count_int;
                    count_int = 0;
                    curr_max = curr_int;        % updeting the best key
                end
                continue;
            end
            if(curr_int ~= max_r_k(j) && (count_int >= 1))
                curr_int = max_r_k(j);
                count_int = 1;
                continue;
            end
        end
    end
end

plot(1:1:100,max_r_i,'+',1:1:100,max_r_i(curr_max,:),'k*',1:1:100,max_r_i(key+1,:),'r*');
xlabel("Number of Traces");
ylabel("Correlation Coefficient");
title(sprintf("Correlation Coefficient Vs. Number of Traces,Best key is %d",curr_max-1));


%% graph d
figure(4)
%CR vs num of traces
%only defenders can calc this cause we know the right key
second_max_corr = f_s_max(1,2:end);         % calc like before find max 
CR = max_r_i(key+1,2:end)./second_max_corr;

plot(2:1:100,CR,'+');
xlabel("Number of Traces");
ylabel("CR");
title(sprintf("CR Vs. Number of Traces, KC = %d",curr_max-1));





