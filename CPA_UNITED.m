function [key_guess] = CPA_UNITED(key,exper,PM,graph,num_t)

% First we load the data

load(strcat("Attack_Course_New\cmos_key"+ key + "_current_100traces_20files_exper" + exper + ".mat")); % Load the measurement file
load("Attack_Course_New\inputs_cmos_100traces_4bits.mat"); % Load the inputs vector

V = zeros(16,100);   % the V matrix
H = zeros(100,16);   % the H matrix

% Key hypotheses
k = 0:15; % 4 bit key
m = 1;    %counter for the matrix

%Power model
%Here we need to create the matrix H from V. Implement HW, HD and HD*HW.
if (PM == 'HW')
    H = zeros(16,100);   % the H matrix
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
            V(m) = SBOX4(p_k);
            H(m) = sum( dec2bin(V(m)).' == '1' );
            m = m +1 ;
    
        end
    end
    
    H = H';
    V = V';

elseif (PM == 'HD')
    for i=1:100
        p = de2bi(plain(i),'left-msb');           % converting the plain to binary for xor
        while length(p)<4
                p = [false,p];
        end
    
        for j=1:16
            k_bin = de2bi(k(j),'left-msb');       % converting the key to binary for xor
            while length(k_bin)<4
                k_bin = [false,k_bin];
            end
    
            p_k_bin = bitxor(p,k_bin);            % the xor between the select key & plain
            p_k = bi2de(p_k_bin,'left-msb');      % back to dec
            V(m) = SBOX4(p_k);
            m = m +1 ;
    
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
            xor_h = bitxor(dec2bin(V(n,i),4)=='1',dec2bin(V(n+1,i),4)=='1');
            xor_d = double(xor_h+48);
            H(n+1,i) = sum( char(xor_d).' == '1');
    
        end
    end

elseif (PM == 'HU')
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
            V(m) = SBOX4(p_k);
            m = m +1 ;
    
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
end
H = H(1:num_t,:);
Imat = Imat(1:num_t,:);
R = corr(H,Imat);  

% graph a
if (graph == 'a')
    [max_row,~] = max(abs(R),[],2);
    [~,max_corr] = max(max_row,[],1);
    key_guess = max_corr - 1;

%     figure
%     max_row = max(abs(R),[],2);      %plot max for each row
%     plot(k,max_row);                
%     hold on
%     xlabel("The Guess key");
%     ylabel("Maximum Correlation Coefficient");

end
% graph b
if (graph == 'b')

    figure
    for i=1:16                        % plot corr coef vs time
        x = (R(i,:));
        hold on
        plot(x,'DisplayName',sprintf('key num %d',i-1));
        hold off
        legend()                     
    end
    xlabel("Sample(time)");
    ylabel("Correlation Coefficient");
    title("Correlation Coefficient Vs. Sample");

end

% graph c
if(graph == 'c' || graph == 'd')
    figure 
    % corr coef vs num of traces
    max_r_i = zeros(16,100);
    max_r_k = zeros(1,100);
    f_s_max = zeros(2,100);
    indx_max = zeros(2,100);
    
    for i = 1:100                  % calc max corr for each key for each num of trace
        R_i = corr(H(1:i,:),Imat(1:i,:));
        
        % max corr for each key with i traces
        max_key_corr_i = max(abs(R_i),[],2);
        max_r_i(:,i) = max_key_corr_i;
        [dummy,max_r_k(i)] = max(max_key_corr_i,[],1);
        if(graph == 'd')
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

end
% graph d
if(graph == 'd')
    figure
    %CR vs num of traces
    %only defenders can calc this cause we know the right key
    second_max_corr = f_s_max(1,2:end);         % calc like before find max 
    CR = max_r_i(key+1,2:end)./second_max_corr;
    
    plot(2:1:100,CR,'+');
    xlabel("Number of Traces");
    ylabel("CR");
    title(sprintf("CR Vs. Number of Traces, KC = %d",curr_max-1));

end

end