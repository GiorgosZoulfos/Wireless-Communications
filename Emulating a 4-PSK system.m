%% MEROS A 

clc; clear; close all

N = 10^5;
iter = 5;


%% Compute std of noise based on SNR

snr = -5:1:10;
std = sqrt(0.5.*(10.^(-snr./10)));

%% Define symbols

s1 = 1 + 0j;
s2 = 0 + 1j;
s3 = -1 + 0j;
s4 = 0 - 1j;

symbols = [s1; s2; s3; s4];

%% Define channel responses

h1 = 1; 
h2 = [0.04; -0.05; 0.07; -0.21; -0.5; 0.72; 0.36; 0; 0.21; 0.03; 0.07];
h3 = [0.407; 0.815; 0.407]; 

%% 10^5 random bits and their modulation

inp = randi([0,1],N,1);
modulatedInp = modulate(inp, symbols);


%% The three channels

meanBer1 = zeros(length(std),1);
meanBer2 = zeros(length(std),1);
meanBer3 = zeros(length(std),1);

for i = 1:length(std)
    
    ber1 = 0; ber2 = 0; ber3 = 0;
    
    for j = 1:iter
        
        
        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);
        
        out1 = channel(h1, std(i), modulatedInp);
        out2 = channel(h2, std(i), modulatedInp);
        out3 = channel(h3, std(i), modulatedInp);

        pred1 = maxLikelihood(out1, symbols);
        pred2 = maxLikelihood(out2, symbols);
        pred3 = maxLikelihood(out3, symbols);

        ber1 = ber1 + calcBer(inp, pred1);
        ber2 = ber2 + calcBer(inp, pred2);
        ber3 = ber3 + calcBer(inp, pred3);
        
    end
    
    meanBer1(i) = ber1/iter;
    meanBer2(i) = ber2/iter;
    meanBer3(i) = ber3/iter;
    
end

%% Plots 
Pe_theoretical = zeros(length(std),1); 
for i = 1:1:length(std)
    Pe_theoretical(i) = qfunc(sqrt(1/std(i)^2));
end

figure(); 
semilogy(snr,meanBer1,'Marker','o','Color','blue'); 
hold on;
semilogy(snr,meanBer2,'Marker','+','Color','red'); 
hold on;
semilogy(snr,meanBer3,'Marker','*','Color','cyan'); 
hold on; 
semilogy(snr,Pe_theoretical,'Marker','>','Color','magenta'); 
xlabel('SNR');
ylabel('BER')
title('Error Probability')
legend('Channel 1','Channel 2', 'Channel 3','Theoritical','Location','southwest')
grid on; 

%% Functions used above

function out = modulate(inp, symbols)

    out = zeros(length(inp)/2,1);

    for i = 1:2:length(inp)
        
        if inp(i:i+1) == [0;0]
            out((i+1)/2) = symbols(1);
        elseif inp(i:i+1) == [0;1]
            out((i+1)/2) = symbols(2);
        elseif inp(i:i+1) == [1;1]
            out((i+1)/2) = symbols(3);
        else 
            out((i+1)/2) = symbols(4);
        end
        
    end
        

end
    

function out = channel(h, std, inp)

    out = zeros(length(inp),1);
    m = floor(length(h)/2);
    
    % Channel
    inp = [zeros(m,1); inp; zeros(m,1)];
    for i = m+1:length(inp)-m
        out(i-m) = h'*inp(i-m:i+m);
    end

    
    % Add Noise
    out = out + std*randn(length(out),1) + std*1j*randn(length(out),1);
    
 
end

function out = maxLikelihood(inp, symbols)

    out = zeros(2*length(inp),1);
    
    for i = 1:length(inp)
        
        dist = zeros(4,1);
        for s = 1:length(symbols)
           dist(s) = norm(symbols(s) - inp(i)); 
        end
        
        [~, ind] = min(dist);
        if ind == 1
            out(2*(i-1)+1:2*i) = [0;0];
        elseif ind == 2
            out(2*(i-1)+1:2*i) = [0;1];
        elseif ind == 3
            out(2*(i-1)+1:2*i) = [1;1];
        else
            out(2*(i-1)+1:2*i) = [1;0];
        end
        
    end
    
end

function out = calcBer(inp1, inp2)

    mask = inp1 == inp2;
    out = 1 - sum(mask)./length(mask);

end


