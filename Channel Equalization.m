%% Meros B

clear; clc; close all;

N = 5000;
iter = 5;

%% Compute std of noise based on SNR

snr = 0:1:20;
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


%% LMS for channel 1

m = 0.035;
meanBer1 = zeros(length(snr),1);

for i = 1:length(snr)
    
    for k = 1:iter

        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);
        
        out1 = channel(h1, std(i), modulatedInp);
        out1 = [0; out1; 0];

        % LMS
        [pred1, error1, coef1] = lmsFun(out1, m, 1, N, symbols, modulatedInp);
        meanBer1(i) = meanBer1(i) + calcBer(inp, pred1);
    
    end
    
    meanBer1(i) = meanBer1(i)/iter;
    
end

%% LMS for channel 2

m = 0.035;
meanBer2 = zeros(length(snr),1);
    
for i = 1:length(snr)

    for k = 1:iter

        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);
        
        out2 = channel(h2, std(i), modulatedInp);
        out2 = [zeros(5,1); out2; zeros(5,1)];

        [pred2, error2, coef2] = lmsFun(out2, m, 5, N, symbols, modulatedInp);
        meanBer2(i) = meanBer2(i) + calcBer(inp, pred2);
        
    end
    
    meanBer2(i) = meanBer2(i)/iter;

end
    
%% LMS for channel 3

m = 0.035;
meanBer3 = zeros(length(snr),1);

for i = 1:length(snr)
    
    for k = 1:iter

        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);
        
        out3 = channel(h3, std(i), modulatedInp);
        out3 = [0; out3; 0];

        [pred3, error3, coef3] = lmsFun(out3, m, 1, N, symbols, modulatedInp);
        meanBer3(i) = meanBer3(i) + calcBer(inp, pred3);

    end
    
    meanBer3(i) = meanBer3(i)/iter;

end
    

%% Plots
Pe_theoretical = zeros(length(std),1); 
for i = 1:length(std)
    Pe_theoretical(i) = qfunc(sqrt(1/std(i)^2));
end

figure 
semilogy(snr,meanBer1,'Marker','o','Color','blue'); 
hold 
semilogy(snr,meanBer2,'Marker','+','Color','red'); 
semilogy(snr,meanBer3,'Marker','*','Color','cyan'); 
semilogy(snr(1:12),Pe_theoretical(1:12),'Marker','>','Color','magenta'); 
xlabel('SNR');
ylabel('BER')
title('Meros B - 2 - LMS, m = '+string(m))
legend('Channel 1','Channel 2', 'Channel 3','Theoritical')
grid on; 

%% DFE for channel 1

m = 0.035;
meanBer1 = zeros(length(snr),1);

len1 = 3;
len2 = 3;

for i = 1:length(snr)
    
    for k = 1:iter

        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);        
        
        out1 = channel(h1, std(i), modulatedInp);
        out1 = [zeros(len1-1,1); out1];

        % LMS
        [pred1, error1, c1, b1] = dfeFun(out1, m, len1, len2, N, symbols, modulatedInp);
        meanBer1(i) = meanBer1(i) + calcBer(inp, pred1);
    
    end
    
    meanBer1(i) = meanBer1(i)/iter;
    
end

%% DFE for channel 2

m = 0.01;
meanBer2 = zeros(length(snr),1);

len1 = 6;
len2 = 5;

for i = 1:length(snr)
    
    for k = 1:iter

        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);
        
        out2 = channel(h2, std(i), modulatedInp);
        out2 = [zeros(len1-1,1); out2];

        % LMS
        [pred2, error2, c2, b2] = dfeFun(out2, m, len1, len2, N, symbols, modulatedInp);
        meanBer2(i) = meanBer2(i) + calcBer(inp, pred2);
    
    end
    
    meanBer2(i) = meanBer2(i)/iter;
    
end

%% DFE for channel 3

m = 0.035;
meanBer3 = zeros(length(snr),1);

len1 = 5;
len2 = 3;

for i = 1:length(snr)
    
    for k = 1:iter

        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);        
        
        out3 = channel(h3, std(i), modulatedInp);
        out3 = [zeros(len1-1,1); out3];

        % LMS
        [pred3, error3, c3, b3] = dfeFun(out3, m, len1, len2, N, symbols, modulatedInp);
        meanBer3(i) = meanBer3(i) + calcBer(inp, pred3);
    
    end
    
    meanBer3(i) = meanBer3(i)/iter;
    
end

%% Plots
Pe_theoretical = zeros(length(std),1); 
for i = 1:length(std)
    Pe_theoretical(i) = qfunc(sqrt(1/std(i)^2));
end

figure 
semilogy(snr,meanBer1,'Marker','o','Color','blue'); 
hold 
semilogy(snr,meanBer2,'Marker','+','Color','red'); 
semilogy(snr,meanBer3,'Marker','*','Color','cyan'); 
semilogy(snr(1:12),Pe_theoretical(1:12),'Marker','>','Color','magenta'); 
xlabel('SNR');
ylabel('BER')
title('Meros B - 2 - DFE')
legend('Channel 1, m = 0.035','Channel 2, m = 0.01', 'Channel 3, m = 0.035', ...
    'Theoritical', 'Location','Southwest')
grid on; 

%% Question 3

iter = 50;

m = [0.1; 0.035; 0.0008];
sqErr = zeros(3,N/2);

figure
hold
for i = 1:length(m)
    
    for k = 1:iter
        
        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);

        out2 = channel(h2, std(21), modulatedInp);
        out2 = [zeros(5,1); out2; zeros(5,1)];
        
        [pred2, error2, coef2] = lmsFun(out2, m(i), 5, N, symbols, modulatedInp);
%         calcBer(inp, pred2)
       
        sqErr(i,:) = sqErr(i,:) + (error2(1:2500).^2)';
    end
    plot(abs(sqErr(i,:))/iter)
end
xlabel('Iterations');
ylabel('MSE')
title('Meros B - 3')
legend('m = ' + string(m(1)),'m = ' + string(m(2)), 'm = ' + string(m(3)))

%% Question 4

iter = 10;

m = 0.01;
meanBer21 = zeros(length(snr),1);
meanBer22 = zeros(length(snr),1);

len1 = 6;
len2 = 5;

for i = 1:length(snr)
    
    for k = 1:iter

        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);        
        
        out2 = channel(h2, std(i), modulatedInp);
        out2 = [zeros(len1-1,1); out2];

        [pred2, error2, c2, b2] = dfeFun(out2, m, len1, len2, N, symbols, modulatedInp);
        meanBer21(i) = meanBer21(i) + calcBer(inp, pred2);
        
        [pred2, error2, c2, b2] = dfeFun2(out2, m, len1, len2, N, symbols, modulatedInp);
        meanBer22(i) = meanBer22(i) + calcBer(inp, pred2);
    
    end
    
    meanBer21(i) = meanBer21(i)/iter;
    meanBer22(i) = meanBer22(i)/iter;
    
end

figure 
semilogy(snr,meanBer21,'Marker','o','Color','blue'); 
hold 
semilogy(snr,meanBer22,'Marker','+','Color','red'); 
xlabel('SNR');
ylabel('BER')
title('Meros B - 4')
legend('Channel 2','Channel 2 - Correct Symbols')
grid on; 

%% Question 5

m = 0.035;

inp = randi([0,1],N,1);
modulatedInp = modulate(inp, symbols);

out2 = channel(h2, std(21), modulatedInp);
out2 = [zeros(5,1); out2; zeros(5,1)];

[pred2, error2, coef2] = lmsFun(out2, m, 5, N, symbols, modulatedInp);

figure
freqz(h2)
title('Channel')

figure
freqz(coef2)
title('Equalizer')

figure
plot(abs(conv(h2, coef2)))
title('Convolution between channel and equalizer')
grid on

%% Question 6

iter = 50;
h2 = [0.04; -0.05; 0.07; -0.21; -0.5; 0.72; 0.36; 0; 0.21; 0.03; 0.07];
h2alt = [0.04; -0.05; 0.07; -0.21; -0.5; -0.72; 0.36; 0; 0.21; 0.03; 0.07];

m = [0.1; 0.035; 0.008];
sqErr = zeros(3,N/2);

figure
hold
for i = 1:length(m)
    
    for k = 1:iter
        
        inp = randi([0,1],N,1);
        modulatedInp = modulate(inp, symbols);

        out2 = channel2(h2, h2alt, std(21), modulatedInp);
        out2 = [zeros(5,1); out2; zeros(5,1)];
        
        [pred2, error2, coef2] = lmsFun(out2, m(i), 5, N, symbols, modulatedInp);
%         calcBer(inp, pred2)
       
        sqErr(i,:) = sqErr(i,:) + (error2(1:2500).^2)';
    end
    plot(abs(sqErr(i,:))/iter)
end
xlabel('Iterations');
ylabel('MSE')
title('Meros B - 6')
legend('m = ' + string(m(1)),'m = ' + string(m(2)), 'm = ' + string(m(3)))

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

function [pred, error, coef] = lmsFun(inp, m, len, N, symbols, modulatedInp) 

    coef = [zeros(len,1); 1; zeros(len,1)];
    pred = zeros(N,1);
    error = zeros(N/2,1);

    for j = len+1:N/2+len

        zk = coef'*inp(j-len:j+len);
        pred(2*(j-len-1)+1:2*(j-len)) = maxLikelihood(zk, symbols);

        if j <= 250
            ak = modulatedInp(j-len);
        else
            ak = modulate(pred(2*(j-len-1)+1:2*(j-len)), symbols);
        end

        error(j) = ak - zk;
        coef = coef + m*inp(j-len:j+len)*conj(error(j));

    end
        
end

function [pred, error, c, b] = dfeFun(inp, m, len1, len2, N, symbols, modulatedInp) 

    c = [zeros(len1-1,1); 1];
    b = zeros(len2,1);
    decisions = zeros(len2,1);
    
    pred = zeros(N,1);
    error = zeros(N/2,1);
    
    for i = len1:N/2+len1-1
        
        zk = c'*inp(i-len1+1:i) + b'*decisions; 
        pred(2*(i-len1)+1:2*(i-len1+1)) = maxLikelihood(zk, symbols);
       
        if i <= 250 + len1
            ak = modulatedInp(i-len1+1);
        else
            ak = modulate(pred(2*(i-len1)+1:2*(i-len1+1)), symbols);
        end
        
        error(i) = ak - zk;
        c = c + m*inp(i-len1+1:i)*conj(error(i));
        b = b + m*decisions*conj(error(i));
        
        decisions = circshift(decisions,-1);
        decisions(end) = ak;
    end

end

function [pred, error, c, b] = dfeFun2(inp, m, len1, len2, N, symbols, modulatedInp) 

    c = [zeros(len1-1,1); 1];
    b = zeros(len2,1);
    decisions = zeros(len2,1);
    
    pred = zeros(N,1);
    error = zeros(N/2,1);
    
    for i = len1:N/2+len1-1
        
        zk = c'*inp(i-len1+1:i) + b'*decisions; 
        pred(2*(i-len1)+1:2*(i-len1+1)) = maxLikelihood(zk, symbols);
       
        ak = modulatedInp(i-len1+1);
        
        error(i) = ak - zk;
        c = c + m*inp(i-len1+1:i)*conj(error(i));
        b = b + m*decisions*conj(error(i));
        
        decisions = circshift(decisions,-1);
        decisions(end) = ak;
    end

end


function out = channel2(h, halt, std, inp)

    out = zeros(length(inp),1);
    m = floor(length(h)/2);
    
    % Channel
    inp = [zeros(m,1); inp; zeros(m,1)];
    for i = m+1:length(inp)-m
        if i <= 1250
            out(i-m) = h'*inp(i-m:i+m);
        else
            out(i-m) = halt'*inp(i-m:i+m);
        end
    end

    
    % Add Noise
    out = out + std*randn(length(out),1) + std*1j*randn(length(out),1);
    
 
end

