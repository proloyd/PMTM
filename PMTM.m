clear all;
close all;

N = 128;  % number of spectrum samples in [0,fs/2) determining the 
          % frequency spacing fs/(2N) between samples  
          
N_max = 128;  % number of desired spectrum samples which tends to be much 
              % lower then N in neural signals because of oversampling   
             
K = 512;  % number of spiking samples per neuron

L = 40;  % number of neurons

iter_Newton = 100;  % number of newton iterations per EM iteration

iter_EM = 100;  % number of EM iterations

tol_Newton = 1e-50; % Stopping criterion for Newton's method

%*************************Toy Example generation***************************
% AR process generation
B = 1; 
multiplier = 0.006;
Coeff = conv( ...
    conv([1 -0.965*exp(-2*1i*pi*0.10)], [1 -0.965*exp(2*1i*pi*0.10)]), ...
    conv([1 -0.975*exp(-2*1i*pi*0.35)], [1 -0.975*exp(2*1i*pi*0.35)]) ...
    );  % AR coeffs

[H_AR,W] = freqz(multiplier, Coeff, 256);  % Frequency response 
                                           % of the filter

b = randn(K+512, 1);                                
x_AR = filter(B, Coeff, multiplier * (b));       
y2 = x_AR(end - K + 1:end);
CIF = (y2 + 0.04);  % baseline spiking rate = 1.2

% Spiking data generation by Thinning
spikes = zeros(length(y2), L);
for i = 1:L
    spikes(:, i) = (CIF > rand(length(y2), 1));
    rng('shuffle');
end

% Plotting the data Set
figure, 
subplot(2, 1, 1), plot(CIF);
xlabel('$$k$$', 'Interpreter', 'Latex');
ylabel('Amplitude','Interpreter','Latex');
xlim([200.5, 350.5]);
xlim([1,K]);
subplot(2, 1, 2), SpikeRasterPlot(spikes(:, 1:10)');
xlabel('$$k$$', 'Interpreter','Latex');
ylabel('Trials','Interpreter','Latex');
xlim([200.5, 350.5]);
%% **************************************************************************

%*****************************PMTM ESTIMATION****************************** 

% Discrete Prolate Spheroidal Sequences
time_halfbandwidth = 5;
dps_seq= dpss(K, time_halfbandwidth);

% Initialization
S_mt = [];                                  
S_pp = [];                                    
S_ss = [];                                  
S_pp_mt = [];                                


A=zeros(K, 2*N_max);
for i=1:K
    for j=1:N_max
        A(i,2*j-1) = cos(i*pi*(j-1)/N);
        A(i,2*j) = -sin(i*pi*(j-1)/N);
    end
end
A = 2*pi*A/N;
A(:, 2) = [];
Ainv = inv(A'*A) * A';


PSTH = mean(spikes, 2); % PSTH calculation

mu_hat = mean(PSTH); % Estimate of true mu

gamma = 1e-3;
   
childPSTH = [];
for k = 1:size(dps_seq, 2) - 5
    % ## Auxiliary spike generation ##
    lamda0 = 1.5 * max(dps_seq(:, k));  % factor 1.5 is arbitrary, can 
                                        % use any factor > 1.
    n_avg = zeros(K, 1);
    for j = 1:L
         n = spikes(:, j);
         n_comp = 1 - spikes(:, j);
         n = n .* dps_seq(:, k) / lamda0;
         n_comp = - n_comp .* dps_seq(:, k) / lamda0;
         n = n .* (n > 0);
         n_comp = n_comp.*(n_comp > 0);
         n_avg  = n_avg + (n + n_comp);
    end
    n_avg = n_avg / L;
    mu_child = mu_hat * (dps_seq(:, k) / lamda0) .* (dps_seq(:, k) >= 0) ...
        + (1 - mu_hat) * (-dps_seq(:, k) / lamda0) .* (dps_seq(:, k) <= 0); 
    ub = max(mu_child) / norm(A(2,:));
    b =  (- Ainv * mu_child);
    
    % Plot Auxiliary Spike trains    
    figure, stem(0:K-1, n_avg);
    str = sprintf('$$%d$$ th child process', k-1);
    title(str, 'Interpreter', 'Latex');
    xlabel('$$k$$', 'Interpreter', 'Latex');
    ylabel('$$\sum_{l=1}^{L-1} n_k^{(l)}$$', 'Interpreter', 'Latex');
    drawnow
    
    % ## EM algorithm ##
    childPSTH = [childPSTH, n_avg];
    mu_v = zeros(2*N_max-1, 1);
    theta = zeros(2*N_max-1, iter_EM);
    PSD_est = zeros(N_max, iter_EM);
    
    % initialize theta for EM
    theta(:, 1) = 1e-0*ones(2*N_max-1, 1);  % 1e-1*ones(2*N_max-1, 1); before
    PSD_est(2:end, 1) = theta(2:2:end, 1) + theta(3:2:end, 1);

    % initialize Newton
    x = mu_child + 0.5 * zeros(K, 1);
    % x = mu + SS_latent_1(n_avg);
    g = L * A' * (childPSTH(:, k) ./ x - (1 - childPSTH(:, k)) ./ (1 - x));
    H = -L * A' * diag(childPSTH(:,k)./(x.^2) ...
        + (1 - childPSTH(:, k)) ./ ((1 - x) .^ 2)) * A ...
        - diag(1 ./ theta(:, 1));
%     d = H\g;
    d = -g;
    
    lambda1 = zeros(size(A,1),1);
    lambda2 = zeros(size(A,1),1);
    tau = 1e-1;
    
    % EM iterations
    for i=2:iter_EM
        for m = 1:10
        % Newton iterations
        for j=1:iter_Newton
%             tau = 1;
%             while (sum(x - tau*A*d + 0.02 <= 0)~=0 || ...
%                     sum(x - tau*A*d - 0.00 >= 1)~=0)
%                 tau = tau / 2;
%                 if tau < 1e-30
%                     break
%                 end
%             end
%             mu_v = mu_v - tau * d;
            mu_v = mu_v - d;
%             mu_v(mu_v < b) = b(mu_v < b);
%             mu_v(mu_v > ub) = ub;
%             mu_v = mu_v / norm(mu_v) * min(norm(mu_v), ub);
        
            x = mu_child + A * mu_v;
%             g = L * A' *(childPSTH(:, k) ./ x - ...
%                 (1 - childPSTH(:, k)) ./ (1 - x)) - mu_v ./ theta(:, i - 1);             
            g = L * A' *(childPSTH(:, k) ./ x - ...
                (1 - childPSTH(:, k)) ./ (1 - x)) - mu_v ./ theta(:, i - 1)...
                -A'*(lambda2-lambda1);
            H = - L * A' * diag(childPSTH(:, k) ./ (x .^ 2) ...
                + (1 - childPSTH(:, k)) ./ ((1 - x) .^ 2)) * A ...
                - diag(1 ./ theta(:, i - 1));
            d = H \ g;
%             d = -g;
%             if - (g' * d) * tau < tol_Newton 
%                 break
%             end
            
            if - (g' * d) < tol_Newton 
                break
            end
        end
        update = A * mu_v;
        lambda1 = lambda1 - tau * update; 
        lambda2 = lambda2 + tau * update;
        lambda1(lambda1 < 0) = 0;
        lambda2(lambda2 < 0) = 0;
        end
 
        E = diag(-H \ eye(2 * N_max - 1)) + mu_v .^ 2;
%         theta(:,i) = E;  % Without regularization works!
        theta(:, i) = (-1 + sqrt(1 + 8 * gamma * E)) / (4 * gamma);  % But 
        % a little regularization helps fast convergence
        PSD_est(2:end, i) = theta(2:2:end, i) + theta(3:2:end, i);
    end
    
    % % Estimation progress
    figure, pcolor(PSD_est);
    shading flat;
    xlabel('Number of Iterations$$\rightarrow$$', 'Interpreter', 'Latex');
    ylabel('Frequency$$\rightarrow$$', 'Interpreter', 'Latex');
    str = [sprintf('$$%d$$th tapered auxiliary process:' , k-1) ...
           'estimation progress'];
    title(str, 'Interpreter', 'Latex');

    % Eigen-Spectra
    S_pp_mt = [S_pp_mt, PSD_est(:, iter_EM) * ((2* pi * lamda0) ^ 2)];
    
end
pp_mt_est = mean(S_pp_mt, 2); % PMTM PSD

%**************************************************************************

%%
%************************Oracle,PSTH,SS PSD********************************

% SS CIF estimation
y_ss = SS_latent_estimation(spikes);

for k = 1:size(dps_seq,2) - 2
    % Oracle PSD
    S_mt = [S_mt, abs(fft((y2 - mean(y2)) .* dps_seq(:, k))) .^ 2];
    
    % PSTH PSD
    S_pp = [S_pp, (abs(fft((PSTH - mean(PSTH)) .* dps_seq(:, k))) .^ 2)]; 
    
    % SS PSD
    S_ss = [S_ss, abs(fft((y_ss - mean(y_ss)) .* dps_seq(:, k))) .^ 2]; 
end

mt_est = mean(S_mt, 2); % Oracle PSD
pp_est = mean(S_pp, 2); % PSTH PSD
pp_ss_est = mean(S_ss, 2); % SS PSD

%**************************************************************************

%%
%********************************Overlay plot******************************
figure, 
plot((1:N_max) / (2 * N), 10 * log10(pp_mt_est * K / N), ...
    'k', 'LineWidth', 1.2);
hold on; 
plot((0:(length(pp_est) - 1)) / length(pp_est), 10 * log10(pp_est), ...
    '--g', 'LineWidth', 1.2);
plot((0:(length(mt_est) - 1)) / length(mt_est), 10 * log10(mt_est), ...
    '--r', 'LineWidth', 1.2);
plot((0:(length(pp_ss_est) - 1)) / length(pp_ss_est), ...
    10 * log10(pp_ss_est), '-.c', 'LineWidth', 1.2);
plot(W / (2 * pi), 20 * log10(abs(H_AR)), 'b', 'LineWidth', 1.2);
legend('Proposed Method','PSTH PSD','Oracle PSD','SS PSD','True PSD');
xlim([0 0.5]);
grid on
xlabel('Normalized Frequency', 'Interpreter', 'Latex');
ylabel('Power/Frequency $$(dB/rad/sample)$$', 'Interpreter', 'Latex');

%*********************************MSE calculation**************************
% lin_MSE_pp = sum(((pp_est(2:256) - abs(H_AR(2:end)) .^ 2) .^2 ) ...
%     ./ (abs(H_AR(2:end)) .^ 2));
% lin_MSE_pp_mt = sum(((pp_mt_est(2:end) - abs(H_AR(2:end)) .^ 2) .^ 2) ...
%     ./ (abs(H_AR(2:end)) .^ 2));
% lin_MSE_pp_ss = sum(((pp_ss_est(2:256) - abs(H_AR(2:end)) .^ 2) .^ 2) ...
%     ./ (abs(H_AR(2:end)) .^ 2));
% 
% fprintf('Method \t \t MSE\n');
% fprintf('PSTH-PSD \t %f\n', lin_MSE_pp);
% fprintf('PMTM-PSD \t %f\n', lin_MSE_pp_mt);
% fprintf('SS-PSD \t \t %f\n', lin_MSE_pp_ss);