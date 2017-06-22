function [ x_SS ] = SS_latent_estimation( spikes ) 
% SS_x returns the estimated latent process under SS modelling

N = length(spikes);
n = mean(spikes,2);

EM_iterations = 200;
sig_EM_init = 0.5;
x_init = 0.0;
sig_x_init = 1;
mu = mean(n);

x_k_k = zeros(N+1,1);
x_k_k1 = zeros(N+1,1);
x_k_K = zeros(N+1,1);
sig_k_k = zeros(N+1,1);
sig_k_k1 = zeros(N+1,1);
sig_k_K = zeros(N+1,1);
sig_k_k1_K = zeros(N+1,1);

x_k_k(1) = x_init;
sig_k_k(1) = sig_x_init;

sig = zeros(EM_iterations,1);
sig(1) = sig_EM_init;

alpha = 0.95;

for i = 1:EM_iterations-1
    
    for j = 2:N+1
        x_k_k1(j) = alpha*x_k_k(j-1);
        sig_k_k1(j) = alpha^2*sig_k_k(j-1)+sig(i);
        r = 0;
        f = @(x) x - x_k_k1(j) - 10*sig_k_k1(j)*(  n(j-1)/(mu + x) - (1-n(j-1))/(1 - mu - x)  ); 
        df = @(x) 1 + 10*sig_k_k1(j)*(  n(j-1)/(mu + x)^2 + (1-n(j-1))/(1 - mu - x)^2  ); 
        d = f(r)/df(r);
        beta = 1;
        for itr = 1:10
            while r - beta*d > 1-mu || r - beta*d  < -mu
                beta = beta/2;
            end
            r1 = r - beta * d;
            if r1 - r  <  10^-5
                break;
            end
            r = r1;
            d = f(r)/df(r);
            beta = 1;
        end   
        x_k_k(j) = r;
        sig_k_k(j) = 1/( 10*n(j-1)/(mu + x_k_k(j))^2 + 10*(1-n(j-1))/(1- mu - x_k_k(j))^2 + 1/sig_k_k1(j) );
    end
    
    x_k_K(N+1) = x_k_k(N+1);
    sig_k_K(N+1) = sig_k_k(N+1);
    
    for z = N:-1:1
        A = alpha*sig_k_k(z)/sig_k_k1(z+1);
        x_k_K(z) = x_k_k(z)+A*(x_k_K(z+1)-x_k_k1(z+1));
        sig_k_K(z) = sig_k_k(z)+(A^2)*(sig_k_K(z+1)-sig_k_k1(z+1));
        sig_k_k1_K(z) = A*sig_k_K(z+1);
    end
    
    x_k_k(1) = x_k_K(1);
    sig_k_k(1) = sig_k_K(1);
    
    
    B = sum(sig_k_K(2:end))+alpha^2*sum(sig_k_K(1:end-1))+sum((x_k_K(2:end)).^2)+alpha^2*sum((x_k_K(1:end-1)).^2)-2*alpha*sum(x_k_K(2:end).*x_k_K(1:end-1)+sig_k_k1_K(1:end-1));
    sig(i+1) = B/N;

end

x_SS = x_k_K(2:N+1);

end