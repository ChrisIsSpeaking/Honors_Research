%% Part A
numsamps = 10000;
rho = 0.5;
cov_mat = [1 rho; rho 1];
R = chol(cov_mat);
C = zeros(numsamps,2);
for c = 1:numsamps
    C(c,:) = (R'*randn(2,1))';
end
cov_matrix = cov(C);
mean_vec = mean(C);

%% Part B
numsamps = 100000;
y0samples = zeros(numsamps,1);
y1samples = zeros(numsamps,1);
%intiializing the first sample
y0samples(1) = 0; y1samples (1) = 0;
%amples(1,:) = 0;
%from derivation, the mean of either depends on one half the mean of the
%other while the standard deviation is fixed sigma squared = 0.75
sigma = sqrt(0.75);
for k = 1:numsamps-1
    y0samples(k+1) = 0.5*y1samples(k) + sigma*randn(1);
    y1samples(k+1) = 0.5*y0samples(k+1) + sigma*randn(1);
end
C = [y0samples, y1samples];
cov_matrix = cov(C);
mean_vec = mean(C);

%% Part C
numsamps= 100000;
%Initializing the first sample
y0 = [0,0]';
C = [numsamps,2]; C(1,:) = y0';
cov_mat_inv = [4/3 -2/3; -2/3 4/3];

%Distribution we want to estimate/one that we know we are proportional to
probfunc = @(y) 1/(2*pi*sqrt(0.75)).*exp(-0.5*(y'*cov_mat_inv*y));

for k = 1:numsamps-1
    %Since the covariance matrix is the identity matrix, the sigma is 1
    yk = C(k,:)';
    %Drawing 
    yk1 = yk + randn(2,1);
    accprob = probfunc(yk1)/probfunc(yk);
    alpha = min([1,accprob]);
    u = rand;
    if u <= alpha
        C(k+1,:)= yk1';
    else
        C(k+1,:) = yk';
    end
end
cov_matrix = cov(C);
mean_vec = mean(C);