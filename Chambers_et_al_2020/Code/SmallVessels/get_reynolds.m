%% Plot the reynolds number as a function of radius.
r_root = 0.02;
aa  = 0.8837;
bb = 0.6666;
r_min = 0.005;

% Determine the number of plots generated
alpha_r = r_root.*aa.^(1:100);
beta_r = r_root.*bb.^(1:100);
n = find(alpha_r<r_min,1);% alpha generations
m = find(beta_r<r_min,1); % beta generations

Palpha = zeros(1024,n);
Qalpha = zeros(1024,n);
Time = zeros(1024,n);
Pbeta = zeros(1024,m);
Qbeta = zeros(1024,m);


%%
ra_d = [];
rb_d = [];

for i = 1:n
    s = num2str(i-1);
    s = strcat('p',s,'_alpha.2d');
    data = load(s);
    Palpha(:,i) = data(:,1);
    Qalpha(:,i) = data(:,2)*10/T; % Dimensionalize correctly!
    Time_a(:,i) = t;
    ra_d = [ra_d r_root*(aa^(i-1))];
end

for i = 1:m
    s = num2str(i-1);
    s = strcat('p',s,'_beta.2d');
    data = load(s);
    Pbeta(:,i) = data(:,1);
    Qbeta(:,i) = data(:,2)*10/T;  % Dimensionalize correctly!
    Time_b(:,i) = t;
    rb_d = [rb_d r_root*(bb^(i-1))];
end


La = lrr(ra_d);
Lb = lrr(rb_d); % vessel's length

Lac = cumsum(La);
Lbc = cumsum(Lb); % cumulative pathe length