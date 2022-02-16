% Load in and read ST resistance values
function ST = plot_ST_resistance(num_term)
name = 'STnetwork_';

% num_term = 45;
vessels = zeros(num_term,1);
alpha   = zeros(num_term,1);
beta    = zeros(num_term,1);
% Z       = zeros(num_term,1);
Z = [];
for i=1:num_term
   term=i-1;
   if term<10
       filename = strcat(name,num2str(term));
   else
      filename = strcat(name,num2str(term)); 
   end
   temp = load(filename);
   vessels(i) = temp(1);
   alpha(i)   = temp(2);
   beta(i)    = temp(3);
   if isempty(Z)
       num_t = length(temp(4:end));
   end
       Z(end+1,1:num_t) =temp(4:end); 
end
disp(sum(vessels))
figure(919); hold on;
plot(Z')

ST = {vessels, alpha, beta, Z};
end