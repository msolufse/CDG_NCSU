%% Filter the signal
function out = filter_signal(input,cutoff,num_pts)

temp = fft(input);
temp(2+cutoff:end-cutoff) = 0;
out = ifft(temp);

%% Interpolate to more points if need be
out = interp1(1:length(out),out,linspace(1,length(out),num_pts));

end