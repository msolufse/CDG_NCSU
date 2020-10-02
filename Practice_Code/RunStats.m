%% Running calculate_statistics code and saving for all 8 residuals 
function Stats = RunStats(data)
for i = 1:1
    for j = 1:8
        Stats{i,j} = calculate_statisticsWORKING(i,j,data);
    end 
end 

end