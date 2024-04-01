function [R_data,L_data] = get_dim_data
%% Do cross validation on the bandwidth parameter

%% Step one: Maximize the bandwidth for the whole set of density
% clear; clc; close all;
conn_name = 'Connectivity_';
dim_name  = 'Dimensions_';
theta1 = [22 25 26  26   26   27  28  29.9 30  30  30  31  31  32  33   33    34  34  35  35  35  36  36  37  43.6];
theta2 = [5   6 4.7 4.8  5.1  5.8 6   4.6 5.7 6.5  8  5.6 6.1 4.1 4.2  5.1   3.3 3.4 3.6 4.8 6.8  4  4.1 3.9 7.62];
num_par = length(theta1);
%%%%%%%%%% LOAD SCALING IF YOU WANT TO CONVERT TO PIXELS
load('ScalingValues.mat');
%% 1 or 0, if you want to scale;
scale_flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vessels to compare
name1 = 'M1P4_lower';
name2 = '_smooth';
name3 = '_Data.mat';
%% BRAND NEW 34 VESSEL SYSTEM
% vessels_included{1}  =  [0 1   31 151 71 32 90    89 109 91  179 219 182 224 220 225 244 231 226 193 197 196 110 111 112      135 186 216 243 237 53 66 63 94] + 2;
% vessels_included{2}  =  [0 1   97 2   155 98 193  156 194 247 282 330 283 338 333 339 387 347 340 298 299 301 195 197 196     226 296 327 382 365 130 147 142 275] + 2;
% vessels_included{3}  =  [0 1   107 2  108 228 109 195 110 168 283 327 284 328 429 390 329 398 391 296 300 299 165 166 167     148 289 308 400 419 275 253 264 190] + 2;
% vessels_included{4}  =  [0 170 171 358 172 299 208 173 209 271 1   50   2  57  53 58  100   59 93  44  48  47  268 269 270    263 9   29   89 76  335 314 325 292] + 2;
% vessels_included{5}  =  [0 169 170 361 171 298 205 172 206 268 1   50   2  52 164 53  94    54 89  17  21  19  265 267 266    216 13  36   78  72 352 332 343 291] + 2;
% vessels_included{6}  =  [0 181 182 373 183 312 221 184 251 222 1   54   2  57 177 58  103   59 96  17  21  19  308 310 309    279  9  41   85  77 348 326 337 224] + 2;
% vessels_included{7}  =  [0 134 198 135 242 199 270 243 298 271 1   38   2  44  39 45  77    53 46  6   7    9  299 300 301    345 35  29   61  65 237 227 234 274] + 2;
% vessels_included{8}  =  [0 1   100 2   158 101 198 159 230 199 295 339 296 345 340 346 387  354 347 311 312 314 231 232 233   261 309 336 384 381 155 146 141 204] + 2;
% vessels_included{9}  =  [0 181 182 388 183 321 224 184 256 225 1   54  2   60  56  61  110  69  62  48  52  51  257 258 259   291 9   30  109 87  365 336 331 227] + 2;
% vessels_included{10} =  [0 183 184 393 185 326 186 285 187 253 1   56  2   58  178 59  112  60  103 17  21  20  250 252 251   227 9   39  90  78  366 341 352 255] + 2;
% vessels_included{11} =  [0 1   98  2   152 99  188 153 222 189 281 323 282 329 324 330 369 340 331 297 298  300 223 224 225   244 288 320 366 365 149 137 146 194] + 2;
% vessels_included{12} =  [0 183 286 184 348 287 349 448 381 350 1   52  2   59  55  60  109 61  100 46  50   49  382 383 384   420 8   32  71  76  294 313 324 353] + 2;
% vessels_included{13}  = [0 127 128 259 129 220 153 130 154 197 1   38  2   39  122 40  71  41  64  32  36   35  155 157 156   188 8   23  63  57  243 231 236 199] + 2;
% vessels_included{14}  = [0 166 261 167 262 382 263 349 295 264 1   48  2   55  51  123 56  131 124 42  46   45  346 348 347   319 8   33  161 145 415 397 392 288] + 2;
% vessels_included{15}  = [0 175 176 375 177 312 217 178 218 283 1   50  2   56  52  130 57  131 168 17  21   20  279 281 280   262 9   36  165 146 365 338 349 306] + 2;
% vessels_included{16}  = [0 198 199 418 200 347 242 201 243 315 1   56  2   63  59  144 64  145 186 50  51   53  310 312 311   262 10  38  147 152 355 373 384 338] + 2;
% vessels_included{17}  = [0 1   117 2   118 261 120 217 121 185 332 333 384 385 387 391 442 402 393 379 383  381 182 184 183   151 342 365 405 413 270 286 297 210] + 2;
% vessels_included{18}  = [0 193 194 388 195 330 232 196 233 297 1   56  2   63  59  64  110 65 103  9   13   11  294 296 295   244 49  33  77  82  350 358 369 323] + 2;
% vessels_included{19}  = [0 1   105 2   169 106 170 268 206 171 305 356 307 357 475 358 404 359 395 339 343  341 207 208 209   240 352 318 391 366 115 131 126 174] + 2;
% vessels_included{20}  = [0 184 185 382 186 313 187 277 188 246 1   56  2   64  59  135 65  145 136 8   12   10  243 245 244   217 48  26  180 149 320 339 334 249] + 2;
% vessels_included{21}  = [0 149 150 330 203 151 243 204 271 244 1   46  2   52  48  113 53  121 114 15  19   17  272 273 274   304 9   27  146 141 184 168 163 246] + 2;
% vessels_included{22}  = [0 1  102  2   160 103 200 161 232 201 297 343 298 349 344 350 385 358 351 302 306  305 233 234 235   259 342 329 382 371 156 147 144 206] + 2;
% vessels_included{23} = [0 149 150  305 151 260 179 152 207 180 1   42  2   48  44  49  86  85  80  13  17   16  256 258 257   242 10  30  77  73  300 282 293 201] + 2;
% vessels_included{24} = [0 1   94   2   150 95  188 151 222 189 285 331 286 337 332  338 373 346 339 301 302 304 223 224 225   253 299 323 370 357 147 123 120 194] + 2;
% vessels_included{25} = [0 1   38   2   62  39  78  63  90  79  120 138 121 140 139  141 162 147 142 128 130 129 91  92  93    119 125 137 148 156 50  59  58  81] + 2;

%% only 32 vessels 

vessels_included{1}  =  [0 1   31 151 71 32 90    89 109 91  179 219 182 224 220 225 244 231 226 193 196 110 111 112      135 186 216 243 237 66 63 94] + 2;
vessels_included{2}  =  [0 1   97 2   155 98 193  156 194 247 282 330 283 338 333 339 387 347 340 298 301 195 197 196     226 296 327 382 365 147 142 275] + 2;
vessels_included{3}  =  [0 1   107 2  108 228 109 195 110 168 283 327 284 328 429 390 329 398 391 296 299 165 166 167     148 289 308 400 419 253 264 190] + 2;
vessels_included{4}  =  [0 170 171 358 172 299 208 173 209 271 1   50   2  57  53 58  100   59 93  44  47  268 269 270    263 9   29   89 76  314 325 292] + 2;
vessels_included{5}  =  [0 169 170 361 171 298 205 172 206 268 1   50   2  52 164 53  94    54 89  17  19  265 267 266    216 13  36   78  72 332 343 291] + 2;
vessels_included{6}  =  [0 181 182 373 183 312 221 184 251 222 1   54   2  57 177 58  103   59 96  17  19  308 310 309    279  9  41   85  77 326 337 224] + 2;
vessels_included{7}  =  [0 134 198 135 242 199 270 243 298 271 1   38   2  44  39 45  77    53 46  6    9  299 300 301    345 35  29   61  65 227 234 274] + 2;
vessels_included{8}  =  [0 1   100 2   158 101 198 159 230 199 295 339 296 345 340 346 387  354 347 311 314 231 232 233   261 309 336 384 381 146 141 204] + 2;
vessels_included{9}  =  [0 181 182 388 183 321 224 184 256 225 1   54  2   60  56  61  110  69  62  48  51  257 258 259   291 9   30  109 87  336 331 227] + 2;
vessels_included{10} =  [0 183 184 393 185 326 186 285 187 253 1   56  2   58  178 59  112  60  103 17  20  250 252 251   227 9   39  90  78  341 352 255] + 2;
vessels_included{11} =  [0 1   98  2   152 99  188 153 222 189 281 323 282 329 324 330 369 340 331 297  300 223 224 225   244 288 320 366 365 137 146 194] + 2;
vessels_included{12} =  [0 183 286 184 348 287 349 448 381 350 1   52  2   59  55  60  109 61  100 46   49  382 383 384   420 8   32  71  76  313 324 353] + 2;
vessels_included{13}  = [0 127 128 259 129 220 153 130 154 197 1   38  2   39  122 40  71  41  64  32   35  155 157 156   188 8   23  63  57  231 236 199] + 2;
vessels_included{14}  = [0 166 261 167 262 382 263 349 295 264 1   48  2   55  51  123 56  131 124 42   45  346 348 347   319 8   33  161 145 397 392 288] + 2;
vessels_included{15}  = [0 175 176 375 177 312 217 178 218 283 1   50  2   56  52  130 57  131 168 17   20  279 281 280   262 9   36  165 146 338 349 306] + 2;
vessels_included{16}  = [0 198 199 418 200 347 242 201 243 315 1   56  2   63  59  144 64  145 186 50   53  310 312 311   262 10  38  147 152 373 384 338] + 2;
vessels_included{17}  = [0 1   117 2   118 261 120 217 121 185 332 333 384 385 387 391 442 402 393 379  381 182 184 183   151 342 365 405 413 286 297 210] + 2;
vessels_included{18}  = [0 193 194 388 195 330 232 196 233 297 1   56  2   63  59  64  110 65 103  9    11  294 296 295   244 49  33  77  82  358 369 323] + 2;
vessels_included{19}  = [0 1   105 2   169 106 170 268 206 171 305 356 307 357 475 358 404 359 395 339  341 207 208 209   240 352 318 391 366 131 126 174] + 2;
vessels_included{20}  = [0 184 185 382 186 313 187 277 188 246 1   56  2   64  59  135 65  145 136 8    10  243 245 244   217 48  26  180 149 339 334 249] + 2;
vessels_included{21}  = [0 149 150 330 203 151 243 204 271 244 1   46  2   52  48  113 53  121 114 15   17  272 273 274   304 9   27  146 141 168 163 246] + 2;
vessels_included{22}  = [0 1  102  2   160 103 200 161 232 201 297 343 298 349 344 350 385 358 351 302  305 233 234 235   259 342 329 382 371 147 144 206] + 2;
vessels_included{23} = [0 149 150  305 151 260 179 152 207 180 1   42  2   48  44  49  86  85  80  13   16  256 258 257   242 10  30  77  73  282 293 201] + 2;
vessels_included{24} = [0 1   94   2   150 95  188 151 222 189 285 331 286 337 332  338 373 346 339 301 304 223 224 225   253 299 323 370 357 123 120 194] + 2;
vessels_included{25} = [0 1   38   2   62  39  78  63  90  79  120 138 121 140 139  141 162 147 142 128 129 91  92  93    119 125 137 148 156  59  58  81] + 2;




%%


num_vessels = length(vessels_included{1});
plottingcolors = linspecer(num_vessels);
legendnames = {};
for i=1:num_vessels
    legendnames{end+1} = strcat('V',num2str(i));
end
%%%%%%%%%% LOAD SCALING IF YOU WANT TO CONVERT TO PIXELS
load('ScalingValues.mat');
%% 1 or 0, if you want to scale;
scale_flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_vessels = length(vessels_included{1});

% This is a vector identifying the parameter values that give crap output
% bad = [7 10 20];
bad = [];
%%
L = zeros(num_vessels,num_par - length(bad));
R = zeros(num_vessels,num_par - length(bad));
% Rmin = zeros(num_vessels,num_par - length(bad));
% Rmed = zeros(num_vessels,num_par - length(bad));
% Rmax = zeros(num_vessels,num_par - length(bad));
R_vessel = cell(num_vessels,num_par - length(bad));
for i=1:num_par - length(bad)
    if ~any(i == bad)
        if scale_flag == 1
            scale = scaling(i);
        else
            scale = 1;
        end
        network = load(strcat(name1,num2str(theta1(i)),name2,num2str(theta2(i)),name3));
        connectivity = network.connectivity;
        dets = network.vessel_details;
        for j=1:num_vessels
            if i==1 && j==4
            elseif i==23 && j==18
            else
            %% New Line here to do 25,50, and 80% of the vessels
            length_vessel = length(dets{vessels_included{i}(j),2});
            c_20 = round([0.4 0.6].*length_vessel);
            c_50 = round([0.25 0.75].*length_vessel);
            f_50 = [1 round(0.5.*length_vessel)];
            c_80 = round([0.1 0.9].*length_vessel);
%             R(j,i) = mean(dets{vessels_included{i}(j),2}(c_20(1):c_20(2),4)./scale);
%             R(j,i) = mean(dets{vessels_included{i}(j),2}(c_50(1):c_50(2),4)./scale);
            R(j,i) = mean(dets{vessels_included{i}(j),2}(f_50(1):f_50(2),4)./scale);
%             R(j,i) = mean(dets{vessels_included{i}(j),2}(c_80(1):c_80(2),4)./scale);
%             R(j,i) = dets{vessels_included{i}(j),4}./scale;
            L(j,i) = dets{vessels_included{i}(j),3}./scale;
            %             R(j,i) = min(dets{vessels_included{i}(j),2}(:,4));
            %             R(j,i) = median(dets{vessels_included{i}(j),2}(:,4));
            %             R(j,i) = max(dets{vessels_included{i}(j),2}(:,4));
            R_vessel{j,i} = dets{vessels_included{i}(j),2}(:,4);
            end
        end
    end
end
%% Try and analyze the data via standardization
L_bin = zeros(num_vessels*(num_par),1);
R_bin = zeros(num_vessels*(num_par),1);
L_bin2 = L_bin;
R_bin2 = R_bin;
LRR_bin = zeros(num_vessels*(num_par),1);
L_data = [];
R_data = [];
R_bin_vessel = [];
%%

for i=1:num_vessels%[1 2 9 3 4 10 11 5 6 12 13 7 8 14 15]
    
    L_curr = L(i,:);
    R_curr = R(i,:);

    R_curr_vessel = [];
    for j=1:num_par
            %% NEW MJC 9/7/18
%             figure(100+i); hold on;
%             plot(R_vessel{i,j},'o','Color',plottingcolors(j,:));
            %%
        R_vessel_curr = R_vessel{i,j};
        size_VC = length(R_vessel_curr);
        R_curr_vessel(end+1:end+size_VC) = R_vessel_curr;
    end
    legendnames = {};
    for leg=1:25
        legendnames{end+1} = strcat('seg',num2str(leg));
    end
    legend(legendnames)
    L_curr = L_curr(L_curr~=0);
    R_curr = R_curr(R_curr~=0);
    mu_L = mean(L_curr);
    mu_R = mean(R_curr);
    mu_R_vessel = mean(R_curr_vessel);
    sigma_L = std(L_curr);
    sigma_R = std(R_curr);
    sigma_R_vessel = std(R_curr_vessel);
    med_L = median(L_curr);
    med_R = median(R_curr);
    med_R_vessel = median(R_curr_vessel);
    mad_L = mad(L_curr);
    mad_R = mad(R_curr);
    mad_R_vessel = mad(R_curr_vessel);
    L_std = (L_curr-mu_L)./sigma_L;
    R_std = (R_curr-mu_R)./sigma_R;
    R_std_vessel = (R_curr_vessel-mu_R_vessel)./sigma_R_vessel;
    
    L_std2 = (L_curr-med_L)./(1.4826*mad_L);
    R_std2 = (R_curr-med_R)./(1.4826*mad_R);
    
    L_data(end+1,1)   = mu_L;
    L_data(end,2)   = sigma_L./mu_L;
    L_data(end,3)   = sigma_L;
    
    
    R_data(end+1,1)   = mu_R;
    R_data(end,2)   = sigma_R./mu_R;
    R_data(end,3)   = sigma_R;
    
    
    startID = (length(R_curr))*(i-1)+1;
    endID   = (length(R_curr))*(i);
    L_bin(startID:endID) = L_std';
    R_bin(startID:endID) = R_std';
    L_bin2(startID:endID) = L_std2';
    R_bin2(startID:endID) = R_std2';
    LRR_bin(startID:endID) = L_curr./R_curr;
    
    % For the whole vessel
    size_rves = length(R_std_vessel);
    R_bin_vessel(end+1:end+size_rves) = R_std_vessel;
    
end

end