clear
close
sep = filesep;
%% load data
dataset_name = {'box1','box2','walking1','walking2','walking_ts'};
for i = 1:length(dataset_name)
data_set = dataset_name{i} ;
data_file = strcat('data',sep,data_set,'.txt');
data_file_output = strcat(data_set,'.mat');
% A struct of TD events with format:
% TD(:,1) = event timestamps in microseconds
% TD(:,2) = pixel X locations
% TD(:,3) = pixel Y locations
% TD(:,4) = event polarity
% TD(:,5) = left or right label left is 0 right is 1
% TD(:,6) = disparity
TD = load(data_file);
% show the input stereo data
ShowTD_stereo(TD);
%% stereo depth estimation
tic
stereo_TD = stereo_matching_bp(TD,1);
time_cost = toc;

%% evaluate and show
ShowTD_depth(stereo_TD);
end

% the histogram of the estimated disparity
figure
hist(stereo_TD.p,1:50);
ylabel('Number of events'),xlabel('Disparity of EMP'),axis ([0,50,0,5000]);

% the detection rate
detection_rate = length(stereo_TD.ts)/length(TD(TD(:,5)==0));

% the matching rate
stereo_TD_error = stereo_TD.p - stereo_TD.disparity_gt;
match_rate = sum(abs(stereo_TD_error)<=1)/length(stereo_TD.p);

% distance matching rate
error_l = 1:-0.01:0;
stereo_TD.depth = 250 * 0.12 ./stereo_TD.p(:);
stereo_TD.depth_gt = 250 * 0.12 ./stereo_TD.disparity_gt(:);
stereo_TD_error = abs(stereo_TD.depth(stereo_TD.depth~=inf) - stereo_TD.depth_gt(stereo_TD.depth~=inf));
stereo_TD_error_norm = stereo_TD_error./stereo_TD.depth_gt(stereo_TD.depth~=inf);

% depth accuracy and error tolarance
for index_error = 1:length(error_l)
    match_rate_d(index_error) = sum(stereo_TD_error_norm <= error_l(index_error))/length(stereo_TD.depth(stereo_TD.depth~=inf));
end

figure(4)
plot(error_l,match_rate_d,'LineWidth',2);
xlabel('Error tolerance','FontSize',14)
ylabel('Depth accuracy','FontSize',14)
legend({'EMP'},'FontSize',14,'Location','southeast')


