%% Compute the mean values per id. 
% They are all the same so in fact is just to get a unique value.
id = unique(semdata(:,1));
n_id = numel(id);
n_features = size(semdata,2);
data = NaN(n_id, n_features);
for i=1:numel(id)
    rows = find(semdata(:,1)==id(i));
    mean_data = mean(semdata(rows,:));
    data(i,:) = mean_data;
end

%% Update data removing features without models and renaming ids
CHANGE_ID = [9 112; 10 113; 16 114; 32 115];
for i=1:size(CHANGE_ID,1)
    remove_row = find(data(:,1)==CHANGE_ID(i,1));
    if ~isempty(remove_row)
        data(remove_row,:) = [];
    end
    update_row = find(data(:,1)==CHANGE_ID(i,2));
    if ~isempty(update_row)
        data(update_row,1) = CHANGE_ID(i,1);
    end
end
n_id = size(data,1);

%% Sort by first row (id)
data = sortrows(data,1);

%% Swap Features
SWAP_ID = [3 4];
ordered_labels = semdata_labels;
for i=1:size(SWAP_ID,1)
    data(:,[SWAP_ID(i,1) SWAP_ID(i,2)]) = data(:,[SWAP_ID(i,2) SWAP_ID(i,1)]);
    aux_label = ordered_labels{SWAP_ID(i,1)};
    ordered_labels{SWAP_ID(i,1)} = ordered_labels{SWAP_ID(i,2)};
    ordered_labels{SWAP_ID(i,2)} = aux_label;
end

%% Remove undesired features: zero, zero 2, zero 3 and zero 4
BAD_FEATURE_IDX = [16, 15, 14, 13, 5, 1];
n_bad_features = numel(BAD_FEATURE_IDX);
n_good_features = n_features - n_bad_features;
for i=1:n_bad_features
    data(:, BAD_FEATURE_IDX(i)) = [];
end

%% Copy desired labels
feature_label = cell(n_good_features, 1);
last_i = 1;
for i=1:n_features
    if ~ismember(i, BAD_FEATURE_IDX)
        feature_label{last_i} = ordered_labels{i};
        last_i = last_i + 1;
    end
end

%% Compute data format
INT_TYPE = 0;
FLOAT_TYPE = 1;
OPTIONS_TYPE = 2;
OPTIONS_LABEL = cell(2,1);
OPTIONS_LABEL{1} = 'Female';
OPTIONS_LABEL{2} = 'Male';
feature_type = NaN(n_good_features,1); % 0-int, 1-float, 2-options
data_format = '';
for i=1:n_good_features
    digit_decimal = regexp(num2str(max(data(:,i))),'\.','split');
    if numel(digit_decimal)>1 % float
        data_format = strcat(data_format, '%10.4f');
        feature_type(i) = FLOAT_TYPE;
    else % int THIS CAN BE IMPROVED
        data_format = strcat(data_format, '%5.0f');
        if numel(unique(data(:,i)))==2
            feature_type(i) = OPTIONS_TYPE;
        else
            feature_type(i) = INT_TYPE;
        end
    end
    if i==n_good_features
        data_format = strcat(data_format, '\n');
    else
%         data_format = strcat(data_format, {' '});
    end
end

%% MINIMUM AND MAXIMUM VALUES BY HAND!!!
min_values = zeros(n_good_features);
max_values = 10.^ceil(log10(ceil(max(data))));

%% Write results in file
FILE_NAME = 'body_data_features.txt';
fileID = fopen(FILE_NAME, 'w');
fprintf(fileID,'%d\n',n_good_features);
fprintf(fileID,'%d\n',n_id);
for i=1:n_good_features
    fprintf(fileID,['%' num2str(numel(feature_label{i})) 's\n'], ...
        feature_label{i});
end
for i=1:n_good_features
    fprintf(fileID,'%5.0f\n', min_values(i));
end
for i=1:n_good_features
    fprintf(fileID,'%5.0f\n', max_values(i));
end
% for i=1:n_good_features
%     fprintf(fileID,'%d\n',feature_type(i));
%     if feature_type(i) == 2 % Options
%         num_options = numel(unique(data(:,i)));
%         fprintf(fileID,'%d\n',num_options);
%         for j=1:num_options
%             fprintf(fileID,['%' num2str(length(OPTIONS_LABEL{j})) 's\n'], ...
%                 OPTIONS_LABEL{j});
%         end
%     end
% end
fprintf(fileID, data_format, data');
fclose(fileID);
