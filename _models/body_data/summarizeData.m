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

%% Remove undesired features: zero, zero 2, zero 3 and zero 4
BAD_FEATURE_IDX = [16, 15, 14, 13, 5];
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
        feature_label{last_i} = semdata_labels{i};
        last_i = last_i + 1;
    end
end

%% Write results in file
FILE_NAME = 'body_data_features.txt';
fileID = fopen(FILE_NAME, 'w');
fprintf(fileID,'%d\n',n_good_features);
fprintf(fileID,'%d\n',n_id);
for i=1:n_good_features
    fprintf(fileID,['%' num2str(numel(feature_label{i})) 's\n'], ...
        feature_label{i});
end
data_format = '';
for i=1:n_good_features
    digit_decimal = regexp(num2str(max(data(:,i))),'\.','split');
    if numel(digit_decimal)>1 % float
        data_format = strcat(data_format, '%10.4f');
    else % int
        data_format = strcat(data_format, '%5.0f');
    end
    if i==n_good_features
        data_format = strcat(data_format, '\n');
    else
%         data_format = strcat(data_format, {' '});
    end
end
fprintf(fileID, data_format, data');
fclose(fileID);
