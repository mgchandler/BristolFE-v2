% Extract the desired portion of the data (750 x 1500:3000)
data_subset = I2(1:50, 2900:3000);

% Flatten the data into a single column vector
data_flat = data_subset(:);

% Calculate the unique values and their counts
unique_values = unique(data_flat);
counts = histcounts(data_flat, [unique_values; max(unique_values)+1]);

% Initialize variables for counting repeated numbers and their occurrences
total_repeats = 0;
total_occurrences = 0;

% Iterate through unique values to count repeats
for i = 1:numel(unique_values)
    if counts(i) > 1 % If the value is repeated
        total_repeats = total_repeats + 1; % Count the repeated value
        total_occurrences = total_occurrences + counts(i); % Count occurrences
    end
end

% Calculate average repetitions
average_repetitions = total_occurrences / total_repeats;

disp(['Average number of repetitions over the specified area: ', num2str(average_repetitions)]);
