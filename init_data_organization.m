% loads, organizes, amd saves large data into structs for data analysis and plotting 
data_folder = [pwd '\data\'];
data_types = {'cue_traces'; 'ITI_traces'};
treatments = {'Saline'; 'MDMA'};
mice_per_group = 8;
num_days = 5;

for run_data_types = 1:2
    curr_data_type = data_types{run_data_types};
    current_struct = struct();
    for treatment = 1:2
        current_struct.(treatments{treatment}) = cell(mice_per_group, num_days);
    end
    for treatment = 1:2
        for run_mice = 1:mice_per_group
            curr_mouse_name = [treatments{treatment} '_mouse_' num2str(run_mice)];
            for run_days = 1:num_days
                filepath = [data_folder '\' curr_data_type '\' curr_mouse_name '_day_' num2str(run_days) '.mat'];
                loaded = load(filepath, 'day_traces_mat');
                current_struct.(treatments{treatment}){run_mice, run_days} = double(loaded.day_traces_mat);
            end
        end
    end
    eval([curr_data_type ' = current_struct;']);
    save([data_folder curr_data_type '.mat'], curr_data_type);
end