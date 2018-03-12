function clear_hi_aggr(current_struct)
labels = fieldnames(current_struct);

% Sum total power of the workers
if ismember('workers',labels)
    for i=1:length(current_struct.workers)
        current_struct.workers(i).Pm_opt = [];
        current_struct.workers(i).e_opt = [];
        current_struct.workers(i).U = [];
    end
end

q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
if ~isempty(q_labels) % if there are Queens, go on recursively
    for i=1:length(q_labels)
    clear_hi_aggr(current_struct.(q_labels{i}));
    end
end


end