function dispstruct(x,indent)

if nargin<2
    indent='';
end

if isstruct(x)
    labels = fieldnames(x);
    
    for i=1:length(labels)
        if isstruct(x.(labels{i}))
            fprintf('\n%s%s',[indent,labels{i}])
            dispstruct(x.(labels{i}),strcat(indent,'--'))
        else
            fprintf('\n%s%s',[indent,labels{i}])
            fprintf('\t[%s] ',class(x.(labels{i})))
        end
    end
end

fprintf('\n')