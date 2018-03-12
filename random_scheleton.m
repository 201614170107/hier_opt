function scheleton = random_scheleton(n_l,n_b,Pnoms,Vnoms)
% Create a random scheleton with 
% Inputs: -n_l: maximum number of levels. Assume uniform random
%         -n_b: number of branches per level. This is fixed
%         -Pnoms,Vnoms, vectors of nominal powers and nominal voltages,
%         these have length n_l. 


N_l = randperm(n_l,1);
scheleton = rec_scheleton([],N_l,n_b,Pnoms,Vnoms);
disp_struct(scheleton)
end

function scheleton = rec_scheleton(scheleton,n_l,n_b,Pnoms,Vnoms)
scheleton.V_ref = Vnoms(1);
scheleton.P_ref = Pnoms(1);
if n_l>1
    for i=1:randperm(n_b,1)
        scheleton.(strcat('q',num2str(i)))  = [];
        scheleton.(strcat('q',num2str(i))) =  rec_scheleton(scheleton.(strcat('q',num2str(i))),n_l-1,n_b,Pnoms(2:end),Vnoms(2:end));
    end
end
end

function disp_struct(x,indent)
    % Custom function for displaying a structure. Use recursion.
    % Input: -x: a nested structure
    %        -indent: do not specify this input. For recursion
    %        purposes only.
    if nargin<2
        indent='';
    end
    if isstruct(x)
        labels = fieldnames(x);
        for i=1:length(labels)
            if isstruct(x.(labels{i}))
                fprintf('\n%s%s',[indent,labels{i}])
                disp_struct(x.(labels{i}),strcat(indent,'---'))
            else
                fprintf('\n%s%s%s',[indent,labels{i},':'])
                fprintf('\t[%ix%i %s] ',[size(x.(labels{i}),1),size((x.(labels{i})),2),class(x.(labels{i}))])
            end
        end
    end
    fprintf('\n')
end