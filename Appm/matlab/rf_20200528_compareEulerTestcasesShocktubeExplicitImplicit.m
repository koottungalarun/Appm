clear
clc

impl = importdata('testcase_shocktube_implicit.csv');
expl = importdata('testcase_shocktube_explicit.csv');

colheaders = [...
    """Points:2""", 
    """Neutrals number density""",
    """Neutrals pressure""",
    """Neutrals velocity:2"""    
    ];

for i = 1 : length(colheaders)
    idx(i) = find(strcmp(impl.colheaders, colheaders{i}));
end

implData = impl.data(:,idx);
explData = expl.data(:,idx);

plot(implData(:,1), implData(:,2:end))
hold on
set(gca, 'ColorOrderIndex', 1)
plot(explData(:,1), explData(:,2:end), '--')
hold off
grid on
shg




