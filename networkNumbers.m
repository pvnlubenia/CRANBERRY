% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    networkNumbers                                                         %
%                                                                           %
%                                                                           %
% OUTPUT: Returns the values of the following network numbers of a chemical %
%    reaction network:                                                      %
%       - Species (m)                                                       %
%       - Complexes (n)                                                     %
%       - Reactant complexes (n_r)                                          %
%       - Reversible reactions (r_rev)                                      %
%       - Irreversible reactions (r_irrev)                                  %
%       - Reactions (r) = r_irrev + 2 r_rev                                 %
%       - Linkage classes (l)                                               %
%       - Strong linkage classes (sl)                                       %
%       - Terminal linkage classes (t)                                      %
%       - Rank (s)                                                          %
%       - Reactant rank (q)                                                 %
%       - Deficiency (delta) = n - l - s                                    %
%       - Reactant deficiency (delta_p) = n_r - q                           %
%    The output variable 'model' allows the user to view the complete       %
%       network with all the species listed in the 'species' field of       %
%       the structure 'model'.                                              %
%                                                                           %
% INPUT: model: a structure, representing the CRN (see README.txt for       %
%    details on how to fill out the structure)                              %
%                                                                           %
% References:                                                               %
%    [1] Arceo C, Jose E, Lao A, Mendoza E (2017) Reactant subspaces and    %
%           kinetics of chemical reaction networks. J Math Chem             %
%           56(5):395–422. https://doi.org/10.1007/s10910-017-0809-x        %
%    [2] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical        %
%           reaction network theory. Bioinform 25(21):2853–2854.            %
%           https://doi.org/10.1093/bioinformatics/btp513                   %
%                                                                           %
% Created: 3 June 2022                                                      %
% Last Modified: 29 August 2024                                             %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [model] = networkNumbers(model)
    
%
% Create a list of all species indicated in the reactions
%

% Initialize list of species
model.species = { };

% Get all species from reactants
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).reactant)
        model.species{end+1} = model.reaction(i).reactant(j).species;
    end
end

% Get species from products
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).product)
        model.species{end+1} = model.reaction(i).product(j).species;
    end
end

% Get only unique species
model.species = unique(model.species);



%
% Species
%

% Count the number of species
m = numel(model.species);



%
% Complexes
%

% Initialize the matrix of reactant complexes
reactant_complexes = [ ];

% Initialize the matrix of product complexes
product_complexes = [ ];

% Initialize the stoichiometric matrix
N = [ ];

% For each reaction in the model
for i = 1:numel(model.reaction)
    
    % Initialize the vector for the reaction's reactant complex
    reactant_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the reactant complex
    for j = 1:numel(model.reaction(i).reactant)
        reactant_complexes(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
    end
    
    % Initialize the vector for the reaction's product complex
    product_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the product complex
    for j = 1:numel(model.reaction(i).product)
        product_complexes(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
    end
    
    % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
    N(:, end+1) = product_complexes(:, end) - reactant_complexes(:, end); % N %
    
    % If the reaction is reversible
    if model.reaction(i).reversible
        
        % Insert a new vector for the reactant complex: make it same as the product complex
        reactant_complexes(:, end+1) = product_complexes(:, end);
        
        % Insert a new vector for the product complex: make it the same as the reactant complex
        product_complexes(:, end+1) = reactant_complexes(:, end-1);
        
        % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
        N(:, end+1) = -N(:, end);
    end
end

% Get just the unique complexes
% index(i) is the index in all_complex of the reactant complex in reaction i
[all_complex, ~, index] = unique([reactant_complexes product_complexes]', 'rows');

% Construct the matrix of complexes
all_complex = all_complex';

% Count the number of complexes
n = size(all_complex, 2);



%
% Reactant complexes
%

% Get just the unique reactant complexes
reactant_complexes_unique = unique([reactant_complexes]', 'rows')';

% Count the number of unique reactant complexes
n_r = size(reactant_complexes_unique, 2);



%
% Reversible, irreversible, and total reactions
%

% Count the number of reversible, irreversible, and total reactions
r_rev = 0;
r_irrev = 0;
for i = 1:numel(model.reaction)
    if (model.reaction(i).reversible == 1)
        r_rev = r_rev + 1;
     else
        r_irrev = r_irrev + 1;
    end
end
r = r_irrev + 2*r_rev;



%
% Linkage classes
%

% Initialize a matrix (complexes x complexes) for the reacts_to relation
% This is for testing reversibility of the network
reacts_to = false(n, n);

% Initialize matrix (complexes x total reactions) for the reacts_in relation
% This is the incidence matrix I_a
reacts_in = zeros(n, r);

% Fill out the entries of the matrices
for i = 1:r
    
    % reacts_to(i, j) = true iff there is a reaction r: y_i -> y_j
    reacts_to(index(i), index(i + r)) = true;
    
    % reacts_in(i, r) = -1 and reacts_in(j, r) = 1) iff there is a reaction r: y_i -> y_j
    reacts_in(index(i), i) = -1;
    reacts_in(index(i+r), i) = 1;
end

% Linkage classes
% Count number of connected components of an undirected graph
linkage_class = conncomp(graph(reacts_to | reacts_to'));

% Count the number of linkage classes
l = max(linkage_class);



%
% Strong linkage classes
%

% Check if the network is reversible
is_reversible = isequal(reacts_to, reacts_to');

% Strong linkage classes
% Count number of strongly connected components of an directed graph
if is_reversible
    strong_linkage_class = linkage_class;
else
    % Count number of connected components of a directed graph
    strong_linkage_class = conncomp(digraph(reacts_to));
end

% Count the number of strong linkage classes
sl = max(strong_linkage_class);



%
% Terminal linkage classes
%

% Count the number of terminal strong linkage classes
% Initialize the count
t = 0;
is_nontrivial_terminal_slc = false(sl, 1);
is_terminal_complex = false(n, 1);
for i = 1:sl
    
    % Locate the indexes in Y of the complexes present in strong-linkage class i
    complexes_i = find(strong_linkage_class == i);
    
    % Locate the indexes in Y of the complexes which in some reactions are products of complexes_i
    products_of_complexes_i = find(any(reacts_to(complexes_i, :), 1));
    
    % Products_of_complexes_i is a subset of complexes_i, so the strong-linkage class i is terminal
    if all(ismember(products_of_complexes_i, complexes_i))
        t = t + 1;
        is_terminal_complex(complexes_i) = true;
        if numel(complexes_i) > 1
            is_nontrivial_terminal_slc(i) = true;
        end
    end
end



%
% Rank
%

% Get the rank of the reaction network
% S = Im N
% dim S = dim (Im N) = rank(N)
% Note: We talk of "dimension of a linear transformation" and "rank of a matrix"
s = rank(N);



%
% Reactant rank
%

% Construct the incidence matrix
% We can decompose this into I_a = I_a^+ - I_a^-
I_a = reacts_in;

% Construct I_a^-: for reactant complexes only
I_a_minus = I_a;
I_a_minus(I_a_minus > 0) = 0;
I_a_minus(I_a_minus < 0) = 1;

% Construct N^-: "reactant subspace" matrix
N_minus = all_complex*I_a_minus;

% Get the reactant rank
% R = Im N^-
% dim R = dim (Im N^-) = rank(N^-)
q = rank(N_minus);


%
% Deficiency
%

% Compute the deficiency of the reaction network
delta = n - l - s;



%
% Reactant deficiency
%

% Compute the reactant deficiency
delta_p = n_r - q;



%
% Prepare table columns
%

% List of characteristics
characteristics = {'Species'; 'Complexes'; 'Reactant complexes'; 'Reversible reactions'; 'Irreversible reactions'; 'Reactions'; 'Linkage classes'; 'Strong linkage classes'; 'Terminal linkage classes'; 'Rank'; 'Reactant rank'; 'Deficiency'; 'Reactant deficiency'};

% List of notations
notation = {'m'; 'n'; 'n_r'; 'r_rev'; 'r_irrev'; 'r'; 'l'; 'sl'; 't'; 's'; 'q'; 'delta'; 'delta_p'};

% List of network numbers
network = {m; n; n_r; r_rev; r_irrev; r; l; sl; t; s; q; delta; delta_p};



%
% Display results
%

% Convert the list to a table
headers = {'Characteristics', 'Notation', 'Value'};
results = cell2table([characteristics, notation, network], 'VariableNames', headers);
results = convertvars(results, results.Properties.VariableNames, 'categorical');

% Print header
% Use 'fprintf' instead of 'disp' to interpret '\n' as 'newline'
fprintf('Network Numbers - %s\n\n', model.id);

% Display the table
disp(results);

end