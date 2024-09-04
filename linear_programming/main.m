addpath(genpath("cvx"))
addpath(genpath("/Users/luismbvarona/GitHub/" + ...
    "LaplacianBandwidth-Julia/linear_programming"))

% load("data/non_WHD_adjacencies_orders3to6.mat")
% load("data/norm_weak_hads_orders1to6.mat", ...
%     "norm_weak_hads_order3", "norm_weak_hads_order4", ...
%     "norm_weak_hads_order5", "norm_weak_hads_order6")

load("data/non_WHD_adjacencies_orders3to6.mat", "non_WHD_adjacencies_order6")
load("data/norm_weak_hads_orders1to6.mat", "norm_weak_hads_order6")
%%
% A3 = double(non_WHD_adjacencies_order3);
% A4 = double(non_WHD_adjacencies_order4);
% A5 = double(non_WHD_adjacencies_order5);
A6 = double(non_WHD_adjacencies_order6);

% WH3 = double(norm_weak_hads_order3);
% WH4 = double(norm_weak_hads_order4);
% WH5 = double(norm_weak_hads_order5);
WH6 = double(norm_weak_hads_order6);
%%
% data = {{A3, WH3}, {A4, WH4}, {A5, WH5}, {A6, WH6}};
% weighted_graphs_orders3to6 = {{}, {}, {}, {}};
% weighted_graphs_orders3to5 = {{}, {}, {}};
%%
% load("data/WHD_weighted_graphs_order6")
weighted_graphs_order6 = {};
%%
% for n = 1:3 % for n = 1:4, but order 6 may take a while
% adjacencies = data{n}{1};
% weak_hads = data{n}{2};

for j = 2:size(A6, 3) % Start at 1 next time
    A = A6(:, :, j);
%     A = adjacencies(:, :, j);
    L = adjacency_to_laplacian(A);
    num_edges = sum(A, "all") / 2;
    
%     for k = 1:length(WH6)
    for k = 1:2000
        W = WH6(:, :, k);
        if cond(W) > 50
            W
        end
%         W = weak_hads(:, :, k);
        if mod(k, 100) == 0
            [j, k]
        end
        
        cvx_begin quiet
            variable weights(num_edges)
            A_weighted = zeros(6) .* weights(1);
%             A_weighted = zeros(n + 2) .* weights(1);
            A_weighted(A == 1) = repmat(weights, [2, 1]);
            L_weighted = adjacency_to_laplacian(A_weighted);
            
            D = inv(W) * L_weighted * W;
            minimize max(weights);
            subject to
                weights >= 1;
                D == diag(diag(D));
        cvx_end
        
        if(cvx_optval < Inf)
            graph.adjacency = A;
            graph.laplacian = L;
            graph.eigvecs = W;
            graph.weights = weights;
            graph.adjacency_weighted = A_weighted;
            graph.laplacian_weighted = L_weighted;
            
            graph
            weighted_graphs_order6 = [weighted_graphs_order6, graph];
            break
        end
    end
end
% end
%%
% save("/Users/luismbvarona/GitHub/LaplacianBandwidth-Julia/" + ...
%     "linear_programming/data/WHD_weighted_graphs_orders3to5", ...
%     "weighted_graphs_orders3to5")

save("/Users/luismbvarona/GitHub/LaplacianBandwidth-Julia/" + ...
    "linear_programming/data/WHD_weighted_graphs_order6", ...
    "weighted_graphs_order6")
%%
% WHD_weighted_graphs_order5 = [];
% for b = weighted_graphs_orders3to5{3}
%     WHD_weighted_graphs_order5 = [WHD_weighted_graphs_order5, b{1}];
% end
% 
% save("/Users/luismbvarona/GitHub/LaplacianBandwidth-Julia/" + ...
%     "linear_programming/data/WHD_weighted_graphs_order5", ...
%     "WHD_weighted_graphs_order5")
%%
function L = adjacency_to_laplacian(A)
    L = diag(sum(A)) - A;
end