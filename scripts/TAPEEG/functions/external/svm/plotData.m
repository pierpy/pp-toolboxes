function plotData(X, y)
%PLOTDATA Plots the data points X and y into a new figure 
%   PLOTDATA(x,y) plots the data points with + for the positive examples
%   and o for the negative examples. X is assumed to be a Mx2 matrix.
%
% Note: This was slightly modified such that it expects y = 1 or y = 0

% Find Indices of Positive and Negative Examples
pos = find(y == 1); neg = find(y == 0);

% Plot Examples
if size(X,2) == 2
    plot(X(pos, 1), X(pos, 2), 'k+','LineWidth', 1, 'MarkerSize', 7)
    hold on;
    plot(X(neg, 1), X(neg, 2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7)
    hold off;
elseif size(X,2) == 3
    plot3(X(pos, 1), X(pos, 2), X(pos, 3), 'k+','LineWidth', 1, 'MarkerSize', 7)
    hold on;
    plot(X(neg, 1), X(neg, 2), X(pos, 3), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7)
    hold off;
end
