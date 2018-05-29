% -------------------------------------------------------------------------
% [Ben] 05/28/18
% Produces a bar graph, where the leftmost bar corresponds to the postoral
% row, and where the height of each bar indicates the no. of BBs
% corresponding to that CR. We assume that the cell is oriented such that
% its anteroposterior axis is vertical, with the anterior end "higher" than
% the posterior end. We move from the postoral row through successive rows
% by going anti-clockwise about the anteroposterior axis.
% -------------------------------------------------------------------------

function visualize_CR_lengths(traceback)

row_lengths = sum((traceback > 0), 2);
mu = mean(row_lengths);
sigma = std(row_lengths);
bar(row_lengths, 'FaceColor', [0.6 0.6 0.6]);
hold on
% draw mean line 
pl_mu = line(xlim, [mu mu], 'Color', 'r', 'LineWidth', 2);
% draw lines demarcating +/- 1 standard deviation
pl_mu_plus = line(xlim, [mu mu]+sigma, 'Color', 'b', 'LineWidth', 1.5);
line(xlim, [mu mu]-sigma, 'Color', 'b', 'LineWidth', 1.5);
legend([pl_mu pl_mu_plus], {'mean', 'mean +/- 1 sd'});
xlabel('Row no.');
ylabel('Row length (no. of BBs)');
title('Graph of ciliary row lengths');
hold off
end