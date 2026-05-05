%% Steady-state U and UP values across buffering implementations
clear; clc; close all;

% Labels
labels = {
    'All buffering'
    'DNA + spRNA'
    'DNA + Dichotomous'
    'spRNA + Dichotomous'
    'DNA only'
    'spRNA only'
    'Dichotomous only'
    'No buffering'
};

% Steady-state values from table
U = [
    454.2475
    371.9712
    245.4750
    454.2475
    107.2355
    371.9712
    245.4750
    107.2355
];

UP = [
    1.1809
    2.3618
    1.1809
    1.1809
    2.3618
    2.3618
    1.1809
    2.3618
];

x = 1:numel(labels);

%% Plot
figure('Color', 'w', 'Position', [100 100 1100 700]);
t = tiledlayout(2,1);
t.TileSpacing = 'tight';   % reduces space between plots
t.Padding = 'tight';       % reduces outer margins

% Panel A: U
nexttile;
bar(x, U, 'FaceColor', [255 145 53]./255, 'LineWidth', 1.5);
ylabel('[U^*] (\muM)');
title('A. Steady-State Unphosphorylated Response Regulator, U^*', 'FontSize',18);
ylim([0 500]);
grid on;
box off;
set(gca, 'FontSize', 12, 'XTick', x, 'XTickLabel', []);

% Add values above bars
for i = 1:numel(U)
    text(x(i), U(i) + 22, sprintf('%.1f', U(i)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% Panel B: UP
nexttile;
bar(x, UP, 'FaceColor', [230 91 0]./255, 'LineWidth', 1.5);
ylabel('[U_P^*] (\muM)');
title('B. Steady-State Phosphorylated Response Regulator, U_P^*', 'FontSize',18);
ylim([0 2.8]);
grid on;
box off;
set(gca, 'FontSize', 12, ...
    'XTick', x, ...
    'XTickLabel', labels, ...
    'XTickLabelRotation', 35);

% Add values above bars
for i = 1:numel(UP)
    text(x(i), UP(i) + 0.12, sprintf('%.2f', UP(i)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

%% Overall title
% sgtitle('Steady-state comparison across buffering implementations', ...
%     'FontSize', 15, ...
%     'FontWeight', 'bold');

%% Optional: save high-resolution figure
% exportgraphics(gcf, 'U_UP_buffering_subplots.png', 'Resolution', 600);