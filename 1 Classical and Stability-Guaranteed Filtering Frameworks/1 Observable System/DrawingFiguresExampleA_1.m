%%% Drawing Figures

close all;
clear;
clc;

load('experiment20210930')


%%  Figure 2
kSequence_selected = [0 1 2 3 6];

for k = kSequence_selected
    k
    
    %   True ranges
    [vertex_k_prior, nr_prior] = con2vert(G_k_prior_total{k+kIndexC}, theta_k_prior_total{k+kIndexC});
    [CH_k_prior, x_prior_volume(k+kIndexC)] = convhull(vertex_k_prior);
    
    [vertex_k_posterior, nr_posterior] = con2vert(G_k_posterior_total{k+kIndexC}, theta_k_posterior_total{k+kIndexC});
    [CH_k_posterior, x_posterior_volume(k+kIndexC)] = convhull(vertex_k_posterior);
    
    x_posterior_diameter(k+kIndexC) = diameter_conv(vertex_k_posterior);
    
    %   Set-membership filter A
    [vertex_k_prior_A, nr_prior_A] = con2vert(G_k_prior_total_A{k+kIndexC}, theta_k_prior_total_A{k+kIndexC});
    [CH_k_prior_A, x_prior_volume_A(k+kIndexC)] = convhull(vertex_k_prior_A);
    
    [vertex_k_posterior_A, nr_posterior_A] = con2vert(G_k_posterior_total_A{k+kIndexC}, theta_k_posterior_total_A{k+kIndexC});
    [CH_k_posterior_A, x_posterior_volume_A(k+kIndexC)] = convhull(vertex_k_posterior_A);
    
    x_posterior_diameter_A(k+kIndexC) = diameter_conv(vertex_k_posterior_A);
    x_posterior_estimation_gap_A(k+kIndexC) = estimation_gap(vertex_k_posterior_A, vertex_k_posterior);
    
    %   Set-membership filter B
    [vertex_k_prior_B, nr_prior_B] = con2vert(G_k_prior_total_B{k+kIndexC}, theta_k_prior_total_B{k+kIndexC});
    [CH_k_prior_B, x_prior_volume_B(k+kIndexC)] = convhull(vertex_k_prior_B);
    
    [vertex_k_posterior_B, nr_posterior_B] = con2vert(G_k_posterior_total_B{k+kIndexC}, theta_k_posterior_total_B{k+kIndexC});
    [CH_k_posterior_B, x_posterior_volume_B(k+kIndexC)] = convhull(vertex_k_posterior_B);
    
    x_posterior_diameter_B(k+kIndexC) = diameter_conv(vertex_k_posterior_B);
    x_posterior_estimation_gap_B(k+kIndexC) = estimation_gap(vertex_k_posterior_B, vertex_k_posterior);
    
    %   Set-membership filter C
    if k == 0
        [vertex_k_prior_C, nr_prior_C] = con2vert(G_k_prior_total_C{k+kIndexC}, theta_k_prior_total_C{k+kIndexC});
        [CH_k_prior_C, x_prior_volume_C(k+kIndexC)] = convhull(vertex_k_prior_C);
    end
    
    [vertex_k_posterior_C, nr_posterior_C] = con2vert(G_k_posterior_total_C{k+kIndexC}, theta_k_posterior_total_C{k+kIndexC});
    [CH_k_posterior_C, x_posterior_volume_C(k+kIndexC)] = convhull(vertex_k_posterior_C);
    
    x_posterior_diameter_C(k+kIndexC) = diameter_conv(vertex_k_posterior_C);
    x_posterior_estimation_gap_C(k+kIndexC) = estimation_gap(vertex_k_posterior_C, vertex_k_posterior);
    
    if k == 0
        figure,
        plot(vertex_k_prior(CH_k_prior, 1), vertex_k_prior(CH_k_prior, 2), '-s', vertex_k_prior_A(CH_k_prior_A, 1), vertex_k_prior_A(CH_k_prior_A, 2), '-o',...
            vertex_k_prior_B(CH_k_prior_B, 1), vertex_k_prior_B(CH_k_prior_B, 2), '-*', vertex_k_prior_C(CH_k_prior_C, 1), vertex_k_prior_C(CH_k_prior_C, 2), '-^', 'LineWidth', 1.2, 'MarkerSize', 8)
        legend('True', 'Alice', 'Bob', 'Carol')
        grid on;
        set(gca,'FontSize',12);
        xlabel('x^{(1)}-axis')
        ylabel('x^{(2)}-axis')
        axis([-4.5 4.5 -4.5 4.5])
%         axis([-5 5 -5 5])
    end
    
    if k > 0
        if k >= delta
            [vertex_k_OIT, nr_OIT] = con2vert(G_k_OIT_total{k+kIndexC}, theta_k_OIT_total{k+kIndexC});
            [CH_k_OIT, x_OIT_volume(k+kIndexC)] = convhull(vertex_k_OIT);
            figure,
            plot(vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2), '-s', vertex_k_posterior_A(CH_k_posterior_A, 1), vertex_k_posterior_A(CH_k_posterior_A, 2), '-o',...
                vertex_k_posterior_B(CH_k_posterior_B, 1), vertex_k_posterior_B(CH_k_posterior_B, 2), '-*', vertex_k_posterior_C(CH_k_posterior_C, 1), vertex_k_posterior_C(CH_k_posterior_C, 2), '-^',...
                vertex_k_OIT(CH_k_OIT, 1), vertex_k_OIT(CH_k_OIT, 2), '--', 'LineWidth', 1.2, 'MarkerSize', 8)
            legend('True', 'Alice', 'Bob', 'Carol', 'OIT')
            grid on;
            set(gca,'FontSize',12);
            xlabel('x^{(1)}-axis')
            ylabel('x^{(2)}-axis')
            if k == 6
                axis([4.5 7.5 -2.5 2])
            end
        else
            figure,
            plot(vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2), '-s', vertex_k_posterior_A(CH_k_posterior_A, 1), vertex_k_posterior_A(CH_k_posterior_A, 2), '-o',...
                vertex_k_posterior_B(CH_k_posterior_B, 1), vertex_k_posterior_B(CH_k_posterior_B, 2), '-*', vertex_k_posterior_C(CH_k_posterior_C, 1), vertex_k_posterior_C(CH_k_posterior_C, 2), '-^',...
                'LineWidth', 1.2, 'MarkerSize', 8)
            legend('True', 'Alice', 'Bob', 'Carol')
            grid on;
            set(gca,'FontSize',12);
            xlabel('x^{(1)}-axis')
            ylabel('x^{(2)}-axis')
        end
    end
end


%%  Maximum Diameters
MD = max(x_posterior_diameter(delta_optimal+kIndexC:end))
MD_A = max(x_posterior_diameter_A(delta_optimal+kIndexC:end))
MD_B = max(x_posterior_diameter_B(delta_optimal+kIndexC:end))
MD_C = max(x_posterior_diameter_C(delta_optimal+kIndexC:end))
diameter_upper_bound

k_draw_max = 50+kIndexC;

figure,
plot(kSequence(1:k_draw_max), x_posterior_diameter(1:k_draw_max), '-s', kSequence(1:k_draw_max), x_posterior_diameter_A(1:k_draw_max), '-o', kSequence(1:k_draw_max), x_posterior_diameter_B(1:k_draw_max), '-*', kSequence(1:k_draw_max), x_posterior_diameter_C(1:k_draw_max), '-^', [min(kSequence(2:k_draw_max)) max(kSequence(2:k_draw_max))], [diameter_upper_bound, diameter_upper_bound], '-', 'LineWidth', 1.2, 'MarkerSize', 6)
hold on;
plot([1 1], [0 9], '--')
legend('True', 'Alice', 'Bob', 'Carol', 'Upper Bound (k \geq 1)')
grid on;
set(gca,'FontSize',10);
ylim([0 9])
xlabel('Time Step')
ylabel('Diameter')