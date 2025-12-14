function MZI_Patel_Analysis()
% MZI_Patel_Analysis.m
% Replica as análises do paper Patel & Hamde (2017)
% Fig 2: Variação do ângulo vs Monitor Output
% Fig 6: Phase Shift vs Monitor Output

    clc; close all;
    
    % Menu de seleção
    fprintf('=== ANÁLISE MZI - PATEL & HAMDE (2017) ===\n');
    fprintf('1. Replicar Figura 2 (Varredura de Ângulo)\n');
    fprintf('2. Replicar Figura 6 (Phase Shift vs Output)\n');
    fprintf('3. Executar ambas as análises\n');
    choice = input('Escolha (1-3): ');
    
    switch choice
        case 1
            analyze_angle_sweep();
        case 2
            analyze_phase_sweep();
        case 3
            analyze_angle_sweep();
            analyze_phase_sweep();
        otherwise
            fprintf('Opção inválida!\n');
    end
end

%% ========== FIGURA 2: VARREDURA DE ÂNGULO ==========
function analyze_angle_sweep()
    fprintf('\n=== ANÁLISE: VARREDURA DE ÂNGULO (Figura 2) ===\n');
    
    % Varredura de ângulos (0° a 70° conforme paper)
    angles = linspace(0,70,100);
    n_angles = length(angles);
    
    monitor_output = zeros(1, n_angles);
    splitting_ratio = zeros(1, n_angles);
    
    hw = waitbar(0, 'Simulando varredura de ângulos...');
    
    for k = 1:n_angles
        angle = angles(k);
        [P_out, ratio] = simulate_mzi_angle(angle);
        monitor_output(k) = P_out;
        splitting_ratio(k) = ratio;
        waitbar(k/n_angles, hw);
    end
    close(hw);
    
    % Normalização
    monitor_norm = monitor_output / max(monitor_output);
    
    % Plotagem - Figura 2
    figure('Color', 'w', 'Position', [100 100 900 600]);
    
    subplot(2,1,1);
    plot(angles, monitor_norm, '-', 'LineWidth', 2.5, 'Color', 'b');
    hold on;
    
    % Destacar ângulo ótimo (18°)
    [~, idx_opt] = max(monitor_norm);
    plot(angles(idx_opt), monitor_norm(idx_opt), 'r*', ...
         'MarkerSize', 15, 'LineWidth', 2);
    
    grid on;
    xlabel('Y-Branch Angle (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Normalized Monitor Output', 'FontSize', 12, 'FontWeight', 'bold');
    title('Figura 2: Análise de Ângulo do Y-Branch', 'FontSize', 14);
    legend({'Monitor Output', sprintf('Ótimo: %.1f°', angles(idx_opt))}, ...
           'Location', 'northeast');
    xlim([0 75]);
    ylim([0 1.1]);
    
    % Adicionar linha de referência em 50%
    plot([0 75], [0.5 0.5], 'k--', 'LineWidth', 1);
    text(35, 0.52, '50% threshold', 'FontSize', 10);
    
    subplot(2,1,2);
    plot(angles, splitting_ratio*100, '-s', 'LineWidth', 2, ...
         'MarkerSize', 7, 'MarkerFaceColor', 'r', 'Color', 'r');
    hold on;
    plot([0 75], [50 50], 'k--', 'LineWidth', 1.5);
    grid on;
    xlabel('Y-Branch Angle (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Power Splitting Ratio (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Razão de Divisão de Potência', 'FontSize', 13);
    xlim([0 75]);
    ylim([0 100]);
    
    % Resultados
    fprintf('\n--- RESULTADOS DA VARREDURA DE ÂNGULO ---\n');
    fprintf('Ângulo Ótimo: %.1f° (Monitor Output = %.4f)\n', ...
            angles(idx_opt), monitor_norm(idx_opt));
    fprintf('Splitting Ratio no ótimo: %.2f%%\n', splitting_ratio(idx_opt)*100);
    fprintf('\nTabela de Resultados:\n');
    fprintf('Ângulo(°) | Monitor Output | Splitting(%%)\n');
    fprintf('----------|----------------|------------\n');
    for k = 1:n_angles
        fprintf('%6.1f    |    %.4f      |   %.2f\n', ...
                angles(k), monitor_norm(k), splitting_ratio(k)*100);
    end
end

%% ========== FIGURA 6: PHASE SHIFT vs MONITOR OUTPUT ==========
function analyze_phase_sweep()
    fprintf('\n=== ANÁLISE: PHASE SHIFT (Figura 6) ===\n');
    
    % Parâmetros
    lambda = 1.55;          % Comprimento de onda (μm)
    L_sensing = 600;        % Comprimento do braço de medição (μm)
    angle_opt = 1.8191;         % Ângulo ótimo identificado
    
    % Varredura de índice de refração
    % Para cobrir 0° a 180° de fase
    n_points = 25;
    phase_deg_target = linspace(0, 180, n_points);
    
    % Calcula Δn necessário para cada fase
    % Δφ = (2π * L * Δn) / λ  =>  Δn = (Δφ * λ) / (2π * L)
    phase_rad = deg2rad(phase_deg_target);
    dn_values = (phase_rad * lambda) / (2*pi * L_sensing);
    
    monitor_output = zeros(1, n_points);
    
    hw = waitbar(0, 'Simulando varredura de fase...');
    
    for k = 1:n_points
        [P_out, ~] = simulate_mzi_phase(dn_values(k), angle_opt);
        monitor_output(k) = P_out;
        waitbar(k/n_points, hw);
    end
    close(hw);
    
    % Normalização
    monitor_norm = monitor_output / max(monitor_output);
    
    % Teoria: P_out ∝ cos²(Δφ/2)
    % CORREÇÃO: phase_deg_target já está em graus, então usar cosd diretamente
    theory_output = cosd(phase_deg_target/2).^2;
    
    % Plotagem - Figura 6
    figure('Color', 'w', 'Position', [150 150 1000 700]);
    
    % Plot principal
    subplot(2,1,1);
    plot(phase_deg_target, monitor_norm, '-o', 'LineWidth', 2.5, ...
         'MarkerSize', 7, 'MarkerFaceColor', 'b', 'Color', 'b');
    hold on;
    plot(phase_deg_target, theory_output, 'r--', 'LineWidth', 2);
    
    % Identificar região linear (aproximadamente de 70° a 110°)
    linear_region = (phase_deg_target >= 70) & (phase_deg_target <= 110);
    plot(phase_deg_target(linear_region), monitor_norm(linear_region), ...
         'g-', 'LineWidth', 3.5, 'Color', [0 0.7 0]);
    
    grid on;
    xlabel('Phase Shift Δφ (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Normalized Monitor Output', 'FontSize', 12, 'FontWeight', 'bold');
    title('Figura 6: Phase Shift vs Monitor Output', 'FontSize', 14);
    legend({'Simulação BPM', 'Teórico: cos²(Δφ/2)', 'Região Linear'}, ...
           'Location', 'northeast', 'FontSize', 11);
    xlim([0 185]);
    ylim([-0.05 1.05]);
    
    % Adicionar marcações de região linear
    plot([70 70], [0 1], 'g--', 'LineWidth', 1);
    plot([110 110], [0 1], 'g--', 'LineWidth', 1);
    text(90, 0.1, 'Região Linear', 'FontSize', 11, 'Color', [0 0.5 0], ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    % Análise de linearidade
    subplot(2,1,2);
    
    % Calcular derivada (sensibilidade)
    dP_dphi = gradient(monitor_norm, phase_deg_target);
    plot(phase_deg_target, abs(dP_dphi), '-', 'LineWidth', 2, 'Color', [0.8 0.3 0.1]);
    hold on;
    plot(phase_deg_target(linear_region), abs(dP_dphi(linear_region)), ...
         'g-', 'LineWidth', 3, 'Color', [0 0.7 0]);
    
    grid on;
    xlabel('Phase Shift Δφ (degrees)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Sensitivity |dP/dφ|', 'FontSize', 12, 'FontWeight', 'bold');
    title('Sensibilidade do Sensor', 'FontSize', 13);
    xlim([0 185]);
    
    % Fit linear na região identificada
    p_linear = polyfit(phase_deg_target(linear_region), ...
                       monitor_norm(linear_region), 1);
    fit_y = polyval(p_linear, phase_deg_target(linear_region));
    R2 = 1 - sum((monitor_norm(linear_region) - fit_y).^2) / ...
             sum((monitor_norm(linear_region) - mean(monitor_norm(linear_region))).^2);
    
    % Resultados
    fprintf('\n--- RESULTADOS DA ANÁLISE DE FASE ---\n');
    fprintf('Região Linear Identificada: %.0f° a %.0f°\n', 70, 110);
    fprintf('Coeficiente de Linearidade (R²): %.4f\n', R2);
    fprintf('Sensibilidade na região linear: %.6f /degree\n', abs(p_linear(1)));
    fprintf('\nΔn correspondente à região linear:\n');
    fprintf('  Δn_min = %.6f (70°)\n', dn_values(find(linear_region, 1, 'first')));
    fprintf('  Δn_max = %.6f (110°)\n', dn_values(find(linear_region, 1, 'last')));
    
    % Erro RMS entre simulação e teoria
    rms_error = sqrt(mean((monitor_norm - theory_output).^2));
    fprintf('\nErro RMS (Simulação vs Teoria): %.6f\n', rms_error);
    
    %% Adicionar visualização de campo para um ponto específico
    fprintf('\nGerando mapa de campo para Δφ = 90° (máxima sensibilidade)...\n');
    idx_90 = find(abs(phase_deg_target - 90) < 5, 1);
    if ~isempty(idx_90)
        simulate_mzi_phase(dn_values(idx_90), angle_opt, true);
    end
end

%% ========== SIMULAÇÃO BPM: VARREDURA DE ÂNGULO ==========
function [P_out, splitting_ratio] = simulate_mzi_angle(angle_deg)
    % Simula MZI com ângulo específico (Δn = 0 para testar splitting)
    
    % Parâmetros físicos
    lambda = 1.55;
    k0 = 2*pi/lambda;
    n_core = 3.3847;
    n_clad = 3.378;
    n0 = n_clad;
    w = 3.0;
    
    % Conversão de ângulo
    angle_rad = deg2rad(angle_deg);
    
    % Grid
    dx = 0.15;
    X_width = 60;
    x = -X_width/2 : dx : X_width/2;
    Nx = length(x);
    
    dz = 0.5;
    L_split = 200;
    L_arm = 600;            
    L_combine = 200;
    L_tail = 0; 
    L_total = L_split + L_arm + L_combine + L_tail;
    
    z = 0 : dz : L_total;
    Nz = length(z);
    
    % Matrizes CN
    e = ones(Nx, 1);
    D2 = spdiags([e -2*e e], -1:1, Nx, Nx) / (dx^2);
    alpha = 1i * dz / (4 * k0 * n0);
    I_mat = speye(Nx);
    M_kin = alpha * D2;
    
    % Campo inicial (modo fundamental)
    sigma = w / 2;
    psi = exp(-(x.^2)/(sigma^2)).';
    psi = psi / sqrt(sum(abs(psi).^2) * dx);
    
    % PML
    mask_abc = ones(Nx, 1);
    border_px = round(5.0/dx);
    mask_abc(1:border_px) = linspace(0, 1, border_px);
    mask_abc(end-border_px+1:end) = linspace(1, 0, border_px);
    
    separation_max = L_split * tan(angle_rad);
    
    % Armazenar potência nos braços após o split
    P_left = 0;
    P_right = 0;
    sample_taken = false;
    
    % Loop BPM
    for i = 1:Nz
        z_curr = z(i);
        n_profile = ones(Nx, 1) * n_clad;
        
        if z_curr <= L_split
            offset = z_curr * tan(angle_rad);
            idx_core = (abs(x - offset) <= w/2) | (abs(x + offset) <= w/2);
            n_profile(idx_core) = n_core;
            
        elseif z_curr <= (L_split + L_arm)
            offset = separation_max;
            idx_left = abs(x + offset) <= w/2;
            idx_right = abs(x - offset) <= w/2;
            n_profile(idx_left) = n_core;
            n_profile(idx_right) = n_core;
            
            % Amostrar potência no meio do braço
            if ~sample_taken && z_curr > (L_split + L_arm/2)
                P_left = sum(abs(psi(idx_left)).^2) * dx;
                P_right = sum(abs(psi(idx_right)).^2) * dx;
                sample_taken = true;
            end
            
        elseif z_curr <= (L_split + L_arm + L_combine)
            z_local = z_curr - (L_split + L_arm);
            offset = separation_max - z_local * tan(angle_rad);
            if offset < 0, offset = 0; end
            idx_core = (abs(x - offset) <= w/2) | (abs(x + offset) <= w/2);
            n_profile(idx_core) = n_core;
            
        else
            idx_core = abs(x) <= w/2;
            n_profile(idx_core) = n_core;
        end
        
        V_vec = k0^2 * (n_profile.^2 - n0^2);
        M_pot = spdiags(alpha * V_vec, 0, Nx, Nx);
        
        A = I_mat - M_kin - M_pot;
        B = I_mat + M_kin + M_pot;
        psi = A \ (B * psi);
        psi = psi .* mask_abc;
    end
    
    % Overlap com modo fundamental na saída
    psi_mode = exp(-(x.^2)/(sigma^2)).';
    psi_mode = psi_mode / sqrt(sum(abs(psi_mode).^2) * dx);
    overlap = sum(psi .* conj(psi_mode)) * dx;
    P_out = abs(overlap)^2;
    
    % Splitting ratio
    P_total_arms = P_left + P_right;
    if P_total_arms > 0
        splitting_ratio = min(P_left, P_right) / P_total_arms;
    else
        splitting_ratio = 0;
    end
end

%% ========== SIMULAÇÃO BPM: VARREDURA DE FASE ==========
function [P_out, phase_actual] = simulate_mzi_phase(d_n_sensing, angle_deg, do_plot)
    % Simula MZI com diferença de índice específica
    
    if nargin < 3
        do_plot = false;
    end
    
    lambda = 1.55;
    k0 = 2*pi/lambda;
    n_core = 3.3847;
    n_clad = 3.378;
    n0 = n_clad;
    w = 3.0;
    angle_rad = deg2rad(angle_deg);
    
    dx = 0.12;
    X_width = 55;
    x = -X_width/2 : dx : X_width/2;
    Nx = length(x);
    
    dz = 0.5;
    L_split = 200;
    L_arm = 600;            
    L_combine = 200;
    L_tail = 0; 
    L_total = L_split + L_arm + L_combine + L_tail;
    
    z = 0 : dz : L_total;
    Nz = length(z);
    
    e = ones(Nx, 1);
    D2 = spdiags([e -2*e e], -1:1, Nx, Nx) / (dx^2);
    alpha = 1i * dz / (4 * k0 * n0);
    I_mat = speye(Nx);
    M_kin = alpha * D2;
    
    sigma = w / 2;
    psi = exp(-(x.^2)/(sigma^2)).';
    psi = psi / sqrt(sum(abs(psi).^2) * dx);
    psi_mode = psi;
    
    mask_abc = ones(Nx, 1);
    border_px = round(5.0/dx);
    mask_abc(1:border_px) = linspace(0, 1, border_px);
    mask_abc(end-border_px+1:end) = linspace(1, 0, border_px);
    
    separation_max = L_split * tan(angle_rad);
    
    % Matriz de intensidade para plot (se solicitado)
    if do_plot
        Intensity_Map = zeros(Nx, Nz);
    end
    
    for i = 1:Nz
        z_curr = z(i);
        n_profile = ones(Nx, 1) * n_clad;
        
        if z_curr <= L_split
            offset = z_curr * tan(angle_rad);
            idx_core = (abs(x - offset) <= w/2) | (abs(x + offset) <= w/2);
            n_profile(idx_core) = n_core;
            
        elseif z_curr <= (L_split + L_arm)
            offset = separation_max;
            idx_ref = abs(x + offset) <= w/2;
            idx_meas = abs(x - offset) <= w/2;
            n_profile(idx_ref) = n_core;
            n_profile(idx_meas) = n_core + d_n_sensing;
            
        elseif z_curr <= (L_split + L_arm + L_combine)
            z_local = z_curr - (L_split + L_arm);
            offset = separation_max - z_local * tan(angle_rad);
            if offset < 0, offset = 0; end
            idx_core = (abs(x - offset) <= w/2) | (abs(x + offset) <= w/2);
            n_profile(idx_core) = n_core;
            
        else
            idx_core = abs(x) <= w/2;
            n_profile(idx_core) = n_core;
        end
        
        V_vec = k0^2 * (n_profile.^2 - n0^2);
        M_pot = spdiags(alpha * V_vec, 0, Nx, Nx);
        
        A = I_mat - M_kin - M_pot;
        B = I_mat + M_kin + M_pot;
        psi = A \ (B * psi);
        psi = psi .* mask_abc;
        
        if do_plot
            Intensity_Map(:, i) = abs(psi).^2;
        end
    end
    
    overlap = sum(psi .* conj(psi_mode)) * dx;
    P_out = abs(overlap)^2;
    
    % Fase real acumulada (em radianos, depois convertida)
    phase_actual = (2*pi * L_arm * d_n_sensing) / lambda;
    
    % Plotagem do campo (se solicitado)
    if do_plot
        phase_deg = rad2deg(phase_actual);
        
        figure('Name', 'Visualização de Campo BPM', 'Color', 'w', ...
               'Position', [200 100 1000 600]);
        
        imagesc(z, x, Intensity_Map);
        colormap('hot'); 
        hcb = colorbar; 
        title(hcb, 'Intensidade', 'FontSize', 10);
        set(gca, 'YDir', 'normal'); 
        
        xlabel('Propagação Z (\mum)', 'FontSize', 12, 'FontWeight', 'bold'); 
        ylabel('Posição Transversal X (\mum)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('MZI Completo (\\Deltan = %.6f, \\Delta\\phi = %.1f°)', ...
              d_n_sensing, phase_deg), 'FontSize', 13);
        
        % Linhas guia visuais para as seções
        hold on;
        y_lim = get(gca, 'YLim');
        
        % Linha vertical: fim do split
        plot([L_split, L_split], y_lim, '--', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);
        text(L_split, y_lim(2)*0.95, ' Split', 'Color', 'white', ...
             'FontSize', 9, 'FontWeight', 'bold');
        
        % Linha vertical: fim dos braços
        plot([L_split+L_arm, L_split+L_arm], y_lim, '--', ...
             'Color', [0.9 0.9 0.9], 'LineWidth', 2);
        text(L_split+L_arm, y_lim(2)*0.95, ' Combine', 'Color', 'white', ...
             'FontSize', 9, 'FontWeight', 'bold');
        
        % Linha vertical: fim do combine
        plot([L_split+L_arm+L_combine, L_split+L_arm+L_combine], y_lim, '--', ...
             'Color', [0.9 0.9 0.9], 'LineWidth', 2);
        text(L_split+L_arm+L_combine, y_lim(2)*0.95, ' Output', 'Color', 'white', ...
             'FontSize', 9, 'FontWeight', 'bold');
        
        % Anotações das regiões
        text(L_split/2, y_lim(1)*0.9, 'Y-Splitter', ...
             'Color', 'yellow', 'FontSize', 11, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center');
        
        text(L_split + L_arm/2, separation_max + 3, 'Sensing Arm', ...
             'Color', 'cyan', 'FontSize', 10, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center');
        
        text(L_split + L_arm/2, -separation_max - 3, 'Reference Arm', ...
             'Color', 'cyan', 'FontSize', 10, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center');
        
        text(L_split + L_arm + L_combine/2, y_lim(1)*0.9, 'Y-Combiner', ...
             'Color', 'yellow', 'FontSize', 11, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center');
        
        grid on;
        set(gca, 'GridColor', [0.3 0.3 0.3], 'GridAlpha', 0.3);
    end
end