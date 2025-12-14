function MZI_Optimize_Angle()
% MZI_Optimize_Angle_V2.m
% Estratégia Otimizada: Foco na região de baixo ângulo onde o guia funciona.
    clc; close all;
    fprintf('=== INICIANDO ROTINA DE OTIMIZAÇÃO E MAPEAMENTO DE ERRO ===\n');

    % --- 1. Mapeamento da Função de Custo para Análise ---
    fprintf('\n[1/4] Mapeando a função de custo (Angle vs. Erro) de 0 a 20 graus...\n');
    [angles_scan, errors_scan] = plot_cost_landscape(100); 
    
    % Encontra o melhor ponto da varredura (Scan)
    [min_error_scan, idx_scan_best] = min(errors_scan);
    angle_scan_best = angles_scan(idx_scan_best);

    % --- 2. Definição do Espaço de Busca Otimizado ---
    lower_bound = 0.001; 
    upper_bound = 2.5; 
    
    % Define a função objetivo (Cost Function)
    objective_fcn = @(angle_deg) MZI_CostFunction(angle_deg);
    
    % Otimizador fminbnd: Busca fina e precisa
    options = optimset('Display', 'iter', 'TolX', 1e-5); 
    
    t_start_opt = tic;
    fprintf('\n[2/4] Executando Otimização Fina fminbnd (%.3f a %.1f graus)...\n', lower_bound, upper_bound);
    [optimized_angle_fminbnd, min_error_fminbnd] = fminbnd(objective_fcn, lower_bound, upper_bound, options);
    t_total_opt = toc(t_start_opt);
    
    % --- 3. Seleção do Melhor Resultado Global ---
    if min_error_fminbnd < min_error_scan
        final_best_angle = optimized_angle_fminbnd;
        final_min_error = min_error_fminbnd;
        fprintf('\n[3/4] Melhor resultado encontrado pelo fminbnd.\n');
    else
        final_best_angle = angle_scan_best;
        final_min_error = min_error_scan;
        fprintf('\n[3/4] Melhor resultado encontrado na varredura inicial.\n');
    end

    fprintf('\n--- RESULTADO FINAL DA OTIMIZAÇÃO ---\n');
    fprintf('Melhor Ângulo Global (angle_deg): %.5f graus\n', final_best_angle);
    fprintf('Erro Mínimo Global (SSD): %.4e\n', final_min_error);
    fprintf('Otimização concluída em %.2f segundos.\n', t_total_opt);
    
    % --- 4. Geração do Resultado Final ---
    fprintf('\n[4/4] Gerando resultado final com o ângulo OTIMIZADO (%.5f graus)...\n', final_best_angle);
    MZI_Full_Simulation(final_best_angle);
end

%% --- FUNÇÃO ADICIONAL: MAPEAMENTO DA FUNÇÃO DE CUSTO ---
function [angles_scan, errors_scan] = plot_cost_landscape(N_points)
    % Mapeia a função de custo em um intervalo amplo para identificar a região de interesse.
    angles_scan = linspace(0.001, 20, N_points);
    errors_scan = zeros(1, N_points);
    
    for k = 1:N_points
        errors_scan(k) = MZI_CostFunction(angles_scan(k));
    end
    
    figure('Color', 'w', 'Name', 'Topologia da Função de Custo (Erro vs. Ângulo)');
    plot(angles_scan, errors_scan, 'b-o', 'LineWidth', 2, 'MarkerSize', 5);
    hold on;
    
    % Plota a região de "Guia Quebrado" (erros altos)
    idx_high_error = errors_scan > 5;
    plot(angles_scan(idx_high_error), errors_scan(idx_high_error), 'r.', 'MarkerSize', 15);
    
    grid on;
    xlabel('Ângulo do Y-Branch (\theta, Degrees)', 'FontSize', 12);
    ylabel('Erro (SSD) entre Simulação e Teoria', 'FontSize', 12);
    title('Função de Custo do Ângulo de Acoplamento', 'FontSize', 14);
    
    % Limita o eixo Y para visualizar melhor os erros pequenos
    if max(errors_scan) > 100
        % Adiciona uma linha pontilhada para mostrar o corte
        y_cut_lim = max(errors_scan(errors_scan < 10)) * 1.5 + 0.1;
        ylim([0 y_cut_lim]); 
        plot(xlim, [1e10/2, 1e10/2], 'k:', 'DisplayName', 'Erro de Corte'); % Linha indicando o erro muito alto
    else
        ylim([0 max(errors_scan) * 1.5]);
    end
end

%% --- FUNÇÃO DE CUSTO (OBJETIVO) ---
function error = MZI_CostFunction(angle_deg)
% Roda a varredura BPM para um dado angle_deg e retorna o erro (SSD)
    
    % 1. Configuração da Varredura (Parâmetros fixos)
    n_points = 15;
    dn_max = 0.00135; 
    dn_values = linspace(0, dn_max, n_points); 
    output_powers = zeros(1, n_points);
    
    % 2. Loop de Simulação
    for k = 1:n_points
        output_powers(k) = run_mzi_bpm(dn_values(k), false, angle_deg); 
    end

    % 3. Cálculo da Curva e Erro
    L_sensing = 600; 
    lambda = 1.55;
    
    % Cálculo da fase (Graus)
    phase_shift_rad = (2*pi * L_sensing * dn_values) / lambda;
    phase_shift_deg = rad2deg(phase_shift_rad);
    
    % Normalização (Essencial para comparação)
    try
        max_power = max(output_powers);
        % Se o max power for muito baixo (guia quebrado), retorna erro alto
        if max_power < 1e-6
            error = 1e10;
            return;
        end
        output_norm = output_powers / max_power;
    catch
        error = 1e10; 
        return;
    end
    
    % Curva Teórica (Lei do Cosseno) avaliada APENAS nos pontos de fase simulados
    power_theory = cosd(phase_shift_deg/2).^2;
    
    % Cálculo do Erro: Soma do Quadrado das Diferenças (SSD)
    error = sum((output_norm - power_theory).^2);
end

%% --- FUNÇÕES LOCAIS (MZI_Full_Simulation e run_mzi_bpm) ---

function MZI_Full_Simulation(optimized_angle)
% Simulação Final Otimizada: Agora aceita o 'optimized_angle'
    fprintf('=== INICIANDO SIMULAÇÃO FINAL COM ÂNGULO DE %.4f GRAUS ===\n', optimized_angle);
    %% 1. Configuração da Varredura
    n_points = 15;
    dn_max = 0.00135; 
    dn_values = linspace(0, dn_max, n_points); 
    output_powers = zeros(1, n_points);
    
    %% 2. Loop de Simulação
    fprintf('Simulando varredura (Filtragem Modal Ativa)...\n');
    hw = waitbar(0, 'Processando BPM...');
    t_start = tic;
    
    for k = 1:n_points
        output_powers(k) = run_mzi_bpm(dn_values(k), false, optimized_angle);
        waitbar(k/n_points, hw);
    end
    close(hw);
    fprintf('Simulação concluída em %.2f segundos.\n', toc(t_start));
    %% 3. Plotagem e Comparação
    L_sensing = 600; 
    lambda = 1.55;
    
    phase_shift_rad = (2*pi * L_sensing * dn_values) / lambda;
    phase_shift_deg = rad2deg(phase_shift_rad);
    
    output_norm = output_powers / max(output_powers);
    
    figure('Color', 'w', 'Name', 'Resultado Final MZI Otimizado');
    
    plot(phase_shift_deg, output_norm, '-o', 'LineWidth', 2, ...
         'MarkerFaceColor', 'b', 'Color', 'b', 'MarkerSize', 6);
    hold on;
    
    phase_theory = linspace(0, 180, 100);
    power_theory = cosd(phase_theory/2).^2;
    plot(phase_theory, power_theory, 'r--', 'LineWidth', 2);
    
    grid on;
    xlabel('Phase Shift \Delta\phi (Degrees)', 'FontSize', 12);
    ylabel('Normalized Output Power', 'FontSize', 12);
    title('Resposta do Interferômetro Mach-Zehnder Otimizado', 'FontSize', 14);
    subtitle(['Ângulo Otimizado: ' num2str(optimized_angle, '%.4f') ' graus'], 'FontSize', 10);
    legend({'Simulação BPM (Otimizada)', 'Teórico'}, 'Location', 'SouthWest', 'FontSize', 11);
    ylim([-0.05 1.05]);
    xlim([0 185]);
    
    fprintf('Gerando mapa de campo para visualização (com ângulo de %.4f graus)...\n', optimized_angle);
    run_mzi_bpm(0, true, optimized_angle); 
end

function [P_out, Intensity_Map] = run_mzi_bpm(d_n_sensing, do_plot, angle_deg)
% run_mzi_bpm: Agora aceita o 'angle_deg' como entrada.
    % Parâmetros Físicos
    lambda = 1.55;          
    k0 = 2*pi/lambda;       
    n_core = 3.3847;        
    n_clad = 3.378;         
    n0 = n_clad;            
    
    w = 3.0;                
    angle_rad = deg2rad(angle_deg);
    % Grid Refinado (dx menor para precisão de fase)
    dx = 0.12;              
    X_width = 50;           
    x = -X_width/2 : dx : X_width/2;
    Nx = length(x);
    
    % Passo Z refinado para evitar erro numérico acumulado
    dz = 0.5;               
    
    % Geometria
    L_split = 200;
    L_arm = 600;            
    L_combine = 200;
    L_tail = 0; 
    L_total = L_split + L_arm + L_combine + L_tail;
    
    z = 0 : dz : L_total;
    Nz = length(z);
    % Matrizes Crank-Nicolson
    e = ones(Nx, 1);
    D2 = spdiags([e -2*e e], -1:1, Nx, Nx) / (dx^2);
    alpha = 1i * dz / (4 * k0 * n0);
    
    I_mat = speye(Nx);
    M_kin = alpha * D2;
    
    % --- CAMPO INICIAL & MODO FUNDAMENTAL ---
    sigma = w / 2; 
    psi_mode = exp(-(x.^2)/(sigma^2)).'; 
    psi_mode = psi_mode / sqrt(sum(abs(psi_mode).^2) * dx); 
    
    psi = psi_mode; 
    
    % Absorção nas bordas
    mask_abc = ones(Nx, 1);
    border_px = round(4.0/dx); 
    mask_abc(1:border_px) = linspace(0, 1, border_px);
    mask_abc(end-border_px+1:end) = linspace(1, 0, border_px);
    separation_max = L_split * tan(angle_rad); 
    
    if do_plot, Intensity_Map = zeros(Nx, Nz); else, Intensity_Map = []; end
    % --- LOOP BPM ---
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
            % Tail
            idx_core = (abs(x) <= w/2);
            n_profile(idx_core) = n_core;
        end
        
        V_vec = k0^2 * (n_profile.^2 - n0^2);
        M_pot = spdiags(alpha * V_vec, 0, Nx, Nx);
        
        A = I_mat - M_kin - M_pot;
        B = I_mat + M_kin + M_pot;
        psi = A \ (B * psi);
        psi = psi .* mask_abc;
        
        if do_plot, Intensity_Map(:, i) = abs(psi).^2; end
    end
    
    % --- CÁLCULO DE POTÊNCIA (OVERLAP INTEGRAL) ---
    overlap = sum(psi .* conj(psi_mode)) * dx;
    P_out = abs(overlap)^2;
    
    % Plotagem (Se solicitado)
    if do_plot
        figure('Name', 'Visualização de Campo BPM', 'Color', 'w');
        imagesc(z, x, Intensity_Map);
        colormap('hot'); 
        hcb = colorbar; title(hcb, 'Intensidade');
        set(gca, 'YDir', 'normal'); 
        xlabel('Propagação Z (\mu m)'); 
        ylabel('Posição Transversal X (\mu m)');
        title(['MZI Completo (\Delta n = ' num2str(d_n_sensing) ', \theta = ' num2str(angle_deg, '%.4f') '^{\circ})']);
        
        % Linhas guia visuais
        hold on;
        y_lim = get(gca, 'YLim');
        plot([L_split, L_split], y_lim, '--', 'Color', [0.8 0.8 0.8]);
        plot([L_split+L_arm, L_split+L_arm], y_lim, '--', 'Color', [0.8 0.8 0.8]);
        plot([L_split+L_arm+L_combine, L_split+L_arm+L_combine], y_lim, '--', 'Color', [0.8 0.8 0.8]);
    end
end