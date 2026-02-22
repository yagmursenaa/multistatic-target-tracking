T = 20;             % Period
c = 1500;           % Ses hızı (m/s)
total_time = 1000;  % Toplam süresi (saniye)
steps = total_time / T;

% helikopter konumu
helicopter_pos = [0, 0]; 

% hedefin gerçek konumu
X_initial = [-8000; 10; 8000; -5];
X_true = zeros(4, steps);
X_true(:,1) = X_initial;

F = [1, T, 0, 0;
     0, 1, 0, 0;
     0, 0, 1, T;
     0, 0, 0, 1];

for k = 1:steps-1
    X_true(:, k+1) = F * X_true(:, k);
end

% SENARYOLAR
scenario = struct();

% üçgensel sonoboy konumları
pos_triangle = [0, 10000; 4000, 2000; -5000, 0];

% doğrusal konumlar
pos_linear = [-5000, 0; 500, 0; 5000, 0];

% 1. Senaryo: üçgensel konumlanma ve her sonoboy için aynı sigma
scenario(1).name = 'S1: Üçgensel Konumlanma ve Aynı Sigma';
scenario(1).sonoboys = pos_triangle;
scenario(1).sigma = [0.1; 0.1; 0.1];

% 2. Senaryo: üçgensel konumlanma ve her sonoboy için farklı sigma
scenario(2).name = 'S2: Üçgensel Konumlanma ve Farklı Sigma';
scenario(2).sonoboys = pos_triangle;
scenario(2).sigma = [0.2; 0.1; 0.05];  

% 3. Senaryo: Doğrusal Konumlanma ve aynı sigma
scenario(3).name = 'S3: Doğrusal Konum ve Aynı Sigma';
scenario(3).sonoboys = pos_linear;
scenario(3).sigma = [0.1; 0.1; 0.1];

% 4. Senaryo: Doğrusal Konumlanma ve farklı sigma
scenario(4).name = 'S4: Doğrusal Konum ve Farklı Sigma';
scenario(4).sonoboys = pos_linear;
scenario(4).sigma = [0.2; 0.1; 0.05];

% Sonuçları saklamak için hücre dizisi
all_results = cell(1, 4);

% Seneryolar için döngü
for s_idx = 1:4
    
    sonoboys = scenario(s_idx).sonoboys;
    sigma_vec = scenario(s_idx).sigma; 
    num_sonoboys = size(sonoboys, 1);
    
    % Ölçüm verileri için
    Z_measurements = zeros(num_sonoboys, steps); 

    % R_cov Matrisi
    R_cov = diag(sigma_vec.^2);  

    for k = 1:steps
        x_t = X_true(1, k);
        y_t = X_true(3, k);
        
        for i = 1:num_sonoboys
            x_s = sonoboys(i, 1);
            y_s = sonoboys(i, 2);
            x_H = helicopter_pos(1);
            y_H = helicopter_pos(2);
            
            % Mesafeler 
            d_direct = sqrt((x_H - x_s)^2 + (y_H - y_s)^2); 
            d_HT    = sqrt((x_H - x_t)^2 + (y_H - y_t)^2); 
            d_TS    = sqrt((x_t - x_s)^2 + (y_t - y_s)^2); 
            
            % Gerçek Zaman Farkı
            delta_t = (d_HT + d_TS - d_direct) / c;
            
            % Ölçüm Gürültüsü Ekleme 
            Z_measurements(i, k) = delta_t + sigma_vec(i) * randn;
        end
    end 

    X_est = zeros(4, steps);
    X_est(:, 1) = [-8200; 0; 8200; 0]; % Başlangıç Önerisi

    P = diag([1e6, 500, 1e6, 500]); 

    Q_var = 0.1; 
    Q = [Q_var, 0; 
         0, Q_var]; 

    gamma = [0.5*T^2, 0;
             T,       0;
             0,       0.5*T^2;
             0,       T];

    for k = 1:steps-1
        X_pred = F * X_est(:, k); 
        P_pred = F * P * F.' + gamma * Q *gamma.' ;
        
        x_T = X_pred(1);
        y_T = X_pred(3);
        
        h_func = zeros(num_sonoboys, 1);
        H = zeros(num_sonoboys, 4);        
        
        x_H = helicopter_pos(1); 
        y_H = helicopter_pos(2);
        
        for i = 1:num_sonoboys
            x_s = sonoboys(i, 1);
            y_s = sonoboys(i, 2);
            
            R_HT = sqrt((x_T - x_H)^2 + (y_T - y_H)^2);
            R_TS = sqrt((x_T - x_s)^2 + (y_T - y_s)^2);
            d_direct = sqrt((x_H - x_s)^2 + (y_H - y_s)^2);
            
            h_func(i) = (R_HT + R_TS - d_direct) / c;
           
            %Jacobian matrisi H
            H(i, 1) = (1/c) * ((x_T - x_H)/R_HT + (x_T - x_s)/R_TS); 
            H(i, 2) = 0; 
            H(i, 3) = (1/c) * ((y_T - y_H)/R_HT + (y_T - y_s)/R_TS); 
            H(i, 4) = 0; 
        end
        
        S = H * P_pred * H.' + R_cov;
        W = P_pred * H.' / S; 
        
        inovation = Z_measurements(:, k+1) - h_func;
        
        X_est(:, k+1) = X_pred + W * inovation;       
        P = P_pred - W*S* W.' ;              
    end
    
    % Sonuçları kaydet
    pos_error = sqrt((X_true(1,:) - X_est(1,:)).^2 + (X_true(3,:) - X_est(3,:)).^2);
    result_struct.X_est = X_est;
    result_struct.pos_error = pos_error;
    result_struct.sonoboys = sonoboys;
    result_struct.rmse = sqrt(mean(pos_error.^2));
    all_results{s_idx} = result_struct;
end

% GRAFİKLER
colors = {'c', 'm'};

% 1. GRAFİK: S1, S2 YÖRÜNGE VE KONUMLAR
figure('Name', 'Grup 1: Üçgen Geometri Yörüngeleri (S1, S2)');
plot(X_true(1,:), X_true(3,:), 'k-', 'LineWidth', 2); hold on;
for i = 1:2
    plot(all_results{i}.X_est(1,:), all_results{i}.X_est(3,:), '--', 'Color', colors{i}, 'LineWidth', 1.5);
end
plot(pos_triangle(:,1), pos_triangle(:,2), 's', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
plot(helicopter_pos(1), helicopter_pos(2), 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
legend('Gerçek Yörünge', 'S1-aynı sigma', 'S2-farklı sigma', 'Sonoboylar', 'Helikopter');
title('Üçgensel geometri için kestirimler - S1 vs S2');

grid on; axis equal;

% 2. GRAFİK: S1, S2 HATA KARŞILAŞTIRMASI
figure('Name', 'Grup 1: Hata Analizi (S1, S2)');
hold on;
for i = 1:2
    plot(1:steps, all_results{i}.pos_error, 'Color', colors{i}, 'LineWidth', 1.5);
end
legend(['S1 (RMSE: ' num2str(all_results{1}.rmse, '%.1f') ')'], ...
       ['S2 (RMSE: ' num2str(all_results{2}.rmse, '%.1f') ')']);
title('Hata karşılaştırması üçgensel geometri- S1 vs S2');
grid on;

% 3. GRAFİK: S3, S4 YÖRÜNGE VE KONUMLAR (LİNEER)
figure('Name', 'Grup 2: Lineer Geometri Yörüngeleri (S3, S4)');
plot(X_true(1,:), X_true(3,:), 'k-', 'LineWidth', 2); hold on;
for i = 3:4
    plot(all_results{i}.X_est(1,:), all_results{i}.X_est(3,:), '--', 'Color', colors{i-2}, 'LineWidth', 1.5);
end
plot(pos_linear(:,1), pos_linear(:,2), 's', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
plot(helicopter_pos(1), helicopter_pos(2), 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
legend('Gerçek Yörünge', 'S3', 'S4', 'Sonoboylar', 'Helikopter');
title('Lineer geometri için kestirimler - S3 vs S4');

grid on; axis equal;

% 4. GRAFİK: S3, S4 HATA KARŞILAŞTIRMASI
figure('Name', 'Grup 2: Hata Analizi (S3, S4)');
hold on;
for i = 3:4
    plot(1:steps, all_results{i}.pos_error, 'Color', colors{i-2}, 'LineWidth', 1.5);
end
legend(['S3 (RMSE: ' num2str(all_results{3}.rmse, '%.1f') ')'], ...
       ['S4 (RMSE: ' num2str(all_results{4}.rmse, '%.1f') ')']);
title('Hata karşılaştırması lineer geometri- S3 vs S4');

grid on;

% 5. GRAFİK: S2 vs S4 Karşılaştırması
figure('Name', 'Karşılaştırma: S2 vs S4');
hold on;
plot(1:steps, all_results{2}.pos_error, 'Color', 'k', 'LineWidth', 2);
plot(1:steps, all_results{4}.pos_error, 'Color', 'b', 'LineWidth', 2);
legend(['S2: Üçgen (RMSE: ' num2str(all_results{2}.rmse, '%.1f') ')'], ...
       ['S4: Lineer (RMSE: ' num2str(all_results{4}.rmse, '%.1f') ')']);
title('Hata karşılaştırması - S2 vs S4 farklı sigmalar');
grid on;

% 6. GRAFİK: S1 vs S3 Karşılaştırması
figure('Name', 'Karşılaştırma: S1 vs S3');
hold on;
plot(1:steps, all_results{1}.pos_error, 'Color', 'k', 'LineWidth', 2);
plot(1:steps, all_results{3}.pos_error, 'Color', 'b', 'LineWidth', 2);
legend(['S1: Üçgen (RMSE: ' num2str(all_results{1}.rmse, '%.1f') ')'], ...
       ['S3: Lineer (RMSE: ' num2str(all_results{3}.rmse, '%.1f') ')']);
title('Hata karşılaştırması - S1 vs S3 aynı sigma');
grid on;