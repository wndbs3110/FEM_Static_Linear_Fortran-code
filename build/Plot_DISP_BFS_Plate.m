clear all; clc; close all;
% =========================================================================
E2=1.0e9;
thick=0.06;
load=75.0;
width=20.0;
% =========================================================================
% =========================================================================
% 1. 파일 열기 및 헤더(Header) 읽기
% =========================================================================
fid = fopen('SATURN_pos.txt','r');
if fid == -1
    error('SATURN_pos.txt 파일을 찾을 수 없습니다.');
end
n            = fscanf(fid, '%i', 1);   % 요소 개수 (size(elem))
nUniqueNodes = fscanf(fid, '%i', 1);   % 전체 노드 개수 (size(node))
nLYR         = fscanf(fid, '%i', 1);   % 적층판 레이어 개수 (nLYR)
scalefactor  = fscanf(fid, '%g', 1);   % 스케일 팩터
fprintf('>> 데이터 로드 완료: 요소=%d개, 노드=%d개, 적층판=%d Layers\n', n, nUniqueNodes, nLYR);
% 개별 응력(SXX, SYY, SXY)을 상세히 확인할 기준 레이어 선택
plot_layer = 1; 
if plot_layer > nLYR || plot_layer < 1
    error('plot_layer는 1과 %d 사이의 값이어야 합니다.', nLYR);
end
% =========================================================================
% 2. 요소별 데이터 읽기 (좌표/변위 + 다중 레이어 응력)
% =========================================================================
nNodes_patch = 4 * n;
Xall    = zeros(nNodes_patch, 1);
Yall    = zeros(nNodes_patch, 1);
Wall    = zeros(nNodes_patch, 1);
Dispall = zeros(nNodes_patch, 1);
sxx_top = zeros(nLYR, nNodes_patch);
sxx_bot = zeros(nLYR, nNodes_patch);
syy_top = zeros(nLYR, nNodes_patch);
syy_bot = zeros(nLYR, nNodes_patch);
sxy_top = zeros(nLYR, nNodes_patch);
sxy_bot = zeros(nLYR, nNodes_patch);
svm_top = zeros(nLYR, nNodes_patch);
svm_bot = zeros(nLYR, nNodes_patch);
for j = 1:n
    % 1) 좌표 + 변위 (32개)
    data32 = fscanf(fid, '%g', 32);
    for i = 1:4
        idx = i + 4*(j-1);
        offset = 8 * (i - 1);
        Xall(idx)    = data32(offset + 1);
        Yall(idx)    = data32(offset + 2);
        Wall(idx)    = data32(offset + 5);
        Dispall(idx) = abs(data32(offset + 5));
    end
    
    % 2) 레이어별 응력 읽기
    for m = 1:nLYR
        stress12 = fscanf(fid, '%g', 24); 
        for i = 1:4
            idx = i + 4*(j-1);
            
            stress_offset = 6 * (i - 1);
            
            sxx_val_top = stress12(stress_offset + 1);
            syy_val_top = stress12(stress_offset + 2);
            sxy_val_top = stress12(stress_offset + 3);
            
            sxx_val_bot = stress12(stress_offset + 4);
            syy_val_bot = stress12(stress_offset + 5);
            sxy_val_bot = stress12(stress_offset + 6);
            
            sxx_top(m, idx) = sxx_val_top;
            syy_top(m, idx) = syy_val_top;
            sxy_top(m, idx) = sxy_val_top;
            svm_top(m, idx) = sqrt(sxx_val_top^2 - sxx_val_top*syy_val_top + syy_val_top^2 + 3*sxy_val_top^2);
            
            sxx_bot(m, idx) = sxx_val_bot;
            syy_bot(m, idx) = syy_val_bot;
            sxy_bot(m, idx) = sxy_val_bot;
            svm_bot(m, idx) = sqrt(sxx_val_bot^2 - sxx_val_bot*syy_val_bot + syy_val_bot^2 + 3*sxy_val_bot^2);
        end
    end
end
% =========================================================================
% 3. 고유 노드 변위 데이터 읽기
% =========================================================================
nodal_data_raw = fscanf(fid, '%g', [7, nUniqueNodes]);
fclose(fid);
nodal_data = nodal_data_raw'; 
exact_U = nodal_data(:, 2);
exact_V = nodal_data(:, 3);
exact_W = nodal_data(:, 4);
exact_Mag = sqrt(exact_W.^2+exact_V.^2+exact_U.^2);
% =========================================================================
% 4. 플롯 데이터 및 Figure 이름 동적 구성
% =========================================================================
fig_names = {};
patchData = {};
fig_idx = 1;
for m = 1:nLYR
    fig_names{fig_idx} = sprintf('Von Mises Stress (Layer %d - Top)', m);
    patchData{fig_idx} = svm_top(m, :)';
    fig_idx = fig_idx + 1;
    
    fig_names{fig_idx} = sprintf('Von Mises Stress (Layer %d - Bottom)', m);
    patchData{fig_idx} = svm_bot(m, :)';
    fig_idx = fig_idx + 1;
end
fig_names{fig_idx} = sprintf('Normal Stress SXX (Layer %d - Top)', plot_layer);
patchData{fig_idx} = sxx_top(plot_layer, :)';
fig_idx = fig_idx + 1;
fig_names{fig_idx} = sprintf('Normal Stress SYY (Layer %d - Top)', plot_layer);
patchData{fig_idx} = syy_top(plot_layer, :)';
fig_idx = fig_idx + 1;
fig_names{fig_idx} = sprintf('Shear Stress SXY (Layer %d - Top)', plot_layer);
patchData{fig_idx} = sxy_top(plot_layer, :)';
fig_idx = fig_idx + 1;
fig_names{fig_idx} = 'Deflection W (Element)';
patchData{fig_idx} = Wall;
fig_idx = fig_idx + 1;
fig_names{fig_idx} = 'Total Displacement |W|';
patchData{fig_idx} = Dispall;
% =========================================================================
% 5. Figure 설정 및 렌더링 (3D 수정 부분)
% =========================================================================
nFig = numel(fig_names);
hFig = gobjects(nFig, 1);
hAx  = gobjects(nFig, 1);
hCb  = gobjects(nFig, 1);
X_patch = reshape(Xall, 4, n);
Y_patch = reshape(Yall, 4, n);
W_patch = reshape(Wall, 4, n); % 💡 Z축으로 사용할 처짐(W) 데이터 생성

for f = 1:nFig
    hFig(f) = figure(f);
    set(hFig(f), 'Name', fig_names{f});
    hAx(f)  = gca;
    
    hCb(f)  = colorbar;
    colormap(hFig(f), jet(100));
    xlabel('X'); ylabel('Y'); zlabel('W (Deflection)'); % 💡 Z라벨 추가
    hold(hAx(f), 'on');
    
    C_patch = reshape(patchData{f}, 4, n);
    
    % 💡 patch에 ZData 추가하여 3D로 플롯
    % (변위가 너무 작아 눈에 안 띈다면 W_patch 대신 W_patch * scalefactor 적용 가능)
    patch(hAx(f), 'XData', X_patch, 'YData', Y_patch, 'ZData', W_patch, 'CData', C_patch, ...
          'FaceColor', 'interp', 'EdgeColor','none'); 
    
    view(3);             % 💡 3차원 뷰로 설정
    grid(hAx(f), 'on');  % 💡 3D 공간을 잘 볼 수 있게 그리드 켬
    axis(hAx(f), 'tight');
    
    % 처짐량이 X, Y 스케일에 비해 매우 작으면 평면처럼 보일 수 있습니다.
    % 3D 왜곡을 보려면 아래 axis equal 설정은 주석 처리하거나 배율을 조정하는 것이 좋습니다.
    % axis(hAx(f), 'equal'); 
    
    dmin = min(patchData{f}); 
    dmax = max(patchData{f});
    if dmax == dmin
        dmax = dmin + eps; 
    end
    clim(hAx(f), [dmin, dmax]);
    hCb(f).Ticks = linspace(dmin, dmax, 10);
    hCb(f).FontSize = 11;
    
    title(hAx(f), fig_names{f}, 'Interpreter', 'none');
end
% =========================================================================
% 6. 통계 계산 및 결과 출력 (터미널)
% =========================================================================
fprintf('\n==========================================================\n');
fprintf('  Post-Processing Results  (Scale=%.4g)\n', scalefactor);
fprintf('==========================================================\n\n');
fprintf(' [ Von Mises Stress by Layer ]\n');
for m = 1:nLYR
    d_top = svm_top(m, :)';
    fprintf(' Layer %d (Top) : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', m, max(d_top), min(d_top), mean(d_top));
    
    d_bot = svm_bot(m, :)';
    fprintf(' Layer %d (Bot) : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', m, max(d_bot), min(d_bot), mean(d_bot));
end
fprintf('\n');
fprintf(' [ Directional Stress (Layer %d - Top) ]\n', plot_layer);
d_sxx = sxx_top(plot_layer, :)'; fprintf(' %-10s : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', 'SXX', max(d_sxx), min(d_sxx), mean(d_sxx));
d_syy = syy_top(plot_layer, :)'; fprintf(' %-10s : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', 'SYY', max(d_syy), min(d_syy), mean(d_syy));
d_sxy = sxy_top(plot_layer, :)'; fprintf(' %-10s : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', 'SXY', max(d_sxy), min(d_sxy), mean(d_sxy));
fprintf('\n');
fprintf(' [ Displacement Information (Exact Nodal Average) ]\n');
fprintf(' %-10s : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', 'U (In-plane)', max(exact_U), min(exact_U), mean(exact_U));
fprintf(' %-10s : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', 'V (In-plane)', max(exact_V), min(exact_V), mean(exact_V));
fprintf(' %-10s : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', 'W (Deflect)', max(exact_W), min(exact_W), mean(exact_W));
fprintf(' %-10s : Max = %+.6e, Min = %+.6e, Avg = %+.6e\n', 'Displacement (Mag)', max(exact_Mag), min(exact_Mag), mean(exact_Mag));
fprintf('==========================================================\n');
% ========================================================================
% w_bar 계산
% ========================================================================
w_bar = max(abs(exact_W))*E2*thick^3/(load*width^4)*100;
fprintf('\nw_bar = %.5f\n', w_bar);
% ========================================================================
% 특정 좌표계의 노드 응력 무차원화 계산
% ========================================================================
a = width;
x_center = a / 2;
y_center = a / 2;
x_edge = a;
y_edge = a;
[~, center_idx] = min((Xall - x_center).^2 + (Yall - y_center).^2);
[~, edge_idx]   = min((Xall - x_edge).^2 + (Yall - y_edge).^2);
sxx_center_val = sxx_top(3, center_idx);
syy_center_val = syy_top(3, center_idx);
sxy_edge_val   = sxy_bot(1, edge_idx);
sxx_bar_center = sxx_center_val * thick^2 / (width^2 * load);
syy_bar_center = syy_center_val * thick^2 / (width^2 * load);
sxy_bar_edge   = sxy_edge_val * thick^2 / (width^2 * load);
% ========================================================================
% 결과 출력
% ========================================================================
fprintf('\n[ Dimensionless Stresses at Specific Coordinates ]\n');
fprintf(' 1. 정 가운데 (a/2, a/2) 노드\n');
fprintf('    - sxx_bar = %.5f\n', sxx_bar_center);
fprintf('    - syy_bar = %.5f\n', syy_bar_center);
fprintf('\n 2. 하단 중앙 (a/2, 0) 노드\n');
fprintf('    - sxy_bar = %.5f\n', sxy_bar_edge);