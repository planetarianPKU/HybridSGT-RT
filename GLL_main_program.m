% Author: Jingnan Sun
% Affiliation: Peking University
% Date: 2025.11

% ============================================================
% Main Program: Synthesize velocity waveforms using the
%               displacement representation theorem.
%
% Workflow:
%   1. Read SGT binary files from OUT_MESH and IN_MESH
%   2. Use prepare_data() to format SGT required by the
%      representation theorem
%   3. Use get_GLLintegral() to perform GLL integration
%   4. Save results and plot the synthetic velocity waveform
% ============================================================

% -------- Root directory (!!!Change here!!!) --------
curdir = fileparts(mfilename('fullpath'));
dir_root = fileparts(fileparts(curdir));


% -------- List of station names --------
snamelst = ["ROC", "XWE"];

% -------- Domain boundaries (xmin, xmax, ymin, ymax, zmin) --------
Bxmin = 469000;
Bxmax = 514000;
Bymin = 3095000;
Bymax = 3140000;
Bz    = -1e4;     % Depth (negative downward)

Boundarys_pos = [Bxmin Bxmax Bymin Bymax Bz];

% ============================================================
%  Loop structure:
%     iEVENT      : event index for SSGT
%     icomp_sta   : source force component (1=X, 2=Y, 3=Z)
%     isname      : station index for RSGT
% ============================================================

for iEVENT = 1                 % Loop over event index
    for icomp_sta = 1:3        % Loop over source-force components
        for isname = 1:2       % Loop over station names

            sname = char(snamelst(isname)); 

            % ===== OUT_MESH directory (receiver-side SGT binfiles) =====
            directoryPath_OUT_MESH = [dir_root '\SGT_upload\binfiles\MESH_OUT\'];

            % Directory where RSGT files are stored for this station and component
            savedir_now_OUT = [dir_root '\SGT_upload\RSGT\' sname '\dir' char(num2str(icomp_sta)) '\'];
            disp(savedir_now_OUT);

            % ===== IN_MESH directory (source-side SGT binfiles) =====
            directoryPath_IN_MESH = [dir_root '\SGT_upload\binfiles\MESH_IN\'];

            % Directory of SSGT input files for the current event
            savedir_now_IN = [dir_root '\SGT_upload\SSGT\event' char(num2str(iEVENT,'%03d')) '\'];
            disp(directoryPath_IN_MESH);

            % ====================================================
            %     Step 1: Build boundary dataset using SSGT and RSGT
            % ====================================================
            [D] = prepare_data(directoryPath_OUT_MESH, savedir_now_OUT, ...
                         directoryPath_IN_MESH, savedir_now_IN, Boundarys_pos);

            % ====================================================
            %     Step 2: Perform GLL integration to obtain the
            %             synthetic velocity waveform
            % ====================================================
            [RP_result, GLL_result] = get_GLLintegral(D);

            % ====================================================
            %     Step 3: Save integration results
            % ====================================================
            save_result_path = [dir_root '\SGT_upload\OUTPUT\' sname];
            mkdir(save_result_path);    % Create directory if not existing

            save_name = ['/dir_' char(num2str(icomp_sta)) ...
                         '_event' char(num2str(iEVENT,'%03d')) '_data_integral.mat'];

            save_path = [save_result_path save_name];

            % Save only GLL_result (contains the velocity time series)
            save(save_path,'GLL_result','-v7.3');

        end
    end
end

% ============================================================
%     Step 4: Plot the synthesized velocity waveform
% ============================================================

figure;
hold on;

% Load last computed result (save_path from the loop)
data = load(save_path);

% Time vector: dt = 0.05 s, total 3000 samples (0–149.95 s)
t = 0:0.05:2999*0.05-0.05;

% Plot waveform
plot(t, data.GLL_result.total_integral);

xlim([0 70]);           % Show first 70 seconds
xlabel('Time (s)');
title('Synthetic Velocity');
grid on;


%%
% D is your structure from get_RP()

faces = {'pxmin','pxmax','pymin','pymax','pzbottom'};
colors = lines(length(faces));  % automatic distinguishable colors

figure;
hold on;
grid on;
axis equal;

points = D.points2_reordered;   % N × 3 coordinates

for i = 1:length(faces)

    face_name = faces{i};
    face_struct = D.face_results.(face_name);

    face_idx = face_struct.faceidx;  % indices for this surface
    coords = points(face_idx, :);    % extract coordinates

    scatter3(coords(:,1), coords(:,2), coords(:,3), ...
             20, colors(i,:), 'filled');

    % add label in plot
    text(mean(coords(:,1)), mean(coords(:,2)), mean(coords(:,3)), ...
         ['  ' face_name], 'Color', colors(i,:), 'FontSize', 10);
end

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Boundary Face Node Locations');
view(3);


