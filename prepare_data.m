function [results]=prepare_data(directoryPath_OUT_MESH,savedir_now_OUT,directoryPath_IN_MESH,savedir_now_IN,Boundarys_pos)


% =========================================================================
% prepare_data
%
% Purpose:
%   Build all physical and geometric quantities required by the
%   Representation Theorem (RP) on the five boundary surfaces of the
%   computational domain. The routine reads mesh information and
%   Green’s tensors from both OUT_MESH (receiver side)
%   and IN_MESH (source side), extracts boundary nodes, computes material
%   parameters and geometric operators, aligns node ordering, and evaluates
%   the boundary integrals for the RP formulation.
%
% Input:
%   directoryPath_OUT_MESH : directory of OUT_MESH binary files (x,y,z coords)
%   savedir_now_OUT        : directory containing OUT-side SGT .mat files
%   directoryPath_IN_MESH  : directory of IN_MESH binary files (geometry/material)
%   savedir_now_IN         : directory containing IN-side SGT .mat files
%   Boundarys_pos          : [Bxmin Bxmax Bymin Bymax Bz]
%                            domain boundary coordinates
%
% Files read from OUT_MESH:
%   proc*x_boundary_size_is*.bin  → x coordinates on the boundary
%   proc*y_boundary_size_is*.bin  → y coordinates on the boundary
%   proc*z_boundary_size_is*.bin  → z coordinates on the boundary
%
% Files read from IN_MESH (geometry, GLL operators, material):
%   Coordinates:
%       x_boundary_size, y_boundary_size, z_boundary_size
%   Jacobian:
%       jacobianl_boundary_size
%   Derivatives of physical coords w.r.t reference coords (xi, eta, gamma):
%       xixl, xiyl, xizl
%       etaxl, etayl, etazl
%       gammaxl, gammayl, gammazl
%   GLL weight functions:
%       wxgll, wygll, wzgll
%   Material properties:
%       kappa  (bulk modulus)
%       mu     (shear modulus)
%       rho    (density)
%
% SGT .mat files (IN and OUT):
%   For each boundary point, SGT files contain:
%       Ivx, Ivy, Ivz      → displacement Green’s tensor components
%       Iv1x1 ... Iv3x3    → spatial derivatives of SGT w.r.t x,y,z
%   OUT-side files contain the same, denoted as Ov*
%
% Physical meaning of SGT variables:
%   Ivx(t,i)  = G_x(t; boundary_i, source)     → displacement response
%   Iv1x1     = ∂G_x / ∂x
%   Iv2x1     = ∂G_x / ∂y
%   Iv3x1     = ∂G_x / ∂z
%   (similarly for y,z components)
%
%
% Output:
%   results : structure containing
%       results.lambda_all, mu_all, rho_all, kappa_all
%       results.xixl_all, xiyl_all, xizl_all
%       results.etaxl_all, etayl_all, etazl_all
%       results.gammaxl_all, gammayl_all, gammazl_all
%       results.wxgll_all, wygll_all, wzgll_all
%       results.jacob_all
%       results.points1, results.points2_reordered
%       results.face_results   → RP integrals for each surface
%
% =========================================================================



%%

Bxmin=Boundarys_pos(1);
Bxmax=Boundarys_pos(2);
Bymin=Boundarys_pos(3);
Bymax=Boundarys_pos(4);
Bz=Boundarys_pos(5);

%%
X_all_2_check = [];
Y_all_2_check = [];
Z_all_2_check = [];
pall_uni_OUT_all = [];
mat_2d_all = [];


% 获取相关文件
filesx = dir(fullfile(directoryPath_OUT_MESH, 'proc*x_boundary_size_is*.bin'));
filesy = dir(fullfile(directoryPath_OUT_MESH, 'proc*y_boundary_size_is*.bin'));
filesz = dir(fullfile(directoryPath_OUT_MESH, 'proc*z_boundary_size_is*.bin'));

for iname = 1:length(filesx)
    % 处理 x 文件
    fileNamex = filesx(iname).name;
    tmp = split(fileNamex, '_');
    tmp2 = char(tmp(6));
    tmp3 = split(tmp2,'.');
    sizee = str2num(char(tmp3(1)));
    
    filePathx = fullfile(directoryPath_OUT_MESH, fileNamex);
    fileId1 = fopen(filePathx, 'rb');
    datax = fread(fileId1, sizee + 2, 'single');
    datax = datax(2:end-1);
    fclose(fileId1);

    % 处理 y 文件
    fileNamey = filesy(iname).name;
    filePathy = fullfile(directoryPath_OUT_MESH, fileNamey);
    fileId2 = fopen(filePathy, 'rb');
    datay = fread(fileId2, sizee + 2, 'single');
    datay = datay(2:end-1);
    fclose(fileId2);

    % 处理 z 文件
    fileNamez = filesz(iname).name;
    filePathz = fullfile(directoryPath_OUT_MESH, fileNamez);
    fileId3 = fopen(filePathz, 'rb');
    dataz = fread(fileId3, sizee + 2, 'single');
    dataz = dataz(2:end-1);
    fclose(fileId3);

    % 筛选坐标
    pxmin_out = find(abs(datax - Bxmin) == 0);
    pxmax_out = find(abs(datax - Bxmax) == 0);
    pymin_out = find(abs(datay - Bymin) == 0);
    pymax_out = find(abs(datay - Bymax) == 0);
    pzbottom_out = find(abs(dataz - Bz) == 0);

    % 逐步筛除不符合条件的点
    pymax_out = pymax_out(datax(pymax_out) >= Bxmin & datax(pymax_out) <= Bxmax & dataz(pymax_out) >= Bz);
    pymin_out = pymin_out(datax(pymin_out) >= Bxmin & datax(pymin_out) <= Bxmax & dataz(pymin_out) >= Bz);
    pxmin_out = pxmin_out(datay(pxmin_out) >= Bymin & datay(pxmin_out) <= Bymax & dataz(pxmin_out) >= Bz);
    pxmax_out = pxmax_out(datay(pxmax_out) >= Bymin & datay(pxmax_out) <= Bymax & dataz(pxmax_out) >= Bz);
    pzbottom_out = pzbottom_out(datay(pzbottom_out) >= Bymin & datay(pzbottom_out) <= Bymax & datax(pzbottom_out) >= Bxmin & datax(pzbottom_out) <= Bxmax);

    pall_uni_OUT_single = unique([pxmin_out; pxmax_out; pymin_out; pymax_out; pzbottom_out]);
    pall_uni_OUT_all = [pall_uni_OUT_all; pall_uni_OUT_single];

    X_all_2_check = [X_all_2_check; datax(pall_uni_OUT_single)];
    Y_all_2_check = [Y_all_2_check; datay(pall_uni_OUT_single)];
    Z_all_2_check = [Z_all_2_check; dataz(pall_uni_OUT_single)];
    
    datax_save = datax(pall_uni_OUT_single);
    datay_save = datay(pall_uni_OUT_single);
    dataz_save = dataz(pall_uni_OUT_single);
    
end

%%
X_all_IN_check = [];
Y_all_IN_check = [];
Z_all_IN_check = [];
pall_uni_IN_all = [];
mat_2d_IN_all = [];

jacob_all_IN_check = [];
xixl_all_IN_check = [];
xiyl_all_IN_check = [];
xizl_all_IN_check = [];

etaxl_all_IN_check = [];
etayl_all_IN_check = [];
etazl_all_IN_check = [];

gammaxl_all_IN_check = [];
gammayl_all_IN_check = [];
gammazl_all_IN_check = [];

wxgll_all_IN_check = [];
wygll_all_IN_check = [];
wzgll_all_IN_check = [];

kappa_all_IN_check = [];
mu_all_IN_check = [];
rho_all_IN_check = [];

% 获取相关文件
filesx = dir(fullfile(directoryPath_IN_MESH, 'proc*x_boundary_size_is*.bin'));
filesy = dir(fullfile(directoryPath_IN_MESH, 'proc*y_boundary_size_is*.bin'));
filesz = dir(fullfile(directoryPath_IN_MESH, 'proc*z_boundary_size_is*.bin'));

files_jacob = dir(fullfile(directoryPath_IN_MESH, 'proc*jacobianl_boundary_size_is*.bin'));

files_xixl = dir(fullfile(directoryPath_IN_MESH, 'proc*xixl_boundary_size_is*.bin'));
files_xiyl = dir(fullfile(directoryPath_IN_MESH, 'proc*xiyl_boundary_size_is*.bin'));
files_xizl = dir(fullfile(directoryPath_IN_MESH, 'proc*xizl_boundary_size_is*.bin'));

files_etaxl = dir(fullfile(directoryPath_IN_MESH, 'proc*etaxl_boundary_size_is*.bin'));
files_etayl = dir(fullfile(directoryPath_IN_MESH, 'proc*etayl_boundary_size_is*.bin'));
files_etazl = dir(fullfile(directoryPath_IN_MESH, 'proc*etazl_boundary_size_is*.bin'));

files_gammaxl = dir(fullfile(directoryPath_IN_MESH, 'proc*gammaxl_boundary_size_is*.bin'));
files_gammayl = dir(fullfile(directoryPath_IN_MESH, 'proc*gammayl_boundary_size_is*.bin'));
files_gammazl = dir(fullfile(directoryPath_IN_MESH, 'proc*gammazl_boundary_size_is*.bin'));

files_wxgll = dir(fullfile(directoryPath_IN_MESH, 'proc*wxgll_boundary_size_is*.bin'));
files_wygll = dir(fullfile(directoryPath_IN_MESH, 'proc*wygll_boundary_size_is*.bin'));
files_wzgll = dir(fullfile(directoryPath_IN_MESH, 'proc*wzgll_boundary_size_is*.bin'));


files_kappa = dir(fullfile(directoryPath_IN_MESH, 'proc*kappa_boundary_size_is*.bin'));
files_mu = dir(fullfile(directoryPath_IN_MESH, 'proc*mu_boundary_size_is*.bin'));
files_rho = dir(fullfile(directoryPath_IN_MESH, 'proc*rho_boundary_size_is*.bin'));




for iname = 1:length(filesx)
    % 处理 x 文件
    fileNamex = filesx(iname).name;
    tmp = split(fileNamex, '_');
    tmp2 = char(tmp(6));
    tmp3 = split(tmp2,'.');
    sizee = str2num(char(tmp3(1)));
    
    filePathx = fullfile(directoryPath_IN_MESH, fileNamex);
    fileId1 = fopen(filePathx, 'rb');
    datax = fread(fileId1, sizee + 2, 'single');
    datax = datax(2:end-1);
    fclose(fileId1);
        
    % 处理 y 文件
    fileNamey = filesy(iname).name;
    filePathy = fullfile(directoryPath_IN_MESH, fileNamey);
    fileId2 = fopen(filePathy, 'rb');
    datay = fread(fileId2, sizee + 2, 'single');
    datay = datay(2:end-1);
    fclose(fileId2);

    % 处理 z 文件
    fileNamez = filesz(iname).name;
    filePathz = fullfile(directoryPath_IN_MESH, fileNamez);
    fileId3 = fopen(filePathz, 'rb');
    dataz = fread(fileId3, sizee + 2, 'single');
    dataz = dataz(2:end-1);
    fclose(fileId3);

    % 处理 jacob 文件
    fileName = files_jacob(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_jacob = fread(fileId, sizee + 2, 'single');
    data_jacob = data_jacob(2:end-1);
    fclose(fileId);

    % 处理 偏导数 文件
    %xi
    fileName = files_xixl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_xixl = fread(fileId, sizee + 2, 'single');
    data_xixl = data_xixl(2:end-1);
    fclose(fileId);

    fileName = files_xiyl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_xiyl = fread(fileId, sizee + 2, 'single');
    data_xiyl = data_xiyl(2:end-1);
    fclose(fileId);

    fileName = files_xizl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_xizl = fread(fileId, sizee + 2, 'single');
    data_xizl = data_xizl(2:end-1);
    fclose(fileId);
    %eta
    fileName = files_etaxl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_etaxl = fread(fileId, sizee + 2, 'single');
    data_etaxl = data_etaxl(2:end-1);
    fclose(fileId);

    fileName = files_etayl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_etayl = fread(fileId, sizee + 2, 'single');
    data_etayl = data_etayl(2:end-1);
    fclose(fileId);

    fileName = files_etazl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_etazl = fread(fileId, sizee + 2, 'single');
    data_etazl = data_etazl(2:end-1);
    fclose(fileId);
    %gamma
    fileName = files_gammaxl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_gammaxl = fread(fileId, sizee + 2, 'single');
    data_gammaxl = data_gammaxl(2:end-1);
    fclose(fileId);
    fileName = files_gammayl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_gammayl = fread(fileId, sizee + 2, 'single');
    data_gammayl = data_gammayl(2:end-1);
    fclose(fileId);
    fileName = files_gammazl(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_gammazl = fread(fileId, sizee + 2, 'single');
    data_gammazl = data_gammazl(2:end-1);
    fclose(fileId);

    %wxgll
    fileName = files_wxgll(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_wxgll = fread(fileId, sizee + 2, 'single');
    data_wxgll = data_wxgll(2:end-1);
    fclose(fileId);
    %wygll
    fileName = files_wygll(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_wygll = fread(fileId, sizee + 2, 'single');
    data_wygll = data_wygll(2:end-1);
    fclose(fileId);
    %wzgll
    fileName = files_wzgll(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_wzgll = fread(fileId, sizee + 2, 'single');
    data_wzgll = data_wzgll(2:end-1);
    fclose(fileId);    

    %kappa
    fileName = files_kappa(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_kappa = fread(fileId, sizee + 2, 'single');
    data_kappa = data_kappa(2:end-1);
    fclose(fileId);       
    %mu
    fileName = files_mu(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_mu = fread(fileId, sizee + 2, 'single');
    data_mu = data_mu(2:end-1);
    fclose(fileId);  
    %rho
    fileName = files_rho(iname).name;
    filePath = fullfile(directoryPath_IN_MESH, fileName);
    fileId = fopen(filePath, 'rb');
    data_rho = fread(fileId, sizee + 2, 'single');
    data_rho = data_rho(2:end-1);
    fclose(fileId);     



    % 筛选坐标
    pxmin_IN = find(abs(datax - Bxmin) == 0);
    pxmax_IN = find(abs(datax - Bxmax) == 0);
    pymin_IN = find(abs(datay - Bymin) == 0);
    pymax_IN = find(abs(datay - Bymax) == 0);
    pzbottom_IN = find(abs(dataz - Bz) == 0);

    % 逐步筛除不符合条件的点
    pymax_IN = pymax_IN(datax(pymax_IN) >= Bxmin & datax(pymax_IN) <= Bxmax & dataz(pymax_IN) >= Bz);
    pymin_IN = pymin_IN(datax(pymin_IN) >= Bxmin & datax(pymin_IN) <= Bxmax & dataz(pymin_IN) >= Bz);
    pxmin_IN = pxmin_IN(datay(pxmin_IN) >= Bymin & datay(pxmin_IN) <= Bymax & dataz(pxmin_IN) >= Bz);
    pxmax_IN = pxmax_IN(datay(pxmax_IN) >= Bymin & datay(pxmax_IN) <= Bymax & dataz(pxmax_IN) >= Bz);
    pzbottom_IN = pzbottom_IN(datay(pzbottom_IN) >= Bymin & datay(pzbottom_IN) <= Bymax & datax(pzbottom_IN) >= Bxmin & datax(pzbottom_IN) <= Bxmax);

    pall_uni_IN_single = unique([pxmin_IN; pxmax_IN; pymin_IN; pymax_IN; pzbottom_IN]);
    pall_uni_IN_all = [pall_uni_IN_all; pall_uni_IN_single];

    X_all_IN_check = [X_all_IN_check; datax(pall_uni_IN_single)];
    Y_all_IN_check = [Y_all_IN_check; datay(pall_uni_IN_single)];
    Z_all_IN_check = [Z_all_IN_check; dataz(pall_uni_IN_single)];
    
jacob_all_IN_check = [jacob_all_IN_check; data_jacob(pall_uni_IN_single)];
xixl_all_IN_check = [xixl_all_IN_check;data_xixl(pall_uni_IN_single)];
xiyl_all_IN_check = [xiyl_all_IN_check;data_xiyl(pall_uni_IN_single)];
xizl_all_IN_check = [xizl_all_IN_check;data_xizl(pall_uni_IN_single)];

etaxl_all_IN_check = [etaxl_all_IN_check;data_etaxl(pall_uni_IN_single)];
etayl_all_IN_check = [etayl_all_IN_check;data_etayl(pall_uni_IN_single)];
etazl_all_IN_check = [etazl_all_IN_check;data_etazl(pall_uni_IN_single)];

gammaxl_all_IN_check = [gammaxl_all_IN_check;data_gammaxl(pall_uni_IN_single)];
gammayl_all_IN_check = [gammayl_all_IN_check;data_gammayl(pall_uni_IN_single)];
gammazl_all_IN_check = [gammazl_all_IN_check;data_gammazl(pall_uni_IN_single)];

wxgll_all_IN_check = [wxgll_all_IN_check;data_wxgll(pall_uni_IN_single)];
wygll_all_IN_check = [wygll_all_IN_check;data_wygll(pall_uni_IN_single)];
wzgll_all_IN_check = [wzgll_all_IN_check;data_wzgll(pall_uni_IN_single)];

kappa_all_IN_check = [kappa_all_IN_check;data_kappa(pall_uni_IN_single)];
mu_all_IN_check = [mu_all_IN_check;data_mu(pall_uni_IN_single)];
rho_all_IN_check = [rho_all_IN_check;data_rho(pall_uni_IN_single)];


    datax_save = datax(pall_uni_IN_single);
    datay_save = datay(pall_uni_IN_single);
    dataz_save = dataz(pall_uni_IN_single);
    

end



kappa_all=kappa_all_IN_check;
mu_all=mu_all_IN_check;
rho_all=rho_all_IN_check;

lambda_all = kappa_all - 2/3 * mu_all;


xixl_all = xixl_all_IN_check;
xiyl_all = xiyl_all_IN_check;
xizl_all = xizl_all_IN_check;

etaxl_all = etaxl_all_IN_check;
etayl_all = etayl_all_IN_check;
etazl_all = etazl_all_IN_check;

gammaxl_all = gammaxl_all_IN_check;
gammayl_all = gammayl_all_IN_check;
gammazl_all = gammazl_all_IN_check;

wxgll_all = wxgll_all_IN_check;
wygll_all = wygll_all_IN_check;
wzgll_all = wzgll_all_IN_check;
jacob_all = jacob_all_IN_check;


% save_result_path = [dir_root '/MESH_OUT/RP_result/' sname ];
% mkdir(save_result_path);
% 
% 
% save_name = ['/dir_' char(num2str(icomp_sta)) '_event' char(num2str(iEVENT,'%03d')) '_data_5surfaces_Jacob_correct_0.mat'];
% save_path = [save_result_path save_name];
% %save(save_path,'Bxixl_all','Bxiyl_all','Bxizl_all','Betaxl_all','Betayl_all','Betazl_all','Bgammaxl_all','Bgammayl_all','Bgammazl_all','Bwxgll_all','Bwygll_all','Bwzgll_all','Bjacob_all','-v7.3');
% save(save_path,'lambda_all','mu_all','rho_all','kappa_all','xixl_all','xiyl_all','xizl_all','etaxl_all','etayl_all','etazl_all','gammaxl_all','gammayl_all','gammazl_all','wxgll_all','wygll_all','wzgll_all','jacob_all','X_all_IN_check','Y_all_IN_check','Z_all_IN_check','-v7.3');



%%
disp('check isequal:')
isequal(X_all_2_check,X_all_IN_check)

% 假设 points1 是一个 N×3 的矩阵
points1 = [X_all_2_check, Y_all_2_check, Z_all_2_check]; % 第一组三维点

% 找到 points1 中的唯一点和每个点第一次出现的索引
[uniquePoints, ~, idxUnique] = unique(points1, 'rows', 'stable');
[~, selectIdx] = ismember(uniquePoints, points1, 'rows');

points1 = points1(selectIdx,:);

% 假设 points1 和 points2 是两个 N×3 的矩阵
points2 = [X_all_IN_check, Y_all_IN_check, Z_all_IN_check]; % 第二组三维点

% 使用 ismember 函数找到 points1 中每个点在 points2 中的索引
[~, reorderIdx] = ismember(points1, points2, 'rows');

% 检查是否有未找到匹配的点
if any(reorderIdx == 0)
    warning('有一些点在 points2 中没有找到匹配。');
end

% 使用排序索引重新排列 points2
points2_reordered = points2(reorderIdx, :);

% 检查 points2_reordered 是否与 points1 相同
assert(isequal(points1, points2_reordered), 'Reordered points2 does not match points1.');

if isequal(points1, points2_reordered)
    disp('Reordered points2 matches points1.');
else
    assert(isequal(points1, points2_reordered), 'Reordered points2 does not match points1.');
end


size(points1)
size(points2_reordered)

dist_now = sqrt( (points1(:,1)-points2_reordered(:,1)).^2  + (points1(:,2)-points2_reordered(:,2)).^2  + (points1(:,3)-points2_reordered(:,3)).^2);

kappa_all = kappa_all(reorderIdx);
rho_all = rho_all(reorderIdx);

lambda_all = lambda_all(reorderIdx);
mu_all = mu_all(reorderIdx);

xixl_all = xixl_all(reorderIdx);
xiyl_all = xiyl_all(reorderIdx);
xizl_all = xizl_all(reorderIdx);

etaxl_all = etaxl_all(reorderIdx);
etayl_all = etayl_all(reorderIdx);
etazl_all = etazl_all(reorderIdx);

gammaxl_all = gammaxl_all(reorderIdx);
gammayl_all = gammayl_all(reorderIdx);
gammazl_all = gammazl_all(reorderIdx);

wxgll_all = wxgll_all(reorderIdx);
wygll_all = wygll_all(reorderIdx);
wzgll_all = wzgll_all(reorderIdx);
jacob_all = jacob_all(reorderIdx);
%%

disp('load SGT mat')
files_mat = dir(fullfile(savedir_now_IN, sprintf('proc*_SGT.mat')));
Iux = [];
Iuy = [];
Iuz = [];
Iu1x1 = [];
Iu1x2 = [];
Iu1x3 = [];
Iu2x1 = [];
Iu2x2 = [];
Iu2x3 = [];
Iu3x1 = [];
Iu3x2 = [];
Iu3x3 = [];
for ifile = 1:length(files_mat)
    fileName = files_mat(ifile).name;
    filePath = fullfile(savedir_now_IN, fileName);
    d = load(filePath);
    Iux = [Iux d.Ivx];
    Iuy = [Iuy d.Ivy];
    Iuz = [Iuz d.Ivz];
    Iu1x1 = [Iu1x1 d.Iv1x1];
    Iu1x2 = [Iu1x2 d.Iv1x2];
    Iu1x3 = [Iu1x3 d.Iv1x3];
    Iu2x1 = [Iu2x1 d.Iv2x1];
    Iu2x2 = [Iu2x2 d.Iv2x2];
    Iu2x3 = [Iu2x3 d.Iv2x3];
    Iu3x1 = [Iu3x1 d.Iv3x1];
    Iu3x2 = [Iu3x2 d.Iv3x2];
    Iu3x3 = [Iu3x3 d.Iv3x3];
end 
Iux = Iux(:,reorderIdx);
Iuy = Iuy(:,reorderIdx);
Iuz = Iuz(:,reorderIdx);

Iu1x1 = Iu1x1(:,reorderIdx);
Iu1x2 = Iu1x2(:,reorderIdx);
Iu1x3 = Iu1x3(:,reorderIdx);

Iu2x1 = Iu2x1(:,reorderIdx);
Iu2x2 = Iu2x2(:,reorderIdx);
Iu2x3 = Iu2x3(:,reorderIdx);

Iu3x1 = Iu3x1(:,reorderIdx);
Iu3x2 = Iu3x2(:,reorderIdx);
Iu3x3 = Iu3x3(:,reorderIdx);
files_mat = dir(fullfile(savedir_now_OUT, sprintf('proc*_SGT.mat')));
Oux = [];
Ouy = [];
Ouz = [];
Ou1x1 = [];
Ou1x2 = [];
Ou1x3 = [];
Ou2x1 = [];
Ou2x2 = [];
Ou2x3 = [];
Ou3x1 = [];
Ou3x2 = [];
Ou3x3 = [];
for ifile = 1:length(files_mat)
    fileName = files_mat(ifile).name;
    filePath = fullfile(savedir_now_OUT, fileName);
    d = load(filePath);

    Oux = [Oux d.Ovx];
    Ouy = [Ouy d.Ovy];
    Ouz = [Ouz d.Ovz];
    Ou1x1 = [Ou1x1 d.Ov1x1];
    Ou1x2 = [Ou1x2 d.Ov1x2];
    Ou1x3 = [Ou1x3 d.Ov1x3];
    Ou2x1 = [Ou2x1 d.Ov2x1];
    Ou2x2 = [Ou2x2 d.Ov2x2];
    Ou2x3 = [Ou2x3 d.Ov2x3];
    Ou3x1 = [Ou3x1 d.Ov3x1];
    Ou3x2 = [Ou3x2 d.Ov3x2];
    Ou3x3 = [Ou3x3 d.Ov3x3];
end 
Oux = Oux(:,selectIdx);
Ouy = Ouy(:,selectIdx);
Ouz = Ouz(:,selectIdx);

Ou1x1 = Ou1x1(:,selectIdx);
Ou1x2 = Ou1x2(:,selectIdx);
Ou1x3 = Ou1x3(:,selectIdx);

Ou2x1 = Ou2x1(:,selectIdx);
Ou2x2 = Ou2x2(:,selectIdx);
Ou2x3 = Ou2x3(:,selectIdx);

Ou3x1 = Ou3x1(:,selectIdx);
Ou3x2 = Ou3x2(:,selectIdx);
Ou3x3 = Ou3x3(:,selectIdx);


Oux = cumsum(Oux);
Ouy = cumsum(Ouy);
Ouz = cumsum(Ouz);

Ou1x1 = cumsum(Ou1x1);
Ou1x2 = cumsum(Ou1x2);
Ou1x3 = cumsum(Ou1x3);

Ou2x1 = cumsum(Ou2x1);
Ou2x2 = cumsum(Ou2x2);
Ou2x3 = cumsum(Ou2x3);

Ou3x1 = cumsum(Ou3x1);
Ou3x2 = cumsum(Ou3x2);
Ou3x3 = cumsum(Ou3x3);

disp(size(Ou3x3))
disp(size(Iu3x3))

NP = size(Ou3x3,2);
disp(NP);



    pxmin_out = find(abs(points1(:,1) - Bxmin) == 0);
    pxmax_out = find(abs(points1(:,1) - Bxmax) == 0);
    pymin_out = find(abs(points1(:,2) - Bymin) == 0);
    pymax_out = find(abs(points1(:,2) - Bymax) == 0);
    pzbottom_out = find(abs(points1(:,3) - Bz) == 0);
    
    disp(size(pxmin_out))
    disp(size(pxmax_out))
    disp(size(pymin_out))
    disp(size(pymax_out))
    disp(size(pzbottom_out))

% 给点编号
num_points1_total = size(points1, 1);
point1_indices = (1:num_points1_total)';
% 将点的编号加入到点数据中
points1_with_idx = [points1, point1_indices];    
Bcoord = points1;
Blambda = lambda_all;
Bmu = mu_all;



Bxixl_all = xixl_all(:);
Bxiyl_all = xiyl_all(:);
Bxizl_all = xizl_all(:);

Betaxl_all = etaxl_all(:);
Betayl_all = etayl_all(:);
Betazl_all = etazl_all(:);

Bgammaxl_all = gammaxl_all(:);
Bgammayl_all = gammayl_all(:);
Bgammazl_all = gammazl_all(:);

Bwxgll_all = wxgll_all(:);
Bwygll_all = wygll_all(:);
Bwzgll_all = wzgll_all(:);
Bjacob_all = jacob_all(:);



% 定义五个面的索引数组
face_names = {'pxmin', 'pxmax', 'pymin', 'pymax', 'pzbottom'};
face_indices = {pxmin_out, pxmax_out, pymin_out, pymax_out, pzbottom_out};

% 结构体存储计算结果
face_results = struct();
    xp = Bcoord(:,1);
    yp = Bcoord(:,2);
    zp = Bcoord(:,3); 

    
figure
hold on;
for f = 1:length(face_names)
    % 当前面名称和索引
    face_name = face_names{f};
    face_out = face_indices{f};
    
    % 获取当前面的数据
    NP_B = size(Ouz(:, face_out), 2);  % 计算当前面上的点数
    Nlst1 = zeros(NP_B, 3);
    
        % 设置法向矢量
    switch face_name
        case 'pxmin', Nlst1(:, 1) = -1;
        case 'pxmax', Nlst1(:, 1) = 1;
        case 'pymin', Nlst1(:, 2) = -1;
        case 'pymax', Nlst1(:, 2) = 1;
        case 'pzbottom', Nlst1(:, 3) = -1;
    end
    
 

% 预分配数组
U_RP_all = zeros(2999,NP_B);
U_RP_term2_all = zeros(2999,NP_B);
U_RP_term3_all = zeros(2999,NP_B);

U_RP=0;
for idxxx = 1:NP_B
    %plot3(xp(face_out),yp(face_out),zp(face_out),'o')
    idx=face_out(idxxx);
     if mod(idxxx,1000)==0
    end
    %坐标
  
    %lambda & mu
    lambda = Blambda(idx);
    mu = Bmu(idx);
    nvec = Nlst1(idxxx,:);


    %法向量，向外大区域为正
    n1 = nvec(1);
    n2 = nvec(2);
    n3 = nvec(3);
    %内部传到边界
     ux   = Iux(1:1:end,idx);
     uy   = Iuy(1:1:end,idx);
     uz   = Iuz(1:1:end,idx);
     u1x1 = Iu1x1(1:1:end,idx);
     u1x2 = Iu1x2(1:1:end,idx);
     u1x3 = Iu1x3(1:1:end,idx);
     u2x1 = Iu2x1(1:1:end,idx);
     u2x2 = Iu2x2(1:1:end,idx);
     u2x3 = Iu2x3(1:1:end,idx);
     u3x1 = Iu3x1(1:1:end,idx);
     u3x2 = Iu3x2(1:1:end,idx);
     u3x3 = Iu3x3(1:1:end,idx);

%     %外部传到边界
     Gux   = Oux(1:1:end,idx);
     Guy   = Ouy(1:1:end,idx);
     Guz   = Ouz(1:1:end,idx);
     Gu1x1 = Ou1x1(1:1:end,idx);
     Gu1x2 = Ou1x2(1:1:end,idx);
     Gu1x3 = Ou1x3(1:1:end,idx);   
     Gu2x1 = Ou2x1(1:1:end,idx);
     Gu2x2 = Ou2x2(1:1:end,idx);
     Gu2x3 = Ou2x3(1:1:end,idx);
     Gu3x1 = Ou3x1(1:1:end,idx);
     Gu3x2 = Ou3x2(1:1:end,idx);
     Gu3x3 = Ou3x3(1:1:end,idx); 
    %for Term 2
    item11 = lambda*n1*(u1x1+u2x2+u3x3);
    item12 = mu*( n1*(u1x1+u1x1) + n2*(u2x1+u1x2) + n3*(u3x1+u1x3) );
    
    item21 = lambda*n2*(u1x1+u2x2+u3x3);
    item22 = mu*( n1*(u1x2+u2x1) + n2*(u2x2+u2x2) + n3*(u3x2+u2x3) );
    
    item31 = lambda*n3*(u1x1+u2x2+u3x3);
    item32 = mu*( n1*(u1x3+u3x1) + n2*(u2x3+u3x2) + n3*(u3x3+u3x3) );
    
    Uterm2 = conv(Gux,item11+item12)+conv(Guy,item21+item22)+conv(Guz,item31+item32);
    
    %for Term 3
    item11 = lambda*n1*(Gu1x1+Gu2x2+Gu3x3);
    item12 = mu*( n1*(Gu1x1+Gu1x1) + n2*(Gu2x1+Gu1x2) + n3*(Gu3x1+Gu1x3) );
    item21 = lambda*n2*(Gu1x1+Gu2x2+Gu3x3);
    item22 = mu*( n1*(Gu1x2+Gu2x1) + n2*(Gu2x2+Gu2x2) + n3*(Gu3x2+Gu2x3) );
    item31 = lambda*n3*(Gu1x1+Gu2x2+Gu3x3);
    item32 = mu*( n1*(Gu1x3+Gu3x1) + n2*(Gu2x3+Gu3x2) + n3*(Gu3x3+Gu3x3) );

    Uterm3 = conv(ux,item11+item12)+conv(uy,item21+item22)+conv(uz,item31+item32);
    
    U = Uterm2 - Uterm3;
    U_RP = U_RP + U;
    U_RP_all(:,idxxx) = U;
    U_RP_term2_all(:,idxxx) = Uterm2;
    U_RP_term3_all(:,idxxx) = Uterm3;

    
    
end

    face_results.(face_name).U_RP_all = U_RP_all;
    face_results.(face_name).U_RP_term2_all = U_RP_term2_all;
    face_results.(face_name).U_RP_term3_all = U_RP_term3_all;
    face_results.(face_name).Nlst1 = Nlst1;
    face_results.(face_name).faceidx=face_out;
    
end

%%

    % 一堆变量的计算
    results.lambda_all = lambda_all;
    results.mu_all     = mu_all;
    results.rho_all    = rho_all;
    results.kappa_all  = kappa_all;
    results.xixl_all   = xixl_all;
    results.xiyl_all   = xiyl_all;
    results.xizl_all   = xizl_all;
    results.etaxl_all  = etaxl_all;
    results.etayl_all  = etayl_all;
    results.etazl_all  = etazl_all;
    results.gammaxl_all = gammaxl_all;
    results.gammayl_all = gammayl_all;
    results.gammazl_all = gammazl_all;
    results.wxgll_all  = wxgll_all;
    results.wygll_all  = wygll_all;
    results.wzgll_all  = wzgll_all;
    results.jacob_all  = jacob_all;
    results.points1    = points1;
    results.points2_reordered = points2_reordered;
    results.face_results = face_results;


disp('save finish')


%% Save extracted datasets (split into META / IN_SGT / OUT_SGT)
disp('Saving extracted datasets...');

%% ------------------------------- 
% META dataset (coordinates + materials + geometry)
% -------------------------------
META = struct();

% Coordinates
META.points_boundary = points1;
META.points_IN_reordered = points2_reordered;

% Material parameters
META.lambda = lambda_all;
META.mu = mu_all;
META.rho = rho_all;
META.kappa = kappa_all;

% Geometry derivatives
META.xixl = xixl_all;
META.xiyl = xiyl_all;
META.xizl = xizl_all;

META.etaxl = etaxl_all;
META.etayl = etayl_all;
META.etazl = etazl_all;

META.gammaxl = gammaxl_all;
META.gammayl = gammayl_all;
META.gammazl = gammazl_all;

% GLL weights
META.wxgll = wxgll_all;
META.wygll = wygll_all;
META.wzgll = wzgll_all;

% Jacobian
META.jacobian = jacob_all;

save(fullfile(savedir_now_IN, 'META_Dataset.mat'), ...
     'META', '-v7.3');
disp('Saved META_Dataset.mat');


%% ------------------------------- 
% IN-side SGT dataset (source → boundary)
% -------------------------------
IN_SGT = struct();

IN_SGT.Iux = Iux;
IN_SGT.Iuy = Iuy;
IN_SGT.Iuz = Iuz;

IN_SGT.Iu1x1 = Iu1x1;   IN_SGT.Iu1x2 = Iu1x2;   IN_SGT.Iu1x3 = Iu1x3;
IN_SGT.Iu2x1 = Iu2x1;   IN_SGT.Iu2x2 = Iu2x2;   IN_SGT.Iu2x3 = Iu2x3;
IN_SGT.Iu3x1 = Iu3x1;   IN_SGT.Iu3x2 = Iu3x2;   IN_SGT.Iu3x3 = Iu3x3;

save(fullfile(savedir_now_IN, 'IN_SGT_Dataset.mat'), ...
     'IN_SGT','-v7.3');
disp('Saved IN_SGT_Dataset.mat');


%% ------------------------------- 
% OUT-side SGT dataset (boundary → receiver)
% -------------------------------
OUT_SGT = struct();

OUT_SGT.Oux = Oux;
OUT_SGT.Ouy = Ouy;
OUT_SGT.Ouz = Ouz;

OUT_SGT.Ou1x1 = Ou1x1;   OUT_SGT.Ou1x2 = Ou1x2;   OUT_SGT.Ou1x3 = Ou1x3;
OUT_SGT.Ou2x1 = Ou2x1;   OUT_SGT.Ou2x2 = Ou2x2;   OUT_SGT.Ou2x3 = Ou2x3;
OUT_SGT.Ou3x1 = Ou3x1;   OUT_SGT.Ou3x2 = Ou3x2;   OUT_SGT.Ou3x3 = Ou3x3;

save(fullfile(savedir_now_OUT, 'OUT_SGT_Dataset.mat'), ...
     'OUT_SGT','-v7.3');
disp('Saved OUT_SGT_Dataset.mat');

disp('All datasets saved successfully.');


end
