% =========================================================================
% @file    calculate_specific_1d_axes.m
% @brief   Calculate anterior-posterior axis and optimal 1D axis for a specific angle pair
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key workflow:
%          1. Calculate anterior-posterior axis (i=51, ii=1: aa=π, bb=0.5π)
%          2. Calculate optimal 1D axis (same angle pair as anterior-posterior axis)
%          3. Compute correlations with lp(:,19) using corr_lp
%          4. Extract final axis values (anterior_posterior_axis/optimal_1d_axis)
% =========================================================================

%% cortest: Row 51 = anterior-posterior axis results, Row 52 = optimal 1D axis results (single angle case)
% Load neurophysiological data (lp(:,19) is target data for correlation analysis)
load('lp99.mat')

% Load 3D coordinate data of 91 brain regions (loc91(:,1-3) = x/y/z axis coordinates respectively)
load('loc91.mat');

% ==================== Calculate Anterior-Posterior Axis ====================
% Fixed angle index i=51 (only this single angle is calculated)
for i =51;
   % Fixed angle index ii=1 (only this single angle is calculated)
   for ii=1;
        % Calculate first rotation angle aa (radians): 
        % aa = (i-1)*0.01π + 0.5π = (51-1)*0.01π + 0.5π = π (180°)
        aa=(i-1)*0.01*pi+0.5*pi;
        clear x1;  % Clear x1 to avoid residual values
        % First rotation: Generate x1 via rotation based on aa (loc91(:,2)/(:,3) = y/z axis coordinates)
        x1=cos(aa)*loc91(:,2)+sin(aa)*loc91(:,3);
        
        % Calculate correlation between x1 and lp(:,19) 
        % (corr_lp = custom correlation function, returns correlation coefficient + p-value)
        [cor_test{i,ii}(1,1),cor_test{i,ii}(2,1)]=corr_lp(x1,lp(:,19));
        cor_test{i,ii}(3,1)=aa;          % Store aa (radians)
        cor_test{i,ii}(4,1)=cos(aa);     % Store cos(aa)
        cor_test{i,ii}(5,1)=sin(aa);     % Store sin(aa)
        cor_test{i,ii}(6,1)=aa/2/pi*360; % Store aa (converted to degrees for interpretability)
        
        % Calculate second rotation angle bb (radians): 
        % bb = (ii-1)*0.01π + 0.5π = 0 + 0.5π = π/2 (90°)
        bb=(ii-1)*0.01*pi+0.5*pi;
        clear x2;  % Clear x2 to avoid residual values
        % Second rotation: Generate x2 via rotation based on bb 
        % (x1 = result of first rotation, loc91(:,1) = x axis coordinate)
        x2=cos(bb)*x1+sin(bb)*loc91(:,1);
        
        % Calculate correlation between x2 and lp(:,19) 
        % (Note: Original code typo - cor_test{i}(2,2) should be cor_test{i,ii}(2,2); kept as original)
        [cor_test{i,ii}(1,2),cor_test{i}(2,2)]=corr_lp(x2,lp(:,19));
        cor_test{i,ii}(3,2)=bb;          % Store bb (radians)
        cor_test{i,ii}(4,2)=cos(bb);     % Store cos(bb)
        cor_test{i,ii}(5,2)=sin(bb);     % Store sin(bb)
        cor_test{i,ii}(6,2)=bb/2/pi*360; % Store bb (converted to degrees for interpretability)
        
        % Store correlation coefficient of x2 vs lp(:,19) in cor matrix
        cor(i,ii)=cor_test{i,ii}(1,2);
   end;
end;

% Extract anterior-posterior axis result (equal to x1 after first rotation)
anterior_posterior_axis=x1;


% ==================== Calculate Optimal 1D Axis ====================
% Fixed angle index i=51 (same angle as anterior-posterior axis)
for i =51;
   % Fixed angle index ii=1 (same angle as anterior-posterior axis)
   for ii=1;
        % Re-calculate first rotation angle aa (identical logic to anterior-posterior axis)
        aa=(i-1)*0.01*pi+0.5*pi;
        clear x1;
        x1=cos(aa)*loc91(:,2)+sin(aa)*loc91(:,3);
        
        % Re-calculate correlation between x1 and lp(:,19)
        [cor_test{i,ii}(1,1),cor_test{i,ii}(2,1)]=corr_lp(x1,lp(:,19));
        cor_test{i,ii}(3,1)=aa;
        cor_test{i,ii}(4,1)=cos(aa);
        cor_test{i,ii}(5,1)=sin(aa);
        cor_test{i,ii}(6,1)=aa/2/pi*360;
        
        % Re-calculate second rotation angle bb (identical logic to anterior-posterior axis)
        bb=(ii-1)*0.01*pi+0.5*pi;
        clear x2;
        x2=cos(bb)*x1+sin(bb)*loc91(:,1);
        
        % Re-calculate correlation between x2 and lp(:,19) (keep original typo)
        [cor_test{i,ii}(1,2),cor_test{i}(2,2)]=corr_lp(x2,lp(:,19));
        cor_test{i,ii}(3,2)=bb;
        cor_test{i,ii}(4,2)=cos(bb);
        cor_test{i,ii}(5,2)=sin(bb);
        cor_test{i,ii}(6,2)=bb/2/pi*360;
        
        % Re-store correlation coefficient of second rotation
        cor(i,ii)=cor_test{i,ii}(1,2);
   end;
end;

% Extract optimal 1D axis result (equal to x1 after first rotation, same as anterior-posterior axis)
optimal_1d_axis=x1;