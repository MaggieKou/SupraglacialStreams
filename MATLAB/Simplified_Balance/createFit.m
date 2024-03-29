function [fitresult, gof] = createFit(t_plot, LHm_plot)
%CREATEFIT(T_PLOT,LHM_PLOT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: t_plot
%      Y Output: LHm_plot
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 27-Jul-2023 10:48:29


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t_plot, LHm_plot );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.372409740055537 0.198118402542975];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'LHm_plot vs. t_plot', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 't_plot', 'Interpreter', 'none' );
ylabel( 'LHm_plot', 'Interpreter', 'none' );
grid on


