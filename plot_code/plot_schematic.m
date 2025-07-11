% Bingchen Liu July 2, 2025
% This code plot the schematic for funwaveC simulation including wave maker
% and sponge

%% fig setup 
close all 
xsize=12;ysize=6;
x0=0.15; 
y0=0.25;  
dx=0.12;
xw=0.75; yw=0.6;
x1= x0+xw+dx; 
pos = [x0 y0 xw yw];

fig_fontsz= 12;
text_fontsz= 12;
subfiglabel_fontsz = 15;
cbar_fontsize =10;
cbar_ticksz= 6;


%%
%--- PARAMETERS ----------------------------------------------------------
h0         = 9;            % flat depth [m]
x_flat1    = -550;           % start of flat region [m]
x_slope_st = -225;           % where slope begins [m]
x_wm       = [-420 -400];    % wavemaker region [m]
gap        = 20;             % white gap width between sponge & wavemaker [m]
gap_left =0;

% y‐limits for patches
y_max = 0;
y_min = -h0;

% compute sponge shading extents (left of the gap)
%x_sponge = [x_flat1+gap_left, x_wm(1)-gap];
x_sponge = [x_flat1+gap_left,-450];

%--- PLOTTING ------------------------------------------------------------
%figure('Color','w')
figure()
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(1,:));
hold on

% 1) sponge patch (dark grey)
patch( [x_sponge, fliplr(x_sponge)], ...
       [y_max, y_max, y_min, y_min], ...
       [0.6 0.6 0.6], 'EdgeColor','none' )

% 2) wavemaker patch (light grey)
patch( [x_wm, fliplr(x_wm)], ...
       [y_max, y_max, y_min, y_min], ...
       [0.85 0.85 0.85], 'EdgeColor','none' )

% 3) water surface
plot(linspace(x_flat1+gap_left,0,200), zeros(1,200), 'k', 'LineWidth', 2)

% 4) bathymetry: flat + slope
plot(linspace(x_flat1+gap_left, x_slope_st, 100), -h0*ones(1,100), 'k', 'LineWidth', 2)
plot(linspace(x_slope_st, 24.25,100), linspace(-h0,0.97,100),   'k', 'LineWidth', 2)

% 5) labels with corrected property names
text(mean(x_sponge),  0.1,    'Sponge',    ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'FontSize',text_fontsz)

text(-392,      -4.5,    'Wavemaker', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'Rotation',90, ...
     'FontSize',text_fontsz)

% …
% 6) axes formatting
ax = gca;
ax.FontSize         = fig_fontsz;
ax.LineWidth        = 1;
ax.TickDir          = 'out';

% draw a full box and show ticks on top & right too
ax.Box              = 'on';

ax.XLim             = [-600, 50];
ax.YLim             = [-10, 2];

% ticks every 100 from -500 to 0
ax.XTick            =[-600,-500,-400,-300,-200,-100,0,50]; %x_flat1:100:50;
ax.YTick            = [-10,-8,-6,-4,-2,0,2];
ax.XMinorTick       = 'on';
ax.YMinorTick       = 'on';
%ax.Label.Rotation = 0;
xtickangle(0);

grid on
ax.GridLineStyle      = ':';
%ax.MinorGridLineStyle = ':';
ax.GridAlpha          = 0.2;
ax.MinorGridAlpha     = 0.15;

xlabel('$x$ (m)','Interpreter','latex','FontSize',fig_fontsz)
ylabel('$z$ (m)','Interpreter','latex','FontSize',fig_fontsz)

%xlim([600,25])
hold off 

print -dpdf -painters '/data1/bliu/vortpaper/fig_schematic.pdf'