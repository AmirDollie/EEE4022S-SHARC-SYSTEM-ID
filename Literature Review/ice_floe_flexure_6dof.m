%% ICE-FLOE FLEXURE AND SIX RIGID-BODY MOTIONS
% Self-contained MATLAB figure. No toolboxes or external files required.
% Run this script to save a PNG and PDF in the current directory.
%
% The flexure is deliberately exaggerated for illustration; it is NOT
% the deflection obtained from the hydroelastic EMM.
% The six arrows describe rigid-body kinematic directions (surge, sway,
% heave, roll, pitch, yaw); they do not imply that every small patch has
% six independent degrees of freedom.
%
% Compatible with MATLAB R2018b+ (local functions in scripts).

clearvars; close all; clc;

%% ----------------------- ADJUST THESE -------------------------------
R        = 2.2;         % illustrated floe radius (arbitrary scale)
plateH   = 0.22;        % visible plate thickness
bendAmp  = 0.49;        % centre-to-rim bending contrast; increase to bend more
warpAmp  = 0.05;        % small nonaxisymmetric deformation
viewMain = [-37, 24];   % 3-D camera azimuth and elevation, degrees
viewZoom = [-38, 22];

% Ice/water palette: muted blues rather than saturated primary colours
iceTop   = [0.89 0.97 0.99];
iceEdge  = [0.62 0.83 0.91];
iceLine  = [0.29 0.51 0.62];
water    = [0.48 0.73 0.83];
ink      = [0.16 0.28 0.35];
surgeCol = [0.16 0.52 0.43];  % sea green: x
swayCol  = [0.84 0.49 0.18];  % ochre: y
heaveCol = [0.16 0.43 0.72];  % ocean blue: z

% Illustrative vertical deflection field (quadratic bowl + mild angular term)
% The centre goes downward and rim upward. Not an EMM mode shape.
w = @(x,y) bendAmp*((x.^2+y.^2)/R^2-0.52) ...
           + warpAmp*(x.^2-y.^2)/R^2;

%% --------------------- FIGURE AND LAYOUT ----------------------------
fig = figure('Color','w','Position',[80 60 1450 900]);
set(fig,'Renderer','opengl');

% Soft panel behind the local rigid element
annotation(fig,'rectangle',[0.555 0.565 0.415 0.386], ...
    'Color',[0.57 0.71 0.77],'LineStyle','--','LineWidth',1.15, ...
    'FaceColor',[0.974 0.990 0.996]);
annotation(fig,'textbox',[0.565 0.911 0.39 0.035], ...
    'String','LOCAL RIGID-BODY MOTIONS', ...
    'EdgeColor','none','Color',ink,'FontWeight','bold', ...
    'FontName','Times New Roman','FontSize',14,'HorizontalAlignment','center');
annotation(fig,'textbox',[0.065 0.916 0.41 0.045], ...
    'String','FLEXURAL DEFORMATION OF AN ICE FLOE', ...
    'EdgeColor','none','Color',ink,'FontWeight','bold', ...
    'FontName','Times New Roman','FontSize',16);

%% ---------------------- LARGE FLEXING FLOE -------------------------
ax1 = axes('Parent',fig,'Position',[0.045 0.055 0.90 0.64]);
hold(ax1,'on');

% Annular water surface around the ice disc (so water does not hide the ice)
th = linspace(0,2*pi,241);
rr = linspace(1.015*R,1.47*R,38);
[Tw,Rw] = meshgrid(th,rr);
Xw = Rw.*cos(Tw); Yw = Rw.*sin(Tw);
Zw = -0.40 + 0.025*sin(4*Xw+3*Yw) + 0.018*cos(7*Yw-2*Xw);
surf(ax1,Xw,Yw,Zw,'FaceColor',water,'FaceAlpha',0.46, ...
    'EdgeColor','none');

% Concentric water ripples, subtle and not distracting
for factor = [1.08,1.19,1.32]
    cx = factor*R*cos(th); cy = factor*R*sin(th);
    cz = -0.40 + 0.025*sin(4*cx+3*cy) + 0.018*cos(7*cy-2*cx)+0.006;
    plot3(ax1,cx,cy,cz,'Color',[0.62 0.85 0.91],'LineWidth',0.8);
end

% Dense polar surface mesh for a smooth continuously bent ice plate
[Theta,Rho] = meshgrid(th,linspace(0,R,66));
X = Rho.*cos(Theta); Y = Rho.*sin(Theta);
Z = w(X,Y);
texture = 0.55 + 0.15*sin(7*Theta).*sin(3*Rho/R) ...
               + 0.07*cos(11*Theta+5*Rho/R);
texture = max(0,min(1,texture));
iceMap = [linspace(0.72,iceTop(1),128)', ...
          linspace(0.88,iceTop(2),128)', ...
          linspace(0.95,iceTop(3),128)'];
colormap(ax1,iceMap); caxis(ax1,[0 1]);

% Lower surface first, then the curved side wall and top surface
surf(ax1,X,Y,Z-plateH,'FaceColor',iceEdge,'EdgeColor','none');
xe = R*cos(th); ye = R*sin(th); ze = w(xe,ye);
surf(ax1,[xe;xe],[ye;ye],[ze;ze-plateH], ...
    'FaceColor',iceEdge,'EdgeColor','none');
surf(ax1,X,Y,Z,texture,'FaceColor','interp','EdgeColor','none');
plot3(ax1,xe,ye,ze,'Color',iceLine,'LineWidth',1.6);

% One meridional profile emphasises the curvature. A faint horizontal
% reference line makes its displacement relative to a flat plate obvious.
xx = linspace(-R,R,240); yy = -0.24*R*ones(size(xx));
zz = w(xx,yy);
plot3(ax1,xx,yy,zz+0.022,'--','Color',ink,'LineWidth',1.8);
plot3(ax1,xx,yy,0.015*ones(size(xx)),':', ...
    'Color',[0.47 0.60 0.65],'LineWidth',1.0);

% A highlighted physical patch, linked schematically to the inset
px = -0.27*R + [-0.24 0.24 0.24 -0.24];
py =  0.29*R + [-0.18 -0.18 0.18 0.18];
pz = w(px,py)+0.027;
patch(ax1,px,py,pz,[0.85 0.94 0.95],'FaceAlpha',0.96, ...
    'EdgeColor',ink,'LineWidth',1.6);

% Bent-plate annotations (keep text in figure margin, not on ice)
text(ax1,-2.18,-1.42,0.54,'Flexure (exaggerated)', ...
    'FontName','Times New Roman','FontSize',14,'FontWeight','bold', ...
    'Color',ink,'Clipping','off');
arrow3(ax1,[-1.87,-1.31,0.34],[-1.47,-0.9,w(-1.47,-0.9)+0.07],ink,1.8);
arrow3(ax1,[0,-0.55,0.05],[0,-0.55,w(0,-0.55)-0.13],heaveCol,2.2);

axis(ax1,'equal'); axis(ax1,'off');
xlim(ax1,[-1.50*R 1.50*R]); ylim(ax1,[-1.48*R 1.48*R]);
zlim(ax1,[-0.9 0.9]);
daspect(ax1,[1 1 0.64]);  % accentuates vertical deformation
view(ax1,viewMain);
lighting(ax1,'gouraud'); camlight(ax1,'headlight');

%% ---------------------- ZOOMED RIGID ELEMENT -----------------------
ax2 = axes('Parent',fig,'Position',[0.577 0.585 0.360 0.292], ...
    'Color','none');
hold(ax2,'on');

% Cuboidal representative ice element (the local coordinate frame)
V = [-0.30 -0.26 -0.10; 0.30 -0.26 -0.10; ...
      0.30  0.26 -0.10;-0.30  0.26 -0.10; ...
     -0.30 -0.26  0.18; 0.30 -0.26  0.18; ...
      0.30  0.26  0.18;-0.30  0.26  0.18];
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; ...
         2 3 7 6; 3 4 8 7; 4 1 5 8];
faceColors = [iceEdge;iceTop; 0.68 0.86 0.93; ...
              0.62 0.81 0.89; 0.78 0.92 0.96; 0.66 0.84 0.91];
for f = 1:6
    patch(ax2,'Vertices',V,'Faces',faces(f,:), ...
        'FaceColor',faceColors(f,:),'EdgeColor',iceLine, ...
        'LineWidth',0.95,'FaceAlpha',0.99);
end

% The translation arrows are EXACTLY parallel to +x,+y,+z respectively.
% In particular, sway is +y, NOT -x as in the early illustration.
arrow3(ax2,[0.30,-0.05,0.02],[1.12,-0.05,0.02],surgeCol,2.7);
arrow3(ax2,[-0.08,0.26,0.04],[-0.08,1.13,0.04],swayCol,2.7);
arrow3(ax2,[-0.08,-0.08,0.18],[-0.08,-0.08,1.14],heaveCol,2.7);

% Curved rotational arrows: roll about x, pitch about y, yaw about z.
rotationArrow(ax2,[0.73, 0.00,0.05],'x',0.29,surgeCol);
rotationArrow(ax2,[0.00,-0.77,0.06],'y',0.29,swayCol);
rotationArrow(ax2,[0.00, 0.00,0.48],'z',0.38,heaveCol);

% DOF labels, in matching ice/ocean-related colours
lab(ax2,[1.24,-0.05,0.02],'Surge  (+x)',surgeCol);
lab(ax2,[-0.08,1.29,0.06],'Sway  (+y)',swayCol);
lab(ax2,[-0.08,-0.08,1.33],'Heave  (+z)',heaveCol);
lab(ax2,[0.86,-0.40,-0.32],'Roll  (about x)',surgeCol);
lab(ax2,[-0.34,-1.07,0.23],'Pitch  (about y)',swayCol);
lab(ax2,[0.40,0.38,0.65],'Yaw  (about z)',heaveCol);

% Small coordinate triad: direction and label agree with the DOF arrows
O = [-0.92,-0.85,-0.32];
len = 0.35;
arrow3(ax2,O,O+[len 0 0],ink,1.5);
arrow3(ax2,O,O+[0 len 0],ink,1.5);
arrow3(ax2,O,O+[0 0 len],ink,1.5);
lab(ax2,O+[1.18*len 0 0],'x',ink);
lab(ax2,O+[0 1.18*len 0],'y',ink);
lab(ax2,O+[0 0 1.18*len],'z',ink);

axis(ax2,'equal'); axis(ax2,'off');
xlim(ax2,[-1.48 1.75]); ylim(ax2,[-1.42 1.61]);
zlim(ax2,[-0.55 1.53]);
daspect(ax2,[1 1 0.9]);
view(ax2,viewZoom);

% Schematic connector (edit these two coordinate pairs if repositioning)
annotation(fig,'line',[0.43 0.56],[0.52 0.64], ...
    'LineStyle','--','Color',[0.40 0.57 0.64],'LineWidth',1.25);

annotation(fig,'textbox',[0.51 0.029 0.46 0.026], ...
    'String','Illustrative flexure, not a computed hydroelastic eigenmode', ...
    'EdgeColor','none','HorizontalAlignment','right', ...
    'Color',[0.42 0.53 0.60],'FontName','Times New Roman','FontSize',10);

%% -------------------------- SAVE -----------------------------------
% Native MATLAB .fig is editable; PDF is convenient for LaTeX.
savefig(fig,'ice_floe_flexure_6dof.fig');
if exist('exportgraphics','file')
    exportgraphics(fig,'ice_floe_flexure_6dof.png','Resolution',300);
    exportgraphics(fig,'ice_floe_flexure_6dof.pdf','ContentType','image', ...
        'Resolution',300);
else
    print(fig,'ice_floe_flexure_6dof','-dpng','-r300');
    print(fig,'ice_floe_flexure_6dof','-dpdf','-r300');
end

%% -------------------- LOCAL DRAWING FUNCTIONS ----------------------
function arrow3(ax,p0,p1,col,lw)
    d = p1-p0;
    quiver3(ax,p0(1),p0(2),p0(3),d(1),d(2),d(3),0, ...
        'LineWidth',lw,'Color',col,'MaxHeadSize',0.55, ...
        'AutoScale','off');
end

function rotationArrow(ax,c,which,r,col)
    a = linspace(-0.75*pi,0.67*pi,90);
    switch lower(which)
        case 'x'  % yz plane, axis of rotation = x
            P = [c(1)*ones(size(a)); c(2)+r*cos(a); c(3)+r*sin(a)];
        case 'y'  % xz plane, axis of rotation = y
            P = [c(1)+r*cos(a); c(2)*ones(size(a)); c(3)+r*sin(a)];
        case 'z'  % xy plane, axis of rotation = z
            P = [c(1)+r*cos(a); c(2)+r*sin(a); c(3)*ones(size(a))];
    end
    plot3(ax,P(1,:),P(2,:),P(3,:),'Color',col,'LineWidth',2.3);
    tang = P(:,end)-P(:,end-8);
    tang = 0.18*tang/norm(tang);
    tip = P(:,end);
    arrow3(ax,tip-tang,tip,col,2.4);
end

function lab(ax,xyz,str,col)
    text(ax,xyz(1),xyz(2),xyz(3),str, ...
        'FontName','Times New Roman','FontSize',11, ...
        'FontWeight','bold','Color',col, ...
        'HorizontalAlignment','center','Clipping','off');
end
