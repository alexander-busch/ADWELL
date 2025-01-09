ax_bottom = gca; % current axes
ax_bottom_pos = ax_bottom.Position; % position of first axes
ax_top = axes('Position',ax_bottom_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

Z=5;
y_cell = h;
alpha_s_cell = 0.55;

R = y/H;
R_b = h/H;

alpha_s = alpha_s_cell.*( (1-R).*R_b./(1-R_b)./R).^Z;
line(alpha_s,y,'Parent',ax_top)

xlim([0 0.63]);

plot(alpha_s.*u,y)
xlim([0 0.03]);