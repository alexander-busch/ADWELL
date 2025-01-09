% Plot Re_pipe vs. Sh_pipe

figure(fig_Re_Sh);
subplot(fig_Re_Sh_sub(counter)); hold on;
set(gca,'ColorOrder', ColorSet,...
    'XScale','log','YScale','log',...
    'xlim',[1e2 1e6],'ylim',[0.01 100]);

sub_title = ['h_b = ' num2str(h_b(ii)) ' r_o, d_p = ' num2str(d_p(jj)*1000) ' mm'];
grid on;
box on;
title(sub_title,'FontSize',10);
xlabel('Re_i_s','FontSize',10);
ylabel('Sh_i_s / Sh_c_r_i_t','FontSize',10);



%%

% Size for symbols representing particle diameters    
% markersizelist = [10,35,80]; 
% ,'SizeData', markersizelist(jj)

% Loop all fluids
% plot_h = zeros(length(PV)+3,1);
% legendentries = cell(length(PV)+3,1);
% xx=1; xx=2; xx=3; xx=4; xx=5; xx=6; xx=7; xx=8; xx=9; xx=10;
set(gca,'ColorOrderIndex',1);
for xx=1:length(PV)
    % Color index
    COI = get(gca,'ColorOrderIndex');
    % Loop d_h and plot pointwise in order to express range with endpoints
    % and a line in between
    % yy=1; yy=2; yy=3; yy=4; yy=5; yy=6;
    for yy=1:max(size(d_h))
        % Pipe
%         plot_h(xx) = scatter([d_h(yy,1) d_h(yy,1)],Re_PL(yy,1:2,xx),'o','filled');
%         set(gca,'ColorOrderIndex',COI);
%         plot([d_h(yy,1) d_h(yy,1)],Re_PL(yy,1:2,xx));
%         set(gca,'ColorOrderIndex',COI);
        % Annulus        
        scatter([Re_PL(yy,3,xx) Re_PL(yy,6,xx)],[Sh_PL(yy,3,xx) Sh_PL(yy,6,xx)],'o','filled','SizeData', 10);
        set(gca,'ColorOrderIndex',COI); 
%         scatter(xdata,[Re_pipe_PL(yy,3,xx) Re_pipe_PL(yy,6,xx)],'wo','MarkerFaceColor', 'w'/2);
%         set(gca,'ColorOrderIndex',COI);
        plot([Re_PL(yy,3,xx) Re_PL(yy,6,xx)],[Sh_PL(yy,3,xx) Sh_PL(yy,6,xx)]);
        if yy==max(size(d_h))
            
        else           
            set(gca,'ColorOrderIndex',COI);
        end
    end
    % Critical Re
%     plot(d_h(:,2),ones(length(d_h),1)*Re_cr_PL(xx,1),'-');
%     set(gca,'ColorOrderIndex',COI);
%     plot(d_h(:,2),ones(length(d_h),1)*Re_cr_PL(xx,2),'--');
%     set(gca,'ColorOrderIndex',COI);
    % Legend (Fluid colors)
%     plot_h(xx) = scatter(100*[d_h(yy,1) d_h(yy,1)],Re_PL(yy,1:2,xx),'s','filled');
%     legendentries{xx} = ['Fluid ' num2str(xx), ' (YP = ', num2str(YP(xx),1), ', PV = ', num2str(PV(xx)),')'];
end
% Legend
% plot_h(xx+1) = scatter(10,10,'ko','filled'); % Pipe
% legendentries(xx+1)={'Pipe'};
% plot_h(xx+2) = scatter(10,10,'ko'); % Annulus
% legendentries(xx+2)={'Annulus'};
% plot_h(xx+3) = plot(10,10,'k-'); % Re_cr_PL = f(n_P)
% legendentries(xx+3)={'Re_c_r = f(n_P)'};
% plot_h(xx+4) = plot(10,10,'k--'); % Re_cr_PL = f(n_A)
% legendentries(xx+4)={'Re_c_r = f(n_A)'};
% legend(plot_h,legendentries,'Location','northwest');
% Format
% set(gca,...
%     'Xdir', 'reverse',...
%     'xlim',[0 1],...
%     'ylim',[1e2 2e6]);




