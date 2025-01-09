
CreateFigure('Flowcurves', 'Shear Rate [1/s]','Apparent viscosity [Pa.s]');
Legend = cell(2*length(fluidlist),1);

SR = logspace(-2,6);

for i = 1:length(fluidlist)

    % Plot Cross fit
    plot(SR,Cross.mu_inf{i}+(Cross.mu_0{i}-Cross.mu_inf{i})./(1+(Cross.lambda{i}.*SR).^Cross.n{i}),'--','Color',colorlist{i});

    % Plot Carreau fit
    plot(SR,Carreau.mu_inf{i}+(Carreau.mu_0{i}-Carreau.mu_inf{i}).*((1+(Carreau.lambda{i}.*SR).^2).^((Carreau.n{i}-1)./2)),'-','Color',colorlist{i});

end
