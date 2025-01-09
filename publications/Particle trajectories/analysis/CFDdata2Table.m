% Assemble pathes for current case 
path_CFD = [ resultspath '\' fluidlist{i} '-cross_U' num2str(U(k)) '_d_p' num2str(d_p(j)) ];
path_CFD_asis = [ path_CFD '_asis' ];
path_CFD_new = [ path_CFD '_new' ];

% Read asis CFD exports if existing and write time, x, y, rep vectors to table
if exist([path_CFD_asis '_x.txt'], 'file')
    data = importCFDdata( [path_CFD_asis '_x.txt'] );
    CFD_asis.t(row,1) = {data(:, 1)};
    CFD_asis.X(row,1) = {data(:, 2)};
else
    CFD_asis.t(row,1) = {0};
    CFD_asis.X(row,1) = {0};
end

if exist([path_CFD_asis '_y.txt'], 'file')
    data = importCFDdata( [path_CFD_asis '_y.txt'] );
    CFD_asis.Y(row,1) = {data(:, 2)};
else
    CFD_asis.Y(row,1) = {0};
end

if exist([path_CFD_asis '_rep.txt'], 'file')
    data = importCFDdata( [path_CFD_asis '_rep.txt'] );
    CFD_asis.rep(row,1) = {data(:, 2)};
else
    CFD_asis.rep(row,1) = {0};
end

if exist([path_CFD_asis '_vx_xy.txt'], 'file')
    data = importCFDdata( [path_CFD_asis '_vx_xy.txt'] );
    CFD_asis.u_xy(row,1) = {data(:, 1)};
    CFD_asis.y(row,1) = {data(:, 2)};
else
    CFD_asis.u_xy(row,1) = {0};
    CFD_asis.y(row,1) = {0};
end


if exist([path_CFD_asis '_vx_xy_262mm.txt'], 'file')
    data = importCFDdata( [path_CFD_asis '_vx_xy_262mm.txt'] );
    CFD_asis.u_xy_262mm(row,1) = {data(:, 1)};
    CFD_asis.y_262mm(row,1) = {data(:, 2)};
else
    CFD_asis.u_xy_262mm(row,1) = {0};
    CFD_asis.y_262mm(row,1) = {0};
end


if exist([path_CFD_asis '_vx_xy_285mm.txt'], 'file')
    data = importCFDdata( [path_CFD_asis '_vx_xy_285mm.txt'] );
    CFD_asis.u_xy_285mm(row,1) = {data(:, 1)};
    CFD_asis.y_285mm(row,1) = {data(:, 2)};
else
    CFD_asis.u_xy_285mm(row,1) = {0};
    CFD_asis.y_285mm(row,1) = {0};
end

if exist([path_CFD_asis '_vx_xz.txt'], 'file')
    data = importCFDdata( [path_CFD_asis '_vx_xz.txt'] );
    CFD_asis.u_xz(row,1) = {data(:, 1)};
    CFD_asis.z(row,1) = {data(:, 2)};
else
    CFD_asis.u_xz(row,1) = {0};
    CFD_asis.z(row,1) = {0};
end

% Read new CFD exports if existing and write time, x, y, rep vectors to table
if exist([path_CFD_new '_x.txt'], 'file')
    data = importCFDdata( [path_CFD_new '_x.txt'] );
    CFD_new.t(row,1) = {data(:, 1)};
    CFD_new.X(row,1) = {data(:, 2)};
else
    CFD_new.t(row,1) = {0};
    CFD_new.X(row,1) = {0};
end

if exist([path_CFD_new '_y.txt'], 'file')
    data = importCFDdata( [path_CFD_new '_y.txt'] );
    CFD_new.Y(row,1) = {data(:, 2)};
else
    CFD_new.Y(row,1) = {0};
end

if exist([path_CFD_new '_rep.txt'], 'file')
    data = importCFDdata( [path_CFD_new '_rep.txt'] );
    CFD_new.rep(row,1) = {data(:, 2)};
else
    CFD_new.rep(row,1) = {0};
    CFD_new.delete(row,1) = 1;    
end

if exist([path_CFD_new '_vx_xy.txt'], 'file')
    data = importCFDdata( [path_CFD_new '_vx_xy.txt'] );
    CFD_new.u_xy(row,1) = {data(:, 1)};
    CFD_new.y(row,1) = {data(:, 2)};
else
    CFD_new.u_xy(row,1) = {0};
    CFD_new.y(row,1) = {0};
end



if exist([path_CFD_new '_vx_xy_262mm.txt'], 'file')
    data = importCFDdata( [path_CFD_new '_vx_xy_262mm.txt'] );
    CFD_new.u_xy_262mm(row,1) = {data(:, 1)};
    CFD_new.y_262mm(row,1) = {data(:, 2)};
else
    CFD_new.u_xy_262mm(row,1) = {0};
    CFD_new.y_262mm(row,1) = {0};
end


if exist([path_CFD_new '_vx_xy_285mm.txt'], 'file')
    data = importCFDdata( [path_CFD_new '_vx_xy_285mm.txt'] );
    CFD_new.u_xy_285mm(row,1) = {data(:, 1)};
    CFD_new.y_285mm(row,1) = {data(:, 2)};
else
    CFD_new.u_xy_285mm(row,1) = {0};
    CFD_new.y_285mm(row,1) = {0};
end




if exist([path_CFD_new '_vx_xz.txt'], 'file')
    data = importCFDdata( [path_CFD_new '_vx_xz.txt'] );
    CFD_new.u_xz(row,1) = {data(:, 1)};
    CFD_new.z(row,1) = {data(:, 2)};
else
    CFD_new.u_xz(row,1) = {0};
    CFD_new.z(row,1) = {0};
end

