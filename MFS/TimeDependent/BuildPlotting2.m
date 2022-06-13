# run after build problem
function BuildPlotting2()
    global problem;
    
    problem.plotting.plotPnts = BuildPlotPnts(); 
    problem.plotting.configure = configurePlotting();
    problem.plotting.plotDomain = @() plotDomain();
    problem.plotting.plotSolutions = @() plotSolutions();

    problem.plotting.plotAll = @() plotAll();
end

function plotPnts = BuildPlotPnts()
    global problem;

    s = BuildPlotPntsS();
    plotPnts.s = s;
    plotPnts.x = BuildPntsX(s);
    plotPnts.nu = BuildPntsNu(s);
    plotPnts.t = BuildPlotPntsT();
end

function configure = configurePlotting()
    configure.fontsize = 18;
    configure.linewidth = 3;
end

function plotAll()
    plotDomain();
    plotSolutions();
end

function plotDomain()
    global problem;

    if ~problem.results.plotDomain
        return;
    end

    configure = problem.plotting.configure;

    fig = figure(1);
    
    n = 40; h = pi / n; k = 0 : 2 * n; sk = k .* h;
    xk1 = problem.example.gamma1(sk); xk2 = problem.example.gamma2(sk);

    plot(xk1(1, :), xk1(2, :), 'displayname', '\Gamma_1', 'linewidth', configure.linewidth, ...
        xk2(1, :), xk2(2, :), 'displayname', '\Gamma_2', 'linewidth', configure.linewidth);
    
    leg = legend;   
    set(get(fig, "currentaxes"), "fontsize", configure.fontsize);
    set(leg, "fontsize", configure.fontsize);
    axis('equal');

    fileName = strcat(problem.results.folder, '\', 'domain');
    saveToFile(fig, fileName);
    
endfunction

# plot approximate \ exact solutions
function plotSolutions()
    global problem;

    s = problem.plotting.plotPnts.s;
    x = problem.plotting.plotPnts.x;
    nu = problem.plotting.plotPnts.nu;
    t = problem.plotting.plotPnts.t;

    if problem.results.plotSln && problem.results.plotSlnNd

        [app, appd] = AppSln2(x, t, nu);
        [ex, exd] = problem.example.exsln(x, t, nu);
    
    elseif problem.results.plotSln
    
        [app, ~] = AppSln2(x, t);
        [ex, ~] = problem.example.exsln(x, t);
    
    elseif problem.results.plotSlnNd
    
        [~, appd] = AppSln2(x, t, nu);
        [~, exd] = problem.example.exsln(x, t, nu);
    
    end

    [s, t] = meshgrid(s, t);

    if problem.results.plotSln
        plotAs3D(app', s, t, figure(2), BuildPlot3DObjConfig(true, false));
        plotAs3D(ex', s, t, figure(3), BuildPlot3DObjConfig(true, true));
    end
    if problem.results.plotSlnNd
        plotAs3D(app', s, t, figure(4), BuildPlot3DObjConfig(false, false));
        plotAs3D(ex', s, t, figure(5), BuildPlot3DObjConfig(false, true));
    end
end

function config = BuildPlot3DObjConfig(isSln, isEx)
    global problem;

    if isSln
        config.zLabel = 'u';
    else
        config.zLabel = '\partial u \ \partial\nu';
    end

    fName = [problem.results.folder '\u'];
    if ~isSln
        fName = [fName 'd'];
    end
    fName = [fName, '_', num2str(problem.model.Nt), '_' num2str(problem.model.Nmfs), '_'];

    if isEx
        fName = [fName, 'ex'];
    else
        fName = [fName, 'app'];
    end

    config.fName = fName;
end

# plot u(x(s), t) function as 3D object
function plotAs3D(u, s, t, fig, objConfig)
    global problem;
    configure = problem.plotting.configure;

    surf(s, t, u, 'FaceColor','interp','EdgeColor','black','FaceLighting','none');
    xlabel('s'); ylabel('t'); zlabel(objConfig.zLabel);
    colormap gray;
    set(get(fig, "currentaxes"), "fontsize", configure.fontsize);

    saveToFile(fig, objConfig.fName);
end

function saveToFile(fig, fname)
    global problem
    
    set(fig, 'PaperOrientation', 'landscape');

    exts = problem.results.plotExtensions;
    for extIdx = 1:length(exts)
        
        ext = char(exts(extIdx));
        ffname = strcat(fname, strcat('.', ext));
        
        if ext ~= 'fig'
            print(fig, ffname, strcat('-d', ext));
        else
            savefig(fig, ffname);
        end

    end
endfunction

function s = BuildPlotPntsS()
    global problem;

    n = problem.results.plotN;
    h = pi / n;
    s = h * (0 : 2 * n - 1);
endfunction

function x = BuildPntsX(s)
    global problem;
    x = problem.example.gamma1(s); #todo: hardcoded gm1
endfunction

function nu = BuildPntsNu(s)
    global problem;
    nu = problem.example.nu1(s); #todo: hardcoded gm1
endfunction

function t = BuildPlotPntsT()
    global problem;

    nt = problem.results.plotT;
    ht = problem.example.maxT / nt;
    t = ht * (1 : nt); # starts from 1
endfunction