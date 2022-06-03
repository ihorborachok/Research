# run after build problem
function BuildPlotting2()
    global problem;

    problem.plotting.configure = configurePlotting();
    problem.plotting.plotDomain = @() plotDomain(problem.plotting.configure);

endfunction

function configure = configurePlotting()
    configure.fontsize = 18;
    configure.linewidth = 3;
    configure.saveDomainToFile = true;
endfunction

function plotDomain(configure)
    global problem;

    fig = figure(1);
    
    n = 40; h = pi / n; k = 0 : 2 * n; sk = k .* h;
    xk1 = problem.example.gamma1(sk); xk2 = problem.example.gamma2(sk);

    plot(xk1(1, :), xk1(2, :), 'displayname', '\Gamma_1', 'linewidth', configure.linewidth, ...
        xk2(1, :), xk2(2, :), 'displayname', '\Gamma_2', 'linewidth', configure.linewidth);
    
    leg = legend;   
    set(get(fig, "currentaxes"), "fontsize", configure.fontsize);
    set(leg, "fontsize", configure.fontsize);
    axis('equal');

    if problem.results.plotDomain
        fileName = strcat(problem.results.folder, '\', 'domain');
        saveToFile(fig, fileName);
    endif

endfunction

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