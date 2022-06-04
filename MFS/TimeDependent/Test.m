# Temp file with some test algorithm parts
function Test()
    #TestFS2();
    #TestAppSln2();
    #TestPlot();
    #TestBatchPlot();
    TestBatchNorm();
endfunction

function TestPlot()
    config.Nt = 10;
    config.Nmfs = 8;
    BuldProblem2(config);

    global problem;
    MFS2();
    problem.plotting.plotAll();
end

function TestBatchPlot()
    for Nt = 10:10:20
        for Nmfs = 2.^ [3:4]
            clear -global problem

            config.Nt = Nt;
            config.Nmfs = Nmfs;
            BuldProblem2(config);

            global problem;
            problem.helper.log(['Nt = ', num2str(Nt), ',     Nmfs = ', num2str(Nmfs)]);
            MFS2();
            problem.plotting.plotAll();
        end
    end
end

function TestBatchNorm()
    
    NtLst = 10*[5:7]; NmfsLst = 2.^ [7:9];
    err = [NmfsLst]; errd = [NmfsLst];

    for Nt = NtLst
        errR = []; errdR = [];
        for Nmfs = NmfsLst
            clear -global problem

            config.Nt = Nt;
            config.Nmfs = Nmfs;
            BuldProblem2(config);

            global problem;
            problem.helper.log(['Nt = ', num2str(Nt), ',     Nmfs = ', num2str(Nmfs)]);
            MFS2();
            
            s = problem.results.plotPnts.s;
            x = problem.results.plotPnts.x;
            nu = problem.results.plotPnts.nu;
            t = problem.results.plotPnts.t;

            [app, appd] = AppSln2(x, t, nu);
            [ex, exd] = problem.example.exsln(x, t, nu);

            errdPrnt = norm(appd-exd)
            
            errR = [errR norm(app-ex)];
            errdR = [errdR norm(appd-exd)];
        end
        err = [err; errR]; errd = [errd; errdR];
    end
    frstCol = [0 NtLst]';
    err = [frstCol err];
    errd = [frstCol errd];
    err
    errd
end

function TestAppSln2()
    config.Nt = 10;
    config.Nmfs = 8;
    BuldProblem2(config);

    global problem;
    #x = [0; 1]; t = 1:5;
    x = [[0.7; 0.6] [0.8; 0.7]]; t = 1:3; nu = x;
    MFS2();
    [app, appd] = AppSln2(x, t, nu);
    [ex, exd] = problem.example.exsln(x, t, nu);
    err = norm(app-ex)
    errd = norm(appd-exd)
endfunction

function TestFS2()
    config.Nt = 10;
    config.Nmfs = 8;
    BuldProblem2(config);

    x = [[0; 1] [0; 2] [0; 3]]; y = [[1; 2] [3; 4]]; p = 0:1; nu = [[0; 1] [0; 2] [0; 3]];
    #x = [zeros(1, 50); 1:50]; y = [1:20; 2:21]; p = 1:50;

    [fs, fsd] = FS2(p, x, y, nu);
    #fs
    fsd
endfunction