function BuldProblem2(config)
    global problem;

    problem.model = BuildModel2(config);
    problem.example = BuldExample2();
    problem.helper = BuildHelper();
    problem.results = BuildResults();

    BuildPlotting2();

endfunction

function example = BuldExample2()
    
    example.gamma1 = @(s) ComputeGamma1(s);
    example.gamma2 = @(s) ComputeGamma2(s);

    example.gamma1d = @(s) ComputeGamma1Derivative(s);
    example.gamma2d = @(s) ComputeGamma2Derivative(s);
    
    example.PntFs = [0; 4];
    example.RghtSd = ComputeExampleRightSide_FromFsNarrowing(example.PntFs); # right side of the problem: f1, f2
    example.exsln = @(x, t, nu) ComputeExSln(x, t, nu);
    example.maxT = 5;
end

function [u, ud] = ComputeExSln(x, t, nu)
    global problem;

    if isargout(1)
        u = [];
        for xEl = x
            u = [u; arrayfun(@(t) 1 / (4 * pi * t) * exp(-norm(xEl - problem.example.PntFs)^2 / (4 * t)) * 100, t)];
        end
    end

    if isargout(2) # compute exact notrmal derivative
        ud = [];
        for xPointIdx = 1 : columns(x)
            xEl = x(:, xPointIdx); nuEl = nu(:, xPointIdx);
            ud = [ud; arrayfun(@(t) 1 / (8 * pi * t^2) * exp(-norm(xEl - problem.example.PntFs)^2 / (4 * t)) * dot((xEl - problem.example.PntFs), nuEl) * 100, t)];
        end
    end
end

# extarnal boundary
function gm = ComputeGamma1(s)
    gm = [0.6*cos(s); 0.5*sin(s)];
endfunction

# internal boundary
function gm = ComputeGamma2(s)
    gm = [cos(s); sin(s)-0.5*sin(s).^2+0.5];
endfunction

function gmd = ComputeGamma1Derivative(s)
    gmd = [-0.6*sin(s); 0.5*cos(s)];
endfunction

# internal boundary
function gmd = ComputeGamma2Derivative(s)
    gmd = [-sin(s); -0.5*sin(2.*s)+cos(s)];
endfunction

function f = ComputeExampleRightSide_FromFsNarrowing(pntFs)
    global problem;

    p = 0 : problem.model.Nt;
    x = problem.model.collpnts; # for MFS system x = coll points
    y = pntFs; # y = pntFs; u_ex(x) = f(x) = FS(x, pntFs), x \in \gamma_{1, 2}

    nrvfs = FS2(p, x, y);

    mult = 1 / (2 * pi) * 100;
    f(:, :) = nrvfs(:, 1, :) * mult;
endfunction

function model = BuildModel2(config)

    model.kappa = 1;
    model.Nt = config.Nt;
    model.Nmfs = config.Nmfs;
    model.beta = model.kappa * ones(model.Nt + 1, 1);
    model.gamma = sqrt(model.beta(1));
    model.collpnts = ComputeCollocationPoints(model.Nmfs);
    model.srcpnts = ComputeSourcePoints(model.Nmfs);

endfunction

function helper = BuildHelper()
    helper.log = @(msg) Log(msg, true); # Log enabled
endfunction

# x - 2x(Nmfs/2) matrix, each column is 2D collocation point
function x = ComputeCollocationPoints(Nmfs)
    idx = 1 : (Nmfs / 2);
    s = idx * 4 * pi / (Nmfs + 1);

    x = [ComputeGamma1(s) ComputeGamma2(s)];
endfunction

# y - 2xNmfs matrix, each column is 2D source point
function y = ComputeSourcePoints(Nmfs)
    step = 2 * pi / Nmfs;
    sGm1 = (1 : 2 : (Nmfs - 1)) * step;
    sGm2 = (2 : 2 : Nmfs) * step;

    y = [0.5 * ComputeGamma1(sGm1) 2 * ComputeGamma2(sGm2)];
endfunction

function Log(msg, enabled)
    if (enabled)
        disp(msg);
    end
end

function results = BuildResults()
    results.folder = 'D:\Private\LNU\Research\MFS\TimeDependent\results';

    results.plotExtensions = {'pdf', 'fig'};

    results.plotDomain = false;
    
    results.plotSln = true; # plot solution
    results.plotSlnNd = true; # plot normal derivate

    results.plotPnts.s = BuildPlotPntsS();
    results.plotPnts.x = BuildPlotPntsX();
    results.plotPnts.nu = BuildPlotPntsNu();
    results.plotPnts.t = BuildPlotPntsT();

end

function s = BuildPlotPntsS()
    global problem;

    n = 10;
    h = pi / n;
    s = h * (0 : 2 * n - 1);
endfunction

function x = BuildPlotPntsX()
    global problem;
    x = problem.example.gamma1(BuildPlotPntsS()); #todo: hardcoded gm1
endfunction

function nu = BuildPlotPntsNu()
    global problem;
    nu = problem.example.gamma1d(BuildPlotPntsS()); #todo: hardcoded gm1
endfunction

function t = BuildPlotPntsT()
    global problem;

    nt = 10;
    ht = problem.example.maxT / nt;
    t = ht * (1 : nt); # starts from 0 or 1?
endfunction