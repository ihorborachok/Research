function BuldProblem2(config)
    global problem;

    problem.model = BuildModel2(config);
    problem.example = BuldExample2();
    problem.helper = BuildHelper();
    problem.results = BuildResults();

    BuildPlotting2();
    BuildNorms2();

endfunction

function example = BuldExample2()
    
    example.gamma1 = @(s) ComputeGamma1(s);
    example.gamma2 = @(s) ComputeGamma2(s);

    example.nu1 = @(s) ComputeNu1(s);
    example.nu2 = @(s) ComputeNu2(s);
    
    example.PntFs = [0; 4];
    example.RghtSd = ComputeRightSide(example.PntFs);
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

function nu = ComputeNu1(s)
    gmd = [-0.6*sin(s); 0.5*cos(s)];
    nu = ComputeNu(gmd);
end

function nu = ComputeNu2(s)
    gmd = [-sin(s); -0.5*sin(2.*s)+cos(s)];
    nu = ComputeNu(gmd);
end

function nu = ComputeNu(gmd)
    nrm = sqrt(gmd(1, :).^2 + gmd(2, :).^2);
    nu = [-gmd(2, :) ; gmd(1, :)] ./ nrm;
end

function f = ComputeRightSide(pntFs)
    global main;

    if main.type.DD
        f = ComputeExampleRightSide_FromFsNarrowing(pntFs);
    elseif main.type.C
        f = [ComputeExampleRightSide_FromFsNarrowing(pntFs); ComputeExampleRightSide_FromFsDNarrowing(pntFs)];
    end
end

function f = ComputeExampleRightSide_FromFsNarrowing(pntFs) # right side of the problem: f1, f2
    global problem;

    p = 0 : problem.model.Nt;
    x = problem.model.collpnts; # for MFS system x = coll points
    y = pntFs; # y = pntFs; u_ex(x) = f(x) = FS(x, pntFs), x \in \gamma_{1, 2}

    nrvfs = FS2(p, x, y);

    mult = 1 / (2 * pi) * 100;
    f(:, :) = nrvfs(:, 1, :) * mult;
end

function f = ComputeExampleRightSide_FromFsDNarrowing(pntFs)  # right side of the problem: f2, g2
    global problem;

    p = 0 : problem.model.Nt;
    x = problem.model.collpnts;
    nu = problem.model.collnu;
    y = pntFs;

    [~, nrvfs] = FS2(p, x, y, nu);

    mult = 1 / (2 * pi) * 100;
    f(:, :) = nrvfs(:, 1, :) * mult;
end

function model = BuildModel2(config)

    model.kappa = 1;
    model.Nt = config.Nt;
    model.Nmfs = config.Nmfs;
    model.beta = model.kappa * ones(model.Nt + 1, 1);
    model.gamma = sqrt(model.beta(1));
    model.srcpnts = ComputeSourcePoints(model.Nmfs);
    [model.collpnts, model.collnu] = ComputeCollocationPoints(model.Nmfs);

endfunction

function helper = BuildHelper()
    helper.log = @(msg) Log(msg, true); # Log enabled
endfunction

# y - 2xNmfs matrix, each column is 2D source point
function y = ComputeSourcePoints(Nmfs)

    step = 2 * pi / Nmfs;
    sGm1 = (1 : 2 : (Nmfs - 1)) * step;
    sGm2 = (2 : 2 : Nmfs) * step;

    y = [0.3 * ComputeGamma1(sGm1) 3 * ComputeGamma2(sGm2)];


    #step = 2 * pi / Nmfs;
    #sGm11 = (1 : 4 : Nmfs) * step;
    #sGm12 = (2 : 4 : Nmfs) * step;
    #sGm21 = (3 : 4 : Nmfs) * step;
    #sGm22 = (4 : 4 : Nmfs) * step;

    #y = [0.4 * ComputeGamma1(sGm11) 3 * ComputeGamma2(sGm21) 4 * ComputeGamma2(sGm12) 2 * ComputeGamma2(sGm22)];

endfunction

# x - 2 x Nmfs matrix, each column is 2D collocation point
function [x, nu] = ComputeCollocationPoints(Nmfs)
    global main;
    
    if main.type.DD # Dirichlet problem; Nmfs / 2 points are on gm1, and Nmfs / 2 points - on gm2

        idx = 1 : (Nmfs / 2);
        s = idx * 4 * pi / (Nmfs + 1);

        x = [ComputeGamma1(s) ComputeGamma2(s)];
        nu = [];

    elseif main.type.C # Cauchy problem; Nmfs / 2 points are on gm2

        idx = 1 : (Nmfs / 2);
        s = idx * 4 * pi / (Nmfs + 1);

        x = ComputeGamma2(s);
        nu = ComputeNu2(s);

    end

end

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

    results.plotT = 10;
    results.plotN = 10;
    
    results.normT = 20;
    results.normN = 32;
end