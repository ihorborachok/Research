function BuldProblem2()
    global problem;

    problem.model = BuildModel2();
    problem.example = BuldExample2();
    problem.helper = BuildHelper();

endfunction

function example = BuldExample2()
    
    example.gamma1 = @(s) gamma1(s);
    example.gamma2 = @(s) gamma2(s);
    example.PntFs = [0; 4];
    example.RghtSd = ComputeExampleRightSide_FromFsNarrowing(example.PntFs); # right side of the problem: f1, f2
    example.exsln = @(x, t) ComputeExactSolution_FsNarrowing(x, t, example.PntFs); # exact solution
end

# extarnal boundary
function gm = ComputeGamma1(s)
    gm = [0.6*cos(s); 0.5*sin(s)];
endfunction

# internal boundary
function gm = ComputeGamma2(s)
    gm = [cos(s); sin(s)-0.5*sin(s).^2+0.5];
endfunction

function u = ComputeExactSolution_FsNarrowing(x, t, pntFs)
    u = 1 / (4 * pi * t) * exp(-norm(x - pntFs)^2 / (4 * t)) * 100;
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

function model = BuildModel2()

    model.kappa = 1;
    model.Nt = 20;
    model.Nmfs = 32;
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