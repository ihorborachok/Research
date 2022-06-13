function BuildNorms2()
    global problem;

    problem.results.normPnts = BuildNormPnts();
    problem.results.computeNorm = @() ComputeNorm();
end

function normPnts = BuildNormPnts()
    global problem;

    s = BuildNormPntsS();

    normPnts.s = s;
    normPnts.x = BuildPntsX(s);
    normPnts.nu = BuildPntsNu(s);
    normPnts.t = BuildNormPntsT();
end

function [unrm, udnrm] = ComputeNorm()
    global problem;

    x = problem.results.normPnts.x;
    nu = problem.results.normPnts.nu;
    t = problem.results.normPnts.t;

    [app, appd] = AppSln2(x, t, nu);
    [ex, exd] = problem.example.exsln(x, t, nu);

    unrm = norm(app-ex);
    udnrm = norm(appd-exd);
end

function x = BuildPntsX(s)
    global problem;
    x = problem.example.gamma1(s); #todo: hardcoded gm1
endfunction

function nu = BuildPntsNu(s)
    global problem;
    nu = problem.example.nu1(s); #todo: hardcoded gm1
endfunction

function t = BuildNormPntsT()
    global problem;

    nt = problem.results.normT;
    ht = problem.example.maxT / nt;
    t = ht * (1 : nt); # starts from 1
endfunction

function s = BuildNormPntsS()
    global problem;

    n = problem.results.normN;
    h = pi / n;
    s = h * (0 : 2 * n - 1);
endfunction