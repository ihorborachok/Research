# Compute Approximate Solution(u) or\and Normal Derivative of Approximate Solution (ud)
# x - each column of x should be a vector (2D point); [[0; 1] [0; 2] [0; 3]].
# nu - normal derivative at x points, only used when computing normal derivative
# t - 1D vector
function [u, ud] = AppSln2(x, t, nu)
    global problem;
    
    problem.helper.log('App Sln computing started');

    pkg load symbolic # for Laguarre polynoms
    lgrPlns = ComputeLaguerrePlns(t); # matrix lgrPlns(t, p)

    computeU = @(up) problem.model.kappa * up * lgrPlns';

    if (isargout(1) && isargout(2))
        problem.helper.log('Computing Solution and Normal Derivative');
        [up, upd] = ComputeAppSlnMfs(x, nu); # matrix up(x, p)
        u = computeU(up);
        ud = computeU(upd);
    elseif (isargout(1))
        problem.helper.log('Computing Solution');
        up = ComputeAppSlnMfs(x);
        u = computeU(up);
    elseif (isargout(2))
        problem.helper.log('Computing Normal Derivative');
        [~, upd] = ComputeAppSlnMfs(x);
        ud = computeU(upd);
    end

    problem.helper.log('App Sln computed');
end

function lgrPlns = ComputeLaguerrePlns(t)
    global problem;

    for p = 0 : problem.model.Nt
        lgrPlns(:, p + 1) = laguerreL(p, problem.model.kappa * t);
    end
end

function [up, upd] = ComputeAppSlnMfs(x, nu)
    global problem;
    p = 0 : problem.model.Nt;
    y = problem.model.srcpnts;

    if (isargout(1) && isargout(2)) #compute only that part, we need
        [fs, fsd] = FS2(p, x, y, nu);
        up = ComputeAppSlnMfsFromFs(fs, x);
        upd = ComputeAppSlnMfsFromFs(fsd, x);
    elseif isargout(1)
        fs = FS2(p, x, y);
        up = ComputeAppSlnMfsFromFs(fs, x);
    elseif isargout(2)
        [~, fsd] = FS2(p, x, y, nu);
        upd = ComputeAppSlnMfsFromFs(fsd, x);
    end
end

function up = ComputeAppSlnMfsFromFs(fs, x)
    global problem;
    
    for p = 0 : problem.model.Nt;
        temp = zeros(size(x)(2), 1); # number of columns in x
        for m = 0 : p
            temp = temp + fs(:, :, p - m + 1) * problem.temp.mfs.alpha(:, m + 1);
        end
        up(:, p + 1) = temp;
    end
end