# Compute Approximate Solution
# x - each column of x should be a vector (2D point); [[0; 1] [0; 2] [0; 3]].
# t - 1D vector
function u = AppSln2(x, t)
    global problem;
    
    problem.helper.log('App Sln computing started');

    pkg load symbolic # for Laguarre polynoms
    
    up = ComputeAppSlnMfs(x); # matrix up(x, p)
    lgrPlns = ComputeLaguerrePlns(t); # matrix lgrPlns(t, p)

    u = problem.model.kappa * up * lgrPlns';

    problem.helper.log('App Sln computed');
end

function lgrPlns = ComputeLaguerrePlns(t)
    global problem;

    for p = 0 : problem.model.Nt
        lgrPlns(:, p + 1) = laguerreL(p, problem.model.kappa * t);
    end
end

function up = ComputeAppSlnMfs(x)
    global problem;
    p = 0 : problem.model.Nt;
    y = problem.model.srcpnts;
    fs = FS2(p, x, y);
    
    for p = 0 : problem.model.Nt;
        temp = zeros(size(x)(2), 1); # number of columns in x
        for m = 0 : p
            temp = temp + fs(:, :, p - m + 1) * problem.temp.mfs.alpha(:, m + 1);
        end
        up(:, p + 1) = temp;
    end
endfunction