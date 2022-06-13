function MFS2()
    global problem;

    problem.helper.log('MFS started');

    A = ComputeA();

    for p = 0 : problem.model.Nt;
        b = ComputeB(p);
        problem.temp.mfs.alpha(:, p + 1) = Solve(A, b);
    end

    problem.helper.log('MFS finished');
end

function b = ComputeB(p)
    global main;
    global problem;

    b = problem.example.RghtSd(:, p + 1); # p starts from 0

    if(p == 0)
        return;
    end

    for m = 0 : (p - 1)
        if main.type.DD

            b = b - problem.temp.mfs.fs(:, :, p - m + 1) * problem.temp.mfs.alpha(:, m + 1);

        elseif main.type.C

            b = b - [problem.temp.mfs.fs(:, :, p - m + 1); problem.temp.mfs.fsd(:, :, p - m + 1)] * problem.temp.mfs.alpha(:, m + 1);

        end
    end
end

function A = ComputeA()
    global main;
    global problem;

    p = 0 : problem.model.Nt;
    x = problem.model.collpnts;
    y = problem.model.srcpnts;

    if main.type.DD # Dirichlet-Dirichlet
        
        problem.temp.mfs.fs = FS2(p, x, y);
        A = problem.temp.mfs.fs(:, :, 1); # p starts from 0;

    elseif main.type.C # Cauchy

        nu = problem.model.collnu;
        [problem.temp.mfs.fs, problem.temp.mfs.fsd] = FS2(p, x, y, nu);
        A = [problem.temp.mfs.fs(:, :, 1); problem.temp.mfs.fsd(:, :, 1)];

    end

end

function x = Solve(A, b)
    global main;

    if main.type.C # ill posed
        x = SolveIllPosed(A, b);
    else
        x = SolveWellPosed(A, b);
    end
end

function x = SolveIllPosed(A, b)
    lambda = 1e-10;
    x = (A' * A + lambda * eye(size(A))) / A' * b;
end

function x = SolveWellPosed(A, b)
    x = A \ b;
end