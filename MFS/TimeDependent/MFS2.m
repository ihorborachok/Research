function MFS2()
    global problem;

    problem.helper.log('MFS started');

    p = 0 : problem.model.Nt;
    x = problem.model.collpnts;
    y = problem.model.srcpnts;
    problem.temp.mfs.fs = FS2(p, x, y);

    A = problem.temp.mfs.fs(:, :, 1); # p starts from 0

    for p = 0 : problem.model.Nt;
        b = ComputeB(p);
        problem.temp.mfs.alpha(:, p + 1) = A \ b;
    end

    problem.helper.log('MFS finished');
end

function b = ComputeB(p)
    global problem;

    b = problem.example.RghtSd(:, p + 1); # p starts from 0

    if(p == 0)
        return;
    end

    for m = 0 : (p - 1)
        b = b - problem.temp.mfs.fs(:, :, p - m + 1) * problem.temp.mfs.alpha(:, m + 1);
    end
end