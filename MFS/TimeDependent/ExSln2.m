function u = ExSln2(x, t)
    global problem;
    u = [];
    for xEl = x
        u = [u; arrayfun(@(t) 1 / (4 * pi * t) * exp(-norm(xEl - problem.example.PntFs)^2 / (4 * t)) * 100, t)];
    end
end