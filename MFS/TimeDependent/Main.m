function Main
    ConfigureMain();
    global main;

    NtTstLst = [20]; NmfsTstLst = [32]; # lists of test parameters
    err = [NmfsTstLst]; errd = [NmfsTstLst];

    for Nt = NtTstLst
        errR = []; errdR = [];
        for Nmfs = NmfsTstLst
            clear -global problem

            config.Nt = Nt;
            config.Nmfs = Nmfs;
            BuldProblem2(config);

            global problem;
            problem.helper.log(['Nt = ', num2str(Nt), ', Nmfs = ', num2str(Nmfs)]);
            
            MFS2();
            [unrm, udnrm] = problem.results.computeNorm();

            if main.plot
                problem.plotting.plotAll();
            end

            problem.helper.log(['e_u = ', num2str(unrm), ', e_ud = ', num2str(udnrm)]);

            errR = [errR unrm]; errdR = [errdR udnrm];
        end
        err = [err; errR]; errd = [errd; errdR];
    end
    frstCol = [0 NtTstLst]';
    err = [frstCol err]; errd = [frstCol errd];
    err
    errd
end

function ConfigureMain
    global main

    # Type of the problem
    main.type.DD = false; # Dirichlet cond on gm1, Dirichlet cond on gm2
    main.type.C = true; # Cauchy problem on gm2

    main.plot = false;
end