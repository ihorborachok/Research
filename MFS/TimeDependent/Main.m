function Main
    NtTstLst = 10*[1:2]; NmfsTstLst = 2.^ [2:3]; # lists of test points
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
            problem.plotting.plotAll();

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