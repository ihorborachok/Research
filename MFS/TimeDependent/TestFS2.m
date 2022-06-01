function TestFS2()
    global gamma;
    global kappa;
    global beta;
    global N;

    kappa = 1; N = 50;
    beta = ones(N + 1, 1) * kappa;
    gamma = sqrt(beta(1));

    x = [[0; 1] [0; 2] [0; 3]]; y = [[1; 2] [3; 4]]; p = 1:10;
    #x = [zeros(1, 50); 1:50]; y = [1:20; 2:21]; p = 1:50;

    FS2(p, x, y)
endfunction