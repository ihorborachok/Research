# Temp file with some test algorithm parts
function Test()
    BuldProblem2();
    TestAppSln2();
endfunction

function TestAppSln2()
    #x = [0; 1]; t = 1:5;
    x = [[0.7; 0.6] [0.8; 0.7]]; t = 1:3;
    MFS2();
    app = AppSln2(x, t)
    ex = ExSln2(x, t)
    err = norm(app-ex)
endfunction

function TestFS2()
    BuldProblem2();

    x = [[0; 1] [0; 2] [0; 3]]; y = [[1; 2] [3; 4]]; p = 1:10;
    #x = [zeros(1, 50); 1:50]; y = [1:20; 2:21]; p = 1:50;

    FS2(p, x, y)
endfunction