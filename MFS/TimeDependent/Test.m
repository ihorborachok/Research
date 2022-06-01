# Temp file with some test algorithm parts
function Test()
    BuldProblem2();
    global problem;

    MFS2()
endfunction

function TestFS2()
    BuldProblem2();

    x = [[0; 1] [0; 2] [0; 3]]; y = [[1; 2] [3; 4]]; p = 1:10;
    #x = [zeros(1, 50); 1:50]; y = [1:20; 2:21]; p = 1:50;

    FS2(p, x, y)
endfunction