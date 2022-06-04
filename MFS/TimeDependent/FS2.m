## p - 1D vector of indexes, could be 0, in a result matrix p index will be p + 1; p = [0 1 2].
## x - each column of x should be a vector (2D point); [[0; 1] [0; 2] [0; 3]].
## y - each column of y should be a vector (2D point); [[1; 2] [3; 4]].
## fs - matrix; fs(2, 1, 3) - fundamental solution for 2-nd point in x-vector and 1-st point in y-vector and 3-rd index from p
## fsd - normal derivative of fs
function [fs, fsd] = FS2(p, x, y, nu)
    global problem;
    problem.temp.a = CalculateA();

    if isargout(2)
        [R, Rd] = CalculateR(x, y, nu); %vector-columns
    else
        R = CalculateR(x, y); %vector-column
    end

    besselX = problem.model.gamma * R;
    besselAlpha = [0 1]; % must be a vector-row
    besselRes = besselk(besselAlpha, besselX); %firs column: K0; second column: K1

    
    for pEl = p
        vRes = arrayfun(@(r) v(pEl, r), R); # v_p
        wRes = arrayfun(@(r) w(pEl, r), R); # w_p

        if isargout(1)
            fsVector = sum(besselRes .*  [vRes wRes], 2); %multiply by elements: [K0*v  K1*w] and sum by row => result saved in vector => transform to matrix
            fs(:, :, pEl + 1) = reshape(fsVector, columns(y), columns(x))'; % p starts from 0 => index: p + 1
        end

        if isargout(2)  # compute normal derivative
            vdRes = arrayfun(@(r) vd(pEl, r), R); # v'_p
            wdRes = arrayfun(@(r) wd(pEl, r), R); # w'_p
            
            fsdVector = Rd .* sum(besselRes .*  [vdRes - problem.model.gamma * wRes, -problem.model.gamma * vRes - wRes ./ R + wdRes], 2);
            fsd(:, :, pEl + 1) = reshape(fsdVector, columns(y), columns(x))';
        end
    end

endfunction

function [R, Rd] = CalculateR(x, y, nu)
    rCounter = 1;
    for xPointIdx = 1 : columns(x)
        for yPointIdx = 1 : columns(y)
            
            R(rCounter) = norm(x(:, xPointIdx) - y(:, yPointIdx));

            if isargout(2) # calculate ((x-y), nu) / |x-y|
                Rd(rCounter) = -dot(x(:, xPointIdx) - y(:, yPointIdx), nu(:, xPointIdx)) / R(rCounter);
            end
            
            rCounter++;
        endfor
    endfor

    R = R'; % must be a vector-column
    if isargout(2)
        Rd = Rd';
    end
endfunction

function resV = v(p, r)
  global problem;
  a = problem.temp.a;
  
  m = 0 : floor(p / 2);
  resV = sum(a(p + 1, 2 .* m + 1) .* (r .^ (2 * m)));
endfunction

function res = w(p, r)
    global problem;
    a = problem.temp.a;

    if(p == 0)
    res = 0;
    return;
    endif

    m = 0 : floor((p - 1) / 2);
    res = sum(a(p + 1, 2 .* m + 2) .* (r .^ (2 * m + 1)));
endfunction

function res = vd(p, r)
    global problem;
    a = problem.temp.a;

    m = 0 : floor(p / 2);
    res = sum(2 .* m .* a(p + 1, 2 .* m + 1) .* (r .^ (2 .* m - 1)));
endfunction

function res = wd(p, r)
    global problem;
    a = problem.temp.a;

    if(p == 0)
        res = 0;
        return;
    endif

    m = 0 : floor((p - 1) / 2);
    res = sum((2 * m + 1) .* a(p + 1, 2 .* m + 2) .* (r .^ (2 * m)));
endfunction

function a = CalculateA()
    global problem;
    Nt = problem.model.Nt;
    gamma = problem.model.gamma;
    beta = problem.model.beta;

    a = zeros(Nt + 1, Nt + 1);
    for p = 0 : Nt
    
        a(p + 1, 1) = 1;
        
        if p ~= 0
            a(p + 1, p + 1) = - 1 / (2 * gamma * p) * beta(2) * a(p, p);
        endif
        
        for k = p - 1 : -1 : 1
            m = k - 1 : 1 : p - 1;
            a(p + 1, k + 1) = 1 / (2 * gamma * k) * (4 * floor((k + 1) / 2) ^ 2 * a(p + 1, k + 2) - sum(beta(p - m + 1) .* a(m + 1, k)));
        endfor  
    
    endfor  
endfunction