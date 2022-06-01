function BuldExample2()
    global example;
    
    example.gamma1 = @(s) gamma1(s);
    example.gamma2 = @(s) gamma2(s);
end

# extarnal boundary
function gm = ComputeGamma1(s)
    gm = [0.6*cos(s); 0.5*sin(s)];
endfunction

# internal boundary
function gm = ComputeGamma2(s)
    gm = [cos(s); sin(s)-0.5*sin(s).^2+0.5];
endfunction