function H0(x)
    1.0
end
function H1(x)
    2.0*x
end
function H2(x)
    -2.0 + 4.0*x^2
end
function H3(x)
    -12.0*x + 8.0*x^3
end
function H4(x)
    12.0 - 48.0*x^2 + 16.0*x^4
end
function H5(x)
    120.0*x - 160.0*x^3 + 32.0*x^5
end
function H6(x)
    -120.0 + 720.0*x^2 - 480.0*x^4 + 64.0*x^6
end
function H7(x)
    -1680.0*x + 3360.0*x^3 - 1344.0*x^5 + 128.0*x^7
end
function H8(x)
    1680.0 - 13440.0*x^2 + 13440.0*x^4 - 3584.0*x^6 + 256.0*x^8
end
function H9(x)
    30240.0*x - 80640.0*x^3 + 48384.0*x^5 - 9216.0*x^7 + 512.0*x^9
end
function H10(x)
    -30240.0 + 302400.0*x^2 - 403200.0*x^4 + 161280.0*x^6 - 23040.0*x^8 + 1024.0*x^10
end

function hermite(x, i)
    if i==0
        return H0(x)
    end
    if i==1
        return H1(x)
    end
    if i==2
        return H2(x)
    end
    if i==3
        return H3(x)
    end
    if i==4
        return H4(x)
    end
    if i==5
        return H5(x)
    end
    if i==6
        return H6(x)
    end
    if i==7
        return H7(x)
    end
    if i==8
        return H8(x)
    end
    if i==9
        return H9(x)
    end
    if i==10
        return H10(x)
    end
end