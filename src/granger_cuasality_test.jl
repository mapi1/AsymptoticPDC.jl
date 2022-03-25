function granger_causality_test(u, A, pf; α = 0.01, verbose::Bool = true)
    nChannels, _, order = size(A)
    Z = get_Z(u, order)
    gamma = Z * Z'
    b = reshape(A, nChannels * nChannels * order)

    n, m = size(pf)
    Va = zeros(n,m)
    Tr = zeros(n,m)
    CO = zeros(n,m)
    pValue = zeros(n,m)
    for i in 1:n
        for j in 1:n
           if i != j
              CO[i,j]=1
              Tr[i,j],Va[i,j],v,th,pValue[i,j]=grangt(CO,b,gamma,pf, 1 - α);
              CO[i,j]=0;
           end
        end
     end
     if verbose
        println("                  GRANGER CAUSALITY TEST")
        println("======================================================================")
        println("Connectivity matrix:")
        Tr[I(nChannels)] .= NaN
        display(Tr)
        println("Granger causality test p-values:")
        pValue[I(nChannels)] .= NaN
        display(pValue)
        println("======================================================================")
     end
     return Tr, pValue
end

function grangt(CO,b,G,SU,significance)
    n,m = size(CO)
    lb = length(b)
    p0 = round(Int, lb / (m*n))
    Ct = repeat(reshape(CO,m*n), p0)

    v = sum(Int, Ct);
    C = zeros(v,m*n*p0);
    for i = 1:v
        p ,q = findmax(Ct)
        C[i,q] = 1;
        Ct[q] = 0;
    end

    value = (C*b)' * inv(C*kron(inv(G),SU)*C')*C*b

    test_dist = Chisq(v)
    th = invlogcdf(test_dist, log(significance))
    y = value >= th;

    pValue = 1 - cdf(test_dist, value)
    return y, value, v, th, pValue
end