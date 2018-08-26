include("lapl.jl")
include("graph.jl")
include("er.jl")
using Laplacians
function VRComp(v, A, cand, eps) #v- vertex chosed; G- graph; cand- candidate edges; eps- epsilon
    n = G.n
    M = round(Int64,(log2(2*n)/eps)/eps);

    Z = zeros(n,M);
    s = [1, -1]
    for i = 1 : n
        for j = 1 : M
            Z[i,j] = rand(s);
        end
    end

    A = sparse(A)
    f = approxCholLap(A,tol=1e-5);

    Y = zeros(n,M);
    for i = 1: M
        Y[:,i] .= f(Z[:,i]);
    end

    ev = zeros(n,1);
    x = zeros(n,1);
    ev[v] = 1;
    x[:,1] = f(ev[:,1]);

    dRv = zeros(n,1);
    r = ApproxiER(A, v, cand);
    for i = 1:n
        if cand[i] > 0
            dRv[i] = n*(x[v]-x[i])*(x[v]-x[i]);
            t = 0;
            for j = 1:M
                t = t + (Y[v,j]-Y[i,j])*(Y[v,j]-Y[i,j]);
            end
            dRv[i] = (dRv[i] + t / M) / (1 + r[i]);
        end
    end

    return dRv
end
