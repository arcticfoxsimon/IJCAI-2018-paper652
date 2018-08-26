include("lapl.jl")
include("graph.jl")
include("er2.jl")
using Laplacians
function VRComp2(v, A, cand, eps) #v- vertex chosed; G- graph; cand- candidate edges; eps- epsilon
    n = G.n
    M = round(Int64,(log2(2*n)/eps)/eps);

#=
    Z = zeros(n,M);
    s = [1, -1]
    for i = 1 : n
        for j = 1 : M
            Z[i,j] = rand(s);
        end
    end
=#
    #A = sparse(A)
    f = approxCholLap(A,tol=1e-5);

    #Y = zeros(n,M);
    dis = zeros(n,1)
    for i = 1: M
        q = zeros(n,1);
        s = [1, -1]
        for i = 1 : n
            q[i] = rand(s);
        end
        z=zeros(n,1)
        z.=f(q[:])
        for j = 1:n
            dis[j] = dis[j] + (z[j]-z[v])^2
        end
#        Y[:,i] .= f(Z[:,i]);
    end

    ev = zeros(n,1);
    x = zeros(n,1);
    ev[v] = 1;
    x[:,1] = f(ev[:,1]);

    dRv = zeros(n,1);
    r = ApproxiER2(A, v, cand);

    for i = 1:n
        if cand[i] > 0
            dRv[i] = n*(x[v]-x[i])*(x[v]-x[i]);
            t = 0;
            for j = 1:M
                #t = t + (Y[v,j]-Y[i,j])*(Y[v,j]-Y[i,j]);
                t = t + dis[i]
            end
            dRv[i] = (dRv[i] + t / M) / (1 + r[i]);
        end
    end

    return dRv
end
