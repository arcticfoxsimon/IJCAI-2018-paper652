include("graph.jl")
include("lapl.jl")
include("VRComp.jl")

function Appro_SM3(v, G, k) #v- the vertex chosed; G- graph; k- the number of edge added

    n = G.n;
    cand = zeros(n) #in general situation; candidate edges are given as a parameter.
    L = lapl(G)
    for i = 1 : n
        if L[v,i] == 0
            cand[i] = 1
        end
    end    # candidate edges are edges which are not existant in the original graph

    ans = zeros(Int64,k)
    A = adja(G)

#return dR = VRComp(v,A,cand,0.5)
Ldag = mppinv(L)
R0 = n*Ldag[v,v]
for i = 1:n
    R0 = R0 + Ldag[i,i]
end
#print("Original Resistant:")
#println(R0)

#@time begin

    for dep = 1:k
        max = 0.0
        u = 0
        dR = VRComp(v, A, cand, 0.5)
        for i = 1: n
            if cand[i] == 1
                if dR[i] > max
                    max = dR[i];
                    u = i
                end
            end
        end
        ans[dep] = u
        cand[u] = 0
        A[v,u] = 1
        A[u,v] = 1
        L[v,u] = -1
        L[u,v] = -1
        L[v,v] = L[v,v]+1
        L[u,u] = L[u,u]+1
    end

#end

    Ldag = mppinv(L)
    R0 = n*Ldag[v,v]
    for i = 1:n
        R0 = R0 + Ldag[i,i]
    end
    #print("Final Resistant:")
    #println(R0)

    return n/R0
end
