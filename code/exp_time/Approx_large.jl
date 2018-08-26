include("graph2.jl")
include("lapl.jl")
include("VRComp2.jl")
include("er3.jl")

function Appro_large(v, G, A, k) #v- the vertex chosed; G- graph; k- the number of edge added

    n = G.n;
    cand = ones(n) #in general situation; candidate edges are given as a parameter.
#=
    L = lapl(G)
    for i = 1 : n
        if L[v,i] == 0
            cand[i] = 1
        end
    end    # candidate edges are edges which are not existant in the original graph
=#
    p = size(G.nbr[v],1)
    for i = 1:p
        c = G.nbr[v][i]
        cand[c] = 0
    end

    ans = zeros(Int64,k)
#    A = sparseAdja2(G)

#return dR = VRComp(v,A,cand,0.5)
#=
Ldag = mppinv(L)
R0 = n*Ldag[v,v]
for i = 1:n
    R0 = R0 + Ldag[i,i]
end
print("Original Resistant:")
println(R0)
=#

start_time = time()

    for dep = 1:k
        println("edges:", dep)
        max = 0.0
        u = 0
        dR = VRComp2(v, A, cand, 0.5)
        for i = 1: n
            if cand[i] == 1
                if dR[i] > max
                    max = dR[i];
                    u = i
                end
            end
        end
#        ans[dep] = u
        cand[u] = 0
        A[v,u] = 1
        A[u,v] = 1
#        L[v,u] = -1
#        L[u,v] = -1
#        L[v,v] = L[v,v]+1
#        L[u,u] = L[u,u]+1
    end

time_elapse = time() - start_time
#println("time: ",time_elapse)

#=
    Ldag = mppinv(L)
    R0 = n*Ldag[v,v]
    for i = 1:n
        R0 = R0 + Ldag[i,i]
    end
=#
    tmp = ones(n)
    r = ApproxiER3(A, v,tmp)
    #return r

    R0 = 0
    for i = 1:n
        if i != v
            R0 = R0 + r[i]
        end
    end

    #print("Final Resistant:")
#    println("Final centrality:", n/R0)

    res = zeros(2)
    res[1] = time_elapse
    res[2] = n/R0

#println("time: ",time()-start_time)

    return res
end
