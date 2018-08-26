include("graph.jl")
include("lapl.jl")
include("VRComp.jl")

function Appro_SM2(v, G, k) #v- the vertex chosed; G- graph; k- the number of edge added

    n = G.n;
    cand = zeros(n) #in general situation; candidate edges are given as a parameter.
    L = lapl(G)
    for i = 1 : n
        if L[v,i] == 0
            cand[i] = 1
        end
    end    # candidate edges are edges which are not existant in the original graph

    ans = zeros(k+1)
    A = adja(G)

#return dR = VRComp(v,A,cand,0.5)
#print("Original Resistant:")
#println(R0)

#@time begin
Ld = mppinv(L)
R0 = n*Ld[v,v]
for i = 1:n
    R0 = R0 + Ld[i,i]
end

#println(f1,n/R0)
ans[1] = n/R0

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
        #ans[dep] = u
        cand[u] = 0
        A[v,u] = 1
        A[u,v] = 1
        L[v,u] = -1
        L[u,v] = -1
        L[v,v] = L[v,v]+1
        L[u,u] = L[u,u]+1

        c = zeros(n,1)
        for i = 1:n
            c[i] = Ld[i,v] - Ld[i,u]
        end
        B = c * c'
        B = B / (1 + Ld[v,v] + Ld[u,u] - Ld[u,v] - Ld[v,u])
        Ld = Ld - B

        R0 = n*Ld[v,v]
        for i = 1:n
            R0 = R0 + Ld[i,i]
        end
    #    print("Final Resistant:")
        #println(f1,n/R0)
        ans[dep+1] = n/R0
    end

#end

    #print("Final Resistant:")
#    println(R0)

    return ans
end
