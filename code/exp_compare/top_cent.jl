include("graph.jl")
include("lapl.jl")

function Top_cent(v, G, k) #v- the vertex chosed; G- graph; k- the number of edge added

    n = G.n;
    cand = zeros(n) #in general situation; candidate edges are given as a parameter.
    L = lapl(G)
    for i = 1 : n
        if L[v,i] == 0
            cand[i] = 1
        end
    end    # candidate edges are edges which are not existant in the original graph

    ans = zeros(k+1);

    Ldag = mppinv(L)
    vcent = zeros(n)
    for i = 1:n
        vcent[i] = n*Ldag[i,i]
    end
    p = sortperm(vcent)

#    Ldag = mppinv(L)
#    R0 = n*Ldag[v,v]
#    for i = 1:n
#        R0 = R0 + Ldag[i,i]
#    end
#    print("Original Resistant:")
#    println(R0)


#    @time begin
    Ld = mppinv(L)

    #fname = "experiment\\exp5.txt"
    #f1 = open(fname,"w")

    #println(f1,"Information")

    R0 = n*Ld[v,v]
    for i = 1:n
        R0 = R0 + Ld[i,i]
    end
#    print("Final Resistant:")
    #println(f1,n/R0)
    ans[1] = n/R0

    t = 1
    for dep = 1:k
        #Ld = mppinv(L)

        u = p[t]
        while cand[u] == 0
            t = t + 1
            u = p[t]
        end

        #ans[dep] = u
        cand[u] = 0
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
    #    println(f1,n/R0)
        ans[dep+1] = n/R0
    end

#close(f1)
#end

    return ans
end
