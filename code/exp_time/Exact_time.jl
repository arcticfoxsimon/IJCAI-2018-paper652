include("graph.jl")
include("lapl.jl")

function Exact_time(v, G, k) #v- the vertex chosed; G- graph; k- the number of edge added

    n = G.n;
    cand = zeros(n) #in general situation; candidate edges are given as a parameter.
    L = lapl(G)
    for i = 1 : n
        if L[v,i] == 0
            cand[i] = 1
        end
    end    # candidate edges are edges which are not existant in the original graph

    ans = zeros(Int64, k);
#=
    Ldag = mppinv(L)
    R0 = n*Ldag[v,v]
    for i = 1:n
        R0 = R0 + Ldag[i,i]
    end
=#
#    print("Original Resistant:")
#    println(R0)


#    @time begin
start_time = time()

    Ld = mppinv(L)

    #fname = "experiment\\exp2.txt"
    #f1 = open(fname,"w")

    #println(f1,"Information")

#    R0 = n*Ld[v,v]
#    for i = 1:n
#        R0 = R0 + Ld[i,i]
#    end
#    print("Final Resistant:")
    #println(f1,n/R0)

#    ans[1] = n/R0

    for dep = 1:k
        #Ld = mppinv(L)
        Ld2 = Ld*Ld
        max = 0
        u = 0
        for i = 1: n
            if cand[i] == 1
                dR = n*(Ld[v,v] - Ld[v,i])*(Ld[v,v]-Ld[i,v])
                dR = dR + Ld2[v,v] + Ld2[i,i] - Ld2[v,i] - Ld2[i,v]
                dR = dR / (1 + Ld[v,v] + Ld[i,i] - Ld[v,i] - Ld[i,v])
                if dR > max
                    max = dR;
                    u = i
                end
            end
        end
        ans[dep] = u
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

#        R0 = n*Ld[v,v]
#        for i = 1:n
#            R0 = R0 + Ld[i,i]
#        end
    #    print("Final Resistant:")
        #println(f1,n/R0)
#        ans[dep+1] = n/R0
    end

    time_elapse = time() - start_time
#    println("time: ",time_elapse)

        Ldag = mppinv(L)
        R0 = n*Ldag[v,v]
        for i = 1:n
            R0 = R0 + Ldag[i,i]
        end
        #print("Final Resistant:")
#        println("Final centrality:", n/R0)

        res = zeros(2)
        res[1] = time_elapse
        res[2] = n/R0

        return res
end
