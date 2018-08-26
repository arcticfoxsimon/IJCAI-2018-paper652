include("graph.jl")
include("lapl.jl")

function Plain_SM2(G, k) #v- the vertex chosed; G- graph; k- the number of edge added

    n = G.n

    cand = zeros(n,n) #in general situation; candidate edges are given as a parameter.

    L = lapl(G)

    ans = zeros(k+1);

    Ld = mppinv(L)
    R0 = 0
    for i = 1:n
        R0 = R0 + Ld[i,i]
    end
    ans[1] = n/R0
    #print("Original Resistant:")
    #println(R0)


#    @time begin
    for dep = 1:k
        #Ld = mppinv(L)
        Ld2 = Ld*Ld
        max = 0
        u1 = 0
        u2 = 0
        for i = 1: n
            for j = 1: n
                if cand[i,j] == 0
                    dR = Ld2[j,j] + Ld2[i,i] - Ld2[j,i] - Ld2[i,j]
                    dR = dR / (1 + Ld[j,j] + Ld[i,i] - Ld[j,i] - Ld[i,j])
                    if dR > max
                        max = dR
                        u1 = i
                        u2 = j
                    end
                end
            end
        end
        #ans[dep] = u
        cand[u1,u2] = 1
        L[u1,u2] = L[u1,u2]-1
        L[u2,u1] = L[u2,u1]-1
        L[u1,u1] = L[u1,u1]+1
        L[u2,u2] = L[u2,u2]+1

        c = zeros(n,1)
        for i = 1:n
            c[i] = Ld[i,u1] - Ld[i,u2]
        end
        B = c * c'
        B = B / (1 + Ld[u1,u1] + Ld[u2,u2] - Ld[u1,u2] - Ld[u2,u1])
        Ld = Ld - B

        R0 = 0
        for i = 1:n
            R0 = R0 + Ld[i,i]
        end
        #print("Final Resistant:")
        #println(f1,n/R0)
        ans[dep+1] = n/R0
    end

#end
#=
    Ldag = mppinv(L)
    R0 = 0
    for i = 1:n
        R0 = R0 + Ldag[i,i]
    end
    print("Final Resistant:")
    println(R0)
=#
    return ans
end
