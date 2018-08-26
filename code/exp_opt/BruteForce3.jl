include("graph.jl")
include("lapl.jl")

function DFS(Ldag,cand,cur,k,min,dep)
    Ld = Ldag
    n = size(Ld,1)
    if dep > k
        R = 0
        for i = 1:n
            R = R + Ld[i,i]
        end

        #println(R)

        if min == 0
            return R
        end
        if R < min
            return R
        end
        return min
    end

    for i = cur:n*n
        if cand[i] == 0
            c = zeros(n,1)
            x = div(i-1,n) + 1
            y = (i-1)%n + 1
            for j = 1:n
                c[j] = Ld[j,x] - Ld[j,y]
            end
            B = c * c'
            B = B / (1 + Ld[x,x] + Ld[y,y] - Ld[x,y] - Ld[y,x])
            Ld = Ld - B

            min = DFS(Ld,cand,i+1,k,min,dep+1)

            Ld = Ld + B
        end
    end

    return min
end

function BruteF3(G, k)

    n = G.n

    cand = zeros(n*n) #in general situation; candidate edges are given as a parameter.

    L = lapl(G)

    ans = zeros(Int64,k);

    Ld = mppinv(L)
    R0 = 0
    for i = 1:n
        R0 = R0 + Ld[i,i]
    end
    print("Original Resistant:")
    println(R0)

    min = DFS(Ld,cand,1,k,0,1)

        #println(min)

        #Ldag = mppinv(L)
        #R0 = n*Ldag[v,v]
        #for i = 1:n
        #    R0 = R0 + Ldag[i,i]
        #end
        #print("Final Resistant:")
        #println(R0)

    return min
end
