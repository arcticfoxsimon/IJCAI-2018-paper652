include("graph.jl")
include("lapl.jl")
include("BruteForce3.jl")
include("appro_SM3.jl")
include("plain_SM2.jl")
include("random2.jl")

G = read_file("data/windsufer.edges")

n = G.n;
cand = zeros(n) #in general situation; candidate edges are given as a parameter.
L = lapl(G)

Ldag = mppinv(L)
vcent = zeros(n)
for i = 1:n
    vcent[i] = n*Ldag[i,i]
end
p = sortperm(vcent)

k = 6
ans1 = zeros(k)
ans2 = zeros(k)
ans3 = zeros(k)
ans4 = zeros(k)

tmp1 = zeros(k)
tmp2 = zeros(k)
tmp3 = zeros(k)
tmp4 = zeros(k)

snum = 20
#group = 3
for dep = 1: snum
    group = ceil(Int64, dep / 5)
    t = rand(round(Int64,n/4*(group-1))+1: round(Int64,n/4*group))
    #v = rand(1:n)

    #println(v)
    println(p[t])

    v = p[t]

    for i = 1:k
        tmp1[i] = Appro_SM3(v,G,i)
        tmp2[i] = Plain_SM2(v,G,i)
        tmp3[i] = BruteF3(v,G,i)
        tmp4[i] = Random2(v,G,i)
    end

    ans1 = ans1 + tmp1
    ans2 = ans2 + tmp2
    ans3 = ans3 + tmp3
    ans4 = ans4 + tmp4

end

fname = "exp2\\exp1.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k
    println(f1,ans1[i]/snum)
end
close(f1)

fname = "exp2\\exp2.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k
    println(f1,ans2[i]/snum)
end
close(f1)


fname = "exp2\\exp3.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k
    println(f1,ans3[i]/snum)
end
close(f1)


fname = "exp2\\exp4.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k
    println(f1,ans4[i]/snum)
end
close(f1)

println("End")
