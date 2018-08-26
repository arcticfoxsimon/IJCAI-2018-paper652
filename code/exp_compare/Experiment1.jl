include("graph.jl")
include("plain_SM2.jl")
include("random2.jl")

G = read_file("data/euroroad.edges")

n = G.n;

k = 20
ans1 = zeros(k+1)
ans2 = zeros(k+1)

ans1 = Plain_SM2(G,k)
ans2 = Random2(G,k)


fname = "exp4/exp1.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k+1
    println(f1,ans1[i])
end
close(f1)

fname = "exp4/exp2.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k+1
    println(f1,ans2[i])
end
close(f1)

#=
snum = 10
#group = 3
for dep = 1: snum
    #group = ceil(Int64, dep / 5)
    #t = rand(round(Int64,n/4*(group-1))+1: round(Int64,n/4*group))
    v = rand(1:n)

    println(v)
    #println(p[t])

    #v = p[t]

    tmp1 = Appro_SM2(v,G,k)
    tmp2 = Plain_SM(v,G,k)
    tmp3 = Top_cent(v,G,k)
    tmp4 = Top_degree(v,G,k)
    tmp5 = Random(v,G,k)


    ans1 = ans1 + tmp1
    ans2 = ans2 + tmp2
    ans3 = ans3 + tmp3
    ans4 = ans4 + tmp4
    ans5 = ans5 + tmp5

end

fname = "exp3\\exp1.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k+1
    println(f1,ans1[i]/snum)
end
close(f1)

fname = "exp3\\exp2.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k+1
    println(f1,ans2[i]/snum)
end
close(f1)


fname = "exp3\\exp3.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k+1
    println(f1,ans3[i]/snum)
end
close(f1)


fname = "exp3\\exp4.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k+1
    println(f1,ans4[i]/snum)
end
close(f1)

fname = "exp3\\exp5.txt"
f1 = open(fname,"w")
println(f1,"Information")

for i = 1 : k+1
    println(f1,ans5[i]/snum)
end
close(f1)

println("End")
=#
