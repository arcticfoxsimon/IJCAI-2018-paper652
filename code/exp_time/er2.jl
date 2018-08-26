function ApproxiER2(a, v, cand; ep=0.3, matrixConcConst=4.0, JLfac=20.0)

  if ep > 1
    warn("Calling sparsify with ep > 1 can produce a disconnected graph.")
  end

  #@time
  f = approxCholLap(a,tol=1e-5);

  n = size(a,1)
  k = round(Int, JLfac*log(n)) # number of dims for JL

  U = wtedEdgeVertexMat(a)
  m = size(U,1)
  er = zeros(n,1)
  #R = randn(m,k)
  #UR = U'*R;

  #V = zeros(n,k)
  #@time begin
  for i in 1:k
      q=randn(m, 1)
      q = U'*q
      z=zeros(n,1)
      z.=f(q[:])
      for j=1: n
#          u=edges[e,1]
#          v=edges[e,2]
          er[j]=er[j]+(z[j]-z[v])^2
      end
    #V[:,i] .= f(UR[:,i])
  end
#end
  #(ai,aj,av) = findnz(triu(a))
  #prs = zeros(av)
  #er = Dict{ Tuple{Int,Int}, Float64 }()

#  er = zeros(n,1)
  for i = 1:n
      if cand[i] == 0
          er[i] = 0
      end
      er[i] = er[i]/k
      #prs[h] = av[h]* (norm(V[i,:]-V[j,:])^2/k)#min(1,av[h]* (norm(V[i,:]-V[j,:])^2/k) * matrixConcConst *log(n)/ep^2)
  end

  return er;

end
