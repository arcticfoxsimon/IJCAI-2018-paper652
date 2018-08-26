function ApproxiER(a, v, cand; ep=0.3, matrixConcConst=4.0, JLfac=20.0)

  if ep > 1
    warn("Calling sparsify with ep > 1 can produce a disconnected graph.")
  end

  #@time
  f = approxCholLap(a,tol=1e-5);

  n = size(a,1)
  k = round(Int, JLfac*log(n)) # number of dims for JL

  U = wtedEdgeVertexMat(a)
  m = size(U,1)
  R = randn(m,k)
  UR = U'*R;

  V = zeros(n,k)
  #@time begin
  for i in 1:k
    V[:,i] .= f(UR[:,i])
  end
#end
  #(ai,aj,av) = findnz(triu(a))
  #prs = zeros(av)
  #er = Dict{ Tuple{Int,Int}, Float64 }()

  er = zeros(n,1)
  for i = 1:n
      if cand[i] > 0
          er[i] =  (norm(V[v,:]-V[i,:])^2/k)
      end
      #prs[h] = av[h]* (norm(V[i,:]-V[j,:])^2/k)#min(1,av[h]* (norm(V[i,:]-V[j,:])^2/k) * matrixConcConst *log(n)/ep^2)
  end

  return er;

end
