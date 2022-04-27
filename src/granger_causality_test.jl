"""
   granger_causality_test(model, u; α = 0.01, verbose::Bool = true)
"""
function granger_causality_test(model, u; α=0.01, verbose::Bool=true)
   nChannels, _, order = size(model.A)
   Z = get_Z(u, order)
   gamma = Z * Z'
   b = reshape(model.A, nChannels * nChannels * order)

   n, m = size(model.pf)
   Va = zeros(n, m)
   Tr = zeros(n, m)
   CO = zeros(n, m)
   pValue = zeros(n, m)
   for i in 1:n
      for j in 1:n
         if i != j
            CO[i, j] = 1
            Tr[i, j], Va[i, j], v, th, pValue[i, j] = grangt(CO, b, gamma, model.pf, 1 - α)
            CO[i, j] = 0
         end
      end
   end
   if verbose
      println("                  GRANGER CAUSALITY TEST")
      println("======================================================================")
      println("Connectivity matrix:")
      Tr[I(nChannels)] .= NaN
      show(stdout, "text/plain", Tr)
      println()

      println("Granger causality test p-values:")
      pValue[I(nChannels)] .= NaN
      show(stdout, "text/plain", pValue)
      println()
      println("======================================================================")
   end
   return Tr, pValue
end

function grangt(CO, b, G, SU, significance)
   n, m = size(CO)
   lb = length(b)
   p0 = round(Int, lb / (m * n))
   Ct = repeat(reshape(CO, m * n), p0)

   v = sum(Int, Ct)
   C = zeros(v, m * n * p0)
   for i = 1:v
      p, q = findmax(Ct)
      C[i, q] = 1
      Ct[q] = 0
   end

   value = (C * b)' * inv(C * kron(inv(G), SU) * C') * C * b

   test_dist = Chisq(v)
   th = invlogcdf(test_dist, log(significance))
   y = value >= th

   pValue = 1 - cdf(test_dist, value)
   return y, value, v, th, pValue
end

"""
   instantaneous_granger_causality_test(model, u; α = 0.01, verbose::Bool = true)
"""
function instantaneous_granger_causality_test(model, u; α=0.01, verbose::Bool=true)
   nChannels, _, IP = size(model.A)
   N = size(u, 1)
   n, m = size(model.pf)
   Va = zeros(n, m)
   Tr = zeros(n, m)
   CO = zeros(n, m)
   pValue = zeros(n, m)
   for i in 1:n
      for j in 1:n
         if i > j
            CO[i, j] = 1
            Tr[i, j], Va[i, j], v, th, pValue[i, j] = instata(CO, model.pf, N, 1 - α)
            CO[i, j] = 0
            Tr[j, i] = Tr[i, j]
            pValue[j, i] = pValue[i, j]
         end
      end
   end

   if verbose
      println("                  INSTANTANEOUS GRANGER CAUSALITY TEST")
      println("======================================================================")
      println("Instantaneous connectivity matrix:")
      Tr[I(nChannels)] .= NaN
      display((Tr))
      println()
      println("Instantaneous Granger causality test p-values:")
      pValue[I(nChannels)] .= NaN
      display((pValue))
      println()
      println("======================================================================")
   end
   return Tr, pValue
end

function instata(CO, pf, N, significance)
   si = vech(pf)
   CO = LowerTriangular(CO)
   C = vech(CO)
   ln = size(pf, 1)
   D = pinv(dmatrix(ln))
   temp = 1 / (2 .* C' * D * kron(pf, pf) * D' * C)
   value = (N*(C' * si)'*temp*(C'*si))[1]
   v = 1
   test_dist = Chisq(v)
   th = invlogcdf(test_dist, log(significance))
   y = value >= th
   pValue = 1 - cdf(test_dist, value)
   return y, value, v, th, pValue
end