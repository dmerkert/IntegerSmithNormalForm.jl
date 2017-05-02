module IntegerSmithNormalForm

#export SNF

"""
    (S,B,T) = negateRow(A,i)

  Negate the i-th row of the integer matrix A such that B = SAT holds true
"""
function negateRow(A::Array{Int,2},i::Int)
  S = eye(Int,size(A,1))
  T = eye(Int,size(A,2))
  S[i,i] = -1
  B = S*A*T
  (S,B,T)
end


"""
    (S,B,T) = negateCol(A,i)

  Negate the i-th column of the integer matrix A such that B = SAT holds true
"""
function negateCol(A::Array{Int,2},i::Int)
  (S,B,T) = negateRow(A',i)
  (T',B',S')
end


"""
    (S,B,T) = swapRows(A,i,j)

  Negate the i-th and j-th row of the integer matrix A such that B = SAT holds true
"""
function swapRows(A::Array{Int,2},i::Int,j::Int)
  S = eye(Int,size(A,1))
  T = eye(Int,size(A,2))

  S[i,i] = 0
  S[j,j] = 0
  S[i,j] = 1
  S[j,i] = 1

  B = S*A*T
  (S,B,T)
end

"""
    (S,B,T) = swapCols(A,i,j)

  Negate the i-th and j-th column of the integer matrix A such that B = SAT holds true
"""
function swapCols(A::Array{Int,2},i::Int,j::Int)
  (S,B,T) = swapRows(A',i,j)
  (T',B',S')
end


"""
    (S,B,T) = addRow(A,i,j,a)

  Add the j-th row a times to the i-th row of the integer matrix A such that B = SAT holds true
"""
function addRow(A::Array{Int,2},i::Int,j::Int,a::Int)
  S = eye(Int,size(A,1))
  T = eye(Int,size(A,2))

  S[i,j] = a

  B = S*A*T
  (S,B,T)
end

"""
    (S,B,T) = addRow(A,i,j,a)

  Add the j-th column a times to the i-th column of the integer matrix A such that B = SAT holds true
"""

function addCol(A::Array{Int,2},i::Int,j::Int,a::Int)
  (S,B,T) = addRow(A',i,j,a)
  (T',B',S')
end

#
#function SNF(M::Array{Int,2})
#  #Computes the Smith normal form according to
#  #https://www.charite.de/sysbio/people/hoppe/Diplomarbeit_Hoppe.pdf 
#  #Algorithmus 1
#  size(M,1) == size(M,2) || error("SNF not implemented for rectangular
#                                  matrices")
#  d = size(M,1)
#  d = 1 && return(M,eye(1),eye(1))
#
#  A = M
#  S = eye(d)
#  T = eye(d)
#  #Ziel: A = S*M*T
#
#  #Find pivot unequal to zero (arbitrary element in matrix)
#  pivot = (0,0)
#  for i in 1:d
#    for j in 1:d
#      if A(i,j) != 0
#        pivot = (i,j)
#        break
#      end
#    end
#  end
#
#  #Put pivot to (1,1)-position
#  (A,S_tmp) = SNF_swapRows(1,pivot[1])
#  (A,T_tmp) = SNF_swapCols(1,pivot[2])
#  S = S_tmp*S
#  T = T*T_tmp
#
#end

end # module
