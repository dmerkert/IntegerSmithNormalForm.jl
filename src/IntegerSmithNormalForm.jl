module IntegerSmithNormalForm

export SNF

"""
    (S,B,T) = negateRow(A,i)

  Negate the i-th row of the integer matrix A such that B = SAT holds true
"""
function negateRow{I}(A::Array{I,2},i)
  S = eye(I,size(A,1))
  T = eye(I,size(A,2))
  S[i,i] = -1
  B = S*A*T
  (S,B,T)
end


"""
    (S,B,T) = negateCol(A,i)

  Negate the i-th column of the integer matrix A such that B = SAT holds true
"""
function negateCol{I}(A::Array{I,2},i)
  (S,B,T) = negateRow(A',i)
  (T',B',S')
end


"""
    (S,B,T) = swapRows(A,i,j)

  Negate the i-th and j-th row of the integer matrix A such that B = SAT holds true
"""
function swapRows{I}(A::Array{I,2},i,j)
  S = eye(I,size(A,1))
  T = eye(I,size(A,2))

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
function swapCols{I}(A::Array{I,2},i,j)
  (S,B,T) = swapRows(A',i,j)
  (T',B',S')
end


"""
    (S,B,T) = addRow(A,i,j,a)

  Add the j-th row a times to the i-th row of the integer matrix A such that B = SAT holds true
"""
function addRow{I}(A::Array{I,2},i,j,a::I)
  S = eye(I,size(A,1))
  T = eye(I,size(A,2))

  S[i,j] = a

  B = S*A*T
  (S,B,T)
end

"""
    (S,B,T) = addRow(A,i,j,a)

  Add the j-th column a times to the i-th column of the integer matrix A such that B = SAT holds true
"""

function addCol{I}(A::Array{I,2},i,j,a::I)
  (S,B,T) = addRow(A',i,j,a)
  (T',B',S')
end


"""
    (S,E,T) = SNF(A)

  Computes the Smith normal form according to
  https://www.charite.de/sysbio/people/hoppe/Diplomarbeit_Hoppe.pdf 
  Algorithmus 1
"""
function SNF{I<:Integer}(A::Array{I,2})
  d1 = size(A,1)
  d2 = size(A,2)

  E = A
  S = eye(I,d1)
  T = eye(I,d2)

  for submatrixIndex in 1:min(d1,d2)
    #Step 1: Find pivot unequal to zero (arbitrary element in sub-matrix)
    pivot = (0,0)
    for i in submatrixIndex:d1
      for j in submatrixIndex:d2
        if A[i,j] != 0
          pivot = (i,j)
          break
        end
      end
    end

    #Step 2: Put pivot to (submatrixIndex,submatrixIndex)-position
    (S_,E,T_) = swapRows(E,submatrixIndex,pivot[1])
    S = S_ * S
    T = T * T_
    (S_,E,T_) = swapCols(E,submatrixIndex,pivot[2])
    S = S_ * S
    T = T * T_

    #check if E[submatrixIndex,submatrixIndex] divides all elements in the
    #column and row
    while
      sum(
          mod(
              E[submatrixIndex:end,submatrixIndex],
              E[submatrixIndex,submatrixIndex]
             )
         ) > 0 |
      sum(
          mod(
              E[submatrixIndex,submatrixIndex:end],
              E[submatrixIndex,submatrixIndex]
             )
         ) > 0

    end








  end
end

end # module
