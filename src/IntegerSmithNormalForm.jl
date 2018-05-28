module IntegerSmithNormalForm

export SNF!, SNF, SNFWithoutTransform

"""
    (S,B,T) = negateRow(A,i)

  Negate the i-th row of the integer matrix A such that B = SAT holds true
"""
function negateRow(A::AbstractArray{I,2},i) where I
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
function negateCol(A::AbstractArray{I,2},i) where I
  (S,B,T) = negateRow(A',i)
  (T',B',S')
end


"""
    (S,B,T) = swapRows(A,i,j)

  Negate the i-th and j-th row of the integer matrix A such that B = SAT holds true
"""
function swapRows(A::AbstractArray{I,2},i,j) where I
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
function swapCols(A::AbstractArray{I,2},i,j) where I
  (S,B,T) = swapRows(A',i,j)
  (T',B',S')
end


"""
    (S,B,T) = addRow(A,i,j,a)

  Add the j-th row a times to the i-th row of the integer matrix A such that B = SAT holds true
"""
function addRow(A::AbstractArray{I,2},i,j,a::I) where I
  S = eye(I,size(A,1))
  T = eye(I,size(A,2))

  S[i,j] = a

  B = S*A*T
  (S,B,T)
end

"""
    (S,B,T) = addCol(A,i,j,a)

  Add the j-th column a times to the i-th column of the integer matrix A such that B = SAT holds true
"""
function addCol(A::AbstractArray{I,2},i,j,a::I) where I
  (S,B,T) = addRow(A',i,j,a)
  (T',B',S')
end


"""
    (S,B,T) = SNF!(A)

Computes the Smith normal form according to
https://www.charite.de/sysbio/people/hoppe/Diplomarbeit_Hoppe.pdf 
Algorithmus 1

The matrices fulfill:
* B is a diagonal matrix
* B[i,i] >= 0 for all i
* B[i,i] divides all B[j,j] for all j > i
"""
function SNF!(A::AbstractArray{I,2}) where {I<:Integer}
  d1 = size(A,1)
  d2 = size(A,2)

  S = eye(I,d1)
  T = eye(I,d2)

  for subMatrixInd in 1:(min(d1,d2))
    ASub = A[subMatrixInd:end,subMatrixInd:end]
    d1Sub = size(ASub,1)
    d2Sub = size(ASub,2)

    SSub = eye(I,d1Sub)
    TSub = eye(I,d2Sub)

    nextPivot = findnext(ASub,1)

    if nextPivot != 0

      pivot = ind2sub(size(ASub),nextPivot)

      #Step 2: Put pivot to (1,1)-position
      (S_,ASub,T_) = swapRows(ASub,1,pivot[1])
      SSub = S_ * SSub
      TSub = TSub * T_
      (S_,ASub,T_) = swapCols(ASub,1,pivot[2])
      SSub = S_ * SSub
      TSub = TSub * T_

      divideMatrix = false
      while !divideMatrix
        #check if ASub[1,1] divides all elements in the
        #column and row
        divideCol = false
        divideRow = false
        while !(divideRow & divideCol)
          divideCol=true
          #check if ASub[1,1] divides all elements in its column
          if d1Sub > 1
            for i in 2:d1Sub
              #does not divide
              if mod.(ASub[i,1],ASub[1,1]) != 0
                divideCol=false
                q = div(ASub[i,1],ASub[1,1])
                #Add row 1 q-times to row i
                (S_,ASub,T_) = addRow(ASub,i,1,-q)
                SSub = S_ * SSub
                TSub = TSub * T_
                #Make (1,i) the new pivot on position (1,1)
                (S_,ASub,T_) = swapRows(ASub,1,i)
                SSub = S_ * SSub
                TSub = TSub * T_
                break
              end
            end
            divideCol || break
          end

          divideRow=true
          #check if ASub[1,1] divides all elements in its row
          if d2Sub > 1
            for i in 1:d2Sub
              #does not divide
              if mod.(ASub[1,i],ASub[1,1]) != 0
                divideRow=false
                q = div(ASub[1,i],ASub[1,1])
                #Add col 1 q-times to col i
                (S_,ASub,T_) = addCol(ASub,i,1,-q)
                SSub = S_ * SSub
                TSub = TSub * T_
                #Make (1,i) the new pivot on position (1,1)
                (S_,ASub,T_) = swapCols(ASub,1,i)
                SSub = S_ * SSub
                TSub = TSub * T_
                break
              end
            end
          end
        end

        #eliminate all entries in the pivot row and column
        if d1Sub > 1
          for i = 2:d1Sub
            q = div(ASub[i,1],ASub[1,1])
            (S_,ASub,T_) = addRow(ASub,i,1,-q)
            SSub = S_ * SSub
            TSub = TSub * T_
          end
        end

        if d2Sub > 1
          for i = 2:d2Sub
            q = div(ASub[1,i],ASub[1,1])
            (S_,ASub,T_) = addCol(ASub,i,1,-q)
            SSub = S_ * SSub
            TSub = TSub * T_
          end
        end

        #check if the pivot element divides all elements of the matrix
        divideMatrix = mod.(ASub,ASub[1,1]) == zeros(size(ASub))

        if !divideMatrix
          #Add row of non-dividing index to pivot row
          nonDividing = ind2sub(size(ASub),findnext(mod.(ASub,ASub[1,1]),1))
          (S_,ASub,T_) = addRow(ASub,1,nonDividing[1],1)
          SSub = S_ * SSub
          TSub = TSub * T_
        end
      end


      S_ = eye(I,d1)
      T_ = eye(I,d2)
      S_[subMatrixInd:end,subMatrixInd:end] = SSub
      T_[subMatrixInd:end,subMatrixInd:end] = TSub
      S = S_ * S
      T = T * T_
      A[subMatrixInd:end,subMatrixInd:end] = ASub
    end
  end

  #make diagonal element positive
  for i in 1:min(d1,d2)
    if A[i,i] < 0
      (S_,A,T_) = negateRow(A,i)
      S = S_ * S
      T = T * T_
    end
  end


  (S,A,T)
end

"""
    (S,B,T) = SNF(A)

Computes the Smith normal form according to
https://www.charite.de/sysbio/people/hoppe/Diplomarbeit_Hoppe.pdf 
Algorithmus 1

The matrices fulfill:
* B = SAT
* B is a diagonal matrix
* B[i,i] >= 0 for all i
* B[i,i] divides all B[j,j] for all j > i
"""
function SNF(A)
  B = copy(A)
  SNF!(B)
end


"""
    B = SNFWithoutTransform(A)

Computes the Smith normal form according to
https://www.charite.de/sysbio/people/hoppe/Diplomarbeit_Hoppe.pdf 
Algorithmus 1

The matrix fulfills:
* There are matrices S,T with |det(S)| = |det(T)| = 1 such that B = SAT
* B is a diagonal matrix
* B[i,i] >= 0 for all i
* B[i,i] divides all B[j,j] for all j > i
"""
function SNFWithoutTransform(A)
  (S,B,T) = SNF(A)
  B
end

end # module
