module IntegerSmithNormalForm
using LinearAlgebra: I

export snf!, snf, snfWithoutTransform

"""
    (S,B,T) = negateRow(A,i)

Negate the i-th row of the integer matrix A such that B = SAT holds true
"""
function negateRow(A::AbstractMatrix{E}, i) where {E}
    m, n = size(A)
    S = Matrix{E}(I, m, m)
    T = Matrix{E}(I, n, n)
    S[i, i] = -one(E)
    B = S * A * T
    (S, B, T)
end


"""
    (S,B,T) = negateCol(A,i)

  Negate the i-th column of the integer matrix A such that B = SAT holds true
"""
function negateCol(A::AbstractMatrix, i)
    (S, B, T) = negateRow(A', i)
    (T', B', S')
end


"""
    (S,B,T) = swapRows(A,i,j)

  Negate the i-th and j-th row of the integer matrix A such that B = SAT holds true
"""
function swapRows(A::AbstractMatrix{E}, i, j) where {E}
    m, n = size(A)
    S = Matrix{E}(I, m, m)
    T = Matrix{E}(I, n, n)

    S[i, i] = zero(E)
    S[j, j] = zero(E)
    S[i, j] = one(E)
    S[j, i] = one(E)

    B = S * A * T
    (S, B, T)
end

"""
    (S,B,T) = swapCols(A,i,j)

  Negate the i-th and j-th column of the integer matrix A such that B = SAT holds true
"""
function swapCols(A::AbstractMatrix, i, j)
    (S, B, T) = swapRows(A', i, j)
    (T', B', S')
end


"""
    (S,B,T) = addRow(A,i,j,a)

  Add the j-th row a times to the i-th row of the integer matrix A such that B = SAT holds true
"""
function addRow(A::AbstractMatrix{E}, i, j, a::E) where {E}
    m, n = size(A)
    S = Matrix{E}(I, m, m)
    T = Matrix{E}(I, n, n)

    S[i, j] = a

    B = S * A * T
    (S, B, T)
end

"""
    (S,B,T) = addCol(A,i,j,a)

  Add the j-th column a times to the i-th column of the integer matrix A such that B = SAT holds true
"""
function addCol(A::AbstractMatrix{E}, i, j, a::E) where {E}
    (S, B, T) = addRow(A', i, j, a)
    (T', B', S')
end


"""
    (S,B,T) = snf!(A)

Computes the Smith normal form according to
https://www.charite.de/sysbio/people/hoppe/Diplomarbeit_Hoppe.pdf
Algorithmus 1

The matrices fulfill:
* B is a diagonal matrix
* B[i,i] >= 0 for all i
* B[i,i] divides all B[j,j] for all j > i
"""
function snf!(A::AbstractMatrix{E}) where {E<:Integer}
    m, n = size(A)
    S = Matrix{E}(I, m, m)
    T = Matrix{E}(I, n, n)

    for k = 1:(min(m, n))
        As = A[k:end, k:end]
        ms, ns = size(As)
        Ss = Matrix{E}(I, ms, ms)
        Ts = Matrix{E}(I, ns, ns)
        pivot = findnext(As .!= 0, CartesianIndex(1, 1))
        if !isnothing(pivot)
            #Step 2: Put pivot to (1,1)-position
            (S_, As, T_) = swapRows(As, 1, pivot[1])
            Ss = S_ * Ss
            Ts = Ts * T_
            (S_, As, T_) = swapCols(As, 1, pivot[2])
            Ss = S_ * Ss
            Ts = Ts * T_

            divideMatrix = false
            while !divideMatrix
                #check if ASub[1,1] divides all elements in the
                #column and row
                divideCol = false
                divideRow = false
                while !(divideRow & divideCol)
                    divideCol = true
                    #check if ASub[1,1] divides all elements in its column
                    if ms > 1
                        for i = 2:ms
                            #does not divide
                            if mod.(As[i, 1], As[1, 1]) != 0
                                divideCol = false
                                q = div(As[i, 1], As[1, 1])
                                #Add row 1 q-times to row i
                                (S_, As, T_) = addRow(As, i, 1, -q)
                                Ss = S_ * Ss
                                Ts = Ts * T_
                                #Make (1,i) the new pivot on position (1,1)
                                (S_, As, T_) = swapRows(As, 1, i)
                                Ss = S_ * Ss
                                Ts = Ts * T_
                                break
                            end
                        end
                        divideCol || break
                    end

                    divideRow = true
                    #check if ASub[1,1] divides all elements in its row
                    if ns > 1
                        for i = 1:ns
                            #does not divide
                            if mod.(As[1, i], As[1, 1]) != 0
                                divideRow = false
                                q = div(As[1, i], As[1, 1])
                                #Add col 1 q-times to col i
                                (S_, As, T_) = addCol(As, i, 1, -q)
                                Ss = S_ * Ss
                                Ts = Ts * T_
                                #Make (1,i) the new pivot on position (1,1)
                                (S_, As, T_) = swapCols(As, 1, i)
                                Ss = S_ * Ss
                                Ts = Ts * T_
                                break
                            end
                        end
                    end
                end

                #eliminate all entries in the pivot row and column
                if ms > 1
                    for i = 2:ms
                        q = div(As[i, 1], As[1, 1])
                        (S_, As, T_) = addRow(As, i, 1, -q)
                        Ss = S_ * Ss
                        Ts = Ts * T_
                    end
                end

                if ns > 1
                    for i = 2:ns
                        q = div(As[1, i], As[1, 1])
                        (S_, As, T_) = addCol(As, i, 1, -q)
                        Ss = S_ * Ss
                        Ts = Ts * T_
                    end
                end

                #check if the pivot element divides all elements of the matrix
                divideMatrix = mod.(As, As[1, 1]) == zeros(size(As))

                if !divideMatrix
                    #Add row of non-dividing index to pivot row

                    non_div = findnext(mod.(As, As[1, 1]) .!= zero(E), CartesianIndex(1, 1))
                    (S_, As, T_) = addRow(As, 1, non_div[1], 1)
                    Ss = S_ * Ss
                    Ts = Ts * T_
                end
            end


            S_ = Matrix{E}(I, m, m)
            T_ = Matrix{E}(I, n, n)
            S_[k:end, k:end] = Ss
            T_[k:end, k:end] = Ts
            S = S_ * S
            T = T * T_
            A[k:end, k:end] = As
        end
    end

    #make diagonal element positive
    for i = 1:min(m, n)
        if A[i, i] < 0
            (S_, A, T_) = negateRow(A, i)
            S = S_ * S
            T = T * T_
        end
    end
    return (S, A, T)
end

"""
    (S,B,T) = snf(A)

Computes the Smith normal form according to
https://www.charite.de/sysbio/people/hoppe/Diplomarbeit_Hoppe.pdf
Algorithmus 1

The matrices fulfill:
* B = SAT
* B is a diagonal matrix
* B[i,i] >= 0 for all i
* B[i,i] divides all B[j,j] for all j > i
"""
function snf(A)
    B = copy(A)
    return snf!(B)
end


"""
    B = snfWithoutTransform(A)

Computes the Smith normal form according to
https://www.charite.de/sysbio/people/hoppe/Diplomarbeit_Hoppe.pdf
Algorithmus 1

The matrix fulfills:
* There are matrices S,T with |det(S)| = |det(T)| = 1 such that B = SAT
* B is a diagonal matrix
* B[i,i] >= 0 for all i
* B[i,i] divides all B[j,j] for all j > i
"""
function snfWithoutTransform(A)
    (S, B, T) = snf(A)
    B
end

end # module
