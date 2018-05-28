using IntegerSmithNormalForm: negateRow, negateCol, swapRows, swapCols, addRow,
addCol, SNF!, SNF, SNFWithoutTransform
using Base.Test

@testset "Negating rows and columns" begin
    A = [[1 2];[3 4]]
    B = [[1 2 3];[4 5 6];[7 8 9]]
    C = [[1 2];[3 4];[5 6]]

    (S,Tmp,T) = negateRow(A,1)
    @test Tmp == [[-1 -2];[3 4]]
    @test S*A*T == Tmp
    @inferred negateRow(A,1)

    (S,Tmp,T) = negateRow(A,2)
    @test Tmp == [[1 2];[-3 -4]]
    @test S*A*T == Tmp
    @inferred negateRow(A,2)

    (S,Tmp,T) = negateRow(B,3)
    @test Tmp == [[1 2 3];[4 5 6];[-7 -8 -9]]
    @test S*B*T == Tmp
    @inferred negateRow(B,3)

    (S,Tmp,T) = negateRow(C,2)
    @test Tmp == [[1 2];[-3 -4];[5 6]]
    @test S*C*T == Tmp
    @inferred negateRow(C,2)


    (S,Tmp,T) = negateCol(A,1)
    @test Tmp == [[-1 2];[-3 4]]
    @test S*A*T == Tmp
    @inferred negateCol(A,1)

    (S,Tmp,T) = negateCol(A,2)
    @test Tmp == [[1 -2];[3 -4]]
    @test S*A*T == Tmp
    @inferred negateCol(A,2)

    (S,Tmp,T) = negateCol(B,3)
    @test Tmp == [[1 2 -3];[4 5 -6];[7 8 -9]]
    @test S*B*T == Tmp
    @inferred negateCol(B,3)

    (S,Tmp,T) = negateCol(C,2)
    @test Tmp == [[1 -2];[3 -4];[5 -6]]
    @test S*C*T == Tmp
    @inferred negateCol(C,2)
end

@testset "Swapping rows and columns" begin
    A = [[1 2];[3 4];[5 6]]

    (S,B,T) = swapRows(A,1,2)
    @test B == [[3 4];[1 2];[5 6]]
    @test S*A*T == B
    @inferred swapRows(A,1,2)

    (S,B,T) = swapRows(A,3,2)
    @test B == [[1 2];[5 6];[3 4]]
    @test S*A*T == B
    @inferred swapRows(A,3,2)

    A = [[1 2 3];[4 5 6]]
    (S,B,T) = swapCols(A,1,2)
    @test B == [[2 1 3];[5 4 6]]
    @test S*A*T == B
    @inferred swapCols(A,1,2)

    (S,B,T) = swapCols(A,3,2)
    @test B == [[1 3 2];[4 6 5]]
    @test S*A*T == B
    @inferred swapCols(A,3,2)
end

@testset "Adding rows and columns" begin
    A = [[1 2];[3 4];[5 6]]

    (S,B,T) = addRow(A,1,2,10)
    @test B == [[31 42];[3 4];[5 6]]
    @test S*A*T == B
    @inferred addRow(A,1,2,10)

    A = [[1 2 3];[4 5 6]]

    (S,B,T) = addCol(A,2,3,10)
    @test B == [[1 32 3];[4 65 6]]
    @test S*A*T == B
    @inferred addCol(A,2,3,10)

end

function getMatrix(
                   k1,k2,k3, a
                  )

  b = 32

  if a == 0
    return [b k1 k2;
            0 b k3;
            0 0 b]
  elseif a == 1
    return [b 0 k2;
            k1 b k3;
            0 0 b]
  elseif a == 2
    return [b k1 0;
            0 b k3;
            k2 0 b]
  elseif a == 3
    return [b k1 k2;
            0 b 0;
            0 k3 b]
  elseif a == 4
    return [b 0 0;
            k1 b k3;
            k2 0 b]
  elseif a == 5
    return [b 0 k2;
            k1 b 0;
            0 k3 b]
  elseif a == 6
    return [b k1 0;
            0 b 0;
            k2 k3 b]
  elseif a == 7
    return [b 0 0;
            k1 b 0;
            k2 k3 b]
  end
end

function testSNF(A)
    B = copy(A)
    @inferred SNF!(B)
    B = copy(A)
    (S,B,T) = SNF!(B)
    @test B == S*A*T
    @test isdiag(B)
    d = diag(B)
    for i in 1:length(d)
        @test d[i] >= 0
        if d[i] != 0
            for j in i:length(d)
                @test mod(d[j],d[i]) == 0
            end
        end
    end
end

@testset "Test Smith Normal Form" begin

    testSNF([16 -16 -16; -16 16 -2; 9 -9 16])
    testSNF([[1 2 3];[4 5 6];[1 4 9]])
    testSNF([[1 2];[4 5];[1 4]])

    aRange = 0:7
    k1Range = [-24:1:24..., 99,100,101]

    simulations = length(aRange)*length(k1Range)^3

    for a in aRange
        for k1 in k1Range
            for k2 in k1Range
                for k3 in k1Range
                    M = getMatrix(k1,k2,k3,a)
                    testSNF(M)
                end
            end
        end
    end



    (S,B,T) = SNF([[1 2 3];[4 5 6];[1 4 9]])
    @inferred SNF([[1 2 3];[4 5 6];[1 4 9]])
    @test B == S*[[1 2 3];[4 5 6];[1 4 9]]*T
    @test isdiag(B)
    d = diag(B)
    for i in 1:length(d)
        @test d[i] >= 0
        for j in i:length(d)
            @test mod(d[j],d[i]) == 0
        end
    end

    B = SNFWithoutTransform([[1 2 3];[4 5 6];[1 4 9]])
    @inferred SNFWithoutTransform([[1 2 3];[4 5 6];[1 4 9]])
    @test isdiag(B)
    d = diag(B)
    for i in 1:length(d)
        @test d[i] >= 0
        for j in i:length(d)
            @test mod(d[j],d[i]) == 0
        end
    end
end
