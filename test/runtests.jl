using IntegerSmithNormalForm: negateRow, negateCol, swapRows, swapCols, addRow, addCol
using Base.Test

@testset "Negating rows and columns" begin
    A = [[1 2];[3 4]]
    B = [[1 2 3];[4 5 6];[7 8 9]]
    C = [[1 2];[3 4];[5 6]]

    (S,Tmp,T) = negateRow(A,1)
    @test Tmp == [[-1 -2];[3 4]]
    @test S*A*T == Tmp
    @test typeof(S) ==   typeof(A)
    @test typeof(Tmp) == typeof(A)
    @test typeof(T) ==   typeof(A)

    (S,Tmp,T) = negateRow(A,2)
    @test Tmp == [[1 2];[-3 -4]]
    @test S*A*T == Tmp

    (S,Tmp,T) = negateRow(B,3)
    @test Tmp == [[1 2 3];[4 5 6];[-7 -8 -9]]
    @test S*B*T == Tmp

    (S,Tmp,T) = negateRow(C,2)
    @test Tmp == [[1 2];[-3 -4];[5 6]]
    @test S*C*T == Tmp


    (S,Tmp,T) = negateCol(A,1)
    @test Tmp == [[-1 2];[-3 4]]
    @test S*A*T == Tmp
    @test typeof(S) ==   typeof(A)
    @test typeof(Tmp) == typeof(A)
    @test typeof(T) ==   typeof(A)

    (S,Tmp,T) = negateCol(A,2)
    @test Tmp == [[1 -2];[3 -4]]
    @test S*A*T == Tmp

    (S,Tmp,T) = negateCol(B,3)
    @test Tmp == [[1 2 -3];[4 5 -6];[7 8 -9]]
    @test S*B*T == Tmp

    (S,Tmp,T) = negateCol(C,2)
    @test Tmp == [[1 -2];[3 -4];[5 -6]]
    @test S*C*T == Tmp
end

@testset "Swapping rows and columns" begin
    A = [[1 2];[3 4];[5 6]]

    (S,B,T) = swapRows(A,1,2)
    @test B == [[3 4];[1 2];[5 6]]
    @test S*A*T == B
    @test typeof(S) == typeof(A)
    @test typeof(B) == typeof(A)
    @test typeof(T) == typeof(A)

    (S,B,T) = swapRows(A,3,2)
    @test B == [[1 2];[5 6];[3 4]]
    @test S*A*T == B

    A = [[1 2 3];[4 5 6]]
    (S,B,T) = swapCols(A,1,2)
    @test B == [[2 1 3];[5 4 6]]
    @test S*A*T == B
    @test typeof(S) == typeof(A)
    @test typeof(B) == typeof(A)
    @test typeof(T) == typeof(A)

    (S,B,T) = swapCols(A,3,2)
    @test B == [[1 3 2];[4 6 5]]
    @test S*A*T == B
end

@testset "Adding rows and columns" begin
    A = [[1 2];[3 4];[5 6]]

    (S,B,T) = addRow(A,1,2,10)
    @test B == [[31 42];[3 4];[5 6]]
    @test S*A*T == B
    @test typeof(S) == typeof(A)
    @test typeof(B) == typeof(A)
    @test typeof(T) == typeof(A)

    A = [[1 2 3];[4 5 6]]

    (S,B,T) = addCol(A,2,3,10)
    @test B == [[1 32 3];[4 65 6]]
    @test S*A*T == B
    @test typeof(S) == typeof(A)
    @test typeof(B) == typeof(A)
    @test typeof(T) == typeof(A)

end
