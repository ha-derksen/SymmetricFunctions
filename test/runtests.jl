using SymmetricFunctions
using Test

@testset "SymmetricFunctions.jl" begin

    f=3*s_(3,1,1)+s_(5)*7-s_(3,2)
    g=2*p_(3,2)-p_(4,1)*3
    h=4*c_(2,2)-c_(3,1)*5


    @testset "SymP,SymS,SymC,Sym" begin
        @test typeof(f)==SymS
        @test typeof(g)==SymP
        @test typeof(h)==SymC
        @test f isa Sym
        @test coeff(f,[5])==7
        @test coeff(g,[4,1])==-3
        @test coeff(h,[2,2])==4
    end

    @testset "conversions" begin
        @test SymS(f)==f
        @test SymP(g)==g
        @test SymC(h)==h
        @test SymS(SymP(f))==f
        @test SymS(SymC(f))==f
        @test SymP(SymS(g))==g 
        @test SymP(SymC(g))==g
        @test SymC(SymP(h))==h 
        @test SymC(SymS(h))==h
    end

    @testset "GL_n tensor product" begin
        @test s_(2,1)*s_(3,1)==s_(3,2,1,1)+s_(3,2,2)+s_(3,3,1)+s_(4,1,1,1)+2*s_(4,2,1)+s_(4,3)+s_(5,1,1)+s_(5,2)
        @test p_(5)*p_(3)*p_(1)*p_(1)==p_(5,3,1,1)
        @test c_(2,2,1,1,1)*c_(2,2,1)==24*c_(2,2,2,2,1,1,1,1)
    end

    @testset "S_n tensor product" begin
        @test s_(3,2,1)&s_(4,1,1)==s_(2,1,1,1,1)+2*s_(2,2,1,1)+s_(2,2,2)+2*s_(3,1,1,1)+4*s_(3,2,1)+s_(3,3)+2*s_(4,1,1)+2*s_(4,2)+s_(5,1)
        @test g&p_(4,1)==-12*p_(4,1)
        @test h&c_(3,1)==-5c_(3,1) 
    end

    @testset "plethysm" begin
        @test plethysm(s_(4),s_(2))==s_(2,2,2,2)+s_(4,2,2)+s_(4,4)+s_(6,2)+s_(8)
        @test plethysm(p_(2,1),p_(3,1)+p_(2,2))==p_(4,4,2,2)+p_(4,4,3,1)+p_(6,2,2,2)+p_(6,3,2,1)
    end
end
        








