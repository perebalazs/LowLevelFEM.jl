cd(@__DIR__)
using SparseArrays
using StaticArrays
using LinearAlgebra
using Gtk
using BenchmarkTools
using Arpack
using NearestNeighbors
include("utilities.jl")


mesh  = Mesh("rect.geo")
mat1 = Material(E=208000)
sec1  = Solid2D(:PlaneStrain,mat1)
#disp1 = Disp((:x,0),:fixed)  
#load1 = Load((:l,:y,20),fy=-100,type=:distributed)

D = dofMatrix(mesh,[sec1])
@btime K = stiffnessMatrix(mesh,[sec1],D)

#kbc = constraint(mesh,K,disp1,D)
#dbc = loadVector(mesh,load1,D)
#q,t = statSolDirect(K,[dbc,kbc]) 
#addView(mesh,D,q,t,scf=200)
#startGMSH()




