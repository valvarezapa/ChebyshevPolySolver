# Last update: 26/03/24 ---------------------------------------------------------------------------# 
#--------------------------------------------------------------------------------------------------#
#------------------ CHEBYSHEV POLYNOMIAL NUMERICAL ROOT-FINDER MODULE -----------------------------#
#--------------------------------------------------------------------------------------------------#
#-------- This module contains algorithms for the specific case of finding the --------------------#
#-------- roots of a Chebyshev polynomial, also known as Chebyshev's nodes. -----------------------#
#-------- Since an analytic formula is known for this nodes, the purpose of this  -----------------#
#-------- module is restricted to numerical procedures, in order to test the ----------------------#
#-------- behaviour of the considered monodromy-based algorithms. ---------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

module ChebyshevPolySolver

include("geometry-lift.jl")
include("combinatorics-lift.jl")

using Polynomials
using Colors
using PyPlot


# Recursively constructs a list with the first n Chebyshev polynomials.
# Note that in the resulting list, the i-th element corresponds to the
# i-th Chebyshev polynomial; that is, the Chebyshev polynomial of degree i-1.
# Therefore, pol_list(10) is the Chebyshev polynomial of degree 9.
function Chebyshev_pol_list(n::Int64)::Array{Array{Complex{Float64},1},1}
    # Pre: n>=0
    # Post: Retuns the list of coefficients of every Chebyshev polynomial of degrees 0 up to n.
    pol_list=[[complex(1.0)]]
    if n>=1
        push!(pol_list,complex([0.0,1.0]))
        i=3
        while i in 3:n+1
            pol=vcat(complex([0.0]),2*pol_list[i-1])
            pol=pol-vcat(pol_list[i-2],complex([0.0,0.0]))
            push!(pol_list,pol)
            i+=1
        end
    end
    return pol_list
  end
  
  
#Example
#Chebyshev_pol_list(10)


# Normalizes a given polynomial, assuming that the norm of the polynomial
# is defined as the sum in absolute value of its coefficients.
function normalize_pol(coeffs::Array{Complex{Float64},1})::Array{Complex{Float64},1}
  # Pre: 'coeffs' does not represent the null polynomial.
  # Post: Retuns the list of coefficients of the normalized polynomial.
  norm=sum(abs,coeffs)
  return coeffs*(1.0/norm)
end


#Example
#=
pol=complex([1.0,0.0,0.0,2.0]) # f(z)=z^3-1
norm_pol=normalize_pol(pol)
=#


# 
function new_x(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64})::Complex{Float64}
  # Pre: The point 'x' is not a critical point of the given polynomial.
  # Post: Computes the next point.
  p=Polynomial(coeffs)
  der_p=derivative(p)       # IMPLEMENT A FASTER DERIVATIVE METHOD AND A FASTER EVALUATION METHOD
  @assert der_p(x)!=0
  return ((new_y-y)/der_p(x))+x
end


#Example
#x_new(Polynomial([1.0-3.0*im,complex(2.5),50.0+6.2*im]),y,new_y,x)


#
function new_precise_x(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64};precision::Float64=12.0)::Complex{Float64}
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post: Computes the last point.
  p=Polynomial(coeffs)
  #print("Until the image of the approximation is equal to ");println(new_y)
  #print("Point: ");print(x);print(" with image ");println(p(x))
  i=0
  while abs(p(new_x(coeffs,y,new_y,x))-new_y)>1.0/10^precision && i<500
    x=new_x(coeffs,y,new_y,x)
    #print("Point: ");print(x);print(" with image ");println(p(x))
    y=p(x)
    i=i+1
  end
  return new_x(coeffs,y,new_y,x)
end


#Example
#new_precise_x([1.0-3.0*im,complex(2.5),50.0+6.2*im],y,new_y,x;precision=8)


# 
function new_precise_x_list(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64};pre::Float64=12.0)::Array{Complex{Float64},1}
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post:
  return [x,new_precise_x(coeffs,y,new_y,x;precision=pre)]
end


#Example
#new_precise_x_list([1.0-3.0*im,complex(2.5),50.0+6.2*im],y,new_y,x;precision=8)


# 
function precise_lifting(coeffs::Array{Complex{Float64},1},y_points::Array{Complex{Float64},1},x::Complex{Float64};pre::Float64=12.0)::Complex{Float64}
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post:
  p=Polynomial(coeffs)
  l=length(y_points)
  @inbounds y_list=[y_points[i] for i in 1:l]
  @assert l>=2
  if l==2 
    return new_precise_x(coeffs,y_points[1],y_points[2],x;precision=pre)
  else
    y=y_list[1]
    new_y=y_list[2]
    new_x=new_precise_x(coeffs,y,new_y,x;precision=pre)
    #println(new_x)
    while length(y_list)>2
      popfirst!(y_list)
      y=y_list[1]
      new_y=y_list[2]
      #=
      println("start")
      println(y)
      println(new_y)
      =#
      new_x=new_precise_x(coeffs,y,new_y,new_x;precision=pre)
      #=
      println("end")
      println(new_x)
      =#
    end
    return new_precise_x(coeffs,y,new_y,new_x;precision=pre)
  end
end


#Example
#precise_lifting([1.0-3.0*im,complex(2.5),50.0+6.2*im],y,new_y,x;precision=8)


# 
function precise_lifting_list(coeffs::Array{Complex{Float64},1},y_points::Array{Complex{Float64},1},x::Complex{Float64};prec::Float64=12.0)::Array{Complex{Float64},1}
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post:
  p=Polynomial(coeffs)
  l=length(y_points)
  @inbounds y_list=[y_points[i] for i in 1:l]
  @assert l>=2
  if l==2 
    return new_precise_x_list(coeffs,y_points[1],y_points[2],x;pre=prec)
  else
    y=y_list[1]
    new_y=y_list[2]
    new_x_list=new_precise_x_list(coeffs,y,new_y,x;pre=prec)
    while length(y_list)>2
      popfirst!(y_list)
      y=y_list[1]
      new_y=y_list[2]
      new_x=new_precise_x(coeffs,y,new_y,last(new_x_list);precision=prec)
      push!(new_x_list,new_x)
    end
    push!(new_x_list,new_precise_x(coeffs,y,new_y,last(new_x_list);precision=prec))
    return new_x_list
  end
end


#Example
#precise_lifting([1.0-3.0*im,complex(2.5),50.0+6.2*im],y,new_y,x;precision=8)


# 
function precise_lifting_plot(coeffs::Array{Complex{Float64},1},y_points::Array{Complex{Float64},1},x::Complex{Float64};preci::Float64=12.0,plotcolor::AbstractString="blue")
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post: 
  p=Polynomial(coeffs)
  l=length(y_points)
  coveringlist=precise_lifting_list(coeffs,y_points,x;prec=preci)
  @inbounds for i in 1:l
    x1=real(coveringlist[i])
    y1=imag(coveringlist[i])
    x2=real(coveringlist[i+1])
    y2=imag(coveringlist[i+1])
    x_values=[x1,x2]
    y_values=[y1,y2]
    PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
  end
  @inbounds for i in 1:l
    x1=real(y_points[i])
    y1=imag(y_points[i])
    x2=real(y_points[i+1])
    y2=imag(y_points[i+1])
    x_values=[x1,x2]
    y_values=[y1,y2]
    PyPlot.plot(x_values,y_values,color="red",linestyle="-")
  end
end


#Example
#precise_lifting_plot(Polynomial(complex([1.0,0.0,1.0])),polygonal_path([1.0+1.0*im],complex(0.0),2.0+2.0*im),3.0*im)


# 
function new_taylor_x(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64})
  # Pre: The polynomial given by 'coeffs' is not zero.
  # Post: 
  degree=length(coeffs)-1
  p=Polynomial(coeffs)
  derivatives_p=[p]
  @inbounds for k in 1:degree
    push!(derivatives_p,derivative(last(derivatives_p)))
  end
  popfirst!(derivatives_p)
  @inbounds taylor_coeffs=[derivatives_p[k](x)/factorial(k) for k in 1:degree]
  @inbounds potential=[abs(taylor_coeffs[k])^(1.0/k) for k in 1:degree]
  @inbounds potenitialy=[abs(new_y-y)^(1.0/k) for k in 1:degree]
  @inbounds comparation=[potential[k]>potenitialy[k] for k in 1:degree]
  first=com.positions(comparation, true)[1]
  @assert abs(taylor_coeffs[first])!=0.0
  return ((new_y-y)/taylor_coeffs[first])^(1.0/first)+x
end


#Example
#new_taylor_x()


# 
function new_precise_taylor_x(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64};precision::Float64=8.0)
  # Pre: The polynomial given by 'coeffs' is not zero.
  # Post: 
  p=Polynomial(coeffs)
  while abs(p(new_taylor_x(coeffs,y,new_y,x))-new_y)>1.0/10^precision
    x=new_taylor_x(coeffs,y,new_y,x)
    y=p(x)
  end
  return new_taylor_x(coeffs,y,new_y,x)
end


#Example
#new_precise_taylor_x()


# Constructs a (discrete) path based on the origin that loops around the closed interval [-I,I], which contains the critical values -I and I of the considered rotated Chebyshev polynomial.
# Equivalently, constructs a representative of the homotopy class of loops around the infinity point, whose lift induces a transitive action over the fiber of an ordinary point,
# in particular over the fiber of the origin, which are the roots of the polynomial. This will be used to compute every root (starting at a given one).
function getLoop()::Array{Complex{Float64},1}
  # Pre: dist>2
  # Post:
  dist=3.0 # Constructs a rectangle of size dist x dist containing the closed interval [-1,1].
  num_points=200 # Preferibly take a number of points divisible by 4 and 10.
  s1=[t*dist/2.0 for t in 0:4.0/num_points:1] # 1/4 of the considered number of points accumulate in the beginning and ending of the loop.
  s2=[(1.0-t)*dist/2+t*(dist/2)*(1.0+im) for t in 0:10.0/num_points:1]
  s3=[(1.0-t)*(dist/2)*(1.0+im)+t*(dist/2)*(-1.0+im) for t in 0:10.0/num_points:1]
  s4=[(1.0-t)*(dist/2)*(-1.0+im)+t*(dist/2)*(-1.0-im) for t in 0:10.0/num_points:1]
  s5=[(1.0-t)*(dist/2)*(-1.0-im)+t*(dist/2)*(1.0-im) for t in 0:10.0/num_points:1]
  s6=[(1.0-t)*(dist/2)*(1.0-im)+t*(dist/2) for t in 0:10.0/num_points:1]
  return [s1;s2;s3;s4;s5;s6;reverse(s1)]
end


#Example
#=
path=getLoop()
geo.plot_polygonal_path(path)
gcf()
=#


# Constructs the fiber T_n^{-1}(0), where n is the degree of T_n and n>0.
function getFiber(coeffs::Array{Complex{Float64},1};preci::Float64=12.0)::Array{Float64,1}
  # Pre: length(coeffs)>1
  # Post:
  deg=length(coeffs)-1
  gen=getLoop()
  if coeffs[1]==0
    solution_list=[complex(0.0)]
  else
    p=Polynomial(coeffs)
    if deg<10
      init=-1.0+1.0*im
    elseif deg<100
      init=-0.1+0.1*im
    else
      init=-0.01+0.01*im
    end
    val=p(init)
    path=[(1-t)*val for t in 0:0.05:1]
    solution_list=[real(precise_lifting(coeffs,path,init;pre=preci))]
  end
  for i in 1:deg-1
    #=
    print("Lifting ");println(i)
    #print("    Coefs: ");println(coeffs)
    #print("    Gen: ");println(gen)
    print("    Init. point: ");println(solution_list[i])
    =#
    push!(solution_list,real(precise_lifting(coeffs,gen,complex(solution_list[i]);pre=preci)))
    #print("    Resulting list: ");println(solution_list)
  end    
  return solution_list
end


#Example
#


# 
function nsolve(coefs::Array{Complex{Float64},1};pre::Float64=12.0)::Array{Float64,1}
  # Pre: The polynomial is not constant.
  # Post:
  #=
  even_deg=coefs[1]!=complex(0.0)
  if even_deg==true # If the degree of the given Chebyshev polynomial is even, we multiply it by the identity polynomial in order to get 0 to be a root.
    coefs=vcat([0.0],coefs) # Note that if the degree of the polynomial is odd, 0 is a root and this product is not necessary.
  end
  =#
  coefs=im*coefs
  #coefs=normalize_pol(coefs)   # ELIMINAR (?)
  nodes=getFiber(coefs;preci=pre)
  #=
  if even_deg==true
    popfirst!(nodes)
  end
  =#
  return nodes
end


#Example
#nsolve()


# 
function nsolve(n::Int64;preci::Float64=12.0)::Array{Array{Float64,1},1}
    # Pre: n>0.
    # Post:
    sol_list=[[]]
    poly_list=Chebyshev_pol_list(n)
    for i in 2:n
      push!(sol_list,nsolve(poly_list[i];pre=preci))
      print("Roots of ");print(i);print(": ");println(last(sol_list))
    end
    return sol_list
end
  
  
#Example
#nsolve()


# 
function Chebyshev_nodes_exact(n::Int64)::Array{Float64,1}
  # Pre: n>0.
  # Post:
  nodes=[]
  for k in 0:n-1
    push!(nodes,cos(pi*(k+1.0/2.0)/n))
  end
  return nodes
end


#Example
#nsolve()


end